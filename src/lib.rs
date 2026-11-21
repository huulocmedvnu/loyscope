//! shared ZINB detection-efficiency estimator (same code as loyfit)
pub const KGRID: &[f64] = &[0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 64.0, 1e6];

pub fn ln_gamma(x: f64) -> f64 {
    // Lanczos g=7
    const C: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.5203681218851,
        -1259.1392167224028,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507343278686905,
        -0.13857109526572012,
        9.984_369_578_019_572e-6,
        1.5056327351493116e-7,
    ];
    if x < 0.5 {
        return (std::f64::consts::PI / (std::f64::consts::PI * x).sin()).ln() - ln_gamma(1.0 - x);
    }
    let x = x - 1.0;
    let mut s = C[0];
    for i in 1..9 {
        s += C[i] / (x + i as f64);
    }
    let t = x + 7.5;
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + s.ln()
}
/// EM for fixed k. returns (f, mu, loglik)
pub fn em(n: &[f64], y: &[f64], k: f64) -> (f64, f64, f64) {
    let m = n.len() as f64;
    let mut f = (y.iter().filter(|&&v| v == 0.0).count() as f64 / m * 0.5).max(1e-4);
    let mut mu = y.iter().sum::<f64>() / n.iter().sum::<f64>().max(1.0) / (1.0 - f);
    if mu <= 0.0 {
        return (1.0, 0.0, 0.0);
    }
    let mut tau = vec![0f64; n.len()];
    for _ in 0..500 {
        let fo = f;
        let muo = mu;
        let mut st = 0.0;
        for i in 0..n.len() {
            if y[i] == 0.0 {
                let p0 = (k / (k + mu * n[i])).powf(k);
                tau[i] = f / (f + (1.0 - f) * p0);
                st += tau[i];
            } else {
                tau[i] = 0.0;
            }
        }
        f = (st / m).clamp(0.0, 1.0);
        for _ in 0..5 {
            // fixed-point for mu
            let (mut a, mut b) = (0.0, 0.0);
            for i in 0..n.len() {
                let w = 1.0 - tau[i];
                a += w * y[i];
                b += w * (y[i] + k) * n[i] / (k + mu * n[i]);
            }
            if b > 0.0 {
                mu = a / b;
            }
        }
        if (f - fo).abs() < 1e-7 && ((mu - muo) / muo).abs() < 1e-7 {
            break;
        }
    }
    let lgk = ln_gamma(k);
    let ll: f64 = (0..n.len())
        .map(|i| {
            let mn = mu * n[i];
            let lnb = ln_gamma(y[i] + k) - lgk - ln_gamma(y[i] + 1.0)
                + k * (k / (k + mn)).ln()
                + y[i] * (mn / (k + mn)).ln();
            if y[i] == 0.0 {
                (f + (1.0 - f) * lnb.exp()).ln()
            } else {
                (1.0 - f).ln() + lnb
            }
        })
        .sum();
    (f, mu, ll)
}
pub fn fit(n: &[f64], y: &[f64]) -> (f64, f64, f64) {
    // (f, mu, k)
    if y.iter().all(|&v| v == 0.0) {
        return (1.0, 0.0, f64::NAN);
    }
    let mut best = (f64::NAN, 0.0, 0.0, f64::NEG_INFINITY);
    for &k in KGRID {
        let (f, mu, ll) = em(n, y, k);
        if ll > best.3 {
            best = (f, mu, k, ll);
        }
    }
    (best.0, best.1, best.2)
}

/// EM with an ambient floor: LOY cells y ~ Pois(rho * mu * N), intact y ~ NB(mu N, k). Returns (f, mu, loglik)
pub fn em_amb(n: &[f64], y: &[f64], k: f64, rho: f64) -> (f64, f64, f64) {
    if rho <= 0.0 {
        return em(n, y, k);
    }
    let m = n.len() as f64;
    let mut f = (y.iter().filter(|&&v| v == 0.0).count() as f64 / m * 0.5).max(1e-4);
    let mut mu = y.iter().sum::<f64>() / n.iter().sum::<f64>().max(1.0) / (1.0 - f);
    if mu <= 0.0 {
        return (1.0, 0.0, 0.0);
    }
    let lgk = ln_gamma(k);
    let lg1: Vec<f64> = y.iter().map(|&v| ln_gamma(v + 1.0)).collect();
    let lgyk: Vec<f64> = y.iter().map(|&v| ln_gamma(v + k)).collect();
    let comp = |mu: f64, i: usize| -> (f64, f64) {
        let lam = rho * mu * n[i];
        let lp = y[i] * lam.ln() - lam - lg1[i];
        let mn = mu * n[i];
        let lnb = lgyk[i] - lgk - lg1[i] + k * (k / (k + mn)).ln() + y[i] * (mn / (k + mn)).ln();
        (lp, lnb)
    };
    let mut tau = vec![0f64; n.len()];
    for _ in 0..500 {
        let (fo, muo) = (f, mu);
        let mut st = 0.0;
        for i in 0..n.len() {
            let (lp, lnb) = comp(mu, i);
            let a1 = f.ln() + lp;
            let a0 = (1.0 - f).ln() + lnb;
            let mx = a1.max(a0);
            tau[i] = (a1 - mx).exp() / ((a1 - mx).exp() + (a0 - mx).exp());
            st += tau[i];
        }
        f = (st / m).clamp(1e-9, 1.0 - 1e-9);
        for _ in 0..5 {
            // joint fixed point: NB part (weight 1-tau) + Poisson part (weight tau, rate rho*mu*N)
            let (mut a, mut b) = (0.0, 0.0);
            for i in 0..n.len() {
                let w = 1.0 - tau[i];
                a += w * y[i] + tau[i] * y[i];
                b += w * (y[i] + k) * n[i] / (k + mu * n[i]) + tau[i] * rho * n[i];
            }
            if b > 0.0 {
                mu = a / b;
            }
        }
        if (f - fo).abs() < 1e-7 && ((mu - muo) / muo).abs() < 1e-7 {
            break;
        }
    }
    let ll: f64 = (0..n.len())
        .map(|i| {
            let (lp, lnb) = comp(mu, i);
            let a1 = f.ln() + lp;
            let a0 = (1.0 - f).ln() + lnb;
            let mx = a1.max(a0);
            mx + ((a1 - mx).exp() + (a0 - mx).exp()).ln()
        })
        .sum();
    (f, mu, ll)
}
pub fn fit_amb(n: &[f64], y: &[f64], rho: f64) -> (f64, f64, f64) {
    if y.iter().all(|&v| v == 0.0) {
        return (1.0, 0.0, f64::NAN);
    }
    let mut best = (f64::NAN, 0.0, 0.0, f64::NEG_INFINITY);
    for &k in KGRID {
        let (f, mu, ll) = em_amb(n, y, k, rho);
        if ll > best.3 {
            best = (f, mu, k, ll);
        }
    }
    (best.0, best.1, best.2)
}

/// EM with an absolute ambient rate: cells without Y have y ~ Pois(bg * N), intact cells y ~ NB(mu N, k).
/// `bg` is the ambient rate per UMI of the gene set (contamination fraction times the share of the set in the soup),
/// fixed per dataset and independent of the expression rate of the group. Returns (f, mu, loglik).
pub fn em_bg(n: &[f64], y: &[f64], k: f64, bg: f64) -> (f64, f64, f64) {
    if bg <= 0.0 {
        return em(n, y, k);
    }
    let m = n.len() as f64;
    let mut f = (y.iter().filter(|&&v| v == 0.0).count() as f64 / m * 0.5).max(1e-4);
    let mut mu = y.iter().sum::<f64>() / n.iter().sum::<f64>().max(1.0) / (1.0 - f);
    if mu <= 0.0 {
        return (1.0, 0.0, 0.0);
    }
    let lgk = ln_gamma(k);
    let lg1: Vec<f64> = y.iter().map(|&v| ln_gamma(v + 1.0)).collect();
    let lgyk: Vec<f64> = y.iter().map(|&v| ln_gamma(v + k)).collect();
    let lp: Vec<f64> = (0..n.len())
        .map(|i| {
            let lam = bg * n[i];
            y[i] * lam.ln() - lam - lg1[i]
        })
        .collect();
    let lnb = |mu: f64, i: usize| -> f64 {
        let mn = mu * n[i];
        lgyk[i] - lgk - lg1[i] + k * (k / (k + mn)).ln() + y[i] * (mn / (k + mn)).ln()
    };
    let mut tau = vec![0f64; n.len()];
    for _ in 0..500 {
        let (fo, muo) = (f, mu);
        let mut st = 0.0;
        for i in 0..n.len() {
            let a1 = f.ln() + lp[i];
            let a0 = (1.0 - f).ln() + lnb(mu, i);
            let mx = a1.max(a0);
            tau[i] = (a1 - mx).exp() / ((a1 - mx).exp() + (a0 - mx).exp());
            st += tau[i];
        }
        f = (st / m).clamp(1e-9, 1.0 - 1e-9);
        for _ in 0..5 {
            // NB part only: the ambient rate does not depend on mu
            let (mut a, mut b) = (0.0, 0.0);
            for i in 0..n.len() {
                let w = 1.0 - tau[i];
                a += w * y[i];
                b += w * (y[i] + k) * n[i] / (k + mu * n[i]);
            }
            if b > 0.0 && a > 0.0 {
                mu = a / b;
            }
        }
        if (f - fo).abs() < 1e-7 && ((mu - muo) / muo).abs() < 1e-7 {
            break;
        }
    }
    let ll: f64 = (0..n.len())
        .map(|i| {
            let a1 = f.ln() + lp[i];
            let a0 = (1.0 - f).ln() + lnb(mu, i);
            let mx = a1.max(a0);
            mx + ((a1 - mx).exp() + (a0 - mx).exp()).ln()
        })
        .sum();
    (f, mu, ll)
}
pub fn fit_bg(n: &[f64], y: &[f64], bg: f64) -> (f64, f64, f64) {
    if y.iter().all(|&v| v == 0.0) {
        return (1.0, 0.0, f64::NAN);
    }
    let mut best = (f64::NAN, 0.0, 0.0, f64::NEG_INFINITY);
    for &k in KGRID {
        let (f, mu, ll) = em_bg(n, y, k, bg);
        if ll > best.3 {
            best = (f, mu, k, ll);
        }
    }
    (best.0, best.1, best.2)
}

/// A block of CSR rows: (first row, local indptr, column indices, values).
pub type Chunk = (usize, Vec<usize>, Vec<u32>, Vec<f32>);
/// Prefetching CSR reader: a dedicated thread reads and decompresses the requested row ranges from the HDF5 file
/// (HDF5 itself is serialized), while the caller processes the previous block in parallel. Two blocks are kept in flight.
pub fn prefetch_csr(
    path: &str,
    mpath: &str,
    indptr: std::sync::Arc<Vec<i64>>,
    ranges: Vec<(usize, usize)>,
) -> std::sync::mpsc::Receiver<anyhow::Result<Chunk>> {
    let (tx, rx) = std::sync::mpsc::sync_channel(2);
    let (path, mpath) = (path.to_string(), mpath.to_string());
    std::thread::spawn(move || {
        let run = || -> anyhow::Result<()> {
            let f = hdf5::File::open(&path)?;
            let m = f.group(&mpath)?;
            let (data, idx) = (m.dataset("data")?, m.dataset("indices")?);
            for (r0, r1) in ranges {
                let (lo, hi) = (indptr[r0] as usize, indptr[r1] as usize);
                let v: Vec<f32> = data.read_slice_1d::<f32, _>(lo..hi)?.to_vec();
                let i: Vec<u32> = idx
                    .read_slice_1d::<i64, _>(lo..hi)?
                    .iter()
                    .map(|&x| x as u32)
                    .collect();
                let p: Vec<usize> = (r0..=r1).map(|r| indptr[r] as usize - lo).collect();
                if tx.send(Ok((r0, p, i, v))).is_err() {
                    break;
                }
            }
            Ok(())
        };
        if let Err(e) = run() {
            let _ = tx.send(Err(e));
        }
    });
    rx
}
pub fn row_ranges(nc: usize, chunk: usize, every: usize) -> Vec<(usize, usize)> {
    (0..nc)
        .step_by(chunk)
        .enumerate()
        .filter(|(i, _)| i % every.max(1) == 0)
        .map(|(_, r0)| (r0, (r0 + chunk).min(nc)))
        .collect()
}
