//! loyfit: per group (dataset x cell type x donor ...), zero-inflated negative binomial fit of Y-set counts given
//! total UMI: y ~ f*delta0 + (1-f)*NB(mu*N, k), EM with k chosen on a grid. Same fit for 20 pseudo-Y sets.
//! usage: loyfit <prefix> [min_cells=50] [n_boot=200]   (reads <prefix>.groups.tsv, .cells.bin, .info.txt)
//! writes <prefix>.fit.tsv
use anyhow::Result;
use loyscope::fit_amb;
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;
use rayon::prelude::*;
use std::io::Write;

const KGRID: &[f64] = &[0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 64.0, 1e6];

fn ln_gamma(x: f64) -> f64 {
    // Lanczos g=7
    const C: [f64; 9] = [
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
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
fn em(n: &[f64], y: &[f64], k: f64) -> (f64, f64, f64) {
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
fn fit(n: &[f64], y: &[f64]) -> (f64, f64, f64) {
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
fn median(v: &mut Vec<f64>) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[v.len() / 2]
}

fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().collect();
    let pre = &a[1];
    let minc: usize = a.get(2).map(|x| x.parse().unwrap()).unwrap_or(50);
    let nboot: usize = a.get(3).map(|x| x.parse().unwrap()).unwrap_or(200);
    let info = std::fs::read_to_string(format!("{pre}.info.txt"))?;
    let rl: usize = info
        .lines()
        .find(|l| l.starts_with("record_u32"))
        .unwrap()
        .split(' ')
        .nth(1)
        .unwrap()
        .parse()?;
    let raw = std::fs::read(format!("{pre}.cells.bin"))?;
    let v: Vec<u32> = raw
        .chunks_exact(4)
        .map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]]))
        .collect();
    let nc = v.len() / rl;
    let np = rl - 3;
    let gl: Vec<Vec<String>> = std::fs::read_to_string(format!("{pre}.groups.tsv"))?
        .lines()
        .skip(1)
        .map(|l| l.split('\t').map(String::from).collect())
        .collect();
    let mut cells: Vec<Vec<usize>> = vec![Vec::new(); gl.len()];
    for c in 0..nc {
        if v[c * rl + 1] > 0 {
            cells[v[c * rl] as usize].push(c);
        }
    }
    let todo: Vec<usize> = (0..gl.len())
        .filter(|&g| {
            cells[g].len() >= minc
                && (std::env::var("LOY_ALL_DISEASE").is_ok() || gl[g][4] == "normal")
                && (gl[g][3] == "male" || gl[g][3] == "female")
        })
        .collect();
    // ambient (LOY_AMBIENT=1): rho = female Y rate per UMI / median male intact rate from a first-pass fit
    let ambient = std::env::var("LOY_AMBIENT").is_ok();
    let mut rho = 0.0;
    if ambient {
        let (mut fy, mut fnn) = (0f64, 0f64);
        for c in 0..nc {
            if gl[v[c * rl] as usize][3] == "female" {
                fy += v[c * rl + 2] as f64;
                fnn += v[c * rl + 1] as f64;
            }
        }
        let mut mus: Vec<f64> = todo
            .par_iter()
            .filter(|&&g| gl[g][3] == "male")
            .filter_map(|&g| {
                let cs = &cells[g];
                let n: Vec<f64> = cs.iter().map(|&c| v[c * rl + 1] as f64).collect();
                let y: Vec<f64> = cs.iter().map(|&c| v[c * rl + 2] as f64).collect();
                let (_, mu, k) = fit(&n, &y);
                if k.is_nan() {
                    None
                } else {
                    Some(mu)
                }
            })
            .collect();
        mus.sort_by(|a, b| a.partial_cmp(b).unwrap());
        if !mus.is_empty() && fnn > 0.0 {
            rho = (fy / fnn) / mus[mus.len() / 2];
        }
        eprintln!("ambient rho {rho:.4e} (female rate {:.3e})", fy / fnn.max(1.0));
    }
    let fitx = |n: &[f64], y: &[f64]| -> (f64, f64, f64) {
        if ambient {
            fit_amb(n, y, rho)
        } else {
            fit(n, y)
        }
    };
    let res: Vec<String> = todo.par_iter().map(|&g| {
        let cs = &cells[g];
        let n: Vec<f64> = cs.iter().map(|&c| v[c * rl + 1] as f64).collect();
        let col = |j: usize| -> Vec<f64> { cs.iter().map(|&c| v[c * rl + 2 + j] as f64).collect() };
        let y = col(0);
        let naive = y.iter().filter(|&&x| x == 0.0).count() as f64 / y.len() as f64;
        let (f, mu, k) = fitx(&n, &y);
        let male = gl[g][3] == "male";
        let (mut lo, mut hi) = (f64::NAN, f64::NAN);
        let (mut pf_med, mut pf_max, mut pn_med) = (f64::NAN, f64::NAN, f64::NAN);
        if male {
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(g as u64 * 104729 + 1);
            let kk = if k.is_nan() { 1.0 } else { k };
            let mut bs: Vec<f64> = (0..nboot).map(|_| {
                let idx: Vec<usize> = (0..n.len()).map(|_| rng.gen_range(0..n.len())).collect();
                let nn: Vec<f64> = idx.iter().map(|&i| n[i]).collect(); let yy: Vec<f64> = idx.iter().map(|&i| y[i]).collect();
                if yy.iter().all(|&x| x == 0.0) { 1.0 } else if ambient { loyscope::em_amb(&nn, &yy, kk, rho).0 } else { em(&nn, &yy, kk).0 }
            }).collect();
            bs.sort_by(|a, b| a.partial_cmp(b).unwrap()); lo = bs[(nboot as f64 * 0.025) as usize]; hi = bs[((nboot as f64 * 0.975) as usize).min(nboot - 1)];
            let mut pf: Vec<f64> = Vec::new(); let mut pn: Vec<f64> = Vec::new();
            for j in 1..=np { let p = col(j); pn.push(p.iter().filter(|&&x| x == 0.0).count() as f64 / p.len() as f64); pf.push(fitx(&n, &p).0); }
            pf_max = pf.iter().cloned().fold(f64::NEG_INFINITY, f64::max); pf_med = median(&mut pf); pn_med = median(&mut pn);
        }
        let mean_n = n.iter().sum::<f64>() / n.len() as f64;
        format!("{}\t{}\t{mean_n:.0}\t{naive:.4}\t{f:.4}\t{lo:.4}\t{hi:.4}\t{mu:.3e}\t{k}\t{pn_med:.4}\t{pf_med:.4}\t{pf_max:.4}", gl[g][..9].join("\t"), cs.len())
    }).collect();
    let mut w = std::fs::File::create(format!("{pre}.{}.tsv", if ambient { "fitamb" } else { "fit" }))?;
    writeln!(w, "gid\tcell_type\ttissue\tsex\tdisease\tdonor\tstage\tassay\tsuspension\tn\tmean_umi\tnaive_zero\tf\tf_lo\tf_hi\tmu\tk\tpseudo_naive_med\tpseudo_f_med\tpseudo_f_max")?;
    for r in res {
        writeln!(w, "{r}")?;
    }
    eprintln!("{pre}: fitted {} groups", todo.len());
    Ok(())
}
