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

/// Element type of a stored 1-D array.
#[derive(Clone, Copy, PartialEq)]
enum Elem {
    F32,
    F64,
    I32,
    I64,
    U32,
    U64,
}
/// Location of every chunk of a chunked 1-D dataset compressed with deflate (optionally shuffled), so that chunks can be
/// read from the file and decompressed in parallel outside the HDF5 library, which serialises all calls in a process.
struct ChunkTable {
    elem: Elem,
    esize: usize,
    celems: usize,
    len: usize,
    // filter pipeline in write order: true = shuffle, false = deflate
    pipeline: Vec<bool>,
    id: hdf5_sys::h5i::hid_t,
    chunks: std::sync::Mutex<std::collections::HashMap<usize, (u64, u64, u32)>>, // file address, stored size, filter mask
}
impl ChunkTable {
    fn open(ds: &hdf5::Dataset) -> Option<ChunkTable> {
        use hdf5::filters::Filter;
        use hdf5::types::{FloatSize, IntSize, TypeDescriptor};
        let elem = match ds.dtype().ok()?.to_descriptor().ok()? {
            TypeDescriptor::Float(FloatSize::U4) => Elem::F32,
            TypeDescriptor::Float(FloatSize::U8) => Elem::F64,
            TypeDescriptor::Integer(IntSize::U4) => Elem::I32,
            TypeDescriptor::Integer(IntSize::U8) => Elem::I64,
            TypeDescriptor::Unsigned(IntSize::U4) => Elem::U32,
            TypeDescriptor::Unsigned(IntSize::U8) => Elem::U64,
            _ => return None,
        };
        let esize = ds.dtype().ok()?.size();
        let celems = *ds.chunk()?.first()?;
        let mut pipeline = Vec::new();
        for f in ds.filters() {
            match f {
                Filter::Shuffle => pipeline.push(true),
                Filter::Deflate(_) => pipeline.push(false),
                _ => return None,
            }
        }
        let len = ds.size();
        let t = ChunkTable {
            elem,
            esize,
            celems,
            len,
            pipeline,
            id: ds.id(),
            chunks: Default::default(),
        };
        t.locate(0, 1).ok()?;
        Some(t)
    }
    /// look up the file location of chunks k0..k1 that are not yet known (one HDF5 call per chunk, under the library lock)
    fn locate(&self, k0: usize, k1: usize) -> anyhow::Result<()> {
        let missing: Vec<usize> = {
            let c = self.chunks.lock().unwrap();
            (k0..k1).filter(|k| !c.contains_key(k)).collect()
        };
        if missing.is_empty() {
            return Ok(());
        }
        let found: Vec<(usize, (u64, u64, u32))> = hdf5::sync::sync(|| {
            missing
                .iter()
                .map(|&k| {
                    let off: [hdf5_sys::h5::hsize_t; 1] = [(k * self.celems) as _];
                    let (mut mask, mut addr, mut size) = (0u32, 0u64, 0u64);
                    let r = unsafe {
                        hdf5_sys::h5d::H5Dget_chunk_info_by_coord(
                            self.id,
                            off.as_ptr(),
                            &mut mask,
                            &mut addr,
                            &mut size,
                        )
                    };
                    if r < 0 || size == 0 {
                        anyhow::bail!("chunk {k} not found")
                    }
                    Ok((k, (addr, size, mask)))
                })
                .collect::<anyhow::Result<_>>()
        })?;
        self.chunks.lock().unwrap().extend(found);
        Ok(())
    }
    /// decoded bytes of chunk k
    fn chunk_bytes(&self, file: &std::fs::File, k: usize) -> anyhow::Result<Vec<u8>> {
        use std::os::unix::fs::FileExt;
        let (addr, size, mask) = self.chunks.lock().unwrap()[&k];
        let mut buf = vec![0u8; size as usize];
        file.read_exact_at(&mut buf, addr)?;
        let full = self.celems * self.esize;
        for (pos, &shuffle) in self.pipeline.iter().enumerate().rev() {
            if mask & (1 << pos) != 0 {
                continue;
            }
            buf = if shuffle {
                let n = buf.len() / self.esize;
                let mut out = vec![0u8; buf.len()];
                for b in 0..self.esize {
                    for i in 0..n {
                        out[i * self.esize + b] = buf[b * n + i];
                    }
                }
                out
            } else {
                let mut out = vec![0u8; full];
                let got = libdeflater::Decompressor::new()
                    .zlib_decompress(&buf, &mut out)
                    .map_err(|e| anyhow::anyhow!("deflate: {e:?}"))?;
                out.truncate(got);
                out
            };
        }
        Ok(buf)
    }
    /// elements lo..hi converted with `conv`, chunks decoded in parallel
    fn read<T: Send + Copy + Default>(
        &self,
        file: &std::fs::File,
        lo: usize,
        hi: usize,
        conv: impl Fn(&[u8], Elem) -> T + Sync,
    ) -> anyhow::Result<Vec<T>> {
        use rayon::prelude::*;
        let hi = hi.min(self.len);
        if lo >= hi {
            return Ok(Vec::new());
        }
        let (k0, k1) = (lo / self.celems, (hi - 1) / self.celems);
        self.locate(k0, k1 + 1)?;
        let parts: Vec<anyhow::Result<Vec<T>>> = (k0..=k1)
            .into_par_iter()
            .map(|k| {
                let bytes = self.chunk_bytes(file, k)?;
                let base = k * self.celems;
                let (a, b) = (lo.max(base) - base, hi.min(base + self.celems) - base);
                Ok((a..b)
                    .map(|i| conv(&bytes[i * self.esize..(i + 1) * self.esize], self.elem))
                    .collect())
            })
            .collect();
        let mut out = Vec::with_capacity(hi - lo);
        for p in parts {
            out.extend(p?);
        }
        Ok(out)
    }
}
fn to_f32(b: &[u8], e: Elem) -> f32 {
    match e {
        Elem::F32 => f32::from_le_bytes(b.try_into().unwrap()),
        Elem::F64 => f64::from_le_bytes(b.try_into().unwrap()) as f32,
        Elem::I32 => i32::from_le_bytes(b.try_into().unwrap()) as f32,
        Elem::I64 => i64::from_le_bytes(b.try_into().unwrap()) as f32,
        Elem::U32 => u32::from_le_bytes(b.try_into().unwrap()) as f32,
        Elem::U64 => u64::from_le_bytes(b.try_into().unwrap()) as f32,
    }
}
fn to_u32(b: &[u8], e: Elem) -> u32 {
    match e {
        Elem::I32 => i32::from_le_bytes(b.try_into().unwrap()) as u32,
        Elem::I64 => i64::from_le_bytes(b.try_into().unwrap()) as u32,
        Elem::U32 => u32::from_le_bytes(b.try_into().unwrap()),
        Elem::U64 => u64::from_le_bytes(b.try_into().unwrap()) as u32,
        Elem::F32 => f32::from_le_bytes(b.try_into().unwrap()) as u32,
        Elem::F64 => f64::from_le_bytes(b.try_into().unwrap()) as u32,
    }
}

/// Prefetching CSR reader. When the data and index arrays are chunked with deflate (and optionally shuffle), chunks are
/// located once and then read and decompressed in parallel by the rayon pool, bypassing the global HDF5 lock; the first
/// block is checked against the HDF5 library and any other storage layout falls back to serial reads through HDF5.
/// One block is kept in flight while the caller processes the previous one. LOY_SERIAL_READ=1 forces the HDF5 path.
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
            let serial = |lo: usize, hi: usize| -> anyhow::Result<(Vec<u32>, Vec<f32>)> {
                let v: Vec<f32> = data.read_slice_1d::<f32, _>(lo..hi)?.to_vec();
                let i: Vec<u32> = idx
                    .read_slice_1d::<i64, _>(lo..hi)?
                    .iter()
                    .map(|&x| x as u32)
                    .collect();
                Ok((i, v))
            };
            let raw = std::fs::File::open(&path)?;
            let mut fast = if std::env::var("LOY_SERIAL_READ").is_ok() {
                None
            } else {
                ChunkTable::open(&data).zip(ChunkTable::open(&idx))
            };
            if let Some((td, ti)) = &fast {
                // agreement check on the first elements
                let n = td.len.min(200_000);
                let (si, sv) = serial(0, n)?;
                let fv = td.read(&raw, 0, n, to_f32)?;
                let fi = ti.read(&raw, 0, n, to_u32)?;
                if fv != sv || fi != si {
                    eprintln!("prefetch_csr: parallel chunk reader disagrees with HDF5 on {path}; using serial reads");
                    fast = None;
                }
            }
            for (r0, r1) in ranges {
                let (lo, hi) = (indptr[r0] as usize, indptr[r1] as usize);
                let (i, v) = match &fast {
                    Some((td, ti)) => (ti.read(&raw, lo, hi, to_u32)?, td.read(&raw, lo, hi, to_f32)?),
                    None => serial(lo, hi)?,
                };
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
