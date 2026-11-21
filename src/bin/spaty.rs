//! spaty: detection-aware search for patches of cells without Y in Visium HD 8 um bins.
//! Per-cell Y counts are far below one expected read, so bins are pooled into square tiles of T x T bins and tested
//! for a deficit of MSY reads relative to their X gametologs, which escape X inactivation and keep their dosage when Y
//! is lost. Within a tile, Y | (Y + X) ~ beta-binomial(Y + X, p) with p the median tile ratio of the section and
//! overdispersion estimated by moments. Twenty autosomal gene-set pairs matched to the Y genes and X gametologs on
//! expression are processed identically as negative controls.
//! usage: spaty <filtered_feature_bc_matrix.h5> <tss.tsv: gene chrom pos> <tile T> <out_prefix> [min_n=20]
//! env: BG_YX = Y/X gametolog ratio of female sections (ambient and off-target background), default 0
use anyhow::{Context, Result};
use hdf5::types::{FixedAscii, VarLenAscii, VarLenUnicode};
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;
use std::collections::HashMap;
use std::io::{BufWriter, Write};

// TMSB4X/TMSB4Y excluded: TMSB4X is expressed two orders of magnitude above TMSB4Y and would dominate the denominator
const YG: &[&str] = &[
    "RPS4Y1", "ZFY", "USP9Y", "DDX3Y", "UTY", "KDM5D", "EIF1AY", "TXLNGY", "NLGN4Y", "PRKY",
];
const XG: &[&str] = &[
    "RPS4X", "ZFX", "USP9X", "DDX3X", "KDM6A", "KDM5C", "EIF1AX", "TXLNG", "NLGN4X", "PRKX",
];
const NP: usize = 20;
const NEAR: usize = 20; // controls drawn among the 20 autosomal genes nearest in total count
const CHUNK: usize = 20_000;

fn strings(ds: &hdf5::Dataset) -> Result<Vec<String>> {
    if let Ok(v) = ds.read_raw::<VarLenUnicode>() {
        return Ok(v.iter().map(|s| s.as_str().to_string()).collect());
    }
    if let Ok(v) = ds.read_raw::<VarLenAscii>() {
        return Ok(v.iter().map(|s| s.as_str().to_string()).collect());
    }
    let v = ds.read_raw::<FixedAscii<256>>().context("strings")?;
    Ok(v.iter().map(|s| s.as_str().to_string()).collect())
}

/// "s_008um_00301_00127-1" -> (301, 127)
fn coords(bc: &str) -> Option<(u32, u32)> {
    let p: Vec<&str> = bc.split('_').collect();
    if p.len() < 4 {
        return None;
    }
    Some((p[2].parse().ok()?, p[3].split('-').next()?.parse().ok()?))
}

struct SetStats {
    tiles: usize,
    p: f64,
    rho: f64,
    n_sig: usize,
    moran: f64,
}

/// tile-level beta-binomial deficit test and Moran's I of the z scores; returns stats and per-tile (z, q)
fn analyse(y: &[f64], x: &[f64], keep: &[bool], nbr: &[Vec<usize>]) -> (SetStats, Vec<f64>, Vec<f64>) {
    let idx: Vec<usize> = (0..y.len()).filter(|&i| keep[i] && y[i] + x[i] > 0.0).collect();
    let mut ratios: Vec<f64> = idx.iter().map(|&i| y[i] / (y[i] + x[i])).collect();
    ratios.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let p = if ratios.is_empty() {
        f64::NAN
    } else {
        ratios[ratios.len() / 2]
    };
    let (mut num, mut den) = (0.0, 0.0);
    for &i in &idx {
        let n = y[i] + x[i];
        num += (y[i] - n * p).powi(2) - n * p * (1.0 - p);
        den += n * (n - 1.0) * p * (1.0 - p);
    }
    let rho = if den > 0.0 { (num / den).max(0.0) } else { 0.0 };
    let mut z = vec![f64::NAN; y.len()];
    let mut pv: Vec<(usize, f64)> = Vec::new();
    for &i in &idx {
        let n = y[i] + x[i];
        let v = n * p * (1.0 - p) * (1.0 + (n - 1.0) * rho);
        let zi = (y[i] - n * p) / v.sqrt().max(1e-12);
        z[i] = zi;
        pv.push((i, 0.5 * erfc(-zi / std::f64::consts::SQRT_2))); // lower tail
    }
    // Benjamini-Hochberg
    let mut q = vec![f64::NAN; y.len()];
    pv.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
    let m = pv.len() as f64;
    let mut run = 1.0f64;
    for (r, &(i, pval)) in pv.iter().enumerate().rev() {
        run = run.min(pval * m / (r as f64 + 1.0));
        q[i] = run;
    }
    let n_sig = q.iter().filter(|v| v.is_finite() && **v < 0.05).count();
    // Moran's I with rook neighbours among tested tiles
    let zm = idx.iter().map(|&i| z[i]).sum::<f64>() / idx.len().max(1) as f64;
    let (mut s0, mut nm, mut dn) = (0.0, 0.0, 0.0);
    for &i in &idx {
        dn += (z[i] - zm).powi(2);
        for &j in &nbr[i] {
            if z[j].is_finite() {
                s0 += 1.0;
                nm += (z[i] - zm) * (z[j] - zm);
            }
        }
    }
    let moran = if s0 > 0.0 && dn > 0.0 {
        idx.len() as f64 / s0 * nm / dn
    } else {
        f64::NAN
    };
    (
        SetStats {
            tiles: idx.len(),
            p,
            rho,
            n_sig,
            moran,
        },
        z,
        q,
    )
}

fn erfc(x: f64) -> f64 {
    // Numerical Recipes erfc, relative error < 1.2e-7
    let z = x.abs();
    let t = 1.0 / (1.0 + 0.5 * z);
    let r = t
        * (-z * z - 1.26551223
            + t * (1.00002368
                + t * (0.37409196
                    + t * (0.09678418
                        + t * (-0.18628806
                            + t * (0.27886807
                                + t * (-1.13520398
                                    + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))))
            .exp();
    if x >= 0.0 {
        r
    } else {
        2.0 - r
    }
}

fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().collect();
    let (h5, tss, tsz, out) = (&a[1], &a[2], a[3].parse::<u32>()?, &a[4]);
    let min_n: f64 = a.get(5).map(|s| s.parse()).transpose()?.unwrap_or(20.0);
    let bg: f64 = std::env::var("BG_YX")
        .ok()
        .map(|s| s.parse())
        .transpose()?
        .unwrap_or(0.0);
    let chrom: HashMap<String, String> = std::fs::read_to_string(tss)?
        .lines()
        .filter_map(|l| {
            let x: Vec<&str> = l.split('\t').collect();
            (x.len() >= 2).then(|| (x[0].to_string(), x[1].to_string()))
        })
        .collect();
    let f = hdf5::File::open(h5).with_context(|| h5.clone())?;
    let m = f.group("matrix")?;
    let names = strings(&m.dataset("features/name")?)?;
    let bcs = strings(&m.dataset("barcodes")?)?;
    let indptr: Vec<i64> = m.dataset("indptr")?.read_raw()?;
    let (data, indices) = (m.dataset("data")?, m.dataset("indices")?);
    let (ng, nb) = (names.len(), bcs.len());
    // pass 1: gene totals and bin totals
    let mut gtot = vec![0f64; ng];
    let mut btot = vec![0f64; nb];
    for c0 in (0..nb).step_by(CHUNK) {
        let c1 = (c0 + CHUNK).min(nb);
        let (lo, hi) = (indptr[c0] as usize, indptr[c1] as usize);
        let d: Vec<i32> = data.read_slice_1d::<i32, _>(lo..hi)?.to_vec();
        let ix: Vec<i64> = indices.read_slice_1d::<i64, _>(lo..hi)?.to_vec();
        for c in c0..c1 {
            for k in (indptr[c] as usize - lo)..(indptr[c + 1] as usize - lo) {
                gtot[ix[k] as usize] += d[k] as f64;
                btot[c] += d[k] as f64;
            }
        }
    }
    let gi: HashMap<&str, usize> = names.iter().enumerate().map(|(i, n)| (n.as_str(), i)).collect();
    // keep gametolog pairs whose Y and X members are both on the probe set
    let pairs: Vec<(usize, usize)> = YG
        .iter()
        .zip(XG)
        .filter_map(|(y, x)| Some((*gi.get(y)?, *gi.get(x)?)))
        .collect();
    eprintln!(
        "{h5}: {nb} bins, {ng} genes; gametolog pairs on panel: {}",
        pairs.len()
    );
    // expression bins of autosomal genes for matched controls
    let auto: Vec<usize> = (0..ng)
        .filter(|&g| {
            gtot[g] > 0.0
                && chrom
                    .get(&names[g])
                    .is_some_and(|c| c != "chrX" && c != "chrY" && c != "chrM")
        })
        .collect();
    let mut auto = auto;
    auto.sort_by(|&p, &q| gtot[p].partial_cmp(&gtot[q]).unwrap());
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(7);
    let mut pick = |g: usize| -> usize {
        let pos = auto.partition_point(|&h| gtot[h] < gtot[g]);
        let lo = pos.saturating_sub(NEAR / 2).min(auto.len().saturating_sub(NEAR));
        auto[lo + rng.gen_range(0..NEAR.min(auto.len()))]
    };
    // set membership: gene -> list of (set id, side) ; set 0 = Y/X, sets 1..=NP controls; side 0 = numerator, 1 = denominator
    let mut member: HashMap<usize, Vec<(usize, usize)>> = HashMap::new();
    for &(yg, xg) in &pairs {
        member.entry(yg).or_default().push((0, 0));
        member.entry(xg).or_default().push((0, 1));
    }
    for s in 1..=NP {
        for &(yg, xg) in &pairs {
            member.entry(pick(yg)).or_default().push((s, 0));
            member.entry(pick(xg)).or_default().push((s, 1));
        }
    }
    // tiles
    let mut tid: HashMap<(u32, u32), usize> = HashMap::new();
    let mut tkey: Vec<(u32, u32)> = Vec::new();
    let mut btile = vec![usize::MAX; nb];
    for (c, b) in bcs.iter().enumerate() {
        if let Some((r, col)) = coords(b) {
            let k = (r / tsz, col / tsz);
            let n = tkey.len();
            let t = *tid.entry(k).or_insert_with(|| {
                tkey.push(k);
                n
            });
            btile[c] = t;
        }
    }
    let nt = tkey.len();
    let mut cnt = vec![[0f64; 2]; nt * (NP + 1)];
    let mut tn = vec![0f64; nt];
    let mut tbins = vec![0u32; nt];
    for c in 0..nb {
        if btile[c] != usize::MAX {
            tn[btile[c]] += btot[c];
            tbins[btile[c]] += 1;
        }
    }
    // pass 2: set counts per tile
    for c0 in (0..nb).step_by(CHUNK) {
        let c1 = (c0 + CHUNK).min(nb);
        let (lo, hi) = (indptr[c0] as usize, indptr[c1] as usize);
        let d: Vec<i32> = data.read_slice_1d::<i32, _>(lo..hi)?.to_vec();
        let ix: Vec<i64> = indices.read_slice_1d::<i64, _>(lo..hi)?.to_vec();
        for c in c0..c1 {
            let t = btile[c];
            if t == usize::MAX {
                continue;
            }
            for k in (indptr[c] as usize - lo)..(indptr[c + 1] as usize - lo) {
                if let Some(v) = member.get(&(ix[k] as usize)) {
                    for &(s, side) in v {
                        cnt[s * nt + t][side] += d[k] as f64;
                    }
                }
            }
        }
    }
    let mut nbr: Vec<Vec<usize>> = vec![Vec::new(); nt];
    for (t, &(r, c)) in tkey.iter().enumerate() {
        for (dr, dc) in [(0i64, 1i64), (1, 0), (0, -1), (-1, 0)] {
            let (rr, cc) = (r as i64 + dr, c as i64 + dc);
            if rr >= 0 && cc >= 0 {
                if let Some(&u) = tid.get(&(rr as u32, cc as u32)) {
                    nbr[t].push(u);
                }
            }
        }
    }
    let mut summ = BufWriter::new(std::fs::File::create(format!("{out}.summary.tsv"))?);
    writeln!(
        summ,
        "set\ttiles\tp_ref\trho\tn_sig_q05\tmoran_I\tsum_num\tsum_den"
    )?;
    let mut y_res = None;
    for s in 0..=NP {
        let y: Vec<f64> = (0..nt).map(|t| cnt[s * nt + t][0]).collect();
        let x: Vec<f64> = (0..nt).map(|t| cnt[s * nt + t][1]).collect();
        let keep: Vec<bool> = (0..nt).map(|t| y[t] + x[t] >= min_n).collect();
        let (st, z, q) = analyse(&y, &x, &keep, &nbr);
        writeln!(
            summ,
            "{}\t{}\t{:.5}\t{:.5}\t{}\t{:.4}\t{}\t{}",
            if s == 0 {
                "Y".to_string()
            } else {
                format!("control{s}")
            },
            st.tiles,
            st.p,
            st.rho,
            st.n_sig,
            st.moran,
            y.iter().sum::<f64>(),
            x.iter().sum::<f64>()
        )?;
        if s == 0 {
            y_res = Some((y, x, z, q, st.p));
        }
    }
    let (y, x, z, q, p) = y_res.unwrap();
    let r_ref = p / (1.0 - p);
    let mut w = BufWriter::new(std::fs::File::create(format!("{out}.tiles.tsv"))?);
    writeln!(w, "tile_row\ttile_col\tbins\tumi\ty\tx\tf_loss\tz\tq")?;
    for t in 0..nt {
        let f_loss = if x[t] > 0.0 {
            1.0 - (y[t] / x[t] - bg) / (r_ref - bg)
        } else {
            f64::NAN
        };
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.3}\t{:.4}",
            tkey[t].0, tkey[t].1, tbins[t], tn[t], y[t], x[t], f_loss, z[t], q[t]
        )?;
    }
    let (sy, sx, sn) = (
        y.iter().sum::<f64>(),
        x.iter().sum::<f64>(),
        btot.iter().sum::<f64>(),
    );
    eprintln!(
        "{out}: {nt} tiles of {tsz}x{tsz} bins; Y per UMI {:.3e}, X gametologs per UMI {:.3e}, Y/X {:.4}",
        sy / sn,
        sx / sn,
        sy / sx.max(1.0)
    );
    Ok(())
}
