//! loyde: genome-wide expression of LOY vs Y+ cells from the saved pseudobulk (male normal groups, cells split by
//! zero vs non-zero Y). Only groups where zero-Y cells are mostly true LOY are used: f / naive_zero >= purity
//! (f from the ambient-aware fit). Per group: log2 CP10k ratio (y0 / ypos) per gene; per cell type: median over groups,
//! share positive, two-sided sign-test P, BH q over genes.
//! usage: loyde <prefix> [purity=0.8] [min_y0_cells=10] [min_cp10k=1]   needs .pb.f32.gz .pb_slots.tsv .genes.txt .groups.tsv .fitamb.tsv
//! writes <prefix>.de.tsv
use anyhow::Result;
use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::{Read, Write};

fn sign_p(k: usize, n: usize) -> f64 { // two-sided binomial(0.5)
    let lc = |n: usize, k: usize| -> f64 { (1..=k).map(|i| ((n - k + i) as f64 / i as f64).ln()).sum() };
    let m = k.min(n - k); let mut s = 0.0; for i in 0..=m { s += (lc(n, i) - n as f64 * 2f64.ln()).exp(); } (2.0 * s).min(1.0)
}
fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().collect();
    let pre = &a[1];
    let purity: f64 = a.get(2).map(|x| x.parse().unwrap()).unwrap_or(0.8);
    let min_y0: f64 = a.get(3).map(|x| x.parse().unwrap()).unwrap_or(10.0);
    let min_cp: f64 = a.get(4).map(|x| x.parse().unwrap()).unwrap_or(1.0);
    let genes: Vec<String> = std::fs::read_to_string(format!("{pre}.genes.txt"))?.lines().map(String::from).collect();
    let ng = genes.len();
    let gl: Vec<Vec<String>> = std::fs::read_to_string(format!("{pre}.groups.tsv"))?.lines().skip(1).map(|l| l.split('\t').map(String::from).collect()).collect();
    let mut fitf: HashMap<usize, (f64, f64)> = HashMap::new();
    for l in std::fs::read_to_string(format!("{pre}.fitamb.tsv"))?.lines().skip(1) { let f: Vec<&str> = l.split('\t').collect(); if f[16] != "NaN" { fitf.insert(f[0].parse()?, (f[12].parse()?, f[11].parse()?)); } }
    let slots: Vec<(usize, usize, f64, f64, f64, f64)> = std::fs::read_to_string(format!("{pre}.pb_slots.tsv"))?.lines().skip(1).map(|l| { let f: Vec<&str> = l.split('\t').collect();
        (f[0].parse().unwrap(), f[1].parse().unwrap(), f[2].parse().unwrap(), f[3].parse().unwrap(), f[4].parse().unwrap(), f[5].parse().unwrap()) }).collect();
    let mut raw = Vec::new(); MultiGzDecoder::new(std::fs::File::open(format!("{pre}.pb.f32.gz"))?).read_to_end(&mut raw)?;
    let pb: Vec<f32> = raw.chunks_exact(4).map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect(); drop(raw);
    // usable slots
    let use_: Vec<&(usize, usize, f64, f64, f64, f64)> = slots.iter().filter(|s| { let Some(&(f, nz)) = fitf.get(&s.1) else { return false }; s.2 >= min_y0 && s.3 >= 20.0 && nz > 0.0 && f / nz >= purity }).collect();
    eprintln!("{pre}: slots {} usable {}", slots.len(), use_.len());
    let mut by_ct: HashMap<String, Vec<usize>> = HashMap::new();
    for (i, s) in use_.iter().enumerate() { by_ct.entry(gl[s.1][1].clone()).or_default().push(i); }
    // per usable slot log2 ratio per gene (NaN if both low)
    let lfc: Vec<Vec<f32>> = use_.par_iter().map(|s| { let (sl, _, _, _, u0, u1) = **s; let b0 = sl * 2 * ng; let b1 = b0 + ng;
        (0..ng).map(|g| { let (c0, c1) = (pb[b0 + g] as f64 / u0 * 1e4, pb[b1 + g] as f64 / u1 * 1e4); if c1 < min_cp { f32::NAN } else { ((c0 + 0.05) / (c1 + 0.05)).log2() as f32 } }).collect() }).collect();
    let mut w = std::fs::File::create(format!("{pre}.de.tsv"))?;
    writeln!(w, "cell_type\tgene\tn_groups\tmedian_log2fc\tshare_up\tsign_p\tq_bh")?;
    let mut cts: Vec<&String> = by_ct.keys().collect(); cts.sort();
    for ct in cts {
        let idx = &by_ct[ct]; if idx.len() < 5 { continue; }
        let mut res: Vec<(usize, usize, f64, f64, f64)> = (0..ng).into_par_iter().filter_map(|g| {
            let mut v: Vec<f32> = idx.iter().map(|&i| lfc[i][g]).filter(|x| x.is_finite()).collect(); if v.len() < 5 { return None; }
            v.sort_by(|a, b| a.partial_cmp(b).unwrap()); let up = v.iter().filter(|&&x| x > 0.0).count();
            Some((g, v.len(), v[v.len() / 2] as f64, up as f64 / v.len() as f64, sign_p(up, v.len()))) }).collect();
        res.sort_by(|a, b| a.4.partial_cmp(&b.4).unwrap()); let m = res.len() as f64;
        let mut q = vec![0f64; res.len()]; let mut run = 1.0f64; for i in (0..res.len()).rev() { run = run.min(res[i].4 * m / (i + 1) as f64); q[i] = run; }
        for (i, r) in res.iter().enumerate() { writeln!(w, "{ct}\t{}\t{}\t{:.4}\t{:.3}\t{:.3e}\t{:.3e}", genes[r.0], r.1, r.2, r.3, r.4, q[i])?; }
    }
    Ok(())
}
