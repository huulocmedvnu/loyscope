//! loydepth: genome-wide expression of cells without Y reads versus cells with Y reads, matched on library size.
//! Within each selected group (male healthy tissue, model estimate >= purity x zero-Y fraction), cells are assigned to
//! deciles of total UMIs N; per decile, normalised expression (counts per 10,000 UMIs) is computed for zero-Y and Y-positive
//! cells, and the two are averaged over deciles with weights equal to the number of zero-Y cells in each decile, as in the
//! pseudoautosomal analysis. The same contrast with cells split by zero versus non-zero counts of the first autosomal control
//! set (which cannot be lost) is computed in the same groups as a control for residual depth effects.
//! usage: loydepth <h5ad> <prefix> <fit.tsv> <out> [purity=0.8] [min_y0=10] [min_ypos=20] [min_cp10k=1]
//!   prefix: loyx outputs (.cells.bin, .info.txt, .groups.tsv); out: <out>.de_depth.tsv and <out>.de_depth_control.tsv
use anyhow::{Context, Result};
use hdf5::types::{VarLenAscii, VarLenUnicode};
use loyscope::{prefetch_csr, row_ranges};
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::Write;
use std::sync::Arc;

const ND: usize = 10;

fn read_strings(ds: &hdf5::Dataset) -> Result<Vec<String>> {
    if let Ok(v) = ds.read_raw::<VarLenUnicode>() {
        return Ok(v.iter().map(|s| s.as_str().to_string()).collect());
    }
    Ok(ds
        .read_raw::<VarLenAscii>()?
        .iter()
        .map(|s| s.as_str().to_string())
        .collect())
}
fn labels(f: &hdf5::File, path: &str) -> Result<Vec<String>> {
    if let Ok(g) = f.group(path) {
        let cats = read_strings(&g.dataset("categories")?)?;
        let codes: Vec<i64> = g.dataset("codes")?.read_raw()?;
        return Ok(codes.iter().map(|&c| cats[c as usize].clone()).collect());
    }
    read_strings(&f.dataset(path)?)
}
fn sign_p(k: usize, n: usize) -> f64 {
    // two-sided exact binomial(0.5)
    let ln_choose =
        |n: usize, k: usize| -> f64 { (1..=k).map(|i| ((n - k + i) as f64 / i as f64).ln()).sum() };
    let lo = k.min(n - k);
    let mut s = 0.0;
    for i in 0..=lo {
        s += (ln_choose(n, i) - n as f64 * 2f64.ln()).exp();
    }
    (2.0 * s).min(1.0)
}

fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().collect();
    let (h5, pre, fitp, out) = (&a[1], &a[2], &a[3], &a[4]);
    let arg = |i: usize, d: f64| -> f64 { a.get(i).and_then(|s| s.parse().ok()).unwrap_or(d) };
    let (purity, min_y0, min_ypos, min_cp) = (arg(5, 0.8), arg(6, 10.0), arg(7, 20.0), arg(8, 1.0));
    let info = std::fs::read_to_string(format!("{pre}.info.txt"))?;
    let rl: usize = info
        .lines()
        .find(|l| l.starts_with("record_u32"))
        .context("record_u32")?
        .split(' ')
        .nth(1)
        .unwrap()
        .parse()?;
    let raw = std::fs::read(format!("{pre}.cells.bin"))?;
    let v: Vec<u32> = raw
        .chunks_exact(4)
        .map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]]))
        .collect();
    drop(raw);
    let nc = v.len() / rl;
    let gl: Vec<Vec<String>> = std::fs::read_to_string(format!("{pre}.groups.tsv"))?
        .lines()
        .skip(1)
        .map(|l| l.split('\t').map(String::from).collect())
        .collect();
    let mut fitf: HashMap<usize, (f64, f64)> = HashMap::new();
    for l in std::fs::read_to_string(fitp)?.lines().skip(1) {
        let f: Vec<&str> = l.split('\t').collect();
        if f[16] != "NaN" {
            fitf.insert(f[0].parse()?, (f[12].parse()?, f[11].parse()?));
        }
    }
    // cells per group
    let mut members: Vec<Vec<usize>> = vec![Vec::new(); gl.len()];
    for c in 0..nc {
        let g = v[c * rl] as usize;
        if g < gl.len() && v[c * rl + 1] > 0 {
            members[g].push(c);
        }
    }
    // selected groups, cell types with at least five
    let ok = |g: usize| -> bool {
        let x = &gl[g];
        if x[3] != "male" || x[4] != "normal" {
            return false;
        }
        let Some(&(f, nz)) = fitf.get(&g) else {
            return false;
        };
        let y0 = members[g].iter().filter(|&&c| v[c * rl + 2] == 0).count() as f64;
        let yp = members[g].len() as f64 - y0;
        y0 >= min_y0 && yp >= min_ypos && nz > 0.0 && f / nz >= purity
    };
    let mut by_ct: HashMap<String, Vec<usize>> = HashMap::new();
    for g in 0..gl.len() {
        if ok(g) {
            by_ct.entry(gl[g][1].clone()).or_default().push(g);
        }
    }
    by_ct.retain(|_, gs| gs.len() >= 5);
    let sel: Vec<usize> = {
        let mut s: Vec<usize> = by_ct.values().flatten().cloned().collect();
        s.sort();
        s
    };
    let slot_of: HashMap<usize, usize> = sel.iter().enumerate().map(|(i, &g)| (g, i)).collect();
    // decile of N within group; cell -> (slot, decile)
    let mut cell_slot = vec![u32::MAX; nc];
    let mut cell_dec = vec![0u8; nc];
    for &g in &sel {
        let mut ns: Vec<u32> = members[g].iter().map(|&c| v[c * rl + 1]).collect();
        ns.sort();
        for &c in &members[g] {
            let r = ns.partition_point(|&x| x < v[c * rl + 1]);
            cell_dec[c] = ((r * ND) / ns.len()).min(ND - 1) as u8;
            cell_slot[c] = slot_of[&g] as u32;
        }
    }
    let f = hdf5::File::open(h5).with_context(|| h5.clone())?;
    let mpath = if f.group("raw/X").is_ok() { "raw/X" } else { "X" };
    let genes = labels(&f, "raw/var/feature_name").or_else(|_| labels(&f, "var/feature_name"))?;
    let ng = genes.len();
    let indptr: Vec<i64> = f.group(mpath)?.dataset("indptr")?.read_raw()?;
    if indptr.len() - 1 != nc {
        anyhow::bail!(
            "cell count differs between {h5} ({}) and {pre}.cells.bin ({nc})",
            indptr.len() - 1
        );
    }
    eprintln!(
        "{pre}: {} groups in {} cell types, {ng} genes",
        sel.len(),
        by_ct.len()
    );
    // accumulators: [slot][split: 0 Y, 1 control][side: 0 zero, 1 positive][decile] x genes; plus cell and UMI counts
    let nslot = sel.len() * 2 * 2 * ND;
    let mut acc = vec![0f32; nslot * ng];
    let mut cnt = vec![[0f64; 2]; nslot]; // cells, UMIs
    let idx =
        |s: usize, split: usize, side: usize, d: usize| -> usize { ((s * 2 + split) * 2 + side) * ND + d };
    for item in prefetch_csr(h5, mpath, Arc::new(indptr.clone()), row_ranges(nc, 50_000, 1)) {
        let (r0, p, ix, val) = item?;
        for c in 0..p.len() - 1 {
            let cell = r0 + c;
            let s = cell_slot[cell];
            if s == u32::MAX {
                continue;
            }
            let (s, d) = (s as usize, cell_dec[cell] as usize);
            let n = v[cell * rl + 1] as f64;
            for split in 0..2 {
                let side = if v[cell * rl + 2 + split] == 0 { 0 } else { 1 };
                let k = idx(s, split, side, d);
                cnt[k][0] += 1.0;
                cnt[k][1] += n;
                let base = k * ng;
                for q in p[c]..p[c + 1] {
                    acc[base + ix[q] as usize] += val[q];
                }
            }
        }
    }
    // per group and gene: decile-matched log2 ratio
    let lfc = |split: usize| -> Vec<Vec<f32>> {
        (0..sel.len())
            .into_par_iter()
            .map(|s| {
                (0..ng)
                    .map(|g| {
                        let (mut w, mut m0, mut m1) = (0.0, 0.0, 0.0);
                        for d in 0..ND {
                            let (k0, k1) = (idx(s, split, 0, d), idx(s, split, 1, d));
                            if cnt[k0][0] == 0.0 || cnt[k1][0] == 0.0 {
                                continue;
                            }
                            let wt = cnt[k0][0];
                            m0 += wt * acc[k0 * ng + g] as f64 / cnt[k0][1] * 1e4;
                            m1 += wt * acc[k1 * ng + g] as f64 / cnt[k1][1] * 1e4;
                            w += wt;
                        }
                        if w == 0.0 || m1 / w < min_cp {
                            f32::NAN
                        } else {
                            ((m0 / w + 0.05) / (m1 / w + 0.05)).log2() as f32
                        }
                    })
                    .collect()
            })
            .collect()
    };
    for (split, suffix) in [(0usize, "de_depth"), (1usize, "de_depth_control")] {
        let l = lfc(split);
        let mut w = std::fs::File::create(format!("{out}.{suffix}.tsv"))?;
        writeln!(
            w,
            "cell_type\tgene\tn_groups\tmedian_log2fc\tshare_up\tsign_p\tq_bh"
        )?;
        let mut cts: Vec<&String> = by_ct.keys().collect();
        cts.sort();
        for ct in cts {
            let ss: Vec<usize> = by_ct[ct].iter().map(|g| slot_of[g]).collect();
            let mut res: Vec<(usize, usize, f64, f64, f64)> = (0..ng)
                .into_par_iter()
                .filter_map(|g| {
                    let mut x: Vec<f32> = ss.iter().map(|&s| l[s][g]).filter(|z| z.is_finite()).collect();
                    if x.len() < 5 {
                        return None;
                    }
                    x.sort_by(|a, b| a.partial_cmp(b).unwrap());
                    let up = x.iter().filter(|&&z| z > 0.0).count();
                    Some((
                        g,
                        x.len(),
                        x[x.len() / 2] as f64,
                        up as f64 / x.len() as f64,
                        sign_p(up, x.len()),
                    ))
                })
                .collect();
            res.sort_by(|a, b| a.4.partial_cmp(&b.4).unwrap());
            let m = res.len() as f64;
            let mut q = vec![0f64; res.len()];
            let mut run = 1.0f64;
            for i in (0..res.len()).rev() {
                run = run.min(res[i].4 * m / (i + 1) as f64);
                q[i] = run;
            }
            for (i, r) in res.iter().enumerate() {
                writeln!(
                    w,
                    "{ct}\t{}\t{}\t{:.4}\t{:.3}\t{:.3e}\t{:.3e}",
                    genes[r.0], r.1, r.2, r.3, r.4, q[i]
                )?;
            }
        }
    }
    Ok(())
}
