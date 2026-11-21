//! loybench: Random Forest LOY caller (as in Nature 2025: 9 YchrS genes, male = WTY vs female = LOY, 70/30 split,
//! 100 trees, sqrt features, bootstrap) vs naive zero-Y calls vs the ZINB estimator, on (i) observed male groups,
//! (ii) spike-in LOY with known fraction (Y counts of chosen cells replaced by ambient-level Poisson counts estimated
//! from female cells), (iii) binomial thinning of depth to one half.
//! RF is histogram-based (64 quantile bins per feature), trees built in parallel.
//! usage: loybench <prefix> [n_train_max=300000] [min_cells=50]
//! needs <prefix>.cells.bin/.groups.tsv/.info.txt/.panel.bin/.panel_genes.txt (loyx v2)
//! writes <prefix>.rf_groups.tsv, <prefix>.bench.tsv, <prefix>.rf_info.txt
use anyhow::{bail, Result};
use loyscope::fit;
use rand::seq::SliceRandom;
use rand::{Rng, SeedableRng};
use rand_distr::{Binomial, Distribution, Poisson};
use rand_xoshiro::Xoshiro256PlusPlus;
use rayon::prelude::*;
use std::io::Write;

const YCHRS: &[&str] = &[
    "DDX3Y", "UTY", "KDM5D", "USP9Y", "ZFY", "RPS4Y1", "TMSB4Y", "EIF1AY", "NLGN4Y",
];
const NB: usize = 64;
const NTREE: usize = 100;

struct Bins {
    cuts: Vec<Vec<f32>>,
}
impl Bins {
    fn fit(x: &[Vec<f32>], idx: &[usize], nf: usize) -> Self {
        let cuts = (0..nf)
            .map(|j| {
                let mut v: Vec<f32> = idx.iter().map(|&i| x[i][j]).filter(|&a| a > 0.0).collect();
                v.sort_by(|a, b| a.partial_cmp(b).unwrap());
                let mut c: Vec<f32> = Vec::new();
                if !v.is_empty() {
                    for q in 1..NB - 1 {
                        let t = v[(q * v.len() / (NB - 1)).min(v.len() - 1)];
                        if c.last().is_none_or(|&l| t > l) {
                            c.push(t);
                        }
                    }
                }
                c
            })
            .collect();
        Bins { cuts }
    }
    #[inline]
    fn bin(&self, j: usize, a: f32) -> u8 {
        if a <= 0.0 {
            0
        } else {
            (1 + self.cuts[j].partition_point(|&c| c < a)).min(NB - 1) as u8
        }
    }
}
#[derive(Clone)]
struct Node {
    feat: u8,
    thr: u8,
    left: u32,
    right: u32,
    prob: f32,
}
fn build_tree(xb: &[Vec<u8>], lab: &[u8], nf: usize, mtry: usize, seed: u64) -> Vec<Node> {
    let n = xb.len();
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    let boot: Vec<u32> = (0..n).map(|_| rng.gen_range(0..n) as u32).collect();
    let mut nodes: Vec<Node> = Vec::new();
    // stack of (node index, sample indices, depth)
    nodes.push(Node {
        feat: 255,
        thr: 0,
        left: 0,
        right: 0,
        prob: 0.0,
    });
    let mut stack: Vec<(usize, Vec<u32>, usize)> = vec![(0, boot, 0)];
    let feats: Vec<usize> = (0..nf).collect();
    while let Some((ni, s, depth)) = stack.pop() {
        let pos = s.iter().filter(|&&i| lab[i as usize] == 1).count();
        let tot = s.len();
        nodes[ni].prob = pos as f32 / tot as f32;
        if pos == 0 || pos == tot || tot < 2 || depth >= 40 {
            continue;
        }
        let parent = 1.0 - ((pos as f64 / tot as f64).powi(2) + ((tot - pos) as f64 / tot as f64).powi(2));
        let mut best = (0.0f64, 255u8, 0u8);
        let mut fs = feats.clone();
        fs.shuffle(&mut rng);
        for &j in fs.iter().take(mtry) {
            let mut h = [[0u32; 2]; NB];
            for &i in &s {
                h[xb[i as usize][j] as usize][lab[i as usize] as usize] += 1;
            }
            let (mut l0, mut l1) = (0u32, 0u32);
            let (t0, t1) = ((tot - pos) as u32, pos as u32);
            for b in 0..NB - 1 {
                l0 += h[b][0];
                l1 += h[b][1];
                let nl = (l0 + l1) as f64;
                let nr = (tot as f64) - nl;
                if nl < 1.0 || nr < 1.0 {
                    continue;
                }
                let gl = 1.0 - ((l0 as f64 / nl).powi(2) + (l1 as f64 / nl).powi(2));
                let (r0, r1) = ((t0 - l0) as f64, (t1 - l1) as f64);
                let gr = 1.0 - ((r0 / nr).powi(2) + (r1 / nr).powi(2));
                let gain = parent - (nl * gl + nr * gr) / tot as f64;
                if gain > best.0 + 1e-12 {
                    best = (gain, j as u8, b as u8);
                }
            }
        }
        if best.1 == 255 {
            continue;
        }
        let (j, t) = (best.1 as usize, best.2);
        let (ls, rs): (Vec<u32>, Vec<u32>) = s.iter().partition(|&&i| xb[i as usize][j] <= t);
        let (li, ri) = (nodes.len(), nodes.len() + 1);
        nodes.push(Node {
            feat: 255,
            thr: 0,
            left: 0,
            right: 0,
            prob: 0.0,
        });
        nodes.push(Node {
            feat: 255,
            thr: 0,
            left: 0,
            right: 0,
            prob: 0.0,
        });
        nodes[ni].feat = j as u8;
        nodes[ni].thr = t;
        nodes[ni].left = li as u32;
        nodes[ni].right = ri as u32;
        stack.push((li, ls, depth + 1));
        stack.push((ri, rs, depth + 1));
    }
    nodes
}
fn predict(forest: &[Vec<Node>], xb: &[u8]) -> f32 {
    let mut s = 0.0;
    for t in forest {
        let mut i = 0;
        loop {
            let nd = &t[i];
            if nd.feat == 255 {
                s += nd.prob;
                break;
            }
            i = if xb[nd.feat as usize] <= nd.thr {
                nd.left
            } else {
                nd.right
            } as usize;
        }
    }
    s / forest.len() as f32
}

fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().collect();
    let pre = &a[1];
    let ntrain_max: usize = a.get(2).map(|x| x.parse().unwrap()).unwrap_or(300_000);
    let minc: usize = a.get(3).map(|x| x.parse().unwrap()).unwrap_or(50);
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
    drop(raw);
    let nc = v.len() / rl;
    let gl: Vec<Vec<String>> = std::fs::read_to_string(format!("{pre}.groups.tsv"))?
        .lines()
        .skip(1)
        .map(|l| l.split('\t').map(String::from).collect())
        .collect();
    let pg: Vec<String> = std::fs::read_to_string(format!("{pre}.panel_genes.txt"))?
        .lines()
        .map(String::from)
        .collect();
    let fmap: Vec<i32> = pg
        .iter()
        .map(|g| YCHRS.iter().position(|y| y == g).map_or(-1, |p| p as i32))
        .collect();
    let present: Vec<usize> = (0..YCHRS.len()).filter(|&j| fmap.contains(&(j as i32))).collect();
    let nf = YCHRS.len();
    if present.is_empty() {
        bail!("no YchrS genes in panel");
    }
    // per cell counts of the 9 genes
    let mut yc = vec![0u16; nc * nf];
    let p = std::fs::read(format!("{pre}.panel.bin"))?;
    for e in p.chunks_exact(8) {
        let g = u16::from_le_bytes([e[4], e[5]]) as usize;
        let j = fmap[g];
        if j >= 0 {
            let cell = u32::from_le_bytes([e[0], e[1], e[2], e[3]]) as usize;
            yc[cell * nf + j as usize] = u16::from_le_bytes([e[6], e[7]]);
        }
    }
    drop(p);
    let ncell = |c: usize| v[c * rl + 1] as f32;
    let gid = |c: usize| v[c * rl] as usize;
    let is_normal = |c: usize| gl[gid(c)][4] == "normal";
    let sex = |c: usize| gl[gid(c)][3].as_str();
    // ambient rate per gene from female normal cells
    let (mut fy, mut fnn) = (vec![0f64; nf], 0f64);
    for c in 0..nc {
        if is_normal(c) && sex(c) == "female" {
            fnn += ncell(c) as f64;
            for j in 0..nf {
                fy[j] += yc[c * nf + j] as f64;
            }
        }
    }
    let amb: Vec<f64> = fy.iter().map(|x| if fnn > 0.0 { x / fnn } else { 0.0 }).collect();
    let feat = |cnt: &[u16], n: f32| -> Vec<f32> {
        cnt.iter()
            .map(|&x| (1.0 + x as f32 / n.max(1.0) * 1e4).ln())
            .collect()
    };
    // RF training data
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(20261121);
    let mut pool: Vec<usize> = (0..nc)
        .filter(|&c| is_normal(c) && (sex(c) == "male" || sex(c) == "female") && ncell(c) > 0.0)
        .collect();
    pool.shuffle(&mut rng);
    pool.truncate(ntrain_max);
    let ntr = pool.len() * 7 / 10;
    let x: Vec<Vec<f32>> = pool
        .iter()
        .map(|&c| feat(&yc[c * nf..(c + 1) * nf], ncell(c)))
        .collect();
    let lab: Vec<u8> = pool.iter().map(|&c| (sex(c) == "female") as u8).collect();
    let tr: Vec<usize> = (0..ntr).collect();
    let bins = Bins::fit(&x, &tr, nf);
    let xb_all: Vec<Vec<u8>> = x
        .iter()
        .map(|r| (0..nf).map(|j| bins.bin(j, r[j])).collect())
        .collect();
    let xtr: Vec<Vec<u8>> = xb_all[..ntr].to_vec();
    let ltr: Vec<u8> = lab[..ntr].to_vec();
    let mtry = ((nf as f64).sqrt().floor() as usize).max(1);
    let forest: Vec<Vec<Node>> = (0..NTREE)
        .into_par_iter()
        .map(|t| build_tree(&xtr, &ltr, nf, mtry, 1000 + t as u64))
        .collect();
    let acc = (ntr..pool.len())
        .filter(|&i| ((predict(&forest, &xb_all[i]) > 0.5) as u8) == lab[i])
        .count() as f64
        / (pool.len() - ntr) as f64;
    let nfem = lab.iter().filter(|&&l| l == 1).count();
    std::fs::write(format!("{pre}.rf_info.txt"), format!("train {ntr} test {} female_share {:.3} test_accuracy {acc:.4} genes_present {}/{} trees {NTREE} mtry {mtry}\nambient per UMI {:?}\n", pool.len() - ntr, nfem as f64 / lab.len() as f64, present.len(), nf, amb))?;
    eprintln!("RF test accuracy {acc:.4}");
    let rf_call = |cnt: &[u16], n: f32| -> bool {
        let f = feat(cnt, n);
        let b: Vec<u8> = (0..nf).map(|j| bins.bin(j, f[j])).collect();
        predict(&forest, &b) > 0.5
    };
    // groups: male normal with >= minc cells
    let mut cells_of: Vec<Vec<usize>> = vec![Vec::new(); gl.len()];
    for c in 0..nc {
        if ncell(c) > 0.0 {
            cells_of[gid(c)].push(c);
        }
    }
    let todo: Vec<usize> = (0..gl.len())
        .filter(|&g| gl[g][3] == "male" && gl[g][4] == "normal" && cells_of[g].len() >= minc)
        .collect();
    eprintln!("male normal groups {}", todo.len());
    let estimate = |ns: &[f32], ys: &[Vec<u16>]| -> (f64, f64, f64, f64) {
        // naive, rf, zinb f, expected Y per intact cell
        let tot: Vec<f64> = ys.iter().map(|r| r.iter().map(|&x| x as f64).sum()).collect();
        let naive = tot.iter().filter(|&&t| t == 0.0).count() as f64 / tot.len() as f64;
        let rf = (0..ns.len()).filter(|&i| rf_call(&ys[i], ns[i])).count() as f64 / ns.len() as f64;
        let nn: Vec<f64> = ns.iter().map(|&x| x as f64).collect();
        let (f, mu, _) = fit(&nn, &tot);
        let mn = nn.iter().sum::<f64>() / nn.len() as f64;
        (naive, rf, f, mu * mn)
    };
    let rows: Vec<String> = todo
        .par_iter()
        .map(|&g| {
            let cs = &cells_of[g];
            let mut rg = Xoshiro256PlusPlus::seed_from_u64(g as u64 * 7 + 3);
            let ns: Vec<f32> = cs.iter().map(|&c| ncell(c)).collect();
            let ys: Vec<Vec<u16>> = cs.iter().map(|&c| yc[c * nf..(c + 1) * nf].to_vec()).collect();
            let mean_umi = ns.iter().map(|&x| x as f64).sum::<f64>() / ns.len() as f64;
            let mut out = Vec::new();
            let head = format!(
                "{g}\t{}\t{}\t{}\t{}\t{mean_umi:.0}",
                gl[g][1],
                gl[g][5],
                gl[g][6],
                cs.len()
            );
            let base = estimate(&ns, &ys);
            out.push(format!(
                "{head}\t{:.3}\tobserved\t0\t{:.4}\t{:.4}\t{:.4}",
                base.3, base.0, base.1, base.2
            ));
            for &ft in &[0.02f64, 0.05, 0.10, 0.20] {
                let mut ys2 = ys.clone();
                let mut idx: Vec<usize> = (0..cs.len()).collect();
                idx.shuffle(&mut rg);
                for &i in idx.iter().take((ft * cs.len() as f64).round() as usize) {
                    for j in 0..nf {
                        let lam = amb[j] * ns[i] as f64;
                        ys2[i][j] = if lam > 0.0 {
                            Poisson::new(lam).unwrap().sample(&mut rg) as u16
                        } else {
                            0
                        };
                    }
                }
                let e = estimate(&ns, &ys2);
                out.push(format!(
                    "{head}\t{:.3}\tspike\t{ft}\t{:.4}\t{:.4}\t{:.4}",
                    base.3, e.0, e.1, e.2
                ));
            }
            // thinning to half depth
            let ns2: Vec<f32> = ns
                .iter()
                .map(|&n| Binomial::new(n as u64, 0.5).unwrap().sample(&mut rg).max(1) as f32)
                .collect();
            let ys3: Vec<Vec<u16>> = ys
                .iter()
                .map(|r| {
                    r.iter()
                        .map(|&x| Binomial::new(x as u64, 0.5).unwrap().sample(&mut rg) as u16)
                        .collect()
                })
                .collect();
            let e = estimate(&ns2, &ys3);
            out.push(format!(
                "{head}\t{:.3}\tthin0.5\t0\t{:.4}\t{:.4}\t{:.4}",
                base.3, e.0, e.1, e.2
            ));
            out.join("\n")
        })
        .collect();
    let mut w = std::fs::File::create(format!("{pre}.bench.tsv"))?;
    writeln!(
        w,
        "gid\tcell_type\tdonor\tstage\tn\tmean_umi\texpY\tscenario\tf_true\tnaive\trf\tzinb"
    )?;
    for r in rows {
        writeln!(w, "{r}")?;
    }
    eprintln!("{pre}: bench done");
    Ok(())
}
