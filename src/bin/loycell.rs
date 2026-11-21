//! loycell: molecular consequences of LOY within cell types. For each male normal group with a fit, every cell gets a
//! posterior P(LOY | y, N) under f*Pois(y; a N) + (1-f)*NB(y; mu N, k), a = ambient Y rate per UMI from female cells.
//! Confident LOY (P >= 0.95) vs confident Y+ (P <= 0.01) cells are compared within depth strata (deciles of N in the
//! group) for every panel gene (chrX/chrY, PAR, XIST, CD99, TCL1A). Control contrast: the same comparison with cells split
//! by pseudo-Y set 1 (zero vs non-zero), which carries the selection-by-depth artifact but no chromosome loss.
//! usage: loycell <prefix> [min_loy_cells=10]   needs .cells.bin .groups.tsv .info.txt .panel.bin .panel_genes.txt .fit.tsv
//! writes <prefix>.cellde.tsv (per group x gene) and <prefix>.loyfrac_ambient.tsv (per group f with ambient-aware posterior mean)
use anyhow::Result;
use loyscope::ln_gamma;
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::Write;

fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().collect();
    let pre = &a[1];
    let min_loy: usize = a.get(2).map(|x| x.parse().unwrap()).unwrap_or(10);
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
    let npg = pg.len();
    // fit parameters per gid
    let mut par: HashMap<usize, (f64, f64, f64)> = HashMap::new();
    for l in std::fs::read_to_string(format!("{pre}.fitamb.tsv"))?
        .lines()
        .skip(1)
    {
        let f: Vec<&str> = l.split('\t').collect();
        if f[16] == "NaN" {
            continue;
        }
        par.insert(f[0].parse()?, (f[12].parse()?, f[15].parse()?, f[16].parse()?));
    }
    // ambient rate from female normal cells
    let (mut fy, mut fnn) = (0f64, 0f64);
    for c in 0..nc {
        let g = &gl[v[c * rl] as usize];
        if g[3] == "female" {
            fy += v[c * rl + 2] as f64;
            fnn += v[c * rl + 1] as f64;
        }
    }
    let amb = if fnn > 0.0 { fy / fnn } else { 0.0 };
    // per cell panel counts, sparse rows
    let mut rows: Vec<Vec<(u16, u16)>> = vec![Vec::new(); nc];
    let p = std::fs::read(format!("{pre}.panel.bin"))?;
    for e in p.chunks_exact(8) {
        let cell = u32::from_le_bytes([e[0], e[1], e[2], e[3]]) as usize;
        rows[cell].push((u16::from_le_bytes([e[4], e[5]]), u16::from_le_bytes([e[6], e[7]])));
    }
    drop(p);
    let mut cells_of: Vec<Vec<usize>> = vec![Vec::new(); gl.len()];
    for c in 0..nc {
        if v[c * rl + 1] > 0 {
            cells_of[v[c * rl] as usize].push(c);
        }
    }
    let todo: Vec<usize> = (0..gl.len())
        .filter(|&g| gl[g][3] == "male" && par.contains_key(&g) && cells_of[g].len() >= 50)
        .collect();
    let res: Vec<(Vec<String>, String)> = todo
        .par_iter()
        .map(|&g| {
            let (f, mu, k) = par[&g];
            let cs = &cells_of[g];
            let n: Vec<f64> = cs.iter().map(|&c| v[c * rl + 1] as f64).collect();
            let y: Vec<f64> = cs.iter().map(|&c| v[c * rl + 2] as f64).collect();
            let pz: Vec<bool> = cs.iter().map(|&c| v[c * rl + 3] == 0).collect();
            let lgk = ln_gamma(k);
            let post: Vec<f64> = (0..cs.len())
                .map(|i| {
                    let lam = amb * n[i];
                    let l1 = if lam > 0.0 {
                        y[i] * lam.ln() - lam - ln_gamma(y[i] + 1.0)
                    } else if y[i] == 0.0 {
                        0.0
                    } else {
                        f64::NEG_INFINITY
                    };
                    let mn = mu * n[i];
                    let l0 = ln_gamma(y[i] + k) - lgk - ln_gamma(y[i] + 1.0)
                        + k * (k / (k + mn)).ln()
                        + y[i] * (mn / (k + mn)).ln();
                    let (a1, a0) = (f.max(1e-12).ln() + l1, (1.0 - f).max(1e-12).ln() + l0);
                    let m = a1.max(a0);
                    ((a1 - m).exp()) / ((a1 - m).exp() + (a0 - m).exp())
                })
                .collect();
            let fpost = post.iter().sum::<f64>() / post.len() as f64;
            let mean_umi = n.iter().sum::<f64>() / n.len() as f64;
            let line = format!(
                "{g}\t{}\t{}\t{}\t{}\t{}\t{mean_umi:.0}\t{f:.4}\t{fpost:.4}\t{amb:.3e}",
                gl[g][1],
                gl[g][4],
                gl[g][5],
                gl[g][6],
                cs.len()
            );
            let nl = post.iter().filter(|&&q| q >= 0.95).count();
            if nl < min_loy {
                return (Vec::new(), line);
            }
            // depth deciles
            let mut ord: Vec<usize> = (0..cs.len()).collect();
            ord.sort_by(|&i, &j| n[i].partial_cmp(&n[j]).unwrap());
            let mut dec = vec![0usize; cs.len()];
            for (r, &i) in ord.iter().enumerate() {
                dec[i] = r * 10 / cs.len();
            }
            // per stratum sums of normalized expression per gene, for 4 sets: LOY, Y+, pseudo-zero, pseudo-nonzero
            let mut s = vec![[[0f64; 4]; 10]; npg];
            let mut cnt = [[0f64; 4]; 10];
            for i in 0..cs.len() {
                let sets = [post[i] >= 0.95, post[i] <= 0.01, pz[i], !pz[i]];
                for (q, &on) in sets.iter().enumerate() {
                    if on {
                        cnt[dec[i]][q] += 1.0;
                    }
                }
                for &(gi, ct) in &rows[cs[i]] {
                    let x = ct as f64 / n[i] * 1e4;
                    for (q, &on) in sets.iter().enumerate() {
                        if on {
                            s[gi as usize][dec[i]][q] += x;
                        }
                    }
                }
            }
            // depth-standardized means: weight strata by number of LOY cells (for pseudo contrast by pseudo-zero cells)
            let mut out = Vec::new();
            for gi in 0..npg {
                let (mut ml, mut mp, mut w1, mut mz, mut mnz, mut w2) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
                for d in 0..10 {
                    if cnt[d][0] > 0.0 && cnt[d][1] > 0.0 {
                        let w = cnt[d][0];
                        ml += w * s[gi][d][0] / cnt[d][0];
                        mp += w * s[gi][d][1] / cnt[d][1];
                        w1 += w;
                    }
                    if cnt[d][2] > 0.0 && cnt[d][3] > 0.0 {
                        let w = cnt[d][2];
                        mz += w * s[gi][d][2] / cnt[d][2];
                        mnz += w * s[gi][d][3] / cnt[d][3];
                        w2 += w;
                    }
                }
                if w1 == 0.0 || mp <= 0.0 {
                    continue;
                }
                let lfc = ((ml / w1 + 1e-3) / (mp / w1 + 1e-3)).log2();
                let lfz = if w2 > 0.0 {
                    ((mz / w2 + 1e-3) / (mnz / w2 + 1e-3)).log2()
                } else {
                    f64::NAN
                };
                out.push(format!(
                    "{g}\t{}\t{}\t{}\t{}\t{nl}\t{:.4}\t{lfc:.4}\t{lfz:.4}",
                    gl[g][1],
                    gl[g][5],
                    gl[g][6],
                    pg[gi],
                    mp / w1,
                ));
            }
            (out, line)
        })
        .collect();
    let mut w = std::fs::File::create(format!("{pre}.cellde.tsv"))?;
    writeln!(
        w,
        "gid\tcell_type\tdonor\tstage\tgene\tn_loy\tmean_ypos_cp10k\tlog2fc_loy\tlog2fc_pseudo"
    )?;
    let mut w2 = std::fs::File::create(format!("{pre}.loyfrac_ambient.tsv"))?;
    writeln!(
        w2,
        "gid\tcell_type\tdisease\tdonor\tstage\tn\tmean_umi\tf_zinb\tf_posterior_ambient\tambient_rate"
    )?;
    for (o, l) in res {
        for x in o {
            writeln!(w, "{x}")?;
        }
        writeln!(w2, "{l}")?;
    }
    eprintln!("{pre}: ambient {amb:.3e}, groups {}", todo.len());
    Ok(())
}
