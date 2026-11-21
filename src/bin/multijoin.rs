//! multijoin: chromatin ground truth for RNA LOY calls in multiome nuclei.
//! ATAC: per male nucleus, chrY (MSY) fragments T_y given total fragments T follow
//!   f * Pois(a_f * T) + (1 - f) * NB(m * T, k), with a_f the female background rate per fragment (fixed);
//!   fitted per RNA group (donor x cell type); posterior P_atac(LOY) per nucleus.
//! RNA: per nucleus naive call (zero Y UMIs) and model posterior from the ambient-aware fit of the same group.
//! usage: multijoin <rna_prefix> <atac.tsv.gz> <out_prefix> [min_nuclei=30]
//! needs <rna_prefix>.cells.bin .groups.tsv .info.txt .obs.tsv.gz .fitamb.tsv; writes <out>.groups_cmp.tsv, <out>.nuclei.tsv.gz,
//! <out>.metrics.txt
use anyhow::Result;
use flate2::read::MultiGzDecoder;
use flate2::{write::GzEncoder, Compression};
use loyscope::{fit_bg, ln_gamma};
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::{BufRead, BufReader, BufWriter, Write};

fn lines_gz(p: &str) -> Result<std::io::Lines<BufReader<MultiGzDecoder<std::fs::File>>>> {
    Ok(BufReader::new(MultiGzDecoder::new(std::fs::File::open(p)?)).lines())
}
/// posterior P(LOY | y, N) under f*Pois(lam_bg) + (1-f)*NB(mu N, k)
fn post(y: f64, n: f64, f: f64, mu: f64, k: f64, bg_rate: f64) -> f64 {
    let lam = (bg_rate * n).max(1e-12);
    let l1 = y * lam.ln() - lam - ln_gamma(y + 1.0);
    let mn = mu * n;
    let l0 = ln_gamma(y + k) - ln_gamma(k) - ln_gamma(y + 1.0)
        + k * (k / (k + mn)).ln()
        + y * (mn / (k + mn)).ln();
    let (a1, a0) = (f.max(1e-12).ln() + l1, (1.0 - f).max(1e-12).ln() + l0);
    let m = a1.max(a0);
    (a1 - m).exp() / ((a1 - m).exp() + (a0 - m).exp())
}

fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().collect();
    let (rp, ap, out) = (&a[1], &a[2], &a[3]);
    let minn: usize = a.get(4).map(|x| x.parse().unwrap()).unwrap_or(30);
    // ATAC table
    let mut atac: HashMap<String, (f64, f64)> = HashMap::new();
    let mut apseudo: HashMap<String, Vec<f64>> = HashMap::new(); // optional pseudo-Y peak sets (columns 5..)
    for l in lines_gz(ap)?.skip(1) {
        let l = l?;
        let f: Vec<&str> = l.split('\t').collect();
        atac.insert(f[0].to_string(), (f[1].parse()?, f[2].parse()?));
        if f.len() > 4 {
            apseudo.insert(
                f[0].to_string(),
                f[4..].iter().map(|x| x.parse().unwrap_or(0.0)).collect(),
            );
        }
    }
    // RNA cells, groups, obs ids
    let info = std::fs::read_to_string(format!("{rp}.info.txt"))?;
    let rl: usize = info
        .lines()
        .find(|l| l.starts_with("record_u32"))
        .unwrap()
        .split(' ')
        .nth(1)
        .unwrap()
        .parse()?;
    let raw = std::fs::read(format!("{rp}.cells.bin"))?;
    let v: Vec<u32> = raw
        .chunks_exact(4)
        .map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]]))
        .collect();
    drop(raw);
    let nc = v.len() / rl;
    let gl: Vec<Vec<String>> = std::fs::read_to_string(format!("{rp}.groups.tsv"))?
        .lines()
        .skip(1)
        .map(|l| l.split('\t').map(String::from).collect())
        .collect();
    let mut ids: Vec<String> = Vec::with_capacity(nc);
    for l in lines_gz(&format!("{rp}.obs.tsv.gz"))?.skip(1) {
        let l = l?;
        ids.push(l.split('\t').next().unwrap_or("").to_string());
    }
    let mut par: HashMap<usize, (f64, f64, f64)> = HashMap::new();
    for l in std::fs::read_to_string(format!("{rp}.fitamb.tsv"))?
        .lines()
        .skip(1)
    {
        let f: Vec<&str> = l.split('\t').collect();
        if f[16] != "NaN" {
            par.insert(f[0].parse()?, (f[12].parse()?, f[15].parse()?, f[16].parse()?));
        }
    }
    // background rates from female nuclei: RNA Y per UMI, ATAC chrY per fragment
    let (mut ry, mut rn, mut ay, mut an, mut matched) = (0f64, 0f64, 0f64, 0f64, 0usize);
    let mut joined: Vec<(usize, f64, f64, f64, f64)> = Vec::new(); // cell, N, y, T, Ty (male only)
    for c in 0..nc {
        let g = &gl[v[c * rl] as usize];
        let Some(&(t, ty)) = atac.get(&ids[c]) else {
            continue;
        };
        matched += 1;
        let (n, y) = (v[c * rl + 1] as f64, v[c * rl + 2] as f64);
        if g[3] == "female" {
            ry += y;
            rn += n;
            ay += ty;
            an += t;
        } else if g[3] == "male" && t > 0.0 && n > 0.0 {
            joined.push((c, n, y, t, ty));
        }
    }
    let (mut rna_bg, mut atac_bg) = (ry / rn.max(1.0), ay / an.max(1.0));
    // libraries without female nuclei: background rates supplied from female samples of the same study
    if let Ok(x) = std::env::var("ATAC_BG") {
        atac_bg = x.parse()?;
    }
    if let Ok(x) = std::env::var("RNA_BG") {
        rna_bg = x.parse()?;
    }
    eprintln!("matched {matched}/{nc} RNA nuclei to ATAC; male joined {}; female bg RNA {rna_bg:.3e}/UMI, ATAC {atac_bg:.3e}/fragment", joined.len());
    let mut by: HashMap<usize, Vec<usize>> = HashMap::new();
    for (i, x) in joined.iter().enumerate() {
        by.entry(v[x.0 * rl] as usize).or_default().push(i);
    }
    // absolute background rates: MSY fragments at the female rate per fragment; for each autosomal control peak set the
    // same contamination fraction times its pooled rate in male nuclei
    let (tot_t, tot_y): (f64, f64) = joined.iter().fold((0.0, 0.0), |a, x| (a.0 + x.3, a.1 + x.4));
    let cfrac = atac_bg / (tot_y / tot_t).max(1e-300);
    let np_all = joined
        .first()
        .and_then(|x| apseudo.get(&ids[x.0]))
        .map_or(0, |v| v.len());
    let pseudo_bg: Vec<f64> = (0..np_all)
        .map(|j| {
            cfrac
                * joined
                    .iter()
                    .map(|x| apseudo.get(&ids[x.0]).map_or(0.0, |v| v[j]))
                    .sum::<f64>()
                / tot_t
        })
        .collect();
    let rho_atac = cfrac; // reported as the contamination fraction implied by female nuclei
    let gids: Vec<usize> = by.keys().copied().filter(|g| by[g].len() >= minn).collect();
    let res: Vec<(String, Vec<String>, [u64; 8])> = gids.par_iter().map(|&g| {
        let idx = &by[&g];
        let t: Vec<f64> = idx.iter().map(|&i| joined[i].3).collect();
        let ty: Vec<f64> = idx.iter().map(|&i| joined[i].4).collect();
        let (fa, ma, ka) = fit_bg(&t, &ty, atac_bg);
        // ATAC pseudo-Y control: same model on autosomal peak sets matched to chrY peaks (true f = 0)
        let np = idx.first().and_then(|&i| apseudo.get(&ids[joined[i].0])).map_or(0, |v| v.len());
        let mut pf: Vec<f64> = (0..np).map(|j| { let pj: Vec<f64> = idx.iter().map(|&i| apseudo[&ids[joined[i].0]][j]).collect(); fit_bg(&t, &pj, pseudo_bg[j]).0 }).collect();
        pf.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let fa_pseudo = if pf.is_empty() { f64::NAN } else { pf[pf.len() / 2] };
        let (fr, mr, kr) = par.get(&g).copied().unwrap_or((f64::NAN, f64::NAN, f64::NAN));
        let mut cm = [0u64; 8]; // naive: TP FP FN TN ; model: TP FP FN TN (truth = P_atac >= 0.5; ambiguous 0.05-0.95 skipped)
        let mut nuc = Vec::new();
        let (mut naive_n, mut pr_sum, mut pa_sum) = (0f64, 0f64, 0f64);
        for &i in idx {
            let (_, n, y, tt, tyy) = joined[i];
            let pa = if ka.is_nan() { f64::NAN } else { post(tyy, tt, fa, ma, ka, atac_bg) };
            let prn = if kr.is_nan() { f64::NAN } else { post(y, n, fr, mr, kr, rna_bg) };
            let naive = y == 0.0;
            naive_n += naive as u8 as f64; pr_sum += prn; pa_sum += pa;
            if pa.is_finite() && (pa >= 0.95 || pa <= 0.05) {
                let truth = pa >= 0.95;
                let o = match (naive, truth) { (true, true) => 0, (true, false) => 1, (false, true) => 2, (false, false) => 3 }; cm[o] += 1;
                if prn.is_finite() { let call = prn >= 0.5; let o2 = 4 + match (call, truth) { (true, true) => 0, (true, false) => 1, (false, true) => 2, (false, false) => 3 }; cm[o2] += 1; }
            }
            nuc.push(format!("{}\t{g}\t{n}\t{y}\t{tt}\t{tyy}\t{pa:.4}\t{prn:.4}", ids[joined[i].0]));
        }
        let m = idx.len() as f64; let mean_umi = idx.iter().map(|&i| joined[i].1).sum::<f64>() / m;
        let line = format!("{g}\t{}\t{}\t{}\t{}\t{}\t{mean_umi:.0}\t{fa:.4}\t{fa_pseudo:.4}\t{:.4}\t{fr:.4}\t{:.4}\t{:.4}\t{:.2}", gl[g][1], gl[g][4], gl[g][5], gl[g][6], idx.len(),
            pa_sum / m, pr_sum / m, naive_n / m, mr * mean_umi);
        (line, nuc, cm)
    }).collect();
    let mut wg = BufWriter::new(std::fs::File::create(format!("{out}.groups_cmp.tsv"))?);
    writeln!(wg, "gid\tcell_type\tdisease\tdonor\tstage\tn\tmean_umi\tf_atac\tf_atac_pseudo\tpost_atac_mean\tf_rna_model\tpost_rna_mean\tnaive_zero\texpY_rna")?;
    let mut wn = GzEncoder::new(
        BufWriter::new(std::fs::File::create(format!("{out}.nuclei.tsv.gz"))?),
        Compression::fast(),
    );
    writeln!(
        wn,
        "cell_id\tgid\trna_umi\trna_y\tatac_total\tatac_chrY\tp_atac\tp_rna_model"
    )?;
    let mut tot = [0u64; 8];
    for (l, nuc, cm) in &res {
        writeln!(wg, "{l}")?;
        for x in nuc {
            writeln!(wn, "{x}")?;
        }
        for j in 0..8 {
            tot[j] += cm[j];
        }
    }
    wn.finish()?;
    let met = |c: &[u64]| {
        let (tp, fp, fnn, tn) = (c[0] as f64, c[1] as f64, c[2] as f64, c[3] as f64);
        format!(
            "TP {tp} FP {fp} FN {fnn} TN {tn} | sensitivity {:.3} specificity {:.4} PPV {:.3}",
            tp / (tp + fnn),
            tn / (tn + fp),
            tp / (tp + fp)
        )
    };
    let txt = format!("matched {matched} of {nc} RNA nuclei; male nuclei joined {}; groups >= {minn}: {}\nATAC female background {atac_bg:.3e}/fragment, contamination fraction {rho_atac:.4}; RNA female background {rna_bg:.3e}/UMI\nnaive zero-Y call vs ATAC truth: {}\nmodel posterior >= 0.5 vs ATAC truth: {}\n",
        joined.len(), res.len(), met(&tot[0..4]), met(&tot[4..8]));
    std::fs::write(format!("{out}.metrics.txt"), &txt)?;
    eprint!("{txt}");
    Ok(())
}
