//! tenx: 10x Genomics filtered_feature_bc_matrix.h5 files (Gene Expression features only) -> loyx-compatible outputs,
//! with sex inferred per sample (Y-set and XIST CP10k), coarse kidney cell types from marker sets, and the published
//! density-trough LOY call (Gaussian KDE of log1p Y CP10k pooled over male cells; threshold = density minimum between 0
//! and the main mode) for direct comparison.
//! usage: tenx <manifest.tsv: sample  path.h5> <tss.tsv> <out_prefix>
use anyhow::{Context, Result};
use hdf5::types::{VarLenAscii, VarLenUnicode, FixedAscii};
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;
use std::collections::HashMap;
use std::io::{BufWriter, Write};

const YSET: &[&str] = &["RPS4Y1", "ZFY", "USP9Y", "DDX3Y", "UTY", "KDM5D", "EIF1AY", "TXLNGY", "NLGN4Y", "PRKY", "TMSB4Y"];
const NP: usize = 20;
const MARKERS: &[(&str, &[&str])] = &[("proximal tubule", &["LRP2", "CUBN", "SLC34A1"]), ("thick ascending limb", &["UMOD", "SLC12A1"]),
    ("distal convoluted tubule", &["SLC12A3"]), ("connecting tubule", &["CALB1", "SLC8A1"]), ("principal cell", &["AQP2", "SCNN1G"]),
    ("intercalated cell", &["SLC4A1", "SLC26A4", "ATP6V0D2"]), ("podocyte", &["NPHS1", "NPHS2", "PODXL"]), ("endothelial cell", &["PECAM1", "EMCN"]),
    ("fibroblast", &["COL1A2", "DCN", "PDGFRB"]), ("immune cell", &["PTPRC", "CD74"])];

fn strings(ds: &hdf5::Dataset) -> Result<Vec<String>> {
    if let Ok(v) = ds.read_raw::<VarLenUnicode>() { return Ok(v.iter().map(|s| s.as_str().to_string()).collect()); }
    if let Ok(v) = ds.read_raw::<VarLenAscii>() { return Ok(v.iter().map(|s| s.as_str().to_string()).collect()); }
    for n in [256usize] { let _ = n; }
    let v = ds.read_raw::<FixedAscii<256>>().context("strings")?; Ok(v.iter().map(|s| s.as_str().to_string()).collect())
}
struct Sample { name: String, genes_idx: Vec<i64>, cols: Vec<(Vec<u32>, Vec<u32>)> } // per cell: (gene ids in global space, counts)

fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().collect();
    let (man, tss, out) = (&a[1], &a[2], &a[3]);
    let mut gidx: HashMap<String, u32> = HashMap::new(); let mut gnames: Vec<String> = Vec::new();
    let mut samples: Vec<Sample> = Vec::new();
    for l in std::fs::read_to_string(man)?.lines() {
        let x: Vec<&str> = l.split('\t').collect(); let f = hdf5::File::open(x[1]).with_context(|| x[1].to_string())?;
        let m = f.group("matrix")?;
        let names = strings(&m.dataset("features/name")?)?; let ftype = strings(&m.dataset("features/feature_type")?)?;
        let map: Vec<i64> = names.iter().zip(&ftype).map(|(n, t)| if t == "Gene Expression" { let k = gidx.len() as u32; *gidx.entry(n.clone()).or_insert_with(|| { gnames.push(n.clone()); k }) as i64 } else { -1 }).collect();
        let data: Vec<i64> = m.dataset("data")?.read_raw::<i64>()?; let ind: Vec<i64> = m.dataset("indices")?.read_raw::<i64>()?; let ip: Vec<i64> = m.dataset("indptr")?.read_raw::<i64>()?;
        let cols: Vec<(Vec<u32>, Vec<u32>)> = (0..ip.len() - 1).map(|c| { let (lo, hi) = (ip[c] as usize, ip[c + 1] as usize); let mut g = Vec::new(); let mut v = Vec::new();
            for q in lo..hi { let gg = map[ind[q] as usize]; if gg >= 0 { g.push(gg as u32); v.push(data[q] as u32); } } (g, v) }).collect();
        eprintln!("{}: {} cells", x[0], cols.len());
        samples.push(Sample { name: x[0].to_string(), genes_idx: map, cols });
    }
    let ng = gnames.len();
    let find = |s: &str| gidx.get(s).map(|&i| i as usize);
    let yidx: Vec<usize> = YSET.iter().filter_map(|s| find(s)).collect();
    let xist = find("XIST");
    // sex per sample
    let mut male = Vec::new(); let mut sexinfo = Vec::new();
    for s in &samples {
        let (mut ny, mut nx, mut nt) = (0f64, 0f64, 0f64);
        for (g, v) in &s.cols { for (gg, vv) in g.iter().zip(v) { let gg = *gg as usize; nt += *vv as f64; if yidx.contains(&gg) { ny += *vv as f64; } if Some(gg) == xist { nx += *vv as f64; } } }
        let (ycp, xcp) = (ny / nt * 1e4, nx / nt * 1e4); let m = ycp > 0.5 && ycp > xcp; male.push(m);
        sexinfo.push(format!("{}\tY_cp10k {ycp:.3}\tXIST_cp10k {xcp:.3}\t{}", s.name, if m { "male" } else { "female" }));
    }
    // chromosome, pseudo sets (expression bins over male cells)
    let mut chr: HashMap<String, String> = HashMap::new();
    for l in std::fs::read_to_string(tss)?.lines() { let x: Vec<&str> = l.split('\t').collect(); if x.len() >= 2 { chr.insert(x[0].into(), x[1].into()); } }
    let mut gsum = vec![0f64; ng];
    for (si, s) in samples.iter().enumerate() { if !male[si] { continue; } for (g, v) in &s.cols { let t: f64 = v.iter().map(|&x| x as f64).sum(); if t > 0.0 { for (gg, vv) in g.iter().zip(v) { gsum[*gg as usize] += *vv as f64 / t * 1e4; } } } }
    let auto: Vec<bool> = gnames.iter().map(|s| chr.get(s).map_or(false, |c| c != "chrX" && c != "chrY" && c != "chrM")).collect();
    let mut expressed: Vec<usize> = (0..ng).filter(|&g| gsum[g] > 0.0).collect(); expressed.sort_by(|&x, &y| gsum[x].partial_cmp(&gsum[y]).unwrap());
    let nb = 25; let mut bin = vec![usize::MAX; ng]; let mut bins: Vec<Vec<usize>> = vec![Vec::new(); nb];
    for (r, &g) in expressed.iter().enumerate() { let b = r * nb / expressed.len(); bin[g] = b; if auto[g] { bins[b].push(g); } }
    let mut member: Vec<Vec<u8>> = vec![Vec::new(); ng]; for &g in &yidx { member[g].push(0); }
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(20261121);
    for s in 1..=NP { for &g in &yidx { if bin[g] == usize::MAX { continue; } let pool = &bins[bin[g]]; member[pool[rng.gen_range(0..pool.len())]].push(s as u8); } }
    let mk: Vec<(String, Vec<usize>)> = MARKERS.iter().map(|(n, gs)| (n.to_string(), gs.iter().filter_map(|g| find(g)).collect())).collect();
    let mut pidx = vec![u16::MAX; ng]; let mut pnames = Vec::new();
    for g in 0..ng { let c = chr.get(&gnames[g]).map(|s| s.as_str()).unwrap_or(""); if c == "chrX" || c == "chrY" || ["XIST", "TSIX", "CD99", "TCL1A"].contains(&gnames[g].as_str()) || yidx.contains(&g) { pidx[g] = pnames.len() as u16; pnames.push(gnames[g].clone()); } }
    // per cell records
    let k = NP + 1;
    let mut gid: HashMap<String, u32> = HashMap::new(); let mut gkeys: Vec<String> = Vec::new(); let mut gn: Vec<u64> = Vec::new();
    let mut w = BufWriter::new(std::fs::File::create(format!("{out}.cells.bin"))?);
    let mut wp = BufWriter::new(std::fs::File::create(format!("{out}.panel.bin"))?);
    let mut trough_x: Vec<(usize, String, f64)> = Vec::new(); // male cells: sample, type, log1p Y cp10k
    let mut cell = 0u32;
    for (si, s) in samples.iter().enumerate() {
        for (g, v) in &s.cols {
            let t: f64 = v.iter().map(|&x| x as f64).sum(); if t <= 0.0 { continue; }
            let mut sc = vec![0f64; mk.len()]; let mut rec = vec![0u32; k];
            for (gg, vv) in g.iter().zip(v) { let gg = *gg as usize; for &m in &member[gg] { rec[m as usize] += *vv; }
                for (mi, (_, ms)) in mk.iter().enumerate() { if ms.contains(&gg) { sc[mi] += *vv as f64 / t * 1e4 / ms.len() as f64; } }
                if pidx[gg] != u16::MAX { wp.write_all(&cell.to_le_bytes())?; wp.write_all(&pidx[gg].to_le_bytes())?; wp.write_all(&(*vv).min(65535).to_le_bytes()[..2])?; } }
            let (bi, bv) = sc.iter().enumerate().fold((0, 0.0), |acc, (i, &x)| if x > acc.1 { (i, x) } else { acc });
            let ct = if bv >= 0.5 { mk[bi].0.clone() } else { "unassigned".into() };
            let dis = if s.name.to_lowercase().starts_with("control") { "normal" } else { "kidney disease" };
            let key = format!("{ct}\tkidney\t{}\t{dis}\t{}\tunknown\t10x GEO author matrix\tnucleus", if male[si] { "male" } else { "female" }, s.name);
            let n = gid.len() as u32; let gi = *gid.entry(key.clone()).or_insert_with(|| { gkeys.push(key); gn.push(0); n }); gn[gi as usize] += 1;
            w.write_all(&gi.to_le_bytes())?; w.write_all(&(t as u32).to_le_bytes())?; for x in &rec { w.write_all(&x.to_le_bytes())?; }
            if male[si] { trough_x.push((si, ct, (1.0 + rec[0] as f64 / t * 1e4).ln())); }
            cell += 1;
        }
    }
    w.flush()?; wp.flush()?;
    let mut wg = BufWriter::new(std::fs::File::create(format!("{out}.groups.tsv"))?);
    writeln!(wg, "gid\tcell_type\ttissue\tsex\tdisease\tdonor\tstage\tassay\tsuspension\tn")?;
    for (i, kk) in gkeys.iter().enumerate() { writeln!(wg, "{i}\t{kk}\t{}", gn[i])?; }
    std::fs::write(format!("{out}.panel_genes.txt"), pnames.join("\n") + "\n")?;
    std::fs::write(format!("{out}.info.txt"), format!("cells {cell}\nmatrix 10x h5 GEO\nrecord_u32 {}\nY genes {}\nsex\n{}\n", k + 2, yidx.len(), sexinfo.join("\n")))?;
    // density-trough method (pooled over male cells), Gaussian KDE on a grid
    let xs: Vec<f64> = trough_x.iter().map(|x| x.2).collect();
    let nx = xs.len() as f64; let mean = xs.iter().sum::<f64>() / nx; let sd = (xs.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / nx).sqrt();
    let h = 1.06 * sd * nx.powf(-0.2); let gmax = xs.iter().cloned().fold(0.0, f64::max); let ngr = 512;
    let grid: Vec<f64> = (0..ngr).map(|i| gmax * i as f64 / (ngr - 1) as f64).collect();
    let mut hist = vec![0f64; ngr]; for &x in &xs { hist[((x / gmax) * (ngr - 1) as f64).round() as usize] += 1.0; }
    let dens: Vec<f64> = grid.iter().map(|&g| hist.iter().enumerate().map(|(j, &c)| if c > 0.0 { c * (-0.5 * ((g - grid[j]) / h).powi(2)).exp() } else { 0.0 }).sum()).collect();
    let imax = (ngr / 10..ngr).max_by(|&i, &j| dens[i].partial_cmp(&dens[j]).unwrap()).unwrap();
    let itr = (1..imax).min_by(|&i, &j| dens[i].partial_cmp(&dens[j]).unwrap()).unwrap_or(0);
    let thr = grid[itr];
    let mut wt = BufWriter::new(std::fs::File::create(format!("{out}.trough.tsv"))?);
    writeln!(wt, "# bandwidth {h:.4} main_mode {:.3} trough {thr:.3} (log1p Y CP10k)\nsample\tcell_type\tn_male_cells\ttrough_loy_frac", grid[imax])?;
    let mut agg: HashMap<(usize, String), (f64, f64)> = HashMap::new();
    for (si, ct, x) in &trough_x { let e = agg.entry((*si, ct.clone())).or_insert((0.0, 0.0)); e.0 += 1.0; if *x < thr { e.1 += 1.0; } let e2 = agg.entry((*si, "ALL".into())).or_insert((0.0, 0.0)); e2.0 += 1.0; if *x < thr { e2.1 += 1.0; } }
    let mut keys: Vec<&(usize, String)> = agg.keys().collect(); keys.sort();
    for kk in keys { let (n, l) = agg[kk]; writeln!(wt, "{}\t{}\t{n}\t{:.4}", samples[kk.0].name, kk.1, l / n)?; }
    let _ = samples.iter().map(|s| s.genes_idx.len()).sum::<usize>();
    eprintln!("cells {cell}, genes {ng}, trough {thr:.3}, male cells {}", xs.len());
    Ok(())
}
