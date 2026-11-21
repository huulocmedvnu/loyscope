//! loyx: stream a CELLxGENE h5ad and write per-cell total UMI, Y-set UMI and 20 expression-matched autosomal
//! pseudo-Y set UMIs. Memory O(chunk + genes).
//! usage: loyx <file.h5ad> <tss.tsv> <out_prefix> [chunk_cells=50000]
//! writes <out>.groups.tsv (gid, cell_type, tissue, sex, disease, donor, stage, assay, suspension, n),
//!        <out>.cells.bin (per cell u32: gid, N, y, p1..p20), <out>.info.txt
use anyhow::{bail, Context, Result};
use flate2::{write::GzEncoder, Compression};
use hdf5::types::{VarLenAscii, VarLenUnicode};
use loyscope::{prefetch_csr, row_ranges};
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::{BufWriter, Write};
use std::sync::Arc;

const YSET: &[&str] = &[
    "RPS4Y1", "ZFY", "USP9Y", "DDX3Y", "UTY", "KDM5D", "EIF1AY", "TXLNGY", "NLGN4Y", "PRKY", "TMSB4Y",
];
const NP: usize = 20;

fn read_strings(ds: &hdf5::Dataset) -> Result<Vec<String>> {
    if let Ok(v) = ds.read_raw::<VarLenUnicode>() {
        return Ok(v.iter().map(|s| s.as_str().to_string()).collect());
    }
    if let Ok(v) = ds.read_raw::<VarLenAscii>() {
        return Ok(v.iter().map(|s| s.as_str().to_string()).collect());
    }
    bail!("cannot read strings from {}", ds.name())
}
fn read_labels(f: &hdf5::File, path: &str) -> Result<Vec<String>> {
    if let Ok(g) = f.group(path) {
        let cats = read_strings(&g.dataset("categories")?)?;
        let codes: Vec<i64> = g.dataset("codes")?.read_raw::<i64>()?;
        return Ok(codes
            .iter()
            .map(|&c| {
                if c < 0 {
                    "NA".into()
                } else {
                    cats[c as usize].clone()
                }
            })
            .collect());
    }
    read_strings(&f.dataset(path)?)
}
struct Csr {
    data: hdf5::Dataset,
    indptr: Vec<i64>,
}

enum Col {
    Cat(Vec<i64>, Vec<String>),
    Str(Vec<String>),
    Num(Vec<f64>),
}
fn read_col(f: &hdf5::File, path: &str) -> Option<Col> {
    if let Ok(g) = f.group(path) {
        let cats = read_strings(&g.dataset("categories").ok()?).ok().or_else(|| {
            g.dataset("categories")
                .ok()?
                .read_raw::<i64>()
                .ok()
                .map(|v| v.iter().map(|x| x.to_string()).collect())
        })?;
        return Some(Col::Cat(g.dataset("codes").ok()?.read_raw::<i64>().ok()?, cats));
    }
    let d = f.dataset(path).ok()?;
    if let Ok(v) = read_strings(&d) {
        return Some(Col::Str(v));
    }
    if let Ok(v) = d.read_raw::<f64>() {
        return Some(Col::Num(v));
    }
    if let Ok(v) = d.read_raw::<i64>() {
        return Some(Col::Num(v.iter().map(|&x| x as f64).collect()));
    }
    if let Ok(v) = d.read_raw::<bool>() {
        return Some(Col::Num(v.iter().map(|&x| x as u8 as f64).collect()));
    }
    None
}
/// all obs columns, one row per cell in matrix order
fn dump_obs(f: &hdf5::File, nc: usize, out: &str) -> Result<()> {
    let obs = f.group("obs")?;
    let mut names = vec!["_index".to_string()];
    let idx: String = obs
        .attr("_index")
        .ok()
        .and_then(|x| x.read_scalar::<VarLenUnicode>().ok())
        .map(|s| s.as_str().to_string())
        .unwrap_or("_index".into());
    let mut cols: Vec<(String, Col)> = Vec::new();
    if let Some(c) = read_col(f, &format!("obs/{idx}")) {
        cols.push(("cell_id".into(), c));
    }
    for n in obs.member_names()? {
        if n == idx || n == "__categories" {
            continue;
        }
        if let Some(c) = read_col(f, &format!("obs/{n}")) {
            let len = match &c {
                Col::Cat(v, _) => v.len(),
                Col::Str(v) => v.len(),
                Col::Num(v) => v.len(),
            };
            if len == nc {
                cols.push((n, c));
            }
        }
    }
    names.clear();
    let mut w = GzEncoder::new(
        BufWriter::new(std::fs::File::create(format!("{out}.obs.tsv.gz"))?),
        Compression::fast(),
    );
    writeln!(
        w,
        "{}",
        cols.iter().map(|c| c.0.clone()).collect::<Vec<_>>().join("\t")
    )?;
    let mut line = String::new();
    for i in 0..nc {
        line.clear();
        for (j, (_, c)) in cols.iter().enumerate() {
            if j > 0 {
                line.push('\t');
            }
            match c {
                Col::Cat(v, cats) => {
                    let k = v[i];
                    if k >= 0 && (k as usize) < cats.len() {
                        line.push_str(&cats[k as usize].replace('\t', " "))
                    } else {
                        line.push_str("NA")
                    }
                }
                Col::Str(v) => line.push_str(&v[i].replace('\t', " ")),
                Col::Num(v) => line.push_str(&format!("{}", v[i])),
            }
        }
        writeln!(w, "{line}")?;
    }
    w.finish()?;
    Ok(())
}

fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().collect();
    let (path, tss, out) = (&a[1], &a[2], &a[3]);
    let chunk: usize = a.get(4).map(|x| x.parse().unwrap()).unwrap_or(50_000);
    let f = hdf5::File::open(path).with_context(|| path.clone())?;
    let mpath = if f.group("raw/X").is_ok() { "raw/X" } else { "X" };
    let syms = read_labels(&f, "raw/var/feature_name").or_else(|_| read_labels(&f, "var/feature_name"))?;
    let ng = syms.len();
    let m = f.group(mpath)?;
    let csr = Csr {
        data: m.dataset("data")?,
        indptr: m.dataset("indptr")?.read_raw::<i64>()?,
    };
    let indptr = Arc::new(csr.indptr.clone());
    let nc = csr.indptr.len() - 1;
    // integer check on the first values
    let probe = csr
        .data
        .read_slice_1d::<f32, _>(0..(csr.indptr[nc] as usize).min(200_000))?;
    if probe.iter().any(|x| (x - x.round()).abs() > 1e-3) {
        bail!("{path}: {mpath} is not raw integer counts");
    }
    let fields = [
        "cell_type",
        "tissue",
        "sex",
        "disease",
        "donor_id",
        "development_stage",
        "assay",
        "suspension_type",
    ];
    let mut cols: Vec<Vec<String>> = fields
        .iter()
        .map(|k| read_labels(&f, &format!("obs/{k}")))
        .collect::<Result<_>>()?;
    // optional LOY_EXTRA = comma list of candidate obs columns (first present is appended to the donor field, e.g. batch/sample)
    let mut extra_used = String::from("none");
    if let Ok(ex) = std::env::var("LOY_EXTRA") {
        for c in ex.split(',') {
            if let Ok(v) = read_labels(&f, &format!("obs/{c}")) {
                for i in 0..nc {
                    cols[4][i] = format!("{}|{}", cols[4][i], v[i]);
                }
                extra_used = c.to_string();
                break;
            }
        }
    }
    eprintln!(
        "extra grouping column: {extra_used}; obs members: {:?}",
        f.group("obs")?.member_names()?
    );
    let mut gid: HashMap<String, u32> = HashMap::new();
    let mut gkeys: Vec<String> = Vec::new();
    let cg: Vec<u32> = (0..nc)
        .map(|i| {
            let k = cols
                .iter()
                .map(|c| c[i].replace('\t', " "))
                .collect::<Vec<_>>()
                .join("\t");
            let n = gid.len() as u32;
            *gid.entry(k.clone()).or_insert_with(|| {
                gkeys.push(k);
                n
            })
        })
        .collect();
    let male_normal: Vec<bool> = (0..nc)
        .map(|i| cols[2][i] == "male" && cols[3][i] == "normal")
        .collect();
    dump_obs(&f, nc, out)?;
    // chromosome of each gene
    let mut chr: HashMap<String, String> = HashMap::new();
    for l in std::fs::read_to_string(tss)?.lines() {
        let x: Vec<&str> = l.split('\t').collect();
        if x.len() >= 2 {
            chr.insert(x[0].into(), x[1].into());
        }
    }
    let auto: Vec<bool> = syms
        .iter()
        .map(|s| {
            chr.get(s)
                .is_some_and(|c| c != "chrX" && c != "chrY" && c != "chrM")
        })
        .collect();
    // pass 1: mean CP10k per gene over male normal cells
    let mut gsum = vec![0f64; ng];
    // pass 1 reads every 8th block only (expression bins need approximate means)
    let every = std::env::var("LOY_PASS1_EVERY")
        .ok()
        .and_then(|x| x.parse().ok())
        .unwrap_or(8usize);
    for item in prefetch_csr(path, mpath, indptr.clone(), row_ranges(nc, chunk, every)) {
        let (r0, p, ix, v) = item?;
        let r1 = r0 + p.len() - 1;
        let part = (0..r1 - r0)
            .into_par_iter()
            .filter(|&c| male_normal[r0 + c])
            .fold(
                || vec![0f64; ng],
                |mut s, c| {
                    let (lo, hi) = (p[c], p[c + 1]);
                    let tot: f64 = v[lo..hi].iter().map(|&x| x as f64).sum();
                    if tot > 0.0 {
                        for k in lo..hi {
                            s[ix[k] as usize] += v[k] as f64 * 1e4 / tot;
                        }
                    }
                    s
                },
            )
            .reduce(
                || vec![0f64; ng],
                |mut x, y| {
                    for (a, b) in x.iter_mut().zip(y) {
                        *a += b;
                    }
                    x
                },
            );
        for (a, b) in gsum.iter_mut().zip(part) {
            *a += b;
        }
    }
    let yset_env: Vec<String> = std::env::var("LOY_YSET")
        .map(|v| v.split(',').map(String::from).collect())
        .unwrap_or_else(|_| YSET.iter().map(|s| s.to_string()).collect());
    let yidx: Vec<usize> = yset_env
        .iter()
        .filter_map(|s| syms.iter().position(|x| x == s))
        .collect();
    if yidx.is_empty() {
        bail!("no Y genes");
    }
    // 25 expression bins over expressed genes; pseudo sets drawn from autosomal genes of the same bin
    let mut expressed: Vec<usize> = (0..ng).filter(|&g| gsum[g] > 0.0).collect();
    expressed.sort_by(|&x, &y| gsum[x].partial_cmp(&gsum[y]).unwrap());
    let nb = 25;
    let mut bin = vec![usize::MAX; ng];
    let mut bins: Vec<Vec<usize>> = vec![Vec::new(); nb];
    for (r, &g) in expressed.iter().enumerate() {
        let b = r * nb / expressed.len();
        bin[g] = b;
        if auto[g] {
            bins[b].push(g);
        }
    }
    let mut member: Vec<Vec<u8>> = vec![Vec::new(); ng]; // 0 = Y, 1..=NP pseudo
    for &g in &yidx {
        member[g].push(0);
    }
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(20261121);
    for s in 1..=NP {
        for &g in &yidx {
            if bin[g] == usize::MAX {
                continue;
            }
            let pool = &bins[bin[g]];
            member[pool[rng.gen_range(0..pool.len())]].push(s as u8);
        }
    }
    // panel: every chrX/chrY gene + extras; mito genes; pseudobulk slots for male normal groups with >= 50 cells
    const EXTRA: &[&str] = &["XIST", "TSIX", "CD99", "TCL1A"];
    let mut pidx = vec![u16::MAX; ng];
    let mut pnames: Vec<String> = Vec::new();
    for g in 0..ng {
        let c = chr.get(&syms[g]).map(|s| s.as_str()).unwrap_or("");
        if c == "chrX" || c == "chrY" || EXTRA.contains(&syms[g].as_str()) || yidx.contains(&g) {
            pidx[g] = pnames.len() as u16;
            pnames.push(syms[g].clone());
        }
    }
    let is_mt: Vec<bool> = syms.iter().map(|s| s.starts_with("MT-")).collect();
    let mut gcount = vec![0usize; gkeys.len()];
    for &g in &cg {
        gcount[g as usize] += 1;
    }
    let mut slot = vec![u32::MAX; gkeys.len()];
    let mut slots: Vec<usize> = Vec::new();
    for (g, kk) in gkeys.iter().enumerate() {
        let x: Vec<&str> = kk.split('\t').collect();
        if x[2] == "male" && x[3] == "normal" && gcount[g] >= 50 {
            slot[g] = slots.len() as u32;
            slots.push(g);
        }
    }
    let max_slots = (24e9 as usize) / (ng * 2 * 4);
    if slots.len() > max_slots {
        slots.sort_by_key(|&g| std::cmp::Reverse(gcount[g]));
        slots.truncate(max_slots);
        slot = vec![u32::MAX; gkeys.len()];
        for (i, &g) in slots.iter().enumerate() {
            slot[g] = i as u32;
        }
    }
    let mut pb = vec![0f32; slots.len() * 2 * ng];
    let mut pbn = vec![[0u64; 4]; slots.len()]; // n_y0, n_ypos, N_y0, N_ypos
    let mut wp = BufWriter::new(std::fs::File::create(format!("{out}.panel.bin"))?); // u32 cell, u16 panel gene, u16 count
    let mut wq = BufWriter::new(std::fs::File::create(format!("{out}.qc.bin"))?); // u32 n_genes, u32 mt umi
                                                                                  // pass 2: per-cell counts
    let k = NP + 1;
    let mut w = BufWriter::new(std::fs::File::create(format!("{out}.cells.bin"))?);
    let mut gn = vec![0u64; gkeys.len()];
    for item in prefetch_csr(path, mpath, indptr.clone(), row_ranges(nc, chunk, 1)) {
        let (r0, p, ix, v) = item?;
        let r1 = r0 + p.len() - 1;
        let rec: Vec<Vec<u32>> = (0..r1 - r0)
            .into_par_iter()
            .map(|c| {
                let (lo, hi) = (p[c], p[c + 1]);
                let mut s = vec![0u32; k + 2];
                s[0] = cg[r0 + c];
                let mut tot = 0f64;
                for q in lo..hi {
                    tot += v[q] as f64;
                    for &mm in &member[ix[q] as usize] {
                        s[2 + mm as usize] += v[q].round() as u32;
                    }
                }
                s[1] = tot.round() as u32;
                s
            })
            .collect();
        for r in &rec {
            gn[r[0] as usize] += 1;
            for x in r {
                w.write_all(&x.to_le_bytes())?;
            }
        }
        for c in 0..r1 - r0 {
            let (lo, hi) = (p[c], p[c + 1]);
            let cell = (r0 + c) as u32;
            let mut mt = 0u32;
            for q in lo..hi {
                let g = ix[q] as usize;
                let val = v[q].round() as u32;
                if is_mt[g] {
                    mt += val;
                }
                if pidx[g] != u16::MAX {
                    wp.write_all(&cell.to_le_bytes())?;
                    wp.write_all(&pidx[g].to_le_bytes())?;
                    wp.write_all(&(val.min(65535) as u16).to_le_bytes())?;
                }
            }
            wq.write_all(&((hi - lo) as u32).to_le_bytes())?;
            wq.write_all(&mt.to_le_bytes())?;
            let sl = slot[cg[r0 + c] as usize];
            if sl != u32::MAX {
                let side = if rec[c][2] == 0 { 0 } else { 1 };
                let base = (sl as usize * 2 + side) * ng;
                for q in lo..hi {
                    pb[base + ix[q] as usize] += v[q];
                }
                pbn[sl as usize][side] += 1;
                pbn[sl as usize][2 + side] += rec[c][1] as u64;
            }
        }
    }
    wp.flush()?;
    wq.flush()?;
    std::fs::write(format!("{out}.panel_genes.txt"), pnames.join("\n") + "\n")?;
    std::fs::write(format!("{out}.genes.txt"), syms.join("\n") + "\n")?;
    {
        let mut wb = GzEncoder::new(
            BufWriter::new(std::fs::File::create(format!("{out}.pb.f32.gz"))?),
            Compression::fast(),
        );
        for x in &pb {
            wb.write_all(&x.to_le_bytes())?;
        }
        wb.finish()?;
    }
    let mut ws = BufWriter::new(std::fs::File::create(format!("{out}.pb_slots.tsv"))?);
    writeln!(ws, "slot\tgid\tn_y0\tn_ypos\tumi_y0\tumi_ypos")?;
    for (i, &g) in slots.iter().enumerate() {
        writeln!(
            ws,
            "{i}\t{g}\t{}\t{}\t{}\t{}",
            pbn[i][0], pbn[i][1], pbn[i][2], pbn[i][3]
        )?;
    }
    drop(ws);
    let mut wg = BufWriter::new(std::fs::File::create(format!("{out}.groups.tsv"))?);
    writeln!(
        wg,
        "gid\tcell_type\ttissue\tsex\tdisease\tdonor\tstage\tassay\tsuspension\tn"
    )?;
    for (i, kk) in gkeys.iter().enumerate() {
        writeln!(wg, "{i}\t{kk}\t{}", gn[i])?;
    }
    let yfound: Vec<String> = yidx
        .iter()
        .map(|&g| {
            format!(
                "{}:{:.3}",
                syms[g],
                gsum[g] / male_normal.iter().filter(|&&x| x).count().max(1) as f64
            )
        })
        .collect();
    std::fs::write(format!("{out}.info.txt"), format!("cells {nc}\nmatrix {mpath}\nrecord_u32 {}\nY {}\ngenes {ng}\npanel_genes {}\npb_slots {} (layout slot x [y0, ypos] x genes, f32 raw count sums)\n", k + 2, yfound.join(" "), pnames.len(), slots.len()))?;
    eprintln!("{path}: {nc} cells, Y genes {}/{}", yidx.len(), yset_env.len());
    Ok(())
}
