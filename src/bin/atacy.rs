//! atacy: per-nucleus ATAC counts on chrY outside PAR (MSY), chrX and all features from an ATAC h5ad (CSR nuclei x peaks
//! or CSC), plus a dump of all obs columns, for joining with RNA LOY calls of the same multiome nuclei.
//! Reading is prefetched on a dedicated thread (HDF5 is serialized) while blocks are reduced in parallel.
//! usage: atacy <atac.h5ad> <out_prefix> [chunk=50000]   writes <out>.atac.tsv.gz (cell_id, total, chrY_msy, chrX), <out>.obs.tsv.gz
use anyhow::{bail, Context, Result};
use flate2::{write::GzEncoder, Compression};
use hdf5::types::{VarLenAscii, VarLenUnicode};
use loyscope::{prefetch_csr, row_ranges};
use rayon::prelude::*;
use std::io::{BufWriter, Write};
use std::sync::Arc;

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
    let (path, out) = (&a[1], &a[2]);
    let chunk: usize = a.get(3).map(|x| x.parse().unwrap()).unwrap_or(50_000);
    let f = hdf5::File::open(path).with_context(|| path.clone())?;
    eprintln!("top: {:?}", f.member_names()?);
    if let Ok(l) = f.group("layers") {
        eprintln!("layers: {:?}", l.member_names()?);
    }
    let mpath = std::env::var("ATAC_MATRIX").unwrap_or_else(|_| {
        if f.group("raw/X").is_ok() {
            "raw/X".into()
        } else {
            "X".into()
        }
    });
    let vbase = if mpath == "raw/X" { "raw/var" } else { "var" };
    let vg = f.group(vbase)?;
    eprintln!("var members: {:?}", vg.member_names()?);
    let vidx: String = vg
        .attr("_index")
        .ok()
        .and_then(|x| x.read_scalar::<VarLenUnicode>().ok())
        .map(|s| s.as_str().to_string())
        .unwrap_or("_index".into());
    let names = read_labels(&f, &format!("{vbase}/{vidx}"))?;
    eprintln!("features {} e.g. {:?}", names.len(), &names[..3.min(names.len())]);
    let cls: Vec<u8> = names
        .iter()
        .map(|n| {
            let t: Vec<&str> = n.split([':', '-', '_']).collect();
            if t.len() < 2 {
                return 0;
            }
            let ch = t[0];
            let st: f64 = t[1].parse().unwrap_or(-1.0);
            if ch == "chrY" || ch == "Y" {
                if (10_001.0..=2_781_479.0).contains(&st) || (56_887_903.0..=57_217_415.0).contains(&st) {
                    0
                } else {
                    1
                }
            } else if ch == "chrX" || ch == "X" {
                2
            } else {
                0
            }
        })
        .collect();
    eprintln!(
        "chrY MSY features {}, chrX features {}",
        cls.iter().filter(|&&c| c == 1).count(),
        cls.iter().filter(|&&c| c == 2).count()
    );
    let m = f.group(&mpath)?;
    let enc: String = m
        .attr("encoding-type")
        .ok()
        .and_then(|x| x.read_scalar::<VarLenUnicode>().ok())
        .map(|s| s.as_str().to_string())
        .unwrap_or("csr_matrix".into());
    let shape: Vec<i64> = m
        .attr("shape")
        .ok()
        .and_then(|x| x.read_raw::<i64>().ok())
        .unwrap_or_default();
    let indptr = Arc::new(m.dataset("indptr")?.read_raw::<i64>()?);
    let csc = enc == "csc_matrix";
    let nc = if csc { shape[0] as usize } else { indptr.len() - 1 };
    eprintln!("matrix {mpath} {enc} shape {:?}", shape);
    dump_obs(&f, nc, out)?;
    let oidx: String = f
        .group("obs")?
        .attr("_index")
        .ok()
        .and_then(|x| x.read_scalar::<VarLenUnicode>().ok())
        .map(|s| s.as_str().to_string())
        .unwrap_or("_index".into());
    let ids = read_labels(&f, &format!("obs/{oidx}"))?;
    let mut acc = vec![(0f64, 0f64, 0f64); nc];
    if csc {
        let nf = indptr.len() - 1;
        for item in prefetch_csr(path, &mpath, indptr.clone(), row_ranges(nf, 20_000, 1)) {
            let (f0, p, ix, v) = item?;
            for j in 0..p.len() - 1 {
                let c = cls[f0 + j];
                for q in p[j]..p[j + 1] {
                    let e = &mut acc[ix[q] as usize];
                    let val = v[q] as f64;
                    e.0 += val;
                    if c == 1 {
                        e.1 += val
                    } else if c == 2 {
                        e.2 += val
                    }
                }
            }
        }
    } else {
        for item in prefetch_csr(path, &mpath, indptr.clone(), row_ranges(nc, chunk, 1)) {
            let (r0, p, ix, v) = item?;
            let rec: Vec<(f64, f64, f64)> = (0..p.len() - 1)
                .into_par_iter()
                .map(|c| {
                    let (mut t, mut y, mut x) = (0.0, 0.0, 0.0);
                    for q in p[c]..p[c + 1] {
                        let val = v[q] as f64;
                        t += val;
                        match cls[ix[q] as usize] {
                            1 => y += val,
                            2 => x += val,
                            _ => {}
                        }
                    }
                    (t, y, x)
                })
                .collect();
            for (c, r) in rec.into_iter().enumerate() {
                acc[r0 + c] = r;
            }
        }
    }
    let mut w = GzEncoder::new(
        BufWriter::new(std::fs::File::create(format!("{out}.atac.tsv.gz"))?),
        Compression::fast(),
    );
    writeln!(w, "cell_id\ttotal\tchrY_msy\tchrX")?;
    for (c, (t, y, x)) in acc.iter().enumerate() {
        writeln!(w, "{}\t{t}\t{y}\t{x}", ids[c])?;
    }
    w.finish()?;
    eprintln!("{path}: {nc} nuclei done");
    Ok(())
}
