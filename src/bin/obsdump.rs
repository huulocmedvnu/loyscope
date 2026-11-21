//! obsdump: write obs index + selected categorical columns of an h5ad as TSV (streams nothing heavy).
//! usage: obsdump <file.h5ad> <out.tsv> [extra_col ...]
use anyhow::{bail, Result};
use hdf5::types::{VarLenAscii, VarLenUnicode};
use std::io::{BufWriter, Write};
fn read_strings(ds: &hdf5::Dataset) -> Result<Vec<String>> {
    if let Ok(v) = ds.read_raw::<VarLenUnicode>() { return Ok(v.iter().map(|s| s.as_str().to_string()).collect()); }
    if let Ok(v) = ds.read_raw::<VarLenAscii>() { return Ok(v.iter().map(|s| s.as_str().to_string()).collect()); }
    bail!("cannot read strings from {}", ds.name())
}
fn read_labels(f: &hdf5::File, path: &str) -> Result<Vec<String>> {
    if let Ok(g) = f.group(path) {
        let cats = read_strings(&g.dataset("categories")?)?;
        let codes: Vec<i64> = g.dataset("codes")?.read_raw::<i64>()?;
        return Ok(codes.iter().map(|&c| if c < 0 { "NA".into() } else { cats[c as usize].clone() }).collect());
    }
    if let Ok(d) = f.dataset(path) { if let Ok(s) = read_strings(&d) { return Ok(s); } if let Ok(v) = d.read_raw::<i64>() { return Ok(v.iter().map(|x| x.to_string()).collect()); } }
    bail!("cannot read {path}")
}
fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().collect();
    let f = hdf5::File::open(&a[1])?;
    let obs = f.group("obs")?;
    let idx_name: String = obs.attr("_index").ok().and_then(|x| x.read_scalar::<VarLenUnicode>().ok()).map(|s| s.as_str().to_string()).unwrap_or("_index".into());
    let mut cols = vec![idx_name.clone(), "cell_type".into(), "donor_id".into(), "sex".into(), "development_stage".into(), "disease".into()];
    cols.extend(a[3..].iter().cloned());
    eprintln!("obs members: {:?}", obs.member_names()?);
    let data: Vec<Vec<String>> = cols.iter().map(|c| read_labels(&f, &format!("obs/{c}"))).collect::<Result<_>>()?;
    let mut w = BufWriter::new(std::fs::File::create(&a[2])?);
    writeln!(w, "{}", cols.join("\t"))?;
    for i in 0..data[0].len() { writeln!(w, "{}", data.iter().map(|d| d[i].replace('\t', " ")).collect::<Vec<_>>().join("\t"))?; }
    Ok(())
}
