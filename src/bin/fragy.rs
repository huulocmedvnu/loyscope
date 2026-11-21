//! fragy: chromatin ground truth from a 10x ATAC fragments file joined to the RNA matrix of the same multiome library.
//! Streams fragments.tsv.gz once (chrom, start, end, barcode, count): per barcode total fragments, fragments on chrY
//! outside PAR (MSY) and on chrX. Reads the RNA filtered_feature_bc_matrix.h5 for per-barcode total UMI and Y-set UMI.
//! usage: fragy <fragments.tsv.gz> <rna.h5> <out.tsv>   (barcodes absent from the RNA matrix are skipped)
use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use hdf5::types::{FixedAscii, VarLenAscii, VarLenUnicode};
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Write};

const YSET: &[&str] = &[
    "RPS4Y1", "ZFY", "USP9Y", "DDX3Y", "UTY", "KDM5D", "EIF1AY", "TXLNGY", "NLGN4Y", "PRKY", "TMSB4Y",
];

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
fn msy(start: u64) -> bool {
    !((10_001..=2_781_479).contains(&start) || (56_887_903..=57_217_415).contains(&start))
}

fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().collect();
    let (frag, rna, out) = (&a[1], &a[2], &a[3]);
    // RNA per barcode
    let f = hdf5::File::open(rna)?;
    let m = f.group("matrix")?;
    let bcs = strings(&m.dataset("barcodes")?)?;
    let names = strings(&m.dataset("features/name")?)?;
    let ftype = strings(&m.dataset("features/feature_type")?)?;
    let isy: Vec<bool> = names.iter().map(|n| YSET.contains(&n.as_str())).collect();
    let gex: Vec<bool> = ftype.iter().map(|t| t == "Gene Expression").collect();
    let data: Vec<i64> = m.dataset("data")?.read_raw::<i64>()?;
    let ind: Vec<i64> = m.dataset("indices")?.read_raw::<i64>()?;
    let ip: Vec<i64> = m.dataset("indptr")?.read_raw::<i64>()?;
    let mut rec: HashMap<String, [f64; 5]> = HashMap::new(); // rna_total, rna_y, atac_total, atac_y, atac_x
    for c in 0..bcs.len() {
        let (mut t, mut y) = (0.0, 0.0);
        for q in ip[c] as usize..ip[c + 1] as usize {
            let g = ind[q] as usize;
            if gex[g] {
                t += data[q] as f64;
                if isy[g] {
                    y += data[q] as f64;
                }
            }
        }
        rec.insert(bcs[c].clone(), [t, y, 0.0, 0.0, 0.0]);
    }
    eprintln!("RNA barcodes {}", rec.len());
    // fragments
    let r = BufReader::with_capacity(1 << 22, MultiGzDecoder::new(std::fs::File::open(frag)?));
    let (mut nf, mut hit) = (0u64, 0u64);
    for l in r.lines() {
        let l = l?;
        if l.starts_with('#') {
            continue;
        }
        let mut it = l.split('\t');
        let (Some(ch), Some(st), Some(_en), Some(bc)) = (it.next(), it.next(), it.next(), it.next()) else {
            continue;
        };
        nf += 1;
        if let Some(e) = rec.get_mut(bc) {
            hit += 1;
            e[2] += 1.0;
            if ch == "chrY" {
                if msy(st.parse().unwrap_or(0)) {
                    e[3] += 1.0;
                }
            } else if ch == "chrX" {
                e[4] += 1.0;
            }
        }
    }
    eprintln!("fragments {nf}, in RNA barcodes {hit}");
    let mut w = std::io::BufWriter::new(std::fs::File::create(out)?);
    writeln!(w, "barcode\trna_umi\trna_y\tatac_frag\tatac_chrY_msy\tatac_chrX")?;
    for (b, e) in &rec {
        if e[2] > 0.0 {
            writeln!(w, "{b}\t{}\t{}\t{}\t{}\t{}", e[0], e[1], e[2], e[3], e[4])?;
        }
    }
    Ok(())
}
