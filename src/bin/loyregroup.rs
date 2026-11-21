//! loyregroup: rebuild loyfit inputs from saved extraction outputs without the h5ad.
//! Regroups cells by extra obs columns appended to the donor field, and optionally recomputes the Y count from the
//! per-gene panel (e.g. RPS4Y1 only). Pseudo-Y columns are kept as extracted.
//! usage: loyregroup <in_prefix> <out_prefix> <extra_cols comma|-> <yset comma|->
use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use std::collections::HashMap;
use std::io::{BufRead, BufReader, BufWriter, Write};

fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().collect();
    let (inp, out, extra, yset) = (&a[1], &a[2], &a[3], &a[4]);
    let info = std::fs::read_to_string(format!("{inp}.info.txt"))?;
    let rl: usize = info.lines().find(|l| l.starts_with("record_u32")).unwrap().split(' ').nth(1).unwrap().parse()?;
    let raw = std::fs::read(format!("{inp}.cells.bin"))?;
    let mut v: Vec<u32> = raw.chunks_exact(4).map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect();
    let nc = v.len() / rl;
    let groups: Vec<String> = std::fs::read_to_string(format!("{inp}.groups.tsv"))?.lines().skip(1).map(|l| l.split('\t').skip(1).take(8).collect::<Vec<_>>().join("\t")).collect();
    // extra columns from obs
    let mut ex: Vec<String> = vec![String::new(); nc];
    if extra != "-" {
        let r = BufReader::new(MultiGzDecoder::new(std::fs::File::open(format!("{inp}.obs.tsv.gz")).context("obs")?));
        let mut lines = r.lines(); let hdr: Vec<String> = lines.next().unwrap()?.split('\t').map(String::from).collect();
        let ci: Vec<usize> = extra.split(',').map(|c| hdr.iter().position(|h| h == c).with_context(|| format!("no obs column {c}"))).collect::<Result<_>>()?;
        for (i, l) in lines.enumerate() { let l = l?; let f: Vec<&str> = l.split('\t').collect(); ex[i] = ci.iter().map(|&j| f[j]).collect::<Vec<_>>().join("|"); }
    }
    // optional Y from panel
    if yset != "-" {
        let pg: Vec<String> = std::fs::read_to_string(format!("{inp}.panel_genes.txt"))?.lines().map(String::from).collect();
        let want: Vec<bool> = pg.iter().map(|g| yset.split(',').any(|y| y == g)).collect();
        eprintln!("Y genes from panel: {}", pg.iter().zip(&want).filter(|x| *x.1).map(|x| x.0.as_str()).collect::<Vec<_>>().join(","));
        for c in 0..nc { v[c * rl + 2] = 0; }
        let p = std::fs::read(format!("{inp}.panel.bin"))?;
        for e in p.chunks_exact(8) { let cell = u32::from_le_bytes([e[0], e[1], e[2], e[3]]) as usize; let g = u16::from_le_bytes([e[4], e[5]]) as usize; let n = u16::from_le_bytes([e[6], e[7]]) as u32; if want[g] { v[cell * rl + 2] += n; } }
    }
    let mut gid: HashMap<String, u32> = HashMap::new(); let mut keys: Vec<String> = Vec::new(); let mut gn: Vec<u64> = Vec::new();
    for c in 0..nc {
        let base = &groups[v[c * rl] as usize]; let k = if extra == "-" { base.clone() } else { let mut f: Vec<&str> = base.split('\t').collect(); let d = format!("{}|{}", f[4], ex[c]); f[4] = &d; f.join("\t") };
        let n = gid.len() as u32; let g = *gid.entry(k.clone()).or_insert_with(|| { keys.push(k); gn.push(0); n }); gn[g as usize] += 1; v[c * rl] = g;
    }
    let mut w = BufWriter::new(std::fs::File::create(format!("{out}.cells.bin"))?); for x in &v { w.write_all(&x.to_le_bytes())?; }
    let mut wg = BufWriter::new(std::fs::File::create(format!("{out}.groups.tsv"))?);
    writeln!(wg, "gid\tcell_type\ttissue\tsex\tdisease\tdonor\tstage\tassay\tsuspension\tn")?;
    for (i, k) in keys.iter().enumerate() { writeln!(wg, "{i}\t{k}\t{}", gn[i])?; }
    std::fs::write(format!("{out}.info.txt"), format!("{info}regroup extra={extra} yset={yset} from {inp}\n"))?;
    eprintln!("{out}: {nc} cells, {} groups", keys.len());
    Ok(())
}
