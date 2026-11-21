//! geox: author count matrices from GEO (dense gene x cell CSV.gz per pool, 10x barcodes) -> same per-cell
//! records as loyx (total UMI, Y-set UMI, 20 pseudo-Y sets), labels joined from a CELLxGENE obs dump.
//! usage: geox <manifest.tsv: pool  counts.csv.gz  barcodes.csv.gz> <obs.tsv> <tss.tsv> <out_prefix>
//! Pass 1: per-gene male sums (for expression bins). Pass 2: per-cell sums. Pools run in parallel.
use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::{BufRead, BufReader, BufWriter, Write};

const YSET: &[&str] = &[
    "RPS4Y1", "ZFY", "USP9Y", "DDX3Y", "UTY", "KDM5D", "EIF1AY", "TXLNGY", "NLGN4Y", "PRKY", "TMSB4Y",
];
const NP: usize = 20;

fn unq(s: &str) -> &str {
    s.trim().trim_matches('"')
}
fn open_gz(p: &str) -> Result<Box<dyn BufRead>> {
    Ok(Box::new(BufReader::with_capacity(
        1 << 22,
        MultiGzDecoder::new(std::fs::File::open(p).with_context(|| p.to_string())?),
    )))
}

/// parse a dense CSV row "GENE",v1,v2,...; calls f(cell_index, value) for non-zero values
fn parse_row(line: &str, mut f: impl FnMut(usize, u32)) -> String {
    let b = line.as_bytes();
    let mut i = 0;
    let mut q = false;
    while i < b.len() {
        if b[i] == b'"' {
            q = !q;
        } else if b[i] == b',' && !q {
            break;
        }
        i += 1;
    }
    let gene = unq(&line[..i]).to_string();
    let (mut c, mut v, mut any, mut frac) = (0usize, 0u64, false, false);
    for &x in &b[(i + 1).min(b.len())..] {
        match x {
            b',' => {
                if any && v > 0 {
                    f(c, v as u32);
                }
                c += 1;
                v = 0;
                any = false;
                frac = false;
            }
            b'0'..=b'9' => {
                if !frac {
                    v = v * 10 + (x - b'0') as u64;
                }
                any = true;
            }
            b'.' => {
                frac = true;
            }
            _ => {}
        }
    }
    if any && v > 0 {
        f(c, v as u32);
    }
    gene
}

fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().collect();
    let (man, obsf, tss, out) = (&a[1], &a[2], &a[3], &a[4]);
    let pools: Vec<(String, String, String)> = std::fs::read_to_string(man)?
        .lines()
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (f[0].into(), f[1].into(), f[2].into())
        })
        .collect();
    // obs: index, cell_type, donor_id, sex, development_stage, disease
    // join key: barcode core (suffix after '-' dropped) + donor id; donor per barcode from the GEO Individual_Barcodes map
    let core = |b: &str| b.split('-').next().unwrap_or(b).to_string();
    let obs_txt = std::fs::read_to_string(obsf)?;
    let mut obs: HashMap<String, Vec<String>> = HashMap::new();
    for l in obs_txt.lines().skip(1) {
        let f: Vec<&str> = l.split('\t').collect();
        obs.insert(
            format!("{}|{}", core(f[0]), f[2]),
            f[1..6].iter().map(|s| s.to_string()).collect(),
        );
    }
    let header = |p: &str| -> Result<Vec<String>> {
        let mut r = open_gz(p)?;
        let mut h = String::new();
        r.read_line(&mut h)?;
        Ok(h.trim_end()
            .split(',')
            .skip(1)
            .map(|s| unq(s).to_string())
            .collect())
    };
    let mut gid: HashMap<String, u32> = HashMap::new();
    let mut gkeys: Vec<String> = Vec::new();
    let mut cg: Vec<Vec<u32>> = Vec::new();
    let mut male: Vec<Vec<bool>> = Vec::new();
    let (mut nmatch, mut ntot) = (0, 0);
    for (_, cf, bf) in &pools {
        let mut ind: HashMap<String, String> = HashMap::new();
        let mut r = open_gz(bf)?;
        let mut l = String::new();
        r.read_line(&mut l)?;
        loop {
            l.clear();
            if r.read_line(&mut l)? == 0 {
                break;
            }
            let x: Vec<&str> = l.trim_end().split(',').collect();
            if x.len() >= 2 {
                ind.insert(core(unq(x[1])), unq(x[0]).to_string());
            }
        }
        let h = header(cf)?;
        let mut g = Vec::with_capacity(h.len());
        let mut m = Vec::with_capacity(h.len());
        for b in &h {
            ntot += 1;
            let o = ind
                .get(&core(b))
                .and_then(|d| obs.get(&format!("{}|{}", core(b), d)));
            match o {
                Some(o) => {
                    nmatch += 1;
                    let k = format!(
                        "{}\tblood\t{}\t{}\t{}\t{}\t10x 3' v2 (GEO author matrix)\tcell",
                        o[0], o[2], o[4], o[1], o[3]
                    );
                    let n = gid.len() as u32;
                    g.push(*gid.entry(k.clone()).or_insert_with(|| {
                        gkeys.push(k);
                        n
                    }));
                    m.push(o[2] == "male" && o[4] == "normal");
                }
                None => {
                    g.push(u32::MAX);
                    m.push(false);
                }
            }
        }
        cg.push(g);
        male.push(m);
    }
    eprintln!("matched {nmatch} / {ntot} matrix barcodes to obs");
    // pass 1: per gene male sums and male total
    let p1: Vec<(HashMap<String, f64>, f64)> = pools
        .par_iter()
        .enumerate()
        .map(|(pi, (_, cf, _))| -> Result<(HashMap<String, f64>, f64)> {
            let mut r = open_gz(cf)?;
            let mut line = String::new();
            r.read_line(&mut line)?;
            let mut gs: HashMap<String, f64> = HashMap::new();
            let mut tot = 0.0;
            loop {
                line.clear();
                if r.read_line(&mut line)? == 0 {
                    break;
                }
                let mut s = 0.0;
                let gene = parse_row(&line, |c, v| {
                    if male[pi][c] {
                        s += v as f64;
                    }
                });
                tot += s;
                *gs.entry(gene).or_default() += s;
            }
            Ok((gs, tot))
        })
        .collect::<Result<_>>()?;
    let mut gsum: HashMap<String, f64> = HashMap::new();
    let mut mtot = 0.0;
    for (g, t) in p1 {
        mtot += t;
        for (k, v) in g {
            *gsum.entry(k).or_default() += v;
        }
    }
    let mut chr: HashMap<String, String> = HashMap::new();
    for l in std::fs::read_to_string(tss)?.lines() {
        let x: Vec<&str> = l.split('\t').collect();
        if x.len() >= 2 {
            chr.insert(x[0].into(), x[1].into());
        }
    }
    let mut genes: Vec<String> = gsum.keys().cloned().collect();
    genes.sort();
    let yg: Vec<String> = YSET
        .iter()
        .filter(|s| gsum.contains_key(**s))
        .map(|s| s.to_string())
        .collect();
    let mut expressed: Vec<&String> = genes.iter().filter(|g| gsum[*g] > 0.0).collect();
    expressed.sort_by(|x, y| gsum[*x].partial_cmp(&gsum[*y]).unwrap());
    let nb = 25;
    let mut bin: HashMap<&String, usize> = HashMap::new();
    let mut bins: Vec<Vec<&String>> = vec![Vec::new(); nb];
    for (r, g) in expressed.iter().enumerate() {
        let b = r * nb / expressed.len();
        bin.insert(g, b);
        if chr
            .get(*g)
            .map_or(false, |c| c != "chrX" && c != "chrY" && c != "chrM")
        {
            bins[b].push(g);
        }
    }
    let mut member: HashMap<String, Vec<u8>> = HashMap::new();
    for g in &yg {
        member.entry(g.clone()).or_default().push(0);
    }
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(20261121);
    for s in 1..=NP {
        for g in &yg {
            if let Some(&b) = bin.get(g) {
                let pool = &bins[b];
                member
                    .entry(pool[rng.gen_range(0..pool.len())].to_string())
                    .or_default()
                    .push(s as u8);
            }
        }
    }
    // pass 2: per cell totals and set sums
    let k = NP + 1;
    let recs: Vec<Vec<u32>> = pools
        .par_iter()
        .enumerate()
        .map(|(pi, (_, cf, _))| -> Result<Vec<u32>> {
            let nc = cg[pi].len();
            let mut s = vec![0u32; nc * (k + 1)];
            let mut r = open_gz(cf)?;
            let mut line = String::new();
            r.read_line(&mut line)?;
            loop {
                line.clear();
                if r.read_line(&mut line)? == 0 {
                    break;
                }
                let b = line.as_bytes();
                let mut i = 0;
                let mut q = false;
                while i < b.len() {
                    if b[i] == b'"' {
                        q = !q;
                    } else if b[i] == b',' && !q {
                        break;
                    }
                    i += 1;
                }
                let mem = member.get(unq(&line[..i])).cloned().unwrap_or_default();
                parse_row(&line, |c, v| {
                    s[c * (k + 1)] += v;
                    for &m in &mem {
                        s[c * (k + 1) + 1 + m as usize] += v;
                    }
                });
            }
            let mut o = Vec::with_capacity(nc * (k + 2));
            for c in 0..nc {
                if cg[pi][c] == u32::MAX {
                    continue;
                }
                o.push(cg[pi][c]);
                o.extend_from_slice(&s[c * (k + 1)..(c + 1) * (k + 1)]);
            }
            Ok(o)
        })
        .collect::<Result<_>>()?;
    let mut w = BufWriter::new(std::fs::File::create(format!("{out}.cells.bin"))?);
    let mut gn = vec![0u64; gkeys.len()];
    for r in &recs {
        for rec in r.chunks_exact(k + 2) {
            gn[rec[0] as usize] += 1;
            for x in rec {
                w.write_all(&x.to_le_bytes())?;
            }
        }
    }
    let mut wg = BufWriter::new(std::fs::File::create(format!("{out}.groups.tsv"))?);
    writeln!(
        wg,
        "gid\tcell_type\ttissue\tsex\tdisease\tdonor\tstage\tassay\tsuspension\tn"
    )?;
    for (i, kk) in gkeys.iter().enumerate() {
        writeln!(wg, "{i}\t{kk}\t{}", gn[i])?;
    }
    let yf: Vec<String> = yg
        .iter()
        .map(|g| format!("{g}:{:.3}", gsum[g] / mtot * 1e4))
        .collect();
    std::fs::write(
        format!("{out}.info.txt"),
        format!(
            "cells {nmatch}\nmatrix GEO author CSV ({} genes)\nrecord_u32 {}\nY {}\n",
            genes.len(),
            k + 2,
            yf.join(" ")
        ),
    )?;
    eprintln!("done: genes {}, Y genes {}/{}", genes.len(), yg.len(), YSET.len());
    Ok(())
}
