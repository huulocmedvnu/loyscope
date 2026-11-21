//! encmulti: one ENCODE 10x multiome library -> RNA records for loyfit (Y and 20 pseudo-Y sets) and chromatin counts for the
//! same nuclei. ATAC barcodes are translated to GEX barcodes with the paired 737K-arc-v1 whitelists (ATAC read is the reverse
//! complement of the whitelist entry). Fragments are streamed once; per nucleus: total, chrY outside PAR (MSY), chrX and
//! counts in 1 Mb autosomal bins, from which 20 control sets of bins matched to the MSY rate are drawn.
//! usage: encmulti <fragments.tsv.gz> <rna_mtx_dir> <atac_wl.gz> <gex_wl.gz> <tss.tsv> <sample> <tissue> <stage> <out_prefix>
//!        encmulti --female <fragments.tsv.gz> <out.tsv>    (background: per-barcode total and MSY fragments, barcodes >= 1000 fragments)
use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use flate2::{write::GzEncoder, Compression};
use rand::{seq::SliceRandom, Rng, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;
use std::collections::HashMap;
use std::io::{BufRead, BufReader, BufWriter, Write};

const YSET: &[&str] = &["RPS4Y1", "ZFY", "USP9Y", "DDX3Y", "UTY", "KDM5D", "EIF1AY", "TXLNGY", "NLGN4Y", "PRKY", "TMSB4Y"];
const NP: usize = 20;
const BIN: u64 = 1_000_000;
const MARKERS: &[(&str, &[&str])] = &[
    ("cardiomyocyte", &["TTN", "RYR2", "MYH7", "TNNT2"]), ("fibroblast", &["DCN", "COL1A2", "LAMA2"]), ("endothelial cell", &["VWF", "PECAM1", "EGFL7"]),
    ("pericyte", &["RGS5", "PDGFRB", "ABCC9"]), ("smooth muscle cell", &["MYH11", "ACTA2"]), ("myeloid cell", &["CD163", "F13A1", "MRC1"]),
    ("lymphocyte", &["PTPRC", "CD247", "IL7R"]), ("hepatocyte", &["ALB", "APOB", "CYP3A4"]), ("adipocyte", &["ADIPOQ", "PLIN1"]),
    ("neuronal", &["NRXN1", "RBFOX3"]), ("epithelial cell", &["EPCAM", "KRT8", "KRT18"]), ("trophoblast", &["KRT7", "GATA3", "PAPPA"])];

fn msy(st: u64) -> bool { !((10_001..=2_781_479).contains(&st) || (56_887_903..=57_217_415).contains(&st)) }
fn rc(s: &str) -> String { s.bytes().rev().map(|b| match b { b'A' => 'T', b'C' => 'G', b'G' => 'C', b'T' => 'A', _ => 'N' }).collect() }
fn gz_lines(p: &str) -> Result<std::io::Lines<BufReader<MultiGzDecoder<std::fs::File>>>> {
    Ok(BufReader::with_capacity(1 << 22, MultiGzDecoder::new(std::fs::File::open(p).with_context(|| p.to_string())?)).lines())
}
fn autosome(ch: &str) -> Option<usize> { ch.strip_prefix("chr")?.parse::<usize>().ok().filter(|&c| (1..=22).contains(&c)) }

fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().collect();
    if a[1] == "--female" {
        let mut m: HashMap<String, (u64, u64)> = HashMap::new();
        for l in gz_lines(&a[2])? { let l = l?; if l.starts_with('#') { continue; } let f: Vec<&str> = l.split('\t').collect(); if f.len() < 4 { continue; }
            let e = m.entry(f[3].to_string()).or_default(); e.0 += 1; if f[0] == "chrY" && msy(f[1].parse().unwrap_or(0)) { e.1 += 1; } }
        let mut w = BufWriter::new(std::fs::File::create(&a[3])?); writeln!(w, "barcode\ttotal\tchrY_msy")?;
        for (b, (t, y)) in m { if t >= 1000 { writeln!(w, "{b}\t{t}\t{y}")?; } }
        return Ok(());
    }
    let (frag, mtx, awl, gwl, tss, sample, tissue, stage, out) = (&a[1], &a[2], &a[3], &a[4], &a[5], &a[6], &a[7], &a[8], &a[9]);
    // whitelist translation ATAC (as read, reverse complement of whitelist) -> GEX
    let wa: Vec<String> = gz_lines(awl)?.map(|l| l.unwrap()).collect(); let wg: Vec<String> = gz_lines(gwl)?.map(|l| l.unwrap()).collect();
    let tr: HashMap<String, usize> = wa.iter().enumerate().map(|(i, s)| (rc(s), i)).collect();
    // RNA matrix (MatrixMarket genes x cells)
    let bcs: Vec<String> = std::fs::read_to_string(format!("{mtx}/barcodes.tsv"))?.lines().map(String::from).collect();
    let feats: Vec<String> = std::fs::read_to_string(format!("{mtx}/features.tsv"))?.lines().map(|l| l.split('\t').nth(1).unwrap_or("").to_string()).collect();
    let ng = feats.len(); let nc = bcs.len();
    let mut cols: Vec<Vec<(u32, u32)>> = vec![Vec::new(); nc];
    let mut hdr = true;
    for l in BufReader::new(std::fs::File::open(format!("{mtx}/matrix.mtx"))?).lines() {
        let l = l?; if l.starts_with('%') { continue; } if hdr { hdr = false; continue; }
        let mut it = l.split_whitespace(); let g: usize = it.next().unwrap().parse()?; let c: usize = it.next().unwrap().parse()?; let v: f64 = it.next().unwrap().parse()?;
        if v > 0.0 { cols[c - 1].push(((g - 1) as u32, v.round() as u32)); }
    }
    let gex_idx: HashMap<&str, usize> = bcs.iter().enumerate().map(|(i, b)| (b.as_str(), i)).collect();
    let find = |s: &str| feats.iter().position(|f| f == s);
    let yidx: Vec<usize> = YSET.iter().filter_map(|s| find(s)).collect();
    // pseudo-Y RNA sets matched on mean CP10k (25 bins) among autosomal genes
    let mut chr: HashMap<String, String> = HashMap::new();
    for l in std::fs::read_to_string(tss)?.lines() { let x: Vec<&str> = l.split('\t').collect(); if x.len() >= 2 { chr.insert(x[0].into(), x[1].into()); } }
    let mut gsum = vec![0f64; ng];
    for c in &cols { let t: f64 = c.iter().map(|x| x.1 as f64).sum(); if t > 0.0 { for &(g, v) in c { gsum[g as usize] += v as f64 / t * 1e4; } } }
    let auto: Vec<bool> = feats.iter().map(|s| chr.get(s).map_or(false, |c| c != "chrX" && c != "chrY" && c != "chrM")).collect();
    let mut expr: Vec<usize> = (0..ng).filter(|&g| gsum[g] > 0.0).collect(); expr.sort_by(|&x, &y| gsum[x].partial_cmp(&gsum[y]).unwrap());
    let nb = 25; let mut bin = vec![usize::MAX; ng]; let mut bins: Vec<Vec<usize>> = vec![Vec::new(); nb];
    for (r, &g) in expr.iter().enumerate() { let b = r * nb / expr.len(); bin[g] = b; if auto[g] { bins[b].push(g); } }
    let mut member: Vec<Vec<u8>> = vec![Vec::new(); ng]; for &g in &yidx { member[g].push(0); }
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(20261121);
    for s in 1..=NP { for &g in &yidx { if bin[g] == usize::MAX { continue; } let p = &bins[bin[g]]; member[p[rng.gen_range(0..p.len())]].push(s as u8); } }
    let mk: Vec<(&str, Vec<usize>)> = MARKERS.iter().map(|(n, gs)| (*n, gs.iter().filter_map(|g| find(g)).collect())).collect();
    // chromatin: one pass over fragments
    let nbins = 3200usize; // 1 Mb bins over chr1..22 (offset table)
    let chrlen_mb = [249, 243, 199, 191, 182, 171, 160, 146, 139, 134, 136, 134, 115, 108, 102, 91, 84, 81, 59, 65, 47, 51];
    let mut off = vec![0usize; 23]; for i in 0..22 { off[i + 1] = off[i] + chrlen_mb[i] + 1; }
    let mut abin = vec![0u16; nc * nbins]; let mut atot = vec![0u32; nc]; let mut ay = vec![0u32; nc]; let mut ax = vec![0u32; nc];
    let mut cache: HashMap<String, Option<usize>> = HashMap::new();
    for l in gz_lines(frag)? {
        let l = l?; if l.starts_with('#') { continue; }
        let mut it = l.split('\t'); let (Some(ch), Some(st), Some(_), Some(bc)) = (it.next(), it.next(), it.next(), it.next()) else { continue };
        let ci = *cache.entry(bc.to_string()).or_insert_with(|| tr.get(bc).and_then(|&i| gex_idx.get(wg[i].as_str()).copied()));
        let Some(c) = ci else { continue };
        atot[c] += 1; let s: u64 = st.parse().unwrap_or(0);
        if ch == "chrY" { if msy(s) { ay[c] += 1; } } else if ch == "chrX" { ax[c] += 1; }
        else if let Some(k) = autosome(ch) { let b = off[k - 1] + (s / BIN) as usize; if b < nbins { let e = &mut abin[c * nbins + b]; *e = e.saturating_add(1); } }
    }
    // pseudo-Y chromatin sets: random bins whose summed mean count matches the MSY mean count
    let yt: f64 = ay.iter().map(|&x| x as f64).sum();
    let bmean: Vec<f64> = (0..nbins).map(|b| (0..nc).map(|c| abin[c * nbins + b] as f64).sum::<f64>()).collect();
    let usable: Vec<usize> = (0..nbins).filter(|&b| bmean[b] > 0.0).collect();
    let mut psets: Vec<Vec<usize>> = Vec::new();
    for _ in 0..NP { let mut order = usable.clone(); order.shuffle(&mut rng); let (mut acc, mut set) = (0.0, Vec::new());
        for b in order { if acc >= yt { break; } acc += bmean[b]; set.push(b); } psets.push(set); }
    // outputs: loyfit-compatible RNA records + chromatin table keyed by GEX barcode
    let mut gid: HashMap<String, u32> = HashMap::new(); let mut gkeys: Vec<String> = Vec::new(); let mut gn: Vec<u64> = Vec::new();
    let mut w = BufWriter::new(std::fs::File::create(format!("{out}.cells.bin"))?);
    let mut wo = GzEncoder::new(BufWriter::new(std::fs::File::create(format!("{out}.obs.tsv.gz"))?), Compression::fast());
    writeln!(wo, "cell_id\tcell_type")?;
    let mut wa_ = GzEncoder::new(BufWriter::new(std::fs::File::create(format!("{out}.atac.tsv.gz"))?), Compression::fast());
    writeln!(wa_, "cell_id\ttotal\tchrY_msy\tchrX\t{}", (1..=NP).map(|i| format!("pseudo{i}")).collect::<Vec<_>>().join("\t"))?;
    for c in 0..nc {
        let t: f64 = cols[c].iter().map(|x| x.1 as f64).sum(); if t <= 0.0 { continue; }
        let mut rec = vec![0u32; NP + 1]; let mut sc = vec![0f64; mk.len()];
        for &(g, v) in &cols[c] { let g = g as usize; for &m in &member[g] { rec[m as usize] += v; } for (mi, (_, ms)) in mk.iter().enumerate() { if ms.contains(&g) { sc[mi] += v as f64 / t * 1e4 / ms.len() as f64; } } }
        let (bi, bv) = sc.iter().enumerate().fold((0, 0.0), |acc, (i, &x)| if x > acc.1 { (i, x) } else { acc });
        let ct = if bv >= 0.5 { mk[bi].0 } else { "unassigned" };
        let key = format!("{ct}\t{tissue}\tmale\tnormal\t{sample}\t{stage}\t10x multiome (ENCODE)\tnucleus");
        let n = gid.len() as u32; let gi = *gid.entry(key.clone()).or_insert_with(|| { gkeys.push(key); gn.push(0); n }); gn[gi as usize] += 1;
        w.write_all(&gi.to_le_bytes())?; w.write_all(&(t as u32).to_le_bytes())?; for x in &rec { w.write_all(&x.to_le_bytes())?; }
        let id = format!("{sample}:{}", bcs[c]); writeln!(wo, "{id}\t{ct}")?;
        let ps: Vec<String> = psets.iter().map(|s| s.iter().map(|&b| abin[c * nbins + b] as u32).sum::<u32>().to_string()).collect();
        writeln!(wa_, "{id}\t{}\t{}\t{}\t{}", atot[c], ay[c], ax[c], ps.join("\t"))?;
    }
    w.flush()?; wo.finish()?; wa_.finish()?;
    let mut wg = BufWriter::new(std::fs::File::create(format!("{out}.groups.tsv"))?);
    writeln!(wg, "gid\tcell_type\ttissue\tsex\tdisease\tdonor\tstage\tassay\tsuspension\tn")?;
    for (i, k) in gkeys.iter().enumerate() { writeln!(wg, "{i}\t{k}\t{}", gn[i])?; }
    let matched = atot.iter().filter(|&&x| x > 0).count();
    std::fs::write(format!("{out}.info.txt"), format!("cells {nc}\nmatrix ENCODE STARsolo GeneFull_Ex50pAS\nrecord_u32 {}\nY genes {}\natac_matched {matched}\n", NP + 3, yidx.len()))?;
    eprintln!("{sample}: RNA cells {nc}, with ATAC fragments {matched}, MSY fragments {yt}, Y genes {}", yidx.len());
    Ok(())
}
