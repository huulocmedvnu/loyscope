//! atacmix: fraction of nuclei without Y from chromatin fragment counts of one library (fragy output).
//! MSY fragments T_y of a nucleus with T fragments follow f * Pois(b T) + (1 - f) * NB(m T, k), with b the MSY fragment
//! rate per fragment of nuclei from a female library of the same study (fixed); fitted with the same estimator as RNA.
//! usage: atacmix <female_fragy.tsv> <male_fragy.tsv> <out.json> [min_fragments=1000]
use anyhow::{Context, Result};
use loyscope::fit_bg;

fn load(path: &str, min: f64) -> Result<(Vec<f64>, Vec<f64>)> {
    let text = std::fs::read_to_string(path).with_context(|| path.to_string())?;
    let mut lines = text.lines();
    let hdr: Vec<&str> = lines.next().context("empty")?.split('\t').collect();
    let col = |n: &str| {
        hdr.iter()
            .position(|h| *h == n)
            .with_context(|| format!("{path}: no column {n}"))
    };
    let (it, iy) = (col("atac_frag")?, col("atac_chrY_msy")?);
    let (mut t, mut y) = (Vec::new(), Vec::new());
    for l in lines {
        let f: Vec<&str> = l.split('\t').collect();
        let tt: f64 = f[it].parse()?;
        if tt >= min {
            t.push(tt);
            y.push(f[iy].parse()?);
        }
    }
    Ok((t, y))
}

fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().collect();
    let min: f64 = a.get(4).and_then(|s| s.parse().ok()).unwrap_or(1000.0);
    let (tf, yf) = load(&a[1], min)?;
    let (tm, ym) = load(&a[2], min)?;
    let bg = yf.iter().sum::<f64>() / tf.iter().sum::<f64>();
    let (f, m, k) = fit_bg(&tm, &ym, bg);
    let zero = ym.iter().filter(|&&v| v == 0.0).count() as f64 / ym.len() as f64;
    let json = format!(
        "{{\"female_nuclei\": {}, \"male_nuclei\": {}, \"background_per_fragment\": {bg:.4e}, \"intact_rate_per_fragment\": {m:.4e}, \"k\": {k}, \"loy_fraction\": {f:.5}, \"zero_msy_fraction\": {zero:.5}}}\n",
        tf.len(),
        tm.len()
    );
    std::fs::write(&a[3], &json)?;
    eprint!("{json}");
    Ok(())
}
