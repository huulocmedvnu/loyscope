//! csrcheck: reads the first rows of an h5ad matrix with prefetch_csr and prints nnz, checksums and time; run once with
//! LOY_SERIAL_READ=1 and once without to compare the parallel chunk reader with serial HDF5 reads.
//! usage: cargo run --release --example csrcheck -- <h5ad> [rows=2000000]
use loyscope::{prefetch_csr, row_ranges};
use std::sync::Arc;
fn main() -> anyhow::Result<()> {
    let a: Vec<String> = std::env::args().collect();
    let f = hdf5::File::open(&a[1])?;
    let mpath = if f.group("raw/X").is_ok() { "raw/X" } else { "X" };
    let indptr: Vec<i64> = f.group(mpath)?.dataset("indptr")?.read_raw()?;
    let nr = a
        .get(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(2_000_000usize)
        .min(indptr.len() - 1);
    drop(f);
    let t = std::time::Instant::now();
    let (mut nnz, mut sv, mut si) = (0u64, 0f64, 0u64);
    for item in prefetch_csr(&a[1], mpath, Arc::new(indptr), row_ranges(nr, 50_000, 1)) {
        let (_, _, i, v) = item?;
        nnz += v.len() as u64;
        sv += v.iter().map(|&x| x as f64).sum::<f64>();
        si += i.iter().map(|&x| x as u64).sum::<u64>();
    }
    println!(
        "rows {nr} nnz {nnz} sum_values {sv:.1} sum_indices {si} seconds {:.1}",
        t.elapsed().as_secs_f64()
    );
    Ok(())
}
