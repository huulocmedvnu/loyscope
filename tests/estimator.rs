//! Recovery of the LOY fraction from simulated detection-limited counts.
use loyscope::{fit, fit_amb, ln_gamma, row_ranges};
use rand::{Rng, SeedableRng};
use rand_distr::{Distribution, Gamma, Poisson};
use rand_xoshiro::Xoshiro256PlusPlus;

/// cells with library sizes 2,000-8,000; intact Y ~ NB(mu N, k), lost Y ~ Poisson(rho mu N)
fn simulate(n_cells: usize, f: f64, mu: f64, k: f64, rho: f64, seed: u64) -> (Vec<f64>, Vec<f64>) {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    let (mut n, mut y) = (Vec::with_capacity(n_cells), Vec::with_capacity(n_cells));
    for _ in 0..n_cells {
        let lib = rng.gen_range(2000.0..8000.0);
        let lam = if rng.gen::<f64>() < f {
            rho * mu * lib
        } else {
            Gamma::new(k, mu * lib / k).unwrap().sample(&mut rng)
        };
        let c = if lam > 0.0 {
            Poisson::new(lam).unwrap().sample(&mut rng)
        } else {
            0.0
        };
        n.push(lib);
        y.push(c);
    }
    (n, y)
}

#[test]
fn ln_gamma_known_values() {
    assert!((ln_gamma(5.0) - 24f64.ln()).abs() < 1e-10);
    assert!((ln_gamma(0.5) - std::f64::consts::PI.sqrt().ln()).abs() < 1e-10);
    assert!((ln_gamma(1.0)).abs() < 1e-10);
}

#[test]
fn recovers_fraction_without_ambient() {
    for &f in &[0.05, 0.2, 0.5] {
        let (n, y) = simulate(20_000, f, 1e-3, 4.0, 0.0, 11);
        let (fh, muh, _) = fit(&n, &y);
        assert!((fh - f).abs() < 0.03, "f {f} estimated {fh}");
        assert!((muh / 1e-3 - 1.0).abs() < 0.1, "mu estimated {muh}");
    }
}

#[test]
fn intact_cells_give_zero() {
    let (n, y) = simulate(20_000, 0.0, 1e-3, 4.0, 0.0, 12);
    let (fh, _, _) = fit(&n, &y);
    assert!(fh < 0.02, "estimated {fh} for intact cells");
}

#[test]
fn zero_read_fraction_overstates_loss_at_low_depth() {
    // about two expected reads per intact cell: many intact cells have no read
    let (n, y) = simulate(20_000, 0.1, 4e-4, 4.0, 0.0, 13);
    let zero = y.iter().filter(|&&v| v == 0.0).count() as f64 / y.len() as f64;
    let (fh, _, _) = fit(&n, &y);
    assert!(zero > 0.25, "zero fraction {zero}");
    assert!((fh - 0.1).abs() < 0.04, "estimated {fh}");
}

#[test]
fn recovers_fraction_with_ambient_floor() {
    let (n, y) = simulate(20_000, 0.2, 1e-3, 4.0, 0.1, 14);
    let (fh, _, _) = fit_amb(&n, &y, 0.1);
    assert!((fh - 0.2).abs() < 0.03, "estimated {fh}");
}

#[test]
fn all_zero_counts_return_full_loss() {
    let n = vec![3000.0; 100];
    let y = vec![0.0; 100];
    assert_eq!(fit(&n, &y).0, 1.0);
}

#[test]
fn row_ranges_cover_all_rows() {
    let r = row_ranges(1003, 100, 1);
    assert_eq!(r.first(), Some(&(0, 100)));
    assert_eq!(r.last(), Some(&(1000, 1003)));
    assert_eq!(r.iter().map(|(a, b)| b - a).sum::<usize>(), 1003);
    assert_eq!(row_ranges(1000, 100, 3).len(), 4);
}
