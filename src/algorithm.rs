extern crate nalgebra as na;

use crate::problem::Problem;
use na::{DMatrix, DVector};
use nalgebra::SymmetricEigen;
use rand_distr::{
    num_traits::{Float, Pow},
    Distribution, Normal,
};

// Algorithm from: https://en.wikipedia.org/wiki/CMA-ES
pub fn cma_es(problem: &mut Problem, max_budget: usize) {
    let n = problem.dimension();
    let mut xmean = random_normal_vector(n, 0.0, 1.0);
    let mut sigma = 0.3;
    let stopfitness = 1e-10;
    let stopeval = max_budget;

    let lambda = 4 + (3.0 * f64::ln(n as f64)).floor() as usize;
    let mu = lambda / 2;
    let weights = calculate_weights(mu);
    let mueff: f64 =
        weights.iter().sum::<f64>().pow(2) / weights.iter().map(|w| w.pow(2)).sum::<f64>();

    let cc: f64 = (4.0 + mueff / n as f64) / (n as f64 + 4.0 + 2.0 * mueff / n as f64);
    let cs: f64 = (mueff + 2.0) / (n as f64 + mueff + 5.0);
    let c1: f64 = 2.0 / ((n as f64 + 1.3).pow(2) + mueff);
    let cmu = f64::min(
        1.0 - c1,
        2.0 * (mueff - 2.0 + 1.0 / mueff) / ((n as f64 + 2.0).pow(2) + mueff),
    );
    let damps = 1.0 + 2.0 * f64::max(0.0, f64::sqrt((mueff - 1.0) / (n as f64 + 1.0)) - 1.0) + cs;

    let mut pc = DVector::<f64>::zeros(n);
    let mut ps = DVector::<f64>::zeros(n);
    let mut b = DMatrix::<f64>::identity(n, n);
    let mut d = DVector::<f64>::from_element(n, 1.0);
    let mut c = &b * DMatrix::<f64>::from_diagonal(&d).map(|d| d.pow(2)) * b.transpose(); // TODO: check if this is correct
    let mut invsqrt_c =
        &b * DMatrix::<f64>::from_diagonal(&d.map(|d| 1.0 / f64::sqrt(d))) * b.transpose();
    let mut eigeneval = 0;
    let chi_n =
        f64::sqrt(n as f64) * (1.0 - 1.0 / (4.0 * n as f64) + 1.0 / (21.0 * (n as f64).pow(2)));

    let mut counteval = 0;
    while counteval < stopeval {
        let mut arx = DMatrix::zeros(n, lambda);
        let mut arfitness = DVector::<f64>::zeros(lambda);

        for k in 0..lambda {
            // sigma moved to the end of the line here
            let col = &xmean + &b * (d.component_mul(&random_normal_vector(n, 0.0, 1.0))) * sigma;
            arx.set_column(k, &col);
            let mut eval = vec![0.0; 1];
            evaluate(problem, col.as_slice(), &mut eval);
            // problem.evaluate_function(col.as_slice(), &mut eval);
            arfitness[k] = eval.iter().sum();
            counteval += 1;
        }

        let mut arindex = vec![0; lambda];
        (arfitness, arindex) = sort_vec_and_indices(arfitness);
        let xold = xmean.clone();
        // TODO: check this line
        let index_vector = DVector::from_vec(arindex.as_slice()[0..mu].to_vec());
        xmean = select_rows_from_matrix(&arx, &index_vector) * weights.clone();

        ps = ps * (1.0 - cs)
            + f64::sqrt(cs * (2.0 - cs) * mueff) * &invsqrt_c * (&xmean - &xold) / sigma;
        let hsig_bool = ps.norm()
            / f64::sqrt(1.0 - (1.0 - cs).pow(2.0 * counteval as f64 / lambda as f64))
            / chi_n
            < 1.4 + 2.0 / (n as f64 + 1.0);
        let hsig = if hsig_bool { 1.0 } else { 0.0 };
        pc = (1.0 - cc) * pc + hsig * f64::sqrt(cc * (2.0 - cc) * mueff) * (&xmean - &xold) / sigma;

        // TODO: check this line
        let artmp =
            (1.0 / sigma) * (select_rows_from_matrix(&arx, &index_vector) - repmap(xold, mu));
        c = (1.0 - c1 - cmu) * &c
            + (&pc * &pc.transpose() * c1 + &c * (1.0 - hsig) * cc * (2.0 - cc))
            + cmu * &artmp * DMatrix::from_diagonal(&weights) * artmp.transpose();

        sigma = sigma * f64::exp((cs / damps) * (ps.norm() / chi_n - 1.0));

        if (counteval - eigeneval) as f64 > lambda as f64 / (c1 + cmu) as f64 / n as f64 / 10.0 {
            eigeneval = counteval;
            // TODO: check this code
            c = triu(&c, 0) + triu(&c, 1).transpose();
            let eig = SymmetricEigen::new(c.clone());
            b = eig.eigenvectors;
            d = eig.eigenvalues.map(|d| d.sqrt());
            let d_inv_sqrt = DMatrix::from_diagonal(&d.map(|x| 1.0 / x));
            invsqrt_c = &b * d_inv_sqrt * b.transpose();
        }

        if arfitness[0] <= stopfitness || d.max() > 1e7 * d.min() {
            break;
        }
    }

    // let xmin = arx.column(arindex[0]).clone();
    // evaluate(problem, &xmin, &mut vec![0.0; 1]);
}

fn calculate_weights(mu: usize) -> DVector<f64> {
    let mut weights = DVector::zeros(mu);
    for i in 0..mu {
        weights[i] = f64::ln(mu as f64 + 0.5) - f64::ln(i as f64 + 1.0);
    }
    let sum = weights.iter().sum::<f64>();
    weights /= sum;
    weights
}

// TODO: Check this function
fn random_normal_vector(size: usize, mean: f64, std_dev: f64) -> DVector<f64> {
    let normal = Normal::new(mean, std_dev).unwrap();
    let mut rng = rand::thread_rng();
    let vec = (0..size).map(|_| normal.sample(&mut rng)).collect();
    DVector::from_vec(vec)
}

fn sort_vec_and_indices(vec: DVector<f64>) -> (DVector<f64>, Vec<usize>) {
    let mut indexed_vec: Vec<(usize, &f64)> = vec
        .iter()
        .enumerate()
        .map(|(index, value)| (index, value))
        .collect();
    indexed_vec.sort_by(|a, b| {
        a.1.partial_cmp(&b.1)
            .unwrap_or_else(|| a.1.is_nan().cmp(&b.1.is_nan()))
    });
    let (indices, sorted_values): (Vec<usize>, Vec<f64>) = indexed_vec.into_iter().unzip();

    (DVector::from_vec(sorted_values), indices)
}

// TODO: check this function
fn select_rows_from_matrix(matrix: &DMatrix<f64>, indices: &DVector<usize>) -> DMatrix<f64> {
    let mut selected_rows = DMatrix::zeros(matrix.nrows(), indices.len());
    for (i, index) in indices.iter().enumerate() {
        selected_rows.set_column(i, &matrix.column(*index));
    }
    selected_rows
}

// TODO: check this function
fn repmap(vector: DVector<f64>, cols: usize) -> DMatrix<f64> {
    let mut matrix = DMatrix::zeros(vector.len(), cols);
    for i in 0..cols {
        matrix.set_column(i, &vector);
    }
    matrix
}

// TODO: check this function
fn triu(matrix: &DMatrix<f64>, offset: usize) -> DMatrix<f64> {
    let mut result = DMatrix::zeros(matrix.nrows(), matrix.ncols());
    for i in 0..matrix.nrows() {
        for j in 0..matrix.ncols() {
            if j >= i + offset {
                result[(i, j)] = matrix[(i, j)];
            }
        }
    }
    result
}

fn evaluate(problem: &mut Problem, x: &[f64], y: &mut [f64]) {
    // let bounds = problem.get_ranges_of_interest();
    // let number_of_integer_variables = problem.number_of_integer_variables();
    // let mut result = vec![0.0; x.len()];
    // for (i, _) in x.iter().enumerate() {
    //     // let (lower, upper) = bounds[i].clone().into_inner();
    //     // result[i] = lower + x[i] * (upper - lower);

    //     if i < number_of_integer_variables {
    //         result[i] = result[i].round();
    //     }
    // }
    problem.evaluate_function(&x, y);
}
