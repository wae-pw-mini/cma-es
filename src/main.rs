mod algorithm;
mod problem;

use algorithm::cma_es;
use problem::Problem;

fn main() {
    let mut problem = Problem::new(vec![|x1| (x1 - 100.0) * (x1 - 100.0), |x2| x2 * x2]);
    let x = &vec![1.0, 2.0];
    let y = &mut vec![0.0; 2];
    problem.evaluate_function(x, y);
    let max_budget = 10000;
    cma_es(&mut problem, max_budget);
    println!("The best known fitness is {}", problem.best_known_fitness);
    println!(
        "The best known solution is {:?}",
        problem.best_known_solution
    );
}
