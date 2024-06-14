pub struct Problem {
    pub dimension: usize,
    pub functions: Vec<fn(f64) -> f64>,
    pub best_known_fitness: f64,
    pub best_known_solution: Vec<f64>,
}

impl Problem {
    pub fn new(functions: Vec<fn(f64) -> f64>) -> Self {
        Self {
            dimension: functions.len(),
            functions,
            best_known_fitness: f64::INFINITY,
            best_known_solution: vec![],
        }
    }

    pub fn initial_solution(&self, x: &mut [f64]) {
        for i in 0..self.dimension {
            x[i] = 5.0;
        }
    }

    pub fn evaluate_function(&mut self, x: &[f64], y: &mut [f64]) {
        let mut sum = 0.0;
        for i in 0..self.dimension {
            sum += self.functions[i](x[i]);
        }
        y[0] = sum;
        if sum < self.best_known_fitness {
            self.best_known_fitness = sum;
            self.best_known_solution = x.to_vec();
        }
    }

    // pub fn get_ranges_of_interest(&self) -> Vec<RangeInclusive<f64>> {
    //     let mut ranges = Vec::new();
    //     for i in 0..self.dimension {
    //         ranges.push(-10000.0..=10000.0);
    //     }
    //     ranges
    // }

    pub fn number_of_integer_variables(&self) -> usize {
        0
    }

    pub fn dimension(&self) -> usize {
        self.dimension
    }

    pub fn number_of_objectives(&self) -> usize {
        1
    }
}
