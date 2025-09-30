use std::f64;
use std::cmp;
use std::rc::Rc;
use std::vec;

use nalgebra as na;
use na::{DMatrix};

use crate::utils::SignedLog;

use crate::tree::{Tree, OptionStyle, OptionType, CalcMethod, FiniteDifferenceType};

// KRL: Karatzas-Ruf-Laplace model for option pricing
// Define the KRL struct
pub struct KRL {
    tree: Rc<Tree>, // reference to the Tree struct
    u: f64, // up factor
    d: f64, // down factor
    p_u: f64, // probability of up move
    p_d: f64, // probability of down move
    p_m: f64, // probability of middle move
    delta_t: f64, // time step
    pub n: u32, // number of time steps
    payoff: Box<dyn Fn(f64, f64) -> f64>, // payoff function
}

// Implement methods for KRL
impl KRL {
    pub fn new(tree: Rc<Tree>, n: u32, payoff: Box<dyn Fn(f64, f64) -> f64>) -> Self {
        // time step
        let delta_t: f64 = tree.t / n as f64;
        // up and down factors and risk-neutral probabilities
        let lambda: f64 = (1.5f64).sqrt();
        let u: f64 = (lambda * tree.sigma * (delta_t).sqrt()).exp();
        let d: f64 = (-lambda * tree.sigma * (delta_t).sqrt()).exp();
        let p_u: f64 = 0.5 / (lambda * lambda) + 0.5 * (tree.r - tree.q - 0.5 * tree.sigma * tree.sigma) * (delta_t).sqrt() / (lambda * tree.sigma);
        let p_d: f64 = 0.5 / (lambda * lambda) - 0.5 * (tree.r - tree.q - 0.5 * tree.sigma * tree.sigma) * (delta_t).sqrt() / (lambda * tree.sigma);
        let p_m: f64 = 1.0 - p_u - p_d;

        KRL {
            tree,
            u,
            d,
            p_u,
            p_d,
            p_m,
            delta_t,
            n,
            payoff,
        }
    }
    pub fn price(&self, calc_method: &CalcMethod) -> f64 {
        let payoff = match calc_method {
            CalcMethod::BackwardInduction => self.price_backward_induction(),
            CalcMethod::Combinatorial => self.price_combinatorial(),
            CalcMethod::FiniteDifference{ kind, n_s, n_t, s_max, s_min } => {
                let delta_s: f64 = (*s_max - *s_min) / *n_s as f64;
                let delta_t: f64 = self.tree.t / *n_t as f64;
                let sigma_sqr: f64 = self.tree.sigma.powi(2);
                let s0_index: i32 = ((self.tree.s - *s_min) / delta_s).floor() as i32;
                if delta_s * s0_index as f64 + *s_min != self.tree.s {
                    panic!("Finite Difference method requires the initial stock price to be on the grid points");
                } else {
                    match kind {
                        FiniteDifferenceType::Explicit => self.price_finite_difference_explicit(*n_s, *n_t, *s_max, *s_min, delta_s, delta_t, sigma_sqr, s0_index),
                        FiniteDifferenceType::Implicit => self.price_finite_difference_implicit(*n_s, *n_t, *s_max, *s_min, delta_s, delta_t, sigma_sqr, s0_index),
                    }
                }
            },
        };
        payoff
    }
    fn price_backward_induction(&self) -> f64 {
        // Calculate the payoff at maturity
        let mut curr: f64 = self.tree.s * self.u.powi(self.n as i32);
        let mut payoffs: Vec<f64> = vec![0.0; (2*self.n + 1) as usize];
        for i in 0..=2*self.n {
            payoffs[i as usize] = f64::max((self.payoff)(curr, self.tree.x), 0.0);
            curr *= self.d;
        }
        // Backward induction to get the option price at time 0
        for j in (0..=self.n-1).rev() {
            curr = self.tree.s * self.u.powi(j as i32);
            for i in 0..=2*j {
                // European option
                if self.tree.option_spec.style == OptionStyle::European {
                    payoffs[i as usize] = (-self.tree.r * self.delta_t).exp() * (self.p_u * payoffs[i as usize] + self.p_m * payoffs[(i + 1) as usize] + self.p_d * payoffs[(i + 2) as usize]);
                } else if self.tree.option_spec.style == OptionStyle::American {
                    // American option
                    let exercise_value = f64::max((self.payoff)(curr, self.tree.x), 0.0);
                    let hold_value = (-self.tree.r * self.delta_t).exp() * (self.p_u * payoffs[i as usize] + self.p_m * payoffs[(i + 1) as usize] + self.p_d * payoffs[(i + 2) as usize]);
                    payoffs[i as usize] = f64::max(exercise_value, hold_value);
                }
                curr *= self.d;
            }
        }
        payoffs[0]
    }

    fn log_c_n_k(&self, n: u32, k: u32) -> f64 {
        if n < k { return 0.0; } // for initialization for invalid cases
        let k = if k > n - k { n - k } else { k }; // to make k smaller than half of n by C(n, k) == C(n, n-k)
        let mut log_c: f64 = 0.0;
        for i in 0..k { log_c += ((n - i) as f64).ln() - ((i + 1) as f64).ln(); }
        log_c
    }

    fn lower(&self, k: i32, c: i32) -> i32 { ((k - c) as f64 / 2.0).floor() as i32 }

    fn upper(&self, k: i32, c: i32) -> i32 { ((k - c) as f64 / 2.0).ceil() as i32 }

    // Helper function for combinatorial pricing
    fn summation<F>(&self, from: i32, to: i32, big_m: F, a: f64, b: f64, v: f64) -> f64
    where
        F: Fn(i32) -> i32,
    {
        if from > to { return 0.0; }
        // Precompute logs
        let log_v = v.ln();
        let log_p_m = self.p_m.ln();
        let log_b = b.ln();
        let log_a = a.ln();

        // Calculation for V(c)
        let mut log_v_k: f64 = log_v + log_p_m * (self.n as f64 - from as f64) + self.log_c_n_k(self.n, from as u32);
        let mut omega_k = SignedLog { sign: 0, log_value: f64::NEG_INFINITY };
        let mut delta_z = SignedLog { sign: 1, log_value: log_b * (from as f64) };
        for z in 0..=(big_m(from)) {
            omega_k = SignedLog::log_sum(&omega_k, &delta_z);
            let next_log = ((from - z) as f64).ln() - (1.0 + (z as f64)).ln() + log_a - log_b;
            delta_z = SignedLog { sign: delta_z.sign, log_value: delta_z.log_value + next_log };
        }
        omega_k.log_value += log_v_k;
        // Initialization of g_k
        let mut g_k: SignedLog;
        if big_m(from) == big_m(from + 1) {
            g_k = SignedLog { sign: 1, log_value: log_v_k + (big_m(from) as f64 * log_a + (from - big_m(from)) as f64 * log_b) };
        } else {
            g_k = SignedLog { sign: 1, log_value: log_v_k + ((big_m(from) + 1) as f64 * log_a + (from - big_m(from) - 1) as f64 * log_b) };
        }
        g_k.log_value += self.log_c_n_k((from - 1) as u32, big_m(from) as u32);
        // If big_m(from) == big_m(from + 1), sign should be negative
        g_k.sign = (-1f64).powi(big_m(from + 1) - big_m(from) - 1) as i8;
        // Summation loop for price
        let mut price = SignedLog { sign: 0, log_value: f64::NEG_INFINITY };
        // k: number of Up-Down moves
        for k in from..=to {
            // Calculate the payoff at maturity
            price = SignedLog::log_sum(&price, &omega_k);
            // Update for next i
            let past_log_v_k = log_v_k;
            log_v_k += ((self.n as i32 - k) as f64).ln() - (((k + 1) as f64) * self.p_m).ln();
            let log_f_k = (b + a).ln() + ((self.n as i32 - k) as f64).ln() - ((k + 1) as f64).ln() - self.p_m.ln();
            if big_m(k) == big_m(k + 1) {
                g_k.log_value += log_v_k - past_log_v_k + log_a;
                // g_k = 0 when big_m(k) == k
                g_k.log_value += if k == big_m(k) { 0.0 } else { (k as f64).ln() - ((k - big_m(k)) as f64).ln() };
            } else {
                g_k.log_value += log_v_k - past_log_v_k + log_b;
                // g_k = 0 when big_m(k) == -1
                g_k.log_value += if big_m(k) == -1 { 0.0 } else { (k as f64).ln() - (big_m(k) as f64 + 1.0).ln() };
            }
            if big_m(k) == k && big_m(k + 1) == k + 1 {
                omega_k.log_value += log_f_k;
            } else {
                omega_k = SignedLog::log_sum(&SignedLog { sign: 1, log_value: log_f_k + omega_k.log_value }, &g_k);
            }
            // Update for next k
            g_k.sign *= -1;
        }
        price.sign as f64 * price.log_value.exp()
    }

    // Combinatorial pricing method
    fn price_combinatorial(&self) -> f64 {
        let mut price: f64 = 0.0;
        let a: f64 = self.d * self.p_d;
        let b: f64 = self.u * self.p_u;
        // Calculate the critical bound for positive payoff
        match self.tree.option_spec.kind {
            OptionType::Call => {
                // Here the setting is for European Call option for a moment
                let c_l: i32 = cmp::max(((self.tree.x / self.tree.s).ln() / self.u.ln()).ceil() as i32, -(self.n as i32)) as i32;
                let c_u: i32 = self.n as i32;
                match self.tree.option_spec.style {
                    OptionStyle::European => {
                        if self.tree.s <= self.tree.x {
                            price += self.summation(c_l, self.n as i32, |k| self.lower(k, c_l), a, b, self.tree.s) - self.summation(c_l, self.n as i32, |k| self.lower(k, c_l), self.p_d, self.p_u, self.tree.x);
                        } else {
                            // Decompose the summation into two parts: c_l to -1 and 0 to c_u = n
                            // 0 to c_u
                            price += self.summation(0, c_u, |k| self.lower(k, 0), a, b, self.tree.s) - self.summation(0, c_u, |k| self.lower(k, 0), self.p_d, self.p_u, self.tree.x);
                            // c_l to -1
                            price += self.summation(-c_l, self.n as i32, |k| self.lower(k, c_l), a, b, self.tree.s) - self.summation(-c_l, self.n as i32, |k| self.lower(k, c_l), self.p_d, self.p_u, self.tree.x);
                            price -= self.summation(-c_l, self.n as i32, |k| self.upper(k, -1) - 1, a, b, self.tree.s) - self.summation(-c_l, self.n as i32, |k| self.upper(k, -1) - 1, self.p_d, self.p_u, self.tree.x);
                            price += self.summation(1, -c_l - 1, |k| k, a, b, self.tree.s) - self.summation(1, -c_l - 1, |k| k, self.p_d, self.p_u, self.tree.x);
                            price -= self.summation(1, -c_l - 1, |k| self.upper(k, -1) - 1, a, b, self.tree.s) - self.summation(1, -c_l - 1, |k| self.upper(k, -1) - 1, self.p_d, self.p_u, self.tree.x);
                        }
                        price * (-self.tree.r * self.tree.t).exp()
                    }
                    OptionStyle::American => {
                        panic!("Combinatorial method is not applicable for American options");
                    }
                }
            }
        OptionType::Put => {
                // Here the setting is for European Put option for a moment
                let c_l: i32 = -(self.n as i32);
                let c_u: i32 = cmp::max(((self.tree.x / self.tree.s).ln() / self.u.ln()).floor() as i32, -(self.n as i32)) as i32;
                match self.tree.option_spec.style {
                    OptionStyle::European => {
                        if self.tree.s >= self.tree.x {
                            price += self.summation(-c_u, -c_l, |k| k, self.p_d, self.p_u, self.tree.x) - self.summation(-c_u, -c_l, |k| k, a, b, self.tree.s);
                            price -= self.summation(-c_u, -c_l, |k| self.upper(k, c_u) - 1, self.p_d, self.p_u, self.tree.x) - self.summation(-c_u, -c_l, |k| self.upper(k, c_u) - 1, a, b, self.tree.s);
                        } else {
                            // Decompose the summation into two parts: -n = c_l to -1 and 0 to c_u
                            // 0 to c_u
                            price += self.summation(c_u, self.n as i32, |k| self.lower(k, 0), self.p_d, self.p_u, self.tree.x) - self.summation(c_u, self.n as i32, |k| self.lower(k, 0), a, b, self.tree.s);
                            price -= self.summation(c_u, self.n as i32, |k| self.upper(k, c_u) - 1, self.p_d, self.p_u, self.tree.x) - self.summation(c_u, self.n as i32, |k| self.upper(k, c_u) - 1, a, b, self.tree.s);
                            price += self.summation(0, c_u - 1, |k| self.lower(k, 0), self.p_d, self.p_u, self.tree.x) - self.summation(0, c_u - 1, |k| self.lower(k, 0), a, b, self.tree.s);
                            // c_l to -1
                            price += self.summation(1, -c_l, |k| k, self.p_d, self.p_u, self.tree.x) - self.summation(1, -c_l, |k| k, a, b, self.tree.s);
                            price -= self.summation(1, -c_l, |k| self.upper(k, -1) - 1, self.p_d, self.p_u, self.tree.x) - self.summation(1, -c_l, |k| self.upper(k, -1) - 1, a, b, self.tree.s);
                        }
                        price * (-self.tree.r * self.tree.t).exp()
                    }
                    OptionStyle::American => {
                        panic!("Combinatorial method is not applicable for American options");
                    }
                }
            }
        }
    }

    // Finite Difference pricing method
    // Explicit method
    fn price_finite_difference_explicit(&self, n_s: u32, n_t: u32, s_max: f64, s_min: f64, delta_s: f64, delta_t: f64, sigma_sqr: f64, s0_index: i32) -> f64 {
        let r1: f64 = 1.0 / (1.0 + self.tree.r * delta_t);
        let r2: f64 = delta_t / (1.0 + self.tree.r * delta_t);
        let mut grid: Vec<Vec<f64>> = vec![vec![0.0; n_s as usize]; 3];
        for i in 0..n_s {
            grid[0][i as usize] = r2 * 0.5 * i as f64 * (-(self.tree.r - self.tree.q) + sigma_sqr * i as f64);
            grid[1][i as usize] = r1 * (1.0 - sigma_sqr * i.pow(2) as f64 * delta_t);
            grid[2][i as usize] = r2 * 0.5 * i as f64 * (self.tree.r - self.tree.q + sigma_sqr * i as f64);
        }
        let mut value_next: Vec<f64> = vec![0.0; (n_s + 1) as usize];
        let mut value_now: Vec<f64> = vec![0.0; (n_s + 1) as usize];
        let mut s_t: Vec<f64> = vec![0.0; (n_s + 1) as usize];
        for i in 0..=n_s {
            s_t[i as usize] = delta_s * i as f64 + s_min;
            value_next[i as usize] = f64::max((self.payoff)(s_t[i as usize], self.tree.x), 0.0);
        }
        for _ in 0..n_t {
            value_now[0] = f64::max((self.payoff)(s_min, self.tree.x), 0.0);
            for i in 1..n_s {
                value_now[i as usize] = grid[0][i as usize] * value_next[i as usize - 1] + grid[1][i as usize] * value_next[i as usize] + grid[2][i as usize] * value_next[i as usize + 1];
            }
            value_now[n_s as usize] = f64::max((self.payoff)(s_max, self.tree.x), 0.0);
            for i in 0..=n_s {
                match self.tree.option_spec.style {
                    OptionStyle::European => {
                        value_next[i as usize] = f64::max(value_now[i as usize], 0.0);
                    }
                    OptionStyle::American => {
                        value_next[i as usize] = f64::max(value_now[i as usize], (self.payoff)(s_t[i as usize], self.tree.x));
                    }
                }
            }
        }
        value_now[s0_index as usize]
    }

    // Implicit method
    fn price_finite_difference_implicit(&self, n_s: u32, n_t: u32, s_max: f64, s_min: f64, delta_s: f64, delta_t: f64, sigma_sqr: f64, s0_index: i32) -> f64 {
        match self.tree.option_spec.style {
            OptionStyle::European => {
                let mut matrix: DMatrix<f64> = DMatrix::zeros(n_s as usize + 1, n_s as usize + 1);
                matrix[(0, 0)] = 1.0 + delta_t * (sigma_sqr * (n_s as f64).powi(2) + self.tree.r); // B
                matrix[(0, 1)] = 0.5 * (n_s as f64) * delta_t * (self.tree.r - self.tree.q - sigma_sqr * (n_s as f64)); // A
                for i in 1..n_s {
                    matrix[(i as usize, i as usize - 1)] = 0.5 * ((n_s - i) as f64) * delta_t * (-(self.tree.r - self.tree.q) - sigma_sqr * ((n_s - i) as f64)); // C
                    matrix[(i as usize, i as usize)] = 1.0 + delta_t * (sigma_sqr * ((n_s - i) as f64).powi(2) + self.tree.r); // B
                    matrix[(i as usize, i as usize + 1)] = 0.5 * ((n_s - i) as f64) * delta_t * (self.tree.r - self.tree.q - sigma_sqr * ((n_s - i) as f64)); // A
                }
                matrix[(n_s as usize, n_s as usize - 1)] = 0.5 * delta_t * (-(self.tree.r - self.tree.q)); // C
                matrix[(n_s as usize, n_s as usize)] = 1.0 + delta_t * self.tree.r; // B
                let inv_matrix: DMatrix<f64> = matrix.try_inverse().expect("Matrix inversion failed in Implicit Finite Difference method");

                let mut value: DMatrix<f64> = DMatrix::zeros(n_s as usize + 1, 1);
                for i in 0..=n_s {
                    let s_t: f64 = delta_s * ((n_s - i) as f64) + s_min;
                    value[(i as usize, 0)] = f64::max((self.payoff)(s_t, self.tree.x), 0.0);
                }
                for _ in 0..n_t {
                    value[(0, 0)] = f64::max((self.payoff)(s_max, self.tree.x), 0.0);
                    value[(n_s as usize, 0)] = f64::max((self.payoff)(s_min, self.tree.x), 0.0);
                    value = &inv_matrix * &value;
                    for i in 1..n_s {
                        value[(i as usize, 0)] = f64::max(value[(i as usize, 0)], 0.0);
                    }
                }
                value[((n_s - s0_index as u32) as usize, 0)]
            },
            OptionStyle::American => {
                let mut matrix: DMatrix<f64> = DMatrix::zeros(n_s as usize - 1, n_s as usize - 1);
                matrix[(0, 0)] = 1.0 + delta_t * (sigma_sqr * (n_s as f64 - 1.0).powi(2) + self.tree.r); // B
                matrix[(0, 1)] = 0.5 * (n_s as f64 - 1.0) * delta_t * (self.tree.r - self.tree.q - sigma_sqr * (n_s as f64 - 1.0)); // A
                for i in 1..(n_s - 2) {
                    matrix[(i as usize, i as usize - 1)] = 0.5 * ((n_s - 1 - i) as f64) * delta_t * (-(self.tree.r - self.tree.q) - sigma_sqr * ((n_s - 1 - i) as f64)); // C
                    matrix[(i as usize, i as usize)] = 1.0 + delta_t * (sigma_sqr * ((n_s - 1 - i) as f64).powi(2) + self.tree.r); // B
                    matrix[(i as usize, i as usize + 1)] = 0.5 * ((n_s - 1 - i) as f64) * delta_t * (self.tree.r - self.tree.q - sigma_sqr * ((n_s - 1 - i) as f64)); // A
                }
                matrix[(n_s as usize - 2, n_s as usize- 3)] = 0.5 * delta_t * (-(self.tree.r - self.tree.q) - sigma_sqr); // C
                matrix[(n_s as usize - 2, n_s as usize- 2)] = 1.0 + delta_t * (sigma_sqr + self.tree.r); // B
                let inv_matrix: DMatrix<f64> = matrix.try_inverse().expect("Matrix inversion failed in Implicit Finite Difference method");

                let mut s_t: DMatrix<f64> = DMatrix::zeros(n_s as usize - 1, 1);
                let mut value: DMatrix<f64> = DMatrix::zeros(n_s as usize - 1, 1);
                for i in 0..(n_s - 1) {
                    s_t[(i as usize, 0)] = delta_s * ((n_s - 1 - i) as f64) + s_min;
                    value[(i as usize, 0)] = f64::max((self.payoff)(s_t[(i as usize, 0)], self.tree.x), 0.0);
                }
                for _ in 0..n_t {
                    value[(0, 0)] = f64::max((self.payoff)(s_max, self.tree.x), 0.0);
                    value = &inv_matrix * &value;
                    for i in 1..(n_s - 1) {
                        value[(i as usize, 0)] = f64::max(value[(i as usize, 0)], (self.payoff)(s_t[(i as usize, 0)], self.tree.x));
                        value[(i as usize, 0)] = f64::max(value[(i as usize, 0)], 0.0);
                    }
                }
                value[((n_s - 1 - s0_index as u32) as usize, 0)]
            }
        }
    }
}