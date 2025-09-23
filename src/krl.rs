use std::f64;
use std::cmp;
use std::rc::Rc;

use crate::tree::Tree;
use tree::OptionStyle;
use tree::CalcMethod;

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
                if self.tree.style == OptionStyle::European {
                    payoffs[i as usize] = (-self.tree.r * self.delta_t).exp() * (self.p_u * payoffs[i as usize] + self.p_m * payoffs[(i + 1) as usize] + self.p_d * payoffs[(i + 2) as usize]);
                } else if self.tree.style == OptionStyle::American {
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

    fn lower(&self, k: i32, c: i32) -> i32 {
        return ((k - c) as f64 / 2.0).floor() as i32;
    }

    fn upper(&self, k: i32, c: i32) -> i32 {
        return ((k - c) as f64 / 2.0).ceil() as i32;
    }

    // Helper for log-sum-exp trick
    fn log_sum_exp(&self, a: f64, b: f64, sign_a: i8, sign_b: i8) -> f64 {
        if a > b { a + (1.0 + (sign_b as f64) * (b - a).exp()).ln() }
        else { b + (1.0 + (sign_a as f64) * (a - b).exp()).ln() }
    }

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
        let mut log_omega_k: f64 = f64::NEG_INFINITY; // initialization of omega_k in log scale as log(0)
        let mut log_delta_z: f64 = log_b * (from as f64);
        for z in 0..=(big_m(from)) {
            log_omega_k = self.log_sum_exp(log_omega_k, log_delta_z, 1, 1);
            log_delta_z += ((from - z) as f64).ln() - (1.0 + (z as f64)).ln() + log_a - log_b;
        }
        log_omega_k += log_v_k;
        // Initialization of g_k
        let mut log_g_k: f64;
        if big_m(from) == big_m(from + 1) {
            log_g_k = log_v_k + (big_m(from) as f64 * log_a + (from - big_m(from)) as f64 * log_b);
        } else {
            log_g_k = log_v_k + ((big_m(from) + 1) as f64 * log_a + (from - big_m(from) - 1) as f64 * log_b);
        }
        log_g_k += self.log_c_n_k((from - 1) as u32, big_m(from) as u32) as f64;
        // If big_m(from) == big_m(from + 1), sign should be negative
        let mut sign_g_k: i8 = (-1f64).powi(big_m(from + 1) - big_m(from) - 1) as i8;
        // Summation loop for price
        let mut log_price: f64 = f64::NEG_INFINITY; // initialization of price in log scale as log(0)
        // k: number of Up-Down moves
        for k in from..=to {
            // Calculate the payoff at maturity
            log_price = self.log_sum_exp(log_price, log_omega_k, 1, 1);
            // Update for next i
            let past_log_v_k: f64 = log_v_k;
            log_v_k += ((self.n as i32 - k) as f64).ln() - (((k + 1) as f64) * self.p_m).ln();
            let log_f_k: f64 = (b + a).ln() + ((self.n as i32 - k) as f64).ln() - ((k + 1) as f64).ln() - self.p_m.ln();
            if big_m(k) == big_m(k + 1) {
                log_g_k += log_v_k - past_log_v_k + log_a;
                log_g_k += if k == big_m(k) { 0.0 } else { (k as f64).ln() - ((k - big_m(k)) as f64).ln() };
            } else {
                log_g_k += log_v_k - past_log_v_k + log_b;
                log_g_k += if k == big_m(k) { 0.0 } else { (k as f64).ln() - (big_m(k) as f64 + 1.0).ln() };
            }
            log_omega_k = self.log_sum_exp(log_f_k + log_omega_k, log_g_k, 1, sign_g_k);
            // Update for next k
            sign_g_k *= -1;
        }
        log_price.exp()
    }

    fn price_combinatorial(&self) -> f64 {
        let mut price: f64 = 0.0;
        let a: f64 = self.d * self.p_d;
        let b: f64 = self.u * self.p_u;
        // Calculate the critical bound for positive payoff
        // Here the setting is for European Call option for a moment
        let c_l: i32 = cmp::max(((self.tree.x / self.tree.s).ln() / self.u.ln()).ceil() as i32, -(self.n as i32)) as i32;
        let c_u: i32 = self.n as i32;
        if self.tree.style == OptionStyle::European {
            if self.tree.s <= self.tree.x {
                price += self.summation(c_l, c_u, |k| self.lower(k, c_l), a, b, self.tree.s) - self.summation(c_l, c_u, |k| self.lower(k, c_l), self.p_d, self.p_u, self.tree.x);
            } else {
                // Decompose the summation into two parts: c_l to -1 and 0 to c_u
                // 0 to c_u
                price += self.summation(0, c_u, |k| self.lower(k, 0), a, b, self.tree.s) - self.summation(0, c_u, |k| self.lower(k, 0), self.p_d, self.p_u, self.tree.x);
                // c_l to -1
                price += self.summation(-c_l, self.n as i32, |k| self.lower(k, c_l), a, b, self.tree.s) - self.summation(-c_l, self.n as i32, |k| self.lower(k, c_l), self.p_d, self.p_u, self.tree.x);
                price -= self.summation(-c_l, self.n as i32, |k| self.upper(k, -1) - 1, a, b, self.tree.s) - self.summation(-c_l, self.n as i32, |k| self.upper(k, -1) - 1, self.p_d, self.p_u, self.tree.x);
                price += self.summation(1, -c_l - 1, |k| k, a, b, self.tree.s) - self.summation(1, -c_l - 1, |k| k, self.p_d, self.p_u, self.tree.x);
                price -= self.summation(1, -c_l - 1, |k| self.upper(k, -1) - 1, a, b, self.tree.s) - self.summation(1, -c_l - 1, |k| self.upper(k, -1) - 1, self.p_d, self.p_u, self.tree.x);
            }
            price * (-self.tree.r * self.tree.t).exp()
        } else {
            panic!("Combinatorial method is not applicable for American options");
        }
    }
}