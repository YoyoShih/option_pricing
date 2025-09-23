use std::f64;
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

    fn c_n_k(&self, n: u32, k: u32) -> u128 {
        let k = if k > n - k { n - k } else { k }; // to make k smaller than half of n by C(n, k) == C(n, n-k)
        let mut c: u128 = 1;
        for i in 0..k {
            c = c * (n - i) as u128 / (i + 1) as u128;
        }
        c
    }

    fn lower(&self, k: i32, c: i32) -> i32 {
        return (((k - c) / 2) as f64).floor() as i32;
    }

    fn upper(&self, k: i32, c: i32) -> i32 {
        return (((k - c) / 2) as f64).ceil() as i32;
    }

    fn summation<F>(&self, from: i32, to: i32, big_m: F, a: f64, b: f64, v: f64) -> f64
    where
        F: Fn(i32) -> i32,
    {
        // Calculation for V(c)
        let mut v_k: f64 = v * self.p_m.powi(self.n as i32 - from as i32) * (self.c_n_k(self.n, from as u32) as f64);
        let mut omega_k: f64 = 0.0;
        let mut delta_z: f64 = b.powi(from as i32);
        for z in 0..=(big_m(from)) {
            omega_k += delta_z;
            delta_z *= ((from - z) as f64 / (z + 1) as f64) * (a / b);
        }
        omega_k *= v_k;
        // Initialization of g_k
        let mut g_k: f64;
        // Using the equations to back-calculate g(k - 1) as initial value
        if big_m(from) == big_m(from + 1) {
            g_k = v_k * a.powi(big_m(from)) * b.powi(from - big_m(from));
            // Avoid negative k or M(k) in c(k, M(k))
            if from > 0 && big_m(from) > 0 {
                g_k *= self.c_n_k((from - 1) as u32, big_m(from) as u32) as f64;
            }
        } else {
            g_k = v_k * a.powi(big_m(from) + 1) * b.powi(from - big_m(from) - 1);
            // Avoid negative k or M(k) in c(k, M(k))
            if from > 0 && big_m(from) > 0 {
                g_k *= self.c_n_k((from - 1) as u32, big_m(from) as u32) as f64;
            }
        }
        // Summation loop for price
        let mut price: f64 = 0.0;
        // k: number of Up-Down moves
        for k in from..=to {
            // Calculate the payoff at maturity
            price += omega_k;
            // Update for next i
            let past_v_k: f64 = v_k;
            // FIX: Use floating-point division for all calculations
            v_k *= ((self.n as i32 - k) as f64) / (((k + 1) as f64) * self.p_m);
            let f_k: f64 = (b + a) * ((self.n as i32 - k) as f64 / ((k + 1) as f64 * self.p_m));
            if big_m(k) == big_m(k + 1) {
                g_k *= -(v_k / past_v_k) * a;
                g_k *= if k == big_m(k) { 1.0 } else { (k as f64) / ((k - big_m(k)) as f64) };
            } else {
                g_k *= -(v_k / past_v_k) * b;
                g_k *= if k == big_m(k) { 1.0 } else { (k as f64) / (big_m(k) as f64 + 1.0) };
            }
            omega_k = f_k * omega_k + g_k;
        }
        price
    }

    fn price_combinatorial(&self) -> f64 {
        let mut price: f64 = 0.0;
        if self.tree.style == OptionStyle::European {
            if self.tree.s <= self.tree.x {
                // Calculate the critical lower bound for positive payoff
                // Here the setting is for European Call option for a moment
                let c_l: i32 = ((self.tree.x / self.tree.s).ln() / self.u.ln()).ceil() as i32;
                let c_u: i32 = self.n as i32;
                let a: f64 = self.d * self.p_d;
                let b: f64 = self.u * self.p_u;
                price += self.summation(c_l, c_u - 1, |k| self.lower(k, c_l), a, b, self.tree.s) - self.summation(c_l, c_u - 1, |k| self.lower(k, c_l), self.p_d, self.p_u, self.tree.x);
                price 
                    += self.summation(c_u, self.n as i32, |k| self.lower(k, c_l), a, b, self.tree.s)
                    - self.summation(c_u, self.n as i32, |k| self.upper(k, c_u) - 1, a, b, self.tree.s)
                    - self.summation(c_u, self.n as i32, |k| self.lower(k, c_l), self.p_d, self.p_u, self.tree.x)
                    + self.summation(c_u, self.n as i32, |k| self.upper(k, c_u) - 1, self.p_d, self.p_u, self.tree.x);
            } else {}
            price * (-self.tree.r * self.tree.t).exp()
        } else {
            panic!("Combinatorial method is not applicable for American options");
        }
    }
}