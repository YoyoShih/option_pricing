use std::f64;
use std::rc::Rc;

use crate::tree::Tree;
use tree::OptionStyle;

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
    pub fn price(&self) -> f64 {
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
}