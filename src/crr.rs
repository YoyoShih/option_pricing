use std::f64;
use std::rc::Rc;

use crate::tree::Tree;
use tree::OptionStyle;

// CRR: Cox-Ross-Rubinstein model for option pricing
// Define the CRR struct
pub struct CRR {
    tree: Rc<Tree>, // reference to the Tree struct
    u: f64, // up factor
    d: f64, // down factor
    p_u: f64, // probability of up move
    p_d: f64, // probability of down move
    delta_t: f64, // time step
    pub n: u32, // number of time steps
    payoff: Box<dyn Fn(f64, f64) -> f64>, // payoff function
}

// Implement methods for CRR
impl CRR {
    pub fn new(tree: Rc<Tree>, n: u32, payoff: Box<dyn Fn(f64, f64) -> f64>) -> Self {
        let delta_t: f64 = tree.t / n as f64; // time step
        // up and down factors and risk-neutral probabilities
        let u: f64 = (tree.sigma * (delta_t).sqrt()).exp();
        let d: f64 = (-tree.sigma * (delta_t).sqrt()).exp();
        let p_u: f64 = (((tree.r - tree.q) * delta_t).exp() - d) / (u - d);
        let p_d: f64 = 1.0 - p_u;

        CRR {
            tree,
            u,
            d,
            p_u,
            p_d,
            delta_t,
            n,
            payoff,
        }
    }
    pub fn price(&self) -> f64 {
        // Calculate the payoff at maturity
        let mut curr: f64 = self.tree.s * self.u.powi(self.n as i32);
        let mut payoffs: Vec<f64> = vec![0.0; (self.n + 1) as usize];
        for i in 0..=self.n {
            payoffs[i as usize] = f64::max((self.payoff)(curr, self.tree.x), 0.0);
            curr *= self.d * self.d;
        }
        // Backward induction to get the option price at time 0
        for i in (1..=self.n).rev() {
            curr = self.tree.s * self.u.powi((i - 1) as i32);
            for j in 0..i {
                // European option
                if self.tree.style == OptionStyle::European {
                    payoffs[j as usize] = (-self.tree.r * self.delta_t).exp() * (self.p_u * payoffs[j as usize] + self.p_d * payoffs[(j + 1) as usize]);
                } else if self.tree.style == OptionStyle::American {
                    // American option
                    let exercise_value = f64::max((self.payoff)(curr, self.tree.x), 0.0);
                    let hold_value = (-self.tree.r * self.delta_t).exp() * (self.p_u * payoffs[j as usize] + self.p_d * payoffs[(j + 1) as usize]);
                    payoffs[j as usize] = f64::max(exercise_value, hold_value);
                    curr *= self.d * self.d;
                }
            }
        }
        payoffs[0]
    }
}