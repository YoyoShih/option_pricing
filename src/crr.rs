use std::f64;
use std::rc::Rc;

use crate::tree::Tree;
use tree::OptionStyle;
use tree::CalcMethod;

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
    pub fn price(&self, calc_method: &CalcMethod) -> f64 {
        let payoff = match calc_method {
            CalcMethod::BackwardInduction => self.price_backward_induction(),
            CalcMethod::Combinatorial => self.price_combinatorial(),
            CalcMethod::FiniteDifference { kind, n_s, n_t, s_max, s_min } => panic!("Finite Difference method is not implemented for CRR"),
        };
        payoff
    }
    fn price_backward_induction(&self) -> f64 {
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
                if self.tree.option_spec.style == OptionStyle::European {
                    payoffs[j as usize] = (-self.tree.r * self.delta_t).exp() * (self.p_u * payoffs[j as usize] + self.p_d * payoffs[(j + 1) as usize]);
                } else if self.tree.option_spec.style == OptionStyle::American {
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

    fn price_combinatorial(&self) -> f64 {
        if self.tree.option_spec.style == OptionStyle::European {
            let mut price: f64 = 0.0;
            let mut s_t: f64 = self.tree.s * self.u.powi(self.n as i32); // stock price at node (n, n)
            // Note that here we take log to comb and prob to avoid underflow when n is large
            let mut log_comb: f64 = 0.0; // log of combination nCi
            let mut log_prob: f64 = (self.n as f64) * self.p_u.ln(); // log of probability p_u^i * p_d^(n-i)

            for i in 0..=self.n {
                // Calculate the payoff at maturity
                let payoff: f64 = f64::max((self.payoff)(s_t, self.tree.x), 0.0);
                price += (log_comb + log_prob).exp() * payoff;
                // Update for next i
                s_t *= self.d * self.d; // s_t = S_0 * u^i * d^(n-i), when i decreases by 1, s_t *= d^2
                log_comb += ((self.n - i) as f64).ln() - ((i + 1) as f64).ln(); // nCi = nC(i-1) * (n-i+1)/i
                log_prob += self.p_d.ln() - self.p_u.ln(); // p_u^i * p_d^(n-i) = p_u^(i-1) * p_d^(n-i+1) * (p_d/p_u)
            }
            price * (-self.tree.r * self.tree.t).exp()
        } else {
            panic!("Combinatorial method is not applicable for American options");
        }
    }
}