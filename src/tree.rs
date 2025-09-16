use std::f64;
use std::rc::Rc;

// Define the Tree struct
pub struct Tree {
    pub s: f64,
    pub x: f64,
    pub t: f64,
    pub r: f64,
    pub sigma: f64,
    pub q: f64,
}

// Implement methods for Tree
impl Tree {
    pub fn new(s: f64, x: f64, t: f64, r: f64, sigma: f64, q: f64) -> Self {
        Tree { s, x, t, r, sigma, q }
    }
}

// CRR: Cox-Ross-Rubinstein model for option pricing
// Define the CRR struct
pub struct CRR {
    pub tree: Rc<Tree>,
    pub u: f64,
    pub d: f64,
    pub p_u: f64,
    pub p_d: f64,
    pub delta_t: f64,
    pub n: u32,
    pub payoff: Box<dyn Fn(f64, f64) -> f64>,
}

// Implement methods for CRR
impl CRR {
    pub fn new(tree: Rc<Tree>, n: u32, payoff: Box<dyn Fn(f64, f64) -> f64>) -> Self {
        let delta_t: f64 = tree.t / n as f64;
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
            // curr = self.tree.s * self.u.powi((i - 1) as i32);
            for j in 0..i {
                payoffs[j as usize] = (-self.tree.r * self.delta_t).exp() * (self.p_u * payoffs[j as usize] + self.p_d * payoffs[(j + 1) as usize]);
            }
        }
        payoffs[0]
    }
}

// KRL: Karatzas-Ruf-Laplace model for option pricing
// Define the KRL struct
pub struct KRL {
    pub tree: Rc<Tree>,
    pub u: f64,
    pub d: f64,
    pub p_u: f64,
    pub p_d: f64,
    pub p_m: f64,
    pub delta_t: f64,
    pub n: u32,
    pub payoff: Box<dyn Fn(f64, f64) -> f64>,
}

// Implement methods for KRL
impl KRL {
    pub fn new(tree: Rc<Tree>, n: u32, payoff: Box<dyn Fn(f64, f64) -> f64>) -> Self {
        let delta_t: f64 = tree.t / n as f64;
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
            for i in 0..=2*j {
                payoffs[i as usize] = (-self.tree.r * self.delta_t).exp() * (self.p_u * payoffs[i as usize] + self.p_m * payoffs[(i + 1) as usize] + self.p_d * payoffs[(i + 2) as usize]);
            }
        }
        payoffs[0]
    }
}