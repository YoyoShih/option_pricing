use std::rc::Rc;

// Declare the crr module
mod tree;
// Import the CRR struct from the crr module
use tree::Tree;

fn main() {
    let s: f64 = 5.0;
    let x: f64 = 5.0;
    let t: f64 = 1.0;
    let r: f64 = 0.1;
    let sigma: f64 = 0.25;
    let q: f64 = 0.0;

    let tree = Rc::new(Tree::new(s, x, t, r, sigma, q));
    println!("Tree created with initial stock price: {}", tree.s);

    fn payoff_function1(price: f64, strike: f64) -> f64 { (price - 4.0) * (price - 5.0) * (price - 6.0) * (price - 7.0) + 5.0 - strike }
    let crr1 = tree::CRR::new(Rc::clone(&tree), 10000, Box::new(payoff_function1));
    println!("CRR created with {} time steps", crr1.n);
    let price = crr1.price();
    println!("Option price: {}", price);

    let krl1 = tree::KRL::new(Rc::clone(&tree), 10000, Box::new(payoff_function1));
    println!("KRL created with {} time steps", krl1.n);
    let price = krl1.price();
    println!("Option price: {}", price);
}
