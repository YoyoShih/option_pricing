use std::rc::Rc;

mod tree;
use tree::Tree;
use tree::OptionStyle;

// Declare the crr & krl module
mod crr;
mod krl;
// Import the CRR & KRL struct from the crr & krl module
use crr::CRR;
use krl::KRL;

fn main() {
    let s: f64 = 5.0;
    let x: f64 = 5.0;
    let t: f64 = 1.0;
    let r: f64 = 0.1;
    let sigma: f64 = 0.25;
    let q: f64 = 0.0;

    // Pure payoff functions for Experiments
    // Note that "the comparison to 0" is done in the pricing functions, not here in these functions
    // fn payoff_function1(price: f64, strike: f64) -> f64 { price - strike }
    fn payoff_function2(price: f64, strike: f64) -> f64 { (price - 4.0) * (price - 5.0) * (price - 6.0) * (price - 7.0) + 5.0 - strike }

    // European Option
    let euro_tree = Rc::new(Tree::new(s, x, t, r, sigma, q, OptionStyle::European));
    println!("Tree created with initial stock price: {}", euro_tree.s);

    // CRR & KRL with 10000 time steps
    let crr1 = CRR::new(Rc::clone(&euro_tree), 10000, Box::new(payoff_function2));
    println!("CRR created with {} time steps", crr1.n);
    let price = crr1.price();
    println!("Option price: {}", price);
    let krl1 = KRL::new(Rc::clone(&euro_tree), 10000, Box::new(payoff_function2));
    println!("KRL created with {} time steps", krl1.n);
    let price = krl1.price();
    println!("Option price: {}", price);

    // American Option
    let amer_tree = Rc::new(Tree::new(s, x, t, r, sigma, q, OptionStyle::American));
    println!("Tree created with initial stock price: {}", amer_tree.s);

    // CRR & KRL with 10000 time steps
    let crr1 = CRR::new(Rc::clone(&amer_tree), 10000, Box::new(payoff_function2));
    println!("CRR created with {} time steps", crr1.n);
    let price = crr1.price();
    println!("Option price: {}", price);
    let krl1 = KRL::new(Rc::clone(&amer_tree), 10000, Box::new(payoff_function2));
    println!("KRL created with {} time steps", krl1.n);
    let price = krl1.price();
    println!("Option price: {}", price);
}
