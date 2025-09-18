use std::rc::Rc;

mod tree;
use tree::Tree;
use tree::OptionStyle;
use tree::CalcMethod;

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
    println!("Tree of European option created with initial stock price: {}", euro_tree.s);

    // CRR & KRL with 10000 time steps (Backward Induction)
    let calc_method_bi = CalcMethod::BackwardInduction;
    let euro_crr_bi = CRR::new(Rc::clone(&euro_tree), 10000, Box::new(payoff_function2));
    println!("CRR (Backward Induction) created with {} time steps", euro_crr_bi.n);
    let euro_crr_bi_price = euro_crr_bi.price(&calc_method_bi);
    println!("Option price: {}", euro_crr_bi_price);
    let euro_krl_bi = KRL::new(Rc::clone(&euro_tree), 10000, Box::new(payoff_function2));
    println!("KRL (Backward Induction) created with {} time steps", euro_krl_bi.n);
    let euro_krl_bi_price = euro_krl_bi.price(&calc_method_bi);
    println!("Option price: {}", euro_krl_bi_price);

    // CRR & KRL with 10000 time steps (Combinatorial)
    let calc_method_comb = CalcMethod::Combinatorial;
    let euro_crr_comb = CRR::new(Rc::clone(&euro_tree), 10000, Box::new(payoff_function2));
    println!("CRR (Combinatorial) created with {} time steps", euro_crr_comb.n);
    let euro_crr_comb_price = euro_crr_comb.price(&calc_method_comb);
    println!("Option price: {}", euro_crr_comb_price);
    let euro_krl_comb = KRL::new(Rc::clone(&euro_tree), 10000, Box::new(payoff_function2));
    println!("KRL (Combinatorial) created with {} time steps", euro_krl_comb.n);
    let euro_krl_comb_price = euro_krl_comb.price(&calc_method_comb);
    println!("Option price: {}", euro_krl_comb_price);

    // American Option
    let amer_tree = Rc::new(Tree::new(s, x, t, r, sigma, q, OptionStyle::American));
    println!("Tree of American option created with initial stock price: {}", amer_tree.s);

    // CRR & KRL with 10000 time steps (Backward Induction)
    let amer_crr_bi = CRR::new(Rc::clone(&amer_tree), 10000, Box::new(payoff_function2));
    println!("CRR (Backward Induction) created with {} time steps", amer_crr_bi.n);
    let amer_crr_bi_price = amer_crr_bi.price(&calc_method_bi);
    println!("Option price: {}", amer_crr_bi_price);
    let amer_krl_bi = KRL::new(Rc::clone(&amer_tree), 10000, Box::new(payoff_function2));
    println!("KRL (Backward Induction) created with {} time steps", amer_krl_bi.n);
    let amer_krl_bi_price = amer_krl_bi.price(&calc_method_bi);
    println!("Option price: {}", amer_krl_bi_price);

    // CRR & KRL with 10000 time steps (Combinatorial)
    let amer_crr_comb = CRR::new(Rc::clone(&amer_tree), 10000, Box::new(payoff_function2));
    println!("CRR (Combinatorial) created with {} time steps", amer_crr_comb.n);
    let amer_crr_comb_price = amer_crr_comb.price(&calc_method_comb);
    println!("Option price: {}", amer_crr_comb_price);
    let amer_krl_comb = KRL::new(Rc::clone(&amer_tree), 10000, Box::new(payoff_function2));
    println!("KRL (Combinatorial) created with {} time steps", amer_krl_comb.n);
    let amer_krl_comb_price = amer_krl_comb.price(&calc_method_comb);
    println!("Option price: {}", amer_krl_comb_price);
}
