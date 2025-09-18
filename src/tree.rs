use std::f64;

// Define the OptionStyle enum
#[derive(PartialEq)]
pub enum OptionStyle {
    European, // default, can only be exercised at maturity
    American, // can be exercised at any time before maturity
}
// Provide a default implementation for OptionStyle
impl Default for OptionStyle {
    fn default() -> Self {
        OptionStyle::European
    }
}

// Define the Tree struct
pub struct Tree {
    pub s: f64, // initial stock price
    pub x: f64, // strike price
    pub t: f64, // time to maturity
    pub r: f64, // risk-free interest rate
    pub sigma: f64, // volatility
    pub q: f64, // dividend yield
    pub style: OptionStyle, // option style
}

// Implement methods for Tree
impl Tree {
    pub fn new(s: f64, x: f64, t: f64, r: f64, sigma: f64, q: f64, style: OptionStyle) -> Self {
        Tree { s, x, t, r, sigma, q, style }
    }
}