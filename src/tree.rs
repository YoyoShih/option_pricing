use std::f64;

// Define the OptionStyle enum
#[derive(PartialEq)]
pub enum OptionStyle {
    European, // default, can only be exercised at maturity
    American, // can be exercised at any time before maturity
}
#[derive(PartialEq)]
pub enum OptionType {
    Call,
    Put,
}
pub struct OptionSpec {
    pub style: OptionStyle,
    pub kind: OptionType,
}
// Provide a default implementation for OptionSpec
impl Default for OptionSpec {
    fn default() -> Self {
        OptionSpec {
            style: OptionStyle::European,
            kind: OptionType::Call,
        }
    }
}

#[derive(PartialEq)]
pub enum FiniteDifferenceType {
    Explicit,
    Implicit,
}

// Define the CalcMethod enum
#[derive(PartialEq)]
pub enum CalcMethod {
    BackwardInduction,
    Combinatorial,
    FiniteDifference {
        kind: FiniteDifferenceType,
        n_s: u32,
        n_t: u32,
        s_max: f64,
        s_min: f64,
    },
}

// Provide a default implementation for CalcMethod
impl Default for CalcMethod {
    fn default() -> Self {
        CalcMethod::BackwardInduction
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
    pub option_spec: OptionSpec, // option style
}

// Implement methods for Tree
impl Tree {
    pub fn new(s: f64, x: f64, t: f64, r: f64, sigma: f64, q: f64, option_spec: OptionSpec) -> Self {
        Tree { s, x, t, r, sigma, q, option_spec }
    }
}