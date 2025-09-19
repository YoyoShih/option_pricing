# Option Pricing in Rust

This project implements option pricing models using Rust, focusing on binomial and trinomial lattice methods. It is designed for educational and research purposes, providing efficient and robust algorithms for pricing vanilla and exotic options.

## Features
- **CRR (Cox-Ross-Rubinstein) Binomial Tree Model**
- **KRL (Karatzas-Ruf-Laplace) Trinomial Tree Model**
- (hasn't fully implement) Combinatorial methods for both Tree models targetting for higher running efficiency
- (hasn't fully implement) Flexible payoff functions (European, American, barrier, etc.)
- Modular code structure for easy extension

## Project Structure
- `src/` — Rust source code
  - `main.rs` — Entry point
  - `tree.rs` — General setting for tree model
  - `crr.rs` — Core CRR tree models and pricing logic
  - `krl.rs` — Core KRL tree models and pricing logic
- `ref/` — Reference files for this project
- `Cargo.toml` — Rust project configuration

## Usage
1. **Build the project:**
   ```powershell
   cargo build --release
   ```
2. **Run the project:**
   ```powershell
   cargo run --release
   ```
3. **Modify parameters:**
   Edit `main.rs` to set option parameters, model type, and payoff function.

## Example
```rust
fn payoff_function2(price: f64, strike: f64) -> f64 { price - strike }

let euro_tree = Rc::new(Tree::new(100.0, 100.0, 1.0, 0.05, 0.2, 0.0, OptionStyle::European));
println!("Tree of European option created with initial stock price: {}", euro_tree.s);

let calc_method_bi = CalcMethod::BackwardInduction;
let euro_crr_bi = CRR::new(Rc::clone(&euro_tree), 10000, Box::new(payoff_function2));
println!("CRR (Backward Induction) created with {} time steps", euro_crr_bi.n);
let euro_crr_bi_price = euro_crr_bi.price(&calc_method_bi);
println!("Option price: {}", euro_crr_bi_price);
```

## References
- "Efficient and Robust Combinatorial Option Pricing Algorithms on the Trinomial Lattice for Polynomial and Barrier Options" (PDF in `ref/`)


## Author
- I-Shao Shih

---
Feel free to modify and extend the models for your research or coursework needs.
