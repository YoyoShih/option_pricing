# Option Pricing in Rust

This project implements option pricing models using Rust, focusing on binomial and trinomial lattice methods. It is designed for educational and research purposes, providing efficient and robust algorithms for pricing vanilla and exotic options.

## Features
- **CRR (Cox-Ross-Rubinstein) Binomial Tree Model**
- **KRL (Karatzas-Ruf-Laplace) Trinomial Tree Model**
- (hasn't implement) Combinatorial methods for both Tree models targetting for higher running efficiency
- (hasn't implement) Flexible payoff functions (European, American, barrier, etc.)
- Modular code structure for easy extension

## Project Structure
- `src/` — Rust source code
  - `main.rs` — Entry point
  - `tree.rs` — Core tree models and pricing logic
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
let tree = Tree::new(100.0, 100.0, 1.0, 0.05, 0.2, 0.0);
let crr = CRR::new(Rc::new(tree), 100, Box::new(|s, x| f64::max(s - x, 0.0)));
let price = crr.price();
println!("Option price: {}", price);
```

## References
- "Efficient and Robust Combinatorial Option Pricing Algorithms on the Trinomial Lattice for Polynomial and Barrier Options" (PDF in `ref/`)


## Author
- I-Shao Shih

---
Feel free to modify and extend the models for your research or coursework needs.
