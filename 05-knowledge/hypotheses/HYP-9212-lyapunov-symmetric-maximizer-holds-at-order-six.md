---
id: HYP-9212
title: "The Lyapunov symmetric-maximizer conjecture holds at order six: for every real 6x6 A, ||L_A restricted to Skew_6|| <= ||L_A restricted to Sym_6||, L_A(X) = AX + XA^T, Frobenius-induced norms"
status: >
  OPEN (the only unresolved order: true for n <= 5, Chen-Tian 2015; false for n >= 7, Kressner-Vandereycken
  arXiv:2608.20875 + padding). FINITE-NUMERICAL evidence with working positive controls (below); no proof.
source: mac-mini-2026-10-06-sixseven, 05-knowledge/results/sixes_and_sevens_20261006.md, section 6
related:
  - 07-reflections/octonions-at-the-center-s6-bott-lyapunov-and-class-rank-frontiers-root-20260824.md (section 8)
  - 04-computation/octonion_s6_lyapunov_exact_audit_20260824.py
  - 04-computation/experiments/sixseven_20261006_lyapunov_adam.py
  - 04-computation/experiments/sixseven_20261006_lyapunov_six.py
---

# HYP-9212 — the symmetric maximizer survives at order six

**Conjecture.** For every `A ∈ R^{6x6}`, `||L_A|_{Skew_6}|| <= ||L_A|_{Sym_6}||`.

## Evidence (FINITE-NUMERICAL)

1. **Adam, the discovery recipe of Kressner–Vandereycken, with positive controls.** The objective is the gap `g = σ_skew - σ_sym` on the unit sphere, with tangent-projected singular-vector gradients. Adam uses `β = (0.9, 0.999)`, step 0.02 (800 iterations) then 0.006 with quadratic decay (800 more), as in KV section 4.

   | order | runs with a positive gap | best gap (ratio − 1) |
   |---|---|---|
   | 9 | 4 / 9 | `1.9e-6` |
   | 8 | 25 / 100 | `7.6e-3` |
   | 7 | 37 / 1712 (three step-size settings) | `7.7e-3` |
   | **6** | **0 / 6000** (step sizes 0.01, 0.02, 0.05; two run lengths) | the best runs end on the equality ridge `g = 0` to machine precision (`~2e-15`) |

   If a run succeeded at order 6 with probability `>= 0.1%`, zero hits in 6000 would have probability `< 0.25%`.
2. **The KV basin cannot be stripped to order 6.**
   * Put `μ(A) = λ_min(A^T A + A A^T)/|A|_F^2`. It vanishes iff `A` is orthogonally `A_6 ⊕ 0`, which is a counterexample iff `A_6` is (padding block identity; reflection, section 8.1).
   * Maximise `h = log(σ_skew/σ_sym)` subject to `μ <= m` along the KV branch: `h ≈ 0.146 m^1.47` as `m -> 0`.
   * The stripped 6x6 block always has `h ≈ -4 h_7 < 0`, and local re-optimisation never crosses 0 (26 values of `m`, 7 starts each).
3. **Structure of the counterexamples found (order 7; 20 collected in 1000 further runs).**
   * The strongest ones fall into one recurrent top basin (gap `~6.2e-3`, reached 7 times), with normalised singular values `(1, .85, .85, .74, .69, .14, .05)`. That is five large ones and two small ones, like KV's rank-5 matrix.
   * Weaker counterexamples do **not** share this profile. One has smallest singular value `0.185`, another `(.., .365, .003)`.
   * So "a rank-5 core plus a 2-dimensional near-kernel" describes the top basin only. It is not a necessary condition and gives no mechanism for order 6.
4. **Hostile control (why item 1 matters).**
   * L-BFGS and random structured starts fail even at `n = 7`: 12,845 starts, best `h = -2e-11`, the equality ridge.
   * A rank-constrained Adam (`A = U V^T`) also failed at `n = 7` rank 5 (0/200). Its order-6 runs (0/1200) are therefore uninformative.

## What would settle it

* **Counterexample.** An `n = 6` matrix with `g > 0`, then KV's exact separator (rational `c` with `c D_S - B_S^T D_S B_S` positive definite by an exact `LDL^T`, and a rational skew witness `K` with `||L_A K||^2 > c ||K||^2`).
* **Proof.** A reduction of the order-6 problem to a structured family where the conjecture holds (normal; nonnegative, nonpositive, tridiagonal — Feng–Lam–Yang–Li). Item 3 rules out the naive route "every counterexample has a 2-dimensional near-kernel".
