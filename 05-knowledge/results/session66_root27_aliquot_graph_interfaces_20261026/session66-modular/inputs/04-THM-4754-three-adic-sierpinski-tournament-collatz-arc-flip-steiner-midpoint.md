---
id: THM-4754
title: "The 3-adic Sierpinski tournament and Collatz. (1) S_k = Cay(Z/3^k, {x : unit part of x = 1 mod 3}) is the k-fold lexicographic power C_3^[k]; on Z_3 the same rule defines S_infinity. (2) The maps x -> 3x + b (b = 0, 1, 2) are its self-similarities: they embed S_(k-1) onto the residue classes, and class b dominates class b+1. So '3x+1' is one of the three Sierpinski maps, and its fixed point -1/2 = ...111 (base 3) is the all-ones address. (3) Multiplication by a unit lambda is an automorphism if lambda = 1 mod 3 and reverses every arc if lambda = 2 mod 3; in particular halving reverses S_infinity, because 2 = -1 mod 3. (4) Collatz arc-flip law: if x = y mod 2^t (same parity vector for t steps), the arc between C^t(x) and C^t(y) in S_infinity equals the arc between x and y exactly when the number of halvings is even. (5) The pasted Steiner rule x*y = -x-y is the midpoint (x+y)/2 in F_3; the triangular reduction is binomial averaging, which collapses to the two endpoints at depth 3^k. (6) Jacobsthal split: H(T) = 3 sum_k J_k alpha_k + I(Omega(T), -1), with J_k = (2^k - (-1)^k)/3, so that each OCF weight 2^k = 3 J_k + (-1)^k is a one-step Collatz relation (J_(2m) = (4^m - 1)/3 = 1, 5, 21, 85, 341 are the predecessors of 4^m). (7) By THM-4756 (which completes THM-086), at the orders n = 2*3^k + 1 = 3, 7, 19, 55, 163, ... every tournament has F(T,x) = H (1 + x^(3^k) + x^(2*3^k)) mod 3."
status: >
  PROVED (elementary), audited 2026-10-22 (05-knowledge/results/sierpinski_audit_20261022.md): statements 1-7. Statement 7 follows from THM-4756, which proves THM-086 in full
  (F(T,x) = H(T) A_n(x) mod 3), and is checked exhaustively here at n = 7 (all 456 classes).
  PRIOR ART: C_3^[k] as a Cayley tournament on F_3^k whose affine lines are cyclic triangles
  (05-knowledge/results/tournament_pseudoprimes_camion_fermat_20261021.md, line 253); C_3[C_3] = Z_9{1,3,4,7} (THM-4743,
  HYP-9353); |Aut(C_3^[k])| = 3^((3^k-1)/2) (Session24 PATHCOVERS); H = I(Omega,-1) mod 3 (HYP-1028..1030) and
  I(Omega,-1) = -reduced Euler characteristic of Ind(Omega) (THM-467, T274); the OCF on a path conflict graph gives Jacobsthal
  numbers (THM-1425); the base-4 repunits as Collatz predecessors of 4^m and the trunk (THM-4517, procgen notes). The repo has no
  Pascal-mod-3 Collatz object (survey 2026-10-22), so statements 2-5 as a Collatz dictionary are new here.
  The orders 3, 7, 19, 163 being Heegner numbers is NUMEROLOGY.
session: mac-mini-2026-10-22-sierp (owner prompt: 3^k+1 and k^3+1 against 3n+1, a pasted Steiner/Sierpinski reduction, {25,1} against {7,21}, 27 in Collatz)
source: 05-knowledge/results/sierpinski_steiner_collatz_tournaments_20261022.md
scripts:
  - 04-computation/experiments/sierpinski_ocf_20261022/threeadic_check.py (+ .out): statements 1-4
  - 04-computation/experiments/sierpinski_ocf_20261022/sierp_ocf.c, analyze1.py (+ ocf_n*.txt): statement 6 on every class n <= 9
  - 04-computation/experiments/sierpinski_ocf_20261022/fwd_mod3.py (+ .out): statement 7 at n = 7
related:
  - 01-canon/theorems/THM-4755-blow-up-gains-skew-spectrum-lexicographic-law-ternary-tower.md
  - 01-canon/theorems/THM-086-universal-taylor-zeros-mod3.md
  - 01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md
---

# THM-4754 — The 3-adic Sierpinski tournament and Collatz

## 1–2. The tournament and its self-similarity

For `x` in `Z/3^k` (or `Z_3`), `x ≠ 0`, write `x = 3^v u` with `3 ∤ u`, and call `u mod 3` the **sign** of `x`. Let `S_k` be the tournament with `x → y` iff the sign of `y − x` is 1.

* **It is a tournament.** The sign of `−z` is `−(sign z)`, so exactly one of `x → y`, `y → x` holds.
* **Self-similarity.**
  * For `x, y` in `Z/3^(k−1)`, `(3y + b) − (3x + b) = 3(y − x)` has the same sign as `y − x`. So `x ↦ 3x + b` is an isomorphism of `S_(k−1)` onto the residue class `b`.
  * If `x ≡ b` and `y ≡ b + 1 (mod 3)`, then `y − x ≡ 1`, so class `b` dominates class `b + 1`.
  * Hence `S_k = C_3[S_(k−1), S_(k−1), S_(k−1)] = C_3^[k]`.
* **The map 3x + 1.** It is the self-similarity onto class 1. Its unique fixed point in `Z_3` is `−1/2 = …1111` (base 3), the vertex whose address is all ones. Iterating `x ↦ 3x + 1` from `x` gives `3^k x + (3^k − 1)/2`: the base-3 repunit is the address of the depth-`k` sub-block reached along ones.

## 3. Multipliers and halving

For a unit `λ`, the sign of `λz` is `(λ mod 3)·(sign z)`. So:
* `x ↦ λx` preserves every arc if `λ ≡ 1 (mod 3)`;
* it reverses every arc if `λ ≡ 2 (mod 3)`.

Since `1/2 ≡ 2 (mod 3)`, **halving reverses every arc** of `S_∞`. Multiplication by 3 preserves arcs.

## 4. The Collatz arc-flip law

Let `C(x) = 3x + 1` for odd `x` and `x/2` for even `x`.
* If `x ≡ y (mod 2^t)`, then `x` and `y` follow the same branches for `t` steps.
* With `o` odd steps and `e` halvings, `C^t(y) − C^t(x) = 3^o (y − x)/2^e`. Its sign is `(−1)^e·sign(y − x)`.

So the `S_∞`-arc between the two trajectories flips exactly once per halving. For the shortcut map `T(x) = (3x+1)/2` or `x/2`, every step flips it.

Checked on 2000 random pairs in `Z/3^12` (`threeadic_check.out`).

## 5. The Steiner rule is the midpoint

In `F_3`, `1/2 = −1`, so the pasted rule is `x ⋆ y = −x − y = (x + y)/2`.
* It is the third point on the line through `x, y` in `AG(1,3)`. Lines of `AG(k,3)` are cyclic triangles of `S_k` (prior art, see status).
* The triangular reduction of a string of length `n` is the binomial average `2^(−(n−1)) Σ_i C(n−1, i) x_i` over `F_3`. Up to sign this is the pasted formula.
* At `n − 1 = 3^k`, the Frobenius identity `(1 + z)^(3^k) = 1 + z^(3^k)` leaves only the endpoints.

## 6. The Jacobsthal split of the OCF

With `J_k = (2^k − (−1)^k)/3` (0, 1, 1, 3, 5, 11, 21, 43, 85, …) we have `2^k = 3J_k + (−1)^k`. Summing the OCF `H = Σ_k 2^k α_k` gives

    H(T) = 3 Σ_k J_k α_k + I(Ω(T), −1).

* **Collatz meaning.** For even `k = 2m`, `3J_(2m) + 1 = 4^m`: the base-4 repunits `1, 5, 21, 85, 341` are the odd numbers sent to a power of 4 by `3n + 1`. For odd `k`, `3J_k − 1 = 2^k`.
* **The forbidden values.** `21 = J_6 = (4³ − 1)/3 → 64`.
* **The mod-3 reduction.** `H ≡ I(Ω, −1) (mod 3)`, which is minus the reduced Euler characteristic of the independence complex (prior art).
* **Checked** on every class with `n ≤ 9` (`sierp_ocf.c`, which aborts on any mismatch).

## 7. The Steiner collapse in the forward-edge polynomial (by THM-4756)

THM-4756 proves `F(T, x) ≡ H(T)·A_n(x)`, which is `H·(x − 1)^(n−1) (mod 3)` for odd `n`; this is THM-086's statement, now proved.
* For `n − 1 = 2·3^k`, Lucas gives `(x − 1)^(2·3^k) ≡ 1 + x^(3^k) + x^(2·3^k) (mod 3)`.
* So the numbers of permutations with `i` `T`-forward steps vanish mod 3 except at `i = 0, 3^k, 2·3^k`, where they equal `H`.
* This holds for all 456 classes at `n = 7` (`fwd_mod3.out`).
* The orders `2·3^k + 1` are `3, 7, 19, 55, 163, …`. That four of them are Heegner numbers is NUMEROLOGY.
