---
id: THM-4755
title: "Blow-up gains through the skew spectrum, the lexicographic law, and the ternary tower product. Let T* be a regular tournament on d vertices with skew matrix S = A - A^T, whose nontrivial eigenvalues are +-i sigma_k. (1) Wormald's blow-up gain rho_T(T*) = lim E[H(T(T*))]/g(dn) equals det(dI - S)/|det(idI - S)| = prod_k (d^2 + sigma_k^2)/(d^2 - sigma_k^2). This is the formula claimed by Bhattacharjee (2026) for Wormald's Conjecture 6; it is rederived here from the block walk, whose transition matrix restricted to the mean-zero space is S/d. (2) Lexicographic law: for regular Q and B the nontrivial skew spectrum of Q[B] is |B| sigma(Q) together with |Q| copies of sigma(B); so the gain factorises by levels. (3) Ternary towers S_k = C_3^[k]: rho_T(S_k) = prod_(i<k) ((3^(2i+1)+1)/(3^(2i+1)-1))^(3^i), so rho_T(C_3[C_3]) = 2 (28/26)^3 = 5488/2197 (Wormald's entry for R'_9), and the limit is kappa_S = 2.79150988819834776750... (4) FINITE-EXACT: H(S_k) for k <= 6 (729 vertices); the normalised counts 2, 2.2286, 2.5724, 2.7135, 2.7649, 2.7826 extrapolate to kappa_S to within 3.1e-7; v_3(H(S_k)) = (3^k-1)/2 = v_3|Aut| for k = 1, 3, 4, 5, 6 and one more at k = 2. (5) log rho_T = sum_k 2 artanh(sigma_k^2/d^2) = -sum_(m>=0) tr(S^(4m+2))/((2m+1) d^(4m+2)); doubly regular T* give ((d+1)/(d-1))^((d-1)/2) -> e, and cyclic R_d give the limit cosh(1)/cos(1) (Bhattacharjee)."
status: >
  Audited 2026-10-22 (05-knowledge/results/sierpinski_audit_20261022.md): no mathematical error (mac-mini-2026-10-22-sierp).
  CITED + INDEPENDENTLY REDERIVED: statement 1. Wormald, "Tournaments with many Hamilton cycles", arXiv:2610.16245 (posted
  2026-10-12), Theorems 1-5 and the table on p. 5. D. Bhattacharjee, github.com/creelie/Breakthrough-problems-1, release v1.0.0
  (2026-10-14, MIT, unreviewed): H = d det(d^2 I + S^2)/det(d^2 I - S^2) for Conjecture 6, the closed form for R_d, the limit
  cosh(1)/cos(1) for Conjecture 2, and the doubly regular lower bound. Its Lean development proves Conjecture 6 only given its
  Proposition 3.1, which identifies Wormald's analytic H with a determinant and is not formalised.
  The derivation here (section 1) is independent and different; it is a derivation through the local limit theorem for
  finite Markov chains, not a formal proof. It agrees exactly with every regular-tournament entry of Wormald's table
  (R_3 ... R_13, R'_9, the non-circulant 7-vertex example; R_21 ... R_101 to all printed digits). Wormald lists seven 9-vertex
  'circulants' (p. 6): six are not tournaments as printed (each connection set contains both r and 9 - r) and lie below the
  Jensen bound 2.44; the seventh, D = {1,2,3,5}, is a tournament with rho = 451696/183937 = 2.45571, matching his 2.4557.
  This observation is Bhattacharjee's (Remark on Wormald's circulants in the release preprint), not new here.
  AUDIT: independently audited 2026-10-22 (05-knowledge/results/sierpinski_audit_20261022.md): no mathematical error;
  credits, digits and labels corrected as required there.
  PROVED: statements 2, 3 and 5 (linear algebra on top of statement 1).
  FINITE-EXACT: statement 4. The C_3 substitution recursion is THM-3121 (codex) and Session24's DP; Session24 holds the values
  through 243 vertices. New here are the value at 729 vertices, the normalised sequence and its limit.
session: mac-mini-2026-10-22-sierp
source: 05-knowledge/results/sierpinski_steiner_collatz_tournaments_20261022.md
scripts:
  - 04-computation/experiments/sierpinski_ocf_20261022/wormald_check.gp, wormald_check2.gp (+ .out): statement 1 against Wormald's table
  - 04-computation/experiments/sierpinski_ocf_20261022/tower_rho.gp (+ .out): statement 3, exact for k <= 3
  - 04-computation/experiments/sierpinski_ocf_20261022/ternary_fast.py (+ ternary_fast_k6.out), ternary_recursion.py, ternary_H.out: statement 4
depends_on:
  - Wormald 2026, Theorems 1, 4, 5 (blow-ups with random blocks; regular blocks give e^(1/d))
  - local limit theorem for occupation counts of a finite aperiodic irreducible Markov chain
related:
  - 01-canon/theorems/THM-4754-three-adic-sierpinski-tournament-collatz-arc-flip-steiner-midpoint.md
  - 05-knowledge/hypotheses/HYP-9364-self-blow-up-principle-tower-limit-rotational-cosh-over-cos.md
  - 01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md
  - 01-canon/theorems/THM-4751-circulant-h-landscape-paley-peak-at-class-number-distance.md
---

# THM-4755 — Blow-up gains through the skew spectrum

## 1. The gain is a determinant ratio (Bhattacharjee's formula, rederived)

Blow up each vertex `i` of a regular tournament `T*` (`d` vertices) into a random tournament `T_i` on `n` vertices, with all arcs between blocks following `T*`.

**The weighted block walk.**
* A Hamiltonian path reads a block word of length `dn`. Each step stays in a block (the arc is random: probability 1/2), moves along an arc of `T*` (probability 1), or moves against one (probability 0).
* Against the random benchmark `2^(−(dn−1))`, a word therefore has weight `Π(1 for a stay, 2 for a forward move, 0 for a backward move)`. This is the walk with transition matrix `P = (I + 2A)/d = (J + S)/d`, whose row sums are 1.

**The density ratio.**
* So `E[H]/g'(dn)` is the ratio of the probabilities that the walk and an i.i.d. uniform sequence have exactly `n` letters of each kind.
* By the local limit theorem it is `(det' Σ_iid / det' Σ_P)^(1/2)`, where the `Σ` are the occupation covariances on the mean-zero space `V = 1^⊥`.
* On `V`, `J = 0`, so `P|_V = S/d`, which is skew-symmetric.
* With `Σ_P = (1/d)[(I − Q)^(−1) + (I − Qᵀ)^(−1) − I]` and `Q = S/d`, one gets `Σ_P = (1/d)(I − Q)^(−1)(I − QQᵀ)(I − Qᵀ)^(−1)`.

Hence

    ρ_T(T*)^2 = |det_V(I − S/d)|^2 / det_V(I + S^2/d^2),  i.e.  ρ_T = Π_k (d^2 + σ_k^2)/(d^2 − σ_k^2) = det(dI − S)/|det(idI − S)|.

**Consequences.**
* No normality assumption on `A` is needed.
* Both determinants are integers up to a factor `d`, so `d/H = ρ_T^2` is the square of a rational. This is Wormald's Conjecture 6, proved here modulo the identification of his `H` with this local-limit constant (also Bhattacharjee's Proposition 3.1).
* `σ_k ≤ d − 1 < d` (row sums of `|S|`), so `H ≠ 0`.

**Agreement with Wormald's table** (`wormald_check*.out`):

| `T*` | `ρ_T` here | Wormald |
|---|---|---|
| `R_3` | 2 | 2 |
| `R_5` | 44/19 | 44/19 |
| `R_7` | 3032/1231 | 3032/1231 |
| `R_9` | 437584/171769 | 437584/171769 |
| `R'_9 = C_3[C_3]` | 5488/2197 | 5488/2197 |
| `R_11` | 108859552/41833969 | same |
| `R_13` | 41637757888/15769379963 | same |
| non-circulant 7-vertex | 3088/1287 | 3088/1287 |
| `R_21, R_41, R_81, R_101` | 2.7215188…, 2.7867038…, 2.820801899…, 2.827747088… | identical digits |

## 2. The lexicographic law

For regular `Q` (`q` vertices) and regular `B` (`b` vertices), `S(Q[B]) = S_Q ⊗ J_b + I_q ⊗ S_B`.
* On `u ⊗ 1` it acts as `b S_Q`.
* On `u ⊗ w` with `w ⊥ 1` it acts as `S_B`.

So the nontrivial spectrum is `b·σ(Q)` together with `q` copies of `σ(B)`, and

    ρ_T(Q[B]) = ρ_T(Q) · Π_(σ in σ(B)) ((qb)^2 + σ^2)/((qb)^2 − σ^2))^q.

## 3. The ternary tower

For `S_k = C_3[S_(k−1)]`, the level-`i` modes are `3^i` pairs with `σ^2 = 3^(2(k−i)−1)`, against `d^2 = 3^(2k)`. So

    ρ_T(S_k) = Π_(i=0)^(k−1) ((3^(2i+1) + 1)/(3^(2i+1) − 1))^(3^i)  →  κ_S = 2 · (28/26)^3 · (244/242)^9 · (2188/2186)^27 ⋯ = 2.79150988819834776750…

This is exact for `k ≤ 3` against the determinant (`tower_rho.out`).

* The first factor is `(3 + 1)/(3 − 1)`, i.e. Wormald's Theorem 2.
* The second is `((3³ + 1)/(3³ − 1))^3`. (Reading 28 and 26 as `7 + 21` and `25 + 1` is NUMEROLOGY.)
* Equivalently, `log κ_S = 2 Σ_(i≥0) 3^i artanh(3^(−2i−1)) = 1 + 2 Σ_(m≥3 odd) 3^(−m)/(m(1 − 3^(1−2m)))`. The leading 1 is the factor `e`.
* `κ_S` exceeds Wormald's regular-block value for `C_3[C_3]`, `(5488/2197)·e^(1/9)`, by a relative `4·10^(−7)`.

## 4. Exact tower counts (FINITE-EXACT)

The C_3 substitution recursion (THM-3121, Session24) gives `H(S_k)` exactly for `k ≤ 6` (`ternary_fast_k6.out`). The values for `k ≤ 3` also come from inclusion–exclusion (`ternary_H.out`). The values for `k ≤ 5` were recomputed independently in the audit (`sierpinski_audit_20261022.md`). Since `Aut` acts freely on Hamiltonian paths (LEM-003), `v_3(H) ≥ v_3|Aut| = (3^k − 1)/2`.

| `k` | `n` | `H·2^(n−1)/n!` | `v_3(H)` | `v_3|Aut|` |
|---|---|---|---|---|
| 1 | 3 | 2 | 1 | 1 |
| 2 | 9 | 2.228571 | 5 | 4 |
| 3 | 27 | 2.572358 | 13 | 13 |
| 4 | 81 | 2.713530 | 40 | 40 |
| 5 | 243 | 2.764891 | 121 | 121 |
| 6 | 729 | 2.782564 | 364 | 364 |

* Successive differences shrink by ratios 0.41, 0.36, 0.34, which tend to 1/3.
* Richardson extrapolation in `1/n` gives 2.79150957576 (four terms), against `κ_S = 2.79150988820`.
* This supports the conjecture HYP-9364 (a) that the towers' own normalised counts converge to the limit of their blow-up gains. Richardson through 3 or 5 levels gives −5.2·10^(−6) and +2.3·10^(−7), so the agreement is robust to about 5–6 digits; it is not a proof.

## 5. Spectral moments, Paley and the rotational limit

* `log ρ_T = Σ_k 2 artanh(σ_k^2/d^2) = −Σ_(m≥0) tr(S^(4m+2))/((2m+1) d^(4m+2))`. Every term is positive, and the `m = 0` term is `(d − 1)/d`.
* By Jensen, `ρ_T ≥ ((d+1)/(d−1))^((d−1)/2)`, with equality exactly for doubly regular `T*` (Bhattacharjee). For Paley `T*` this tends to `e`.
* For cyclic `R_d`, `σ_j/d → 2/(πj)` for odd `j`, and `σ_j/d → 0` for even `j`. So `ρ_T(R_d) → Π_(j odd) (1 + 4/(π² j²))/(1 − 4/(π² j²)) = cosh(1)/cos(1) = 2.855957892565…`.
* This is Wormald's fitted constant 2.8559579 (Bhattacharjee's closed form `((d+1)^d + (d−1)^d)/((d+i)^d + (d−i)^d)` has the same limit).
* **Reading.** A flat spectrum gives `e`; the most concentrated spectrum gives `cosh(1)/cos(1)`.

## Addendum: Wormald's Conjecture 1 for regular bases, d ≤ 13 (FINITE-EXACT)

Over every regular tournament on `d` vertices, `ρ = Π (d² + σ_k²)/(d² − σ_k²)` has a unique maximum at the rotational `R_d` (`wormald_conj1.out`, `wormald_conj1_13.out`). The cases `d ≤ 9` are already in Bhattacharjee's release (its table of the regular tournaments on `d ≤ 9` vertices, with `R_d` maximal); `d = 11, 13` are new:

| `d` | regular classes | `ρ(R_d)` | runner-up |
|---|---|---|---|
| 5 | 1 | 2.3157894737 | — |
| 7 | 3 | 2.4630381803 | 2.3993783994 |
| 9 | 15 | 2.5475143943 | 2.4981158008 (443488/177529; `C_3[C_3]` = 5488/2197 = 2.4979517524 is third) |
| 11 | 1223 | 2.6021808258 | 2.5682251025 |
| 13 | 1,495,297 | 2.6404182020 | 2.6170585197 |

**Context.**
* Conjecture 1 asks for the maximum over all tournaments on `d` vertices. Non-regular bases are not covered by the determinant formula: their balanced block counts are exponentially unlikely, so `ρ = 0` in the limit. That reduction is not proved here.
* `ρ` is a Schur-convex function of `(σ_k²)` with the fixed sum `d(d − 1)/2`. But the squared skew spectrum of `R_d` does not majorise those of all regular `d`-tournaments: there is 1 violation at `d = 7`, 7 at `d = 9` and 626 at `d = 11` (`majorization_test.out`).
* So the conjecture is not a Schur-convexity consequence. `R_d` wins through its dominant top value `σ_1² ≈ 4d²/π²`, even where later partial sums are beaten.
