---
id: THM-4757
title: "Loop expansion of Hamiltonian paths in tournaments. (1) For every tournament, kappa(T) = H(T) 2^(n-1)/n! = sum_L N_L(S)/(n)_L, where N_L(S) is the sum of prod S_e over the directed linear forests of K_n with L arcs (S = A - A^T). (2) Moebius inversion over coincidences of run positions writes log of Phi_T(z) = sum_L N_L (z/n)^L as a sum of connected clusters graded by excess e = arcs - blocks; for regular T, clusters with a free run end vanish, so e >= 0 (no tree level). (3) One loop: the e = 0 clusters (single cycles) sum exactly to (1/2) log det(I - z^2 Shat^2) - (1/2) log det(I + z^2 Shat^2), Shat = S/n; at z = 1 this is log rho(T), the Wormald blow-up gain det(nI - S)/|det(inI - S)|. (4) Two loops: the e = 1 clusters (figure-eight, theta, dumbbell) sum to sum_v [2 G_vv H_vv - (3/2) G_vv^2] + sum_(u,v) G_uv (H_uv^2 - O_uv^2 - (2/3) G_uv^2) - d^T G d, with G = -Shat^2 (I + Shat^2)^-1, H = Shat^2 (I - Shat^2)^-1, O = Shat (I - Shat^2)^-1, d = diag(H - G); the dumbbell vanishes for vertex-transitive T."
status: >
  PROVED, pre-audit (mac-mini-2026-10-22-sierp): (1)-(3) are elementary identities, (3) as a formal power series in z,
  convergent at z = 1 because every |x_k| < 1. (4) is proved as a formal power series modulo the kernel computation:
  each kernel is 4-periodic in every path length (rational generating functions with simple poles at fourth roots of
  unity, below), and the periodic values are fixed by exhaustive enumeration (two_loop_kernels.py). Independent check:
  the closed formula's z^4 and z^6 coefficients equal the brute-force excess-1 cluster sums for R_15 and P_19, and the
  brute-force excess-0 sums equal the log rho coefficients for L <= 6 (cluster_check.py, two_loop_series_check.py).
  The asymptotic consequence kappa = rho (1 - y/n + O(n^-2)) is NOT proved here; it is HYP-9365.
  PRIOR ART: (1) is inclusion-exclusion over forced arcs (classical in permanent/Hamiltonicity counting). The blow-up
  gain rho is Wormald's (arXiv:2610.16245) in Bhattacharjee's determinant form (GitHub release v1.0.0, 2026-10-14);
  THM-4755 rederived it for blow-ups via a block walk. The single-tournament reading of rho as the one-loop part of the
  forest gas, the excess grading and the two-loop formula are new here as far as the repo and these two sources go.
session: mac-mini-2026-10-22-sierp
source: 05-knowledge/results/sierpinski_steiner_collatz_tournaments_20261022.md
scripts:
  - 04-computation/experiments/sierpinski_ocf_20261022/cluster_check.py (+ cluster_check_rot15_L6.out, cluster_check_paley19_L6.out)
  - 04-computation/experiments/sierpinski_ocf_20261022/cluster_shapes.py
  - 04-computation/experiments/sierpinski_ocf_20261022/two_loop_kernels.py, theta_table.py, dumbbell_table.py
  - 04-computation/experiments/sierpinski_ocf_20261022/two_loop_series_check.py
related:
  - 05-knowledge/hypotheses/HYP-9365-two-loop-law-for-hamiltonian-paths-in-regular-tournaments.md
  - 05-knowledge/hypotheses/HYP-9364-self-blow-up-principle-tower-limit-rotational-cosh-over-cos.md
  - 01-canon/theorems/THM-4755-blow-up-gains-skew-spectrum-lexicographic-law-ternary-tower.md
  - 01-canon/theorems/THM-4756-forward-step-polynomial-mod-3-is-h-times-eulerian.md
---

# THM-4757 — Loop expansion of Hamiltonian paths

**Notation.** `T` is a tournament on `[n]`, `A` its adjacency matrix, and `S = A − A^T` the skew matrix (entries ±1 off the diagonal). Write `Ŝ = S/n`. Its eigenvalues are `i x_k` with `x_k` real and `|x_k| ≤ (n−1)/n < 1`. A *directed linear forest* is a set of arcs of the complete digraph on `[n]` forming vertex-disjoint directed paths. Set `N_L(S) = Σ_{|F| = L} Π_{e∈F} S_e`.

## 1. The forest identity

For every tournament,

    2^(n−1) H(T) = Σ_F (n − |F|)! Π_(e∈F) S_e,      i.e.      κ(T) := H(T) 2^(n−1)/n! = Σ_(L=0)^(n−1) N_L(S)/(n)_L.

*Proof.* For an ordering `σ` of `[n]`, `Π_i (1 + S_(σ_i σ_(i+1)))` equals `2^(n−1)` if `σ` is a Hamiltonian path of `T` and 0 otherwise. Expand the product and exchange the sums. A linear forest with `L` arcs has `n − L` components, counting isolated vertices, and lies in exactly `(n − L)!` orderings. ∎

The top term is `N_(n−1) = Σ_σ (−1)^(#backward steps) = (−1)^(n−1) F(T, −1)`, the alternating value of the forward-step polynomial of THM-4756.

## 2. Clusters and excess

Put `Φ_T(z) = Σ_L N_L (z/n)^L`.
* Write a forest as `r` labelled runs, i.e. injective vertex sequences of lengths `ℓ_j ≥ 1`, with weight `1/r!`.
* Replace injectivity by Möbius inversion on the partition lattice of all run positions: a partition `π` gets the weight `μ(π) = Π_B (−1)^(|B|−1)(|B|−1)!` and the free homomorphism sum `Z_π = Σ_(f: blocks → [n]) Π_arcs (zŜ)_(f(·) f(·))`. Each `Z_π` factorises over the connected components of the quotient multigraph.
* So, by the exponential formula,

      log Φ_T(z) = Σ_(e ≥ 0) X_e(z),

  where `X_e` sums `μ(π) Z_π / p!` over connected clusters (`p` runs and a partition) of **excess** `e = #arcs − #blocks`, i.e. cyclomatic number `e + 1`.

**Regular `T`.** If `S·1 = 0`, a block that is a single run end gives `Σ_v Ŝ_(vw) = 0`. So only quotient graphs of minimum degree ≥ 2 survive, and `e ≥ 0`. There is no tree level. The grading is a loop expansion with `1/n` as the loop-counting parameter: an excess-`e` cluster is `O(n^(−e))`.

## 3. One loop: the blow-up gain

**Theorem.** For regular `T`,

    X_0(z) = (1/2) Σ_k log((1 + z² x_k²)/(1 − z² x_k²)) = (1/2) log det(I − z²Ŝ²) − (1/2) log det(I + z²Ŝ²).

So `exp X_0(1) = ρ(T) = det(nI − S)/|det(inI − S)|`, Wormald's blow-up gain (THM-4755).

*Proof.* Excess-0 connected quotients have every vertex of degree 2, so they are single cycles. Interior run positions are singleton blocks, and the `2p` run ends are matched in pairs, so `μ = (−1)^p`. Among the matchings of `p` labelled runs, the ones that close into one cycle number `(p−1)! 2^(p−1)`: a cyclic order of runs 2..p relative to run 1, and their orientations. A run traversed against its direction contributes `(−Ŝ)^ℓ`. All factors are powers of `Ŝ`, so the cycle's value is `Σ_k y_k^(ℓ_1) Π_(j≥2) (y_k^(ℓ_j) + (−y_k)^(ℓ_j))` with `y_k = i z x_k`. Summing the lengths gives

    X_0 = Σ_k Σ_(p≥1) ((−1)^p/p) A(y_k) B(y_k)^(p−1) = −Σ_k (A/B) log(1 + B),    A = y/(1−y),  B = 2y²/(1−y²).

The eigenvalues come in pairs `±i x`, so only the even part `y²/(1−y²) = B/2` of `A` survives. Hence `X_0 = −(1/2) Σ_k log(1 + B(y_k))`, and `1 + B = (1 + y²)/(1 − y²)` with `y² = −z²x²`. ∎

**Remarks.**
* The series `log ρ = Σ_k Σ_(m odd) x_k^(2m)/m` has no `x⁴, x⁸, …` terms. The 4-cycles cancel among themselves, which is the source of the "odd `m` only" pattern in THM-4755.
* Brute force over every `(ℓ-tuple, π)` with `L ≤ 6` reproduces `Σx², 0, Σx⁶/3` to 12 digits for `R_15` and `P_19` (`cluster_check_*_L6.out`).

## 4. Two loops: figure-eight, theta, dumbbell

Every excess-1 quotient is a subdivision of one of the three two-loop cores:
* the **figure-eight**: two loops at a degree-4 vertex;
* the **theta**: three paths between two degree-3 vertices;
* the **dumbbell**: two loops joined by a path.

Fix reference orientations of the core paths. The *kernel* `K` of a core with path lengths `(a, b, …)` is the sum over labelled run decompositions of `μ` times `Π_runs ((−1)^α + (−1)^(ℓ−α))`, where `α` is the number of arcs a run traverses against the reference.
* At a degree-2 vertex there is a pass-through (`μ = 1`) or a junction of two run ends (`μ = −1`).
* At a core vertex there is a set of pass-through pairs, with `μ = (−1)^(s−1)(s−1)!` for block size `s`.
* Closed runs are excluded.

Dividing by the automorphism groups of the labelled cores (orders 8, 12, 8), with `ε_a = (−1)^(a/2)`:

| core | kernel (zero unless stated) |
|---|---|
| figure-eight `(a, b)` | `8(ε_a + ε_b) − 12 ε_a ε_b` for `a, b` even |
| theta `(a, b, c)` | `4(ε_a + ε_b + ε_c) − 8 ε_a ε_b ε_c` all even; `−4 ε_c` if `a, b` odd, `c` even |
| dumbbell `(a, b; c)` | `−32 ε_c` if `a ≡ b ≡ 2 (mod 4)`, `c` even |

**Why the kernels are 4-periodic.** Along a core path, a junction has weight `−1`, an internal segment of length `λ` contributes `1 + (−1)^λ`, and an end segment contributes `(sz)^λ` with `s = ±1` its run's orientation sign. The resulting path generating functions are:
* one segment: `sz/(1 − sz)`;
* two or more segments: `−s_1 s_2 z²(1 − z²)/((1 − s_1 z)(1 − s_2 z)(1 + z²))`.

All have simple poles at fourth roots of unity and no polynomial part beyond `z⁰`, so their coefficients are 4-periodic for `λ ≥ 1`. A kernel is a finite sum over core configurations and orientation signs of products of such coefficients, so it is 4-periodic in each length. The table entries were enumerated on representatives of every residue class (`two_loop_kernels.py`, `theta_table.py`, `dumbbell_table.py`; lengths to 8–10, which is more than one period).

**Resummation.** Even paths with sign `ε_a` sum to `G = Σ ε_a Ŝ^a = −Ŝ²(I + Ŝ²)^(−1)`. Plain even paths give `H = Ŝ²(I − Ŝ²)^(−1)`, and odd paths give `O = Ŝ(I − Ŝ²)^(−1)`. Series composition is the matrix product, and parallel composition is the entrywise product. Hence

    X_1 = Σ_v [2 G_vv H_vv − (3/2) G_vv²]  +  Σ_(u,v) G_uv (H_uv² − O_uv² − (2/3) G_uv²)  −  dᵀ G d,     d = diag(H − G),

with `Ŝ → zŜ` for the series in `z`.
* The `z⁴` and `z⁶` coefficients match brute force: −0.19911 against −0.19910 and −0.060250 against −0.060250 for `R_15`; −0.16271 against −0.16270 and −0.0068139 against −0.0068140 for `P_19` (`two_loop_series_check.py`).
* For vertex-transitive `T`, `d` is constant and `G·1 = 0`, so the dumbbell vanishes. The figure-eight becomes `(1/n)[2 tr G tr H − (3/2)(tr G)²]`, a function of the spectrum alone. The theta becomes `(1/n) Σ_(k_1+k_2+k_3 ≡ 0)` of products of `g, h, o` at three frequencies: the first place where the additive structure of the connection set, and not only the spectrum, enters.

**Series–parallel ("least resistance") remark.** All three two-loop cores are series–parallel, and their values come from series and parallel composition of three propagators. The first non-series–parallel vacuum graph, `K_4`, appears at three loops (`e = 2`). That is where, for example, Paley character sums stop collapsing to polynomials in `q`.
