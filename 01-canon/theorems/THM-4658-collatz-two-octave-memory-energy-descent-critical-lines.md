---
id: THM-4658
title: "Two-octave memory and 4-adic energy descent for the Collatz mod-3 discrepancy: H(4N) = H(N) + H_(2,2)(3N) - H_(2,8)(6N) exactly, so E_(4X) = 4E_X + 8C_X + 4D_X + O(X^(3/2)) and sqrt(E_(4X)) <= 2 sqrt(E_X) + 2 sqrt(D_X) + 4 sqrt(3X); the class-sum recursion C(A;N) = -C(A/2; floor(N/2)) - C((3A+1)/2; 3 ceil(N/2)) makes every power-law main term (density or excess) an eigen-distribution of a 3-adic transfer operator (up to a rank-one correction), whose two branches have squared L^2(Z_3)-norms summing to 4^(-sigma)(1 + 3^(2 sigma - 1)) - equal to 1 exactly on Re s = 1/2 and Re s = 1"
status: >
  PROVED (elementary consequences of the exact refinement identity (2) of collatz_wavelet_energy_20261013d.md, here in
  measure form). FINITE-EXACT control: the energy transport (2) checked on every octave 2^10..2^30 of an exact census
  through 2^34 (agreement to 4 decimals; census reproduces S(2^34) = 24953018 and H(2^34) = 26869 of the repository's
  independent 2^34 census). No bound on the actual Collatz discrepancy is claimed.
source: opus-2026-10-15-S37
depends_on: [collatz_wavelet_energy_20261013d.md (identity (2), renewal_shapes_oct11), collatz_amplitude_boundary_20261013c.md (Lemma B1, D(s))]
verification:
  - 04-computation/experiments/collatz_octave_memory_census_20261015.py (output 05-knowledge/results/collatz_octave_memory_census_20261015.json/.npz/.out)
  - 04-computation/experiments/collatz_octave_memory_analysis_20261015.py (output collatz_octave_memory_analysis_20261015.json/.out)
note: 05-knowledge/results/collatz_octave_memory_bias_energy_budgets_20261015.md
related: [HYP-9282, HYP-9283, THM-4629, THM-4630, THM-4581]
---

# THM-4658: two-octave memory, 4-adic energy descent, and the critical lines of the 3-adic transfer operator

## Setting

- `T(n) = n/2` or `(3n+1)/2`. `τ(n)` is the number of `T`-steps from `n` to 1, with `t(n) = (−1)^τ(n)` on the basin of 1 and `t(n) = 0` elsewhere; no convergence is assumed.
- Since the terminal cycle `1 → 2 → 1` has even length, `t(n) = −t(T(n))` for every `n >= 1`.
- For a set `A` of residues modulo `3^q`, put `C(A; N) = Σ_(n <= N, n ∈ A) t(n)` and `S(N) = C(Z; N)`.
- The wavelet channels are `H_(q,r)(N) = 3C(r mod 3^q; N) − C(r mod 3^(q−1); N)`, with `H = H_(1,2) = Σ_(n<=N) w_n t(n)` and `w_n = 3[n ≡ 2 (3)] − 1`.
- `E_X = Σ_(X <= n < 2X) H(n)²`. `D_X = Σ_(X <= N < 2X) Δ(N)²`, `C_X = Σ_(X <= N < 2X) H(N)Δ(N)`, where `Δ(N) := H(4N) − H(N)`.

## Statement

1. **Class-sum recursion (measure form of the wavelet identity (2)).** For every residue set `A` modulo `3^q` and every integer `N >= 0`,

   `C(A; N) = −C(A/2; ⌊N/2⌋) − C((3A+1)/2; 3⌈N/2⌉)`,

   where `A/2 = 2^(−1)A` modulo `3^q` and `(3A+1)/2 = 2^(−1)(3A+1)` modulo `3^(q+1)`. The second set lies in `2 + 3Z`.
2. **Two-octave memory.** For every integer `N >= 0`,

   `H(4N) = H(N) + H_(2,2)(3N) − H_(2,8)(6N)`,

   and, at one octave, `H_(1,2)(2N) = −H_(1,1)(N) − H_(2,8)(3N)`, `H_(1,1)(2N) = −H_(1,2)(N) − H_(2,2)(3N)`, `H_(1,0)(2N) = −H_(1,0)(N) − H_(2,5)(3N)`.
3. **Energy transport and 4-adic descent.** For every `X >= 1`,

   `E_(4X) = 4E_X + 8C_X + 4D_X + R_X`, with `|R_X| <= 24 sqrt(X E_(4X)) + 48X`,

   `sqrt(E_(4X)) <= 2 sqrt(E_X) + 2 sqrt(D_X) + 4 sqrt(3X)`.

   In normalized form, with `e_X = E_X/X²` and `δ_X = D_X/X²`: `sqrt(e_(4X)) <= (sqrt(e_X) + sqrt(δ_X))/2 + sqrt(3/X)`. Hence along every chain `X_j = 4^j X_0`,

   `sqrt(e_(X_J)) <= 2^(−J) sqrt(e_(X_0)) + Σ_(j<J) 2^(−(J−j)) (sqrt(δ_(X_j)) + 2 sqrt(3/X_j))`,

   so a uniform bound on the level-2 increment energy `δ` bounds the top-channel energy with the same constant in the limit.
4. **Main terms are transfer-operator eigen-distributions.** Let `ρ` be complex with `σ = Re ρ > 0`, and let `P_ρ` act on finitely additive signed distributions on `Z_3` by `(P_ρ β)(A) = −2^(−ρ) β(A/2) − (3/2)^ρ β((3A+1)/2)`.
   - **(a) Density level.** `P_ρ Haar = (−2^(−ρ) − (3/2)^ρ/3) Haar`, so Haar measure is an eigen-distribution with eigenvalue 1 exactly at the zeros of `D(ρ) = 1 + 2^(−ρ) + (3/2)^ρ/3`, the Mercator characteristic function of the amplitude note. If `S(N) = c N^ρ + o(N^σ)` with `c ≠ 0`, then `D(ρ) = 0`.
   - **(b) Excess level.** Suppose every cylinder `A` has `C(A; N) = Haar(A) S(N) + β(A) N^ρ + o(N^σ)`. Then `β(Z_3) = 0`, and `β` is fixed by the rank-one modification

     `(P̃_ρ β)(A) = −2^(−ρ) β(A/2) − (3/2)^ρ [β((3A+1)/2) − Haar(A) β(2 + 3Z_3)]`.

     The correction term carries the `H`-forcing: `C(2 + 3Z; N) − S(N)/3 = H(N)/3`.
5. **Critical lines.** On `L²(Z_3, Haar)`, the adjoint `P_ρ^* f(y) = −2^(−ρ) f(2y) − (3/2)^ρ 1_(y ≡ 2 (3)) f((2y−1)/3)` is a sum of two branches.
   - The first is `2^(−ρ)` times a unitary map.
   - The second is `(3/2)^ρ` times an operator `V` with `V^*V = I/3`.
   - Their norms are `2^(−σ)` and `3^(σ−1/2) 2^(−σ)` (`σ = Re ρ`), and their squared norms sum to `φ(σ) = 4^(−σ)(1 + 3^(2σ−1))`.
   - `φ` is strictly convex, `φ(σ) = 1` exactly at `σ = 1/2` and `σ = 1`, `φ < 1` on `(1/2, 1)` and `φ > 1` outside `[1/2, 1]`.
   - Re `σ = 1/2` is also the unique line on which the two branch norms are equal.

## Proof

**1.** The even sources `n = 2m ∈ A` contribute `t(2m) = −t(m)` with `m ∈ A/2` and `m <= ⌊N/2⌋`. The odd sources `n = 2m − 1 ∈ A` contribute `t(2m−1) = −t(3m−1)`.
- With `j = 3m − 1` one has `n = (2j−1)/3`, so `n ∈ A` iff `j ∈ (3A+1)/2` modulo `3^(q+1)`.
- `m <= ⌈N/2⌉` iff `j <= 3⌈N/2⌉ − 1`, and `3⌈N/2⌉ ≡ 0 (3)` lies outside `(3A+1)/2 ⊂ 2 + 3Z`.
- The case `n = 1` is the root identity `t(1) = −t(2)`. ∎

**2.** Statement 1 with `A = {2}`, `{1}`, `{0}` modulo 3 gives the one-octave relations, using `2^(−1) ≡ 2 (3)`, `2^(−1) ≡ 5 (9)`: `o(2) = 5·7 ≡ 8`, `o(1) = 5·4 ≡ 2`, `o(0) = 5 (9)`. The channel form is identity (2) of the wavelet note. Apply the `H_(1,2)` relation at `2N`, then the `H_(1,1)` relation at `N`. ∎

**3.** For `M = 4N + j`, `0 <= j <= 3`, write `H(M) = H(4N) + ε_(N,j)` with `ε_(N,0) = 0` and `ε_(N,j) = Σ_(i<=j) w_(4N+i) t(4N+i)`. Three consecutive `w` are a permutation of `(−1, −1, 2)`, so `|ε_(N,j)| <= 4`.
- *Identity.* Expand `Σ_N Σ_j (H(N) + Δ(N) + ε_(N,j))²`. The cross and square terms in `ε` give `R_X`, bounded by Cauchy–Schwarz with `Σ_N H(4N)² <= E_(4X)` and `Σ ε² <= 48X`.
- *Descent.* Minkowski's inequality in `ℓ²` of the `4X` pairs `(N, j)`: the main vector has norm `2‖H + Δ‖ <= 2(sqrt(E_X) + sqrt(D_X))` and `‖ε‖ <= 4 sqrt(3X)`.
- *Chain bound.* Iterate the normalized inequality. ∎

**4.** Use `⌊N/2⌋^ρ = 2^(−ρ)N^ρ + O(N^(σ − 1))` and `(3⌈N/2⌉)^ρ = (3/2)^ρ N^ρ + O(N^(σ − 1))`.
- *(a)* For Haar measure, `A/2` has the same measure as `A` and `(3A+1)/2` has a third of it. Statement 1 with `A = Z` reads `S(N) = −S(⌊N/2⌋) − C(2 + 3Z; 3⌈N/2⌉)`. If `S` has the stated main term, every class sum shares it with weight Haar at leading order, and comparing coefficients gives `c = −c 2^(−ρ) − c (3/2)^ρ/3`.
- *(b)* Put `Ẽ(A; N) = C(A; N) − Haar(A) S(N)`. Statement 1 for `A` and for `Z`, together with `Haar((3A+1)/2) = Haar(A)/3`, gives exactly

  `Ẽ(A; N) = −Ẽ(A/2; ⌊N/2⌋) − Ẽ((3A+1)/2; 3⌈N/2⌉) + Haar(A) Ẽ(2 + 3Z; 3⌈N/2⌉)`.

  Insert the main terms and compare coefficients of `N^ρ`. `β(Z_3) = 0` because `Ẽ(Z; N) = 0`. ∎

**5.** `f ↦ f(2·)` is measure-preserving and invertible on `Z_3`, hence unitary. For `V f(y) = 1_(y≡2) f((2y−1)/3)`, the substitution `z = (2y−1)/3` maps `2 + 3Z_3` bijectively onto `Z_3` with `dy = dz/3`, so `⟨Vf, Vg⟩ = ⟨f, g⟩/3`. Then `φ(σ) = 4^(−σ) + (9/4)^σ/3` is a sum of two exponentials, hence strictly convex, with `φ(1/2) = φ(1) = 1`. Equal norms means `3^(σ−1/2) = 1`. ∎

## Remarks

- **Readings.**
  - Statement 2 is the "mean-current" transport: the octave mean of `H` at `4X` inherits the octave mean at `X` with weight 1/2 in `sqrt(N)`-normalized units, plus the level-2 increment.
  - Statement 3 is a source-preserving energy rule in the sense of the break-clock note: it keeps the integer clock `N ↦ 4N, 3N, 6N` exactly.
  - Statement 5 explains why a white-sign model of `t` is self-consistent exactly at the `sqrt(N)` scale. The variance budget of the two refinement branches balances on `Re s = 1/2`, while the density and Mercator modes sit on the other balanced line, `Re s = 1`.
- **What it does not give.**
  - Statement 3 transfers the energy problem one level down; it does not close it. The census shows `Δ(N)` strongly correlated with `H(N)`: the cross term has `|κ|` up to 0.97, and the two-octave memory is almost exactly cancelled.
  - Statement 4 is conditional on main terms existing. HYP-9282 records evidence for nonzero real `sqrt(N)` main terms on many residue classes; by statement 4(b) their coefficients would be a fixed point of `P̃_(1/2)`.
- **Credit.**
  - The channel identity (2) and the wavelet hierarchy are the renewal_shapes_oct11 notes'.
  - The clock-mesh theorems M1/M2 and the terminal-controller theorems T1/T2 of clock_duality_oct15 (`collatz_clock_mesh_20261015a.md`, `collatz_terminal_mesh_20261015a.md`) pay the refinement energy; this theorem addresses the coarse/mean side.
