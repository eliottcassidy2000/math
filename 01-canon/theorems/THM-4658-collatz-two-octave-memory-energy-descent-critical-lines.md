---
id: THM-4658
title: "Two-octave memory and 4-adic energy descent for the Collatz mod-3 discrepancy: H(4N) = H(N) + H_(2,2)(3N) - H_(2,8)(6N) exactly, so E_(4X) = 4E_X + 8C_X + 4D_X + R_X with |R_X| <= 24 sqrt(X E_(4X)) + 48X and sqrt(E_(4X)) <= 2 sqrt(E_X) + 2 sqrt(D_X) + 4 sqrt(3X); the class-sum recursion C(A;N) = -C(A/2; floor(N/2)) - C((3A+1)/2; 3 ceil(N/2)) makes Haar measure an eigen-distribution of a 3-adic transfer operator exactly at the zeros of D, puts power-law excess main terms under a rank-one-modified fixed-point equation, gives every halving eigenfunction (Dirichlet characters mod 3^q) an exact one-octave inheritance C_psi(N) = -psi(2) C_psi(floor(N/2)) + O_psi(N), and on densities of any finite 3-adic depth has spectrum {1 - D(s)} together with -2^(-s) times roots of unity (no finite-depth excess mode off Re s = 0); the two operator branches have squared L^2(Z_3)-norms summing to 4^(-sigma)(1 + 3^(2 sigma - 1)), equal to 1 exactly at Re s = 1/2 and Re s = 1"
status: >
  PROVED (elementary consequences of the exact refinement identity (2) of collatz_wavelet_energy_20261013d.md, here in
  measure form). AUDITED 2026-10-15 (statements 1-5: first S37 referee, 04-computation/experiments/s37_audit_20261015/
  REPORT.md, statement 4(a) repaired by clock_duality_oct15; statements 6-7: second S37 referee,
  04-computation/experiments/s37_audit2_20261015/REPORT.md, spectra, multiplicities and proofs confirmed). FINITE-EXACT control: the energy transport of statement 3 checked
  on every octave 2^10..2^30 of an exact census through 2^34 (identity exact; the remainder R_X is visible at the 4th
  decimal for k = 10..13, 16 and below 1e-6 from k = 17); census reproduces S(2^34) = 24953018 and H(2^34) = 26869 of
  the repository's independent 2^34 census. No bound on the actual Collatz discrepancy is claimed.
source: opus-2026-10-15-S37
depends_on: [collatz_wavelet_energy_20261013d.md (identity (2), renewal_shapes_oct11), collatz_amplitude_boundary_20261013c.md (Lemma B1, D(s))]
verification:
  - 04-computation/experiments/collatz_octave_memory_census_20261015.py (output 05-knowledge/results/collatz_octave_memory_census_20261015.json/.npz/.out)
  - 04-computation/experiments/collatz_octave_memory_analysis_20261015.py (output collatz_octave_memory_analysis_20261015.json/.out)
  - 04-computation/experiments/collatz_twisted_split_memory_20261015.py (statement 6 for the twisted splits D_j; output collatz_twisted_split_memory_20261015.json)
  - 04-computation/experiments/s37_audit_20261015/a1_thm4658_identities.py (independent referee checks of statements 1-5)
  - 04-computation/experiments/collatz_excess_operator_resonance_scan_20261015.py (numerical control of statement 7, Q <= 4)
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
   - **(a) Density level.** `P_ρ Haar = (−2^(−ρ) − (3/2)^ρ/3) Haar`, so Haar measure is an eigen-distribution with eigenvalue 1 exactly at the zeros of `D(ρ) = 1 + 2^(−ρ) + (3/2)^ρ/3`, the Mercator characteristic function of the amplitude note. If `S(N) = c N^ρ + o(N^σ)` with `c ≠ 0` **and** `H(N) = o(N^σ)`, then `D(ρ) = 0`. The latter hypothesis is implied by Haar-leading class asymptotics; it is not inferred from the total-sum asymptotic alone.
   - **(b) Excess level.** Suppose every cylinder `A` has `C(A; N) = Haar(A) S(N) + β(A) N^ρ + o(N^σ)`. Then `β(Z_3) = 0`, and `β` is fixed by the rank-one modification

     `(P̃_ρ β)(A) = −2^(−ρ) β(A/2) − (3/2)^ρ [β((3A+1)/2) − Haar(A) β(2 + 3Z_3)]`.

     The correction term carries the `H`-forcing: `C(2 + 3Z; N) − S(N)/3 = H(N)/3`.
5. **Critical lines.** On `L²(Z_3, Haar)`, the adjoint `P_ρ^* f(y) = −2^(−ρ) f(2y) − (3/2)^ρ 1_(y ≡ 2 (3)) f((2y−1)/3)` is a sum of two branches.
   - The first is `2^(−ρ)` times a unitary map.
   - The second is `(3/2)^ρ` times an operator `V` with `V^*V = I/3`.
   - Their norms are `2^(−σ)` and `3^(σ−1/2) 2^(−σ)` (`σ = Re ρ`), and their squared norms sum to `φ(σ) = 4^(−σ)(1 + 3^(2σ−1))`.
   - `φ` is strictly convex, `φ(σ) = 1` exactly at `σ = 1/2` and `σ = 1`, `φ < 1` on `(1/2, 1)` and `φ > 1` outside `[1/2, 1]`.
   - Re `σ = 1/2` is also the unique line on which the two branch norms are equal.
6. **Character inheritance (halving sectors).** Let `φ` be a function on `Z` that depends only on `n` modulo `3^q` and satisfies `φ(2n) = λ φ(n)`; for example a Dirichlet character `ψ` modulo `3^q` with `λ = ψ(2)`, or the twisted split weight `φ_j(n) = [3^j ∥ n] χ(n/3^j)` with `χ` the nontrivial character modulo 3 and `λ = −1`. Put `C_φ(N) = Σ_(n<=N) φ(n) t(n)`. Then for every integer `N >= 0`,

   `C_φ(N) = −λ C_φ(⌊N/2⌋) + O_φ(N)`,  `O_φ(N) = Σ_(n<=N, n odd) φ(n) t(n) = −Σ_(m <= 3⌈N/2⌉, m ≡ 2 (3)) φ((2m−1)/3) t(m)`.

   In particular the twisted splits `D_j = C_(φ_j)` have positive one-octave memory, `D_j(N) = D_j(⌊N/2⌋) + O_j(N)`, while Haar-type sums (`λ = 1`, e.g. `S` and the untwisted multiples of `3^j`) alternate. Since `2` generates `(Z/3^q)^×`, the eigenfunctions of `n ↦ 2n` on the units modulo `3^q` are exactly the multiples of Dirichlet characters, and `λ = −1` occurs only for the character of order 2.
7. **Finite-depth spectrum: no finite-depth excess modes off the imaginary axis.** Fix `Q >= 1` and identify densities constant on the residue classes modulo `3^Q` with functions `f` on `Z/3^Q`. On them `P_ρ` (density form, `β = f dHaar`) acts by

   `(P_ρ f)(a) = −2^(−ρ) f(2^(−1)a) − (1/3)(3/2)^ρ f(2^(−1)(3a+1))`,

   and `P̃_ρ = P_ρ + (3/2)^ρ ⟨f, 1_(2+3Z)⟩` (with `⟨f, g⟩ = 3^(−Q) Σ_a f(a) g(a)`) preserves `Σ_a f(a) = 0`. Then:
   - the spectrum of `P_ρ` on these functions is `{1 − D(ρ)} ∪ {−2^(−ρ) ζ : ζ^(2·3^k) = 1, 0 <= k <= Q − 1}`;
   - the spectrum of `P̃_ρ` on the mean-zero ones is `{−2^(−ρ) ζ : ζ^(2·3^k) = 1, 0 <= k <= Q − 1}`.

   So the finite-depth transfer operator has eigenvalue 1 off the imaginary axis exactly at the zeros of `D`, and the excess operator never does: `P̃_ρ f = f` with `Σ f = 0`, `f ≠ 0`, forces `Re ρ = 0`. Consequently a power-law excess main term `β(A) N^ρ` as in 4(b) with `Re ρ ≠ 0` cannot have a density constant on the classes of any fixed modulus `3^Q`: its 3-adic Fourier transform must reach arbitrarily fine frequencies.

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
- *(a)* For Haar measure, `A/2` has the same measure as `A` and `(3A+1)/2` has a third of it. Statement 1 with `A = Z` reads `S(N) = −S(⌊N/2⌋) − C(2 + 3Z; 3⌈N/2⌉)`. Now use the exact identity `C(2 + 3Z; M) = (S(M) + H(M))/3`. Under both displayed hypotheses, comparing coefficients gives `c = −c 2^(−ρ) − c (3/2)^ρ/3`. Without the small-H assumption, an order-N^σ excess term remains and this comparison does not give `D(ρ)=0`.
- *(b)* Put `Ẽ(A; N) = C(A; N) − Haar(A) S(N)`. Statement 1 for `A` and for `Z`, together with `Haar((3A+1)/2) = Haar(A)/3`, gives exactly

  `Ẽ(A; N) = −Ẽ(A/2; ⌊N/2⌋) − Ẽ((3A+1)/2; 3⌈N/2⌉) + Haar(A) Ẽ(2 + 3Z; 3⌈N/2⌉)`.

  Insert the main terms and compare coefficients of `N^ρ`. `β(Z_3) = 0` because `Ẽ(Z; N) = 0`. ∎

**5.** `f ↦ f(2·)` is measure-preserving and invertible on `Z_3`, hence unitary. For `V f(y) = 1_(y≡2) f((2y−1)/3)`, the substitution `z = (2y−1)/3` maps `2 + 3Z_3` bijectively onto `Z_3` with `dy = dz/3`, so `⟨Vf, Vg⟩ = ⟨f, g⟩/3`. Then `φ(σ) = 4^(−σ) + (9/4)^σ/3` is a sum of two exponentials, hence strictly convex, with `φ(1/2) = φ(1) = 1`. Equal norms means `3^(σ−1/2) = 1`. ∎

**6.** The even sources `n = 2m <= N` give `φ(2m) t(2m) = −λ φ(m) t(m)` with `m <= ⌊N/2⌋`. The odd sources form `O_φ(N)`; the substitution `n = 2m − 1 ↦ j = 3m − 1` of statement 1 and `t(n) = −t(T(n))` give the second expression. For a character, `φ(2n) = ψ(2)φ(n)`. For `φ_j`, `3^j ∥ 2n` iff `3^j ∥ n`, and `χ(2n/3^j) = χ(2)χ(n/3^j) = −χ(n/3^j)`. A function on the cyclic group `(Z/3^q)^× = ⟨2⟩` with `φ(2n) = λφ(n)` is a character times a constant. ∎

**7.** The density form follows from 4(b) as in statement 4: `β(A/2) = 3^(−q) f(a/2)` and `β((3A+1)/2) = 3^(−q−1) f((3a+1)/2)` for a level-`q >= Q` cylinder `A = a + 3^q Z_3`. Expand `f = Σ_m f̂(m) e_m`, `e_m(a) = exp(2πi ma/3^Q)`. Then `f(2^(−1)·)` has coefficient `f̂(2m)` at `e_m`, and `f(2^(−1)(3·+1)) = Σ_k f̂(k) e_k(2^(−1)) e_(3·2^(−1) k)`, whose coefficient at `e_m` collects the `k` with `3·2^(−1)k ≡ m`, all of which have `v_3(k) = v_3(m) − 1` (none if `v_3(m) = 0`); the rank-one term only touches `e_0`. Order the frequencies by `v_3(m) = 0, 1, …, Q − 1` and put `m = 0` last. In this order both operators are block lower triangular. The diagonal block on `{v_3(m) = j}` is `f̂(m) ↦ −2^(−ρ) f̂(2m)`: `−2^(−ρ)` times the permutation `u ↦ 2u` of `(Z/3^(Q−j))^×`, a single cycle of length `2·3^(Q−j−1)` because `2` is a primitive root modulo every power of 3. Its eigenvalues are `−2^(−ρ)ζ` with `ζ^(2·3^(Q−j−1)) = 1`. The diagonal entry at `m = 0` is `−2^(−ρ) − (3/2)^ρ/3 = 1 − D(ρ)` for `P_ρ`; for `P̃_ρ` the frequency `0` is absent on mean-zero functions. Eigenvalue 1 of `−2^(−ρ)ζ` needs `|2^(−ρ)| = 1`. ∎

## Remarks

- **Readings.**
  - Statement 2 is the "mean-current" transport: the octave mean of `H` at `4X` inherits the octave mean at `X` with weight 1/2 in `sqrt(N)`-normalized units, plus the level-2 increment. Since `Δ(N) = H(4N) − H(N)` by definition, the empirical "cancelled memory" (two-octave regression slope 0.006) says only that the two-octave autocorrelation of the normalized profile is near zero. The new content is that `Δ` is an exact level-2 expression; clock_duality_oct15 writes the same increment as the annular source row `Δ(N) = Σ_(N<n<=4N) w_n t(n)`.
  - Statement 3 is a source-preserving energy rule in the sense of the break-clock note: it keeps the integer clock `N ↦ 4N, 3N, 6N` exactly. It is the two-octave analogue of the one-octave energy rule (3) of the wavelet note (renewal_shapes_oct11).
  - Statement 5 is the variance recursion for independent white-sign branches: if the branch contributions were uncorrelated, `sqrt(N)` would be the self-consistent scale. It is not an operator-norm statement. The exact level-6 norms of `P*_ρ` are 1.398, 1.362, 1.352 at `ρ = 1/2, 3/4, 1` (about 1.34 on mean-zero functions; S37 referee, `a9_transfer_operator_norms.py`), so "critical line" is a heuristic label. The Mercator zeros sit just below `Re s = 1` (`s_0 ≈ 0.99821 + 240.248i`).
  - Statement 7 explains the structure: in 3-adic Fourier coordinates the odd branch only moves information to coarser frequencies, so at any finite depth it contributes only strictly lower-triangular blocks, and only the density mode (frequency 0) can resonate, at the zeros of `D`. The eigenvalue-1 points on `Re ρ = 0` are the `ρ` with `2^(−ρ) = −ζ̄` (e.g. `Im ρ` = 4.532, 9.065, 13.597 for `Q = 1`). An argument-principle scan finds no eigenvalue-1 point in `[0.05, 1.5] × (0, 300]` for `Q <= 4`, and the predicted spectra match to `10^(−14)` (`collatz_excess_operator_resonance_scan_20261015.py`). The multiplicity of `−2^(−ρ)ζ` is the number of `k < Q` with `ζ^(2·3^k) = 1` (no Jordan blocks at `ζ = ±1`; second referee). For non-real `ρ`, a main term of a real class sum is a conjugate pair `β N^ρ + β̄ N^(ρ̄)`, and 4(b) and 7 apply to each member. The infinite-depth operator is not covered: there the triangular induction has no finest level, and eigen-distributions with unbounded frequency depth are not excluded.
  - **Infinite depth (remark with proof).** Statement 7 cannot be extended to all finitely additive distributions. In the Fourier coordinates of its proof, a mean-zero distribution is a sequence of blocks `f̂_d` on the frequencies of depth `d >= 1` (`ξ = u/3^d`, `u` a unit modulo `3^d`), and the excess eigen-equation reads `f̂_d = L_d f̂_(d+1)` with `L_d = −(1/3)(3/2)^ρ (I + 2^(−ρ)Π_d)^(−1) In_d`. Here `Π_d` is the permutation `u ↦ 2u`, and `In_d` sums each fibre of three depth-`(d+1)` frequencies over `ξ` with unimodular phases, so `In_d In_d^* = 3I`. Hence, for `Re ρ ≠ 0`, every `L_d` is onto, and its nonzero singular values are `(3/2)^σ 3^(−1/2) |1 + 2^(−ρ)ζ|^(−1)`, `ζ` running over the `2·3^(d−1)`-th roots of unity (`σ = Re ρ`).
    - *Existence of distributions.* Choosing `f̂_1 ≠ 0` and preimages successively gives an infinite-dimensional space of nonzero excess eigen-distributions at every `ρ` with `Re ρ ≠ 0`, all of unbounded depth (the second referee built them upward level by level for every `ρ` tried, including `s_0`, with all cylinder equations at levels 0-6 verified to 5e-15). So statement 4(b) alone never excludes an excess main term; an exclusion needs regularity of `β`.
    - *Square-integrable densities.* If every singular value exceeds 1, that is `(3/2)^σ > sqrt(3)(1 + 2^(−σ))`, i.e. `σ > σ_2 = 1.92963...`, then minimal-norm preimages shrink geometrically and the construction yields an `L²(Z_3)` eigen-density. That range is irrelevant to main terms, since `|C(A; N)| <= N` forces `Re ρ <= 1`. For `0 < σ <= 1` the norms do not decide: for real `ρ` the largest singular value is `(3/2)^σ/(sqrt(3)(1 − 2^(−σ))) > 1` and the smallest is below 1. Whether a square-integrable excess eigen-density exists with `0 < Re ρ <= 1` is OPEN. In the other direction, clock_duality_oct15 proved that the excess operator contracts in the weighted Wiener norm `Σ_ξ |f̂(ξ)| R^(q(ξ))` (`q` the conductor exponent) whenever `R > 3^(σ−1)/(2^σ − 1)`, so no excess eigen-density of finite weighted norm exists there (Proposition W1 of `collatz_conductor_wiener_20261015d.md`; at `σ = 1/2`, `R = sqrt 2` works).
  - Statement 6 is the halving half of the recursion in character coordinates. Averaged over octaves, any exact relation `F(N) = a F(⌊N/2⌋) + G(N)` gives `(1 − a 2^(−1/2)) ⟨F/sqrt(N)⟩ = ⟨G/sqrt(N)⟩` up to end effects of relative size about one over the number of octaves, for every sequence whatsoever. So the ratio `⟨D_j⟩/⟨O_j⟩ ≈ 2 + sqrt(2)` (measured 3.12, 3.48, 3.36 for `j = 0, 1, 2`) and the octave-averaged fixed-point "consistency" of statement 4(b) (residuals 6-12% of the bias size, below their standard errors, levels 1-4) are code checks, not evidence of main terms. The same inheritance holds verbatim for random signs with `t(2m) = −t(m)` (the halving null of the S37 note). This is why channel octave means carry strong octave-to-octave memory (pooled lag-one correlation about 0.65 at every level).
- **What it does not give.**
  - Correction (clock_duality_oct15, 2026-10-15): statement 4(a)'s pre-audit scalar-to-Haar inference omitted `H(N)=o(N^σ)`. The amended implication retains that assumption; the Haar eigenvalue computation and statement 4(b) remain unchanged. The correction mechanism is recorded in `01-canon/MISTAKES.md`.
  - The branch squared-norm sum in statement 5 is not the norm of their sum. A mean-zero depth-two cylinder function gives norm expansion at real `ρ=3/4`, despite `φ(3/4)<1`; see [the exact incoming-work audit](../../05-knowledge/results/collatz_incoming_audit_20261015b.md). No full-operator contraction is asserted.
  - Statement 3 transfers the energy problem one level down; it does not close it. The census shows `Δ(N)` strongly aligned with `H(N)`: the uncentred cosine `κ = C_X/sqrt(E_X D_X)` ranges over −0.97..+0.91 on 2^14..2^30 (the wavelet note reported −0.949..+0.936 for its one-octave cross terms).
  - Statement 4 is conditional on main terms of order exactly `N^σ` with `o(N^σ)` errors. HYP-9282 proposed real `sqrt(N)` main terms on many residue classes; it failed its pre-registered out-of-sample test (S37 census to 2^37) and is REFUTED. A model "main term plus mean-zero fluctuation of the same order" would in any case need a logarithmic-mean version of 4(b), which is not stated here.
- **Credit.**
  - The channel identity (2) and the wavelet hierarchy are the renewal_shapes_oct11 notes'.
  - The clock-mesh theorems M1/M2 of clock_duality_oct15 (`collatz_clock_mesh_20261015a.md`) pay the refinement energy. Its terminal-controller theorems T1/T2 (`collatz_terminal_mesh_20261015a.md`) are coarse-mean certificates. This theorem addresses the coarse/mean side through the exact transport.
  - The repair of statement 4(a) and the branch-cross-term hostile are clock_duality_oct15's (`collatz_incoming_audit_20261015b.md`); the independent S37 referee found the same gap.
