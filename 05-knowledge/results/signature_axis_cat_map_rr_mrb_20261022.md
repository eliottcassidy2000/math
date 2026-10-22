# The signature axis: definite and indefinite, Euclidean and Lorentzian, semisimple and logarithmic

**Session:** opus-2026-10-22-S55. **Status:** synthesis note. Its proved and numerical parts are registered as THM-4762, THM-4763, HYP-9372 and HYP-9373. Status labels are given per item.

## 0. Scope and concurrency

The owner asked for a merge of many sources:
* definite versus indefinite lattices, and Euclidean versus Lorentzian space;
* the Kazhdan–Lusztig correspondence and modular tensor categories (vertex-algebra and quantum-group module categories);
* the Fibonacci periods "15 in base 4, 24 in base 6, 120 in base 30";
* the Bailey mod 9, Dyson mod 27, Göllnitz–Gordon, Jackson–Slater, Rogers mod 14 and Rogers–Selberg identities;
* three MRB-constant PDFs, a Collatz abstract, and lists of "15"s.

A concurrent lane, **Session31** (`05-knowledge/results/session31_mirror_mrb_20261022/`), works the same PDFs through Lambert-W mirror series, rooted forests, CM data and a Collatz audit. This note avoids those routes. Its organising object is one matrix:
```
A = [[2,1],[1,1]] = t s t^-1 s^-1 = Q²      (the Arnold cat map; Q = Fibonacci matrix).
```
`A` is hyperbolic: a Lorentz boost of rapidity `2 log φ` in `SL2(R) ≅ Spin(2,1)`. It is also the commutator of the two elliptic (Euclidean) generators.

## 1. The Fibonacci periods (CHECKED; THM-4762 (1))

The last digit of `F_n` in base `b` has period `π(b)` (Pisano). The values are:
* `π(4) = 6`, **not 15**. In fact `π(N)` is even for every `N ≥ 3`, so the period 15 never occurs in any base.
  * A likely source of the 15 is the 4-bit maximal-length "Fibonacci" linear feedback shift register (`x⁴ + x³ + 1`), whose period is `2⁴ − 1 = 15`. That is a different object.
* `π(6) = 24` ✓.
* `π(30) = 120` ✓, and 30 is the least such base ✓ (`π(30) = lcm(3, 8, 20)`).

The coincidences `π(6) = 24 = |SL2(F_3)|` and `π(30) = 120 = |SL2(F_5)|` (binary icosahedral) have no map behind them: NUMEROLOGY.

In modular data, the cat map has order `π(N)/2` on the ground states of the `Z_N` toric code on a torus:
* order 3 for `Z_4`;
* order 12 for `Z_6`;
* order 60 for `Z_30`.

## 2. The cat map in modular data (THM-4762, HYP-9372)

For modular data `(S, T)`, the matrix `C = T S T⁻¹ S⁻¹` is a normalisation-free invariant.
* **What it does not see.** It is blind to the central charge, and for lattice (pointed) categories to the signature. Definite and Lorentzian lattices with the same discriminant form give the same `C`. Examples: `A_1`, `A_1 ⊕ U`, `A_1(−1)`; `A_2`, `A_2 ⊕ U`; `E_8`, `U`.
* **Pisano bounds.** By the congruence theorems, `ord C` divides a Pisano half-period.
* **Census.** 135 modular data were tested; 88 attain the bound.
* **The Ising exception.** The Ising factor alone is anomalous. There `C` acts through the spin structures, with period `3 = π(2)` (HYP-9372).
* **The golden Lorentzian lattice.** `(Z[φ], 2N)`, with Gram matrix `[[2,1],[1,−2]]`, has the cat map among its own automorphs. Its discriminant form `Z_5` gives `C` of order `10 = π(5)/2`, the same as the Fibonacci category.

## 3. The logarithmic side of Kazhdan–Lusztig (THM-4763, HYP-9373)

* **Semisimple side** (affine VOAs at positive integer level ≃ quantum groups at roots of unity). The cat map has finite order: Fibonacci mod `N`.
* **Logarithmic side** (triplet W(p) ≃ restricted quantum group, FGST). The torus amplitudes carry `R_(p+1) ⊕ C² ⊗ R_(p−1)`, so the cat map has spectral radius `φ²` and traces `Tr R_(p+1)(Aⁿ) + L_(2n) Tr R_(p−1)(Aⁿ)`, with `L` the Lucas numbers. This is verified for `p = 2..6`.
* HYP-9373 conjectures that the growth exponent is the logarithmic degree `K`.

Euclidean versus Lorentzian is literally semisimple versus logarithmic here. Only the logarithmic side retains the defining representation of `Spin(2,1)`.

## 4. The Rogers–Ramanujan family: definite products, indefinite partners

All 18 identities are verified exactly to `q^300` (`rr_family.out`). Products are written with `θ_(a,m) = (q^a, q^(m−a), q^m; q^m)_∞`.

The table uses two derived quantities:
* `c_eff` comes from part densities (`rr_asymptotics.out`). The coefficients grow like `exp(2π√(c_eff·n/6))`.
* The q → 1 ratios of companion functions are the sine ratios, i.e. quantum dimensions `[a]_q` at `q = exp(iπ/m)`.
* Characters `χ` in the table are equal up to an overall power of `q`.

| family | product side | `c_eff` | `6/c_eff` | q → 1 ratios | modular reading |
|---|---|---|---|---|---|
| Rogers–Ramanujan G, H | `θ_(a,5)/(q;q)` | 2/5 | **15** | `1/φ` | Lee–Yang `M(2,5)`, Galois conjugate of Fibonacci; cat map 10 |
| Rogers–Selberg A, B, C | `θ_(a,7)/(q²;q²) = χ^(M(2,7))_a · (q;q²)_∞` | 1/14 | 84 | heptagon `sin(kπ/7)/sin(3π/7)` | `M(2,7)` times a fermionic factor; cat map 4 |
| Rogers mod 14 A, B, C | `θ_(a/2,7)(q²)/(q;q) = χ^(M(2,7))(q²) · (−q;q)_∞` | 11/14 | 84/11 | heptagon | `M(2,7)` at `q²` times the free-fermion factor |
| Göllnitz–Gordon 1, 2 | `1/((q,q⁴,q⁷;q⁸)_∞)` | 3/8 | 16 | `√2 − 1 = tan(π/8)` | octagon; the silver ratio `1 + √2` plays the role of `φ` |
| Bailey mod 9 A, B, C | `θ_(a,9)/(q³;q³)` | 0 | — | nonagon `sin(kπ/9)/sin(4π/9)` | weight-0 modular units (no growth; signed coefficients) |
| Dyson mod 27 A–D | `θ_(b,9)(q³)/(q;q) = Bailey_b(q³) · (q⁹;q⁹)/(q;q)` | 8/9 | 27/4 | nonagon and `sin(π/3)` | the `M(2,27)` characters with `3 \| a`; Dyson B is the lift of the trivial `b = 3` |
| Jackson–Slater | `(q,q⁷,q⁸;q⁸)(q⁶,q¹⁰;q¹⁶)/(q;q)` | 1/2 | 12 | — | **Ising vacuum character** `q^(1/48) χ_0` |

Notes:
* The "15" in the growth `exp(2π√(n/15))` of the Rogers–Ramanujan partitions is exactly `6/c_eff` of Lee–Yang.
* The Jackson–Slater identity is the Ising vacuum character in fermionic form: the even part of Euler's `(−q^(1/2); q)_∞`. This ties it to the anomalous Ising row of §2 and to the Ising rings of §5.

**The indefinite (Lorentzian) partners.**
* Each sum side has a "sign-twisted" companion that is mock modular: a Hecke-type theta series over a cone of a signature-(1,1) lattice.
* For Rogers–Ramanujan it is the fifth-order mock theta function `f_0(q) = Σ q^(n²)/(−q;q)_n`. It satisfies
  ```
  (q;q)_∞ f_0(q) = Σ_(n≥0, |j|≤n) (−1)^j q^(n(5n+1)/2 − j²)(1 − q^(4n+2)),
  ```
  VERIFIED to `q^300`, against Andrews (TAMS 293, 1986).
* The quadratic form `5n²/2 − j²` is indefinite; the cone `|j| ≤ n` makes the sum converge.
* The mock theta conjectures (Hickerson 1988) express `f_0` through Rogers–Ramanujan functions plus Appell–Lerch sums.
* The seventh-order analogues (Andrews 1986) accompany the mod-7 family. They are cited, not reproduced: a simple ansatz was tested and failed.

So the definite product (a one-dimensional positive-definite Jacobi theta) and the indefinite Hecke sum are the Euclidean and Lorentzian faces of the same identity family.

## 5. The supplied MRB documents (data; checked, `mrb_checks.out`)

* **Finite irrationality theorem ("alternating roots").** The trace argument is correct. Each `n^(1/n)`, `n ≥ 2`, has trace 0 in `L_N`, so a rational `S_N` would equal 0 or −1. The bounds hold: `S_5 = 0.0064479 > 1/200`, the minimum odd partial sum for `N ≤ 59` is `S_5`, and the minimum even partial sum is 0.2239. As the document itself says, this does not reach the limit.
* **Prime-defect theorem ("prime persistence").**
  * The defect `D_p` was computed numerically for `p = 3, 5, 7, 11, 101, 211`. `|D_p + 4πi log2/(p−1)|` lies within the stated bound `(30 log p + 661)/p²` in every case.
  * `D_p` has a nonzero real part coming from the remainder (`0.0067` at `p = 101`).
  * The paper's own caveat stands: a Galois conjugate need not be continuous in the real embedding.
  * Typed ANALOGY with modular data: Galois symmetry of modular tensor categories (Coste–Gannon) acts through `Gal(Q(ζ_N)/Q)` on S-matrix entries. The MRB defect uses Kummer automorphisms that fix every root of unity and move radicals. These are complementary directions in `Q(μ_∞, radicals)`. No theorem is claimed.
* **Ising rings.**
  * The identity `C_MRB = Σ[(g_(n,−)/2)^(1/n) − (g_(n,+)/2)^(1/n)]` is exact.
  * The periodic and antiperiodic rings are two spin structures on the circle. Under `SL2(Z)`, the three even torus spin structures are cyclically permuted by the cat map with period `π(2) = 3`. That is the same period the Ising category's `C` has (HYP-9372 (c)).
  * At the special coupling `e^(−2K) = √5`, the transfer-matrix eigenvalue ratio is `(1−√5)/(1+√5) = −1/φ²`, and the ring partition functions become `2^n L_n` and `2^n √5 F_n`: Lucas and Fibonacci.

## 6. The Collatz abstract (data; checked)

The abstract claims a complete proof. Its parametrisation `N = a·2^n − 2` with `a` odd is the trailing-ones (Mersenne-line) coordinate `N/2 = a·2^(n−1) − 1`; see the repo's Mersenne-line work, S21.

Its key mechanism is that the "15 mod 16" run turns into a "huge factor `2^k`" at the vertical maximum. This is false in general:
* starting from `2^k − 1`, `k` rising steps reach `3^k − 1`;
* by lifting the exponent, `v_2(3^k − 1) = 1` for odd `k` and `2 + v_2(k)` for even `k` (`mrb_checks.out`, `k = 2..12`).

So `15 = 2⁴ − 1 → 80 = 2⁴·5` is the unique coincidence `v_2 = k`; already `31 → 242 = 2·121`. The abstract's "100%" blocking of divergence therefore does not follow. Session31 audits a supplied Collatz contraction proof separately.

## 7. The "15"s

* No base has Fibonacci period 15 (§1).
* 15 = `6/c_eff` of Lee–Yang (§4).
* The cat map's projective order is 15 for `SU(2)_3`, `M(3,5)`, `M(4,5)`, `M(4,11)` and `SU(2)_9` (THM-4762 census).
* The supersingular primes (Ogg), Mazur's 15 torsion groups, the 15 pentagon tilings and the Conway–Schneeberger theorem have no map to the objects here. They were not pursued: NUMEROLOGY.

## 8. Next questions

* Is faithfulness modulo the centre (HYP-9372 (b)) a theorem for every modular data without an Ising-type factor?
* Do the logarithmic `W(p,p')` models realise `K = 2` and spectral radius `φ⁴` (HYP-9373)?
* Is there a mock-modular (indefinite) analogue of THM-4763? Completions of mock theta functions carry non-semisimple `SL2(Z)` actions through Eichler integrals. Do hyperbolic elements then grow like a power of the dilatation?
* Is there a twisted (Kummer) analogue of the Galois action on modular data that sees the MRB defect?
