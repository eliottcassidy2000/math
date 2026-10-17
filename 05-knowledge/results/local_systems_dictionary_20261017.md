# One mechanism behind the fast computations: local systems and their moments, across tournaments, unit distances, Collatz and Keller maps

2026-10-17, opus-2026-10-17-S39. **PROVED scoped pieces + FINITE-EXACT checks + CITED inputs + ANALOGY for the
cross-row reading** (each item typed). Canon: [THM-4663](../../01-canon/theorems/THM-4663-polygon-ladder-of-the-finite-unit-distance-graphs-closed-walks-are-kloosterman-moments-pentagons-level-15-hexagons-level-6.md).
Hypothesis: [HYP-9286](../hypotheses/HYP-9286-local-system-dictionary-for-tournaments-unit-distances-collatz-and-keller-maps.md).
No Collatz, planar-Jacobian or Hadwiger–Nelson conclusion follows.

## 1. The question

The owner asked for more fundamental relations between the planar Jacobian conjecture, Collatz and Hadwiger–Nelson,
blended with the mathematics that made the repository's tournament and sequence computations fast. Looking back at
those computations — `H(P_31)`, the Paley wheel and run ladder (S34–S35), the Collatz censuses and the twisted
transfer law (S26–S38) — the common mechanism is this: **each object carries a natural local system on a parameter
space over a finite field, the quantities computed are the moments of that local system, closed forms exist exactly
where the moments are motives of low weight, and the "cores" of each ladder are the places where a new motive
enters.** Below are five rows, one per setting, with what is proved in each.

## 2. The five rows

| row | object | local system | rank / weight / monodromy | moments = repository counts | arithmetic |
|---|---|---|---|---|---|
| T | Paley tournament `P_q` | quadratic character; Legendre family (Jacobi sums, Greene) | 2 / 1 / `SL_2` | wheel and run sums (THM-4652, HYP-9278) | level 8 (`W_4` = `η(2z)^4 η(4z)^4`, `W_6`); first core `K_6 - P_6`, level 32 (HYP-9276) |
| H | unit-distance graph `G_q` | Kloosterman sheaf `Kl_2` (closed chains), four-bar surface (open 3-chains) | 2 / 1 / `SL_2`; elliptic surface of `Γ_1(6)` | polygon and theta counts (THM-4663) | levels 15, 6, 525, 6; first core `θ(2,2,3)`: 30a1 |
| R | the real plane | Bessel functions (`J_0` = Fourier transform of the circle) | de Rham realization of row H | walk densities, Bessel moments | `L`-values of the same forms (THM-4663 (7)) |
| C | Collatz modulo a prime `p` | the affine pair `x/2, (3x+1)/2` on `F_p` (Schreier digraph of `F_p ⋊ <2,3>`) | permutation (weight 0) | traces = periodic points; degenerations at `2^k ≡ 3^j (mod p)` | none: finite-place shadows only (§4) |
| J | a Keller map `F` | the defect `F_* Q_l / Q_l` off the Jelonek set | `d - 1` / 0 / finite (`S_(d+1)` for the sweep family) | fibre statistics; Fourier energy = collisions | value-set constants `c_(d+1) → 1 - 1/e`; spectrum on the dual of the Jelonek curve at infinity (THM-4664) (§5) |

## 3. Rows T, H, R: two elliptic modular surfaces, two second moments (PROVED / CITED / FINITE-EXACT)

The cleanest relation found this session is a precise parallel between rows T and H.

* **Tournaments.** The Paley 4-wheel is governed by `η(2z)^4 η(4z)^4` (S34). In moment form this is the
  `χ(λ)`-twisted second moment of the Legendre family `y^2 = x(x-1)(x-λ)`:
  `Σ_(λ ≠ 0,1) χ(λ) a_λ^2 = -1 - p - c_p`, already in the repo as [HYP-9274](../hypotheses/HYP-9274-paley-rate-three-halves-wheel-w4-apery-form.md) (B)
  (Ahlgren–Ono Lemma 2.2 + Theorem 6, every odd prime); re-checked here at every prime `5 ≤ p < 300`
  ([script](../../04-computation/experiments/legendre_fourbar_second_moments_20261017.py)).
  The untwisted second moment is not of this form.
* **Unit distances.** The closed equilateral hexagons of `G_p` are governed by `(η1 η2 η3 η6)^2` (THM-4663 (3)), and the
  same form is the second moment of the four-bar family `E'_c` (THM-4663 (9), the elliptic modular surface of
  `Γ_1(6)`): `Σ_(c ≠ 0,1,9) a_c^2 = (p+1)(p-4) - b_p` (THM-4663 (10), PROVED by Eichler–Shimura for `X_1(6)`, whose
  weight-4 cusp space is spanned by `b`; exact at every prime `11 ≤ p < 1000`). Its third moment is a weight-5 level-6
  trace (found by the audit).

So in both rows a rank-2 object of weight one (the quadratic character; the unit circle) has, behind its pattern
sums, a rational elliptic modular surface (Legendre; four-bar `Γ_1(6)`), and the first weight-4 form of the ladder is
that surface's (twisted) second moment, the motive of a rigid Calabi–Yau fibre product (Schoen's construction;
Eichler–Shimura for the modular surface, or Gouvêa–Yui). The level-15 pentagon form is produced by the same four-bar
surface: its `η(c(c-4))`-twisted first moment is `-a_p(f_(3,15))` up to Legendre terms (THM-4663 (10)); the base change
along the codegree conic is a K3 whose transcendental Frobenius trace is `+a_p(f_(3,15))`, CM by `Q(sqrt(-15))`. The seventh moment adds level 525. Row R is the same arithmetic at the
real place: the five-step planar walk density `p_5'(0) = 6 L(f_15, 1)/π^2` (NUMERICAL, 25 digits), and the three-step
density is a period of the four-bar curve, `p_3(r) = (2/π^2) r Ω(r^2)` (PROVED, by convolution; 22 digits) — the same curve whose
`F_p`-points count the finite 3-chains (THM-4663 (9)).

## 4. Row C: the Collatz pair modulo `p` (PROVED; [script](../../04-computation/experiments/collatz_digraph_mod_p_trace_formula_20261017.py))

For a prime `p ∤ 6` let `A_p` be the adjacency matrix of the digraph on `F_p` with arcs `x → x/2`, `x → (3x+1)/2`.
The group generated is `F_p ⋊ <2,3>` (this is in Codex's [affine guarded lifts note](collatz_affine_guarded_lifts_20261004.md),
which uses Euler circuits of its Cayley digraph); new here:

* **Trace formula.** `tr(A_p^k) = 2^k + Σ_(w : 3^j ≡ 2^k (mod p)) (p [p | c_w] - 1)`, where `w` runs over parity words of
  length `k` with `j` odd steps and `g_w(x) = (3^j x + c_w)/2^k`. Every word has the rational 2-adic periodic point
  `x_w = c_w/(2^k - 3^j)` (Böhm–Sontacchi; Lagarias); the closed walks count its reductions modulo `p`, and the words
  whose denominator vanishes modulo `p` contribute `p` or `0` according as the numerator does. Reading (ANALOGY): the
  finite-place shadow of THM-4659's parity-matched near-returns `2^k ≈ 3^j` at the real place. Verified for `k ≤ 14`,
  twelve primes (and by the audit for every prime `5 ≤ p < 62`, `k ≤ 20`).
* **Equidistribution of cycle numerators.** The digraph is strongly connected with a loop, so Perron–Frobenius gives
  a simple eigenvalue 2 and `|tr A_p^k - 2^k| ≤ (p - 1) ρ_p^k` with `ρ_p < 2`. Hence, along the lengths `k` with
  `2^k ∈ <3> (mod p)` (otherwise there are no such words, `B_k = 0`: every odd `k` at `p = 11`, `k ≢ 0 (mod 4)` at
  `p = 13`), among the words with `p | 2^k - 3^j` the fraction with `p | c_w` is `1/p + O(p ρ_p^k / B_k)`: the prime `p` cancels
  from the rational cycle point with probability tending to `1/p`, at the exponential rate `ρ_p/2`
  (`ρ_5 = 1.618`, `ρ_7 = 1.428`, `ρ_11 = 1.338`; fractions at `k = 20`: 0.2116, 0.1446, 0.0903 against 1/5, 1/7, 1/11).
* **Arborescences.** The digraph is 2-in 2-out, so by BEST the number of Eulerian circuits (the binary necklaces of
  length `2p` driving the affine walk through every residue-parity pair once) is the arborescence count `τ_p`:
  2, 6, 96, 288, 3582, 5994, ... for `p = 5, 7, 11, 13, 17, 19`. There is no visible closed form (large prime factors);
  the characteristic polynomial on the nonconstant functions often has linear factors (eigenvalues 0, ±1) and at `p = 5`
  the non-cyclotomic factor `z^2 + z - 1`, besides one large factor; at `p = 13` and `p = 23` it is irreducible.

**What row C does not do.** The finite places carry no obstruction: every `A_p` is strongly connected with a
spectral gap, and the trace formula holds for any affine pair. The Collatz conjecture is a statement at the real and
2-adic places; the near-return structure that THM-4659 found at the real place has only this weight-zero shadow at
the primes `p`. This is the honest contrast with rows T and H, where the finite places carry the whole arithmetic.

## 5. Row J: Keller maps over finite fields (PROVED given cited inputs; [script](../../04-computation/experiments/keller_sweep_value_set_law_20261017.py))

* **Planar criterion (assembled from standard results; no novelty claimed).** For `F = (P, Q) ∈ Z[x, y]^2` with
  `det JF = 1` the following are equivalent: (a) `F` is an automorphism; (b) `F mod p` permutes `F_p^2` for all large
  `p`; (c) for all large `p` and all `(a, b) ≠ 0`, `Σ_(x ∈ F_p^2) ψ(aP(x) + bQ(x)) = 0`; (d) for all large `p` there
  are no collisions `F(x) = F(x')`, `x ≠ x'`. (b)⇔(d) by finiteness of `F_p^2`; (c)⇔(d) by Parseval:
  `Σ_((a,b) ≠ 0) |S(a,b)|^2 = p^2 · #collisions`.
  (b)⇒(a): a non-automorphism Keller map has generic degree `d ≥ 2` (Keller: birational Keller maps are
  automorphisms), is finite étale of degree `d` off its Jelonek curve, and for the positive-density set of primes
  splitting completely in the constant field of its Galois closure, Chebotarev and the Cameron–Cohen derangement bound
  give `#F(F_p^2) ≤ (1 - 1/d) p^2 + O(p^(3/2))`. Every linear combination `aP + bQ` of a Keller pair has no affine
  critical point, so a counterexample's exponential sums would be supported entirely at infinity (Milnor number at
  infinity, Suzuki's formula `χ = 1 - μ - λ`; CITED).
* **The dimension-3 family: a value-set law.** For the tangent-sweep Keller maps of THM-4625 (THM-1300 is `d = 2`),
  statement (1) there makes the fibre over `(a, b, c)`, `c ≠ 0`, the set of simple roots of the degree-`(d+1)`
  polynomial `H_(b,c)(w) = c^2 a`. So the image of `F mod p` is a union of polynomial value sets, and with the
  monodromy `S_(d+1)` (THM-4625 (5)) Chebotarev gives fibre-size frequencies equal to the fixed-point law of `S_(d+1)` and
  image density `c_(d+1) = Σ_(k=1)^(d+1) (-1)^(k-1)/k!` (2/3, 5/8, 19/30, 91/144, ... → `1 - 1/e`), the
  Birch–Swinnerton-Dyer / S. D. Cohen value-set constant. THM-4625 (6)'s "fibre sizes never exactly `d`" is the shadow
  of "no permutation of `d + 1` letters fixes exactly `d`". FINITE illustration at `p = 61, 101`: `d = 2, 3` are the actual
  Keller maps (THM-1300, Gao's map; the audit reproduced the histograms by evaluating the maps on `F_p^3`); the `d = 4, 5`
  rows are only the one-variable family `H_(b,c)` (those `p(w)` fail Gao's side conditions, so they are not Keller maps,
  MISTAKE-624); observed gaps below 0.012 from the law, and no fibre of size `d`.
* **The spectrum of the dimension-3 counterexample ([THM-4664](../../01-canon/theorems/THM-4664-fourier-transform-of-the-thm1300-keller-map-mod-p-concentrates-on-the-dual-of-its-asymptotic-jelonek-cubic.md); PROVED decomposition, FINITE-EXACT sizes).**
  The exponential sums `S(α, β, γ) = Σ_(F_p^3) e_p(αF_1 + βF_2 + γF_3)` of THM-1300's map are
  `p e_p(-4β^2/α)(p[α^2γ = 4β^3] - 1)` plus two explicit two-variable sums of size at most `2.5p` (`p ≤ 13`); the large
  coefficients (`|S| > p^2/2`) are exactly the `p^2 - 1` nonzero points of the cone `α^2γ = 4β^3` for every prime
  `11 ≤ p ≤ 41`, carrying 92–98% of the Fourier energy. That cone is the projective dual of the cuspidal cubic
  `27a^2c + b^3 = 0`, which with the line `c = 0` is the curve at infinity of the Jelonek surface (THM-4621): the Keller
  defect's spectrum sits on the dual of the asymptotic non-properness locus. (A tempting guess, that the fibre count
  equals the number of roots of the plane-section cubic `2at^3 + bt^2 + c`, is false: it agrees at only about a third of the
  points; the phases on the cone matter.)
* **Where this leaves planar JC.** The finite places carry the planar obstruction completely (criterion above), but no
  planar candidate exists to test. The dimension-3 counterexamples are families of one-variable polynomials minus
  their critical points over a two-dimensional base; the planar analogue would need such a family over a
  one-dimensional base with total space `A^2`. This is the repo's slice question (THM-4562: THM-1300 is fully
  non-slice), and nothing here decides it.

## 6. Numerology avoided

`H(P_q) 2^(q-1)/q! → e` (THM-438) and the Keller image deficiency `→ 1/e` are unrelated Poisson limits (path counts
in a quasi-random tournament; fixed points of a random permutation). The level-15 form of the pentagons and the
`Q(sqrt(-15))` in the Polymath field `Q(sqrt(-3), sqrt(-11), sqrt(-15))` (THM-4558) have no connection found: the CM field
of the pentagon motive is a property of the closed five-chain, the Polymath field of particular rigid gadgets.
Both are NUMEROLOGY until a mechanism appears.

## 7. Next obligations

1. Row H: the theta moments `θ(2,3,3)`, `θ(3,3,3)` (higher moments of the `Γ_1(6)` surface), the five-bar surfaces
   `R_4(c)`, and an independent check of the seventh rung.
2. Row T/H: is there a single construction producing both elliptic modular surfaces from the respective connection
   sets (the quadratic residues; the unit circle)?
3. Row J: square-root bounds for THM-4664's remainders; the dual-cone law for Gao's map and the whole sweep family;
   and the planar prediction of HYP-9286: a planar counterexample's exponential sums would be large only on the finitely
   many frequency lines dual to its Jelonek curve's points at infinity.
4. Row C: none from this lens beyond the trace formula; the census programme (HYP-9283/9285) continues at the real place.
