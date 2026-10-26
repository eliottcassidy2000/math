# Pythagorean hourglass, busy-beaver analogies, and the Paley census tower

2026-10-25, mac-mini-2026-10-25-pyth. Owner prompt (pasted block plus the owner's own text):
- Cauchy–Schwarz as the triangle inequality; CLT and Tychonoff "analogously".
- "Quadratic reciprocity just describes that a tournament on 4 vertices has 4 isomorphism classes".
- The area under a curve as a Möbius-like band that flips its face below zero.
- FTC with FTA, and "the infinite spiraling Riemann surface … the strip returns".
- Rule 30.
- Past analogies between primitive Pythagorean triples and busy-beaver problems such as Collatz.
- The twisted-hourglass lifting statement of THM-4773: one squareclass; the central conic is the unit circle for the twist `(1,1,−1,1,−1,1,−1)`; `(3,4,5) ↦ (0:3:3:4:4:5:5)`; genus five on the outer lines; the untwisted surface has none.

Canon produced:
- [THM-4780](../../01-canon/theorems/THM-4780-paley-census-tower-to-six-vertices-wheel-octahedron-pgl2-hamiltonian-path-complement.md) (PROVED).
- [HYP-9390](../hypotheses/HYP-9390-the-six-vertex-paley-census-is-the-level-32-adjoint.md) (OPEN, conditional on HYP-9276).
- HYP-9391 was reserved and is released unused.

Scripts are in `04-computation/experiments/pythagorean_hourglass_20261025/`. Every `.py` there has its `.out` beside it.

## 0. Answers in brief

| Clause | Verdict | Where |
|---|---|---|
| THM-4773 (hourglass lifts) | **Verified independently**, every part | §1 |
| Pythagorean triples vs busy beaver / Collatz | A **contrast**, not an analogy: Berggren descent halts by a Lorentzian potential; Collatz has no invariant form | §2 |
| QR "is" the four 4-vertex classes | The numeral four is a coincidence (klein's verdict stands). The four classes nevertheless carry reciprocity-level arithmetic in three exact ways | §3 |
| (new) the transitive census of Paley tournaments | `TT4` = the Gauss identity, `TT5` = the Apéry form, `TT6` = the level-32 adjoint; **THM-4780**, **HYP-9390** | §4 |
| Möbius band, FTC/FTA spiral, CS/triangle, CLT/Tychonoff, Rule 30 | Mostly covered by klein's atlas. New: the square-root band is the hourglass's lifting monodromy, the χ(−1) face flip is the census vanishing law, and FTA = FTC for `log p` on the spiral | §5 |

## 1. THM-4773 checked independently

`verify_thm4773.py` was written from the statement alone. It shares no code with klein's packet.

- **(A) Pointwise criterion.** The criterion and the count `2^{6−r}` of lifts agree with brute force on 9,995 (twist, point) pairs, with 0 disagreements.
- **(B) Central line.** The Hilbert-symbol decision of the conic `w² = d1d5d3²s² + d3d5t²` was compared with a brute-force search (`|x|, |z| ≤ 40`) on 504 twists.
  - 465 agree.
  - The 39 "disagreements" are all twists that are soluble but have no point inside the box. All 39 have explicit lifts of height up to 625 (B'), so the decision is exact on all 504.
  - The displayed twists come out as stated: Pythagorean yes; `(1,1,−1,1,−1,7,−7)` no; `(1,2,−2,3,−3,6,−6)` no; untwisted no.
- **(C) Genus.** Riemann–Hurwitz `2g − 2 = 2^k(−2) + b·2^{k−1}` gives the following.
  - The central line `E = 0`: three branch points, a rank-2 cover, eight components of genus 0.
  - Each of the six outer lines: five branch points, a rank-4 cover, two components of genus 5. Faltings then makes the lifts finite.
- **(D) Untwisted.** No boundary lift with `|x|, |z| ≤ 60`.
- **(E) Pythagorean twist.** The central base points `(0:x:z)` with `0 < x, z < 400` that lift are exactly the 6 of the form `(a², b²)` with `a² + b² = c²`: `(3,4,5)`, `(4,3,5)`, `(5,12,13)`, `(12,5,13)`, `(8,15,17)`, `(15,8,17)`.

So the lifted central points are the primitive Pythagorean triples. `(3,4,5) ↦ (0:3:3:4:4:5:5)` is the root of Berggren's tree, and every other lift has its Berggren address.

## 2. Pythagorean triples versus busy-beaver problems

**Prior repo work.**
- [`odd_square_20260921_synthesis.md`](odd_square_20260921_synthesis.md).
- [`collatz_mod6_20260921_berggren_edge_transport.md`](collatz_mod6_20260921_berggren_edge_transport.md).
- [`ternary_berggren_20260925.md`](ternary_berggren_20260925.md). It contains the separation: `tr²/det` is an integer on every Berggren word and on no Collatz word.
- `BERGGREN_POSITIVITY_NONREGULAR.md`.
- THM-4591 (Pythagorean *tuning*) is a false friend.
- The BB(6) cryptids (THM-4770, THM-4771, HYP-9380) include Antihydra, whose odd-orbit walk has halting probability `φ^{−6}` in the i.i.d. model (THM-4771).

**What is new here** (`berggren_vs_collatz.py`, FINITE-EXACT).

1. **Berggren descent provably halts.**
   - The parent map lowers the hypotenuse strictly. All 31,819 primitive triples with `c ≤ 200,000` reach `(3,4,5)`.
   - The potential `c` is the Lorentz-invariant height: the descent is a walk in `O(2,1; Z)` with a positive potential.
   - Collatz has no such potential: no quadratic form is preserved, which is the separation of `ternary_berggren_20260925.md`.
2. **The two depth laws differ in kind.**

   | range | Berggren depth: mean, sd²/log | Collatz stopping time: mean, sd²/log |
   |---|---|---|
   | `[10³, 2·10³)` | `1.09 log c`, 4.2 | `6.82 log n`, 96.4 |
   | `[10⁴, 2·10⁴)` | `1.34 log c`, 17.2 | `6.85 log n`, 90.8 |
   | `[10⁵, 2·10⁵)` | `1.60 log c`, 60.1 | `6.83 log n`, 90.7 |

   - Collatz stopping times are Gaussian-like at the drift `log n / log(2/√3) = 6.95 log n`.
   - Berggren depths are heavy-tailed and grow like `(log c)²`. The mean tracks `(6/π²)(log √c)²`, the Yao–Knuth law for sums of partial quotients.
   - Depth correlates 0.92 with the partial-quotient sum `S(m/n)` of the generator ratio.
   - The Pythagorean tree is a continued-fraction (Euclid/Stern–Brocot) object. Collatz is a random walk with drift.
3. **The hourglass link is exactly this tree.**
   - Under the Pythagorean twist the central conic becomes the unit circle (THM-4773).
   - Its rational points are the triples, and the twisted lifts inherit Berggren's ternary structure.

**Verdict.** The "Pythagorean triples are busy-beaver-like" analogy fails at one precise point, the invariant quadratic form. With the form, halting is a theorem (descent in a lattice of `O(2,1)`). Without it, halting is the open problem. All the past attempts in the repo that went beyond vocabulary ended at that separation.

## 3. Reciprocity and the four 4-vertex classes

klein's atlas (`analogy_reduction_atlas_20261025/QR_TOURNAMENT.md`, HYP-9375) is right that the numeral four is a coincidence: `A000568(4) = 4` and the four cases `(p, q) mod 4` have no common source. Three exact statements make the classes arithmetic anyway.

1. **Paley level (THM-4780 (7)).**
   - In a regular tournament, `TT4 ≥ n(n−1)(n−3)(n−7)/64`, with equality iff the tournament is doubly regular.
   - `TT4` fixes all four class counts.
   - For `P_q`, equality is `S² = J − qI`, i.e. the Gauss identity `g² = −q`. So the four class frequencies of `P_q` are the Gauss identity.
   - Checks: exhaustive over the regular classes at `n = 7, 11` (3 and 1223 classes); random regular tournaments at 15, 19, 23; `P_q` for `q ≤ 43` (`qr_four_census.out`).
2. **The first supplement.**
   - At `q ≡ 3 (mod 4)` the `K4` sum vanishes because the transposition `(3 4)` reverses one arc (THM-4652).
   - At `q ≡ 1 (mod 4)` the same sum is `a_q(η(4z)⁶)`, the CM form of `Q(i)`, which vanishes exactly at the inert primes `q ≡ 3 (mod 4)`.
   - So "`χ(−1) = −1`", "Paley is a tournament" and "the `Q(i)` form vanishes" are one fact.
   - On the graph side, `K4(G_q) = q(q−1)(q² − 20q + 81 + a_q(η(4z)⁶))/1536`. This is FINITE-EXACT at 59 primes (`paley_graph_k4.out`) and equivalent to Evans–Pulham–Sheehan's `q(q−1)((q−9)² − 4y²)/1536` with `q = x² + y²`, `y` even.
3. **Four primes ≡ 3 (mod 8).**
   - Put `p_i → p_j` iff `(p_j/p_i) = 1`. This is a tournament by reciprocity.
   - The pure 2-Selmer dimension of `y² = x³ − n²x`, `n = p1p2p3p4`, is the Monsky nullity `dim_F2 ker(A² + A + I)` (`bsd_selmer_growth_20261013a.md` §§2–3).
   - It depends only on the isomorphism class: TT4 gives 0, source + C3 gives 2, sink + C3 gives 2, strong C4 gives 0. This is exhaustive over the 64 labelled tournaments.
   - All four classes occur, with the full Monsky matrix as an independent check:
     - `(3, 11, 59, 83)`: TT4, 0;
     - `(3, 11, 19, 107)`: source + C3, 2;
     - `(3, 11, 19, 163)`: sink + C3, 2;
     - `(3, 11, 19, 43)`: strong, 0.
   - The 2 occurs exactly on the classes with one cyclic triple, which is THM-1455's `|Pf| = 3` dichotomy.

## 4. The Paley census tower (THM-4780, HYP-9390)

**Method.** Expand the indicator of a transitive `k`-set into `2^{C(k,2)}` character sums `T_G` over subgraphs of `K_k`.
- THM-4652's calculus eliminates vertices of degree `≤ 2`. Its vanishing laws then kill every dense graph with an odd number of edges, or with an automorphism that reverses an odd number of arcs, since `χ(−1) = −1`.
- A new PGL(2) lemma reduces the dense graphs whose degrees are all even.
- What survives:

| `k` | `TT_k(P_q)` (times `2^{C(k,2)}/(q(q−1))`) | surviving dense sums | arithmetic |
|---|---|---|---|
| 3 | `q − 3` | none | none |
| 4 | `(q−3)(q−7)` | none (`K4` killed) | the Gauss identity (doubly regular) |
| 5 | `(q−3)(q²−21q+105) + S_W4` | the wheel `W4` (multiplicity +1) | `S_W4 = −4 − a_q`, Apéry form `η(2z)⁴η(4z)⁴` (`₄F₃(1)`, Ahlgren–Ono) |
| 6 | `(q−3)(q³−46q²+682q−3239) + 7(q−15)S_W4 − 8S_{K6−P6}` | `K6 − P6` (−8); the octahedron (+1), which reduces to `(q−5)W4` | `S_{K6−P6} = |b(q)|² − q² − 7`, the adjoint of `32.3.c.a` (THM-4653 + HYP-9276) |
| 7 | open | 37 genuine + 1 reducible (HYP-9390 (c)) | open |

**The formulas.** In closed modular form, for primes `p ≡ 3 (mod 4)` (`TT6` conditional on HYP-9276):
```
TT5(P_p) = p(p−1)(p³ − 24p² + 168p − 319 − a_p)/2^10                                   (THM-4780 (6), PROVED)
TT6(P_p) = p(p−1)(p⁴ − 49p³ + 828p² − 5313p + 10193 − 7(p−15)a_p − 8|b(p)|²)/2^15        (HYP-9390)
```

**Novelty.** McCarthy–Springfield (arXiv:2311.02135) give orders 3 and 4 only, for `k`-th power Paley digraphs; we found no order-5 or order-6 count.

**What was proved and how.**
- The census formulas in terms of `S_W4` and `S_{K6−P6}` follow by exact interpolation:
  - the reduction gives a priori degree bounds;
  - the systems (7 unknowns at 41 primes, 10 unknowns at 24 primes) have full rank;
  - the signed multiplicities come out of the solve as `+1` and `−8`, matching the independent count.
- `−8` is the Désarménien–Foata signed Eulerian value `½ Σ_{σ∈S6} sgn(σ)(−1)^{des σ} = ½·2³A_3(−1)`. A relabelling reverses `inv − des` complement arcs of the Hamiltonian path.
  - This is the same tangent-number family as HYP-9366's `F(T, −1)`: `A_3(−1) = −T_3`.
- **The octahedron.** HYP-9275 (opus S35) left `U(K_{2,2,2}) = q(q−1)((5−q)a_q + 3q)` open. The PGL(2) lemma gives `T_oct = (q − 5)T_W4`:
  - every in-degree of the order-oriented octahedron is even;
  - so the six-point sum is invariant on `P¹`;
  - `P¹` has `q + 1` points, and the six classes with one vertex at `∞` are each a wheel.
  - The kernel-partition bookkeeping converts this into HYP-9275's all-maps normalisation (THM-4780, part 4). So the octahedron identity is PROVED at every odd prime.
- **Checks.**
  - The PGL(2) identity at all 27 odd primes below 110, in both classes.
  - HYP-9275's `U` reproduced by brute force at `q = 5, 7, 11`.
  - The census at the prime powers `27, 243, 343, 1331, 2187`, with Frobenius traces.
  - Zeros at `q = 3, 7, 11` (`TT5`) and `3, 7, 11, 19, 23, 27` (`TT6`); `P_27` is ST27, the unique largest `TT6`-free tournament.
  - Each zero is an identity between modular coefficients, e.g. `TT5(P_7) = 0` ⟺ `a_7 = 24`.

**The ladder.** The 2026-10-14 notes found the same hierarchy in the Paley run expansion of Hamiltonian paths: `₂F₁` inside every terminal, then `₄F₃(1)` for the wheel, then `₂F₁·₃F₂` for `K6 − P6` (`paley_k6p6_adjoint_weight3_level32_20261014.md`). The census tower is that ladder seen through transitive subtournaments.
- The adjoint first appears at six vertices, through the Hamiltonian-path complement, as in the runs.
- **Consequence for the loop expansion.** For `P_q`, every five-vertex class count is affine in `a_q`, and every six-vertex class count is affine in `(a_q, |b(q)|²)`. So the four-loop term of HYP-9371 is a rational function of `q` plus a rational function times `a_q`, and the five-loop term adds `|b(q)|²`.

**Negative results (receipts kept).**
- `A(p) + B(p)a_p` alone does not fit `TT6`.
- The weight-6 level-8 newform `8.6.a.a` (`q + 20q³ − 74q⁵ − 24q⁷ + …`), the CM forms of `Q(√−2)`, `Q(√−3)`, `Q(√−7)` (weights 3 and 5), `(2/p)` twists and `a_p²` all have coefficient 0 or fail (`tt6_adjoint_search.out`, `tt6_search.out`).

## 5. The remaining clauses

**Cauchy–Schwarz and the triangle inequality.** For a norm that comes from an inner product they are equivalent:
- CS gives `‖x+y‖² ≤ (‖x‖+‖y‖)²`;
- conversely, expanding `‖x + ty‖ ≤ ‖x‖ + |t|‖y‖` for small `t` gives CS.
- The `ℓ¹` norm is the hostile example: triangle inequality without CS.
- See klein's atlas. Nothing is added here.

**CLT and Tychonoff.** The shared content is compactness plus a fixed point:
- tightness and Prokhorov (a Tychonoff/Banach–Alaoglu compactness) give limit points;
- the Gaussian is the fixed point of the renormalised convolution.
- klein's atlas covers this. THM-4761's `√n(κ−1) ⇒ N(0, 2)` is the repo's own conditional CLT.

**The band that flips its face.**
- As stated, the region under a graph is an orientable strip, and the face flip is just the sign of `f`.
- A genuine Möbius band appears for **square roots**. Around a sign change of `L`, the two branches of `√L` are exchanged, and the double cover is a Möbius band over a loop.
- That is literally the hourglass. A base point lifts through `s_i² = d_iL_i` iff the seven Möbius covers can be trivialised simultaneously, i.e. the nonzero `d_iL_i` share one squareclass (THM-4773 (1)).
- In the census the face flip is the factor `χ(−1) = −1` acquired by reversing an arc:
  - a graph with an odd automorphism (an odd number of forced flips) cannot be untwisted, and its sum vanishes;
  - the surviving sums `W4`, `K6 − P6` and the octahedron are exactly the graphs whose flips can be paired.
- This is a dictionary, not a theorem.

**FTC, FTA and the spiral.** "The infinite spiralling Riemann surface is the only topology that permits every polynomial to have a root" is not correct as stated: FTA is a statement about `C`, proved by the degree of `p` on large circles. The precise version, in the spirit of the prompt:
- the log spiral is the universal cover of `C*`;
- FTA is FTC for `d log p` along a large circle lifted to the spiral;
- the lift climbs `deg p` sheets, i.e. `(1/2πi)∮ p'/p = deg p ≠ 0`;
- a loop that climbs sheets cannot be contracted in `C*`, so `p` has a zero inside.
- "The strip returns" is the passage from the log spiral (`Z` monodromy) to the square-root band (`Z/2` monodromy, the quotient by `2Z`). That band is what the hourglass and the census use.

**Rule 30.** The repo has Rule 30 theorems THM-3456..4210, with their transfer barriers, and Collatz as a base-6 automaton (`collatz_oscillation_gilbreath_20260926.md`). No new claim is made.

## 6. Open problems

1. **Seven vertices** (HYP-9390 (c)). Compute `TT7(P_q)` for `q ≤ 250` and test the span of `{q^i, q^i a_q, q^i |b(q)|²}`. The 37 genuine sums either cancel or bring new arithmetic.
2. **HYP-9276** (`M_q = −b(q²) − 2q − 1`). Proving it makes `TT6` unconditional. The routes are listed in the 2026-10-14 note: Hijikata trace, Katz middle convolution, an elliptic K3 pullback.
3. **The wheel identity at prime powers.** It is FINITE-EXACT at five prime powers. Extend Ahlgren–Ono's threefold count to `F_q`.
4. **The McCarthy–Springfield cross-link.** Is their weight-3 level-32 non-CM form the same `32.3.c.a`? If so, it is linear in the 4-vertex count of the fourth-power Paley digraph and adjoint in the 6-vertex count of the Paley tournament.
5. **The Legendre tournament.** The four-prime Selmer split extends to `n` primes: the dimension is `dim ker(A² + A + I)` of the Legendre tournament. Its distribution over isomorphism classes is a tournament-census question, and klein's HYP-9375 quasirandomness would make it computable.

## 7. Reproduction

Run from `04-computation/experiments/pythagorean_hourglass_20261025/`.

| script | output | what |
|---|---|---|
| `verify_thm4773.py` | `.out` | §1 |
| `berggren_vs_collatz.py` | `.out` | §2 |
| `qr_four_census.py` | `.out` | THM-4780 (7) |
| `legendre_tournament_selmer.py` | `.out` | §3.3 |
| `paley_graph_k4.py` | `.out` | §3.2 |
| `paley_census_tower.py` | `.out` | `TT4`, `TT5` table to 199 |
| `paley_tt5_arith.py 1500` | `paley_tt5_data.out` | `TT5(P_p)`, `p < 1500` |
| `tt5_exact_poly.py`, `tt5_residual_fit.py`, `tt5_scaling.py`, `tt5_modular_search.py` | — | the TT5 search history |
| `tt5_exact_f8.py` | `.out` | exact TT5 = A + B·a_p fit |
| `paley_graph_sums.py`, `graph_sum_modular_search.py` | `.out` | `K4`, `K5`, `W4` sums, all odd `p ≤ 419` |
| `paley_tt6.py 760` | `paley_tt6_data.out` | `TT6(P_p)`, `p ≤ 751` |
| `tt6_adjoint_search.py`, `tt6_search.py` | `.out` | §4 search and negative results |
| `graph_sums6.py 200` | `.out` | the nine dense six-vertex sums |
| `census_interpolation.py` | `.out` | THM-4780 (5) |
| `octahedron_pgl2_check.py` | `.out` | THM-4780 (4) and HYP-9275's normalisation |
| `census_level_forecast.py 7` | `.out` | survivors, `k = 4..7` |
| `paley_prime_powers.py` | `.out` | prime powers |

Sources:
- McCarthy–Springfield, arXiv:2311.02135.
- Ahlgren–Ono (2000), `₄F₃(1)`.
- Evans–Pulham–Sheehan (1981), `K4` in Paley graphs.
- Désarménien–Foata (1992), signed Eulerian numbers.
- Sanchez-Flores (Graphs Combin. 1996) and Neiman–Mackey–Heule (arXiv:2011.00683), ST27.
