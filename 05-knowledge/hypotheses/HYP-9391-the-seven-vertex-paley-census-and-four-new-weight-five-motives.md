---
id: HYP-9391
title: "The seven-vertex Paley census. (a) For primes q = 3 (mod 4), TT7(P_q) = q(q-1)[q^5 - 90q^4 + 3116q^3 - 49252q^2 + 322374q - 673743 - 2(14q^2 - 701q + 7254)a_q + 46 c_q - 83 d_q - 64(q - 31)|b(q)|^2 - 64 Y02 - 48 Y05 - 64 Y09 + 96 Y18]/2^21. Here a_q, c_q, d_q are the coefficients of 8.4.a.a, 8.6.a.a, 4.6.a.a, b = 32.3.c.a, and the Y are four explicit two-parameter Legendre triple products, so TT7 costs O(q^2). (b) The chi-twisted fourth moment of the Legendre family is sum_l chi(l) a_l^4 = -2q^2 - 3q - 1 - 3q a_q - c_q - d_q (the k = 6 case of Frechette-Ono-Papanikolas 2004, Prop. 2.1). (c) Modulo the known motives, the triple products span exactly four new directions. Numerically their new parts behave like four irreducible weight-5 motives, pairwise not quadratic twists and not GL(2) (normalised traces exceed 2), with trace statistics consistent with USp(4) but not distinguishing it from USp(2n), n >= 3: conjecturally the spinor motives of four genus-2 Siegel eigenforms."
status: >
  OPEN.
  (a) FINITE-EXACT: the formula reproduces all 71 directly counted values TT7(P_q), 7 <= q <= 811. It rests on three
  inputs:
  the exact reduced form of 2^21 TT7 (audit of THM-4780; PROVED structure);
  37 graph-sum identities S_H = (known combination) + (at most one triple product, coefficient +-1), FINITE-EXACT at
  43 primes 7..443 (23 primes 7..199 for the four sums without a degree-3 vertex);
  the wheel (HYP-9274 (B), PROVED for primes), the K6 - P6 terminal (HYP-9276, FINITE-EXACT) and (b).
  (b) PROVED for every odd prime (both classes mod 4): the k = 6 case of Frechette-Ono-Papanikolas (IMRN 2004, Prop. 2.1;
  part (1) Ahlgren, Finite Fields Appl. 8 (2002), part (2) Ahlgren-Ono, Monatsh. Math. 129 (2000)); see (b) below. Re-checked at the 428 odd primes 5..2999 (audit).
  A monodromy count predicts the two weight-6 forms; the untwisted companion (Ahlgren's Gamma_0(4) case) carries 4.6.a.a only.
  (c) NUMERICAL: 338 primes q = 3 (mod 4) up to 4999 (statistics over the 307 primes in [307, 4999]). Identification OPEN.
  Audited 2026-10-26 (tt7/motives/audit_hyp9391/): (a) re-derived with independent TT7 counts and Y sums, exact at all
  618 primes q = 3 (mod 4) in 7..9967; all 37 identities re-checked with brute-force sums at 18 primes 7..223
  (211, 223 out of sample, also for the four sums without a degree-3 vertex); Y04 = -Y02 and Y06 = Y05 PROVED; (b) is known (FOP 2004);
  in (c) the projection error at q >= 2000 is about 0.14-0.19, not < 0.05; the Gram and not-GL(2) conclusions survive, but
  the trace statistics do not single out rank 4 (MISTAKE-671).
session: mac-mini-2026-10-25-pyth (continuation, 2026-10-26)
source: 05-knowledge/results/pythagorean_hourglass_census_tower_20261025.md
depends_on: [THM-4780, HYP-9390, HYP-9274, HYP-9276]
scripts:
  - 04-computation/experiments/pythagorean_hourglass_20261025/tt7/motives/dump_genuine.py (-> genuine37.json): the 37 genuine graphs and multiplicities from the audit's engine
  - 04-computation/experiments/pythagorean_hourglass_20261025/tt7/motives/sum7.py, sum7_fast.py: exact O(q^5) and O(q^4) (degree-3 elimination) evaluators, validated against brute force
  - 04-computation/experiments/pythagorean_hourglass_20261025/tt7/motives/run_sums.py (-> sums_7_199.jsonl): all 38 survivor sums, sum m_H S_H = g(q) checked at its 23 primes; run_fast.py (-> fast_*.jsonl): the 33 genuine sums with a degree-3 vertex (no g check possible there)
  - 04-computation/experiments/pythagorean_hourglass_20261025/tt7/motives/classify_sums.py, exact_tt7_formula.py (+ .out, g_formula.json): the identities and the formula
  - 04-computation/experiments/pythagorean_hourglass_20261025/tt7/motives/cone_reduce.py, test_cones.py, closed_forms_check.py, run_cheap_Y.py (-> cheapY_*.txt): the triple products
  - 04-computation/experiments/pythagorean_hourglass_20261025/tt7/motives/legendre_moments.py (+ .out): part (b)
  - 04-computation/experiments/pythagorean_hourglass_20261025/tt7/motives/analyze_new.py, moments_new.py, moments_mod8.py, twist_tests.py, induced_test.py (+ .out): part (c)
  - 04-computation/experiments/pythagorean_hourglass_20261025/tt7/motives/identify_new.py, identify_products.py, identify_odd_levels.py, identify_Y.py, identify_classnumbers.py, candidate_forms_odd.gp: negative identifications (only candidate_forms_odd.out was saved; reruns in audit_hyp9391/rerun_identify_*.out)
  - 04-computation/experiments/pythagorean_hourglass_20261025/tt7/motives/run_cheap_Y4.py (-> cheapY_5000_9999.txt; writes the Y04, Y06 columns as -Y02, Y05), untwisted_moment.py (+ .out)
related:
  - 01-canon/theorems/THM-4780-paley-census-tower-to-six-vertices-wheel-octahedron-pgl2-hamiltonian-path-complement.md
  - 05-knowledge/hypotheses/HYP-9390-the-six-vertex-paley-census-is-the-level-32-adjoint.md
  - 05-knowledge/hypotheses/HYP-9276-k6-minus-p6-terminal-is-the-adjoint-of-a-weight-3-level-32-newform.md
---

# HYP-9391 — The seven-vertex Paley census and four new weight-five motives

Notation:
- `q ≡ 3 (mod 4)` is prime and `φ` is the quadratic character.
- `L(λ) = Σ_w φ(w(w−1)(w−λ)) = −a_λ` is the Legendre-family trace.
- `a_q`, `c_q`, `d_q` are the coefficients of `8.4.a.a = η(2z)⁴η(4z)⁴`, `8.6.a.a` and `4.6.a.a = η(2z)^12`.
- `|b(q)|²` is the square modulus of the coefficient of `32.3.c.a` (HYP-9276).
- The sums below run over `s ≠ t` with `s, t ∉ {0, 1}`.

## (a) The census in closed form

The four triple products are:
```
Y02 = − Σ φ(s) φ(t) φ(t−1) L(s)² L(t/s)          Y05 = − Σ φ(s−1) φ(t) φ(t−1) L(s)² L(t/s)
Y09 =   Σ φ(s) φ(s−1) φ(t) L(s) L(t) L(t/s)       Y18 =   Σ φ(s−1) φ(t−1) L(s) L(t) L(t/s)
```
Then
```
2^21 TT7(P_q) / (q(q−1)) = q⁵ − 90q⁴ + 3116q³ − 49252q² + 322374q − 673743 − 2(14q² − 701q + 7254) a_q
                           + 46 c_q − 83 d_q − 64(q − 31)|b(q)|² − 64 Y02 − 48 Y05 − 64 Y09 + 96 Y18 .
```

**How it was found.**
- By THM-4780's calculus (exact reduced form from its audit), `2^21 TT7 = A7 + B7 T_W4 + C7 T_{K6−P6} + Σ_{37} m_H T_H`.
- The 37 genuine sums were computed exactly at 43 primes (`sum7.py`, `sum7_fast.py`, `run_*.py`). Each was then fitted in the span of `1, q, q², q³, a_q, q a_q, c_q, d_q, |b|², Y02, Y05, Y09, Y18`.
- Every one is exact with integer coefficients:
  - 21 lie in the span of the known motives alone;
  - the other 16 are known terms plus exactly one `±Y` (`exact_tt7_formula.out`).
- The cone reductions are what bring in the triple products:
  - summing out three pairwise non-adjacent degree-3 vertices with the `L`-table leaves a four-point base `{0, 1, s, t}`;
  - the "merge" corrections lie in the known span;
  - `Y04 = −Y02` and `Y06 = Y05` hold for every prime `q ≡ 3 (mod 4)`. Proof (audit): with `A(u) = Σ_{s≠0,1} φ(s)φ(u−s)L(s)²`, one has `−Y02 = Σ_{u≠0,1} φ(u)L(u)A(u)` and `Y04 = Σ_{u≠0,1} φ(u)φ(u−1)L(u)A(u)`; on `{u : φ(u−1) = −1}` the fixed-point-free involution `u ↦ u/(u−1)` fixes `L(u)` and `A(u)` and flips `φ(u)`, so that part of `−Y02` vanishes. The same argument with `φ(s)φ(s−1)` in place of `φ(s)` gives `Y06 = Y05`. (`run_cheap_Y4.py` writes the `Y04`, `Y06` columns as `−Y02`, `Y05`, so in `cheapY_5000_9999.txt` they are not independent data.)

**Check.** The formula gives all 71 counted values of `TT7(P_q)` for `7 ≤ q ≤ 811` (`paley_tt7_data*.out`). It costs `O(q²)`. (Counting costs `O(q⁴)`, not `O(q⁵)`: the maps `x ↦ a²x + b` act simply transitively on the arcs, so `TT7 = (q(q−1)/2)·T5(N⁺(0) ∩ N⁺(1))`. The audit's counter uses this; the formula holds at all 618 primes `q ≡ 3 (mod 4)` up to 9967.)

## (b) The twisted fourth Legendre moment

```
Σ_{λ ≠ 0,1} φ(λ) a_λ⁴ = −2q² − 3q − 1 − 3q·a_q − c_q − d_q        (every odd prime; FOP 2004, Prop. 2.1, k = 6)
```

- The twisted second moment `Σ φ(λ) a_λ² = −q − 1 − a_q` is the wheel (Ahlgren–Ono).
- Prediction of the weight-5 part:
  - the local system `Sym⁴(H¹) ⊗ L_φ` on `P¹ − {0, 1, ∞}` has Euler characteristic `−5`;
  - it has one invariant (at `λ = 1`), so its parabolic `H¹` has dimension 4 and weight 5;
  - the identity finds that part to be `8.6.a.a ⊕ 4.6.a.a`. (For `Sym²` the same count gives `8.4.a.a`.)
- The identity is known. Frechette–Ono–Papanikolas (IMRN 2004, no. 60, 3233–3262, Prop. 2.1: (1) Ahlgren, Finite Fields Appl. 8 (2002), the `η(2z)^12` fivefold; (2) Ahlgren–Ono, Monatsh. Math. 129 (2000); stated for every even `k ≥ 4`) give `Tr_k(Γ₀(4), p) = −3 − Σ_{λ=2}^{p−1} G_k(a_λ, p)` and `Tr_k(Γ₀(8), p) = −4 − Σ_{λ=2}^{p−2} G_k(a_{λ²}, p)`, where `G_6(s, p) = s⁴ − 3ps² + p²` is the `Sym⁴` trace. Since `Σ_{λ=2}^{p−2} f(λ²) = Σ_{μ≠0,1} (1 + φ(μ)) f(μ)`, `Tr_6(Γ₀(8), p) = c_p + 2d_p` and `Tr_6(Γ₀(4), p) = d_p`, this gives `Σ φ(μ) Sym⁴(a_μ) = −1 − c_p − d_p` and `Σ Sym⁴(a_μ) = −3 − d_p`; with `a⁴ = Sym⁴ + 3p·Sym² + 2p²` and the `k = 4` case these are the twisted and untwisted identities. So (b) is PROVED for all odd primes.
- The untwisted companion is `Σ_{λ≠0,1} a_λ⁴ = 2q³ − 4q² − 9q − 3 − d_q` (76 primes; `untwisted_moment.out`). Only `4.6.a.a` appears, the classical `Γ₀(4)` moment, so `8.6.a.a` enters through the `φ(λ)` twist.
- In the census, the `c_q`, `d_q` terms are these fourth moments. For example `K_{3,4}` (four degree-3 vertices on one untwisted triple) has `S = −171 + 50q + 4q² − q³ − d_q`.

## (c) Four new motives

Here `N` is the new part: the residual of a `q^{−5}`-weighted least-squares fit of each `Y` on the nine known columns, normalised by `q^{5/2}`.
- **Number and irreducibility.**
  - Modulo the known span, `Y02, Y05, Y09, Y18` span exactly four new directions (exact rank).
  - The Gram matrix of their new parts is the identity, to within 0.07, in the basis `m1 = N(Y02)`, `m2 = N(Y09)`, `m3 = N(Y18)`, `m4 = N(Y05) + N(Y18)`.
  - `E[m_i²] = 0.96, 0.93, 0.96, 0.96` over the 307 primes in `[307, 4999]`. So, heuristically (integral multiplicities and Rankin–Selberg orthogonality on the class `q ≡ 3 (mod 4)`, where a motive and its `χ₋₄`-twist differ only by sign), there are four distinct irreducible motives of weight 5, each entering one `Y` (or `Y05 + Y18`) once.
- **Not GL(2).**
  - At the 184 primes `≥ 2000`, `max |m_i| = 2.59, 3.16, 3.10, 3.01`.
  - The share of those primes with `|m_i| > 2` is 3–6.5%; USp(4) predicts 5.0%.
  - A GL(2) motive satisfies `|x| ≤ 2`, so the rank is at least 3. An irreducible self-dual motive of odd weight is expected (not known in general) to be symplectic, hence of even rank, so the rank is at least 4.
  - The fitted known part is uncertain by about 0.14 (mean) to 0.19 (max) per prime there (audit, white-noise model; coefficient standard errors about 21 for `q²`, 2.3 for `|b|²`, 100 for `a_q`, 2·10⁴ for `q`), and reasonable refits move single values by up to about 0.3. The conclusion survives: on the 464 primes in `[2000, 9999]` (audit refit on 587 primes) 3–8 values per motive exceed 2.5, against an expected 0.02–0.5 for a GL(2) motive plus that error.
- **USp-type statistics.** Compared with the exact reference values. They do not fix the rank: `E[x⁴] = 3` for every `USp(2n)` with `n ≥ 2`, `E[x⁶] = 14` for USp(4) against 15 for `n ≥ 3`, and on the audit's 587 primes the KS distances to USp(6), USp(8) and `N(0, 1)` are as small as to USp(4).

  | | USp(4) | SU(2) | `m1` | `m2` | `m3` | `m4` |
  |---|---|---|---|---|---|---|
  | `E[x⁴]` | 3.0 | 2.0 | 3.21 | 2.34 | 2.60 | 2.70 |
  | `E[x⁶]` | 14.0 | 5.0 | 17.2 | 9.9 | 10.9 | 11.3 |
  | `P(\|x\| < 0.2)` | 0.169 | 0.127 | 0.179 | 0.150 | 0.169 | 0.156 |

  The standard error of `E[x⁴]` is about 0.5 at 307 primes. A tensor product `SU(2)×SU(2)` would give `E[x⁴] = 4` and `Sym³(SU(2))` gives 4 (with `E[x⁶] = 34`). On the audit's 587 primes `[307, 9999]`: `E[x⁴] = 3.13, 2.96, 2.88, 2.70`, `E[x⁶] = 16.0, 13.0, 13.6, 10.5`, and a KS test rejects `SU(2)×SU(2)` and `Sym³` at the 5% level for all four.
- **No coincidences.**
  - No quadratic twist relates two of them, and none is a self-twist: for characters `(2/p), (±3/p), (5/p), (7/p)`, every `|E[χ m_i m_j]|` (including `i = j`) is at most 0.18, about two standard errors, against 1 for a twist.
  - None is induced from `Q(√D)` for the `D` tested: every `E[χ_D m_i²]` is near 0 rather than 1, and no class mod 8 has `E[x²]` near 0. The share of `|x| < 0.05` (2–6% on the 307 primes) is not a sharp test, because the projection error (about 0.15) is three times the window: an induced motive would show about 15%, not 50%.
- **Negative identifications** (none removes a new direction; tested at 31–71 primes; `identify_*.py`):
  - 97 newform orbits of 2-power level up to 256 (weights 2–6; characters `1, χ8, χ−4, χ−8`);
  - 260 orbits of levels `2^a·{3, 9, 5, 7, 15, 21, 45, 63}` (weights 2, 4, 6);
  - 329 products and cubes of these;
  - `|b|²` and `q²` times weight-2 orbits;
  - twisted polynomials.
- **Conjecture.** The four motives are the spinor motives of four genus-2 Siegel modular eigenforms with 2-power level (paramodular or `Γ₀`-type). Motivic weight 5 with a nonzero `(5,0)` Hodge part means scalar weight 4, or vector-valued weight `(3, 2)` (or the non-cohomological `(2, 4)`). Tate twists such as `q·M` (with `M` of paramodular weight 3) or `q²·A` (an abelian surface) also have weight 5. The audit's mod-`q` test excludes them heuristically: no congruence `Y ≡ c₀ + c₁a_q + c₂c_q + c₃d_q + c₄|b|²` (mod `q`) with small integers `c_i` holds for any of the four, even after adding `(2/q)`, `(3/q)`, `h(−q)`, `(2/q)h(−q)`, `h(−8q)`, while the positive controls are recovered (`modq_hodge_test.out`). So each new trace is nonzero mod `q`. Rank 4, hence 'genus 2', is part of the conjecture.
- **Mechanism.**
  - Heuristic only. Multiplicative convolution of hypergeometric sheaves with disjoint parameter lists concatenates the lists (Katz, ESDE ch. 8). So the bare sum `Σ_s L(s)L(t/s)` is the rank-4 symplectic family of `₄F₃(½,½,½,½; 1,1,1)`.
  - The actual `s`-integrands are `φ(s)L(s)²` (Y02), `φ(s−1)L(s)²` (Y05), `φ(s)φ(s−1)L(s)` (Y09) and `φ(s−1)L(s)` (Y18). `L(s)² = Sym² + q` has rank 4, and the twist by `φ(s−1)` is not a Kummer twist on `G_m` (its local monodromy at `s = 1` is `−J₂`, not a pseudo-reflection). So none of the four is a hypergeometric convolution; the relevant operation is Katz's middle convolution of rigid local systems.
  - For Y02 and Y05 the second Legendre factor sits in `s` (a rank-4 sheaf convolved with the Legendre sheaf, then paired with `φ(t)φ(t−1)`); for Y09 and Y18 a rank-2 sheaf is convolved with the Legendre sheaf and tensored with `L(t)`. The ranks need recounting: with middle extensions the Euler–Poincaré formula (generic rank `2·r_F·r_G − dim F^{I₁}·r_G − r_F·dim G^{I₁}` for the convolution) gives 8 for Y02 but 12 for Y05, Y09 and Y18, whose `s`-integrands have no inertia invariants at `s = 1`. The new motives are pieces of `H¹_!` of these local systems.
- **Caveat.** The split into "known part + new motives" is statistical. The coefficients of `1, q, q², a_q` in the projection are not determined exactly, so the exact integer traces of the four motives are not yet available. The rank and Sato–Tate evidence uses only large primes, where the projection error is about 0.15 per prime; small-value frequencies are therefore unreliable.

## Next steps

1. **Exact separation.** Compute the boundary and degenerate contributions of each `Y` exactly, giving exact traces `t_k(q)`.
2. **Euler factors.** Evaluate the `Y` over `F_{q²}` (and `F_{q³}`) for small `q` to get degree-4 Euler factors.
3. **Conductor and weight.** Test the functional equation with PARI `lfun` (degree 4, motivic weight 5) over 2-power conductors.
4. **Comparison.** Check against the LMFDB / Poor–Yuen paramodular data.
5. **Seven vertices at prime powers.** Extend the formula of (a) to prime powers, as THM-4780 does at six vertices.
