# Independent audit: THM-4780, HYP-9390 and the Pythagorean-hourglass results note

2026-10-25. Worktree `math-wt-fanopsp`, audited at commit `5df2594f12` (mac-mini-2026-10-25-pyth checkpoint).

Audited documents:
- `01-canon/theorems/THM-4780-paley-census-tower-to-six-vertices-wheel-octahedron-pgl2-hamiltonian-path-complement.md`
- `05-knowledge/hypotheses/HYP-9390-the-six-vertex-paley-census-is-the-level-32-adjoint.md`
- `05-knowledge/results/pythagorean_hourglass_census_tower_20261025.md`

Independence of the checks:
- All code in this directory was written for the audit.
- Lane scripts were read, never imported. Lane outputs were loaded only to compare numbers.
- Modular data come from PARI/GP's `mf` package (`modular_coeffs.gp`), not from the eta/theta products the lane used.
- Finite fields `GF(p^r)` use the audit's own log-table implementation (`audit_gf.py`), with a different irreducible polynomial from the lane's.
- No repository file was edited. The lane's untracked `tt7/` directory appeared during the audit and was not touched.

(The auditing agent could not write this file itself; the session saved its report text verbatim.)

## Verdict

No error invalidates a theorem-level claim. Every number in the three documents that I re-derived is correct. The census formulas of THM-4780 (5) were re-derived without interpolation, by an independent symbolic implementation of the elimination calculus.

The fixes needed:
- two false statements in HYP-9390: a sign error in the Sym² reading of `|b(q)|²`, and a wrong coefficient ring in Remark 1;
- wrong credits: the elimination calculus is not THM-4652's;
- a novelty overclaim for the PGL(2) invariance;
- uncited classical results;
- about a dozen minor wording and citation issues.

All 41 recommended edits are in `required_edits.json`. Each `old` string occurs exactly once in the current file, and the edits for each file apply in sequence (`make_edits.py`, `make_edits.out`).

| Document / claim | Verdict |
|---|---|
| THM-4780 (1) expansion | CORRECT |
| THM-4780 (2) reduction lemma, degree bounds | CORRECT WITH FIXES (wording; credit) |
| THM-4780 (3) vanishing laws and survivor table, k = 4..7 | CORRECT |
| THM-4780 (4) PGL(2) lemma, `T_oct = (q−5)T_W4`, the `q ≡ 1 (mod 4)` case | CORRECT WITH FIXES (clarity; the invariance is classical) |
| THM-4780 (4) HYP-9275 normalisation and the octahedron identity at every odd prime | CORRECT WITH FIXES (misattribution of `U(W4)`) |
| THM-4780 (5) `TT5`, `TT6` formulas, multiplicities `+1`, `−8`, Désarménien–Foata, `M_q` form | CORRECT |
| THM-4780 (6) `TT5(P_p)` with the Apéry form | CORRECT (attribution wording) |
| THM-4780 (7) `TT4` lemma, class counts, Gauss identity | CORRECT; classical, uncited |
| THM-4780 data (zeros, ST27, prime powers), status PROVED | CORRECT; credit for earlier prime-power checks missing |
| HYP-9390 (a) formula, conditional proof, FINITE-EXACT at 75 `q` | CORRECT; the prose reading of `|b(q)|²` in title and body is WRONG (sign) |
| HYP-9390 (b) prime-power wheel identity | CORRECT as FINITE-EXACT; earlier checks uncredited |
| HYP-9390 (c) seven-vertex survivors and degree bounds | CORRECT; the multiplicity list is INACCURATE |
| HYP-9390 Remark 1 (coefficient rings) | WRONG |
| HYP-9390 Remarks 2, 3; search history | CORRECT / minor fix / one negative result UNSUPPORTED (`Q(√−3)`) |
| Results note §1 (THM-4773 check) | not audited (out of scope) |
| Results note §2 (Berggren vs Collatz) | numbers CORRECT (reproduced); three statements wrong or unsupported (minor) |
| Results note §3 (Gauss identity, `K4`/EPS, four-prime Selmer split) | CORRECT (one convention note) |
| Results note §4 | CORRECT WITH FIXES (credit; conditional statements unmarked) |
| Results note §5 | CS/FTA/CLT CORRECT; the `ℓ¹` and Möbius-band sentences imprecise |
| Results note §7 and sources | minor errors |

## What was run (all outputs in this directory)

| script | output | check | method | result |
|---|---|---|---|---|
| `audit_tt_census.py` | `.out`, `audit_tt_values.out` | A | `TT5`, `TT6` by the source–sink decomposition `TT_k = q(q−1)/2 · TT_(k−2)(T[M(0,1)])` (the lane used the source decomposition `q·TT_(k−1)(T[QR])`), and by plain DFS enumeration of transitive chains with no symmetry | agree with the lane at all 119 primes (`TT5`), 67 primes (`TT6`), the 5 prime powers, and `q = 3, 7, 11`; DFS agrees at 11 values `q ≤ 71`, including 27 |
| `modular_coeffs.gp` | `modular_coeffs_gp.out` | G | `a_p(8.4.a.a)` and `b(p)` (32.3.c.a, character Conrey 32.31 = `χ₋₄`) from `mfinit` / `mfeigenbasis` for `p ≤ 2179` | `b(3), b(7), b(11), b(19), b(23), b(31) = 4i, −8i, −4i, −12i, 40i, −32i`, matching HYP-9276 |
| `audit_graph_sums.py` | `.out` | B | exact-injectivity brute force; unnormalised `T` at `q = 7, 11, 19` | see B |
| `audit_census_symbolic.py` | `.out` | C, D | own canonical forms, automorphisms and odd-reversal test; completeness by orbit sizes summing to `2^C(k,2)`; multiplicities as sums of signs over labelled copies; THM-4652 vanishing plus the degree-≤2 elimination, implemented symbolically with exact `Z[q]` coefficients and applied to all `2^C(k,2)` labelled subgraphs | see C, D |
| `audit_interpolation_modular.py` | `.out` | D, G | re-solve of the interpolation systems from the audit's own data; THM-4780 (6) and HYP-9390 (a) against PARI data; algebra by sympy; prime powers with PARI's Satake roots | same unique solutions; no failure at 122 primes (`TT5`) or 70 primes (`TT6`); prime powers all OK |
| `audit_primepower_sums.py` | `.out` | D, G | `S_W4`, `S_(K6−P6)` by brute force over the audit's `GF(27)`, `GF(243)`, and `S_W4` over `GF(343)` | THM-4780 (5) holds at 27, 243 (`TT5`, `TT6`) and 343 (`TT5`) with no modular input; the Frobenius readings hold |
| `audit_pgl2.py` | `.out` | E | exhaustive `PGL(2)`-invariance of `F` on `P¹(F_7)`; random configurations at 11, 13; identity tests on further even-degree graphs | see E |
| `audit_kernel.py` | `.out` | F | brute-force all-maps `U` and injective `T` at `q = 5, 7, 11, 13, 17`; Ahlgren–Ono moment at all odd primes `≤ 199` | all identities hold (see F) |
| `audit_selmer_tt4.py` | `.out` | H | Monsky matrix; four-prime witnesses; class table in both orientation conventions; `|Pf|`; PARI `ellrank`; `TT4` lemma on all regular tournaments of orders 7, 9, 11 (`gentourng`) and on `P_q`; `K4` in Paley graphs | all claims hold (see H) |
| `audit_class_counts.py` | `.out` | I (HYP-9390 Remark 1) | exact class counts of `C3`, the regular 5-tournament `R5` and `TT5` from the symbolic engine, against direct counts in `P_19`, `P_23` | counts correct; odd denominators present |
| `audit_tt7_structure.py`, `audit_tt7_explicit.py` (+ `audit_tt7_structure_helpers.py`) | `.out` | extra | exact reduced form of `2^21 TT7`, checked numerically by brute-forcing all 41 surviving dense sums at `q = 7, 11, 19, 23, 31` | exact agreement |
| `audit_berggren.py` | `.out` | extra (§2) | Berggren depths by the `(m, n)` generator recursion; Collatz `T`-stopping times; sampling to `c ≈ 2·10¹⁰` | lane numbers reproduced; Yao–Knuth "tracking" not supported |
| `make_edits.py` | `make_edits.out`, `required_edits.json` | — | builds and validates the edits | 41 edits valid |

## A. `TT5`, `TT6` by different algorithms

**Source–sink decomposition.**
- A transitive `k`-set has a unique source `s` and sink `t`, and `s → t`. Its middle `k − 2` vertices form a transitive set in `N⁺(s) ∩ N⁻(t)`.
- The affine squares act transitively on the `q(q−1)/2` arcs. Hence `TT_k = q(q−1)/2 · TT_(k−2)(T[M])`, with `M = {x : χ(x) = χ(1−x) = 1}`.
- Inside `T[M]`: `TT3 = Σ_(s→t) (A²)_st` and `TT4 = Σ_(s→t) C((A²)_st, 2)`.

**DFS.** Ordered chains `v1 → … → vk` with all pairs `i < j` arcs, over all start vertices, using integer bitsets.

**Results.**
- No mismatches with `paley_tt5_data.out` (119 primes), `paley_tt6_data.out` (67 primes) or `paley_prime_powers.out` (27, 243, 343, 1331, 2187).
- Among primes below 1500, `TT5(P_q) = 0` exactly at `q = 3, 7, 11` (and `TT5(P_27) = 4212`).
- Among primes up to 751, `TT6(P_q) = 0` exactly at `q = 3, 7, 11, 19, 23`; also `TT6(P_27) = 0`.
- DFS values include `TT6(P_59) = 503034`, `TT6(P_71) = 1724590` and `TT5(P_27) = 4212`.

## B. Graph sums by brute force

- Full enumeration at `q = 7, 11` (W4, `K6 − P6`, octahedron) and at 19 (W4) gives `T_G = q(q−1)S_G` in every case. For example `T_(K6−P6)(7) = 336`, `T_(K6−P6)(11) = −12320` and `T_oct(7) = −2352`.
- `S_W4` at all 44 odd primes `5..199` agrees with the `K5−2K2` column of `paley_graph_sums.out`. It equals `−4 − a_p` with PARI's `a_p`, in both classes mod 4.
- The nine dense even six-vertex sums at `q = 7, 11, 19, 23, 31, 43, 47, 59` agree with `graph_sums6.out` entry by entry. The seven killed sums vanish.
- At `q ≡ 1 (mod 4)` (5, 13, 17, 29), `S_(K6−P6) = |b|² − q² + 21`, as HYP-9276 predicts, and `S_oct = (q − 5)S_W4`.
- Side remark on `paley_graph_sums.py`, harmless because unused:
  - its `K5 − e` column is a normalised sum of an odd-edge graph, so it is not `T_(K5−e)/(q(q−1))` (that ratio is identically 0);
  - the column equals the wheel column by the `S3` symmetry of the Legendre family.

## C. Survivor census and multiplicities

**Classes.**
- Own canonical forms give 11, 34, 156 and 1044 classes for `k = 4..7`. The orbit sizes sum to `2^6, 2^10, 2^15, 2^21`, so the class lists are complete. geng supplied representatives for `k = 7` only.
- Dense / even / killed / surviving counts:
  - `k = 4`: 1 / 1 / 1 / 0;
  - `k = 5`: 3 / 2 / 1 / 1 (W4);
  - `k = 6`: 19 / 9 / 7 / 2 (`K6 − P6`, octahedron);
  - `k = 7`: 150 / 75 / 37 / 38.
- Exactly one of the 38 has all degrees even: `K7 − C7`, with `|Aut| = 14`. This is THM-4780's table.

**Multiplicities** (sums of signs over labelled copies, in the documents' labellings).
- W4 (hub 0, non-edges `{1,2},{3,4}`): `+1`, over 15 copies.
- Octahedron (antipodal pairs `{1,2},{3,4},{5,6}`): `+1`, over 15 copies.
- `K6 − P6` (path `1−…−6`): `−8`, over 360 copies.
- `rev_(K6−P6)(σ) = inv(σ) − des(σ)` holds for all 720 `σ`.
- The Désarménien–Foata formula `Σ_(S_n) sgn(σ) t^(des σ) = (1−t)^⌊n/2⌋ A_⌈n/2⌉(t)` holds for `n = 1..9`. At `n = 6` and `t = −1` the sum is `−16`, so `m = −8`. The attribution is correct (Discrete Math. 99 (1992), Thm 1.1).

**Seven vertices.**
- The magnitudes are `|m_H| ∈ {1, 3, 4, 6, 8, 12, 16, 24, 32, 40, 48}`, as HYP-9390 says. The signs are labelling conventions; they differ between my canonical labellings and geng's.
- HYP-9390's list `{±1, ±3, ±4, ±6, ±8, ±12, ±16, ±24, −32, −40, −48}` contains `−1, −3, −6`. These occur in no labelling used, by the lane or by the audit.
- In geng's labellings the set is `{+1, +3, ±4, +6, ±8, ±12, ±16, ±24, −32, −40, −48}`. Minor.

## D. Reduction lemma and interpolation

**Scrutiny of part 2.**
- Toggles never create loops, since `z ≠ u, w`.
- A doubled edge cancels to `+1` when both factors point the same way, and to `χ(−1) = −1` when they point opposite ways. The text says only "`χ² = 1`". The next sentence about `χ(−1)` covers the sign, but the wording should say so (edit).
- No merges occur in the injective calculus. They occur only in the `U ↔ T` conversion of part 4, where every merged double arc is same-direction (checked for both labellings).
- The degree bound is right. Each eliminated vertex contributes a factor of degree `≤ 1`: degree 0 gives `q − k + 1`; degrees 1 and 2 give constants.
- The statement "`H` runs over dense graphs on fewer vertices" is false for dense `G`, where trivially `H = G` (edit).

**Credit (moderate).**
- THM-4652 contains the vanishing laws and series-parallel polynomiality, not the vertex elimination.
- The degree-≤2 elimination, with subtraction of occupied values, is codex-tiling's `wheel_run_interface_20261014h.md`.
- The ten-edge terminal classification (`K6 − P6` survives; the other ten-edge six-vertex terminals vanish) is tiling-modular-atoms' `paley_next_cores_20261014i.md`.
- THM-4780 (title, part 2 heading, depends_on, prior art) and the results note §4 credit THM-4652 instead.

**Interpolation-free derivation.** The audit's symbolic engine applies THM-4652's vanishing and the elimination to every labelled subgraph of `K_k`, with exact integer polynomials. It gives:
- `2^3 TT3 = q(q−1)(q−3)` and `2^6 TT4 = q(q−1)(q−3)(q−7)`;
- `2^10 TT5 = q(q−1)(q−3)(q²−21q+105) + 1·T_W4`;
- `2^15 TT6 = q(q−1)(q−3)(q³−46q²+682q−3239) + (6q − 100)T_W4 + 1·T_oct − 8T_(K6−P6)`. With `T_oct = (q−5)T_W4` the wheel coefficient becomes `7(q−15)`.

These are exactly THM-4780 (5), and they confirm HYP-9390's split of the wheel coefficient into `q − 5` from the octahedron plus `6q − 100`. So part 5 holds even without the interpolation step.

**Re-solved interpolation.**
- `TT5`, from the audit's own data: 24 equations (`q = 3` and 23 primes `7..199`), rank 7, solution `[0, 315, −483, 192, −25, 1 | 1]`.
- `TT6`: 14 equations (`q = 3` and 13 primes `7..103`), rank 10, also without the `q = 3` row. Solution `[0, −9717, 15002, −6105, 869, −50, 1 | −105, 7 | −8]`.
- Both are identical to the lane's.
- The logic is sound: full column rank, plus the true coefficients satisfying the system, makes the unique solution the truth. This needs the a priori form (verified) and exact data.
- The lane's float64 matrix products in `graph_sums6.py` are exact: every partial sum is an integer of size at most `p⁴ < 2^53`.

**Prime powers without modular input.** The brute-force values below satisfy part 5 against the census, confirming directly that the polynomials are the same for all `q ≡ 3 (mod 4)`:
- `S_W4(27) = −264`, `S_(K6−P6)(27) = 1200`;
- `S_W4(243) = 6960`, `S_(K6−P6)(243) = −3360`;
- `S_W4(343) = 10868`.

## E. PGL(2) lemma (part 4)

**Line by line.**
- Affine maps: correct.
- `ι(x) = −1/x` with no vertex at 0 or ∞: correct.
- `u` at 0: `F(ιx) = (−1)^(in u) F(x)`, re-derived.
- `u` at 0 and `w` at ∞: I get `F(ιx)/F(x) = (−1)^(out_R(w) + in_R(u))`, with `R` the other `k − 2` vertices. This equals the document's `(−1)^(deg'(w)+in'(w)+in'(u))`.
- The document says "every one of these exponents is even", but `deg'` and `in'` are undefined. When `u ~ w`, `deg'(w)` and one of `in'(u)`, `in'(w)` are odd; only the total exponent is even. The proof is correct but needs one more sentence (edit).
- Correct as stated:
  - the parity-orientation lemma;
  - the orbit count `Φ = (q+1)T_(G−u)` for every `u`, so all `T_(G−w)` are equal;
  - `T_G = (q+1−k)T_(G−u)`.
- The octahedron in-degrees are `0, 0, 2, 2, 4, 4`, so `ε = +1`. The relabelling described is the one for deleting `u = 6` (edit: say so).

**Exhaustive invariance on `P¹(F_7)`.**
- Checked `F(gx) = F(x)` over all injective configurations for four maps: `g = x+1`, a non-square dilation, `−1/x`, and a random Möbius map.
- Graphs checked:
  - the octahedron, 20160 configurations;
  - `K7 − C7`, 40320 configurations, with an even in-degree orientation found by search (flip the arc (2,4)).
- Result: 0 failures. Random samples gave the same at 11 and 13.
- Control: flip one arc of the octahedron, so that some in-degrees are odd. Invariance under `−1/x` then fails on 8640 of the 20160 configurations at `q = 7`. It survives at `q = 5, 13`, as the lemma's `q ≡ 1 (mod 4)` clause predicts.

**Identities on further even-degree graphs** (all exact).
- `S_K5 = (q−4)S_K4` at `q = 5, 13, 17, 29, 37, 41`, where both are nonzero.
- `T_C4 = (q−3)T_P3 = −q(q−1)(q−3)` at 7, 11, 13, 17.
- `S_(K7−C7) = (q−6)S_(K6−P6)` (`ε = +1`, `u = 6`) at 7, 11, 19, 23, 13, 17. For example `−2912` at 19 and `18088` at 23.
- At `q ≡ 1 (mod 4)` the lemma forces coincidences, confirmed by the data:
  - `S_(K6−e) = S_(K6−K3)`, from `K7 − K3`;
  - `S_(K6−P4) = S_(K6−(C5+K1))`, from `K7 − C5`;
  - `S_(K6−(C4+K2)) = S_(K6−(C3+P3))`, from `K7 − (C3+C4)`.
  
  At `q ≡ 3 (mod 4)` all of these vanish.

**Prior art.**
- In homogeneous coordinates, with an even in-degree orientation, `F(x) = χ(∏_(i→j) det(x_i, x_j))`. Checked: finite arcs give `χ(−1)` each, an arc into ∞ gives `χ(−1)`, and the total exponent `|E| − out(u)` is even.
- The expression is well defined on `P¹` because all degrees are even, and `GL(2)`-invariant because `|E|` is even.
- This is the classical Möbius invariance of such character products. The triangle case is the oriented Paley two-graph on `P¹(F_q)` (e.g. Gunderson–Semeraro, arXiv:2204.10775; Babai–Cameron on switching classes).
- New to the repo: the counting consequence and its application to the octahedron. The claim "New here: the PGL(2) lemma" should be narrowed (edit).
- The homogeneous form also gives a two-line proof of the lemma.

## F. Kernel-partition conversion

**Brute force at `q = 5, 7, 11, 13, 17`.** All of the following hold. The lane checked 5, 7, 11; 13 and 17 are new.
- `U(W4) = T_W4 − q(q−1)(q−4) = −q(q−1)(q + a_q)`.
- `U(oct) = T_oct + 3(q−4)T_C4 + 3(q)_4 + (q)_3 = q(q−1)((5−q)a_q + 3q)`.
- `T_oct = (q−5)T_W4`.

The algebra from `T_oct = (q−5)(−4 − a_q)q(q−1)` to `(5−q)a_q + 3q` was re-done by hand.

**The wheel moment.**
- The Ahlgren–Ono moment `Σ_b χ(b)t(b)² = −a_p − p` holds at all odd primes `3..199`.
- `U(W4) = q(q−1)·(moment)` exactly, with no collision correction.

**Attribution (minor; edits supplied).**
- HYP-9274 (B) states the injective identity `Σ*(W4) = −q(q−1)(a_q + 4)`.
- The all-maps `U(W4) = −q(q−1)(q + a_q)` is recorded in HYP-9275 (A).
- So two sentences of THM-4780 misattribute: "HYP-9274 (B) states `U(W4) = …`", and part 6's "gives `U(W4)` … plus a collision correction".
- The mathematics is consistent.

**Same orientation as HYP-9275.** HYP-9275's own octahedron (`paley_six_vertex_cores_20261014.py`, edges listed low to high) gives the same signed `S` as THM-4780's at 7, 11 and 19. So the proof applies to HYP-9275's statement as written.

**Legendre form.** `S_W4 = Σ_(x≠0,1) χ(x)a_x² + (q−3)` holds at all primes `5..199`, in both classes. I re-derived it by hand.

## G. HYP-9390

**Algebra and numbers.**
- The passage from THM-4780 (5) to HYP-9390 (a), (6) and the `M_q` form was checked with sympy.
- With the audit's `TT6` and PARI's `a_p`, `|b(p)|²`, the formula holds at all 70 primes `3..751` with `q ≡ 3 (mod 4)`. Part (6) holds at all 122 primes `3..1499`.

**Prime powers.** I used PARI's `a_p` and `b(p)` at `p = 3, 7, 11`, with Satake roots `α_b + β_b = b(p)` and `α_bβ_b = −p²`.

| `q` | `A_q` | `h_r²` | Sym² trace `S_q` |
|---|---|---|---|
| 27 | 260 | 1936 | −1207 |
| 243 | −6964 | 55696 | 3353 |
| 343 | −10872 | 440896 | −323247 |
| 1331 | 90508 | 1926544 | −154983 |
| 2187 | 75092 | 16289296 | −11506327 |

The `TT5` and `TT6` formulas hold at all five, and in every case `|b(q)|² = h_r² = q² − S_q`.

**Sign error (moderate).**
- The title says "`|b(q)|²` the Sym² Frobenius trace minus `q²`". The body says "`−b(q²)` here is HYP-9276's Sym² trace minus `q²`". Both are false.
- The first gives `S_q − q² = −|b(q)|²`. The second, at primes, says `−b(q²) = b(q²) − q²`.
- Correct: `|b(q)|² = q² − S_q`, with `S_p = b(p²)`, so `|b(p)|² − p² = −b(p²)` at Paley primes.
- The formulas and numbers in HYP-9390 use `h_r²` and are right; only the prose is wrong.

**Smaller fixes to (a) and (b).**
- "Equivalent at primes to HYP-9276" holds only at `p ≡ 3 (mod 4)`, against HYP-9276 (1) at those primes. Remark 3 should name HYP-9276 (2) at prime powers. Minor.
- Credit for (b). The Frobenius-trace wheel identity at prime powers was already FINITE-EXACT by direct evaluation at 9, 25, 27, 49, 81, 121, 125 and 169 (`wheel_apery_extension_20261014h.md` §4, codex-tiling). The census adds 243, 343, 1331 and 2187. Minor.
- Worth saying in (b): `a_q` here is not the Fourier coefficient. For example `a(27) = 152`, but the trace is 260.

**Remark 1 is false as stated (moderate).**
- Class counts are `(1/|Aut τ|)` times labelled counts. `|Aut τ|` is odd, and often divisible by 3 or 5.
- Already the 3-cycle count is `q(q²−1)/24`.
- The regular 5-tournament count is `(q⁵ − 5q⁴ + 25q² − 21q)/5120 + T_W4/1024`, checked against direct counts in `P_19` and `P_23` (`audit_class_counts.out`).
- Neither of these lies in `Z[1/2][q] + …`.
- The `S_(K6−P6)` coefficient is `c·q(q−1)`, not a constant in `Z[1/2]`.
- Correct: `Q[q] + Q[q]S_W4 (+ Q·q(q−1)S_(K6−P6))`, with `Z[1/2]` only for labelled counts.

**(c): structure and the explicit seven-vertex form.** The symbolic engine confirms the a priori degrees `deg A7 ≤ 7`, `deg B7 ≤ 2`, `deg C7 ≤ 1`, and gives:

```
2^21 TT7 = A7 + B7·T_W4 + C7·T_(K6−P6) + Σ_(37 genuine H) m_H T_H
A7 = q(q−1)(q−3)(q⁴ − 87q³ + 2785q² − 38827q + 198842)
B7 = 28q² − 1288q + 13238      (octahedron reduced)
C7 = −8(8q − 239)              (K7 − C7 reduced: m = −8, ε = +1)
```

- The full 42-term form matches `2^21 TT7` exactly at `q = 7, 11, 19, 23, 31`, with every surviving dense sum brute-forced.
- For example `TT7(P_31) = 465` and `TT7(P_43) = 7224`. Both agree with the lane's concurrent `tt7/` data.
- The genuine part is not zero: `G(7) = 25536` and `G(31) = 861522240`. So "the 37 sums cancel" can only mean "their combination lies in the span" (edit).
- The lane's `tt7/` span tests use exactly this a priori form. `A7`, `B7` and `C7` above can be subtracted directly.

**Search history (minor).**
- `tt6_search.py` takes `Q(√−3)` CM traces from the first representation `4p = a² + 3b²` and ignores units.
- In weight 3 this gives `−13` and `−37` at `p = 7, 19`. The newform `η(2z)³η(6z)³` (PARI `mfinit([12,3,−3])`) has `2` and `26`.
- So "CM forms of `Q(√−3)` give no identity" is not established. The positive identity is unaffected.

## H. Results note §3 (four vertices, Selmer)

**Monsky matrix.** I used the formula from Monsky's appendix to Heath-Brown (Invent. Math. 118 (1994)), as stated in Watkins–Donnelly–Elkies–Fisher–Granville–Rogers §4.1:
- `A_ij = [p_j/p_i]` for `i ≠ j`, and `A_ii = Σ_(j≠i) A_ij`;
- `D_u = diag([u/p_i])`;
- `M = [[A + D_2, D_2], [D_2, A + D_−2]]`, with `dim Sel_2(E_n) − 2 = nullity(M)`.

Sanity check: the parity of the nullity equals the root-number parity (odd iff `n ≡ 5, 7 (mod 8)`) for all 1214 odd squarefree `n < 3000`.

**Witnesses.** All four reproduce:

| primes | class | nullity |
|---|---|---|
| `(3, 11, 59, 83)` | TT4 | 0 |
| `(3, 11, 19, 107)` | source + C3 | 2 |
| `(3, 11, 19, 163)` | sink + C3 | 2 |
| `(3, 11, 19, 43)` | strong | 0 |

- For all 64 labelled 4-tournaments, the nullities of `M` and of `A² + A + I` depend only on the class (`0, 2, 2, 0`), in both orientation conventions.
- `|Pf| = 3` holds exactly on the two classes with one cyclic triple (THM-1455).
- Orientation convention (minor; edit). The note orients the Legendre tournament as `p_i → p_j` iff `(p_j/p_i) = 1`. This is the converse of `bsd_selmer_growth_20261013a.md` §3, which the note cites. The values are invariant, but the labels of the two middle witnesses swap.
- Extra, from PARI 2-descent (`ellrank`):
  - `n = 67089` has `[r1, r2, s] = [0, 0, 2]`: rank 0, with the Selmer excess in Sha;
  - `n = 102201` has `[2, 2, 0]`: rank 2.
  
  So Selmer dimension 2 is not a rank statement. The note correctly says "Selmer".

**`TT4` lemma and class counts.**
- Checked exhaustively on all 3, 15 and 1223 regular tournaments of orders 7, 9 and 11.
- There are no violations. Equality holds iff the tournament is doubly regular, and every class-count relation holds.
- At `n = 9 ≡ 1 (mod 4)` equality is impossible: the minimum is 18 against a bound of 13.5.
- On `P_q` (`q ≤ 43`) the `TT4` value and the class counts follow. For example `P_43` gives 40635 / 16555 / 16555 / 49665.
- The lemma is classical in substance (edit: cite it). For regular `T`:
  - `TT4 = n·C(k,3) − Σ_v c3(N⁺(v))`;
  - `c4 = n(n²−1)(n−3)/48 − Σ_v c3(N⁺(v))`.
  
  So it is the known bound `c4(T) ≥ c4(DR_n)`, with the same equality case (Savchenko, arXiv:2403.07629, §1 and references).

**`K4` and EPS.**
- `K4(G_q) = q(q−1)(q² − 20q + 81 + a_q(η(4z)⁶))/1536`, using `η(4z)⁶` coefficients from PARI's 16.3 form.
- It holds at all 21 primes `q ≡ 1 (mod 4)` up to 200, against brute-force `K4` counts.
- So does EPS's `q(q−1)((q−9)² − 4y²)/1536`, where `q = x² + y²` with `y` even. The equivalence is the identity `a_q = 2x² − 2y²`.
- `S_K4 = a_q(η(4z)⁶)` at `q ≡ 1 (mod 4)` and `S_K4 = 0` at `q ≡ 3 (mod 4)`, for all primes below 200.

## I. Prior art, statuses, and the §2 and §5 statements

**Novelty.**
- No published order-5 or order-6 count of transitive subtournaments of Paley tournaments was found.
- McCarthy–Springfield (arXiv:2311.02135) give orders 3 and 4 only. Cor. 2.3 and 2.7 are confirmed, as is the weight-3 level-32 remark after Cor. 2.4 for `k = 4`, `q ≡ 5 (mod 8)`.
- The `TT3` and `TT4` values themselves hold for all regular, respectively doubly regular, tournaments (classical).
- The PGL(2) mechanism is classical (see E).
- The `k ≤ 6` survivor census was largely in the 2026-10-14 packets and HYP-9275. THM-4780 claims only `k = 7` as new, which is fair.

**Statuses.**
- THM-4780 "PROVED" is honest: part 5 is a computer-assisted exact argument, now confirmed by two independent routes.
- HYP-9390 "OPEN, CONDITIONAL" is honest.
- HYP-9390 (c) "PROVED (exhaustive)" is honest for the list, though the signs are conventions.
- The results note should mark three things conditional on HYP-9276: the summary row "`TT6` = the level-32 adjoint", the `k = 6` table row, and the loop-expansion consequence (edits).

**§2 (minor).**
- "`c` is the Lorentz-invariant height" is wrong: the form `a² + b² − c²` is invariant, not `c`.
- "The mean tracks `(6/π²)(log √c)²`" is not supported. The ratio is 0.98, 0.92 and 0.88 in the note's windows, and 0.86, 0.82 and 0.80 near `c = 1.5·10⁶, 10⁸, 10¹⁰`.
- "Collatz … no quadratic form is preserved" misquotes the cited result:
  - `ternary_berggren_20260925.md` §6 proves only that the integer `J = tr²/det` blocks a fixed projective conjugacy, and says it does not distinguish convergence;
  - a single rescaled Collatz affine word is hyperbolic and does preserve an indefinite binary form.
- The depth and stopping-time numbers themselves are reproduced exactly.

**§5.**
- CS ⟺ triangle for inner-product norms is correct, and any `t ≠ 0` works, not just small `t`.
- "The `ℓ¹` norm … triangle inequality without CS" is false for the dot product, since `|x·y| ≤ ‖x‖₁‖y‖₁`. The intended point is that `ℓ¹` comes from no inner product (edit).
- FTA as FTC for `d log p` on the spiral is correct. So is "the strip returns" (`Z → Z/2`, quotient by `2Z`).
- The Möbius-band sentence conflates the real and complex pictures. On the real line nothing is exchanged at a sign change; the exchange is the monodromy around a simple zero (edit).
- CLT/Tychonoff and the Rule 30 references (THM-3456 and THM-4210 exist) are fine.

**§7 and sources (minor).**
- "Every `.py` there has its `.out` beside it" is false for four TT5 search-history scripts. The reproduction table already marks them "—".
- "Sanchez-Flores (Graphs Combin. 1996)" should be Graphs Combin. 10 (1994) 367–376 (`v(6) = 28`, ST27) and 14 (1998) 181–200 (uniqueness).
- ST27 is the quadratic-residue ("Galois") tournament on `GF(27)` (Neiman–Mackey–Heule), consistent with `TT6(P_27) = 0`.

## Errors and overclaims, by severity

No major errors: none invalidates a stated theorem or numerical claim.

**Moderate.**
1. HYP-9390, title and (a) body: the Sym² reading of `|b(q)|²` has the wrong sign. It should read `|b(q)|² = q² − S_q` (G).
2. HYP-9390 Remark 1: wrong coefficient ring. Class counts have odd denominators `|Aut τ|`, and the `K6 − P6` coefficient is `c·q(q−1)` (G).
3. Credit: the degree-≤2 elimination calculus is attributed to THM-4652 in THM-4780 (title, part 2, depends_on, prior art) and in the note §4. It is from `wheel_run_interface_20261014h.md` and `paley_next_cores_20261014i.md` (D).
4. Novelty overclaim: "New here: the PGL(2) lemma". The Möbius invariance of even-degree χ-products is classical; the counting consequence and its application are new (E).

**Minor.**
- THM-4780:
  - `U(W4)` misattributed to HYP-9274 (B), twice;
  - part 2 wording: dense `G`, and opposite-direction doubled edges;
  - part 4: undefined `deg'` and `in'`, "every exponent even", and which `u` is deleted;
  - part 7 and the `TT3`/`TT4` values are classical and uncited;
  - earlier prime-power wheel checks are uncredited.
- HYP-9390:
  - restrict "equivalent at primes" to `p ≡ 3 (mod 4)`;
  - Remark 3 at prime powers;
  - (b): the credit, and "Frobenius trace, not Fourier coefficient";
  - (c): the multiplicity sign list, "the 37 sums cancel", and `K7 − C7 → ±(q−6)T_(K6−P6)`;
  - the `Q(√−3)` negative result.
- Results note:
  - conditionality unmarked (§0, §4 table, §4 consequence);
  - §2: three statements;
  - §3.3: orientation convention;
  - §4: credit, and "new PGL(2) lemma";
  - §5: the `ℓ¹` and Möbius-band sentences;
  - §7: the reproduction claim, and the Sanchez-Flores year.

## Optional improvements (not in `required_edits.json`)

- Give part 4 the two-line homogeneous-coordinate proof (E).
- The zero lists are complete. Every tournament on at least 14 vertices contains a `TT5` (Reid–Parker 1970), and every tournament on at least 28 contains a `TT6` (Sanchez-Flores 1994). So:
  - `TT5(P_q) = 0` iff `q ∈ {3, 7, 11}`;
  - `TT6(P_q) = 0` iff `q ∈ {3, 7, 11, 19, 23, 27}`, i.e. every Paley tournament below 28 is `TT6`-free.
- Add the exact reduced form of `2^21 TT7` (G(c)) to HYP-9390 (c). The 37 genuine sums are listed, with multiplicities, in `audit_tt7_explicit.out`.
- Record the `q ≡ 1 (mod 4)` coincidences among six-vertex sums forced by the PGL(2) lemma (E).

## Reproduction

Run from this directory, in this order:
1. `gp -q modular_coeffs.gp < /dev/null > modular_coeffs_gp.out`.
2. Any `audit_*.py`; each writes its `.out`.
3. `python3 make_edits.py`.

Total running time is under 3 minutes on the mac-mini. Peak memory is a few hundred MB (the `GF(2187)` difference tables).

Sources:
- McCarthy–Springfield, arXiv:2311.02135.
- Neiman–Mackey–Heule, arXiv:2011.00683.
- Savchenko, arXiv:2403.07629.
- Désarménien–Foata, Discrete Math. 99 (1992); arXiv:1812.09098.
- Gunderson–Semeraro, arXiv:2204.10775.
- Sanchez-Flores, Graphs Combin. 10 (1994) and 14 (1998): doi 10.1007/bf02986687 and 10.1007/s003730050025.
