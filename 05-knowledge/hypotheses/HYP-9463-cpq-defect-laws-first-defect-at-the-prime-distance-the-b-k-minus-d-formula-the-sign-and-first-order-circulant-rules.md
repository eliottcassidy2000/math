---
id: HYP-9463
title: "After THM-4871 (the (p, q, 1^k) defect family of HYP-9340): (1) first-defect law: for coprime odd prime powers p != q and H = C_pq of type (p, q, 1^k), tau_T(H) - tau_E(H) = 0 for k < |p - q| and = 1 at k = |p - q| (62 census rows over thirteen families; PROVED for primes p < q < 2p by THM-4871 statement 7, numerical for q > 2p and for prime powers); (2) the defect formula tau_T(k) - tau_E(k) = b_(k-d)(4) - E_T^F(k - d), d = |p - q|, where b_j(4) counts 4-coloured tournaments on j points and E_T^F(j) the tournament coincidences created by the fixed points (true in every computed row, giving 1, 4, 16, 84, 672 for C15 and 1, 4, 16, 88, 704 for (5,7)); (3) tau_T(H) >= tau_E(H) for every odd-order H (HYP-9344 (1), now with 64 further rows); (4) the uniform first-order rule question for prime circulants (HYP-9344 (3)) in its exact form, with the quantifier-free reduction to the Codex degree bound; (5) the fold reading of the owner's mirror remarks is typed, not conjectured"
status: >
  (1) PROVED for primes p < q < 2p (THM-4871 statement 7: 0 below k = q - p, exactly 1 at k = q - p); NUMERICAL for
  q > 2p and for the prime-power families (all 62 rows of THM-4871 statement 5). (2) NUMERICAL
  (every computed row, including the rows with E_T^F > 0: C15 k = 5, 6 and the prime-power rows where the non-CI
  coincidences are excluded from E_T^F); for primes E_T^F = E_T by THM-4871 statement 3. (3) NUMERICAL. (4) OPEN;
  the quantifier-free homogeneous fragment is excluded by the Codex Boolean-degree theorem (CITED). (5) typed in the
  note (section 5); no claim. (4)'s pre-registered model check is FINITE-EXACT (section 4 (d)). Audit A (2026-11-04,
  PASS WITH CORRECTIONS) applied: statement 1's first-defect description restricted to the proved range, statement
  2's "equivalent form" and "indexing" reading corrected, row counts (MISTAKE-710). Pre-registered decisive tests below.
session: mac-mini-2026-11-04-oddprimes
source: 05-knowledge/results/cpq_defect_family_mirror_folds_20261104.md
scripts:
  - 04-computation/experiments/cpq_defect_20261104/cpq_defect_census.py (+ cpq_<p>_<q>_k<k>.out), witness_examples.py
related:
  - 01-canon/theorems/THM-4871-cpq-defect-family-of-hyp-9340-normaliser-orbit-parametrisation-disjoint-orbit-lemma-for-odd-order-groups-first-defect-at-the-prime-distance-and-the-census-refuting-hyp-9344-2.md
  - 01-canon/theorems/THM-4728-hyp-9340-beyond-p-groups-defect-identity-pronormal-frattini-lemma-and-the-smallest-failure-c15-at-n-10.md
  - 05-knowledge/hypotheses/HYP-9344-after-the-sylow-day-defect-growth-of-c15-tau-t-at-least-tau-e-affine-and-first-order-rules-for-prime-circulants-and-prime-power-circulants.md (statement 2 REFUTED by THM-4871: 84, not 88)
  - 05-knowledge/hypotheses/HYP-9340-refined-equinumerosity-for-every-odd-order-symmetry-and-what-a-tournaments-even-graphs-bijection-can-be.md
  - 05-knowledge/results/affine_circulant_ideal_obstruction_20261020, circulant_boolean_degree_20261020 (Codex, CITED)
---

# HYP-9463 — what the (p, q, 1^k) census leaves open

Notation as in THM-4871: `d = |p - q|`, `b_j(4) = 1, 4, 16, 88, 704, 8912, 188288, ...` (4-coloured tournaments on `j`
points), `E_T^F(j)` = the part of the tournament multiplicity excess at level `j` due to conjugates of `H` whose `p`- or
`q`-orbit lies among the fixed points (all of `E_T` when `p`, `q` are primes, THM-4871 statement 3; for `p = 9` the excess
`2 t_q b_j(4)` coming from the isomorphic non-dilation-equivalent circulant 9-tournaments is excluded).

## Statements

1. **First-defect law.** `tau_T(H) - tau_E(H) = 0` for `k < d` and `= 1` for `k = d`. For primes `p < q < 2p` this is
   THM-4871 statement 7, and there the unique witness at `k = d` is `K_(q,q)` on the larger orbit against the smaller
   orbit plus `d` fixed points, with `E_even(d) = E_wit(d) = E_T(d) = 0` (`(3,5) k = 2`, `(5,7) k = 2`, `(7,11) k = 4`,
   `(11,13) k = 2`). In the families with `d >= p` or a prime-power orbit the same values `0, ..., 0, 1` occur, but the
   first defect `1` is a net value of many witnesses and coincidences: `(3,7) k = 4`: `W = 17` (one `K_(7,7)` and
   sixteen `K_(3,3)`-type configurations), `E_even = 48`, `E_T = 64`; `(9,5) k = 4`: `W = 1409`, `E_T = 1408`;
   `(9,7) k = 2`: `65 / 64`; `(9,11) k = 2`: `129 / 128`; `(9,13) k = 4`: `8449 / 8448`; `(3,11)`, `(3,13)`, `(7,13)`
   are zero through the computed range (`k <= 4 < d`). Open exactly there: `q > 2p` and prime powers.
2. **Defect formula.** `tau_T(k) - tau_E(k) = b_(k-d)(4) - E_T^F(k - d)` for `k >= d` (and `0` below). Rows with a
   nonzero correction: C15 `k = 5`: `88 - E_T(3) = 88 - 4 = 84`; `k = 6`: `704 - E_T(4) = 704 - 32 = 672`. Rows without:
   `(5,7)` `k = 2..6`: `1, 4, 16, 88, 704` (`E_T(j) = 0` for `j <= 4 < p`); `(3,7)` `k = 4, 5, 6`: `1, 4, 16`; `(9,7)`,
   `(9,11)` `k = 2, 3`: `1, 4`; `(9,5)`, `(9,13)` `k = 4`: `1`. Equivalent form: `tau_E(k) = tau_T(k) - R(k - d)/(2 t_p t_q)
   + E_T^F(k - d)`: the defect is the number of 4-coloured tournaments on the `k - d` non-padding fixed points minus the
   fixed-point coincidences at that level, with no factor `t_p t_q` (the `(5,7)` and `(9,q)` families have
   `t_p t_q > 1` and the same `1, 4, 16, ...`). The formula is supported only as a net identity: the witnesses are not
   indexed by those coloured tournaments (`(5,7)`: `W(3) = 2 < b_1(4) = 4`, `W(4) = 6 < 16`, and the defects `4, 16, 88`
   there are `W + E_even - E_T = 2 + 2, 6 + 10, 24 + 72 - 8`). A proof would have to match the `K_(p,p)`-type witnesses,
   the even-graph coincidences and the tournament coincidences bijectively, as the `k = p = 3` case does explicitly
   (`W = E_even = 2 t_q`, `E_T = 4 t_q`) and as THM-4871 statement 7 does for `k <= q - p`, `q < 2p`.
3. **Sign.** `tau_T(H) >= tau_E(H)` for every odd-order `H <= S_n` (HYP-9344 (1)); all 62 rows here satisfy it,
   including the prime-power rows where `E_T` is large (`(9,13)` at `k = 4`: `E_T = 8448`, `W = 8449`).
4. **Uniform first-order rules for prime circulants (OPEN).** Precisely: is there a first-order formula `phi(x)` in
   the language `{+, ., 0, 1, S}` (`S` a unary predicate) such that for every prime `p` and every skew set `S <= F_p^*`
   the set `C_phi(S) = {x : (F_p, S) |= phi(x)}` is a symmetric set and `S |-> C_phi(S)` induces a bijection from the
   dilation classes of skew sets onto the dilation classes of multiplier-even symmetric sets? Observations:
   (a) a quantifier-free `phi` whose atoms are `S(c_i x)` (`c_i in Z`, no constant terms) is dilation-equivariant and each
   output bit is a Boolean function of at most `r = #atoms` input bits, so its Boolean degree is `<= r`; the Codex
   degree theorem (`circulant_boolean_degree_20261020`: every dilation-equivariant type bijection has degree
   `>= ceil((e - 1)/a)`, `e = 2^a` the dyadic part of `(p - 1)/2`, unbounded along the primes dividing Fermat numbers)
   excludes every such `phi` uniformly; (b) atoms `S(c x + d)` with `d != 0`, atoms in `x^2`, and quantifiers break
   equivariance, and then the well-definedness of the induced map on classes is itself a first-order sentence
   `forall a exists b forall x (phi^(aS)(x) <-> phi^S(b^(-1) x))` that every candidate must satisfy for all `p` and all
   `S`, a uniform constraint that the quantifier-free case shows to be severe; (c) the `S`-free part of any candidate
   must send the Paley skew set (`p = 3 mod 4`) to a ring-definable multiplier-even symmetric set. HYP-9341 (c)
   conjectures that no such `phi` exists; the right tool is still missing (definable sets in pseudofinite fields with a
   predicate are not constrained by the Chatzidakis-van den Dries-Macintyre counting).
   (d) **Pre-registered model check (FINITE-EXACT; `fo_rule_modelcheck.py`, `fo_rule_modelcheck3.py`).** Among the
   two-atom quantifier-free rules `B(S(c_1 x + d_1), S(c_2 x + d_2))` with `|c_i|, |d_i| <= 4` (42,048 rules) and
   `<= 8` (594,048 rules), and the three-atom rules with `|c_i|, |d_i| <= 4` (16,594,944 rules), **none** is a type
   bijection at `p = 13` or at `p = 17` (the first primes with `(p - 1)/2` even and greater than `2`); at `p = 7` and
   `11` exactly the XOR rules `S(c_1 x) + S(c_2 x)` with `c_2/c_1` in the right class pass (10 and 12 of the box-4
   two-atom rules, two of them at both primes), as the Codex affine theorem predicts (`m = (p - 1)/2` odd), and at
   `p = 5` (`t_5 = e_5 = 1`) the test is vacuous. So no uniform rule with at most three atoms in these boxes exists.
5. **Fold reading (typed, no claim).** The owner's "same value twice on either side of the peak, small changes on one
   side are large on the other" is the saddle-node mechanism: `x^(1/x)` at `e` (mirror bases `r_n`, derivative
   `m'(x) = f'(x)/f'(m(x)) -> -infty` as `x -> 1+`), the unit-ball volume `V(x) = pi^(x/2)/Gamma(x/2 + 1)` at
   `x* = 5.2569` (mirror dimensions `0 <-> 12.764`, `1 <-> 10.814`, `2 <-> 9.205`, `3 <-> 7.820`, `4 <-> 6.603`,
   `5 <-> 5.517`, no two integer dimensions with equal volume), and the period-3 saddle-node of `x^2 + c` at `c = -7/4`
   (two real 3-cycles merge; the parabolic cycle generates `Q(zeta_7)^+`). The number `5` shared by the ball peak and
   the `8 + 4` self-converse split of the twelve 5-tournaments is a coincidence of two unrelated unimodal profiles
   (NUMEROLOGY). See the note, section 5.

## Pre-registered decisive tests

- Statement 1 at `(7,11)` (`d = 4`): run after pre-registration, `0, 0, 0, 0, 1` with `W = 1` at `k = 4` (now the
  proved case, THM-4871 statement 7); `(5,11)` (`d = 6`, `k <= 3` predicted `0`) did not run (the chain stalled behind
  the `(3,25) k = 2` job) and remains pre-registered; the first open cases of the law are `(3,7) k <= 3` -> `0` (done)
  and the prime-power rows, all consistent.
- Statement 2 at C15 `k = 7` (`n = 15`, `2 b_7(4) = 13,694,400` objects per side): predicted defect
  `b_5(4) - E_T(5) = 8912 - 260 = 8652`; at `(3,7)` `k = 7`: `b_3(4) - E_T(3) = 88 - 8 = 80`; at `(5,7)` `k = 7`:
  `b_5(4) - E_T(5) = 8912 - 8 = 8904`. Memory on the 8 GB machine, not time, is the obstacle (hash the canonical
  forms, stream the coloured types).
- Statement 3: a counterexample needs a type with many `Aut`-classes of `H` and few witnesses; the natural place is
  `k >= q` (alternative `q`-orbits inside the fixed points), first reachable at `(3,5)` `k = 5` where `E_wit` first
  appears (`2`), and the sign still holds there by `84`.
- Statement 4: model-check all quantifier-free rules with two atoms of the forms `S(c x)`, `S(c x + d)`, `S(x^2 + d)`
  at `p = 13, 17` for the well-definedness sentence.
