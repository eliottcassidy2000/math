---
id: HYP-9463
title: "After THM-4871 (the (p, q, 1^k) defect family of HYP-9340): (1) first-defect law: for coprime odd prime powers p != q and H = C_pq of type (p, q, 1^k), tau_T(H) - tau_E(H) = 0 for k < |p - q| and = 1 at k = |p - q| (64 census rows over eleven families; proved only as >= 1 for primes p < q < 2p); (2) the defect formula tau_T(k) - tau_E(k) = b_(k-d)(4) - E_T^F(k - d), d = |p - q|, where b_j(4) counts 4-coloured tournaments on j points and E_T^F(j) the tournament coincidences created by the fixed points (true in every computed row, giving 1, 4, 16, 84, 672 for C15 and 1, 4, 16, 88, 704 for (5,7)); (3) tau_T(H) >= tau_E(H) for every odd-order H (HYP-9344 (1), now with 64 further rows); (4) the uniform first-order rule question for prime circulants (HYP-9344 (3)) in its exact form, with the quantifier-free reduction to the Codex degree bound; (5) the fold reading of the owner's mirror remarks is typed, not conjectured"
status: >
  (1) NUMERICAL (all 64 rows of THM-4871 statement 5; the inequality >= 1 at k = q - p is PROVED for primes
  p < q < 2p, THM-4871 statement 4; = 0 for k < |p - q| is proved only insofar as E_T = 0 for k < p). (2) NUMERICAL
  (every computed row, including the rows with E_T^F > 0: C15 k = 5, 6 and the prime-power rows where the non-CI
  coincidences are excluded from E_T^F); for primes E_T^F = E_T by THM-4871 statement 3. (3) NUMERICAL. (4) OPEN;
  the quantifier-free homogeneous fragment is excluded by the Codex Boolean-degree theorem (CITED). (5) typed in the
  note (section 5); no claim. Pre-registered decisive tests below.
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

1. **First-defect law.** `tau_T(H) - tau_E(H) = 0` for `k < d` and `= 1` for `k = d`; the unique witness at `k = d` is
   `K_(q,q)` on the larger orbit against the smaller orbit plus `d` fixed points (THM-4871 statement 4(a)), with
   `E_even(d) = E_wit(d) = 0` and `E_T(d) = E_T^F(d) = 0`. Evidence: `(3,5)`, `(5,7)`, `(9,7)`, `(9,11)` at `k = 2`;
   `(3,7)`, `(9,5)`, `(9,13)` at `k = 4`; `(3,11)`, `(3,13)`, `(7,13)` zero through the computed range (`k <= 4 < d`).
2. **Defect formula.** `tau_T(k) - tau_E(k) = b_(k-d)(4) - E_T^F(k - d)` for `k >= d` (and `0` below). Rows with a
   nonzero correction: C15 `k = 5`: `88 - E_T(3) = 88 - 4 = 84`; `k = 6`: `704 - E_T(4) = 704 - 32 = 672`. Rows without:
   `(5,7)` `k = 2..6`: `1, 4, 16, 88, 704` (`E_T(j) = 0` for `j <= 4 < p`); `(3,7)` `k = 4, 5, 6`: `1, 4, 16`; `(9,7)`,
   `(9,11)` `k = 2, 3`: `1, 4`; `(9,5)`, `(9,13)` `k = 4`: `1`. Equivalent form: `tau_E(k) = tau_T(k) - tau_T(k - d) + b_(k-d)(4)`
   whenever all coincidences at level `k - d` come from the fixed points. The reading is that the `K_(q,q)` witnesses
   at level `k` are indexed by the `H`-symmetric *tournaments* at level `k - d` with the cross orientation fixed (the
   `d` padding points are forced, the remaining `k - d` points carry a 4-coloured structure), and the `K_(p,p)`-type
   witnesses, the even-graph coincidences and the tournament coincidences cancel among themselves; a proof would have to
   match the three cancelling families bijectively, as the `k = p = 3` case does explicitly (`W = E_even = 2 t_q`,
   `E_T = 4 t_q`).
3. **Sign.** `tau_T(H) >= tau_E(H)` for every odd-order `H <= S_n` (HYP-9344 (1)); all 64 rows here satisfy it,
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
5. **Fold reading (typed, no claim).** The owner's "same value twice on either side of the peak, small changes on one
   side are large on the other" is the saddle-node mechanism: `x^(1/x)` at `e` (mirror bases `r_n`, derivative
   `m'(x) = f'(x)/f'(m(x)) -> -infty` as `x -> 1+`), the unit-ball volume `V(x) = pi^(x/2)/Gamma(x/2 + 1)` at
   `x* = 5.2569` (mirror dimensions `0 <-> 12.764`, `1 <-> 10.814`, `2 <-> 9.205`, `3 <-> 7.820`, `4 <-> 6.603`,
   `5 <-> 5.517`, no two integer dimensions with equal volume), and the period-3 saddle-node of `x^2 + c` at `c = -7/4`
   (two real 3-cycles merge; the parabolic cycle generates `Q(zeta_7)^+`). The number `5` shared by the ball peak and
   the `8 + 4` self-converse split of the twelve 5-tournaments is a coincidence of two unrelated unimodal profiles
   (NUMEROLOGY). See the note, section 5.

## Pre-registered decisive tests

- Statement 1 at `(7,11)` and `(5,11)` (`d = 4`, `6`): the runs at `k <= 3` must show defect `0` (recorded in the note
  when they finish); the first nonzero value `1` at `k = 4` for `(7,11)` is a `2^... `-object run at `n = 22` reachable
  with the same tool.
- Statement 2 at C15 `k = 7` (`n = 15`, `2 b_7(4) = 13,694,400` objects per side): predicted defect
  `b_5(4) - E_T(5) = 8912 - 260 = 8652`; at `(3,7)` `k = 7`: `b_3(4) - E_T(3) = 88 - 8 = 80`; at `(5,7)` `k = 7`:
  `b_5(4) - E_T(5) = 8912 - 8 = 8904`. Memory on the 8 GB machine, not time, is the obstacle (hash the canonical
  forms, stream the coloured types).
- Statement 3: a counterexample needs a type with many `Aut`-classes of `H` and few witnesses; the natural place is
  `k >= q` (alternative `q`-orbits inside the fixed points), first reachable at `(3,5)` `k = 5` where `E_wit` first
  appears (`2`), and the sign still holds there by `84`.
- Statement 4: model-check all quantifier-free rules with two atoms of the forms `S(c x)`, `S(c x + d)`, `S(x^2 + d)`
  at `p = 13, 17` for the well-definedness sentence.
