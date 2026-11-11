---
id: THM-5129
title: "One-vertex sources of the Hamiltonian-path spectrum (rooted tournaments (B, y), offset D = H(B), slope U = U_B(y), closed values D + 2^t U). (A) Finite supply: for every odd D the slopes of offset D, over bases of every order, are Slopes_0(D) * 2^j (j >= 0) with Slopes_0(D) finite and computable from the strong tournaments of order at most n(D). (B) Complete lists for every offset D <= 73 and for D = 129, 257; at the Fermat-type offsets 3, 5, 9, 17, 33, 65, 129, 257 the slopes have 1, 5, 8, 11, 34, 66, 177, 359 odd parts, and at level 6 of the 9-grid the least exponent is 5, above the bound 3 given by U >= D (the KNOWN bound j >= s - 3 is attained by the known sources at levels 3, 4, 5, 7, 8). (C) The two-step / near-join values D + 2^t U, t >= 1, omit 7, 11, 13, 21, 47, 67 at every order and attain every other odd number in [3, 225); 11, 13, 47, 67 are head values D - E + U. (D) No window [h, 3h/2] lies in the slope set at the offsets computed, except h = 1, 5: attained h <= 73 (mean covered share 0.50 on 37 <= h <= 73) and h = 129, 257 (shares 0.75, 0.72); whether windows fill up for large h is OPEN. (E) Head correction E = U_odd - U_even. (F) Line counts: every odd N in [129, 16384) has at least N/64 sources on its line D + 2U = N (rooted order <= 10; pre-registered; a count on a finite range)."
status: >
  (A) PROVED. KNOWN in (A): the product law (A1) is THM-5130 (F2)(b) (kind-pasteur-2026-11-10-S22) for a one-vertex
  set O = {y}, U_B(O) = U_{C_a}(O) prod_{i<a} H(C_i) prod_{i>a} Pi_2(C_i) (the factorisation of a closing coefficient
  over an ordered join goes back to Session216, star/PROOF.md, Theorem 3); it is re-proved here in three lines
  because (A2) and (B)-(D) depend on it. KNOWN ingredients of (A2): the multiplicativity of H over strong
  components, the insertion lemma H(T) >= H(T - v) and the monotonicity of the strong floor f (THM-1370 (1)-(3));
  the finiteness of the strong tournaments with at most D paths (classical: Moon 1972; H(S) >= |S| from Camion's
  Hamiltonian cycle already suffices); the finite inventory of nontrivial strong components at a fixed path count,
  valid at every order (Session 24, sections 1 and 5 of
  05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md, complete for
  H <= 73; THM-4744 (3) for every count below 625; THM-4884 (2)), and all-order classifications obtained from it
  (chanpak_55_105_golden_20261110/targets/PROOF.md, P1-P4; decomposition16_odd_spectrum_20261110, J1; THM-5130, a
  seed never has slope 7); the use of the strong floor to decide a two-step question at every order from the
  classes of order <= 8 (THM-4873 (iii): marks of digit 1 and 2 on the lines 3h + 2 sigma). NEW in (A): the rooted
  form (A2) (the slopes of a fixed offset are finitely many odd parts times powers of two) and its counting
  corollary (A3); the step f(n) >= f(n-1) + 2 from the odd-cycle formula is added only to keep the proof
  self-contained.
  (B), (C), (D): FINITE-EXACT computations over the strong tournaments of order <= 8 (source_supply.py), valid at
  EVERY order of the base by (A) and f(9) = 75; cross-checked against the complete rooted census of order <= 9
  (3,239 pairs with D <= 73, all in the product-law set); the correlations quoted in (D) are NUMERICAL. (B'):
  FINITE-EXACT, from an exhaustive enumeration of the strong tournaments of orders 10 and 11 with at most 257
  paths (source_supply_levels78.py); these tournaments are in the KNOWN inventory of THM-4744 (3), which also has
  f(11) = 225 and f(12) = 375, and were enumerated again because U and Pi_2 are needed. (B') was added after the
  audit was launched and was recomputed by the audit (audit-5129, after seeing the committed output, so not blind).
  KNOWN before: the exceptional set 1, 7, 11, 13, 21, 47, 67 as a FINITE-EXACT fact through base order 8 for the
  near-join family (HYP-9453 (1), klein-2026-11-02) and through base orders 9 and 10 for the two-step family
  (THM-4865 (e); HYP-9454 (6)), and its conjectured permanence (HYP-9453 (1)); the attainment of every odd number in
  [69, 382687] by the near-join values of bases of order <= 9, and in [69, 6171561] of order <= 10 (THM-4865 (e));
  the size of the all-vertex digit sets (about 857 members on average for h in [300, 1000), THM-4865 (e)); the laws
  H(B + z) = D - E + U and D + 2U (THM-4865); the bound U >= D (HYP-9453; THM-5130 (F1)) and its consequence
  j >= s - 3 for a source (2^s + 1, 9 * 2^j) (kind-pasteur, SOURCES_AND_DOMINATION.md, Proposition 4:
  v >= s - floor(log2 d)); the remark that the sets of responses of a fixed h are not intervals for h <= 31
  (HYP-9453 (2)); the multiplicities m(h) of the strong classes (THM-4884 (2): m(h) = 2 exactly at 9, 17, 19, 35,
  39, 47, 63); pure cones from slopes that are powers of two (THM-5109 (C)); the expectation that the number of
  ready pairs on a line 3H + 2 sigma = v grows linearly, for the sole-exit sweep (HYP-9454 (1), klein-2026-11-04;
  a heuristic, OPEN); one-vertex sources of the 9-grid at levels 2-6 (Session261) and at levels 7, 8, 9, 11
  (kind-pasteur-2026-11-10-S22, SOURCES_AND_DOMINATION.md, pushed before this file); the table of levels with a
  one-vertex source of order <= 9 for d <= 129 (ibid.); U_B(o) >= 2^outdeg(o) (ibid.); the word "qualifying",
  defined first for the permutations of a k-vertex extension
  (05-knowledge/results/odd_c_seed_symmetry_20261110/EXTENSIONS.md) and read as a cover count of a one-vertex
  extension by kind-pasteur (kps_s22 README, section 4); (E) uses it in a narrower sense.
  NEW here: permanence of the omissions 11, 13, 47, 67 at every order (C) (7 and 21 are not path counts at all,
  THM-1370); the slope lists (B), (B') and their completeness at every order (of the least exponents of the 9-grid
  only the one at level 6 needs the lists: at levels 3, 4, 5, 7, 8 the known sources attain the KNOWN bound
  j >= s - 3, and at level 2 the known source has j = 0); the exact window counts (D); the line counts (F) with
  four PRE-REGISTERED tests (L1-L4, commit e06b161113 precedes the run), all CONFIRMED (L1, L2, L4 one-sided, L3
  already true in sample; passed with wide margins, see (F) 3). (E) PROVED: a two-line corollary of the
  odd-cycle formula (THM-002) and the pendant law (THM-4865); it is the one-vertex case of Session261's marked
  odd-cycle identity (bridges/EXTENSION4.md, eq. (5)). Independent audit 2026-11-11 (seven auditors; record
  05-knowledge/results/s91_f1_local_densities_20261110/audit/AUDIT_REPORT.md): PASS WITH CORRECTIONS, all applied;
  MISTAKE-761. Chan-Pak remains OPEN: (F) is a lower bound on a finite range.
session: opus-2026-11-10-S91
date: 2026-11-11
source: 05-knowledge/results/s91_f1_local_densities_20261110/README.md
scripts:
  - 04-computation/experiments/s91_f1_local_densities_20261110/rooted_spectrum/source_supply.py (+ .out, slope_lists_le73.txt): (A)-(D); section (6) of the output has the strong classes per offset and the correlations of (D)
  - 04-computation/experiments/s91_f1_local_densities_20261110/rooted_spectrum/source_supply_levels78.py (+ .out): (B'), strong tournaments of orders 10, 11 with at most 257 paths
  - 04-computation/experiments/s91_f1_local_densities_20261110/rooted_spectrum/head_parity_check.py (+ .out): (E), heads of 11, 13, 47, 67
  - 04-computation/experiments/s91_f1_local_densities_20261110/rooted_spectrum/line_density_order10.py, line_density_order9.out, line_density_order10.out; 04-computation/experiments/s91_f1_local_densities_20261110/PREREGISTRATION.md (second registration; the file is in the experiment's root directory): (F)
  - 04-computation/experiments/s91_f1_local_densities_20261110/rooted_spectrum/source_witnesses.py (+ .out), source_atlas.py (+ .out): literal recounts of 9-grid sources at levels 2-9, numbers of rooted tournaments of each order with each 9-grid pair
  - 04-computation/experiments/s91_f1_local_densities_20261110/rooted_spectrum/census.py, core.py, reps_le8.npz: the rooted census of orders <= 9 (rooted_le9.npz and codes9.npy are regenerated by census.py)
  - 04-computation/experiments/s91_f1_local_densities_20261110/rooted_spectrum/census10.py (+ .out): the aggregated census of order 10 of 2026-11-10, which predates the second registration (it writes agg10.npz, not kept in the repository; see (F) 3)
related:
  - 05-knowledge/hypotheses/HYP-9453-* (near-join family, exceptional set, windows), HYP-9454-* (lines 3H + 2 sigma)
  - 01-canon/theorems/THM-4865-* (two-step law and census), THM-5126-* (A) (tail law, offsets), THM-5130-* (F) (U >= H, product over components)
  - 05-knowledge/results/kps_s22_sign_flip_towers_20261110/SOURCES_AND_DOMINATION.md (concurrent: grids, excess types, abundancy index)
  - 05-knowledge/results/session261_creative_chanpak_20261111/ (domination-defect inequality, 9-grid levels 1-6), session262_interface_supply_20261111/
  - 01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md, THM-002 (odd-cycle formula)
  - 01-canon/theorems/THM-4744-* (3) (strong inventories below 625), THM-4884-* (2) (strong multiplicities m(h)), THM-4873-* (iii) (strong floor used at every order), THM-5109-* (C) (dyadic cones)
  - 05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md (Session 24, finite strong cores, inventory complete for H <= 73)
  - 05-knowledge/results/chanpak_55_105_golden_20261110/targets/PROOF.md, 05-knowledge/results/decomposition16_odd_spectrum_20261110/JOINT_CORES.md (all-order classifications from finite strong inventories)
  - 05-knowledge/results/odd_c_seed_symmetry_20261110/EXTENSIONS.md (qualifying permutations, k-vertex extensions), 05-knowledge/results/session216_overnight_constructions_20261107/star/PROOF.md (Theorem 3)
  - 05-knowledge/results/kps_s22_sign_flip_towers_20261110/README.md section 4 (qualifying count read as a cover count), 05-knowledge/results/session257_permutohedron_squares_20261111/history/sources/INCOMING_EXTENSION_CONTENT.md (refinement of the extension compiler)
  - 05-knowledge/hypotheses/HYP-9661-* (this session: line density)
  - 05-knowledge/results/s91_f1_local_densities_20261110/audit/AUDIT_REPORT.md (independent audit; the auditors' outputs cited in the text are in audit/audit-5129 and audit/numbers-vs-outputs; further recomputations of (B)-(D) in audit/prior-art-wording/a3_supply_replication.out and of the strong floors in audit/citations/check_classical_facts.out)
---

# THM-5129 - one-vertex sources: finite supply, complete small lists, permanent exceptions, line counts

## Setting

`H(T)` is the number of Hamiltonian paths of a tournament `T` (`H` of the empty tournament is 1). For a tournament
`B` and a vertex `y`:

- `D = H(B)` (the *offset*); `E = E_B(y)` = number of Hamiltonian paths of `B` ending at `y`;
- `U = U_B(y)` (the *slope*) = number of ordered pairs `(P, R)` of vertex-disjoint directed paths with
  `V(P) + V(R) = V(B)`, `R` starting at `y`, `P` possibly empty. Equivalently `U = sum over directed paths R
  starting at y of H(B - V(R))`; `U_even`, `U_odd` are the parts with `|V(R)|` even, odd;
- `Pi_2(C)` = number of ordered pairs `(P, R)` of vertex-disjoint paths covering `C`, each possibly empty
  (`Pi_2(point) = 2`).

`X => Y` means every vertex of `X` beats every vertex of `Y`. Closing laws (KNOWN: THM-4865; HYP-9453; THM-5126 (A)):

    H(B + z) = D - E + U             (z -> y, every other vertex -> z; the head),
    H(B => TT_t => z) = D + 2^t U    (B => TT_t, TT_t => z, z -> y, every other vertex of B -> z; t >= 1;
                                      t = 1 is the two-step value h + 2B of THM-4865).

A pair `(D, U)` arising from some `(B, y)` is a *one-vertex source*; `U >= D` (cut each Hamiltonian path at `y`).
Appending a sink to `B` keeps `D` and doubles `U`, so `{D + 2^t U : t >= 1} = {D + 2U' : U' = U 2^j}`.
`f(n)` is the least `H` of a strong tournament of order `n` (`n >= 3`; `f(1) = 1`, and there is no strong
tournament of order 2).

## (A) Product law and finite supply

**(A1) (KNOWN: THM-5130 (F2)(b) for the one-vertex set `O = {y}`.)** Let `C_1 => C_2 => ... => C_m` be the strong
components of `B` and `y` in `C_k`. Then

    D = prod_i H(C_i),        U = prod_{i<k} H(C_i) * U_{C_k}(y) * prod_{i>k} Pi_2(C_i).

*Proof.* Arcs between components point forward, so a path meets the components in increasing order and meets each
in one (possibly empty) path. `R` starts in `C_k`, so `P` contains `C_1, ..., C_{k-1}` and meets each in a
Hamiltonian path. In `C_k` the two traces form a cover with second path starting at `y`, first possibly empty; in
`C_i`, `i > k`, an ordered pair of possibly empty paths covering `C_i`. Any such choice glues back. ∎

A trivial component contributes `1` before the root, `2` after it, and `U = 1` if it is the root. The formula for
`D`, the multiplicativity of `H` over strong components, is classical (recorded in Busch 2006, cited below; in the
repository THM-1370 (1)).

**(A2) (Finite supply.)** For every odd `D` put `n(D) = max{n : f(n) <= D}`. Then the set of slopes of offset `D`,
over bases of every order, is

    Slopes(D) = Slopes_0(D) * {2^j : j >= 0},
    Slopes_0(D) = { d_1 u w :  D = d_1 h d_3,  d_1 attained,  u in Root(h),  w in Suf(d_3) },

where `Root(1) = {1}`, `Root(h) = {U_C(y) : C strong of order >= 3, H(C) = h, y in C}`, and `Suf(d)` is the set of
products `prod Pi_2(C_i)` over ordered factorisations `d = prod H(C_i)` into strong tournaments of order `>= 3`
(`Suf(1) = {1}`). All strong tournaments involved have order `<= n(D)`, so `Slopes_0(D)` is finite.

*Proof.* By (A1), with the nontrivial components before the root collected in `d_1`, the root component in `u`, the
later nontrivial ones in `w` and the later trivial ones in `2^j`: every slope of offset `D` lies in
`Slopes_0(D) * 2^j`. Conversely every product `d_1 u w 2^j` is realised, by the ordered join of a tournament with
`d_1` paths, the root component (a single vertex if `h = 1`), the chosen strong tournaments and `j` single
vertices. A strong component `C` of `B` has `H(C) <= D` (it divides `D`), hence order `<= n(D)` once `f` is
nondecreasing and unbounded: a strong tournament of order `n >= 4` has a vertex whose deletion leaves it strong (it
has a cycle of length `n - 1`: Moser, in Harary and Moser, Amer. Math. Monthly 73 (1966); every vertex lies on a
cycle of every length `3..n`: Moon, Canad. Math. Bull. 9 (1966) 297-301, Theorem 1), and `H(T) >= H(T - v)`
(insert `v` into a Hamiltonian path of `T - v` before the first vertex it beats, or at the end; injective), so
`f(n) >= f(n-1)` (KNOWN: THM-1370 (2), (3)). Moreover `v` lies on a 3-cycle of the strong tournament `T` (Moon's
theorem), and by the odd-cycle formula each odd cycle through `v` adds at least 2,
`H(T) - H(T - v) = 2 sum over odd cycles Z through v of H(T - V(Z)) >= 2`, so `f(n) >= f(n-1) + 2` and `f` is
unbounded. ∎

(The true size is exponential: for `n >= 3`, `f(n) = 5^((n-1)/3)`, `3 * 5^((n-3)/3)`, `9 * 5^((n-5)/3)` for
`n = 1, 0, 2 (mod 3)`: A. H. Busch, Electron. J. Combin. 13 (2006) #N3, proving sharp the construction of
J. W. Moon, Publ. Math. Debrecen 19 (1972) 101-104, whose lower bound `6^((n-1)/4)` already gives the finiteness
used here; CITED, not used. The product of `H` over strong components is in the same sources.) For the lists below
only this is used: `f(3), ..., f(9) = 3, 5, 9, 15, 25, 45, 75` (exhaustive in source_supply.py, order 9 on a list
containing every isomorphism class; KNOWN values: THM-1370 (4) through order 8, THM-4744 through order 13), so by
monotonicity a strong component of a base with `D <= 73` has order `<= 8`. That the inventory of strong
tournaments with `H <= 73` is complete at every order is KNOWN (Session 24, PROOF.md, sections 1 and 5; for every
`H < 625`: THM-4744 (3)); what is computed on it here is `U` and `Pi_2`.

**(A3) (Counting.)** For `X >= 1`, `#{U in Slopes(D) : U <= X} <= |Delta(D)| (1 + log_2 X)`, `Delta(D)` the set of
odd parts of `Slopes_0(D)`. In a dyadic interval `[a, 2a)` there are at most `|Delta(D)|` slopes of offset `D`.

*Proof.* By (A2) a slope of offset `D` is `delta 2^i` with `delta` in `Delta(D)` and `i >= j(delta)`, the least
exponent of `delta`. For each `delta` at most `1 + log_2 X` exponents `i >= 0` give a value `<= X`, and at most one
gives a value in `[a, 2a)`. ∎

## (B) The complete lists for `D <= 73` (FINITE-EXACT; valid at every order by (A))

Strong classes of orders 3-8: `1, 1, 6, 35, 353, 6008`, with `1, 4, 20, 156, 1984, 38935` distinct pairs
`(H, U)` and `1, 1, 5, 16, 131, 1567` distinct pairs `(H, Pi_2)`. The attained odd numbers `<= 73` are all except
7 and 21. The lists of all 35 attained offsets `<= 73` (1,066 pairs of an offset and an odd part) are in
rooted_spectrum/slope_lists_le73.txt; the audit's own lists (audit-5129/s3_slope_lists.txt) are identical. At the
six Fermat-type offsets, the odd parts `d` of the slopes with their least exponent `j` (so `d 2^i` is a slope iff
`i >= j`):

| offset | slopes `d : j` |
|---|---|
| 3 | `3:0` |
| 5 | `1:5, 3:1, 5:0, 7:0, 9:0` |
| 9 | `1:4, 5:1, 9:0, 11:0, 15:0, 21:0, 23:0, 25:0` |
| 17 | `5:2, 7:2, 9:1, 11:2, 17:0, 21:0, 31:0, 33:0, 55:0, 61:0, 69:0` |
| 33 | 34 odd parts (`1:6, 3:4, 5:3, 9:2, 13:2, 15:2, 17:1, 19:1, 21:1, 23:1, 31:2, 33:0, ...`; full list in source_supply.out) |
| 65 | 66 odd parts (`3:5, 5:4, 7:6, 9:5, 11:3, 13:3, 15:5, 17:2, ...`; full list in source_supply.out) |

Consequences, in the grid coordinates `N(d, r, s) = 1 + 2^s + d 2^(r+s)` of Session261:

1. A one-vertex source of offset exactly `2^s + 1` supplies the row `(d, s)` for all `r >= 1` iff `d` is listed with
   `j <= s`, and for `r >= j - s + 1` if `j > s`. At level 1 the only row is `d = 3` (slopes `3 * 2^i`).
2. By `U >= D` a source `(2^s + 1, 9 * 2^j)` has `9 * 2^j >= 2^s + 1`, i.e. `j >= s - 3` (KNOWN: kind-pasteur,
   SOURCES_AND_DOMINATION.md, Proposition 4, `v >= s - floor(log2 d)`). The five 9-grid sources of Session261
   (levels 2-6, `j = 0, 0, 1, 2, 5`) attain this bound at levels 3, 4, 5 (and have `j = 0` at level 2), so there
   their exponents are least possible without any list. The only level where the complete lists say more is 6: the
   bound allows `j = 3`, but the least exponent is 5 - neither `(65, 72)` nor `(65, 144)` is a source at any order.
3. Dyadic slopes (pure cones `D + 2^m`; the cone mechanism is KNOWN, THM-5109 (C)) at these six offsets:
   `(5, 32)`, `(9, 16)`, `(33, 64)` and their doubles only; none at 3, 17, 65. Each has `j > s`, as `U >= D`
   forces for a dyadic slope (`2^j >= 2^s + 1` gives `j >= s + 1`): the rows `d = 1` have finite heads at levels
   2, 3, 5. At levels 3 and 5 that bound is attained and the head is `r = 1` only (`N = 25`, `97`); at level 2 the
   least exponent is 5 where the bound allows 3 (`(5, 8)` and `(5, 16)` are not sources at any order), and the
   head (`r = 1, 2, 3`: `13, 21, 37`) contains `21 = 5 + 16 = N(1, 2, 2)`. `7 = 3 + 4 = N(1, 1, 1)` is the first
   point of the row `d = 1` at level 1, where the offset 3 has no dyadic slope at all.

The census of order `<= 9` contains every odd part of the 35 lists with its least exponent, except one: the slope
`(55, 352 = 11 * 2^5)` at offset `55 = 5 * 11` needs a base of order 10 (`C => {y} => C'` with `C` strong of order
5, `H(C) = 11`, and `C'` the strong tournament of order 4, `H(C') = 5`, `Pi_2(C') = 32`; source_supply.out prints
the odd part and its exponent, the audit the least order: numbers-vs-outputs/verify_supply.out (2'),
audit-5129/s9_misc.out).

**(B') Offsets 129 and 257 (FINITE-EXACT; valid at every order; source_supply_levels78.py).** The strong
tournaments with 129 or 257 Hamiltonian paths were enumerated completely: those of order `<= 9` from the class
lists of the census (every class of orders `<= 9`), the strong tournaments of order 10 with `H <= 257` (436
isomorphism classes, 423 with `H <= 255`), the strong one-vertex extensions of the latter with `H <= 257` (31
classes of order 11, least count `225 = f(11)`; by the proof of (A2) a strong tournament of order 11 with
`H <= 257` has a strong vertex-deleted subtournament with `H <= 255`), and the check that no strong tournament
of order 12 has `H <= 257`. (In the labelled lists of source_supply_levels78.py: 8,642 strong among
144,154 of order 10 and 4,368 strong among 21,156 of order 11; these four numbers depend on the generation scheme,
the class counts are the audit's, audit-5129/s10_levels78.out.) Strong classes with `H = 129`: 7, 26, 7 at orders
7, 8, 9; with `H = 257`: 30, 30, 13 at orders 8, 9, 10; none at other orders. KNOWN: these tournaments are in the
inventory of THM-4744 (3) (every strong class with `H < 625`; `f(11) = 225`, `f(12) = 375`), and the counts by
order agree with its `P_129 = 7z^7 + 26z^8 + 7z^9` and `P_257 = 30z^8 + 30z^9 + 13z^10`
(04-computation/experiments/h_spectrum_multiplicities_20261021/strong_ears_375.out, whose `P_d` also sum to the
class counts 436, 423 and 31 above). They were enumerated again because `U` and `Pi_2` are needed. The audit
recomputed (B') with its own class enumeration, after seeing the committed output, so not blind
(audit-5129/s11_supply257.out, s12_compare_out.out): the two lists are identical. With `129 = 3 * 43` and the
strong tournaments with 43 paths (orders 6, 7):

| offset | odd parts of the slopes | least exponent of `d = 9` | dyadic slope | slopes in `[h, 3h/2]` |
|---|---|---|---|---|
| 129 | 177, all with `j <= 7` | 4 | none | 49 of 65 (`0.754`) |
| 257 | 359, all with `j <= 8` | 5 | none | 93 of 129 (`0.721`) |

So at levels 7 and 8 every grid with a one-vertex source of offset `2^s + 1` is supplied for all `r >= 1`, and the
numbers of odd parts at the offsets `3, 5, 9, 17, 33, 65, 129, 257` are `1, 5, 8, 11, 34, 66, 177, 359`. The
sources `(129, 144)` and `(257, 288)` printed below attain the bound `j >= s - 3` (least possible by `U >= D`
alone); the enumeration is not needed for that. What (B') gives is the complete lists of 177 and 359 odd parts, the
absence of a dyadic slope at both offsets, and the window shares.

## (C) The six permanent omissions of the two-step family (PROVED from (A), (B))

For odd `N` let `r(N)` be the number of odd `D <= N/3` with `(N - D)/2` in `Slopes(D)`. Then for odd `N` in
`[3, 225)`: `r(N) = 0` exactly for `N = 7, 11, 13, 21, 47, 67`. So these six numbers are not of the form
`D + 2^t U` (`t >= 1`) for any rooted tournament of any order; with THM-4865 (e) they are the only odd numbers in
`[3, 382687]` with this property (in `[3, 6171561]` with its census of bases of order `<= 10`). `7` and `21` are
not path counts at all (THM-1370). The other four are head values `D - E + U = D + 2 U_even` at a root that is not
a sink:

    11 = 5 - 1 + 7,     13 = 5 - 1 + 9,     47 = 17 - 1 + 31,     67 = 33 - 7 + 41

(bases of order 4, 4, 6, 6). At a sink `E = U = D`, so every attained number is trivially a head value; only
heads with `E < D` carry information. `r(N) >= 5` on `[129, 225)`, mean `8.5`.

*Proof.* `U >= D` gives `D <= N/3 < 75`, where the lists are complete. ∎

## (D) Window sets (FINITE-EXACT, complete counts; the correlations are NUMERICAL)

HYP-9453 (2) asks for bases `h` whose response set contains every integer of a window `[h, c h]`, `c > 3/2`. For
attained `h <= 73` the number of slopes in `[h, 3h/2]` is

    1:1/1 3:1/2 5:3/3 9:3/5 11:3/6 13:4/7 15:4/8 17:4/9 19:5/10 23:8/12 25:7/13 27:5/14 29:9/15 31:9/16 33:13/17
    35:5/18 37:12/19 39:7/20 41:14/21 43:9/22 45:12/23 47:5/24 49:15/25 51:14/26 53:15/27 55:10/28 57:11/29
    59:19/30 61:12/31 63:5/32 65:24/33 67:9/34 69:26/35 71:23/36 73:27/37

(slopes / integers). The whole window is present only for `h = 1` and `h = 5`; on `37 <= h <= 73` the covered share
has mean `0.50`, minimum `0.156` (`h = 63`), maximum `0.743` (`h = 69`). On `23 <= h <= 73` (26 attained offsets)
the four offsets with exactly two strong classes (`35, 39, 47, 63`: shares `0.278, 0.350, 0.208, 0.156`) are among
the five poorest; the third poorest, `67` (`9/34 = 0.265`), has six. The numbers `m(h)` of strong classes with
`H = h` are KNOWN (THM-4884 (2)). The covered share and `m(h)` have Pearson correlation `0.68` on `23 <= h <= 73`,
`0.66` on `15 <= h <= 73` and `0.41` on all attained `3 <= h <= 73` (source_supply.out, section (6); NUMERICAL: the
value depends on the lower cutoff, which is a choice; `9, 17, 19` also have two strong classes, with shares
`0.60, 0.44, 0.50`).

By (A3) a fixed `h` has at most `|Delta(h)|` slopes per octave, hence at most `|Delta(h)|` in `[h, 3h/2]`. This
bound is smaller than the number of integers of the window only at `h = 3, 13, 35, 39, 43, 47, 63, 67` among the
attained `h <= 73` (slope_lists_le73.txt); at every other offset it does not exclude a full window: at
`h = 33, 65, 129, 257` there are `34, 66, 177, 359` odd parts against `17, 33, 65, 129` integers. What is
established is the count: no window `[h, 3h/2]` is full at the offsets computed (attained `h <= 73` other than 1
and 5; 129 and 257 in (B')); the covered share rises from about `0.50` (`37 <= h <= 73`) to `0.75` and `0.72`
(`129`, `257`); whether windows fill up for large `h` is OPEN. (F) counts the bases that an induction on
`N = D + 2U` has available instead.

## (E) The head correction is an alternating sum (PROVED)

    E = U_odd - U_even,        H(B + z) = D + 2 U_even,        H(B => TT_t => z) = D + 2^t (U_even + U_odd).

*Proof.* In `B + z` the odd cycles through `z` are `z, R` with `R` a path of `B` starting at `y` on an even number
of vertices (its last vertex beats `z` because it is not `y`). By the odd-cycle formula (THM-002),
`H(B + z) = H(B) + 2 sum over such R of H(B - V(R)) = D + 2 U_even`. Compare with `H(B + z) = D - E + U`. ∎

So a bare apex closes only the qualifying paths with an even number of vertices, and one sink of a tail removes
the parity condition (a path can pick up any subset of the sinks; half of the subsets have each parity). The head
is not the case `t = 0` of the tail formula: it differs by `E`. Checked on every root of every class of order
`<= 8` (58,653 rooted tournaments), with a separate count of `H(B + z)` through order 6.

Wording. Here a *qualifying path* is a directed path `R` of `B` starting at `y`, counted with the weight
`H(B - V(R))`, so that the total weight is `U`. This is the session's restriction to the one-vertex interface with
a transitive tail, an interpretation and not the general object: *qualifying* was defined first for the
permutations of the `n + k` vertices of a `k`-vertex extension with arbitrary attachments
(05-knowledge/results/odd_c_seed_symmetry_20261110/EXTENSIONS.md: an exact compiler of degree at most `2k`, whose
mean is THM-5126 (C); refined in Session257, INCOMING_EXTENSION_CONTENT.md), and kind-pasteur-2026-11-10-S22
(kps_s22_sign_flip_towers_20261110/README.md, section 4) first read the qualifying count of a one-vertex
extension as a cover count.

## (F) Line counts through rooted order 10 (FINITE-EXACT lower bounds; PRE-REGISTERED)

Let `r_m(N) <= r(N)` count the sources on the line `D + 2U = N` coming from rooted tournaments of order `<= m`
(with appended sinks). A line has about `N/6` lattice points (`D` odd, `D <= N/3`). The *share* of a block is the
mean of `r_m(N)/(N/6)` over its odd `N` (the registered definition; with the exact number of lattice points
`L(N) = floor((floor(N/3) + 1)/2)` in place of `N/6`, the mean of `r_m(N)/L(N)` differs from the entries below by
at most `0.002`, the largest difference being `0.294` against `0.296` in the first block:
audit-5129/s8_lines10.out).

| block | min `r_10` (at `N`) | median | share of the lattice points (`r_9`) | share (`r_10`) |
|---|---|---|---|---|
| `[2^7, 2^8)` | 5 (131, first of six) | 9 | 0.296 | 0.296 |
| `[2^8, 2^9)` | 7 (259) | 22 | 0.333 | 0.334 |
| `[2^9, 2^10)` | 24 (521) | 49 | 0.364 | 0.381 |
| `[2^10, 2^11)` | 63 (1027, 1073) | 120 | 0.414 | 0.464 |
| `[2^11, 2^12)` | 146 (2101) | 263 | 0.417 | 0.518 |
| `[2^12, 2^13)` | 352 (4109) | 606 | 0.373 | 0.588 |
| `[2^13, 2^14)` | 794 (8203) | 1259 | 0.165 | 0.616 |
| `[2^14, 2^15)` | 1660 (16447) | 2480 | 0.039 | 0.602 |

1. Every odd `N` in `[69, 32767]` has `r(N) >= 1` (KNOWN: to 382,687 from bases of order `<= 9` and to 6,171,561
   from order `<= 10`, THM-4865 (e)), and `r(N) >= N/64` for every odd `N` in `[129, 16384)`.
2. For `N` in `[4096, 8192)` the share of lattice points that are sources, by decile of `3D/N`, is
   `0.00 0.04 0.10 0.25 0.42 0.59 0.74 0.88 0.92 0.94` for odd `U` and
   `0.11 0.38 0.62 0.75 0.87 0.95 0.74 0.83 0.84 0.84` for even `U`: near the diagonal `U = D` (the root is close
   to a sink; `N = 3D + 2 (U - D)` is the form `3h + 2 sigma` of HYP-9454) `0.94` (odd `U`) and `0.84` (even `U`)
   of the lattice points are sources (`0.89` with both parities pooled: audit, numbers-vs-outputs/verify_lines.out).
3. Pre-registered before the order-10 run
   (04-computation/experiments/s91_f1_local_densities_20261110/PREREGISTRATION.md, second registration): L1
   `r_10(N) >= N/64` on `[129, 2^14)` (in sample at order 9: 758 violations) - CONFIRMED, no violation; L2 shares
   `>= 0.45`, `>= 0.43` on `[2^12, 2^13)`, `[2^11, 2^12)` - CONFIRMED (`0.588`, `0.518`); L3 shares increasing over
   the blocks `[2^6, 2^7)` to `[2^11, 2^12)` - CONFIRMED (`0.271, 0.296, 0.334, 0.381, 0.464, 0.518`); L4
   top-decile share `>= 0.45`, `>= 0.40` - CONFIRMED (`0.94`, `0.84`). L1, L2 and L4 are one-sided (`r_m` is
   monotone in `m`), L3 already held in sample at order 9, and all four were passed with wide margins. Aggregated
   data of the order-10 census existed before the registration (rooted_spectrum/census10.py and census10.out,
   2026-11-10; the file agg10.npz that it writes is not kept in the repository: value sets, of `U` among others,
   and pendant pairs `(D, D - E + U)`); the pairs `(D, U)` of order 10, from which the registered statistics are
   computed, were not among them. The first two rows of the table are exact values of `r`: `r_10(N) = r(N)` for
   every odd `N <= 543` (the audit's complete lists for every `D <= 257` give the exact `r(N)` for odd
   `N <= 773`: audit-5129/s11_supply257.out, s16_saturation.out). From `[2^9, 2^10)` on the rows are lower bounds
   only (`r_10(545) = 34 < r(545) = 35`: the source `(165, 190)` needs a base of order 11), and the last two
   blocks are far from saturated.

These are counts on a finite range. They do not prove that every line carries a source: that `r(N)` grows
linearly for all `N` is HYP-9661 (A) (in klein's form for the sole-exit sweep, HYP-9454 (1)), OPEN.

## 9-grid sources at levels 7, 8, 9 (FINITE-EXACT; concurrent with kind-pasteur-2026-11-10-S22, which is first)

Out-masks (bit `j` of row `i` set iff `i -> j`), root `y`, recounted literally (`H(B)` by subset DP and by all
permutations, `U` from the definition, `H` of the closed tournaments for `t = 1, 2, 3`):

    s = 7:  B = (444,273,258,6,268,286,319,382,8),  y = 8:  (D, U) = (129, 144 = 9 * 2^4)   [order 9]
    s = 8:  B = (82,36,129,135,142,157,190,3),      y = 7:  (D, U) = (257, 288 = 9 * 2^5)   [order 8]
    s = 9:  B = (402,360,3,261,270,285,317,382,4),  y = 8:  (D, U) = (513, 576 = 9 * 2^6)   [order 9]

Each is the only rooted tournament of its order with its pair (source_witnesses.out prints the numbers by order;
audit: numbers-vs-outputs/verify_supply.out (7), one rooted tournament up to isomorphism in each case), so the
sources at `s = 8, 9` are, up to isomorphism, those of kind-pasteur's table. At `s = 7` his table has the source
`(129, 288)`, of order 8 (also the only rooted tournament of its order with its pair); `(129, 144)` first occurs
at order 9. Each of the three sources above has the least exponent allowed by `U >= D` (`j = s - 3`). No source
`(65, 72)` exists (by (B), at any order), so the family `(2^s + 1, 9 * 2^(s-3))`, present at
`s = 3, 4, 5, 7, 8, 9` and `11` (`s = 11`: `(2049, 2304)`, order 9, kind-pasteur's table), absent at `s = 6` at
every order and at `s = 10` among the orders `<= 9` (kind-pasteur; audit-5129/s9_misc.out), is not uniform. A
uniform family for the 9-grid remains OPEN (Session262 OPEN_QUESTIONS 6).
