# Independent adversarial audit A of the (p, q, 1^k) defect family packet (mac-mini-2026-11-04-oddprimes)

**Auditor:** independent Claude session (audit A), 2026-11-04, worktree `math-wt-oddprimes-20261104`.
**Audited files:** THM-4871, HYP-9463, `05-knowledge/results/cpq_defect_family_mirror_folds_20261104.md` (sections 2,
2.3, 2.4, 5 in detail), `04-computation/experiments/cpq_defect_20261104/{cpq_defect_census.py, witness_examples.py,
cpq_*.out}`. Context read: THM-4721 (Theorem R, sign convention), THM-4728 (defect identity, C15 census to n = 13),
HYP-9344 (2), HYP-9340.
**Scripts and outputs of this audit:** `04-computation/experiments/cpq_defect_20261104/audit_A/` (listed in section 6).
No existing file was modified. Untracked files that appeared in the worktree during the audit from the still-running
session (`cpq_3_5_k7.out`, `cpq_defect_k7.py`, `fo_rule_modelcheck3.*`) are outside the audited set and were not
examined beyond noting that the `k = 7` run was in progress.

## 0. Verdicts

| file | verdict |
|---|---|
| THM-4871 | **PASS WITH CORRECTIONS** — statements 1, 2, 4 and the proofs are correct; statement 3(c)'s constant is wrong (`E_T(p) = 4 t_p^2 t_q`, not `4 t_p t_q`; the proof drops the `S_P` factor; invisible in the census because `t_3 = t_5 = 1`); two unsupported brute-force/row-count claims in the status; one wrong citation. |
| HYP-9463 | **PASS WITH CORRECTIONS** — statement 1's description of the first-defect row (unique witness, zero excesses) is false in five of the eleven families; statement 2's "equivalent form" holds only for `(3,5)`; the "indexing" reading of statement 2 is contradicted by the `(5,7)` rows. |
| note `cpq_defect_family_mirror_folds_20261104.md` | **PASS WITH CORRECTIONS** — same `4 t_p t_q` slip (2.3); the brute-force agreement claim for `(3,25) k <= 1` is not supported by the committed output (interrupted run, `k = 0` only); the characterisation of HYP-9344 (2) in sections 0 and 2.2 is inaccurate; section 5 is correct in every checked numeric and exact fact. |
| scripts and `.out` files | **PASS** — `cpq_defect_census.py` is correct (orbit method and `--brute` mode); every row re-run here or recomputed by two independent implementations agrees; `cpq_3_25_k2.out` is an interrupted run (no `k = 2` row, no summary table, brute only at `k = 0`). |

**Decisive independent confirmations.** (i) The headline refutation of HYP-9344 (2) is confirmed by a fresh
orbit-method implementation (own coloured-type enumeration validated against Burnside, own object builder, nauty
certificates cross-checked by networkx VF2 on all 260 tournament and 288 graph coincidences and 600 sampled
non-coincidences): C15 `(3,5,1^5)`: `tau_T = 17564`, `tau_E = 17480`, defect `84`, `W = 60`, `E_T = 260`,
`E_even = 282`, `E_wit = 2`; and `(3,5,1^6)`: `373804 / 373132`, defect `672`. (ii) A definition-level brute force
over all `2^c` `H`-fixed objects (N-action and N-evenness computed directly from the normaliser, no product
criterion, no orbit parametrisation) reproduces the `(3,5)` rows `k <= 4` (`2/2, 8/8, 32/31, 172/168, 1376/1360`
with `W = 0, 0, 1, 4, 14`, `E_T = 0, 0, 0, 4, 32`, `E_even = 0, 0, 0, 4, 34`), `(5,7) k <= 3`, `(3,7) k = 3`,
`(9,5) k <= 1`. (iii) The disjoint-orbit lemma survives an exhaustive search (17,752 overlapping prime-cycle pairs,
all lengths `<= 7` exhaustive, `11` sampled: every generated group has even order) and the theorem's sanity number
(838 prime cycles in `Aut T`, `n <= 8`, no overlap) is reproduced exactly. (iv) The one substantive error:
`E_T(p) = 4 t_p^2 t_q` (for `(7,11)` at `k = 7` the circulant constant-colour subfamily has 128 N-orbits forming 64
types of multiplicity 2, so `E_T(7) = 64`, not `32`; for `(7,13)`: 96, not 48).

## 1. THM-4871 — itemised

### 1.1 Statement 1 (normaliser-orbit parametrisation): PASS

* `N = Hol(Z_p) x Hol(Z_q) x Sym(F)`: the closure of the stated generators has order `p phi(p) q phi(q) k!` for
  `(3,5,k<=3)`, `(5,7,2)`, `(9,5,0..1)`, `(3,7,3)` (`audit_a_counts.out`), which equals `|C_{S_n}(h)| phi(pq)`
  `= pq k! phi(pq)`, the order of the full normaliser (all generators of `H` are `S_n`-conjugate to `h`); every closure
  element normalises `<h>` (asserted). Also for prime powers (`|N| = 1080` for `(9,5)`, `3000` for `(3,25)`).
* `H`-orbits on pairs and `c = (p+q)/2 + 2k + C(k,2)`: re-derived and recomputed (`c = 4, 6, 9, 13, 18` for `(3,5)`;
  `7, 9` for `(9,5)`; `14` for `(3,25)`), all as stated. The orbit list is right for prime powers too (a translate of
  `{i, i+d}` equals it only for `j = 0` since `2d != 0` mod odd `m`).
* `R = 2 t_p t_q b_k(4)`: the action of `N` on the datum `(S_P, S_Q, cross, kappa)` is the product action, so orbits are
  products. `b_k(4)` by Burnside over odd cycle types (`4^cycles 2^(pair orbits)/z_lambda`) gives
  `1, 4, 16, 88, 704, 8912, 188288, 6847200, 437069312` (k = 0..8), and explicit orbit enumeration (min over `S_k`)
  gives `1, 4, 16, 88, 704` (k <= 4): **matches**. `t_m` by brute-force dilation classes: `1, 1, 2, 4, 4, 6, 208` for
  `m = 3, 5, 7, 9, 11, 13, 25`; prime formula `t_p = (1/(p-1)) sum_{d | p-1 odd} phi(d) 2^((p-1)/(2d))` agrees
  (`t_7 = 2, t_11 = 4, t_13 = 6` confirmed, item (c)); `e_m = t_m` in every case.
* Product criterion for `N`-evenness: the order argument (an automorphism supported on one block reverses only edges
  inside that block under `P < Q < F`) is correct as written; order-independence of the sign (adjacent-transposition
  argument of THM-4721) is correct for automorphisms. Independent confirmation: my brute force computes `N`-evenness
  directly from the definition (every `g in N` fixing the graph word, reversal parity by orbit masks) and finds the
  number of `N`-even `N`-orbits equal to `R` in every row, e.g. `1408` at `(3,5,1^4)` (`audit_e1_brute.out`).
* `ec_k(4) = b_k(4)`: Theorem R with `P = 1` and the Young subgroup is a correct application; numerically
  colour-sign-even 4-coloured graph types `= 1, 4, 16, 88, 704, 8912, 188288 = b_k(4)` for `k <= 6`
  (`audit_e2_orbit.out`).

### 1.2 Statement 2 (disjoint-orbit lemma): PASS

Every step was checked. Transitivity of `<g, g'>` on `Delta` (the two orbits meet). The block argument: the
dichotomy per cycle (either the whole support lies in one invariant block, or every block meets the support in at
most one point) is exhaustive, and the three combinations are all contradictory, including `O' inside O` (then
`Delta = O` and case (i) gives `b = m` directly; prime degree is primitive anyway). Nested case: Burnside's theorem
(odd order excludes 2-transitivity because `r(r-1)` would divide `|K|`), `AGL(1,r)` has no non-identity element with
two fixed points, and `r - s >= 2` holds because `r > s` are odd primes (so `r - s` is even and positive): justified.
`t = m - r >= 3`: Jordan's theorem in exactly Wielandt's form (primitive, degree `m = r + t`, `t >= 3`, a prime
`r`-cycle, so `A_m <= K`), hypothesis `r <= m - 3` satisfied. `t = 1`: `m` even divides `|K|`. `t = 2`: Feit-Thompson
gives solvability, Galois gives the regular elementary abelian normal subgroup `V = F_l^e` with `K <= AGL(e, l)`, the
stabiliser of a point is `K cap GL(e, l)`, conjugating by a translation puts a fixed point of `g` at `0`, and the
fixed-point set of a linear map is a subspace of size `l^j != 2` for odd `l`: correct. The lemma does use `r, s`
odd, which is automatic from `|G|` odd. The alternative via Jones 2014 is also valid (all primitive groups containing
a cycle other than `C_l <= K <= AGL(1,l)` are `A_n`, `S_n`, `PGL(d,q) <= K <= PGamma L(d,q)`, `PSL(2,11)`,
`M_11, M_12, M_23, M_24`, all of even order).

Adversarial search (`audit_b_lemma_search.py`, Part 1): for prime lengths `r <= s` in `{3, 5, 7}` and every overlap,
EVERY `s`-cycle on the second support (exhaustive, up to `Sym`-conjugacy this covers all overlapping pairs) and 300
random `s`-cycles for `s = 11`: 17,752 pairs, every group `<g, g'>` of even order; no counterexample. Part 2 reproduces
the theorem's sanity check exactly: over all tournaments of orders 3-8 (2, 4, 12, 56, 456, 6880 types), the elements
of `Aut T` that are single prime-length cycles number `2, 4, 12, 32, 124, 664`, total **838**, with no two overlapping
unequal supports.

### 1.3 Statement 3 (tournament multiplicities): PASS WITH ONE CORRECTION

(a) and (b) are correct: `Aut T` has odd order; the lemma gives `P' = P` or `P' <= F`, `Q' = Q` or `Q' <= F`; the
`p`-cycles with support `P` generate (on `P`) an odd-order transitive group of prime degree, inside `AGL(1,p)` by
Burnside, with the unique Sylow `p`-subgroup `<c_P|_P>`; CRT gives `H' = H`. (c): `m_T <= 2`, `m_T = 2` iff `kappa` is
circulant with constant colour, and the two classes are not `Aut T`-conjugate (an automorphism conjugating `H` to
`H'` would map `P` onto `F` against the constant direction of the `P`-`F` arcs): all correct, and the datum
transformation `(cross, (x,y), S_Q, kappa) -> (y, (1-x, cross), S_Q, T[P])` was verified to land in the same type for
every type of the subfamily at `(3,5), (3,7), (5,7), (7,11), (7,13)` (`audit_c_ET_count.out`).

**CORRECTION C1 (substantive).** Quote (statement 3(c)): "`E_T(p) = 4 t_p t_q` ... there are `8 t_p t_q` such
`N`-orbits"; proof: "The `N`-orbits with circulant constant-colour `kappa` number `2 (cross) x t_q (S_Q) x 4 (colour)
x t_p (circulant tournaments on `p` points up to isomorphism ...), i.e. `8 t_p t_q`"; title: "equals `4 t_p t_q` at
`k = p`". The count omits the factor `t_p` for the dilation class of `S_P`. Corrected statement: **the `N`-orbits with
circulant constant-colour `kappa` number `t_p (S_P) x t_q (S_Q) x 2 (cross) x 4 (colour) x t_p (kappa) = 8 t_p^2 t_q`,
so `E_T(p) = 4 t_p^2 t_q`.** The census cannot distinguish the two formulas because `t_3 = t_5 = 1` (the `k = p` rows
exist only for `p = 3, 5`). Independent computation of the subfamily (`audit_c_ET_count.out`, nauty certificates
cross-checked by VF2): `(7,11)` at `k = 7`, `n = 25`: 128 N-orbits, 64 types, all of multiplicity 2, `E_T(7) = 64`
(theorem: 32); `(7,13)` at `k = 7`: 192 orbits, 96 types, `E_T(7) = 96` (theorem: 48); `(3,5)`, `(3,7)`, `(5,7)`:
`4, 8, 8` as in the census. The proof's structure (`m_T <= 2`, characterisation of `m_T = 2`) is right, so the
corrected statement is PROVED by the same argument. Downstream: the pre-registered tests of HYP-9463 (2) at
`(3,7) k = 7` and `(5,7) k = 7` use `E_T(3) = 8`, `E_T(5) = 8` and are unaffected; any extension to `p >= 7` must use
`4 t_p^2 t_q`. The same slip is in the note, section 2.3 ("`E_T(p) = 4 t_p t_q` exactly").

### 1.4 Statement 4 (witnesses and the first defect): PASS

(a) `K_{q,q}` with parts `Q` and `P cup F`: verified directly (`audit_d_witness.out`). `(3,5)` at `n = 10`: `H`-fixed,
`N <= Aut X` (`|N| = 240`), sign `+1` on all of `N`; `|Aut X| = 28800 = 2 (5!)^2` (nauty, and by closure), exactly
half of `Aut X` has odd sign, the side swap has sign `-1`; `N cap Aut X = N` all even. `(5,7)` at `n = 14`: `|N| = 1680`
all even on `X`, `|Aut X| = 50803200 = 2 (7!)^2`, one odd nauty generator, side swap odd. The orbit-form sign argument
(a `<g>`-orbit of edges contributes the parity of "some power of `g` swaps the pair") was re-derived and is correct.
The logic of the inequality is right and correctly conditioned: `E_T(q - p) = 0` needs `q - p < p`, i.e. `q < 2p`, and
`E_even, E_wit >= 0` by definition; Bertrand's postulate makes the family infinite. The census rows `(7,11) k = 4` and
`(11,13) k = 2` (in the `.out` files, not in the theorem's table) are the next two members and show exactly `W = 1`,
defect `1` (recomputed here, `audit_e2_orbit_extra.out`).
(b) `(3,7) k = 3`: the four configurations are `H`-fixed, `N`-even (`|N cap Aut X| = 1512, 1512, 504, 504`, all even),
sign-odd (the swap `P <-> F` reverses `9` or `51` edges in the order `P < Q < F`, odd either way), pairwise
non-isomorphic (distinct certificates and VF2); my orbit census at `(3,7,1^3)` finds exactly these four witnesses
(`C_Q in {empty, {2,3,4,5}}` x `cross`). Minor wording: "reversing `9` edges (the `3` swapped edges `{i, sigma i}` and
unflipped pairs)" conflates the direct count (`9` or `9 + 6q`) with the orbit form (`3` self-paired orbits); both are
odd, so the conclusion stands.

### 1.5 Statements 5-6 (census): PASS, with two status corrections

* Re-runs of `cpq_defect_census.py` for `(3,5) k <= 4 --brute 3`, `(5,7) k <= 3`, `(9,5) k <= 2 --brute 1`
  (`audit_A/rerun_*.out`) reproduce the committed `.out` rows and THM-4728's table (`R = 2, 8, 32, 176, 1408`;
  `tau_T = 2, 8, 32, 172, 1376`; `tau_E = 2, 8, 31, 168, 1360`; `W = 0, 0, 1, 4, 14`; `E_T = 0, 0, 0, 4, 32`;
  `E_even = 0, 0, 0, 4, 34`), with brute-force AGREE at `(3,5) k <= 3` and `(9,5) k <= 1`.
* Independent brute force from the definitions (`audit_e1_brute.py`, no shared code; the `N`-action on orbit words,
  `N`-orbit representatives, `N`-evenness, types, signs, `m_T`, `m_X`, `m^ev_X` all computed directly):
  `(3,5) k = 0..4`, `(9,5) k = 0, 1`, `(5,7) k = 2, 3`, `(3,7) k = 3` — every number (`R`, `tau_T`, `tau_E`, `W`,
  `E_T`, `E_even`, `E_wit`, Theorem R, the identity) agrees with the `.out` files. In particular `32/31/1` and
  `172/168/4` (with `W = 4`, `E_T = 4`, `E_even = 4`), and `1376/1360/16` at `k = 4`, which the committed outputs do
  not brute-force.
* Independent orbit method (`audit_e2_orbit.py`): `(3,5) k = 2..6`, `(3,7) k = 3, 4`, `(5,7) k = 5`, `(9,5) k = 0, 1, 4`,
  `(9,7) k = 3`, `(7,11) k = 4`, `(11,13) k = 2`, `(3,25) k = 0`: all rows agree, including the headline
  `17564 / 17480 / 84` and `373804 / 373132 / 672`, `(5,7) k = 5`: `35640 / 35552 / 88`, `(3,7) k = 4`:
  `2752 / 2751 / 1` with `W = 17, E_T = 64, E_even = 48`, `(9,5) k = 4`: `4224 / 4223 / 1` with `W = 1409`,
  `E_T = 1408`, `(3,25) k = 0`: `410 / 410` with `W = 4, E_T = 6, E_wit = 2`. VF2 cross-checks of nauty at `k = 5`:
  all 260 coincident tournament pairs and 288 coincident graph pairs are isomorphic, 600 sampled distinct-certificate
  pairs are non-isomorphic.
* Statement 6: `Z_9`: the 4 dilation classes of circulant 9-tournaments give 3 isomorphism types; `Z_25`: 208 classes
  give 205 types, hence `E_T(0) = 2 (208 - 205) = 6` for `(3,25)` (`audit_a_counts.out`); `E_T(k) = 2 t_q b_k(4)` and
  `W = E_T` below `|p - q|` hold in every `(9,q)` row of the `.out` files. FINITE-EXACT typing is appropriate.

**CORRECTION C2.** Quote (status): "a brute-force pass over all `2^c` `H`-fixed objects with the full normaliser
action agrees for `k <= 3` at `(3,5)` and `k <= 1` at `(9,5), (3,25)`" (also note, section 2: "it agrees with the
orbit method at `(3,5) k <= 3`, `(9,5) k <= 1`, `(3,25) k <= 1`"). The committed `cpq_3_5_k6.out` contains no brute
lines at all, and `cpq_3_25_k2.out` contains a brute line only for `k = 0`; the file stops after the `k = 1` orbit row
without the `k = 2` row or the summary table (interrupted run, `2^16` objects x `|N| = 3000`). Corrected: the `(3,5)
k <= 3` and `(9,5) k <= 1` agreements are reproduced by this audit's re-run and by its own brute force (which also
covers `(3,5) k = 4`); **`(3,25) k = 1` has not been brute-forced by anyone**, only `k = 0`.

**CORRECTION C3.** Quote: "`tau_T >= tau_E` in all 64 rows; THM-4721 Theorem R and THM-4728's identity hold in all 64
rows" (and HYP-9463 "64 census rows over eleven families", the note "64 census rows"). The table of statement 5 has
**53** rows over eleven families; the `.out` files have **62** distinct rows over thirteen families (adding `(7,11)
k <= 4`, `(11,13) k <= 2`, `(3,25) k = 1`). No count gives 64. Replace by the actual count and say which rows are meant.

### 1.6 Citations and typing

* **CORRECTION C4.** Quote: "G. A. Jones, J. Algebraic Combin. 40 (2014) 227-240". The paper "Primitive permutation
  groups containing a cycle" is **Bull. Aust. Math. Soc. 89 (2014), no. 1, 159-165** (online 18 July 2013, DOI
  10.1017/S000497271300049X; arXiv:1209.5169). Its abstract: of the primitive groups containing a cycle, only the
  alternating and symmetric groups contain a cycle fixing at least three points.
* Jordan (Wielandt 13.9) and Burnside (Wielandt 11.7) are used with the right hypotheses; Galois's theorem and
  Feit-Thompson as stated; Alspach-Parsons/Babai for `Z_p` CI and Palfy/Muzychuk for `Z_9`, `Z_25` not CI are
  appropriate (for tournaments the relevant fact is that `Z_9` is not a DCI-group, which the census exhibits
  directly: 4 classes, 3 types). The Burnside section number ("2nd ed., sect. 169") was not checked.
* Typing: statements 1-4 are PROVED after correction C1 (the corrected 3(c) is proved by the same argument); 5-6
  FINITE-EXACT; nothing numerical is typed PROVED. "HYP-9340 (a) fails for every prime pair `p < q < 2p` at degree
  `2q`" is correctly derived from 3(b) + 4(a) and correctly restricted (the `Scope` section says equality and the
  `(3,7)`, `(9,q)` rows are census only). No overclaim found in the theorem about `(3,7)`, `(9,5)`.
* Scope omission (minor): the sanity check "838 prime-length cycles" has no script or output in the experiment
  directory; it is reproduced by `audit_b_lemma_search.py` (Part 2), which can be cited.

## 2. HYP-9463 — itemised

**CORRECTION H1 (statement 1).** Quote: "the unique witness at `k = d` is `K_(q,q)` on the larger orbit against the
smaller orbit plus `d` fixed points (THM-4871 statement 4(a)), with `E_even(d) = E_wit(d) = 0` and
`E_T(d) = E_T^F(d) = 0`." This is true exactly in the rows with `p, q` prime and `d < p`: `(3,5) k=2`, `(5,7) k=2`,
`(7,11) k=4`, `(11,13) k=2` (one witness, the `K_{q,q}`, all excesses zero: `audit_h_witness_anatomy.out`). It is
**false** in the other five families with a computed first defect: `(3,7) k=4`: `W = 17` (one `K_{7,7}` and sixteen
`K_{3,3}`-type configurations), `E_even = 48`, `E_T = 64`; `(9,5) k=4`: `W = 1409`, `E_T = 1408`; `(9,7) k=2`:
`W = 65`, `E_T = 64`; `(9,11) k=2`: `129 / 128`; `(9,13) k=4`: `8449 / 8448`. The note's own section 2.4 states the
`(3,7)` anatomy (`17 + 48 - 64 = 1`) and contradicts this sentence. Corrected: "for primes `p < q < 2p` the census
shows a unique witness `K_{q,q}` with all excesses zero; in the families with `d >= p` or a prime-power orbit the
first defect `1` is a net value of many witnesses and coincidences."

**CORRECTION H2 (statement 2, "equivalent form").** Quote: "Equivalent form: `tau_E(k) = tau_T(k) - tau_T(k - d) +
b_(k-d)(4)` whenever all coincidences at level `k - d` come from the fixed points." Since
`tau_T(k-d) = 2 t_p t_q b_(k-d)(4) - E_T(k-d)`, this is equivalent to the defect formula only when `t_p t_q = 1`, i.e.
for `(3,5)`. Counter-check at `(5,7) k = 3`: `352 - 16 + 4 = 340 != tau_E(3) = 348`. Corrected equivalent form:
`tau_E(k) = tau_T(k) - b_(k-d)(4) + E_T^F(k-d)`, or `defect(k) = R(k-d)/(2 t_p t_q) - E_T^F(k-d)`.

**CORRECTION H3 (statement 2, reading).** Quote: "The reading is that the `K_(q,q)` witnesses at level `k` are
indexed by the `H`-symmetric tournaments at level `k - d` with the cross orientation fixed (... the remaining `k - d`
points carry a 4-coloured structure)". This would give `W(k) >= b_(k-d)(4)`; the census contradicts it: `(5,7)`
`W(3) = 2 < b_1(4) = 4`, `W(4) = 6 < 16`, `W(5) = 24 < 88` (the defect `4, 16, 88` there is `W + E_even`:
`2 + 2`, `6 + 10`, `24 + 72 - 8`). The formula `b_(k-d)(4) - E_T^F(k-d)` is supported only as a net identity.

**H4.** The row count "64" (statements 1, 3 and the title): see C3.

Statements 3-5 and the pre-registered tests are correctly typed (NUMERICAL / OPEN / typed). The `(7,11) k <= 4` test of
statement 1 has in fact been run (`cpq_7_11_k4.out`: `0, 0, 0, 0, 1`) and agrees with the pre-registration, though the
`k = 4` value `>= 1` is the proved case. Statement 4(a)'s reduction (a quantifier-free `phi` in the atoms `S(c_i x)`
has Boolean degree `<= #atoms`, excluded by the Codex bound `degree >= ceil((e-1)/a)`, Corollary 8 of
`circulant_boolean_degree_20261020/BOOLEAN_DEGREE_BOUNDARY.md`) is consistent with the cited packet (read here).
The "E_T^F" convention for `p = 9` is a definition; the equality `E_T = 2 t_q b_k(4)` in the `(9,q)` rows is
numerical, as typed.

## 3. The note — itemised

* **N1** = C1 (section 2.3: "`E_T(p) = 4 t_p t_q` exactly" -> `4 t_p^2 t_q`).
* **N2** = C2 (section 2: brute-force agreement at `(3,25) k <= 1` unsupported; section 7's reproduction command
  `python3 cpq_defect_census.py 3 25 2 --brute 1` did not finish in the committed output).
* **N3.** Quote (section 0): "HYP-9344 (2)'s `88` (it took the witness count `b_(k-2)(4)` for the defect and missed the
  tournament coincidences `E_T(k - 2)`...)"; section 11: "`b_(k-2)(4)` was the right scale and the wrong quantity".
  HYP-9344 (2) did not take a witness count for the defect: it tabulated `W = 1, 4, 14` separately and extrapolated the
  *defect* sequence `1, 4, 16 = b_0, b_1, b_2 = R(k-2)/2`. Corrected: "it extrapolated `defect(k) = b_(k-2)(4)` from
  `k = 2, 3, 4` and missed the correction `-E_T(k-2)`, which first appears at `k = 5`". The refutation itself is fairly
  described: `84` lies outside `{64, 88}` and so refutes both pre-registered alternatives, as HYP-9344 (2) itself
  stipulated; "what it got right is the scale" is a fair reading given that `b_(k-2)(4) - E_T(k-2)` is only NUMERICAL.
* **N4** (section 2, nauty keying). Quote: "keying by the canonical form alone silently merged non-isomorphic
  coloured tournaments (`82` instead of `b_3(4) = 88`)". The phenomenon is real (pynauty's certificate under a
  vertex colouring is not an invariant across partition shapes), but the number `82` was not reproduced: keying the
  512 coloured 3-tournaments by `certificate(graph with vertex_coloring)` alone gives 8 keys here; the figure depends
  on the exact construction. Harmless; state the construction or drop the number. The present keying (cell sizes +
  certificate) is an isomorphism invariant whose number of values equals the Burnside count, which proves it is a
  complete invariant on these sets (so the coloured-type enumeration is correct, not merely plausible).
* **N5** = C3 (section 2.1 and 3: "64 rows").
* **N6** (section 2.4): "at `(3,7) k = 4` the 17 witnesses are the one `K_(7,7)` and sixteen `K_(3,3)`-type
  configurations (2 cross x 2 circulant classes x 4 colours of the extra point), and `17 + 48 - 64 = 1`": confirmed
  (`audit_h_witness_anatomy.out`: 16 configurations with `C_P` empty and three fixed points of colour `(1, cross)`
  pairwise non-adjacent, 1 `K_{7,7}`; the extra point is joined to `P cup F'` or not, and to `Q` or not).
  Section 2.3's description of the lemma's proof ("Jordan (union size `>= r + 3`), parity (`r + 1`), or Galois +
  Feit-Thompson (`r + 2`)") is accurate.
* **N7** (section 5): see item (g) below; every numeric and exact fact checked is correct.

## 4. Item (g): section 5 of the note (`audit_g_numerics.out`)

* MRB `C = 0.18785964246206712024851793405427...` (mpmath, alternating-series acceleration, 40 digits): matches.
  `B = 2.09679112439152612572503067508...` with `r_1 = 1, r_2 = 4, r_n = exp(-W_0(-ln n/n))`: matches.
  `B - C = 1.90893148192945900547651274103...` matches, and equals `4 - sqrt 2 + sum_{n>=3} (-1)^n (r_n - n^(1/n))`
  to all 32 digits computed. `B = 3 - length(E)` and `C = sqrt 2 - 1 - int_E f'` are exact rearrangements (pairing of
  consecutive terms of a convergent alternating series; `int_[r_(2k+2), r_(2k+1)] f' = f(r_(2k+1)) - f(r_(2k+2)) =
  f(2k+1) - f(2k+2)` by the mirror property, checked on `[r_4, r_3]`); the paired partial sums coincide with the
  alternating partial sums (see the output). `(r_n - n^(1/n))/(ln n/n)^2 = 1.0175, 1.0023, 1.00029` at
  `n = 10^3, 10^4, 10^5`: matches (`z^2 + O(z^3)`). `m(1.01) = 658.805`, `m'(1.01) = -76722.4`, asymptotic
  `-m^2/(ln m - 1) = -79051`: matches. Cayley expansion at `n = 3`: 60 terms leave `0.0914`: matches.
* Ball: `V(n) = 1, 2, 3.1416, 4.1888, 4.9348, 5.2638, 5.1677, 4.7248`; `x* = 5.25695`, `V(x*) = 5.27777`; mirror
  dimensions `12.7641, 10.8136, 9.2046, 7.8204, 6.6028, 5.5174` for `n = 0..5`: all match to the digits given.
  `sum (-1)^n V(n) = e^pi erfc(sqrt pi) = 0.282059176175682645184119553197` (both sides computed to 30 digits; the
  generating function `sum V(n) t^n = e^(pi t^2)(1 + erf(sqrt pi t))` was re-derived from `e^(x^2) erf x =
  (2/sqrt pi) sum 2^j x^(2j+1)/(2j+1)!!`). "No two distinct integer dimensions have the same ball volume": the
  parity/`|n - m| >= 2` case is a nonzero power of `pi` times a rational (Lindemann), and for `|n - m| = 1` the ratio
  `sqrt pi Gamma((n+1)/2)/Gamma(n/2 + 1)` is strictly decreasing (`psi(x + 1/2) < psi(x + 1)`), equal to `16/15` at
  `n = 5` and `5 pi/16` at `n = 6` (checked numerically, strictly decreasing for `n <= 14`): the proof is correct.
  (Remark: for `n` even, `m = n - 1` the ratio already contains `pi^1`, so Lindemann covers it too; only `n` odd,
  `m = n - 1`, where the ratio is rational, needs the monotonicity.)
* `-7/4`, `-29/16` (sympy, gp): `(f^3(x) - x)/(f(x) - x) = (8x^3 + 4x^2 - 18x - 1)^2/64` at `c = -7/4`: exact. The
  cubic has roots `-1.74698, -0.054958, 1.30194`, polynomial discriminant `200704 = 2^12 7^2`, **field discriminant
  49**, Galois group `A_3`, and `nfisisom` with `x^3 + x^2 - 2x - 1` succeeds: it generates `Q(zeta_7)^+`. At
  `c = -29/16`: `(f^3 - x)/(f - x) = (4x-5)(4x+1)(4x+7)(64x^3 + 16x^2 - 164x + 23)/4096`, the rational cycle
  `-7/4 -> 5/4 -> -1/4` has multiplier `8 (-7/4)(5/4)(-1/4) = 35/8`, the companion cubic has multiplier
  `8 (-23/64) = -23/8` and field discriminant `961 = 31^2` (cyclic cubic): all exact facts confirmed.
  `c(c+1)^2 = -63/64`, `f^3(0) = -7/256` at `c = -7/4`: confirmed.
* Five-vertex tournaments (gentourng + nauty): 12 types, **8 self-converse + 4 not**, `|Aut|` distribution
  `{1: 7, 3: 4, 5: 1}`, **6 strong + 6 not**: all as stated. Self-converse counts `n = 1..7`:
  `1, 1, 2, 2, 8, 12, 88`; A002785 (OEIS, fetched 2026-11-04) reads `1, 1, 2, 2, 8, 12, 88, 176, 2752, ...`:
  consistent with the note's "match ... as far as computed here". OEIS (JSON API, 2026-11-04): no match for
  `1, 4, 16, 88, 704, 8912, 188288` or for `1, 4, 16, 84, 672`, as the note says.

## 5. Item (f): the reading of HYP-9344 (2)

HYP-9344 (2) pre-registered `defect(5) = 88` against `64` and stated that a value outside `{64, 88}` refutes both.
THM-4871 and the note report `84` and call HYP-9344 (2) refuted: fair. "What HYP-9344 (2) got right is the scale" is
fair (the first term of the NUMERICAL law HYP-9463 (2) is HYP-9344's `b_(k-2)(4)`). What is not fair is the
characterisation in sections 0 and 11 of the note (N3): HYP-9344 (2) never equated witnesses with the defect; its
witness list (`W = 1, 4, 14`) was separate and correct, and its error was extrapolating the defect sequence without the
`-E_T(k-2)` term. THM-4871 statement 5's wording ("pre-registered `defect(5) = 88` (against the alternative `64`)...
the value is `84`") is accurate.

## 6. Scripts and outputs of this audit (`04-computation/experiments/cpq_defect_20261104/audit_A/`)

| file | content |
|---|---|
| `audit_a_counts.py/.out` | `b_k(4)` by Burnside and by orbit enumeration; `t_m`, `e_m`, prime formula; `c`, `|N|` by closure; `Z_9`, `Z_25` isomorphism types of circulant tournaments; nauty keying check |
| `audit_b_lemma_search.py/.out` | exhaustive/sampled search for overlapping prime cycles generating odd-order groups (sympy); prime cycles in `Aut T`, `n <= 8` (838, no overlap) |
| `audit_c_ET_count.py/.out` | the circulant constant-colour subfamily at `k = p` for `(3,5), (3,7), (5,7), (7,11), (7,13)`: N-orbits, types, multiplicities, VF2 cross-check, datum transformation; `E_T(p) = 4 t_p^2 t_q` |
| `audit_d_witness.py/.out` | `K_{q,q}` at `(3,5) n=10`, `(5,7) n=14`: H-fixed, `N`-signs, `Aut X` order and signs, side swap; the four `(3,7) k=3` configurations |
| `audit_e1_brute.py/.out` | definition-level brute force over all `2^c` objects (own N-action/N-evenness), `(3,5) k<=4`, `(9,5) k<=1`, `(5,7) k<=3`, `(3,7) k=3` |
| `audit_e2_orbit.py/.out`, `audit_e2_orbit_extra.out` | own orbit-method census: `(3,5) k=2..6` (VF2 at `k<=3, 5`), `(3,7) k=3,4`, `(5,7) k=5`, `(9,5) k=0,1,4`, `(9,7) k=3`, `(7,11) k=4`, `(11,13) k=2`, `(3,25) k=0` |
| `audit_h_witness_anatomy.py/.out` | witness anatomy at the first-defect rows; `(3,7) k=4`: 1 `K_{7,7}` + 16 `K_{3,3}`-type |
| `audit_g_numerics.py/.out` | section 5 numerics (mpmath), factorisations (sympy), cubic fields (gp), 5-tournaments and self-converse counts (gentourng + nauty) |
| `rerun_cpq_3_5_k4_brute3.out`, `rerun_cpq_5_7_k3.out`, `rerun_cpq_9_5_k2_brute1.out`, `rerun_witness_examples.out` | re-runs of the author's scripts |

## 7. Summary of required corrections (for the author)

1. THM-4871 title, statement 3(c), proof of 3; note 2.3: `E_T(p) = 4 t_p^2 t_q` and `8 t_p^2 t_q` N-orbits (C1).
2. THM-4871 status, note section 2: brute-force agreement claimed for `(3,25) k <= 1` is supported only at `k = 0`;
   the `(3,5) k <= 3` brute lines are not in the committed output (now reproduced here) (C2).
3. "64 rows" -> 53 (table) or 62 (`.out` files, 13 families) in THM-4871, HYP-9463, the note (C3).
4. Jones 2014: Bull. Aust. Math. Soc. 89 (2014) 159-165 (C4).
5. HYP-9463 (1): the unique-witness/zero-excess description holds only for the prime rows with `d < p` (H1).
6. HYP-9463 (2): the "equivalent form" is valid only for `t_p t_q = 1`; the "indexing" reading is contradicted by
   `(5,7)` (H2, H3).
7. Note sections 0 and 11: HYP-9344 (2) extrapolated the defect, not a witness count (N3); the nauty "82" is
   unreproduced (N4); cite a script for the 838-cycle sanity check (Scope).
8. Optional wording: 4(b) "reversing 9 edges" (1.4).
