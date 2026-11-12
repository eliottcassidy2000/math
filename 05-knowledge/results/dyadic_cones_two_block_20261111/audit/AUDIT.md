# Independent audit of `dyadic_cones_two_block_20261111` (THM-5145, HYP-9675)

Auditor: independent audit agent, 2026-11-11. Scope: README.md, PREREGISTRATION.md, programs, `out/`,
THM-5145, HYP-9675, and the inherited statements (THM-5126, THM-4875 (D), THM-5109, THM-5085, HYP-9658).
All audit files are in this directory. The packet's code was not reused. The audit wrote its own census, its own
literal recount and its own one-back-arc check.

## Overall verdict: PASS WITH CORRECTIONS

The three theorems (sections 3, 1 to 3) are **correct**. The census they rest on was reproduced by an independent
program. The audit's order-9 and order-10 type tables match the packet's files entry for entry. Every witness the
audit sampled, and every type that had no witness on file, was recounted literally. The audit found no
mathematical error in Theorems 1 to 3, the window lemma, the one-back-arc formulas, or the template examples.

The majors are about credit, one unsupported equivalence and one overstatement of what was verified:

- **M1:** the order-8 two-parameter template census and its "no independent pair" exclusion are prior work.
- **M2:** cones and rays with negative offsets, and the columns `r = +-1`, are prior work.
- **M3:** the equivalence "quantitative type conjecture <=> Chan-Pak" is not supported, and it is false at `N = 1`.
- **M4:** 293 of the Theorem 2 types had no witness on file.
- **M5:** parts of section 4 are classical or known.

| item | verdict |
|---|---|
| A. Theorem 1 + cone census | PASS WITH CORRECTIONS (minor wording; credit M2) |
| B. Theorems 2, 3 + type census | PASS WITH CORRECTIONS (M4) |
| C. Window lemma; "every path count is H(L o 0)"; density remark | PASS WITH CORRECTIONS (M3) |
| D. One back arc | PASS WITH CORRECTIONS (M5, m3) |
| E. Blow-up templates | PASS WITH CORRECTIONS (M1, m4-m6) |
| F. Section 6 routes | PASS WITH CORRECTIONS (M3) |
| G. Pre-registration | PASS (m7) |
| Prior art / "NEW" claims of THM-5145 | partly REFUTED (M1, M2, M5) |

## What was recomputed, and how

1. **Own census: `acensus.c`.** Offset and response are computed from end- and start-path tables of `B`, using
   `c0 = H(B) - sum_{e != f in O} Q(e,f)` and `M = sum_{y in O} U_B(y)` as defined in THM-5126 (1'). For each
   `(c, r)` with `c` odd `<= 255` and `|r| <= 2047`, the program records the least exponent, plus a witness when
   `c <= 31` and `|r| <= 300`.
   - **Universes.** Roots of orders 2 to 8 use all labelled tournaments of orders 1 to 7 as bases (2,097,152 at
     order 7). Roots of order 9 use the 6,880 classes. Roots of order 10 use all 1,761,280 one-vertex
     extensions. The order-10 run took 7 min 54 s.
   - **Order 9.** The table is identical to `out/blocks_roots_order9.txt` (87,182 types).
   - **Order 10.** The table is identical to `out/blocks_roots_order10.txt` (243,055 types). The union over all
     orders 2 to 10 is still 243,055, so the packet's omission of bases of order `<= 6` changes nothing.
   - **Order 8.** The cone offsets from the labelled universe are identical to `out/cones_roots_order8.txt`.
2. **Class list: `canon8.c`.** The 6,880 lines of `out/classes_order8.txt` are tournaments and pairwise
   non-isomorphic (minimum code over all 8! relabellings). Since A000568(8) = 6,880, the list is complete, and its
   one-vertex extensions contain every tournament of order 9.
3. **Census claims, read from the files and from the audit's own tables (`analyse_mine.py`, output in
   `analyse_mine.out`).**
   - **Cone offsets.** Every odd `r` in `[-791, 647]` has a cone witness line. The least exponent over this range
     is 13, attained at `r = -769, -765, -709, -701, -627, -611, -441, -301, -295, 511`.
   - **Gaps.** The only odd `|r| < 1000` with no cone root at any order `<= 10` (orders 2 to 7 included) are
     `-891, -809, -793, 649, 825`. This was confirmed by exhaustive search, not only by reading the files.
   - **Box `c <= 31`, `|r| <= 275`.** Complete, with maximum exponent 12. It fails first at `(31, 277)`.
   - **Box `c <= 15`, `|r| <= 295`.** Complete, with maximum exponent 13. It fails first at `(3, 297)`.
   - **Box `c <= 63`, `|r| <= 9`.** Complete. It fails at `(35, 11)`.
   - **`(245, 1)`.** Absent.
   - **Exponents for Theorem 3.** For `1 + 2^a`, `a = 1..10`: `5, 5, 4, 7, 6, 8, 10, 10, 10, 11`. For `2^a - 1`,
     `a = 1..10`: `0, 5, 5, 7, 5, 6, 7, 10, 13, 12`. Both lists are as stated.
4. **Literal recounts: `averify.c`.** The program builds `L = B + v` and then `L o 0^t` as in THM-5126 (the chain
   `z_1 -> ... -> z_t`, `B` beats every `z`, every `z` beats `v`). It counts Hamiltonian paths by subset DP on an
   adjacency matrix, for `t = 1, 2, 3`. No failure in any batch:
   - 67 cone offsets: 201 counts. These include `-791`, `647`, `3`, `7`, `21`, `27`, `257`, `513`, `1025`,
     `2^a +- 1`, and 45 random offsets in `[-791, 647]`.
   - 72 types `(c, r)` with `c = 3..31`, spread over the box, taken from the packet's two type-witness files: 216
     counts.
   - All 293 box types that have no witness line on file (see M4), using the audit's own order-10 witnesses: 879
     counts. Inputs and outputs: `nowitness_recount.in` and `nowitness_recount.out`.
5. **One back arc: `aoneback.py`.**
   - The closed forms were checked by the audit's own subset DP for `3 <= m <= 9` (84 cases), with no failure.
   - Among all `|O| = 1` roots over these bases, the power-of-two responses occur exactly at `(k = 2, y < i)`,
     offset 5, and `(k = 3, a = b = 1)`, offset 9.
   - Diophantine checks by brute force: for `k < 200`, `3^k + 2^(k+1) - 1` is a power of two only at `k = 2`. For
     `a, b < 80`, `2^(b+1) + 2^(a+1) 3^b - [a = 0]` is a power of two only at `(1, 1)`.
   - The valuation proof was re-derived by hand and is correct. It needs only that `2^(k+1)` divides `3^k - 1`.
     When `a != b` the bracket is odd and greater than 1. When `a = b`, `3^b + 1` is a power of two only for
     `b <= 1`.
6. **Templates: `ablow.c`.** The five example families of section 5 were rebuilt from the order-7 and order-8
   output lines and blown up literally for `1 <= a, b <= 4` (80 counts, no mismatch). From the output files: 52,
   200 and 59 distinct unit polynomials at orders 7, 8 and 9 (as stated). No unit template has `gamma = 0` and
   `alpha beta != 0`, and none has the shape `1 + 2^j v + 2^i uv` (or its mirror).
7. **Roots with nonempty `O`: `acensus_h.c`.** Every odd `N <= 400` other than 1, 7 and 21 is `H(L o 0)` for some
   root with `O` nonempty and `B` of order 6 or 8 (see M3).
8. **Git, read-only.**
   - The pre-registration commit `18571fc3da` (18:39:09) is on `origin/main` and predates the output commit
     `aaf352fa43` (18:49:29).
   - The prediction text is unchanged between the two commits. Only the Outcomes section was filled in.
   - An untracked `out/types_c9_to_c31_roots_order10.running` was present during the audit.

## MAJOR findings

**M1. Section 5's two-parameter template census and the exclusion "no independent pair" are prior work. THM-5145
types both as NEW.**
- **Where:** THM-5145 status, "NEW ... the census (5) and its two exclusions"; README section 5, first bullet.
- **The prior census.** kind-pasteur S21 (`kps_s21_two_adic_lines_20261110`, README section 4 table row
  "two-module families (`modfam.c`)") is the same blow-up `K[TT_a at x, TT_b at y]`, fitted as
  `al + be 2^a + ga 2^b + de 2^(a+b)`. It covers every class of order `<= 8` and every ordered pair of vertices
  (406,412 triples), finding 3,477 distinct exact tuples. It states that "none is additive (`de = 0` with
  `be, ga != 0`)". That covers every offset, so it is stronger than the packet's unit-offset case.
- **Also prior.** That packet's audit (AUDIT.md lines 179 and 252) proves exactness on the grid and records
  "`de > 0` iff both on cycles". THM-5124 (opus S90, two-module additive rigidity) is credited there as "no additive
  pair".
- **Fix.**
  - Credit `modfam.c`, THM-5126's packet and THM-5124.
  - Restrict the NEW claim to what is new: the count of unit polynomials, the targeted order-9 census, and the
    weight-three-shape exclusion `1 + 2^j v + 2^i uv`. No earlier statement of that exclusion was found.

**M2. Negative-offset cones and rays, and the columns `r = +-1` of the type conjecture, are partly on record.**
- **Where:**
  - README section 1 ("one negative offset (`2^K - 159`, THM-5085)").
  - THM-5145 status, "NEW ... numbers beginning with the binary digits 11 reached through negative offsets".
  - HYP-9675.
- **On record:**
  - **THM-5130 / kps_s22 (README lines 24, 127; THM-5130 (D), (E)).** A pair of order 8 with offset `-1` and
    slope `2^10`, giving `2^(10+t) - 1`. Rays `2^s - 1` from one pair. Depth-3 towers with power-of-two slope
    (`2^N - 1`). Anti-seeds `c 2^s - 1` for every odd `c <= 4097`.
  - **THM-4748 (4).** The numbers `2^k - 1` are known path counts.
  - **Session195 A/K/D rays (cited in kps_s22 README line 225).** `66 * 2^k - 17` and `88 * 2^k - 13`, i.e. types
    `(33, -17)` and `(11, -13)`.
  - **THM-5129 (B).** Complete slope lists, at every order, for one-vertex sources with offsets `D <= 73` and
    `D = 129, 257`, including "absence of a dyadic slope" at 129 and 257. These are cone roots with `|O| = 1`.
  - **HYP-9662 (A).** This is exactly the `r = +-1` columns of HYP-9675 and is not cited there.
  - **Theorem 3, second family.** `2^b + 2^a - 1 = (2^(b-a) + 1) 2^a - 1` lies on the anti-seed line
    `c = 2^(b-a) + 1`, so for `b - a <= 12` its ray part is THM-5130 (E).
  - **Tooling.** kps_s22 `scripts/portrait.c` already computes `(r, M)` for every `(B, O)` and any `|r| <= RMAX`.
    It was run only at `RMAX = 1`.
- **Not found on record:** the complete census of cone roots and of general types `(c, r)` at order 10; Theorems
  1, 2 and 3 as stated; the cone gaps 649, 825, -793, -809, -891.
- **Fix.**
  - Replace "one negative offset" with the list above.
  - Drop "numbers beginning 11 reached through negative offsets" from the NEW list. The phenomenon is THM-5085 and
    THM-5130; what is new is the complete interval.
  - Cite HYP-9662 (A) and THM-5130 (E) in HYP-9675 and in Theorem 3.
  - Cite THM-5129 (B) next to section 4 (b).

**M3. The equivalence "quantitative type conjecture <=> Chan-Pak" is unsupported, and false as literally stated.**
- **Where:**
  - HYP-9675 status: "The quantitative form is equivalent to the Chan-Pak conjecture (every path count is
    H(L o 0) for some root)".
  - README section 6, opening paragraph and the HYP-9675 paragraph ("Chan-Pak ... needs ... for every odd N a
    splitting ...").
- **Why the justification fails.** The argument adds a sink `s` and a vertex beating only `s`. That gives a root
  with `O = {}`, so `M = 0`, and such a root has no type. With `M = 0` allowed, "every odd N is `c0 + M 2^t`" is a
  trivial restatement of Chan-Pak.
- **Counterexample at `N = 1`.** The quantitative form needs `M = c 2^j >= 1`, and it fails at `N = 1`. In
  `L o 0^t`, the vertex `z_t` has out-degree 1. In a transitive tournament the unique vertex of out-degree 1 points
  to the sink, so `v` is a sink and `O` is empty. Hence `H = 1` has no representation with a type.
- **The converse direction.** It needs every path count `N >= 3` to have a realiser `L o 0` with `O` nonempty.
  This is true for all odd `3 <= N <= 400` other than 7 and 21 (audit item 7), but it is not proved.
- **Fix.** State only "the quantitative form implies Chan-Pak". Then either prove the representation with
  `O != {}` for `N >= 3`, or label it as checked to 400 and OPEN. In README section 6, add that the reformulation
  is vacuous when `M = 0` is allowed.

**M4. Verification is overstated for Theorem 2.**
- **Where:**
  - README table row 2 ("PROVED (same)", where row 1 reads "every witness recounted literally").
  - README section 2, "Literal verification".
  - THM-5145 status "Machine checks".
- **What is missing.** The boxes `c <= 31, |r| <= 275` and `c <= 15, |r| <= 295` contain 293 types, all with
  `c = 9..31`, whose only roots are of order 10. The packet has no witness line for them, because witness files
  exist only for `c = 3..31` at order 9 and `c = 3, 5, 7` at order 10. So they were never recounted literally.
  Their breakdown by `c` is in `analyse_mine.out`.
- **The proof stands.** The tail law (1') is audited, and the census is reproduced independently. The audit
  recounted all 293 literally, using its own order-10 witnesses (879 counts, 0 failures,
  `nowitness_recount.in`/`.out`).
- **Fix.** Add a witness file for these types: the packet's `.running` job, or the audit's file. Or reword the
  verification claim. Only after that can the theorem be said to rest on literally recounted witnesses.

**M5. Section 4 is typed NEW as a whole, but parts are classical or known.**
- **Where:** THM-5145 status, "NEW ... (4)".
- **Classical.** `H(B) = 1 + 2^k` for a transitive tournament with one back arc (THM-4748, notation `NT_s`; Moon
  1972; THM-250).
- **Known.** The three-slot origin of the factor `3^k` is THM-4879 (B) / THM-5109 (B): the coefficient
  `2^(...) 3^(z-u-1)`, "each vertex strictly between an earlier entry and a later exit has three admissible
  slots". The two cones of (b), `(5, 0, 32)` and `(9, 2, 16)`, are Session236's seeds (offsets/PROOF.md section 2;
  THM-5109 (C)).
- **New, as far as the search shows:** the closed forms of `U_B(y)` and `W_ij` for this base, and the
  classification (b).
- **Fix.** Credit these sources.

## MINOR findings

- **m1. README section 2.** "So the census is over all rooted tournaments of order at most 10" is inaccurate: bases
  of order `<= 6` (roots of order `<= 7`) were not run. The audit ran them, and they add no type and no gap. Fix:
  add them, or say "orders 8 to 10 (smaller roots add nothing; audit)".
- **m2. Window lemma.** The bound `H(B) <= U_B(y)` is also THM-5130 (F) ("`U_B(o) >= H(B)`") and kps_s22
  SOURCES_AND_DOMINATION Proposition 4. Credit them next to THM-5109. The proof given is correct.
- **m3. Section 4 proof.**
  - The proof is a one-sentence sketch backed by a check for `m <= 11`, yet the status says PROVED for all `m`.
    Write out the count. Example for `y > j`: `P2` is increasing inside `{y, ..., m-1}`. `P1` is a Hamiltonian path
    of `B[{0..y-1} u S]`, and there are `1 + 2^k` of them for each of the `2^(m-1-y)` choices of `S`. Treat the
    other four cases the same way.
  - In the number theory, say that divisibility forces `v_2(3^k - 1) >= k + 1`, which gives `k <= 2`, and that
    `k = 1` fails because `v_2 = 1`.
- **m4. Section 5 hypothesis.**
  - "If no Hamiltonian walk of K visits x (resp. y) three times" is not the right hypothesis for two modules. The
    polynomial needs `w_jk = 0` whenever `max(j, k) >= 3`, and a walk could visit `x` three times only when `y` is
    visited twice.
  - The program's test is nevertheless sufficient: the four fit points, then an exact check at `(3, 3)`. At
    `(3, 3)` the defect is a positive combination of the `w_jk` with `j, k <= 3` and `max = 3`. The lemma of
    THM-5126 (B) (audit A2), applied to `K[TT_b at y]` and `K[TT_3 at x]`, then propagates to `j` or `k >= 4`.
  - Add this argument. The checks at `(4, 2)`, `(2, 4)` and `a = 4` are redundant.
- **m5. Order-9 targets.** The list in `out/order9_target_values.txt` is `1 + 2^p + 2^q` with `1 <= p <= 6`,
  `2 <= q <= 11` only. Ten candidates up to 2113 were not searched: 257, 385, 513, 641, 769, 1025, 1153, 1281,
  1537, 2049. In addition, 49 to 73 are below `f(9) = 75` and so vacuous at order 9. Say so where the 38 values are
  mentioned (README section 5, THM-5145 (5)).
- **m6. Section 5, "every candidate with `H < 45`".** Correct, because `f(8) = 45` and the floor is increasing
  (THM-4810 gives `f(9) = 75`; floor law of Moon and Busch). Cite the floor at orders `>= 9`, not only 3 to 8. The
  "in particular 19, 35, 37" is a subset: 33 and 41 are also candidates, and 37 also arises as `(j, i) = (4, 0)`.
- **m7. PREREGISTRATION line 29.** "at order 9 the positive offsets below 1000 with a cone root are fewer than 181"
  is wrong. All 181 are below 1000 (largest 597). This is harmless to the scoring of the predictions.
- **m8. Pre-registration outcomes.** These were checked against the files and the audit's census. Q1, Q2 and Q4
  are CONFIRMED and Q3 FAILED, as stated. The in-sample statements are exact: 181 positive and 548 negative
  offsets; missing `27, 79, 99, 155, 195, 199` and `-117, -169`; `1 + 2^e` for `e = 1..7, 9`.
- **m9. Head.** Optionally note that the verified head reaches `2^35` (THM-4875 addendum, as recorded in THM-5126).
  Nothing changes, since `2^32` suffices.

## Section-by-section notes (checked, no finding)

- **Theorems 1 to 3: logic.** The tail law needs `e >= j + 1`. When `e <= j`, `N` is at most `2^13 + 647`, or
  `31 * 2^12 + 275`, or `15 * 2^13 + 295`, or `1 + 2^10 + 2^13`, all below `2^32`. Parity: `N` odd and `e >= 1`
  give `r` odd. Positivity is a hypothesis. Rays never contain 7 or 21, by THM-1370 together with the tail law. The
  HYP-9658 (C) remainder `a >= 11`, `b - a >= 13` is right, because the lines `c = 2^d + 1 <= 4097` mean
  `d <= 12`.
- **Window lemma.** The proof is correct. Cutting a Hamiltonian path before `y` is injective, `W >= 0`, and
  `c0 + 2M = H(L o 0) >= 1`.
- **Density remark.** Correct. Theorem 1 has about `720 log_2 X` members up to `X`.
- **Section 4 (a).** `c0 = 1 + 2^k - 2^k = 1` and `M = 2^s (2^(k+1) - 1 + 2^k + 1) = 3 * 2^(s+k)`. Correct.
- **Section 5.** "A non-product unit template lives in a strong tournament" holds in the sense that it lives in one
  strong component. If `x` and `y` are in different components, `H` factors as a product, and `P(0, 0) = 1` forces
  the other components to have `H = 1`. The H formula for the weight-three shape is `P(2, 2)`. The shapes
  `(0,0)`, `(0,1)`, `(1,1)`, `(1,2)` contain 7 or 21, checked.
- **Section 6.**
  - **Route 1.** The argument is correct as worded ("for every `T`"): the `P_j(TT_a)` with `j >= 2` grow with
    `a` while `H = 1`. The remark about 7 and 21 is correct.
  - **Route 2.** The `N = 1 (mod 8)` implication is correct: take `k` equal to the popcount of `(N - 1)/8`.
  - **Routes 3 to 6.** These are honest and scoped. The overclaim is M3.

## Could not check / not rechecked

- The totals of single and pair templates in the section 5 table (1,006 / 1,446; 8,832 / 8,481; 28,859 / 23,605)
  and the 27,315 strong order-9 tournaments were not recomputed. The audit only read the output files: their
  internal counts, the distinct unit polynomials, both exclusions, and the five families recounted literally.
- `sparse_cones.c` (the offsets 3, 17, 33, 65 with two and three back arcs) was not rerun.
- The prior-art search used grep over `01-canon/theorems`, `05-knowledge/hypotheses` and the result folders dated
  2026-11-09 to 2026-11-11. There are about 50 folders on 2026-11-11 alone, and only the hits were read in full
  (THM-5109, THM-5126, HYP-9658, THM-5129, THM-5130, THM-4879, THM-4748, THM-5124, Session236 offsets/PROOF.md,
  kps_s19/s21/s22 READMEs, odd_c_seed_symmetry). A census in a folder that uses none of the searched words could
  have been missed.

## Files

- `acensus.c`: the audit's own census.
- `acensus_h.c`: the same census, with the `H(L o 0)` values for nonempty `O` added.
- `averify.c`: the literal recount.
- `canon8.c`: the completeness check of the order-8 class list.
- `aoneback.py`: the section 4 check.
- `ablow.c`: the template blow-ups.
- `parse_witness.py` and `analyse_mine.py`: the analysis scripts. `analyse_mine.out` holds the analysis output.
- `ord2.txt` ... `ord10.txt`: the type tables, with `W` witness lines.
- `sample_*.in` / `.out`, `nowitness_recount.in` / `.out`: recount inputs and outputs.
- `hvals_Onz.txt`: the `H(L o 0)` values.
