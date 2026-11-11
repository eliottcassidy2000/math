# Independent audit of kind-pasteur-2026-11-11-S23 (THM-5133, HYP-9664, packet `kps_s23_flip_spin_kernels_20261111`)

Audited state: commit `bdb2f8edd8` (README.md, PREREGISTRATION.md, THM-5133, HYP-9664, scripts/, runs/, data/).
Auditor: independent adversarial pass, 2026-11-11. All numbers below were recomputed with programs written for this
audit (directory `audit/`); nothing in `scripts/` was imported or run. The author's scripts and run files were read,
and outputs were compared.

## Verdict: PASS WITH CORRECTIONS

No error was found in the mathematics of the main results: every FINITE-EXACT number of the packet was reproduced
(one exception in prose, minor 10), every proof was checked line by line, and Theorem 3.4 (PG(3,2)) was
confirmed by a method that does not use the author's lemma, his orbit argument or his programs. The corrections
concern one hypothesis that is false (HYP-9664 (C)), three places where known results are presented as specific
to this packet or as new path counts, one statement that is false in a boundary case, and wording.

Five majors, 23 minors.

---

## MAJOR findings

### M1. HYP-9664 (C) is false for n = 5 and n = 6 (its own "first test" refutes it)

*Where:* HYP-9664 title, status (C) and section (C). README section 3.4 (last sentence, "So the doubling is
extremal among sign patterns with quaternionic lines") and THM-5133 (F5) are statements about `n = 4` and are
correct there; their wording must say so.

*What is wrong:* (C) says that for every `n >= 4` the Cayley-Dickson tournament minimises the number of even
zero-sum 4-sets among tournaments on `PG(n-1, 2)` with all lines cyclic (630 at `n = 5`, 6,510 at `n = 6`) and that
at most `[n choose 3]_2 - 7 [n-1 choose 3]_2` planes are odd at once (50 at `n = 5`, 310 at `n = 6`). Type A of
Theorem 3.4 generalises and beats the doubling: fix a subspace `U` of dimension `n - 3` and ask that all
`2^(3(n-3))` planes skew to `U` be odd, with all lines cyclic. The F_2 system is consistent, and the explicit
solutions, recounted by direct counts of cyclic triangles (`a10_hyp_c_refute.py`, `a10_verify_witness.py`), give

| n | odd planes, type-A generalisation | odd planes, CD_n | even 4-sets, type-A generalisation | even 4-sets, CD_n |
|---|---|---|---|---|
| 4 | 8 | 8 | 42 | 42 |
| 5 | 64 | 50 | 546 | 630 |
| 6 | 512 | 310 | 5,298 | 6,510 |

all with every line cyclic and regular scores. The test proposed in the hypothesis ("decide by linear algebra over
F_2 whether more than 50 of the 155 planes can be odd at once") is a few seconds of linear algebra. A random
greedy search already finds 55 odd planes on average and up to 64 (`a10_hyp_c_probe.out`).

After the audited commit the author's working tree gained `scripts/cd_pg42.py`, `runs/cd_pg42.txt` and
`data/pg42_71_odd_planes.txt` (files dated 10:48, not in `bdb2f8edd8`): 71 odd planes at `n = 5` (the 64 planes
skew to a line and the 7 planes containing it), 504 even 4-sets. I recounted that certificate directly
(`a10_verify_author_71.py`): 155/155 cyclic lines, 71 odd planes, 84 planes of defect 6, 504 even 4-sets. So the
refutation is confirmed twice, by two different constructions.

*Correction:* mark HYP-9664 (C) REFUTED for `n = 5, 6` (both clauses), with the certificates; keep only the
`n = 4` statement (already Theorem 3.4 (5): the sedenions tie with type A). Restrict "the doubling is extremal
among sign patterns with quaternionic lines" to `n = 4` everywhere. If a replacement hypothesis is wanted it must
start from the values 71 (`n = 5`) and at least 512 (`n = 6`); the only upper bound I know at `n = 5` is 82
(each of the 31 solids holds at most 8 odd planes and each plane lies in 3 solids).

### M2. `H(CD_4) = 198,335,025`, its 21 automorphisms and the DRT tower are already in canon (THM-447/448/450); listed as "specific here"

*Where:* README section 3.5 item 4 and the closing paragraph of section 3.5 ("What is specific here: ... `H(CD_4)`");
THM-5133 title, (G) and status (G); README section 0 item 1 and section 6 row 1 ("the exact dictionary is
`2^n - 1` vertices").

*What is wrong:* THM-448 (kind-pasteur-2026-06-09-S1, "The DRT/Mersenne doubling tower") already has a doubly regular
tournament `T_15` with `H(T_15) = 198,335,025` and `|Aut| = 21`, proves that every member of its tower on `2^k - 1`
vertices is doubly regular, and records `T_7 = Paley`. THM-450 (5) reads that tower as "a CD tower on tournaments
... this one lives at `2^k - 1`"; THM-585 uses `15 | 198335025`. The audit found an explicit isomorphism
`CD_4 = T_15` (own backtracking and, separately, networkx VF2; `a8_thm448.py`, `a8_thm448e.py`); `CD_4` is not
isomorphic to the converse of `T_15`. So `H(CD_4)`, the automorphism count and double regularity at `n = 4` are
KNOWN facts about a tournament that is already in the author's own canon, and the "`2^n - 1` vertices" reading is
THM-450's. What the packet could claim is the identification `CD_4 = T_15`, which it does not state. Grepping the
constant `198,335,025` finds THM-448.

New information from the audit: the two towers part at `n = 5`. `CD_5` is not isomorphic to `T_31` of THM-448 nor
to its converse (different vertex profiles of triple out-intersections; VF2 agrees), although both have
`|Aut| = 21`, the same global triple-intersection spectrum `{0:28, 1:84, 2:252, 3:3192, 4:812, 5:84, 6:28, 7:15}`
and the same multiset of neighbourhood path counts up to exchanging out and in (`a8_thm448c/d.out`).

*Correction:* in README 3.5, THM-5133 (G) and its status: "`CD_4` is isomorphic to `T_15` of THM-448; `H`, `|Aut|`
and double regularity at `n = 4` are KNOWN (THM-447/448/450, THM-585)". Remove `H(CD_4)` from "what is specific
here". Credit THM-450 (5) for the `2^n - 1` reading. State that `CD_5` and `T_31` are different doubly regular
tournaments (audit). What remains specific in 3.5 is the description of the planes (8 Paley, 7 of defect 6).

### M3. "New path counts with proof: the minus grid of d = 31" is an overclaim; the lines are in THM-5130 (E)

*Where:* README section 5 first bullet, section 0 item 4, Proposition 2.5 (the list of lines); THM-5133 status (C)
("New: the minus carrier d = 31").

*What is wrong:* `N-(31, r, s) = (31 2^r - 1) 2^s - 1`. THM-5130 (E) (the author's S22, audited the day before)
certifies an anti-seed for every odd `c <= 4097`, that is every line `c 2^s - 1` from a level `s_-(c)`. For the
lines quoted in Proposition 2.5:

| line `c 2^s - 1` | first level here | first level in THM-5130 (E) (`certs_minus.txt`) | audit |
|---|---|---|---|
| c = 61 (r = 1) | 5 | 5 | tower certificate, quoted |
| c = 123 (r = 2) | 5 | 3 | pair rebuilt, `H = 983, 1967, 3935` recounted |
| c = 247 (r = 3) | 5 | 2 | recounted |
| c = 495, 991, 1983, 3967 (r = 4..7) | 5 | 1 | recounted |
| c = 27 (row b = 1) | 3 | 1 | recounted (also `27 D - 1` in uniform_return_metagraph_20261111) |
| c = 29 (row b = 2) | 4 | 4 | recounted |
| c = 15 (row b = 3) | 6 | 1 | quoted (also Session261 `15 2^m - 1`, `m >= 1`; carrier 21 of OFFSET_TRANSPORT.md) |

(`a3_prior_lines.py`.) So every line of the minus grid with `r <= 7`, and all three extra rows, were path counts
before this packet, from levels not higher than here. Moreover every member below `2^32` is a path count by
THM-4875 (D). What is new: one kernel of order 8 with two free transitive modules that carries the whole grid, and
the lines `c = 31 2^r - 1` for `r >= 8` (`c >= 7935`), of which only the members above `2^32` are new numbers.

*Correction:* README section 5: "A uniform carrier for the minus grid of `d = 31` (one kernel, two free modules).
As path counts, the lines with `r <= 7` and the rows `c = 27, 29, 15` are in THM-5130 (E) from levels
5, 3, 2, 1, 1, 1, 1 and 1, 4, 1; new lines: `c = 31 2^r - 1`, `r >= 8`." Add the same sentence to Proposition 2.5
and to THM-5133 status (C).

### M4. Cayley-Dickson statements: three missing credits, one of them in the paper the README already cites

*Where:* README Proposition 3.5 items 1, 3, 4, and the KNOWN / CITED paragraph; THM-5133 (G) and status (G);
HYP-9664 status (C) ("For n = 5, 6 only the values of CD_n are known").

(a) **The counts for every n are a theorem in Wilmot, arXiv:2505.11747 (v4).** Read with `pdftotext`: Table 5 lists
8, 50, 310, 2,046, ... octonion and 7, 105, 1,085, 9,765, ... quasi-octonion subalgebras for dimensions 16, 32, 64,
128, ...; Theorem 8: `S_(m+1) = 7 (O_m + S_m)`; Theorem 9: closed forms for `O_m`, `S_m`; Theorem 11: the number of
zero-divisor pairs is `(1/2) N (N - 1)(N - 3)`, `N = 2^(n-1) - 1`, that is 84, 1,260, 13,020, ... (its proof uses
an enumeration from an earlier reference for the 12 pairs in one quasi-octonion algebra). In the README's
language: odd planes of `CD_n` = octonion subalgebras, planes of defect 6 = quasi-octonion subalgebras, and
`(1/2) N (N-1)(N-3) = 84 [n-1 choose 3]_2`. So "50 odd planes, 105 of defect 6" at `n = 5`, "6,510 even 4-sets" at
`n = 6` and the pattern `42 [n-1 choose 3]_2` are CITED for all `n`, not a FINITE-EXACT pattern for `4 <= n <= 6`;
the README cites this paper only for 84 and 1,260.

(b) **The block recursion of the sign table (3.5.1) is in canon:** Session168
`session168_computational_structures_20261103/cayley_dickson/NOTE.md` eq. (1) (2026-11-03), with the same doubling
formula, in Boolean form `f_(u,r)(y,t) = f_u(y) + r delta(y) + t (c_u(y) + r)`; the four cases `r, t = 0, 1` are
the four blocks `G, G^T, G D, -G^T D`. Not cited.

(c) **`G_n` is a Hadamard matrix (3.5.3):** B. Kivunge, "Sign matrices for frames of 2^n-ons under Smith-Conway
and Cayley-Dickson multiplication" (abstract, Mile High Conference on Nonassociative Mathematics, 2013) states it
for the same doubling formula `(a,b)(c,d) = (ac - conj(d) b, da + b conj(c))`. The README says a web search "did
not locate statement 3"; one query locates it. (No novelty was claimed, so this part is a citation repair.)

*Correction:* cite Wilmot Table 5 and Theorems 8, 9, 11 for the counts at every `n` and relabel item 4 "CITED
(all n); recomputed here for n <= 6"; in HYP-9664 replace "only the values of CD_n are known"; cite Session168
eq. (1) as KNOWN for 3.5.1; cite Kivunge for 3.5.3 and delete the "not located" sentence in README and THM status
(G).

### M5. Proposition 2.3 / THM-5133 (B): "elliptic iff r > 0" is false when O is empty

*Where:* README Proposition 2.3 ("Hence `K` is elliptic iff the offset `r` is positive ..."); THM-5133 (B) ("It is
elliptic iff `r > 0` and hyperbolic iff `r < 0`").

*What is wrong:* THM-5130 defines a pair `(B, O)` with `O` "possibly empty". For `O` empty, `M = 0`, `r = H(B) > 0`
and `det G = r M = 0`: the kernel is join-like, not elliptic. Smallest case: `B_0` a point, `O` empty, `K = TT_4`
with the two middle vertices as ports, `G = [[1, 0], [0, 0]]`. In `a1_visit_matrix.py` the claim fails on exactly
the 28 sampled pairs with `O` empty and on no other. `det G = r M`, "at most two visits" and "`r` is odd, never 0"
are correct (140 random pairs with the state computed from the definitions of THM-5130; 0 failures).

*Correction:* add the hypothesis "`O` nonempty (so `M >= |O| H(B) >= 1`)" in both places, or state "`det G = r M`;
for `O` nonempty the kernel is elliptic iff `r > 0` and hyperbolic iff `r < 0`; for `O` empty it is a join".

---

## MINOR findings

1. **THM-5133 title (D), README section 0 item 2.** "`(e_a + e_b)(e_c +- e_d) = 0` iff `a + b = c + d` and the 4-set
   contains an even number of cyclic triangles" is false for the third pairing of an even 4-set (only two of the
   three pairings give a zero product). The body (D2) and Proposition 3.2 are correct. Write: "some pairing of
   `{a, b, c, d}` gives a zero product iff `a + b + c + d = 0` and the 4-set contains an even number of cyclic
   triangles, and then exactly two pairings do".
2. **README section 0 item 2.** "The octonions are the tournaments on the Fano plane in which every line complement
   holds exactly one cyclic triangle": there are 1,024 such tournaments; the octonion tables are the 16 with all
   lines cyclic (the other 1,008 have 1, 2 or 4 cyclic lines: 560, 336, 112). Section 3.3 has it right.
3. **Proposition 3.1, proof.** "replace `e_c` by `-e_c` ...; this reverses the two arcs at `c`": it reverses all
   three arcs of the line (the arc `{a, b}` as well; it is the switching at `c` times the coboundary at `c`, as
   section 3.3 says). Reversing only the two arcs at `c` would turn a cyclic triangle into a transitive one. The
   reduction is valid because all three arcs are reversed. The rest of the proof is correct.
4. **Theorem 2.1.** "`N_ij = 0` unless `1 <= i, j <= m - 1`": walks with `i = 0` or `j = 0` exist (for `C_3`:
   `N_00 = N_01 = N_10 = 1`); they are not used. Say "the sum runs over `1 <= i, j <= m - 1`". The bound `m - 1`
   is correct and attained at every order 3..7 (audit).
5. **Theorem 2.1, status.** The formula is THM-3121 eq. (4) with `m - 2` singleton factors, not only "the
   principle"; the two-port expansion `sum c_rs C_r(A) C_s(B)` and the factorisation for one kernel
   (`37 H = x x' + 329 H H'`, `329 = det`) are in Session261 BILINEAR_KERNEL.md (2026-11-11 02:20, credited);
   `de = N_22` and the two-visit criterion are THM-5126 (E); the S21 audit computed the `c_ij`. Label 2.1 KNOWN and
   list as specific: the general identity `N_22 H = L L' + det G H H'`, the trichotomy, the census.
6. **THM-5133 (A).** "If no walk visits a port three times" should be "three or more times" (README Corollary 2.2
   has it right; the S21 audit flagged the same wording in THM-5126 (E)).
7. **Proposition 2.3, proof of "at most two visits".** One sentence about the closing law is not a proof. Write:
   no arc goes from `X` to `B_0` and `v` is used once, so a Hamiltonian path leaves `X` at most once; it meets `X`
   in at most two segments, each of which runs through `A` and then through `B`.
8. **Table of 2.4 and THM-5133 (C).** The column "parabolic" excludes the join-like kernels, which are parabolic
   by the definition (`det G = 0`). Head it "parabolic, not join-like" (with them: 6, 24, 80, 278, 1,108, 5,610,
   43,110). README section 0 item 3 "parabolic (one product: joins)": 12, 54, 218, ... parabolic kernels are not
   joins.
9. **Proposition 2.5, "exactly".** The statement is about the values of `d`. At order 8 each of the three forms
   occurs in exactly two classes, a class and its converse (audit: codes isomorphic to 17082497, 172278145,
   71647617 and their converses). Say so, and say that counts are over (class representative, ordered port pair).
10. **README section 4, F1.** "one class of order 7 with `H = 129`, `Pi_2 = 3 2^9`": there are four classes of
    order 7 (and 8 and 12 at orders 8 and 9, the numbers expected from adding sources and sinks). The levels
    1, 2, 3, 7 are correct.
11. **F2, THM-5133 (H).** The restriction "chain length at most 22" is vacuous: in an irreducible word of length
    `L >= 3` the pair `(1, L)` alone contributes `2^(L-2)` and the second letter gives another inversion, so
    `H > 2^(L-2) + 1`, and `H = 2^s + 1` forces `L <= s + 1` (attained at every level up to 15). The run therefore
    covers all one-apex words at the levels 3..21. State the stronger fact.
12. **README section 4.** The congruence is step 3 of THM-4756 (the reading modulo 3 of the identity of step 2):
    cite "steps 2-3".
13. **README section 4, S-unit remark.** It needs `O` nonempty (`u_2 != 0`; for `O` empty `Pi_2 = 2 r_2 2^t`) and
    one line on vanishing subsums. The formula `Pi_2 = 2 [u 3^t + r_2 2^t + s]` was checked (25 pairs).
14. **HYP-9664 status (A).** "d = 43 (THM-5130 (F), ...)" should be THM-5130 (H3)(iii). The remark on excluded
    addresses should include `N-(3, 2, 1) = 21`, which lies inside the range `d >= 3` of the hypothesis.
15. **README 2.4.** "the table of S22's `modfam.c`": `modfam.c` is in the S21 packet (section 7 says so).
16. **README 2.4, "the grid of Sessions 260/261".** The unique address `N(d, r, s) = 1 + 2^s + d 2^(r+s)` is R9 of
    `seed_recursion_births_20261111/PROOFS.md` (2026-11-11 00:44), which also closes the grids 3 and 21.
17. **Closing law.** `H = r H(X) + M Pi_2(X)` for an arbitrary module is also eq. (2) of Session193
    `profile_substitution/PROOF.md` (2026-11-05), earlier than Session216. Add it where the law is cited.
18. **Packet.** `runs/ker_7.txt` is listed in section 7 (`m = 3..7`) but is not in the packet. The order-7 row is
    in `ker_census_le8.txt` and was reproduced.
19. **PREREGISTRATION.md caveats** (labels are correct, see below). Add: (a) K2 was derivable by hand exactly like
    K1: working claims 1 and 2 of the file give `H(C_3[1, A => B, 1]) = H H' + Pi Pi'`, the form of discriminant
    -4; (b) C1: the tournament is `T_15` of THM-448, doubly regular in canon since June (M2), and the Hadamard
    property is in the literature (M4 (c)); (c) F2 covers all one-apex words (minor 11).
20. **Concurrency.** Session263 (`session263_spin_flip_20261111`, checkpoints 09:20-09:36, apparently the same
    pasted prompt) and Session264 (09:48-10:04) are not mentioned. I found no mathematical overlap (physical
    source algebras, selectors, Fourier filters; multi-return seams), but the README should name them as
    CONCURRENT and say so.
21. **"What is specific here: the tournament formulations (Propositions 3.1-3.3)".** The sign tournament of the
    octonions, the action of `e_u -> -e_u` on line orientations, the Hamming coset of 16 and the split 16 / 112 are
    in `arithmetic_braids2_20260917_fano_code.md` and S13-S16 of the 2026-09-22 note (credited for the count
    only). The second half of C2 follows from them in two lines (reversing one line of a gauge orientation makes
    the six other line complements even). No literature search is recorded for Theorem 3.4; "specific here" should
    read "not found in the repo; literature not searched".
22. **README section 4.** "the bound `n + 1` ... is attained only by the point and `C_3`" is a statement about
    orders at most 9 (FINITE-EXACT); say so.
23. **README 3.5, citation for 84.** The number 84 is stated in Wilmot (abstract and Theorem 11; "12 zero divisor
    pairs for each quasi-octonion", 7 x 12). The abstract of de Marrais (math/0011260) gives 42 assessors and 168
    primitive unit zero divisors; I did not find 84 there (body not read). Attach "84" to Wilmot (and Cawagas),
    "42 assessors" to de Marrais.

---

## Proof checks

- **Theorem 2.1.** Correct (bijection between Hamiltonian paths and decorated walks; bound by counting the
  non-port vertices between two visits). See minors 4, 5.
- **Corollary 2.2.** The four coefficients, `al de - be ga = det G` and the factorisation were re-derived by hand
  and checked on 282 cases. The reading "sum / difference of two products" needs `N_22 > 0`; no kernel of order at
  most 7 has `N_22 = 0` and `N_12 N_21 > 0` (in agreement with THM-5126 (E)), so it is never vacuous there.
- **Proposition 2.3.** `det G = r M` and the two identities are correct; see M5 and minor 7.
- **Proposition 3.1.** Correct for line-associative algebras; see minor 3. The subgroup statement for split forms
  has a two-line proof (`chi(a) = -e_a^2` is a character: `e_(i+h)^2 = -g e_i^2`), so it need not be FINITE-EXACT.
- **Proposition 3.2.** Correct; the four score sequences were checked by hand.
- **Proposition 3.3.** Invariance proof correct; rescaling = switching times coboundary, verified (21 arcs, the
  three arcs of the seven lines through the point).
- **Theorem 3.4, proof of 4.** The case analysis is complete: `o <= 4` gives `E >= 22` (needs "defect never 1",
  verified by direct enumeration of all `2^21` tournaments); `5 <= o <= 7` gives `E >= 24 - k`; `o = 8` by the two
  types; `o >= 9` impossible. The search of `cd_f5_orbits.py` is correctly designed: the property is invariant
  under linear maps, inconsistency is monotone so the depth-first search is complete, and the 21 vectors of defect
  2 of a plane are exactly the 21 pairs of lines (each once). An independent re-implementation of the same search
  reproduces every count (below). Independently of that lemma, the audit enumerated ALL violation patterns next
  to EVERY one of the 2,835 feasible 5-sets (4.6 x 10^11 patterns, no symmetry used): the least number of even
  4-sets with at least five odd planes is 21, and the values that occur are 21, 25, 26, 27, .... With the bound 22
  for `o <= 4` this proves statement 4 again. By-product: 21 is attained only in type A (`o = 7` forces at least
  42, `o = 6` at least 25, `o = 5` at least 27), so the minimisers are exactly the `15 x 64` patterns of
  statement 3.
- **Theorem 3.4, proof of 5.** Correct.
- **Proposition 3.5, proofs of 1, 2, 3.** Re-derived from the doubling formula; correct (all three kinds of lines
  are covered; the cross term `G D G^T - G^T D G` vanishes because the first row and column are all ones; the
  count `N/4 - 1` follows from orthogonality and regularity).

---

## What was recomputed (audit numbers next to the author's)

| item | author | audit | program |
|---|---|---|---|
| visit-matrix formula, all kernels of order 3-5 and 150 of order 6, random modules | 220 cases, 0 failures | 900 cases (360 on kernels with three or more visits), 0 failures; walks by plain enumeration | `a1_visit_matrix.py` |
| factorisation `N_22 H = L L' + det H H'` | 160 cases | 282 cases, 0 failures | `a1_visit_matrix.py` |
| Gaussian kernel, first hyperbolic kernel | `G = [[5,2],[2,1]]`, `[[11,6],[6,3]]`; 120 cases | same matrices, states (1,1), (-1,3) from the definitions; 300 cases, 0 failures | `a1_visit_matrix.py` |
| census, order 3 | 12 / 12 / 6 / 6 / 0 / 0 | same | `a2_census.c` (plain walk enumeration) |
| order 4 | 48 / 46 / 12 / 22 / 12 / 0 | same | |
| order 5 | 240 / 160 / 26 / 76 / 54 / 4 | same; the four hyperbolic kernels are those of code 209 | |
| order 6 | 1,680 / 630 / 60 / 290 / 218 / 62 | same; all 1,680 rows equal to `runs/ker_6.txt` | |
| order 7 | 19,152 / 2,892 / 162 / 1,242 / 946 / 542 | same | |
| order 8 | 385,280 / 16,962 / 602 / 6,732 / 5,008 / 4,620 | same (own class list, own DP) | `gen_classes.c`, `a8_census_dp.c` |
| order 9 | 13,790,592 / 148,596 / 3,804 / 53,804 / 39,306 / 51,682 | same (191,536 classes generated from scratch) | |
| discriminants of symmetric kernels, order <= 6 | -4, -8, -12, -16, -20, -24, -32, -36, -44, -72; 0; 8, 12, 20, 24, 88, 108, 140, 148 | same | `a2_census_summary.py` |
| carriers through order 9 | plus `d = 3, 21, 35`; minus `d = 31` only | same; at order 8 each in a class and its converse | `a8_census_run.py` |
| minus carrier: `G`, visits, formula | `[[215,124],[108,62]]`, two visits, `31 2^(a+b+1) - 2^(a+4) - 1` (16 cases) | same; 25 direct counts `a, b <= 5` (order up to 16) and 12 non-transitive module pairs, 0 mismatches | `a3_carriers.py` |
| plus carriers 3, 21, 35 | formulas of 2.5 | 25 direct counts each, 0 mismatches; the kernels are the carriers of OFFSET_TRANSPORT.md with a one-vertex tail (3: the converse) | `a3_carriers.py`, `a3_offset_transport.py` |
| grid arithmetic | unique addresses, level rule, lines 27, 29, 15, 61, 123, 247, 495 | confirmed for odd `N <= 200,001` | `a3_carriers.py` |
| zero products, dimension 16 / 32 | 84 / 1,260, the even pairings | 84 / 1,260 unordered pairings, exactly the predicted ones over all quadruples of distinct units and both signs; none outside zero-sum 4-sets; 42 / 294 pairs involved | `a5_cd.py` (own doubling) |
| block recursion, lines, Hadamard, double regularity | `n <= 5` | `n <= 6` (recursion 1 -> 6) | `a5_cd.py` |
| `H(CD_4)`, `|Aut|` | 198,335,025; 21 | same | `a5_cd.py` |
| planes of `CD_4`, `CD_5` | 8 + 7; 50 + 105 | same (`H = 189` / `171`) | `a5_cd.py` |
| even 4-sets of `CD_4, CD_5, CD_6` | 42, 630, 6,510 | same | `a5_cd.py` |
| split forms: line rule, subgroup, cyclic lines, even 4-sets, `H = 45` | 1/7, 7/35, 35/155; 42, 78, 90 | same (all parameter choices, dimensions 4-32) | `a5_cd.py` |
| Fano: all-odd tournaments | 1,024 = one orbit; 16 with cyclic lines | same (direct triangle counts on all `2^21`; orbit by BFS; rank of moves 10) | `a6_fano.c` |
| Fano: distribution of the defect | 1024, 21504, 143360, 322560, 688128, 637952, 282624; never 1 | same | `a6_fano.c` |
| 128 cyclic orientations | 16 (0, 189) + 112 (6, 171) | same | `a6_fano.c` |
| PG(3,2): rank, consistency, kernel | rank 80, inconsistent, 25 = moves | same | `a7_pg32.py` |
| odd plane sets by size | 15, 105, 455, 1365, 2835, 3325, 1395, 135; none of size 9 | same | `a7_pg32.py`, `a7_extra.py` |
| sets of 8 | 15 type A + 120 type B | same; solution dimensions 44 and 25; with cyclic lines 17 and 11 | `a7_pg32.py` |
| type B | 42, defect 6 in each of the seven planes | same, all 120 sets | `a7_extra.py` |
| type A | 21:64, 27:896, 32:5376, 33:448, 36:17920, 37:5376, 39:35840, 40:26880, 41:44352, 42:104448, 43:163968, 44:107520, 45:11200 | same histogram for all 15 sets; the 64 have defect 3 in each plane (checked on one set) | `a7_enum.c`, `a7_extra.py` |
| feasible 5-sets, orbits | 2,835 in orbits 315, 1680, 840 (ranks 55, 53, 52) | same, under the full GL(4,2) (20,160 matrices) | `a7_pg32.py` |
| search next to five odd planes | (98, 1405, 1816, 0), (84, 1038, 3456, 0), (66, 819, 1759, 0); `k <= 3` | same counts with an independent implementation; `k <= 3` also by full enumeration | `a7_f5_dfs.py`, `a7_enum.c` |
| least number of even 4-sets | 21 | 21: full enumeration next to all 2,835 5-sets (minima 38, 36, 21 by orbit); annealing with direct triangle counts, 2,700 restarts, always 21 | `a7_run_all_five.py`, `a7_anneal.c` |
| the minimiser of `runs/cd_minimum.txt` | 21, 6 cyclic lines, `H = 46,700,719` | same by direct recount | `a7_extra.py` |
| sedenions in Theorem 3.4 | type B | type B: the plane `{1..7}` and the planes through the point 8 | `a7_extra.py` |
| `Pi_2 = (-1)^n H (mod 3)` | 198,948 classes, 0 violations | same (own class lists) | `a9_pairs.py` |
| F1 | levels 1, 2, 3, 7; 2 of 9 | same; four classes at level 7, order 7 | `a9_pairs.py` |
| `Pi_2` a power of two | 18 values of `H` | same list (124 classes) | `a9_pairs.py` |
| least `Pi_2 / H`, orders 3..9 | 4, 6.4, 7.2, 9.51, 10.67, 12.98, 14.51 | same | `a9_pairs.py` |
| distinct `H`, distinct pairs | 1,520; 50,173 | same | `a9_pairs.py` |
| one-apex formulas of `words_f2.c` | - | 0 mismatches against direct counts for all 510 words of length <= 8 and 40 longer words | `a9_pairs.py` |
| F2 table | 21 levels | identical, from formulas derived for the audit (subset form) | `a9_words.c` |
| fans, `Pi_2` of closings | - | 60 and 25 cases, 0 failures | `a4_misc.py` |

Class lists: `gen_classes.c` produces 1, 1, 2, 4, 12, 56, 456, 6,880, 191,536 classes (OEIS A000568); the multisets
of `H` agree with `data/cls_m.txt` for `m <= 7`.

---

## Pre-registration

- The first 84 lines of PREREGISTRATION.md are byte-identical to the version of commit `b91157e81d`
  (2026-11-11 09:24:44, contained in `origin/main`); nothing above the outcomes line was changed.
- Labels: F1 FAILED (2 of 9), F2 FAILED (0 of 14, indeed 0 of 19), K1 CONFIRMED (0, 0, 4, 62), K2 CONFIRMED
  (`[[5,2],[2,1]]`), C1 CONFIRMED, C2 CONFIRMED (16 / 112). All six labels are right.
- Caveats 1-4 are candid. Additions needed: minor 19.

---

## Credits opened

| credited source | accurate? | sufficient? |
|---|---|---|
| Session193 `profile_substitution/PROOF.md` Theorem 1 | yes (multi-port transfer on cover polynomials, mirror orientation); its eq. (2) is the closing law | add eq. (2) (minor 17) |
| Session260 `QUASIPERFECT_NORMS.md` eq. (1) | yes: `H(C_3[A, B, 1]) = x x' + 2 y y'` for arbitrary modules | yes |
| Session261 `BILINEAR_KERNEL.md` | yes: two-port kernel, `149 hg + 72 (cg + hd) + 37 cd`, `37 H = x x' + 329 H H'` | it also contains the two-port expansion itself (minor 5) |
| OFFSET_TRANSPORT.md section 4 (6) | yes: carriers (16, 24), (4, 42), (32, 35); the order-8 and order-9 plus kernels are these carriers | yes |
| THM-4756 | yes (step 3) | yes (minor 12) |
| 2026-09-22 octonion note S14-S16 | yes: 16 / 112, 189 / 171 | S13 and the fano_code note should be cited for the sign tournament and the gauge action (minor 21) |
| THM-3121 | yes; eq. (4) is Theorem 2.1 | relabel (minor 5) |
| not credited | THM-447/448/450, THM-585 (M2); THM-5130 (E) for the lines (M3); Wilmot Theorems 8, 9, 11, Session168 eq. (1), Kivunge (M4) | |

Searches made for the items called new or specific: "minus grid", "anti-carrier", families with coefficient 31,
lines `c 2^s - 1` for `c = 15, 27, 29, 61, 123, 247`; "hyperbolic / indefinite / signature" in the packets of
2026-11-10/11 (only s91's signature (2,1) for one-arc joins, credited, and an unrelated cycle-gluing kernel in
Session257); "PG(3,2)", "zero-sum 4-set / quadrangle", "odd plane" in canon and hypotheses (nothing before this
packet); "sign tournament", "sedenion", "Cayley-Dickson" with "Hadamard / doubly regular" (THM-447/448/450); the
constant 198,335,025 (THM-448, THM-585). The minus carrier as a kernel, the trichotomy by the sign of `det G`, and
Theorem 3.4 were not found earlier in the repo.

---

## Not verified

- That the runs were made after the pre-registration was pushed: only the commit order can be checked
  (`b91157e81d` 09:24:44 before `bdb2f8edd8` 10:42:21).
- The literature beyond what was read: Wilmot's paper was read for Table 5 and Theorems 8, 9, 11 only (not its
  proofs); Kivunge's abstract only; de Marrais (abstract only), Moreno and Cawagas not opened. No search of the
  literature for Theorem 3.4 was attempted.
- The two tower certificates of THM-5130 (E) for `c = 61` and `c = 15` were quoted from the audited theorem, not
  rebuilt (the other eight were rebuilt and recounted).
- HYP-9664 (A), (B), (D): open statements, not probed. For (C) the maximum number of odd planes at `n >= 5` is
  open (71 at `n = 5` by the author's later certificate; upper bound 82).
- `H(CD_5)`, `H(CD_6)` (not claimed in the packet).
- The author's scripts were read but, by instruction, not executed; `runs/` was compared with the audit outputs.

---

## External sources read (2026-11-11)

- G. P. Wilmot, "Structure of the Cayley-Dickson algebras", arXiv:2505.11747v4 (https://arxiv.org/abs/2505.11747):
  abstract, introduction, Table 5, Theorems 8, 9, 11 (text extracted with `pdftotext`).
- B. Kivunge, "Sign matrices for frames of 2^n-ons under Smith-Conway and Cayley-Dickson multiplication", abstract,
  Mile High Conference on Nonassociative Mathematics 2013
  (https://www.cs.du.edu/~petr/milehigh/2013/old/abstracts/kivunge.pdf).
- R. P. C. de Marrais, "The 42 Assessors and the Box-Kites they fly", arXiv:math/0011260
  (https://arxiv.org/abs/math/0011260): abstract only (42 assessors; the abstract itself gives 168 primitive unit
  zero divisors, not the number 84).
- Seen in search results, not opened: J. W. Bales, "Cayley-Dickson and Clifford algebras as twisted group
  algebras", arXiv:1107.1375 (https://arxiv.org/pdf/1107.1375).

---

## Reproduce (from `audit/`)

    gcc -O2 -o hp hp.c ; gcc -O2 -o a2_census a2_census.c ; gcc -O2 -o gen_classes gen_classes.c
    gcc -O2 -o a6_fano a6_fano.c ; gcc -O3 -mpopcnt -o a7_enum a7_enum.c ; gcc -O2 -o a7_anneal a7_anneal.c -lm
    gcc -O2 -o a8_census_dp a8_census_dp.c ; gcc -O2 -o a9_words a9_words.c
    ./gen_classes 9 aud_cls                                   # about 1 minute; aud_cls_9.txt (4.4 MB) is not kept
    python a1_visit_matrix.py
    for m in 3 4 5 6 7: ./a2_census m < ../data/cls_m.txt > a2_census_m.out ; python a2_census_summary.py
    python a3_carriers.py ; python a3_prior_lines.py ; python a3_offset_transport.py ; python a4_misc.py
    python a5_cd.py ; ./a6_fano
    python a7_pg32.py                                         # writes a7_problems_*.txt (the 3.6 MB file is not kept)
    ./a7_enum a7_problems_AB.txt 1 1 0 ; ./a7_enum a7_problems_five_reps.txt 1 1 0
    python a7_run_all_five.py                                 # 4.6e11 patterns, about 4.5 minutes on 11 processes
    python a7_extra.py ; python a7_f5_dfs.py ; ./a7_anneal 300 200000 1
    python a8_census_run.py                                   # order 9: about 3 minutes on 11 processes
    python a8_thm448.py ; python a8_thm448c.py ; python a8_thm448d.py ; python a8_thm448e.py
    ./a9_words check 8 > a9_words_check8.txt ; ./a9_words f2 22 ; python a9_pairs.py
    python a10_hyp_c_probe.py ; python a10_hyp_c_refute.py ; python a10_verify_witness.py ; python a10_verify_author_71.py

Outputs are the `.out` files next to the programs. Windows binaries and the two large regenerable intermediates
(`aud_cls_9.txt`, `a7_problems_five.txt`) were removed after the runs. The Python drivers call `hp.exe`,
`a2_census.exe`, `a7_enum.exe`, `a8_census_dp.exe` (MinGW names); on another system compile with `-o name.exe`.
