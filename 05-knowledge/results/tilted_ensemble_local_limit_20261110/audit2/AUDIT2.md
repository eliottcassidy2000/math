# AUDIT 2 of the post-audit additions to the tilted-ensemble packet (THM-5127 / HYP-9659)

**VERDICT: PASS WITH CORRECTIONS** - 1 MAJOR (missed prior art for the two equalities of the parity law (4.5) and
for the exact form of (4.4)), 15 minor; every proof correct, every number reproduced exactly with independent code;
late items L1-L4 confirmed.

**Object audited.** Only what was added after `audit/AUDIT.md`: README.md sections 4.1a and 4.1b, the paragraph
"Explicit polar form (4.4)", and the corresponding sentences of the header, Outcome items 3-4, section 7 and section 9
(README md5 `203ec22ea65d96a8e0166b0b56871f23`, worktree HEAD `994ed3a7f9`; the text of the audited sections is the
one committed in `612fa44c81`); THM-5127 (C) (5), (6), title, status block and FINITE-EXACT bullets (md5 `5bbab7e9...`);
HYP-9659 (B'') (md5 `43308a35...`); PREREGISTRATION.md second block; the programs `biasweights.c`, `merge_bias.py`,
`deletion_pair_corr.py`, `box_norm_check.py`, `oddedge.c`, `oddedge_model.py` and their outputs; and four late items
L1-L4 sent by the coordinator during the audit (`gate_mod4.py`, `affine_span.py`, `digit_degree.py`,
`deletion_rank.py`), which are not yet in the README.
**Date.** 2026-11-10. **Auditor.** Independent adversarial subagent. Every program in this folder was written from
scratch. No packet code was imported or copied; the packet's programs were run only to compare their output with
`out/` (section 7), and the packet's `out/` files were read only for comparison (`x9_compare_author.py`).

## Verdict: PASS WITH CORRECTIONS

- Every proved statement of the audited parts is correct as a statement: Proposition 4.3 (both equalities, the
  Cauchy-Schwarz consequences, the chain), Lemma 4.5, Corollary 4.6 (a)-(c) with its equality cases, and the late
  items L1 and L3. Each proof was checked line by line and each step that can be executed was executed (section 3).
- Every number of the audited parts was recomputed with independent code and agrees exactly, including all order-10
  numbers (section 4): the packet's raw totals at orders 2-10 (signed sum, plain sum, absolute sum, classes and
  labelled weights per rank and per zero / non-zero bias, edge histograms) are identical to the auditor's, field by
  field.
- Corrections required: **one MAJOR** (missed prior art: both equalities of the parity law (4.5) and the exact form
  of (4.4) are in canon since August; only the congruence and its corollaries are new) and **15 minor** (two
  overclaims in wording, one mis-typed "PROVED" in HYP-9659, a missing range, a calibration stated for the wrong
  ensemble, and typing / reproduction items).
- The pre-registration is in order: the P3-P6 text is in commit `23d6c4033f` with the intervals the README quotes, and
  the timestamps are consistent with a run after the push (section 5). P3 and P6 are, however, not independent of P4
  (minor m8).
- Late items: L1, L2, L3, L4 all confirmed; L3 can be sharpened to an exact criterion (section 8).

---

## 1. MAJOR finding

### M1. Lemma 4.5: the two equalities of (4.5) are prior canon, and (4.4) has a known exact integer form (credit)

*Where.* README section 4.1b (Lemma 4.5 "(PROVED)", presented as an addition of this session), section 4 paragraph
"Explicit polar form", Outcome item 3; THM-5127 status block ("NEW here: ... the parity law (4.5) with its
corollary"), (C) (2) and (C) (6); HYP-9659 (B'') ("from the parity law").

*What is wrong.* (4.5) is `sum_{i<j} mu_R(i,j) = sum_i a_i = sum_S H(R[S]) H(R[V \ S]) = 2 (mod 4)`. Both equalities
are in the repository, in files of the ear-response / exposure-capacity lane that the packet already cites in part
(THM-4118) but not for this:

- **First equality** (`sum_i a_i = sum_{i<j} mu`): THM-4118 eq. (2) `F_T(emptyset) = F_T(V) = H` with eq. (4)
  (codex, 2026-08-25).
- **Second equality** (the total of the pair terms is the number of ordered pairs of disjoint non-empty paths
  covering `V`): THM-4111 eq. (6), `sum_(a != b) Q_T(a,b) = (n-1) H(T) + F_1(T)`, `F_1` = orderings with exactly one
  backward adjacency (codex, 2026-08-25; such an ordering cut at its backward adjacency, or a Hamiltonian path cut at
  any of its `n-1` adjacencies, is exactly an ordered two-path cover); restated as THM-4115 eq. (7)
  `W = ((n-1)H + F_1(T))/2` with `W = sum_(i<j) w_ij`, `w = mu/2`; as THM-4181 Theorem 2, eq. (6)
  `U^0 + U^1 = W + H` with `W = sum_(i<j) c_ij` and the proof sentence "those objects [ordered pairs of a non-empty
  path and a non-empty Hamiltonian path of the complement] contribute `W`" (`c_ij` is THM-4167's exposure capacity,
  eq. (2) there: `c_e(T) = F_T({i}) + F_T({j}) - F_T(empty) - F_T({i,j}) = Q(i,j) + Q(j,i)`, i.e. `mu_R(i,j)`);
  as Session 210 `session210/graded-spectrum-gates/PROOF.md` eq. (8) `P_T(2) = (n+1) H(T) + m_(n-2)`; and it is
  `W(O)` at `O = V` in the same-day THM-5126 eq. (1').
- **Half the total.** THM-4184 Theorem 1 (2026-08-26; eq. (5) and the last line of its proof, `P^0 = W/2`), i.e.
  `sum_{i<j} mu_R(i,j) / 2 = sum over the non-empty EVEN subsets S of H(R[S]) H(R[V \ S])` (`H(empty) = 1`); the same
  balance is Session 210 Lemma 3 (`F_T(x) F_T(-x) = 1`) and is used in THM-4892 (F) 3.
- **(4.4).** THM-4177 Theorem 3, eq. (16) (2026-08-26) is the exact integer formula
  `mu_R(i,j) = 2 sum_P H(R - V(P))`, the sum over the directed paths `P` between `i` and `j` (either direction) with an
  odd number of arcs, i.e. an even number of vertices. (4.4) is its reduction mod 2 (by Redei), and its one-arc layer,
  eq. (18) `c^[1]_ij = 2 H(T - {i,j})`, is the bound `mu_R(i,j) >= 2 H(R - {i,j})` of Theorem C (1), two months before
  the Session 136 note credited there.

What is **new** is the congruence `= 2 (mod 4)` and what follows from it (`rank(beta_R) >= 2`, `abs(b_R) <= 1/2`,
Corollary 4.6). No earlier statement of the congruence, of `rank(beta_R) >= 2`, of "the number of even-order paths is
odd", or of `1/4 <= P(H = 1 mod 4) <= 3/4` was found (searches in section 6). It is one line from any of the above:
by THM-4184 (5) the half-total is a sum of `2^(m-1) - 1` odd numbers. The README's proof (pairing `S` with its
complement) is the same line and is correct.

*Checked.* `mu_R(i,j) = 2 sum_P H(R - V(P))` (THM-4177 (16)) on every labelled tournament of orders 2-7 and every class
of orders 2-9, with `mu` taken from second differences of directly counted `H(R + v_X)` (`x1_cube.c`, flag `p`,
counter `fail_4177`: 0); `sum mu = (m-1) H + F_1` and `sum mu / 2 = sum over non-empty even S` on all labelled
tournaments of orders 2-5 by permutation counting (`x3_bruteforce_small.py`, counters `thm4111`, `thm4184`: 0).

*Replacement text.*

- README section 4.1b, first line of the lemma: "**Lemma 4.5 (the two equalities are KNOWN; the congruence and its
  consequences are PROVED here).**" and, after the proof: "*Credits.* The first equality is `F_T(V) = H` (THM-4118
  eqs. (2), (4)). The second is THM-4111 eq. (6), `sum_(a != b) Q_T(a,b) = (n-1) H(T) + F_1(T)` (`F_1` = orderings
  with exactly one backward adjacency), restated as THM-4115 eq. (7), THM-4181 eq. (6) (`W = sum_(i<j) c_ij` counts
  the exposed words; `c_ij = mu_R(i,j)` is THM-4167's exposure capacity) and Session 210 eq. (8)
  (`P_T(2) = (n+1) H(T) + m_(n-2)`). The congruence is not stated in those files; it also follows in one line from
  THM-4184 eq. (5) (`W/2 = P^0`, the sum of `H(R[S]) H(R[V \ S])` over the non-empty even `S`) or from the odd-subset
  gate formula THM-4892 (F) 3. Read in THM-4118's terms it says that the response-lattice gcd `d_T` is `2 mod 4` (eq.
  (9) there has only `d_T` even), and in THM-4115's terms that `F_1(T) = 2 - (n-1) H(T) (mod 4)` (there: `F_1` even for
  odd `n`)."
- README section 4, end of the paragraph "Explicit polar form": "(4.4) is the reduction mod 2 of THM-4177 Theorem 3,
  eq. (16): `mu_R(i,j) = 2 sum_P H(R - V(P))` over the directed paths `P` between `i` and `j`, in either direction, on
  an even number of vertices; its one-arc layer, eq. (18), is the bound `mu_R(i,j) >= 2 H(R - {i,j})` of Theorem C (1)."
- THM-5127 status, "KNOWN, credited": add "- the two equalities of the parity law (4.5): `F_T(V) = H` (THM-4118 (2))
  and `sum_(a != b) Q_T(a,b) = (n-1)H + F_1` = number of ordered two-path covers (THM-4111 (6), THM-4115 (7),
  THM-4181 (6), Session 210 (8)); half the total as a sum over even subsets (THM-4184 (5)); the exact integer form of
  (4.4), `mu_R(i,j) = 2 sum over even-order paths P between i and j of H(R - V(P))` (THM-4177 (16); (18) is
  `mu >= 2 H(R - {i,j})`)". In "NEW here" replace "the parity law (4.5) with its corollary" by "the congruence
  `= 2 mod 4` in (4.5) with its corollaries (`rank(beta_R) >= 2`, `abs(b_R) <= 1/2`,
  `1/4 <= P(H = 1 mod 4) <= 3/4`)". In (C) (6) insert after "(Parity law; added after the audit": "; the two
  equalities are KNOWN, THM-4118 (2) and THM-4111 (6)". Add THM-4111, THM-4115, THM-4167, THM-4177, THM-4181,
  THM-4184, THM-4892 to `related`.

---

## 2. Minor findings

**m1. "independent" (overclaim in wording; the number is right).** *Where.* README 4.1a: "measures how far the mod-4
classes of two *different* one-vertex extensions of the same base are from independent"; P3-P6 outcome paragraph: "So
the mod-4 classes of two different one-vertex extensions of a random 10-vertex base are independent to within 4% of
the twin term". *What is wrong.* The quantity is one signed total: `rho_m - 1` is the sum over all ordered pairs
`X != X'` of `E_R (-1)^(q_R(X) + q_R(X'))`, in units of the twin term; per pair it is `(rho_m - 1)/(2^m - 1)`
(`3.97 * 10^(-5)` at `m = 10`). It bounds no individual pair, and individual pairs are far from uncorrelated: the
pair `(empty, V)` has correlation 1 (both extensions have `H = H(R)`), and the `2^m` antipodal pairs `(X, V \ X)`
have, on average over `X`, correlation exactly `P_R(all beta_R-degrees even)` (`0.002048` at `m = 10`, 52 times the
overall average; exact identity `q_R(X) + q_R(V \ X) = beta_R(X, 1)`, which follows from
`q_R(V) = 0`, so that `E_X (-1)^(q_R(X) + q_R(V \ X))` is 1 if all `beta_R`-degrees are even and 0 otherwise; checked
on every tournament of every run, counter `fail_anti`: 0). Nor are the two deletions independent in the literal
sense: `E[s s'] = 1.04 * 2^(-10)` against `D_11^2 = 2.2 * 10^(-7)`. *Replacement.* (i) "... measures the total
correlation between the mod-4 classes of two *different* one-vertex extensions of the same base, in units of the twin
term (an average over the pair of
out-sets: `(rho_m - 1)/(2^m - 1)` per ordered pair)". (ii) "So at order 10 the pairs of different out-sets add
4.06% to the twin term in total: for independent uniform out-sets `X != X'` of a random 10-vertex base the mod-4
classes of the two extensions have correlation `3.97 * 10^(-5)` on average. Individual pairs are not uncorrelated:
`(empty, V)` has correlation 1, and the antipodal pairs `(X, V \ X)` have correlation
`P_R(all beta_R-degrees even) = 0.00205` on average (5% of the excess). The exact labelled mod-4 discrepancy at
order 11 is `260546929/2^39`."

**m2. Proposition 4.3, chain: the range of `k` is missing, and one clause of the proof.** The definition has
"(`1 <= k <= n`)" and the inequality `B(n,k)^2 <= B(n,k+1)` no range. The step needs an unmarked vertex, i.e.
`n - k >= 1`, and `B(n, n+1)` does not exist. *Replacement.* "Then `B(n,k)^2 <= B(n,k+1)` for `0 <= k <= n-1`, hence
`abs(D_n)^(2^k) <= B(n,k)` for `0 <= k <= n`." (THM-5127 (C) (5): same.) In the proof, after "two independent stars
`u, u'`" add: "(the arc between `u` and `u'` lies on no counted cycle; orient it uniformly at random, so that the
tournament on the `n + k + 1` vertices is uniform)". With these the proof is complete. Remark for the text:
`B(n,k)` is not monotone in `k` and the end of the chain is trivial (`B(n,n) = 1` for even `n`, since a cycle meeting
every pair exactly once has `n` vertices): `B(4, .) = 1/2, 1/4, 1/8, 1/8, 1`.

**m3. The calibration "for a uniformly random quadratic form the left side is exactly `2^(-m)`" holds for one
ensemble and not for the null model of section 4.1b.** It is true when the linear part `ell` is uniform on `F_2^m`
and independent of the polar form (for instance all `2^(m + C(m,2))` forms without constant term equally likely):
then `E (-1)^(q(X) + q(X')) = 0` for `X != X'`. Tournament forms are never of that kind: `q_R(V) = 0`, `ell` has odd
weight and `beta` an odd number of edges (Lemma 4.5). Under the uniform law on that constrained set the value is
`2^(-m) (1 + 2^(-(m-1)))` for `m >= 3` (brute force over all forms, `m <= 5`: `x8_nullmodel.out`), because of the
antipodal pairs of m1. *Replacement* (README 4.1a; THM-5127 (C) (5) and HYP-9659 (B'') analogously): "... is exactly
`2^(-m)` when the linear part is uniform on `F_2^m` and independent of the polar form, because the values at two
distinct points are then uncorrelated. (Under the constraints that every `q_R` satisfies - `q_R(V) = 0`, linear part
of odd weight, odd number of edges, section 4.1b - the uniform law gives `2^(-m)(1 + 2^(-(m-1)))`, `m >= 3`: the
antipodal pair `X, V \ X` has `q(X) + q(V \ X) = beta(X, 1)`. For tournaments this antipodal term is
`2^(-m) P_R(all beta_R-degrees even)`, with `P_R = 1, 1/4, 1/16, 0.0088, 0.0242, 0.0075, 0.0047, 0.0020` at
`m = 3, ..., 10`.)"

**m4. "the statement 'the rank of `beta_R` is large' can be written without forms".** Only one direction holds.
*Replacement.* "and "the rank of `beta_R` is large in probability" implies a statement without forms: the mod-4
classes of two one-vertex deletions of a random tournament decorrelate (not conversely: forms of low rank and zero
bias do not contribute to `E_R[b_R^2]`)".

**m5. "(4.5) is the only linear relation on `beta_R` found here"; HYP-9659 "no second linear relation on `beta_R` is
known".** (4.5) is an affine relation, and with late item L2 the sentence can be a statement. *Replacement.* "(4.5) is
the only affine relation over `F_2` among the entries of `beta_R` that holds for all tournaments of a given order
`4 <= m <= 11` (FINITE-EXACT certificate, `affine_span.py`, audit 2 `x6_affine.c`: the `beta_R` affinely span the
hyperplane `sum = 1`, and the pairs `(ell_R, beta_R)` the flat `sum ell = sum beta = 1`). This does not make `beta_R`
generic: `R -> beta_R` reaches only 22 of the 32, 352 of the 512 and 5,578 of the 16,384 graphs with an odd number of
edges at `m = 4, 5, 6` (`audit2/x10_image.out`)."

**m6. "differs ... by less than `3 * 10^(-5)` in every rank class for `m >= 8`" is stated for all `m >= 8` but
computed for `m <= 11`.** It is true for all `m >= 8`, with an exact formula. Let `S(m, r)` be the sum of
`(-1)^(number of edges)` over the alternating `m x m` matrices of rank `r`. Splitting off the last vertex
(`A = [[A', u], [u^T, 0]]`: the rank is unchanged iff `u` lies in the image of `A'`; the sum of `(-1)^|u|` over the
image is `2^rank` if `A' 1 = 0` and `0` otherwise; the matrices with `A' 1 = 0` are the even graphs, in rank- and
edge-parity-preserving bijection with all graphs on one vertex fewer) gives
`S(m, r) = 2^r S(m-2, r) - 2^(r-2) S(m-2, r-2)`, hence `S(m, 2s) = (-1)^s 4^C(s,2) [n choose s]_4` with
`n = floor(m/2)` (Gaussian binomial at `q = 4`), and
`P(rank = 2s | odd) - P(rank = 2s) = -S(m, 2s) / 2^C(m,2)`, of absolute value at most `1.46 * 2^(-n^2)`. Checked against
brute force over all graphs for `m <= 8` (`x5_altrank.c`, `x5_oddedge_check.out`). The largest difference is
`85/2^22 = 2.03 * 10^(-5)` at `m = 8` (rank 6), `7.9 * 10^(-8)` at `m = 9`; at `m = 7` it is `4.0 * 10^(-5)`, so "8" is
sharp for the stated threshold. *Replacement.* "... differs from the unconditioned law by exactly
`4^C(s,2) [n choose s]_4 / 2^C(m,2)` in the rank class `2s` (`n = floor(m/2)`; audit 2), at most `2.03 * 10^(-5)` for
`m >= 8`, so the comparison table of section 4.1 is unaffected."

**m7. OBSERVATION (rank law equal to the odd-edge law at `m = 4, 5`).** Confirmed exactly (`5/8, 3/8`;
`5/32, 27/32`; at `m = 6`: `549/1024, 475/1024` against `21/1024, 567/1024, 436/1024`). Add: "(trivially also at
`m = 2, 3`, where every rank is 2). It is a coincidence of the two rank laws only: the law of `beta_R` itself is not
the odd-edge law (22 of 32 graphs reached at `m = 4`, with fibres of sizes 2 and 4; no graph with one edge occurs;
`audit2/x10_image.out`)."

**m8. P3, P4 and P6 are one prediction, not three, and the labelled rank law of order 10 was in-sample.** With the
labelled rank fractions `p_r` of order 10 - computed by the first audit before the pre-registration
(`audit/b_rankonly_10_parts.out`, in the same commit `23d6c4033f`) - and `f_c` the non-zero-bias frequency at corank
`c`: `rho_10 = p_10 + 4 p_8 f_2 + 16 p_6 f_4 + ... = 0.41894 + 2.23766 f_2 + 0.34537 f_4 + 0.00365 f_6 + ...` and
`E_R abs(b_R) = 0.013092 + 0.034963 f_2 + 0.002698 f_4 + ...`. The P4 interval `f_2 in [0.255, 0.290]` maps (other
terms at their observed values) to `rho_10 in [1.020, 1.099]` and `E_R abs(b_R) in [0.02225, 0.02347]`, inside the P3
and P6 intervals; with `f_4` anywhere between the uniform-form model `1/16` and the order-8 value `0.2` the images are
`[1.012, 1.137]` and `[0.02218, 0.02377]`, essentially the registered intervals `[1.01, 1.12]` and
`[0.0220, 0.0238]`. So the run confirms two independent predictions, P4 and P5. PREREGISTRATION.md says "In-sample
data: orders 2-9" only. *Replacement* (README 4.1a, after
"All four CONFIRMED (...)"): "(P3 and P6 are not independent of P4: given the labelled rank fractions of order 10,
known from the first audit before this pre-registration, `rho_10` and `E_R abs(b_R)` are affine functions of the
non-zero-bias frequencies, dominated by the corank-2 term; the independent confirmations are P4 and P5.)" Add the
same sentence as a dated annotation to PREREGISTRATION.md (annotate, do not rewrite).

**m9. HYP-9659 (B''), "Relations (PROVED): ... neither of (B') and bounded `rho_m` implies the other".** Both are
open statements about one sequence, so a non-implication between them cannot be typed PROVED (if both hold, each
implies the other). What the parenthesis shows is that the inequalities at hand give neither. The two implications
before it are correct: (B) gives `E_R[b_R^2] <= E_R[2^(-rank)] -> 0`; bounded `rho_m` gives
`D_(m+1)^2 <= E_R[b_R^2] <= C 2^(-m)`. *Replacement.* "Relations (PROVED): (B) implies the qualitative form; (B')
implies it with rate `O(2^(-m/2))`; bounded `rho_m` implies a mod-4 discrepancy `O(2^(-m/2))` at order `m+1`
(Cauchy-Schwarz) and (B') restricted to the forms of non-zero bias (`E_R abs(b_R) <= sqrt(rho_m) 2^(-m/2)`). No
implication between (B') and bounded `rho_m` is known: forms of zero bias count in (B') only, and (B') bounds
`E_R[b_R^2]` only by `C 2^(-m/2)`."

**m10. HYP-9659 (B''), "This is the smallest open statement on the mod-4 route."** `D_n -> 0` is weaker (it follows
from the qualitative form of (B'') by `D_(m+1)^2 <= E_R[b_R^2]`) and is also open. *Replacement.* "The weakest open
statement on the mod-4 route is `D_n -> 0`; the qualitative form of (B'') implies it."

**m11. Three imprecise phrases.** (i) 4.1b "By (4.4) the lemma says that the number of directed paths ... is odd":
(4.4) translates the mod-2 shadow of the lemma. Replace "the lemma says" by "the statement that `beta_R` has an odd
number of entries says". (ii) 4.1a "forms of positive corank whose linear part vanishes on the radical": on the
radical it is `q_R` that is linear; replace by "forms of positive corank for which `q_R` vanishes on the radical".
(iii) Table header "`E_R abs(b_R)` (bound of (4))" and PREREGISTRATION P6 "the bound of Theorem C (4)": Theorem C (4)
as stated has `E_R[2^(-rank/2)]` on the right (`0.05077` at `m = 10`); `E_R abs(b_R)` is the sharper bound inside its
proof. Replace by "(the bound in the proof of (4))".

**m12. Rounding.** README section 4.1 ("Exact values from the audit") has `63103/8388608 = 0.007523`; the 4.1a table
has `0.007522`. The fraction is `0.00752246`, so `0.007522` is right (the first audit's report has the same slip).

**m13. Section 9, new lines.** All of them run (section 7). (i) The shipped `out/biasweights_summary_2_9.out`,
`out/biasweights_summary_10.out` and `out/box_norm_check.out` are the standard output of `merge_bias.py` /
`box_norm_check.py`, but the lines do not redirect; add `> out/...`. (ii) `deletion_pair_corr.py 7` is announced as
20 minutes of pure Python; the two entries it adds over argument `6` (route 1 at `m = 5`, route 2 at `m = 6`) can be
had in seconds from compiled code (here: `x2_box 7 labelled`, `x1_cube 6 labelled`). (iii) No lines yet for the four
late-item scripts.

**m14. Status placeholders.** HYP-9659 "(B'') Deletion-pair form (added after the audit; not audited)" -> "(added
after the first audit; audited, `audit2/AUDIT2.md`)"; README header "second audit (of these additions): pending at
this commit" and THM-5127 "SECOND AUDIT: PENDING at this commit" -> the verdict of this report.

**m15. Connections the packet missed (not errors).** (i) THM-4118 Theorem 1: the response lattice of a tournament is
`d_T Z` with `d_T` "a positive even integer" (eq. (9)); Lemma 4.5 gives `d_T = 2 (mod 4)` for every tournament of
order `>= 2`, and Corollary 4.6 (a) is the quantitative form of "`F_T(S) = H mod d_T`, and not mod 4". (ii) THM-4115
(13a): `F_1 = (n-1) C_1 + 2 C_2` and "`F_1` even for odd `n`"; Lemma 4.5 is `F_1 = 2 - (n-1) H (mod 4)`. (iii) The
gate (late item L1). (iv) THM-478 (4) records that on the arc cube the digits of `H` leave the low-degree zone at once
(digit 2 of degree 3, 6, 7, 11 at `n = 4..7`); on a vertex cube the ladder is the classical `2^k` (late item L3).
(v) (4.3) holds at every frequency, with the same proof: `E_R abs(E_X e(theta H(R + v_X)))^2 = E_T e(theta (H(T - v') -
H(T - v)))`; this is the form in which the third Principle of section 7 states it, and the second moment that Lemma
C4 bounds. (vi) For general `k` the README defines `B(n,k)` through odd cycles only; it has the same "deletion" form
as the middle term of (4.3): `(-1)^(N_k) = prod_D s(W - D)` over the `2^k` transversals `D` of the marked pairs,
`s = +-1` for `H = 1, 3 mod 4`, which is what makes it a box norm (used and checked here as route 2).

---

## 3. Proofs: what was checked and how

| item | result | how |
|---|---|---|
| Proposition 4.3, first equality | correct | by hand: `T - v' = R + v_X`, `T - v = R + v'_X'`, `X, X'` independent and uniform given `R`; `(a-c)/2 + (b-c)/2 - (a-b)/2 = b - c` is even. Executed: both sides by direct enumeration, all labelled tournaments: `m <= 5` (`x2_box N labelled`, `N <= 7`, against `x1_cube m labelled`), DP-free for `m <= 4` (`x3`); from class representatives with weights for `m = 6, 7, 8` (classes of orders 8, 9, 10 against vertex cubes of orders 6, 7, 8) |
| Proposition 4.3, second equality (`N_1`) | correct | cycles avoiding both vertices cancel, cycles through both lie in neither deletion. Executed pointwise: `(-1)^((H(T-v)-H(T-v'))/2) = (-1)^(N_1(T))` with `N_1` by explicit enumeration of odd cycles on every labelled `T`, `n <= 6` (`x3`), and by a Hamiltonian-cycle DP on every labelled `T` with `n <= 7`, on all pairs `v, v'` of every class of order 8 and on a random sample of pairs at orders 9 and 10 (`x2_box`): 0 mismatches |
| `D_(m+1)^2 <= (E abs(b_R))^2 <= E[b_R^2]` | correct | proof of Theorem C (4), then Cauchy-Schwarz; exact values `m = 2..10` (`x1_summary.out`); equality throughout at `m = 2, 3`, first inequality an equality at `m = 4, 5` |
| chain `B(n,k)^2 <= B(n,k+1)` | correct for `0 <= k <= n-1`; range missing (m2) | general step redone by hand: condition on `W - u`, Cauchy-Schwarz, second star `u'`, arc `u u'` free; the enlarged tournament is uniform and has `k+1` disjoint pairs; the exponent is `N_(k+1)`. Executed: every `B(n,k)` with `n + k <= 10`, two routes (below) |
| `abs(D_n)^(2^k) <= B(n,k)` | correct for `0 <= k <= n` | induction; `B(n,k) >= 0` for `k >= 1` because it is a mean of squares |
| calibration "twins contribute `2^(-m)`" | correct | `P(X = X') = 2^(-m)` |
| calibration "uniformly random quadratic form: exactly `2^(-m)`" | correct for `ell` uniform and independent of `beta`; not for the constrained null model (m3) | brute force over all forms, `m <= 5` (`x8_nullmodel.out`): `1/4, 1/8, 1/16, 1/32`; constrained: `1/4, 5/32, 9/128, 17/512` |
| Lemma 4.5, `sum a_i = sum mu` | correct (KNOWN, M1) | `H(R + v_V) = H(R)`; executed with `a`, `mu` from first / second differences of directly counted `H(R + v_X)` |
| Lemma 4.5, `sum a_i` = ordered two-path covers `= sum_S H(R[S]) H(R[V \ S])` | correct (KNOWN, M1) | (4.1) at `X = {i}`; executed with the split sum from directly counted `H` of induced subtournaments |
| Lemma 4.5, congruence, `m = 2` | correct | `2^(m-1) - 1` unordered pairs `{S, V \ S}`, an odd number for `m >= 2`; `m = 2`: `mu = 2`, `a = (2, 0)`, split `= 1 + 1`. (At `m = 1` all three sums are 0.) |
| odd number of `beta`-edges, odd weight of the linear part, `rank >= 2`, `abs(b_R) <= 1/2` | correct | counters `fail_edges`, `fail_ell`, rank tables, `fail_c3` |
| "the number of directed paths on an even number of vertices is odd" | correct | `sum_{i<j} beta_ij = sum_(i != j) P_ij`; directly: `2^(m-1) - 1` non-empty even subsets, Redei on each. Executed by enumeration of all vertex sequences, labelled `m <= 5` (`x3`), and by DP on all classes `<= 10` (`fail_evenpaths`) |
| Corollary 4.6 (a), (b), (c); "equality at `n = 3, 4`" | correct | (a) from `abs(b_R) <= 1/2`; (b) averages (a) over `R = T - v`, needs `n - 1 >= 2`; `D_3 = D_4 = 1/2` (`P(H = 1 mod 4) = 3/4`), `D_5 = 13/32`; (c): equality at `m = 2, 3` |
| HYP-9659 (B''): (B) => qualitative form; bounded `rho_m` => `D = O(2^(-m/2))` | correct | by hand (m9) |
| HYP-9659 (B''): "neither of (B') and bounded `rho_m` implies the other", typed PROVED | not a provable statement as typed (m9) | |
| L1 (gate corollary) | correct | section 8 |
| L3 (digit-degree law), incl. the integrality / degree step | correct; sharpened | section 8 |

Scale of the executed checks (all counters zero; `x1_summary.out`). Direct mode means: `H(R + v_X)` for every
out-set `X` by a fresh subset DP on the `(m+1)`-vertex tournament; `a_i`, `mu_ij`, `beta`, `ell`, `q` are then read
off that table and nothing uses the cut formula.

| universe | mode | tournaments | checks |
|---|---|---|---|
| all labelled, `m = 2..5` | permutations only (`x3`) | 1,098 | (4.2), (4.5), THM-4111 (6), THM-4184 (5), even-path parity, (4.4), Corollary 4.6 (a) |
| all labelled, `m = 2..7` | direct, no automorphism groups | 2,131,018 | (4.2) on every `X`, parity and `mu >= 2 H(R - {i,j})`, Theorem C (3), (4.5) (three counters), edge and weight parity, (4.4), THM-4177 (16), antipodal identity, gate formula, digit degrees |
| all classes, `m = 2..8` | direct, weights `m!/abs(Aut)` | 7,411 | the same |
| all classes, `m = 9` | direct | 191,536 | the same (98,066,432 ten-vertex tournaments counted) |
| all classes, `m = 10` | insertion formula (4.1) with exact integers + one random `X` per class recounted directly | 9,733,056 | Theorem C (3), (4.5), parities, antipodal identity, digit degrees; 9,733,056 direct controls |

Other audited sentences.

- *Paragraph "Explicit polar form (4.4)".* The text is what the first audit verified as m24
  (`q_R(X) = sum_{i in X, j notin X} P_ij`, `beta_R(i,j) = P_ij + P_ji`, linear part `q_R({i}) = sum_{j != i} P_ij`,
  all mod 2). Re-verified here on a larger universe: every labelled tournament of orders 2-7 and every class of orders
  2-9, polar form and linear part (counter `fail_44`: 0). Missing credit: M1.
- *Header and Outcome items 3-4.* The statements and numbers agree with sections 4.1a / 4.1b ("1.041 times the twin
  term `2^(-10)`", "`260546929/2^39 = 0.000474`", "orders 8, 9, 10: `0.0415, 0.00752, 0.00242`": confirmed). Outcome
  item 3 lists the parity law under "(PROVED)" without the credits of M1.
- *Section 7, "One picture".* It has no sentence about 4.1a / 4.1b; the replacement R2 of the first audit is applied
  faithfully, with two added precisions ("of an `h >= 3`", "bases of order `<= 10`") that agree with the first
  audit's data (its M5: the only `h` whose window is filled by one base is `h = 1`).
- *Section 7, third Principle.* "the second moment of the bias is the correlation of two one-vertex deletions (4.3)":
  correct, and true at every frequency (m15 (v)). "two new vertices give a quartic": correct (a Hamiltonian path uses
  at most two arcs at each new vertex).

---

## 4. Recomputation table (claimed / recomputed)

Exact labelled statistics of the vertex cubes (`x1_summary.out`; `m <= 9` from direct counts, `m <= 7` also without
automorphism groups; `E_R[b_R^2]` for `m <= 8` also from the deletion-pair side, `x2_summary.out`).

| `m` | `E_R[b_R^2]` claimed | recomputed | `rho_m` | `E_R abs(b_R)` claimed / recomputed | `D_(m+1)` claimed | recomputed |
|---|---|---|---|---|---|---|
| 2 | 1/4 | 1/4 | 1.000 / 1.000000 | 0.5 / 1/2 | 1/2 | 1/2 |
| 3 | 1/4 | 1/4 | 2.000 / 2.000000 | 0.5 / 1/2 | 1/2 | 1/2 |
| 4 | 23/128 | 23/128 | 2.875 / 2.875000 | 0.40625 / 13/32 | 13/32 | 13/32 |
| 5 | 21/256 | 21/256 | 2.625 / 2.625000 | 0.25 / 1/4 | 1/4 | 1/4 |
| 6 | 1621/65536 | 1621/65536 | 1.583 / 1.583008 | 0.12793 / 131/1024 = 0.127930 | 853/8192 | 853/8192 |
| 7 | 57619/4194304 | 57619/4194304 | 1.758 / 1.758392 | 0.08835 / 46321/524288 = 0.088350 | 1361/32768 | 1361/32768 |
| 8 | 1304557/268435456 | 1304557/268435456 | 1.244 / 1.244123 | 0.05028 / 52721/1048576 = 0.050279 | 63103/8388608 | 63103/8388608 = 0.007522464 |
| 9 | 155129173/68719476736 | 155129173/68719476736 | 1.156 / 1.155802 | 0.03155 / 67761877/2147483648 = 0.031554 | 650471/268435456 | 650471/268435456 = 0.002423193 |
| 10 | 71511414601/70368744177664 | 71511414601/70368744177664 | 1.041 (1.040628) / 1.040628 | 0.02256 (0.022563) / 6202089377/274877906944 = 0.022563 | 260546929/549755813888 | 260546929/549755813888 = 0.000473932 |

| quantity | claimed | recomputed |
|---|---|---|
| order 10: classes; labelled weight | 9,733,056; `2^45` | 9,733,056; 35184372088832 `= 2^45` |
| order 10: labelled rank fractions (rank 10 / 8 / 6 / 4) | 0.4189420 / 0.5594148 / 0.0215859 / 0.0000570 | 28789471611/68719476736 = 0.4189420 / 76885390455/137438953472 = 0.5594148 / 1483371925/68719476736 = 0.0215859 / 489825/8589934592 = 0.0000570 (rank 2: 38745/137438953472) |
| order 10: classes zero-bias + non-zero per rank 2, 4, 6, 8, 10 | 5+2, 523+145, 191763+18799, 4005054+1439603, 0+4077162 | identical |
| order 10: non-zero-bias frequency at rank 8 (P4) | 0.26407 | 0.26407 (rank 6: 0.08753; rank 4: 0.14972; rank 2: 0.14634) |
| P3 / P4 / P5 / P6 observed | 1.0406 / 0.2641 / 0.000474 / 0.02256 | 1.040628 / 0.26407 / 0.000473932 / 0.022563; all inside the registered intervals |
| control at order 10 (Theorem C (1)-(3) on every class) | no failure | no failure (and 9,733,056 direct recounts of a random `H(R + v_X)`: 0 mismatches) |
| non-zero-bias frequency, corank 1 (`m = 3, 5, 7, 9`) | 1, 0.815, 0.636, 0.527 | 1, 0.81481, 0.63612, 0.52705 |
| corank 2 (`m = 4, 6, 8, 10`) | 1, 0.522, 0.326, 0.264 | 1, 0.52186, 0.32637, 0.26407 |
| corank 3 (`m = 5, 7, 9`) | 1, 0.539, 0.195 | 1, 0.53908, 0.19515 |
| corank 4 (`m = 8, 10`) | 0.200, 0.088 | 0.20000, 0.08753 |
| (not listed in the README) corank 5 (`m = 9`), 6 (`m = 8, 10`), 7 (`m = 9`), 8 (`m = 10`) | - | 0.22546; 1, 0.14972; 1; 0.14634 |
| (4.5) on all classes of orders 2-9 and of order 10 | no failure (198,947 + 9,733,056 classes) | no failure (same counts); histograms of the number of `beta`-edges identical to `out/oddedge_*.out` |
| `B(3,1), B(3,2), B(3,3), B(4,1), B(4,2), B(5,1)` | 1/4, 1/8, 1/8, 1/4, 1/8, 23/128 | 1/4, 1/8, 1/8, 1/4, 1/8, 23/128 (two routes each, see below) |
| odd-edge law against unconditioned law, `m >= 8` | difference `< 3 * 10^(-5)` in every rank class | exact maximum `85/4194304 = 2.03 * 10^(-5)` (`m = 8`, rank 6); proved for all `m >= 8` (m6) |
| labelled rank law at `m = 4`, `m = 5` | `5/8, 3/8`; `5/32, 27/32`, equal to the odd-edge law | `{2: 5/8, 4: 3/8}`; `{2: 5/32, 4: 27/32}`; equal |
| `m = 6` | not equal, no rank 2 | `{4: 549/1024, 6: 475/1024}` against `{2: 21/1024, 4: 567/1024, 6: 109/256}` |
| THM-5127: `rho_m`, `m = 2..9` | 1, 2, 2.875, 2.625, 1.583, 1.758, 1.244, 1.156 | 1, 2, 2.875, 2.625, 1.583008, 1.758392, 1.244123, 1.155802 |
| packet raw totals, orders 2-10 (`out/biasweights_2_9.out`, `out/biasweights_10_parts.out`, `out/oddedge_2_9.out`, `out/oddedge_10_parts.out`) | - | identical field by field: classes, weight, signed sum `S`, plain sum `U`, absolute sum `A`, classes and labelled weight per rank and per zero / non-zero bias, edge histogram (`x9_compare_author.out`) |

**Box norms** (`x2_summary.out`). Route 1: `N_k` by explicit counting of the directed Hamiltonian cycles of every odd
vertex set meeting each marked pair once. Route 2: `(-1)^(N_k) = prod over the 2^k transversals D (one vertex of each
pair) of s(W - D)`, `s = +-1` for `H = 1, 3 mod 4` (a cycle is missed by `2^(#pairs it avoids)` transversals, or by
none if it contains a pair). The two routes agree pointwise on all labelled tournaments with `N <= 7`, on all
5,366,730 (class, placement) pairs of orders 2-8, on a random 1/64 of those of order 9 (7,842,119 placements) and on
4,473,488 placements of a slice of 30,144 classes of order 10; the table is route 2 over all placements, classes
weighted by `N!/abs(Aut)`.

| `n` | `B(n,0) = D_n` | `B(n,1)` | `B(n,2)` | `B(n,3)` | `B(n,4)` | `B(n,5)` |
|---|---|---|---|---|---|---|
| 3 | 1/2 | 1/4 | 1/8 | 1/8 | | |
| 4 | 1/2 | 1/4 | 1/8 | 1/8 | 1 | |
| 5 | 13/32 | 23/128 | 85/1024 | 221/4096 | 6971/131072 | 6971/131072 |
| 6 | 1/4 | 21/256 | 161/4096 | 1011/32768 | 6971/131072 | |
| 7 | 853/8192 | 1621/65536 | 19581/1048576 | 237873/16777216 | | |
| 8 | 1361/32768 | 57619/4194304 | 156235/16777216 | | | |
| 9 | 63103/8388608 | 1304557/268435456 | | | | |
| 10 | 650471/268435456 | | | | | |

`B(n,1) = E_R[b_R^2]` at order `m = n - 1`: the column reproduces the table above from the deletion-pair side
(classes of order `n + 1`). Every inequality `B(n,k)^2 <= B(n,k+1)` and `abs(D_n)^(2^k) <= B(n,k)` available in the
table holds; equality `D_n^2 = B(n,1)` at `n = 2, 3, 4`. The repeated value `6971/131072` is not an artefact:
`B(6,4) = B(5,5)` identically (with four pairs and two unmarked vertices a counted cycle has five vertices, so the two
unmarked vertices act as a fifth pair), and `B(n, n-1) = B(n, n)` at `n = 3, 5` is an equality of values (at `n = 3`:
the bias of the star of the unmarked vertex is `0` or `+-1/2` and its sign has mean `1/2`). The packet's range
`n + k <= 6` is the first four rows up to the anti-diagonal.

**Not in the packet (new exact data).** Labelled residues of `H` mod 8 (classes 1, 3, 5, 7) at order 11:
`0.25013, 0.24987, 0.25011, 0.24989`. `P_R(all beta_R-degrees even)` at `m = 3..10`:
`1, 1/4, 1/16, 9/1024, 99/4096, 3915/524288, 637373/134217728, 35183785/17179869184`. Decomposition of the excess at
order 10: `rho_10 - 1 = 0.040628 = 0.031485 (corank 2) + 0.008643 (corank 4) + 0.000489 (corank 6) + 0.000010 (corank
8)`.

---

## 5. Pre-registration P3-P6

- `git show 23d6c4033f:.../PREREGISTRATION.md` contains the complete second block (P3 `[1.01, 1.12]`, P4
  `[0.255, 0.290]`, P5 positive and `[0.0003, 0.0011]`, P6 `[0.0220, 0.0238]`, control, refutation clause). The
  README, HYP-9659 and the outcome paragraph quote the same intervals. The working tree / later commits add only the
  outcome paragraph. At that commit the README table row for `m = 10` reads "pre-registered (P3-P6); run pending at
  this commit", and no order-10 bias output is in the tree.
- Timestamps (read-only `git log`, `git reflog`, file mtimes): commit `23d6c4033f` 13:11:35 -0700; `origin/main`
  "update by push" to `06824d5ce5` (a descendant) 13:13:16; `out/biasweights_10_parts.out` and
  `out/biasweights_summary_10.out` written 13:17:53. The packet's program needs about 16 CPU-minutes for the 9,733,056
  classes (measured on a 1/128 slice: 5.8 s), i.e. 2 to 4.5 minutes on 8 shared cores. So the files are consistent
  with a run started after the push, and in any case the commit precedes the earliest possible start. (Timestamps
  cannot exclude an earlier unrecorded run; nothing suggests one.)
- All four outcomes are confirmed by the independent recomputation (section 4).
- Limits of the evidence: m8 (P3 and P6 follow from P4 and data that were in-sample). The P5 interval spans a factor
  3.7 (the in-sample ratios `D_n / D_(n+1)` are 2.5, 5.5, 3.1); its sign prediction is a genuine one.

---

## 6. Prior art in the repository

Searches: `grep -rlE --include='*.md'` over `01-canon/theorems`, `05-knowledge/hypotheses`, `05-knowledge/results`,
`session1*`, `session2*` for: two-path cover, path cover, exposed word, exposure capacity, capacity tensor, `P_T(2)`,
even number of vertices, even-order path, `H mod 4` and variants, mod-4 discrepancy, box norm, Gowers, uniformity
norm, deletion pair, two deletions, response lattice, `d_T`, gate with mod 4; then line-level reading of the hits in
the cut-quadratic lane (THM-4097, 4111, 4114, 4115, 4118, 4123, 4167, 4177, 4181, 4184, 5122, 5126, Session 136
`positive_cut_geometry.md`), the gate lane (Session 210 `graded-spectrum-gates/PROOF.md`, THM-4892, THM-4900) and the
`H mod 4` lane (THM-466, THM-478, THM-791, THM-865, THM-1980, THM-4748, THM-4930).

| statement | found? | where |
|---|---|---|
| `sum_i a_i = sum_{i<j} mu_R(i,j)` | yes (implicit) | THM-4118 eqs. (2), (4) |
| `sum_{i<j} mu_R(i,j)` = number of ordered two-path covers `= (n-1)H + F_1` | **yes** | THM-4111 eq. (6); THM-4115 eq. (7); THM-4181 eq. (6); Session 210 eq. (8); THM-5126 eq. (1') at `O = V` |
| half of it = sum over non-empty even `S` of `H(R[S]) H(R[V \ S])` | **yes** | THM-4184 eq. (5) (`P^0 = W/2`); Session 210 Lemma 3; THM-4892 (F) 3 |
| the residue `2 mod 4` of that total | no | one line from the previous row and Redei. THM-4757's `N_L` are signed sums over linear forests (weights `prod S_e`), not the two-path-cover count, and that file has no residue of it |
| `rank(beta_R) >= 2`; "odd number of even-order paths" | no | THM-4118 eq. (9) has only "`d_T` is a positive even integer" |
| `1/4 <= P(H = 1 mod 4) <= 3/4`, or any bound on the mod-4 law of `H` over random tournaments | no | THM-466, THM-478, THM-1980, THM-791 concern `H mod 4` of individual tournaments or structured families |
| `mu_R(i,j)/2 = P_ij + P_ji (mod 2)` (4.4) | **yes, in exact form** | THM-4177 eq. (16); Session 136 (`w_ij = R_ij + R_ji`), already noted by the first audit |
| Proposition 4.3 (deletion-pair identity), box-norm chain | no | "box norm" / "Gowers" occur only in the quasirandom-group and LRC packets |
| the gate preserves `H mod 4` (L1), residue of `P_T(2)` | no | Session 210, THM-4892, THM-4900 state the formulas and no residue mod 4 |
| digit-degree bounds for `H` | on the arc cube only | THM-478 (2), (4); its note cites the classical `2^i` law for the Hamming weight |

---

## 7. Reproduction lines of section 9

Run in a scratch copy of the output names (nothing in `out/` was overwritten):

- `cc -O2 -o biasweights biasweights.c`; orders 2-9: output **byte-identical** to `out/biasweights_2_9.out`;
  `merge_bias.py` on it: byte-identical to `out/biasweights_summary_2_9.out`; `merge_bias.py` on the shipped
  `out/biasweights_10_parts.out`: byte-identical to `out/biasweights_summary_10.out`.
- `cc -O2 -o oddedge oddedge.c`; orders 2-9: byte-identical to `out/oddedge_2_9.out`.
- `oddedge_model.py 11`: byte-identical to `out/oddedge_model.out`; `box_norm_check.py`: byte-identical to
  `out/box_norm_check.out`.
- `deletion_pair_corr.py 6` (4 s): identical to the shipped file except for the two entries that need argument 7
  (not rerun: 20 minutes of Python; both entries, `21/256` by route 1 at `m = 5` and `1621/65536` by route 2 at
  `m = 6`, are confirmed by the auditor's own programs).
- The two order-10 loops (`biasweights 10`, `oddedge 10`) were not rerun; the program was run on a 1/128 slice, and
  the shipped order-10 parts agree with the auditor's totals field by field.
- Late-item scripts: `gate_mod4.py`, `affine_span.py 3..11` and `deletion_rank.py 4..7` reproduce
  `out/gate_mod4.out`, `out/affine_span.out`, `out/deletion_rank.out` byte for byte; `digit_degree.py 3..6` reproduces
  the first four lines of `out/digit_degree.out` (the order-7 line was not rerun in Python; its content is confirmed by
  the auditor's direct counts).

---

## 8. Late items

### L1. Gate corollary: `H(Gamma R) = H(R) (mod 4)` for every tournament `R` on `m >= 2` vertices - CORRECT

- *Citations.* Session 210 `graded-spectrum-gates/PROOF.md`, section 4-5: `P_T(q)` = number of ordered `q`-tuples of
  possibly empty directed paths whose vertex sets partition `V(T)`, `P_empty(q) = 1`, and eq. (5)
  `H(Gamma(T)) = H(T) + P_T(2)` "for nonempty `T`" - as quoted. THM-4892: the item is **(F) 3** ("Odd-subset gate
  formula"); its own display is the odd-subset form `H(Gamma(T)) = H(T) + 2 sum_{abs(S) odd} H(T[S]) H(T[V - S])`,
  and the all-subsets form is quoted there as Session 210's. Cite it as "THM-4892 (F) 3", not "statement 3". THM-4900
  lists the gate and `P_Gamma(q) = q[P(q+1) + q P(q) - P(q-1)]` as KNOWN from Session 210. The gate is
  `a -> b -> T -> a` in all three, as in the claim.
- *Proof.* `P_R(2) = 2 H(R) + (split sum)`; `2 H(R) = 2` and, by Lemma 4.5, split `= 2 (mod 4)`. Correct; `m >= 2`
  is needed (`m = 1`: `Gamma` is the 3-cycle, `H = 3`). The corollary does not need Lemma 4.5: by THM-4892 (F) 3,
  `H(Gamma R) - H(R)` is twice a sum of `2^(m-1)` odd numbers. Conversely this is a second proof of the congruence in
  (4.5). So L1 and the congruence of Lemma 4.5 are the same statement read through Session 210 (5), and both are one
  line from formulas of canon; say so.
- *Recomputation.* `H(Gamma R)` by a fresh DP on the `(m+2)`-vertex tournament against `H(R) + sum over all S`, and
  `P_R(2) = 0 mod 4`: every labelled tournament of orders 2-7 (2,131,018) and every class of orders 2-9 (198,947);
  counters `fail_gate`, `fail_gate4`: 0. The packet's `gate_mod4.py` covers the 531 classes of orders 2-7.
- *Residue already noted?* No file found (section 6).
- *Text.* "**Corollary (gate).** For every tournament `R` on `m >= 2` vertices, `P_R(2) = 0 (mod 4)` and
  `H(Gamma R) = H(R) (mod 4)`, where `Gamma(R) = (a -> b -> R -> a)` and `H(Gamma R) = H(R) + P_R(2)` (Session 210
  eq. (5); THM-4892 (F) 3). Proof: `P_R(2) = 2 H(R) +` the total of (4.5); or directly, `P_R(2)` is twice the sum of
  `H(R[S]) H(R[V \ S])` over the `2^(m-1)` odd subsets (THM-4892 (F) 3). Hence every gate tower `Gamma^r(R)`, `r >= 0`,
  keeps `H mod 4`. (At `m = 1` the residue changes: `H(Gamma(point)) = 3`.) FINITE-EXACT check: all classes of orders
  2-9, all labelled tournaments of orders 2-7 (audit 2)."

### L2. "(4.5) is the only affine relation" for `4 <= m <= 11` - CORRECT (FINITE-EXACT certificate)

- Own program `x6_affine.c` ((ell, beta) from the insertion decomposition, validated against direct counts by
  `x1_cube.c`): **every** labelled tournament for `m <= 7` (2,097,152 at `m = 7`), and 400 pseudo-random labelled
  tournaments (own generator and seed) for `m = 8, 9, 10, 11`. Affine dimension of `{beta_R}`:
  `0, 0, 5, 9, 14, 20, 27, 35, 44, 54` at `m = 2, ..., 11`, equal to `C(m,2) - 1` for `m = 2` and `4 <= m <= 11`; of
  `{(ell_R, beta_R)}`: `1, 2, 8, 13, 19, 26, 34, 43, 53, 64`, equal to `m + C(m,2) - 2` for `m = 2` and
  `4 <= m <= 11`. Full dimension is reached after 37, 50, 54, 65 random tournaments at `m = 8, ..., 11`. Independent
  confirmation without the insertion formula for `m <= 6` (`x1_cube m labelled`, from direct counts): same dimensions.
- `m = 3`: every `beta_R` is the triangle (dimension 0 of at most 2), and the pairs `(ell_R, beta_R)` have affine
  dimension 2 of at most 4 - as the packet's `out/affine_span.out` shows; say both.
- Typing: FINITE-EXACT (the lower bound is a certificate; the upper bound is Lemma 4.5). An all-`m` statement is not
  claimed and not proved. Use "affine", not "linear". Replacement text: m5.

### L3. Digit-degree law on a vertex cube - CORRECT (PROVED), and it can be made exact

- *Proof check.* `G'(X) = (H(R + v_X) - 1)/2 = (H(R)-1)/2 + sum (a_i/2) x_i - sum (mu_ij/2) x_i x_j` has integer
  coefficients (parities of Theorem C (1)) and is `>= 0` on the cube. Digit `j` of `N >= 0` is `C(N, 2^j) mod 2`
  (Lucas). `C(G', 2^j)` is a rational polynomial of degree `<= 2 * 2^j` in the `x_i`; reducing `x_i^2 = x_i` does not
  raise the degree; the multilinear polynomial that represents an integer-valued function on the cube is unique and has
  integer coefficients (Moebius inversion); its reduction mod 2 is the algebraic normal form. So digit `j` of `G'`,
  which is digit `j+1` of `H`, has `F_2`-degree `<= 2^(j+1)`. Every step is correct. This is the standard argument
  (digit `j` of an integer-valued polynomial of degree `d` has degree `<= d 2^j`; for `d = 1` it is the classical law
  for the Hamming weight that THM-478's note cites); type it PROVED and cite it as standard.
- *Sharpening (auditor; PROVED).* For `abs(U) = 2^k`, `k >= 1`, the coefficient of `prod_{i in U} x_i` in the algebraic
  normal form of digit `k` of `H(R + v_X)` is the Pfaffian of `beta_R` restricted to `U`. (With `t = 2^(k-1)`: the
  degree-`2t` multilinear part of `C(G', t)` comes only from `(quadratic part)^t / t!`, and each perfect matching of
  `U` arises `t!` times, so the coefficient is the hafnian of `mu/2` on `U`, which is the Pfaffian mod 2.) Hence
  **digit `k` has degree exactly `2^k` if and only if `rank(beta_R) >= 2^k`**. This contains both "exactly" claims:
  digit 1 has degree 2 for every `R` with `m >= 2` because `rank(beta_R) >= 2` (Lemma 4.5) - correct as claimed; and
  digit 2 has degree 4 on every class of order 6 because no class of order 6 has `rank(beta_R) = 2` - correct, also at
  order 7, but **not** an all-order law: it fails exactly on the rank-2 classes (3, 2, 2, 8, 7 classes at orders 4, 5,
  8, 9, 10).
- *Recomputation* (classes; direct counts at orders `<= 9`, insertion formula with one direct control per class at
  order 10): identical to `out/digit_degree.out` at orders 3-7.

  | order | digit 2 | digit 3 |
  |---|---|---|
  | 4 | {2: 2, 3: 1, 4: 1} | {2: 1, 3: 1, 4: 2} |
  | 5 | {2: 1, 3: 1, 4: 10} | {2: 1, 4: 9, 5: 2} |
  | 6 | {4: 56} | {4: 2, 5: 31, 6: 23} |
  | 7 | {4: 456} | {4: 2, 5: 9, 6: 319, 7: 126} |
  | 8 | {3: 2, 4: 6878} | {5: 2, 6: 548, 7: 3452, 8: 2878} |
  | 9 | {2: 2, 3: 6, 4: 191528} | {4: 2, 5: 7, 6: 907, 7: 32063, 8: 158557} |
  | 10 | {3: 7, 4: 9733049} | {5: 5, 6: 885, 7: 210347, 8: 9521819} |

  Digit 1 has degree 2 on every tournament of every run. No violation of `degree <= 2^k` anywhere (the bound for
  digit 3 is not vacuous at orders 9 and 10). The number of classes with digit-3 degree 8 equals the number of classes
  with `rank(beta_R) >= 8` (2,878; 158,557; 5,444,657 + 4,077,162 = 9,521,819), and the classes with digit-2 degree
  below 4 are exactly the rank-2 classes (2, 8, 7 at orders 8, 9, 10), as the criterion requires.

### L4. Deleting a vertex is not a bounded-rank perturbation of `beta_R` (negative probe) - numbers CORRECT

- Recomputed with `beta` obtained by a different route (parity of the even-order paths, (4.4);
  `x7_deletion_rank.py`), all classes of orders 4-7, all vertices: the four tables of `out/deletion_rank.out` are
  reproduced exactly. Order 7: rank of `beta_R|_(V-u) + beta_(R-u)`: `{0: 118, 2: 278, 4: 1864, 6: 932}`;
  `(rank beta_R, rank beta_(R-u))`: `{(4,4): 303, (4,6): 215, (6,4): 1387, (6,6): 1287}`. Order 5:
  `{(2,2): 7, (2,4): 3, (4,2): 31, (4,4): 19}`, so 3 pairs with `rank beta_R = 2 < rank beta_(R-u) = 4`. Order 6:
  `{0: 24, 2: 80, 4: 232}` and `{(4,2): 36, (4,4): 144, (6,2): 16, (6,4): 140}`.
- Typing: FINITE-EXACT, orders 4-7; counts are over (class representative, vertex) pairs, unweighted. "Not a
  bounded-rank perturbation" is an extrapolation; what the data show is that at each order 4-7 the difference attains
  the largest possible rank `2 floor((m-1)/2) = 2, 4, 4, 6`. Suggested sentence: "over all classes of orders 4-7 and
  all vertices `u` the difference `beta_R|_(V-u) + beta_(R-u)` attains the largest possible rank
  `2 floor((m-1)/2)` (at order 7: ranks 0, 2, 4, 6 with counts 118, 278, 1864, 932), and the rank of `beta` can go up
  as well as down under deletion."

### Planned citation of THM-5126 - the identification is right

THM-5126 (line 62 and eq. (1')): `Q_B(e,f)` = number of ordered pairs of disjoint non-empty paths covering `B`, the
first ending at `e`, the second starting at `f`; `W_ef = Q_B(e,f) + Q_B(f,e)`; `W(O)` = sum of `W_ef` over the pairs
of `O`. The README's `M_R[j,i]` is `Q_R(j,i)`, so `W_ef = mu_B(e,f)` and (4.5) reads `W(V(B)) = 2 (mod 4)`. Write
`W(V(B))` (in THM-5126, `V` is not a name for the vertex set of the base). THM-5126 itself notes only that `W(O)` is
even for every `O` (line 107).

---

## 9. Not checked

- `deletion_pair_corr.py 7` and the packet's two order-10 loops were not rerun (section 7).
- No literature lookup (Dickson, Lucas, Delsarte-Goethals for alternating forms; the recursion of m6 is proved here
  and checked by brute force instead). No search outside the repository for the box-norm chain or for the mod-4 law of
  `H` over random tournaments.
- The proofs of Session 210 eq. (5) and THM-4892 (F) 3 were not rechecked; the formula was verified numerically (L1).
- L2 at `m = 8..11` rests on pseudo-random certificates (sufficient for a lower bound); L4 only at orders 4-7.
- Nothing outside the audited parts (first audit); HYP-9659 (A), (C) and the window census were not revisited.
- File timestamps were read, not proved (section 5).

## 10. Rerun

    sh audit2/run_all.sh            # stages a b c d e f; about 1, 5, 10, 9, 12, 0 minutes with 6 processes on a
                                    # heavily loaded 8-core machine (stage e: 41 CPU-minutes)
    sh audit2/run_all.sh a          # small orders only (1 minute)
    sh audit2/run_all.sh f          # summaries, once the raw files of all stages are present
    # binaries and part files go to $W (default /tmp/claude-audit2-ensemble); P=<processes>, GENTOURNG=<path>

The files in this folder come from these stages, each run once with the final sources: a, b, c, e and f through the
script; stage d (order 10) by its command line typed by hand before the script existed, with a binary built from the
final `x1_cube.c`. The script was not run end to end in one go.

Programs: `x1_cube.c` + `x1_merge.py` (vertex cubes: bias, ranks, (4.2), (4.4), (4.5), THM-4177 (16), Theorem C (3),
gate, digit degrees, antipodal identity); `x2_box.c` + `x2_merge.py` (deletion pairs, box norms, two routes);
`x3_bruteforce_small.py` (permutation counting only); `x5_altrank.c` + `x5_oddedge_check.py` (alternating matrices);
`x6_affine.c` (L2); `x7_deletion_rank.py` (L4); `x8_nullmodel.py` (m3); `x9_compare_author.py` (comparison with the
packet's raw totals); `x10_image.py` (image of `R -> beta_R`, m5, m7). Outputs: `x1_raw_*.out`, `x1_summary.out`,
`x2_raw_*.out` (with the order-10 slice check of route 1), `x2_summary.out`,
`x3_bruteforce_small.out`, `x5_raw_altrank_2_8.out`, `x5_oddedge_check.out`, `x6_affine.out`, `x7_deletion_rank.out`,
`x8_nullmodel.out`, `x9_compare_author.out`, `x10_image.out`.
