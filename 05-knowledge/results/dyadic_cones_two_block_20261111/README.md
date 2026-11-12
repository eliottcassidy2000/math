# Dyadic neighbourhoods and two-block families of Hamiltonian-path counts; an attempt on Chan-Pak through rays

**Session:** mac-mini-2026-11-11-rays (follow-on of mac-mini-2026-11-11-gauge). **Canon:** THM-5145.
**Hypotheses:** HYP-9675. **Audit:** PASS WITH CORRECTIONS ([audit/AUDIT.md](audit/AUDIT.md); independent, own
programs). Theorems 1 to 3 and the census were reproduced. Five major findings (three on credit to earlier work,
one unsupported equivalence, one overstated verification) and nine minor ones; all are applied in this version
(section 10; MISTAKE-767).

Owner directive (2026-11-11): continue the session and extend it into another that aims to prove the Chan-Pak
conjecture (every odd number other than 7 and 21 is the number `H(T)` of Hamiltonian paths of a tournament).

**The conjecture is not proved here. It remains OPEN.** This note records what was tried, what was proved on
the way (three theorems giving infinite families of path counts), and where each route stops.

## 0. What is in this note

| | statement | status |
|---|---|---|
| 1 | **Dyadic neighbourhoods.** Every odd `N` with `-791 <= N - 2^e <= 647` for some `e >= 1` is a path count, except 7 and 21. | PROVED (complete census of the rooted tournaments of order at most 10, every witness recounted literally, with the tail law of THM-5126 and the verified head of THM-4875) |
| 2 | **Two blocks.** For every odd `c <= 31`, every odd `r` with `-275 <= r <= 275` and every `e >= 1` with `c 2^e + r >= 1`, the number `c 2^e + r` is a path count, except 7 and 21 (also `c <= 15`, `-295 <= r <= 295`). | PROVED (same; every type of the two boxes has a literally recounted witness on file) |
| 3 | **Weight three, rows.** `1 + 2^a + 2^b` is a path count for `1 <= a <= 10` and every `b > a`, except 7 and 21. With the complete lines `c <= 4097` this leaves `a >= 11`, `b - a >= 13` of HYP-9658 (C). | PROVED (same) |
| 4 | One back arc: closed forms of the responses; with one out-neighbour the only cone offsets are 5 and 9. | PROVED. On record before: the count `1 + 2^k` (classical), the factor `3^k` (THM-4879), the two cones (Session236) |
| 5 | Blow-up templates in two parameters: 52 and 200 distinct unit polynomials at orders 7 and 8; no template of the shape that would prove weight three with `H(K) < 75` at any order, or with `H(K) <= 2113` at order 9. | FINITE-EXACT. On record before (kind-pasteur S21, THM-5124): the census of two-module families to order 8 and the absence of an additive pair |
| 6 | The tail letters and the change of root as three integer matrices on four coordinates; no digit system for words of length at most 6. | FINITE-EXACT; an independent check of a result that opus S93 pushed the same evening (one-mark block law, THM-5150 reserved). Nothing new is claimed |
| 7 | Type conjecture (HYP-9675): every pair of odd integers `(c, r)`, `c >= 1`, is the type of a rooted tournament. At order 10: 92.8 percent of the types with `c <= 255`, `r` of absolute value below 2047; every odd `N` in `[3, 266979]` other than 7 and 21 is `c 2^e + r` for a type of exponent below `e`. | CONJECTURE (the columns `r = 1` and `r = -1` are HYP-9658 (A) and HYP-9662 (A)) |

## 1. Inheritance

- **The tail law and rays (KNOWN).** THM-5126 (A) (kind-pasteur S21, with Session216): for a rooted tournament
  `(L, v)`, `B = L - v`, `O = N+(v)`, put `c0 = H(B) - sum_{pairs {e,f} of O} W_ef` (the offset) and
  `M = sum_{y in O} U_B(y)` (the response); then `H(L o 0^t) = c0 + M 2^t` for every `t >= 1`, where `L o 0^t` adds a
  chain `z_1 -> ... -> z_t` below `B` and above `v` (every vertex of `B` beats every `z`, every `z` beats `v`).
  `U_B(y)` is the number of ordered pairs `(P1, P2)` of disjoint paths covering `B` with `P2` starting at `y` (`P1`
  possibly empty); `W_ef` is the number of ordered two-path covers with the first path ending at one of `e, f` and
  the second starting at the other.
- **Lines (KNOWN).** Every odd `N > 1` is `1 + c 2^s` with `c` odd; THM-5126 (B3) and HYP-9658 (B): every line with
  `c <= 4097` is complete except for 7 and 21 (unit-offset roots, `c0 = 1`). HYP-9658 (A) (every `c` has a seed) and
  (C) (weight three) are open.
- **Cones (KNOWN mechanism).** THM-5109 (kind-pasteur S20), THM-5085, Session236: a root whose response is a power
  of two gives `h + 2^K` for all large `K`; the dyadic ports of order at most 8 (181 triples; offsets that are
  path counts, e.g. 5, 19, 217, 223, 599 and their products); the window bound "a cone value begins with the
  binary digits 10".
- **Offsets of both signs (KNOWN in part; this credit was added after the audit, M2).**
  - `2^K - 159` (THM-5085).
  - THM-5130 (D), (E) (kind-pasteur S22): a pair of order 8 with offset `-1` and slope `2^10`, hence `2^(10+t) - 1`;
    towers of depth three with a power-of-two slope (`2^N - 1` for every `N >= 24`); an anti-seed (offset `-1`,
    slope `c 2^j`) for every odd `c <= 4097`. The numbers `2^k - 1` themselves are path counts by THM-4748 (4).
  - Session195: the rays `66 * 2^k - 17` and `88 * 2^k - 13`, i.e. the types `(33, -17)` and `(11, -13)`.
  - THM-5129 (B) (opus S91): for one out-neighbour, the complete lists of slopes, at every order, for the offsets
    `D <= 73` and `D = 129, 257`.
  - HYP-9662 (A) (kind-pasteur S22): every odd `c` has a seed and an anti-seed. These are the columns `r = 1` and
    `r = -1` of the type conjecture below. kind-pasteur's `portrait.c` computes offset and slope of every
    `(B, O)`; it was run for offsets `+1` and `-1`.
  - Same day, in parallel (opus S93, seeds lane, THM-5149 reserved): the census of unit pairs of both signs
    extended to order 10 and to odd `c` below 53,559.
  - **Added here:** the complete census of cone roots and of general types `(c, r)` among the roots of order at
    most 10, with arbitrary out-sets and offsets of both signs; Theorems 1 to 3 as stated; the gaps
    `649, 825, -793, -809, -891`. The audit's search found no earlier statement of these.
- **Verified head (KNOWN).** THM-4875 (D): every odd number up to `2^32` other than 7 and 21 is a path count
  (audited); an addendum records one further run to `2^35`, not recomputed independently; `2^32` suffices here.
  THM-1370: 7 and 21 are not path counts.
- **Substitution (KNOWN).** THM-3121 (4), THM-5126 (B1): `H(K[M at x]) = sum_j c_j(K, x) P_j(M)`, `c_j` the number of
  Hamiltonian walks of `K` visiting `x` exactly `j` times, `P_j(M)` the number of ordered `j`-tuples of disjoint
  non-empty paths covering `M`; `P_2(TT_a) = 2^a - 2`; if no walk visits `x` exactly three times, none visits it
  more often (THM-5126, audit A2).
- **Data used.** The 6,880 classes of order 8 and the 456 of order 7
  (`../actual_hamilton_multiplicities_20261021/session24-multiplicities/`). THM-4810: the least `H` of a strong
  tournament of order 3, ..., 9 is 3, 5, 9, 15, 25, 45, 75, and it increases with the order (floor law of Moon
  1972 and Busch 2006, as cited there).

## 2. Root types and the census

A rooted tournament with response `M = c 2^j` (`c` odd) and offset `r` has *type* `(c, r)` and *exponent* `j`; a
*cone root* is a root of type `(1, r)`. A root has a type only if `O` is not empty (`M >= 1`). By the tail law a
root of type `(c, r)` and exponent `j` puts `c 2^e + r` into the spectrum for every `e >= j + 1`.

**Lemma (window).** If `O` is not empty then `c0 <= H(B) <= M` and `c0 >= 1 - 2M`.

*Proof.* `c0 <= H(B)` because `W >= 0`. For `y` in `O`, cutting a Hamiltonian path of `B` just before `y` (or not at
all if it starts at `y`) is an injection into the covers counted by `U_B(y)`, so `H(B) <= U_B(y) <= M`. And
`c0 + 2M = H(L o 0) >= 1`. (KNOWN: the first half is the window bound of THM-5109 in this setting, and
`U_B(y) >= H(B)` is THM-5130 (F1).)

**Census** (`cone_census.c`, `block_census.c`). For every tournament `B` among the 456 and 6,880 classes of orders
7 and 8 and all 1,761,280 one-vertex extensions of the classes of order 8 (which contain every tournament of
order 9 up to isomorphism), and for every out-set `O`, the offset and the response are computed from `H(B)`,
`U_B`, `W` (subset dynamic programming). These are the roots of orders 8, 9 and 10. Roots of smaller order add
nothing: a source added to `B` outside `O` changes neither the offset nor the response (it is a source of
`L o 0^t`), so every type of a smaller root recurs with the same exponent at every larger order. The audit also
ran the bases of orders 1 to 6 directly (`audit/acensus.c`): no new type, no change in the gaps. So the census is
over all rooted tournaments of order at most 10.

- Cone roots (`out/summary_cones.txt`): 1,012 positive and 4,165 negative offsets below 20,000 in absolute
  value. **Every odd `r` with `-791 <= r <= 647` is the offset of a cone root**, with exponent at most 13. Missing
  below 1000: `649, 825` and `-793, -809, -891`. Offsets `1 + 2^a`: present for `a = 1, ..., 10` (exponents
  5, 5, 4, 7, 6, 8, 10, 10, 10, 11); offsets `2^a - 1`: present for `a <= 10`.
- Types (`out/summary_blocks.txt`): 243,055 of the 261,888 pairs (`c` odd `<= 255`, `r` odd, `|r| < 2047`), 92.8
  percent. **Every odd `c <= 31` with every odd `|r| <= 275`** is a type (exponent at most 12); every odd `c <= 15`
  with `|r| <= 295` (exponent at most 13). For `c <= 63` the box closes only to `|r| <= 9` (`c = 35`); `c = 245` has
  no root of type `(245, 1)`, in agreement with HYP-9658.
- **Coverage** (`coverage.py`, `out/coverage_by_types_order10.txt`). Every odd `N` with `3 <= N <= 266,979` other
  than 7 and 21 is `c 2^e + r` for a type `(c, r)` of this table whose exponent is below `e`; the first number not
  of this form is 266,981 (31 odd numbers of `[2^18, 2^19)` are not; beyond, the bounds `c <= 255`, `|r| < 2047` of
  the table bind). This gives no new path count - the head is verified to `2^32` - but it says that up to 266,979
  every path count `N >= 3` is `H(L o 0^t)`, `t >= 1`, for a root of order at most 10 with a non-empty out-set
  (used in section 6).
- **Literal verification** (`verify_cones.c`, a subset count written without the offset formula). For each of the
  118 + 729 + 5,177 cone witnesses of orders 8, 9, 10 the tournaments `L o 0^t`, `t = 1, 2, 3`, are built and
  their Hamiltonian paths counted; all 18,072 counts equal `r + 2^(j+t)`. The same for the witnesses of the types
  with `c = 3, ..., 31`: 10,179 at order 9 (30,537 counts), 5,670 at order 10 for `c = 3, 5, 7` (17,010 counts) and
  3,598 at order 10 for `c = 9, ..., 31` (10,794 counts). No failure. The last file was completed after the audit,
  which found that 293 types of the boxes of Theorem 2 (all with `c >= 9` and roots of order 10 only) had no
  witness on file, and recounted them with witnesses of its own (M4). Every cone offset in `[-791, 647]` and every
  type of the two boxes now has a witness line on file that was recounted literally (checked by a script over
  the six witness files).

## 3. Three theorems

**Theorem 1 (dyadic neighbourhoods).** Let `N` be an odd positive integer other than 7 and 21. If
`-791 <= N - 2^e <= 647` for some integer `e >= 1`, then `N = H(T)` for some tournament `T`.

*Proof.* Put `r = N - 2^e`. By the census there is a cone root with offset `r` and exponent `j <= 13`; by the tail
law `2^e + r` is a path count when `e >= j + 1`. If `e <= j` then `N <= 2^13 + 647 < 2^32`, and THM-4875 (D) applies.

**Theorem 2 (two blocks).** Let `c <= 31` be odd, `r` odd with `|r| <= 275`, `e >= 1`, and `N = c 2^e + r >= 1`,
`N` not 7 or 21. Then `N` is a path count. The same holds for `c <= 15`, `|r| <= 295`.

*Proof.* The same, with the root of type `(c, r)`; here `N <= 31 * 2^12 + 275 < 2^32` in the remaining cases.

In binary: `N` is a path count as soon as it consists of a block of at most five digits, a run of zeros, and a
block of at most eight digits (or is 275 or less below such a number). The lines of THM-5126 are the column
`r = 1` and the anti-seed lines of THM-5130 (E) the column `r = -1` (both KNOWN, and known for `c <= 4097`);
Theorem 1 is the row `c = 1`, to a larger bound.

**Theorem 3 (weight three, ten rows).** For `1 <= a <= 10` and `b > a`, the number `1 + 2^a + 2^b` is a path count,
except `7 = 1 + 2 + 4` and `21 = 1 + 4 + 16`. Likewise `2^b + 2^a - 1` for `1 <= a <= 10`, `b > a`.

*Proof.* Cone roots with offsets `1 + 2^a` and `2^a - 1`, `a <= 10`, exponents at most 13; the numbers with
`b <= 13` are below `2^32`.

What is added by Theorem 3 is the range `b - a >= 13`. For `b - a <= 12` the numbers `1 + 2^a + 2^b` lie on the
complete lines `c = 2^(b-a) + 1 <= 4097` (THM-5126, HYP-9658 (B)) and the numbers
`2^b + 2^a - 1 = (2^(b-a) + 1) 2^a - 1` on the anti-seed lines of THM-5130 (E). With those lines the theorem
proves HYP-9658 (C) except when `a >= 11` and `b - a >= 13`.

**None of the three has positive density.** Up to `X` they contain `O(log X)`, `O(log X)` and `O(log X)` numbers
per type; the union over the 243,055 types is still `O(log X)`.

## 4. One back arc

Let `B` be the transitive tournament on `0 < 1 < ... < m-1` with one back arc `j -> i` of gap `k = j - i - 1 >= 1`;
`s = m - 1 - j`.

On record before (audit M5): `H(B) = 1 + 2^k` is classical (Moon 1972; THM-250; THM-4748); the factor `3^k` -
a vertex strictly between an entry and a later exit has three admissible slots - is THM-4879 (B) and
THM-5109 (B); the two cones of (b) below are the seeds `(5, 0, 32)` and `(9, 2, 16)` of Session236 (THM-5109 (C));
THM-5129 (B) lists the slopes of every one-vertex source with these offsets. Added here: the closed forms of
`U_B(y)` and `W_ij` for this base and the classification (b).

**Proposition.** `H(B) = 1 + 2^k`, `W_ij = 2^k`, and

    U_B(y) = 2^(q+s+1) (3^k + 2^(k+1) - 1)              (y < i,  q = i - y - 1)
    U_B(i) = 2^s (2^(k+1) - 1)
    U_B(y) = 2^s (2^(b+1) + 2^(a+1) 3^b - [a = 0])      (i < y < j,  a = y - i - 1,  b = j - y - 1)
    U_B(j) = 2^s (2^k + 1) ,     U_B(y) = 2^(m-1-y) (2^k + 1)   (y > j) .

Consequently: (a) `O = {i, j}` is a unit-offset root with response `3 * 2^(s+k)` for every `k` (the line `c = 3`
from its second point, as in THM-5126 (B2)); (b) with one out-neighbour (`O = {y}`, offset `1 + 2^k`), the response
is a power of two only for `k = 2`, `y < i` (offset 5, response `2^(q+s+5)`) and for `k = 3`, `y` the middle vertex of
the gap (offset 9, response `2^(s+4)`).

*Proof.* A path of `B` is an increasing sequence, or two increasing sequences joined by the step `j -> i`; an
increasing sequence is a path unless it contains `i` and `j` and no vertex of the gap. In a path through the back
arc the part before `j` lies below `j` and the part after `i` lies above `i`. The gap has `k` vertices, `s`
vertices lie above `j`.

`H(B)`: the increasing path, and the paths through the back arc: the vertices below `i` precede `j`, those above
`j` follow `i`, each gap vertex chooses a side. So `H(B) = 1 + 2^k`. `W_ij`: a cover with first path ending at `j`
and second starting at `i` is a Hamiltonian path through the back arc cut at it (`2^k`); a cover with first path
ending at `i` and second starting at `j` leaves the gap uncovered. So `W_ij = 2^k`.

`U_B(y)` counts the pairs `(P1, P2)`, `P2` starting at `y`. Write `N0`, `N1`, `N2` for the numbers of pairs in
which no path, `P1`, `P2` uses the back arc. If `P2` is increasing, every vertex below `y` is in `P1`.

- `y > j`. `N0`: any subset of the `m - 1 - y` vertices above `y` joins `P2`; `P1` contains `i`, `j` and the gap.
  `N1`: the same choice for `P2`, then `2^k` choices for the gap inside `P1`. `N2 = 0` (`P2` would have to pass
  from `y` up to `j`). Total `2^(m-1-y) (1 + 2^k)`.
- `y = j`. `N0 = 2^s`. `N2`: `P2 = j, i, ...` continues with any subset of the gap and of the `s` upper vertices:
  `2^(k+s)`. `N1 = 0`. Total `2^s (1 + 2^k)`.
- `y = i`. `N0`: `P2 = {i}` with a subset of the `k + 1 + s` vertices above `i`, not allowed when it contains `j`
  and no gap vertex: `2^(k+1+s) - 2^s`. `N1 = N2 = 0` (`i` would be preceded by `j`). Total `2^s (2^(k+1) - 1)`.
- `y < i`. `N0`: the `q + k + s + 2` vertices above `y` are distributed over two increasing sequences; not allowed
  when `i`, `j` are in one of them and the whole gap in the other: `2^(q+k+s+2) - 2 * 2^(q+s)`. `N2`: the `q`
  vertices between `y` and `i` are in `P1` or in `P2` before `j`; a gap vertex is in `P1`, in `P2` before `j` or in
  `P2` after `i`; an upper vertex is in `P1` or in `P2`: `2^q 3^k 2^s`. `N1`: the same count with the roles of the
  two paths exchanged. Total `2^(q+s+1) (2^(k+1) - 1 + 3^k)`.
- `i < y < j`, with `a` gap vertices below `y` and `b` above (`a + b + 1 = k`). `N0`: the `b + 1 + s` vertices above
  `y` are distributed; `P1` contains `i` and the `a` lower gap vertices, so the only forbidden case is `a = 0`,
  `j` in `P1`, the upper gap in `P2`: `2^(b+1+s) - [a = 0] 2^s`. `N2`: `P2` runs from `y` to `j`, then `i`; a lower
  gap vertex is in `P1` or after `i`, an upper gap vertex has three slots, an upper vertex two:
  `2^a 3^b 2^s`. `N1`: the lower gap vertices are in `P1` before `j` or after `i`, the upper gap vertices have
  three slots, the upper vertices two: `2^a 3^b 2^s`. Total `2^s (2^(b+1) - [a = 0] + 2^(a+1) 3^b)`.

The formulas are also checked against the subset count for all `m <= 11` (`one_back_arc.py`, 165 cases) and by
the audit's own program for `m <= 9`.

(a): `c0 = 1 + 2^k - W_ij = 1` and `M = U_B(i) + U_B(j) = 3 * 2^(s+k)`.

(b): `2^k + 1` and `2^(k+1) - 1` are odd and larger than 1. If `3^k + 2^(k+1) - 1 = 2^u` then `u > k + 1`, so
`2^(k+1)` divides `3^k - 1`; the 2-adic valuation of `3^k - 1` is 1 for odd `k` and `2 + v_2(k)` for even `k`, and
`2 + v_2(k) >= k + 1` holds only for `k = 2` (`9 + 8 - 1 = 16`). For `a = 0` the bracket
`2^(b+1) + 2 * 3^b - 1` is odd and larger than 1. For `a >= 1`, `a != b`, the bracket is `2^(min(a,b)+1)` times an
odd number larger than 1. For `a = b >= 1` it is `2^(a+1) (3^b + 1)`, and `3^b + 1` is 2 or 4 modulo 8, hence a
power of two only for `b = 1` (`k = 3`, response `2^(s+4)`).

So in the simplest infinite family the cones exist by two coincidences. With two and three back arcs the
cone offsets `3, 17, 33, 65` appear (`sparse_cones.c`); the general census of section 2 has no such restriction.

## 5. Blow-up templates in two parameters

**Formula and hypothesis (KNOWN: THM-3121 (4), THM-5126 (B1); kind-pasteur S21).** Let `w_jk` be the number of
Hamiltonian walks of `K` that visit `x` exactly `j` times, `y` exactly `k` times and every other vertex once. Then
`H(K[TT_a at x, TT_b at y]) = sum_{j,k} w_jk s(a, j) s(b, k)`, where `s(a, j)` is the number of surjections of an
`a`-set onto a `j`-set (the ordered covers of `TT_a` by `j` non-empty paths). If no walk visits `x` or `y` three or
more times (`w_jk = 0` whenever `max(j, k) >= 3`), then

`H(K[TT_a at x, TT_b at y]) = P(2^a, 2^b)`, `P(u, v) = w_11 + w_21 (u - 2) + w_12 (v - 2) + w_22 (u - 2)(v - 2)`.

A *unit pair template* has `P(0, 0) = 1`; it gives the two-parameter family
`1 + alpha 2^a + beta 2^b + gamma 2^(a+b)`, `a, b >= 1`. Products `(1 + alpha u)(1 + beta v)` come from ordinal sums. A
unit pair template that is not a product lives in one strong component of its tournament; the other components
have one Hamiltonian path each, and the component carries the same polynomial.

**On record before (audit M1).** The census of these two-module families over every class of order at most 8 and
every ordered pair of vertices is kind-pasteur S21's (`modfam.c`, `kps_s21_two_adic_lines_20261110`, README
section 4: 406,412 triples, 3,477 distinct exact tuples), with the statement that none is additive (`gamma = 0`,
`alpha beta != 0`) for any offset, in agreement with THM-5124 (opus S90, no additive pair). The first exclusion
below is therefore KNOWN, in a stronger form. Added here: the counts of distinct unit polynomials, the targeted
census at order 9, and the exclusion of the weight-three shape.

**The test of the program is sufficient (audit m4).** `templates.c` fits `P` at `(a, b) = (1,1), (2,1), (1,2), (2,2)`,
which determines `w_11, w_21, w_12, w_22`, and compares with the literal count at `(3, 3)`. There
`H - P = 6 w_31 + 6 w_13 + 36 (w_32 + w_23 + w_33) >= 0`, so equality means `w_jk = 0` for `j, k <= 3`,
`max(j, k) = 3`. The number of walks of `K[TT_b at y]` visiting `x` exactly `j` times is `sum_k w_jk s(b, k)`. By the
three-visit lemma applied to `x` in `K[TT_3 at y]`, `w_jk = 0` for `j >= 3`, `k <= 3`; by symmetry `w_jk = 0` for
`j <= 3`, `k >= 3`. Hence `w_3k = 0` for every `k`, and the lemma applied to `x` in `K[TT_b at y]` for every `b` gives
`w_jk = 0` for all `j >= 3`; the case `j <= 2`, `k >= 3` is already covered. So a pair that passes the test at
`(3, 3)` is exact for all `a, b >= 1`. The further checks of the program (at `(4, 2)`, `(2, 4)` for order 7, and
at `a = 4` for single templates) are redundant.

Census (`templates.c`, `templates9.c`):

| order | tournaments | single templates (unit) | pair templates (unit) | distinct unit polynomials |
|---|---|---|---|---|
| 7 | 456 | 1,006 (122) | 1,446 (388) | 52 |
| 8 | 6,880 | 8,832 (324) | 8,481 (1,304) | 200 |
| 9, strong, `H` one of 38 values (first run) | 27,315 labelled | 28,859 (693) | 23,605 (3,276) | 59 |
| 9, strong, `H` one of 10 further values (run after the audit, m5) | 12,136 labelled | 11,880 (288) | 5,280 (710) | 13 |

(Polynomials are counted up to the exchange of `u` and `v`. A tournament of order below 8 with added sources has
the same templates, so the row of order 8 covers every smaller order.)

Examples of proved two-parameter families (`a, b >= 1`): `1 + 2^a + 2^b + 6 * 2^(a+b)`; `1 + 2^a + 2^b + 14 * 2^(a+b)`;
`1 + 2^(b+1) (3 + 7 * 2^(a-1))`, which contains the lines `c = 3 + 7 * 2^i`, `i >= 1`, from their second point;
`1 + 2^(b+2) (1 + 21 * 2^(a-1))`; `1 + 2^(b+4) (1 + 3 * 2^(a-1))`. (All of these numbers with small parameters are
in the verified head; the families are infinite.)

**What the census excludes (FINITE-EXACT).**

- *No independent pair (KNOWN to order 8, see above).* In no tournament of order at most 8, and in none of the
  39,451 targeted ones of order 9, is there a unit pair template with `gamma = 0` and `alpha beta != 0`. So no
  tournament of these orders carries two rays `1 + alpha 2^a + beta 2^b` that can be moved independently.
- *No template for weight three (added here).* A unit pair template `1 + 2^j v + 2^i u v` would give all numbers
  `1 + 2^(b+j) + 2^(a+b+i)`. It is not a product, so it lives in a strong tournament `K` with
  `H(K) = P(2, 2) = 1 + 2^(j+1) + 2^(i+2)`. A strong tournament of order at least 9 has `H >= 75`, so a candidate
  with `H(K) < 75` has order at most 8, and the census of order 8 is complete: **no such template exists with
  `H(K) < 75`, at any order** (the candidate values are 9, 13, 17, 19, 25, 33, 35, 37, 41, 49, 65, 67, 69, 73; of
  the other values of the form, 7 and 21 are not path counts and 11 belongs only to an impossible shape).
  At order 9 the candidate values up to 2113 are 43 numbers between 81 and 2113; the first run searched 33 of
  them (its list of 38 values, `out/order9_target_values.txt`, also contains 49, 65, 67, 69, 73, which no strong
  tournament of order 9 attains), and the audit pointed out that 257, 385, 513, 641, 769, 1025, 1153, 1281, 1537,
  2049 were missing; the second run (`out/order9_target_values_b.txt`, `out/pair_templates_order9_targets_b.txt`)
  searched these. **None exists at order 9 with `H(K) <= 2113`.** Values of `H` above 2113 at order 9 and all
  orders above 9 are not searched. The shapes `(j, i) = (0,0), (0,1), (1,1), (1,2)` are impossible for every order
  because they would contain 7 or 21.

## 6. What was tried towards the conjecture, and where each route stops

**What the tail formalism can express (corrected after the audit, M3).** A tournament is `L o 0^t` with `t >= 1` for
a root with a non-empty out-set exactly when its last strong component has at least three vertices and contains
a vertex of out-degree one. (In `L o 0^t` the vertex `z_t` beats only `v`, and `v -> y -> z_t -> v` for `y` in `O`;
conversely, if `z` in the last component beats only `v`, delete `z` and root at `v`.) Since `H` is the product over
the strong components and does not depend on their order, the following are equivalent for an odd `N`:
`N = c0 + M 2^t` for a rooted tournament with `M >= 1` and `t >= 1`; `N` has a splitting `N = c 2^e + r` with a root
of type `(c, r)` and exponent below `e`; `N` has a realiser with a strong component containing a vertex whose
out-degree inside the component is one. This holds for every odd `N` in `[3, 266979]` other than 7 and 21
(section 2). Whether every path count `N >= 3` has such a realiser is not proved. (If `M = 0` is allowed the reformulation is empty: add a sink `v` and a vertex beating
only `v` to a realiser; then `O` is empty and `c0 = N`. The number 1 has no representation with `M >= 1`: in a
transitive tournament the vertex of out-degree one beats the sink. The pre-audit version of this paragraph
claimed that the formalism "loses nothing"; that was this empty reformulation.)

1. **Closure under `N -> 2N - 1` and `N -> 2N + 1`.** Every odd number above 3 is of one of the two forms. A ray
   is closed under the first map. The second map starts a new line and needs a root with unit offset and
   response `N`. Substituting `T` for a vertex of a fixed tournament gives `sum_j c_j P_j(T)`, a non-negative
   combination of the cover numbers of `T` with no constant term, so the only maps `H(T) -> alpha H(T) + beta` obtained
   in this way for every `T` are the products `beta = 0`. The step therefore has to use finer data of `T`; no class
   of realisers closed under both maps was found, and 7 and 21 show that none can contain realisers of 3 and of 11.
2. **Independent rays.** If for every `k` some tournament carried `k` independent rays of multiplier 4
   (`H = 1 + 4 sum_i 2^(a_i)`), every `N = 1 (mod 8)` would be a path count. Section 5: no two independent unit
   rays exist in order at most 8 (KNOWN: kind-pasteur S21, THM-5124).
3. **A single two-parameter template for weight three.** Section 5: not of the simple shapes with `H(K) < 75`, nor
   at order 9 with `H(K) <= 2113`.
4. **Sparse families with closed forms.** Section 4: the seeds are coincidences of exponential Diophantine
   equations (`3^k + 2^(k+1) - 1 = 2^u`); with more back arcs more coincidences, no parametric family in the
   offset.
5. **A 2-adic law on the arc cube.** False (the previous note of this machine, HYP-9669): the divisibility of
   the switching differences is special to switching cubes.
6. **A digit system inside the word semigroup (an independent check; nothing new is claimed).** The tail
   letters and the change of root act linearly on the state `(alpha, M, E, d)` of a rooted tournament (`alpha = H(L)`,
   `M` the response, `E` the number of paths ending at the root, `d = H(B)`): tail letter with `z -> v`:
   `(alpha + M + d - E, 2M, d, d)`; tail letter with `v -> z`: `(alpha, 2M + d, 0, d)`; tail letter followed by
   re-rooting at the new vertex: `(alpha + M + d - E, M + d, alpha - E, alpha)` (`semigroup.py`: 300 random words of
   length at most 10 agree with literal counts).
   - *Credit.* These three maps are the one-vertex operations `a`, `c`, `b` of the one-mark block law that opus
     S93 pushed the same evening (`04-computation/experiments/s93_chanpak_wave_20261111/digits/PROOFS_finish.txt`,
     P1; state `(h, e, B, g)` with `B = M + d`; THM-5150 reserved), which treats arbitrary blocks. They were found
     independently here; the priority is S93's. Both are matrix forms of the O/Z/J words with root changes
     (THM-4888, THM-4904). Chan-Pak as cofiniteness of the orbit of a matrix semigroup is HYP-9640 (kind-pasteur
     S18); see also the four-state gadget law THM-5130 (A) and kind-pasteur S25's block semigroup (THM-5148,
     same day).
   - *Search.* If two words acted on a common invariant line as `alpha -> 2 alpha - 1` and `alpha -> 2 alpha + 1`, or
     with any complete system of digits, a full arithmetic progression of path counts would follow. For words of
     length at most 6 (1,092 words) the only plane on which two independent words act as
     `[[lambda, b], [0, 1]]` is `E = d = 0`, which contains no state of a tournament (`planes.py`). Eigen-directions
     shared by several words do exist (45, `lines.py`; for instance `(3, 2, 1, 1)` is an eigenvector of `0r = r1`
     with eigenvalue 3 and of `10r` with eigenvalue 5), but in no case is there a second direction with
     eigenvalue one in a common plane. opus S93's digits lane reports the same absence (no exact pair of
     non-commuting affine maps) in every finite system it searched, with blocks up to order 5, and its thinorbit
     lane that no published local-global theorem applies to these semigroups. The comparison with problems of
     Zaremba type (which integers are first coordinates of an orbit of a finitely generated matrix semigroup;
     the known methods give density one and not every integer) is on record in the repository since
     kind-pasteur S11 (see HYP-9661, Session205 and the notes of opus S90 and S91) and is an ANALOGY; no theorem
     of that kind is proved or used here.
7. **What does work, and its limit.** The census: almost every small type exists already at order 10, and the
   three theorems above are its content. Each type covers one number per exponent; any finite set of types has
   density zero. A proof along these lines needs roots of unbounded order, hence a recursive source of types.

**HYP-9675 (type conjecture).** For every pair of odd integers `(c, r)` with `c >= 1` there is a rooted
tournament of type `(c, r)`. The column `r = 1` is HYP-9658 (A); the columns `r = 1` and `r = -1` together are
HYP-9662 (A) (kind-pasteur S22), verified for every odd `c <= 4097` (THM-5130 (E)). The conjecture implies that
every number `c 2^e + r` with `e` large is a path count. It does not imply Chan-Pak. The *quantitative form* - for
every odd `N >= 3` other than 7 and 21 a splitting `N = c 2^e + r` with a root of type `(c, r)` and exponent below
`e` - implies Chan-Pak. By the first paragraph of this section it is equivalent to Chan-Pak together with the
statement that every path count `N >= 3` has a realiser with a strong component containing a vertex of out-degree
one; that statement is checked up to 266,979 and OPEN. (The pre-audit version asserted the equivalence with
Chan-Pak outright.) At order 10 the least exponent over the box `c <= 31`, `|r| <= 275` is at most 12.

## 7. Pre-registration

`PREREGISTRATION.md`: (Q1), (Q2), (Q4) CONFIRMED, (Q3) FAILED (the offsets 649 and 825). The registration was
written while the two runs were in progress and before any output was read; this is said in the file. The audit
checked the commit order and that the prediction text was not changed afterwards; one in-sample sentence of
the file was wrong and carries an erratum (m7).

## 8. Reproduction

    cc -O2 -o cone_census cone_census.c && ./cone_census 8 0 20000 < out/classes_order8.txt     # roots of order 9
    ./cone_census 8 1 20000 < out/classes_order8.txt                                            # roots of order 10, about 15 minutes
    cc -O2 -o block_census block_census.c && ./block_census 8 1 2047 < out/classes_order8.txt
    for c in 9 11 13 15 17 19 21 23 25 27 29 31; do ./cone_census 8 1 300 $c < out/classes_order8.txt; done   # type witnesses, order 10
    cc -O2 -o verify_cones verify_cones.c && ./verify_cones < out/cones_roots_order10.txt
    ./verify_cones < out/types_c9_to_c31_roots_order10.txt
    python3 analyse_census.py cones out/cones_roots_order8.txt out/cones_roots_order9.txt out/cones_roots_order10.txt
    python3 coverage.py 4194304
    python3 one_back_arc.py
    cc -O2 -o templates templates.c && ./templates 8 0 < out/classes_order8.txt                 # pair templates, order 8
    cc -O2 -o ext9 ext9.c && ./ext9 257 385 513 641 769 1025 1153 1281 1537 2049 < out/classes_order8.txt > t9.txt
    cc -O2 -o templates9 templates9.c && ./templates9 9 0 < t9.txt                              # pair templates, order 9
    python3 semigroup.py; python3 planes.py 6; python3 lines.py 6                               # needs sympy

The audit's own programs and tables are in `audit/` (`acensus.c`, `averify.c`, `canon8.c`, `aoneback.py`,
`ablow.c`).

## 9. What is not claimed

- The Chan-Pak conjecture: OPEN. HYP-9658 (A), (C) and HYP-9675: OPEN.
- The three theorems rest on finite censuses and on two inherited results (the tail law, the verified head to
  `2^32`); they have density zero.
- That every path count `N >= 3` is `H(L o 0^t)` for a root with a non-empty out-set (checked to 266,979 only), and
  hence the equivalence of the quantitative type statement with Chan-Pak.
- Novelty of the template census to order 8, of the absence of an additive pair, of cones and rays with negative
  offsets as such, of the columns `r = 1` and `r = -1`, of the count `1 + 2^k`, of the three-slot factor `3^k`, of
  the cones with offsets 5 and 9, and of the matrix form of the tail letters: all on record elsewhere (sections
  1, 4, 5, 6).
- The census of templates of order 9 is restricted to `H <= 2113`.

## 10. Audit and corrections

Independent audit by a fresh agent with its own census, literal recount and one-back-arc check
([audit/AUDIT.md](audit/AUDIT.md)): **PASS WITH CORRECTIONS.** Its census reproduces the type tables of orders 9
and 10 entry for entry (87,182 and 243,055 types), the class list of order 8 is complete, and every recount
agreed. No error was found in Theorems 1 to 3, the window lemma, the one-back-arc formulas or the template
examples. Findings and what was done:

| | finding | applied |
|---|---|---|
| M1 | The two-module template census to order 8 and "no additive pair" are kind-pasteur S21's (`modfam.c`) and THM-5124's; they were typed as new. | Credited in sections 0, 5, 6, 9; the new part is delimited. |
| M2 | Cones and rays with negative offsets and the columns `r = 1`, `r = -1` are partly on record (THM-5130, THM-4748, Session195, THM-5129, HYP-9662). | Credited in sections 1, 3, 6; the claim about "numbers beginning with the digits 11" is withdrawn as a novelty. |
| M3 | "The quantitative type statement is equivalent to Chan-Pak" was unsupported, and false at `N = 1`. | Replaced by "implies", with the exact missing statement and a check to 266,979 (sections 2, 6, 9; HYP-9675). |
| M4 | 293 types of Theorem 2 had no witness on file, so "every witness recounted" overstated the verification. | Witness file completed and recounted (3,598 witnesses, 10,794 counts, no failure); every type of the boxes now has a recounted witness. |
| M5 | Parts of section 4 are classical or known (`1 + 2^k`; the factor `3^k`; the two cones). | Credited in sections 0, 4, 9. |
| m1 | Bases of order at most 6 were not run. | Covered by the source-vertex remark and by the audit's own run (section 2). |
| m2 | `U_B(y) >= H(B)` is also THM-5130 (F1). | Credited. |
| m3 | The proof of section 4 was a sketch. | Written out, case by case. |
| m4 | The template hypothesis must read "three or more times"; the sufficiency of the `(3, 3)` test was not argued. | Hypothesis restated; argument added (section 5). |
| m5 | Ten candidate values up to 2113 were not searched at order 9; five of the 38 listed values are vacuous there. | Second run on the ten values (no template); the statement now covers every candidate with `H(K) < 75` at any order and `H(K) <= 2113` at order 9. |
| m6 | The floor should be cited at orders 9 and above; "19, 35, 37" was only a subset of the candidates. | Corrected (sections 1, 5). |
| m7 | A wrong in-sample sentence in the pre-registration ("fewer than 181"). | Erratum added to the file; the predictions and their scoring are unchanged. |
| m8, m9 | Outcomes confirmed; the head also reaches `2^35` by one run. | Noted in section 1. |

After the audit, and not seen by the auditor: the coverage computation (`coverage.py`), the second template run
at order 9, the written-out proofs of sections 4 and 5, the strong-component characterisation of section 6, and
item 6 of section 6 with its credit to opus S93 (whose lane folders were pushed while the audit was running).
