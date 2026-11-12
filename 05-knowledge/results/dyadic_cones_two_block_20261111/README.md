# Dyadic neighbourhoods and two-block families of Hamiltonian-path counts; an attempt on Chan-Pak through rays

**Session:** mac-mini-2026-11-11-rays (follow-on of mac-mini-2026-11-11-gauge). **Canon:** THM-5145.
**Hypotheses:** HYP-9675. **Audit:** AUDIT-PENDING (pre-audit checkpoint).

Owner directive (2026-11-11): continue the session and extend it into another that aims to prove the Chan-Pak
conjecture (every odd number other than 7 and 21 is the number `H(T)` of Hamiltonian paths of a tournament).

**The conjecture is not proved here. It remains OPEN.** This note records what was tried, what was proved on
the way (three theorems giving new infinite families of path counts), and where each route stops.

## 0. What is in this note

| | statement | status |
|---|---|---|
| 1 | **Dyadic neighbourhoods.** Every odd `N` with `-791 <= N - 2^e <= 647` for some `e >= 1` is a path count, except 7 and 21. | PROVED (complete census of the rooted tournaments of order at most 10, every witness recounted literally, with the tail law of THM-5126 and the verified head of THM-4875) |
| 2 | **Two blocks.** For every odd `c <= 31`, every odd `r` with `|r| <= 275` and every `e >= 1` with `c 2^e + r >= 1`, the number `c 2^e + r` is a path count, except 7 and 21 (also `c <= 15`, `|r| <= 295`). | PROVED (same) |
| 3 | **Weight three, rows.** `1 + 2^a + 2^b` is a path count for `1 <= a <= 10` and every `b > a`, except 7 and 21. With the complete lines `c <= 4097` this leaves `a >= 11`, `b - a >= 13` of HYP-9658 (C). | PROVED (same) |
| 4 | One back arc: closed forms; with one out-neighbour the only cone offsets are 5 and 9. | PROVED |
| 5 | Blow-up templates in two parameters: census to order 8 (200 unit polynomials), no independent pair, no template of the shape that would prove weight three. | FINITE-EXACT |
| 6 | Type conjecture (HYP-9675): every pair of odd integers `(c, r)`, `c >= 1`, is the type of a rooted tournament. 92.8 percent of the types with `c <= 255`, `|r| < 2047` occur at order 10. | CONJECTURE |

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
  binary digits 10"; one negative offset (`2^K - 159`, THM-5085). What is added here is the complete census of
  cone roots, with arbitrary out-sets and offsets of both signs, to order 10, and its consequences.
- **Verified head (KNOWN).** THM-4875 (D): every odd number up to `2^32` other than 7 and 21 is a path count
  (audited); THM-1370: 7 and 21 are not.
- **Substitution (KNOWN).** THM-3121 (4), THM-5126 (B1): `H(K[M at x]) = sum_j c_j(K, x) P_j(M)`, `c_j` the number of
  Hamiltonian walks of `K` visiting `x` exactly `j` times, `P_j(M)` the number of ordered `j`-tuples of disjoint
  non-empty paths covering `M`; `P_2(TT_a) = 2^a - 2`; if no walk visits `x` three times, none visits it more often.
- **Data used.** The 6,880 classes of order 8 and the 456 of order 7
  (`../actual_hamilton_multiplicities_20261021/session24-multiplicities/`); THM-4810 (the least `H` of a strong
  tournament of order 3..8 is 3, 5, 9, 15, 25, 45).

## 2. Root types and the census

A rooted tournament with response `M = c 2^j` (`c` odd) and offset `r` has *type* `(c, r)` and *exponent* `j`; a
*cone root* is a root of type `(1, r)`. By the tail law a root of type `(c, r)` and exponent `j` puts
`c 2^e + r` into the spectrum for every `e >= j + 1`.

**Lemma (window).** If `O` is not empty then `c0 <= H(B) <= M` and `c0 >= 1 - 2M`.

*Proof.* `c0 <= H(B)` because `W >= 0`. For `y` in `O`, cutting a Hamiltonian path of `B` just before `y` (or not at
all if it starts at `y`) is an injection into the covers counted by `U_B(y)`, so `H(B) <= U_B(y) <= M`. And
`c0 + 2M = H(L o 0) >= 1`. (The first half is the window bound of THM-5109 in this setting.)

**Census** (`cone_census.c`, `block_census.c`). For every tournament `B` of order at most 9 - the 456 and 6,880
classes of orders 7 and 8, and all 1,761,280 one-vertex extensions of the classes of order 8, which contain
every tournament of order 9 up to isomorphism - and every out-set `O`, the offset and the response are computed
from `H(B)`, `U_B`, `W` (subset dynamic programming). So the census is over all rooted tournaments of order at
most 10.

- Cone roots (`out/summary_cones.txt`): 1,012 positive and 4,165 negative offsets below 20,000 in absolute
  value. **Every odd `r` with `-791 <= r <= 647` is the offset of a cone root**, with exponent at most 13. Missing
  below 1000: `649, 825` and `-793, -809, -891`. Offsets `1 + 2^a`: present for `a = 1, ..., 10` (exponents
  5, 5, 4, 7, 6, 8, 10, 10, 10, 11); offsets `2^a - 1`: present for `a <= 10`.
- Types (`out/summary_blocks.txt`): 243,055 of the 261,888 pairs (`c` odd `<= 255`, `r` odd, `|r| < 2047`), 92.8
  percent. **Every odd `c <= 31` with every odd `|r| <= 275`** is a type (exponent at most 12); every odd `c <= 15`
  with `|r| <= 295` (exponent at most 13). For `c <= 63` the box closes only to `|r| <= 9` (`c = 35`); `c = 245` has
  no root of type `(245, 1)`, in agreement with HYP-9658.
- **Literal verification** (`verify_cones.c`, a subset count written without the offset formula): for each of the
  118 + 729 + 5,177 cone witnesses of orders 8, 9, 10 the tournaments `L o 0^t`, `t = 1, 2, 3`, are built and
  their Hamiltonian paths counted; all 18,072 counts equal `r + 2^(j+t)`. The same for 10,179 witnesses of types
  with `c = 3, ..., 31` at order 9 (30,537 counts). Witnesses of order 10 for `c = 3, 5, 7`:
  `out/types_c3_c5_c7_roots_order10.txt`.

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
`r = 1`; Theorem 1 is the row `c = 1`, to a larger bound.

**Theorem 3 (weight three, ten rows).** For `1 <= a <= 10` and `b > a`, the number `1 + 2^a + 2^b` is a path count,
except `7 = 1 + 2 + 4` and `21 = 1 + 4 + 16`. Likewise `2^b + 2^a - 1` for `1 <= a <= 10`, `b > a`.

*Proof.* Cone roots with offsets `1 + 2^a` and `2^a - 1`, `a <= 10`, exponents at most 13; the numbers with
`b <= 13` are below `2^32`.

With the lines `c = 2^d + 1 <= 4097` this proves HYP-9658 (C) except when `a >= 11` and `b - a >= 13`.

**None of the three has positive density.** Up to `X` they contain `O(log X)`, `O(log X)` and `O(log X)` numbers
per type; the union over the 243,055 types is still `O(log X)`.

## 4. One back arc

Let `B` be the transitive tournament on `0 < 1 < ... < m-1` with one back arc `j -> i` of gap `k = j - i - 1 >= 1`;
`s = m - 1 - j`.

**Proposition.** `H(B) = 1 + 2^k`, `W_ij = 2^k`, and

    U_B(y) = 2^(q+s+1) (3^k + 2^(k+1) - 1)              (y < i,  q = i - y - 1)
    U_B(i) = 2^s (2^(k+1) - 1)
    U_B(y) = 2^s (2^(b+1) + 2^(a+1) 3^b - [a = 0])      (i < y < j,  a = y - i - 1,  b = j - y - 1)
    U_B(j) = 2^s (2^k + 1) ,     U_B(y) = 2^(m-1-y) (2^k + 1)   (y > j) .

Consequently: (a) `O = {i, j}` is a unit-offset root with response `3 * 2^(s+k)` for every `k` (the line `c = 3`
from its second point, as in THM-5126 (B2)); (b) with one out-neighbour (`O = {y}`, offset `1 + 2^k`), the response
is a power of two only for `k = 2`, `y < i` (offset 5, response `2^(q+s+5)`) and for `k = 3`, `y` the middle vertex of
the gap (offset 9, response `2^(s+4)`).

*Proof.* A path of `B` is increasing except for at most one step `j -> i`. Counting the covers `(P1, P2)` by
whether `P2`, `P1` use the back arc and by where the free vertices go gives the formulas (the factor `3^k`: a
gap vertex can lie in `P1`, in `P2` before `j`, or in `P2` after `i`); they are checked against the subset count
for all `m <= 11` (`one_back_arc.py`, 165 cases). (b): `3^k + 2^(k+1) - 1 = 2^u` forces `2^(k+1)` to divide `3^k - 1`,
whose 2-adic valuation is 1 for odd `k` and `2 + v_2(k)` for even `k`, so `k = 2`; and
`2^(b+1) + 2^(a+1) 3^b` is a power of two only if `a = b` and `3^b + 1` is one, so `a = b = 1`; the other three
expressions are odd multiples of powers of two.

So in the simplest infinite family the cones exist by two coincidences. With two and three back arcs the
cone offsets `3, 17, 33, 65` appear (`sparse_cones.c`); the general census of section 2 has no such restriction.

## 5. Blow-up templates in two parameters

If no Hamiltonian walk of `K` visits `x` (resp. `y`) three times, then
`H(K[TT_a at x, TT_b at y]) = P(2^a, 2^b)` with `P(u, v) = w_11 + w_21 (u - 2) + w_12 (v - 2) + w_22 (u - 2)(v - 2)`,
`w_jk` the number of Hamiltonian walks visiting `x` `j` times and `y` `k` times. A *unit pair template* has
`P(0, 0) = 1`; it gives the two-parameter family `1 + alpha 2^a + beta 2^b + gamma 2^(a+b)`, `a, b >= 1`. Products
`(1 + alpha u)(1 + beta v)` come from ordinal sums; a unit pair template that is not a product lives in a strong
tournament.

Census (`templates.c`; the polynomial is fitted at `(1,1), (2,1), (1,2), (2,2)` and checked at `(3,3)`, and at
`(4,2), (2,4)` for order 7; single templates are checked at `a = 3, 4`):

| order | tournaments | single templates (unit) | pair templates (unit) | distinct unit polynomials |
|---|---|---|---|---|
| 7 | 456 | 1,006 (122) | 1,446 (388) | 52 |
| 8 | 6,880 | 8,832 (324) | 8,481 (1,304) | 200 |
| 9, strong, `H` in a target list | 27,315 labelled | 28,859 (693) | 23,605 (3,276) | 59 |

Examples of proved two-parameter families (`a, b >= 1`): `1 + 2^a + 2^b + 6 * 2^(a+b)`; `1 + 2^a + 2^b + 14 * 2^(a+b)`;
`1 + 2^(b+1) (3 + 7 * 2^(a-1))`, which contains the lines `c = 3 + 7 * 2^i`, `i >= 1`, from their second point;
`1 + 2^(b+2) (1 + 21 * 2^(a-1))`; `1 + 2^(b+4) (1 + 3 * 2^(a-1))`. (All of these numbers with small parameters are
in the verified head; the families are infinite.)

**What the census excludes (FINITE-EXACT).**

- *No independent pair.* In no tournament of order at most 8, and in none of the 27,315 targeted ones of order
  9, is there a unit pair template with `gamma = 0` and `alpha beta != 0`. So no tournament of these orders carries
  two rays `1 + alpha 2^a + beta 2^b` that can be moved independently.
- *No template for weight three.* A unit pair template `1 + 2^j v + 2^i u v` would give all numbers
  `1 + 2^(b+j) + 2^(a+b+i)`. Its tournament is strong with `H = 1 + 2^(j+1) + 2^(i+2)`. None exists in order at most
  8 - which is every candidate with `H < 45`, in particular `H = 19, 35, 37` (`(j,i) = (0,2), (0,3), (1,3)`) - and
  none among the strong tournaments of order 9 with `H` one of the 38 target values up to 2113
  (`out/order9_target_values.txt`). The shapes `(j,i) = (0,0), (0,1), (1,1), (1,2)` are impossible for every
  order because they would contain 7 or 21.

## 6. What was tried towards the conjecture, and where each route stops

The tail formalism loses nothing: every path count is `H(L o 0)` for some root (add a sink `s` to a realiser and
a vertex beating only `s`). So Chan-Pak is the statement that every odd `N` other than 7 and 21 is
`c0 + M 2^t` for some rooted tournament and some `t >= 1`.

1. **Closure under `N -> 2N - 1` and `N -> 2N + 1`.** Every odd number above 3 is of one of the two forms. A ray
   is closed under the first map. The second map starts a new line and needs a root with unit offset and
   response `N`. Substituting `T` for a vertex of a fixed tournament gives `sum_j c_j P_j(T)`, a non-negative
   combination of the cover numbers of `T` with no constant term, so the only maps `H(T) -> alpha H(T) + beta` obtained
   in this way for every `T` are the products `beta = 0`. The step therefore has to use finer data of `T`; no class
   of realisers closed under both maps was found, and 7 and 21 show that none can contain realisers of 3 and of 11.
2. **Independent rays.** If for every `k` some tournament carried `k` independent rays of multiplier 4
   (`H = 1 + 4 sum_i 2^(a_i)`), every `N = 1 (mod 8)` would be a path count. Section 5: no two independent unit
   rays exist in order at most 8.
3. **A single two-parameter template for weight three.** Section 5: not of the simple shapes, up to order 9.
4. **Sparse families with closed forms.** Section 4: the seeds are coincidences of exponential Diophantine
   equations (`3^k + 2^(k+1) - 1 = 2^u`); with more back arcs more coincidences, no parametric family in the
   offset.
5. **A 2-adic law on the arc cube.** False (the previous note of this machine, HYP-9669): the divisibility of
   the switching differences is special to switching cubes.
6. **What does work, and its limit.** The census: almost every small type exists already at order 10, and the
   three theorems above are its content. Each type covers one number per exponent; any finite set of types has
   density zero. A proof along these lines needs roots of unbounded order, hence a recursive source of types.

**HYP-9675 (type conjecture).** For every pair of odd integers `(c, r)` with `c >= 1` there is a rooted
tournament of type `(c, r)`. (The column `r = 1` is HYP-9658 (A).) It implies that every number
`c 2^e + r` with `e` large is a path count; it does not imply Chan-Pak, which needs the exponent `j(c, r)` to be
small: for every odd `N` a splitting `N = c 2^e + r` with a root of type `(c, r)` and exponent below `e`. At order
10 the least exponent over the box `c <= 31`, `|r| <= 275` is at most 12.

## 7. Pre-registration

`PREREGISTRATION.md`: (Q1), (Q2), (Q4) CONFIRMED, (Q3) FAILED (the offsets 649 and 825). The registration was
written while the two runs were in progress and before any output was read; this is said in the file.

## 8. Reproduction

    cc -O2 -o cone_census cone_census.c && ./cone_census 8 0 20000 < out/classes_order8.txt     # roots of order 9
    ./cone_census 8 1 20000 < out/classes_order8.txt                                            # roots of order 10, about 15 minutes
    cc -O2 -o block_census block_census.c && ./block_census 8 1 2047 < out/classes_order8.txt
    cc -O2 -o verify_cones verify_cones.c && ./verify_cones < out/cones_roots_order10.txt
    python3 analyse_census.py cones out/cones_roots_order8.txt out/cones_roots_order9.txt out/cones_roots_order10.txt
    python3 one_back_arc.py
    cc -O2 -o templates templates.c && ./templates 8 0 < out/classes_order8.txt                 # pair templates, order 8

## 9. What is not claimed

- The Chan-Pak conjecture: OPEN. HYP-9658 (A), (C) and HYP-9675: OPEN.
- The three theorems rest on finite censuses and on two inherited results (the tail law, the verified head to
  `2^32`); they have density zero.
- The census of templates of order 9 is restricted to a list of values of `H`.
