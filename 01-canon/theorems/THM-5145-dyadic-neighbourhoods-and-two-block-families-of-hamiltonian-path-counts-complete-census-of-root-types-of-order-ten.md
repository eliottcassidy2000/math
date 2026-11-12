---
id: THM-5145
title: "Dyadic neighbourhoods and two-block families of Hamiltonian-path counts. A rooted tournament (L, v) with offset c0 = r and response M = c 2^j (c odd) - type (c, r) - puts c 2^e + r into the spectrum for every e >= j + 1 (the tail law, KNOWN). Complete census of the rooted tournaments of order at most 10: every odd r with -791 <= r <= 647 is the offset of a root whose response is a power of two; every odd c <= 31 and odd r with |r| <= 275 is a type. Hence (1) every odd N with -791 <= N - 2^e <= 647 for some e >= 1 is a path count, except 7 and 21; (2) every number c 2^e + r >= 1 with c <= 31 odd, |r| <= 275 odd, e >= 1 is a path count, except 7 and 21; (3) 1 + 2^a + 2^b is a path count for 1 <= a <= 10 and every b > a, except 7 and 21. (4) For the transitive tournament with one back arc, closed forms of U_B, and the only roots with one out-neighbour and a power-of-two response have offsets 5 and 9 (the two cones KNOWN, Session236). (5) Two-parameter blow-up templates (census to order 8 and the absence of an additive pair KNOWN: kind-pasteur S21, THM-5124): 200 distinct unit polynomials at order 8, and no template 1 + 2^j v + 2^i u v with H(K) < 75 at any order or with H(K) <= 2113 at order 9"
status: >
  PROVED: (1)-(3) from a complete FINITE-EXACT census (every base B among the classes of order 7 and 8 and all
  1,761,280 one-vertex extensions of the classes of order 8, every out-set O; smaller roots recur at larger
  orders, and the audit ran them directly) together with two inherited results: the tail law
  H(L o 0^t) = c0 + M 2^t (THM-5126 (A)) and the verified head, every odd number up to 2^32 other than 7 and 21
  (THM-4875 (D)); the largest exponent j in the ranges used is 13. (4) PROVED by a case-by-case count of covers
  (written out in the note), the 2-adic valuation of 3^k - 1, and 3^b + 1 = 2^u only for b <= 1. (5) FINITE-EXACT.
  Machine checks: every cone witness (118, 729, 5,177 at root orders 8, 9, 10) recounted literally at tail lengths
  1, 2, 3 by a subset count written without the offset formula (18,072 counts, no failure); the witnesses of the
  types with c = 3..31: 10,179 at order 9 (30,537 counts), 5,670 (c = 3, 5, 7) and 3,598 (c = 9..31) at order 10
  (17,010 and 10,794 counts); every cone offset of (1) and every type of the boxes of (2) has a recounted witness
  on file (the last file was completed after the audit, finding M4); the closed forms of (4) for all m <= 11
  (165 cases); templates fitted at four points and checked at (3, 3), which is sufficient by the three-visit
  lemma of THM-5126.
  KNOWN, credited: the tail law, rays, unit-offset roots and the complete lines c <= 4097 (THM-5126, HYP-9658;
  kind-pasteur S21; Session216); the cone mechanism, the window bound and the dyadic ports of order <= 8
  (THM-5109, kind-pasteur S20; Session236); fixed-port tails and the negative offset 2^K - 159 (THM-5085);
  negative offsets and anti-seeds: offset -1 with slope 2^10, 2^N - 1 by towers, anti-seeds for every odd
  c <= 4097 (THM-5130 (D), (E), kind-pasteur S22), 2^k - 1 (THM-4748 (4)), the rays 66 * 2^k - 17 and 88 * 2^k - 13
  (Session195), complete slope lists of one-vertex sources for the offsets D <= 73, 129, 257 (THM-5129 (B)), the
  columns r = 1 and r = -1 as a conjecture (HYP-9662 (A)); U_B(y) >= H(B) (THM-5130 (F1)); in (4): H(B) = 1 + 2^k
  (Moon 1972; THM-250; THM-4748), the three-slot factor 3^k (THM-4879 (B), THM-5109 (B)), the cones with offsets
  5 and 9 (Session236 seeds (5,0,32), (9,2,16); THM-5109 (C)); in (5): the census of two-module families over
  every class of order <= 8 and the absence of an additive pair (kind-pasteur S21, modfam.c; THM-5124, opus S90);
  the substitution law and the three-visit lemma (THM-3121, THM-5126 (B1)); the verified ranges (THM-4875,
  THM-4885); 7 and 21 are not path counts (THM-1370); strong floors (THM-4810; Moon 1972, Busch 2006).
  NEW to the repository as far as the author's and the audit's searches found: the complete census of cone
  roots and of types (c, r) to order 10 with offsets of both signs, and the gaps 649, 825, -793, -809, -891;
  (1), (2), (3) as stated (in (3) the range b - a >= 13; in (2) the columns other than r = 1, -1); in (4) the
  closed forms of U_B and W for this base and the classification of the power-of-two responses; in (5) the
  counts of distinct unit polynomials, the targeted census at order 9 and the exclusion of the weight-three
  shape.
  NOT CLAIMED: the Chan-Pak conjecture (OPEN); positive density (each family has O(log X) members up to X);
  completeness of the template census at order 9 (restricted to H <= 2113) or above.
  Pre-registered (PREREGISTRATION.md of the source packet, written while the order-10 runs were in progress and
  before any output was read): Q1, Q2, Q4 CONFIRMED, Q3 FAILED (offsets 649, 825).
  AUDIT: PASS WITH CORRECTIONS (independent census, literal recount and one-back-arc check; (1)-(3) and the
  census confirmed; five major findings - three on credit (M1, M2, M5), one unsupported equivalence in HYP-9675
  (M3), one overstated verification (M4) - and nine minor ones, all applied; MISTAKE-767).
session: mac-mini-2026-11-11-rays
source: 05-knowledge/results/dyadic_cones_two_block_20261111/README.md
related:
  - THM-5126
  - HYP-9658
  - THM-5109
  - THM-5085
  - THM-5130
  - THM-5129
  - THM-5124
  - THM-4748
  - THM-4879
  - THM-4875
  - THM-4810
  - THM-1370
  - THM-3121
  - HYP-9662
  - HYP-9675
  - THM-5139
---

# THM-5145 - dyadic neighbourhoods and two-block families of path counts

Source: [the note](../../05-knowledge/results/dyadic_cones_two_block_20261111/README.md), sections 2-5.
Independent audit: PASS WITH CORRECTIONS
([report](../../05-knowledge/results/dyadic_cones_two_block_20261111/audit/AUDIT.md); the corrections are listed
in section 10 of the note).

## Statement

`H(T)` is the number of Hamiltonian paths of a tournament `T`. For a rooted tournament `(L, v)`, `B = L - v`,
`O = N+(v)`: offset `c0 = H(B) - sum_{pairs of O} W_ef`, response `M = sum_{y in O} U_B(y)` (THM-5126 (A)); type
`(c, r)` and exponent `j` when `M = c 2^j`, `c` odd, `c0 = r` (so `O` is not empty).

1. **Dyadic neighbourhoods.** If `N` is odd, `N` is not 7 or 21, and `-791 <= N - 2^e <= 647` for some `e >= 1`, then
   `N` is a path count.
2. **Two blocks.** If `c <= 31` is odd, `r` is odd with `|r| <= 275`, `e >= 1`, and `N = c 2^e + r >= 1` is not 7 or 21,
   then `N` is a path count. The same for `c <= 15`, `|r| <= 295`. (The columns `r = 1` and `r = -1` are KNOWN, for
   every odd `c <= 4097`: THM-5126, THM-5130 (E).)
3. **Weight three, rows.** For `1 <= a <= 10` and `b > a`, `1 + 2^a + 2^b` is a path count except 7 and 21, and so
   is `2^b + 2^a - 1`. (For `b - a <= 12` these lie on the known lines and anti-seed lines `c = 2^(b-a) + 1`; the
   addition is `b - a >= 13`.)
4. **One back arc.** For the transitive tournament on `0..m-1` with one back arc `j -> i`, `k = j - i - 1 >= 1`,
   `s = m - 1 - j`: `H(B) = 1 + 2^k` (classical), `W_ij = 2^k`, `U_B(y) = 2^(q+s+1)(3^k + 2^(k+1) - 1)` for `y < i`
   (`q = i - y - 1`), `U_B(i) = 2^s(2^(k+1) - 1)`, `U_B(y) = 2^s(2^(b+1) + 2^(a+1) 3^b - [a = 0])` for `i < y < j`
   (`a = y - i - 1`, `b = j - y - 1`), `U_B(j) = 2^s(2^k + 1)`, `U_B(y) = 2^(m-1-y)(2^k + 1)` for `y > j`. With `O = {y}`
   the response is a power of two only for `k = 2`, `y < i` (offset 5) and `k = 3`, `a = b = 1` (offset 9); these two
   cones are Session236's seeds.
5. **Templates.** A unit pair template is a tournament `K` with two vertices `x`, `y` such that
   `H(K[TT_a at x, TT_b at y]) = 1 + alpha 2^a + beta 2^b + gamma 2^(a+b)` for all `a, b >= 1`. There are 52 and 200
   distinct unit polynomials at orders 7 and 8 (up to the exchange of the parameters). No template of the form
   `1 + 2^j v + 2^i u v` exists with `H(K) < 75` at any order, nor at order 9 with `H(K) <= 2113`. (KNOWN: among the
   tournaments of order at most 8 no pair template has `gamma = 0` and `alpha beta != 0`: kind-pasteur S21, THM-5124;
   the same holds for the 39,451 strong tournaments of order 9 examined here.)

## Proof

Sections 2-5 of the note.
