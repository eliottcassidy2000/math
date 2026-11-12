---
id: THM-5145
title: "Dyadic neighbourhoods and two-block families of Hamiltonian-path counts. A rooted tournament (L, v) with offset c0 = r and response M = c 2^j (c odd) - type (c, r) - puts c 2^e + r into the spectrum for every e >= j + 1 (the tail law, KNOWN). Complete census of the rooted tournaments of order at most 10: every odd r with -791 <= r <= 647 is the offset of a root whose response is a power of two; every odd c <= 31 and odd r with |r| <= 275 is a type. Hence (1) every odd N with -791 <= N - 2^e <= 647 for some e >= 1 is a path count, except 7 and 21; (2) every number c 2^e + r >= 1 with c <= 31 odd, |r| <= 275 odd, e >= 1 is a path count, except 7 and 21; (3) 1 + 2^a + 2^b is a path count for 1 <= a <= 10 and every b > a, except 7 and 21. (4) For the transitive tournament with one back arc, closed forms of U_B, and the only roots with one out-neighbour and a power-of-two response have offsets 5 and 9. (5) Census of two-parameter blow-up templates to order 8: no independent unit pair, no template 1 + 2^j v + 2^i u v"
status: >
  PROVED: (1)-(3) from a complete FINITE-EXACT census (every base B among the classes of order 7 and 8 and all
  1,761,280 one-vertex extensions of the classes of order 8, every out-set O) together with two inherited results:
  the tail law H(L o 0^t) = c0 + M 2^t (THM-5126 (A)) and the verified head, every odd number up to 2^32 other than
  7 and 21 (THM-4875 (D)); the largest exponent j in the ranges used is 13. (4) by direct counting of covers, the
  2-adic valuation of 3^k - 1 and 3^b + 1 = 2^u only for b = 1. (5) FINITE-EXACT.
  Machine checks: every cone witness (118, 729, 5,177 at root orders 8, 9, 10) recounted literally at tail lengths
  1, 2, 3 by a subset count written without the offset formula (18,072 counts, no failure); 10,179 witnesses of
  types with c = 3..31 at order 9 (30,537 counts); the closed forms of (4) for all m <= 11 (165 cases); templates
  fitted at four points and checked at further points.
  KNOWN, credited: the tail law, rays, unit-offset roots and the complete lines c <= 4097 (THM-5126, HYP-9658;
  kind-pasteur S21; Session216); the cone mechanism, the window bound and the dyadic ports of order <= 8
  (THM-5109, kind-pasteur S20; Session236); fixed-port tails and one negative offset 2^K - 159 (THM-5085); the
  substitution law (THM-3121, THM-5126 (B1)); the verified ranges (THM-4875, THM-4885); 7 and 21 are not path
  counts (THM-1370); strong floors (THM-4810).
  NEW to the repository (not found by the author; the audit is asked to refute this): the complete census of
  cone roots and of types (c, r) to order 10 with offsets of both signs; (1), (2), (3) as stated (rows a <= 10
  of the weight-three statement HYP-9658 (C); numbers beginning with the binary digits 11 reached through
  negative offsets); (4); the census (5) and its two exclusions.
  NOT CLAIMED: the Chan-Pak conjecture (OPEN); positive density (each family has O(log X) members up to X);
  completeness of the template census at order 9 (restricted to a list of values of H).
  Pre-registered (PREREGISTRATION.md of the source packet, written while the order-10 runs were in progress and
  before any output was read): Q1, Q2, Q4 CONFIRMED, Q3 FAILED (offsets 649, 825).
  AUDIT: AUDIT-PENDING (pre-audit checkpoint)
session: mac-mini-2026-11-11-rays
source: 05-knowledge/results/dyadic_cones_two_block_20261111/README.md
related:
  - THM-5126
  - HYP-9658
  - THM-5109
  - THM-5085
  - THM-4875
  - THM-1370
  - THM-3121
  - HYP-9675
  - THM-5139
---

# THM-5145 - dyadic neighbourhoods and two-block families of path counts

Source: [the note](../../05-knowledge/results/dyadic_cones_two_block_20261111/README.md), sections 2-5.
Independent audit: AUDIT-PENDING (pre-audit checkpoint)

## Statement

`H(T)` is the number of Hamiltonian paths of a tournament `T`. For a rooted tournament `(L, v)`, `B = L - v`,
`O = N+(v)`: offset `c0 = H(B) - sum_{pairs of O} W_ef`, response `M = sum_{y in O} U_B(y)` (THM-5126 (A)); type
`(c, r)` and exponent `j` when `M = c 2^j`, `c` odd, `c0 = r`.

1. **Dyadic neighbourhoods.** If `N` is odd, `N` is not 7 or 21, and `-791 <= N - 2^e <= 647` for some `e >= 1`, then
   `N` is a path count.
2. **Two blocks.** If `c <= 31` is odd, `r` is odd with `|r| <= 275`, `e >= 1`, and `N = c 2^e + r >= 1` is not 7 or 21,
   then `N` is a path count. The same for `c <= 15`, `|r| <= 295`.
3. **Weight three, rows.** For `1 <= a <= 10` and `b > a`, `1 + 2^a + 2^b` is a path count except 7 and 21, and so
   is `2^b + 2^a - 1`.
4. **One back arc.** For the transitive tournament on `0..m-1` with one back arc `j -> i`, `k = j - i - 1 >= 1`,
   `s = m - 1 - j`: `H(B) = 1 + 2^k`, `W_ij = 2^k`, `U_B(y) = 2^(q+s+1)(3^k + 2^(k+1) - 1)` for `y < i` (`q = i - y - 1`),
   `U_B(i) = 2^s(2^(k+1) - 1)`, `U_B(y) = 2^s(2^(b+1) + 2^(a+1) 3^b - [a = 0])` for `i < y < j` (`a = y - i - 1`,
   `b = j - y - 1`), `U_B(j) = 2^s(2^k + 1)`, `U_B(y) = 2^(m-1-y)(2^k + 1)` for `y > j`. With `O = {y}` the response is a
   power of two only for `k = 2`, `y < i` (offset 5) and `k = 3`, `a = b = 1` (offset 9).
5. **Templates.** Among the tournaments of order at most 8 (and 27,315 strong ones of order 9 with `H` in a list
   of 38 values) there is no unit pair template `1 + alpha u + beta v` with `alpha beta != 0`, and none of the form
   `1 + 2^j v + 2^i u v`.

## Proof

Sections 2-5 of the note.
