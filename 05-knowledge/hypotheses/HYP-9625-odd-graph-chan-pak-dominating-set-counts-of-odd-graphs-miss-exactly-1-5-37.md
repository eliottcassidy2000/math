---
id: HYP-9625
title: "Odd-graph Chan-Pak: the numbers of dominating sets of odd graphs (some automorphism reverses an odd number of edges) are exactly the odd numbers other than 1, 5 and 37. Equivalent to: every prime p >= 512 is d of an odd graph. The companion of THM-5055 (even graphs: exactly the odd numbers other than 3, 7, 13)"
status: >
  OPEN. EVIDENCE: exhaustive below 512 (THM-4996). PROVED: the reduction to primes, and the inclusion of every odd
  number whose binary expansion starts with 11. The proof route is the THM-5055 machinery.
session: mac-mini-2026-11-09-digraph721
date: 2026-11-09
related: [THM-5055, THM-4996, HYP-9575]
note: 05-knowledge/results/digraph721_20261109/README.md
---

# HYP-9625: odd-graph Chan-Pak

**Conjecture.** O = {d(G) : G odd} consists of all odd numbers except 1, 5 and 37.

**Evidence (FINITE-EXACT, THM-4996).** The census of all connected graphs with d <= 511 shows that, below 512, O misses
exactly 1, 5 and 37:
- 1 is realised only by edgeless graphs;
- 5 only by P_3;
- 37 only by four connected graphs, all even.

Since 37 is prime, every realiser of 37 is one of those plus isolated vertices. So 1, 5 and 37 are excluded for every
order.

## Proved reductions

1. **(Unions)** A graph is odd iff some component is odd. Since d is multiplicative, O * (any odd number) is
   contained in O.
2. **(Join with K_2)** For every graph G on m vertices, K_2 v G is odd, because the two K_2-vertices are adjacent
   twins. Its value is d = 3 * 2^m + d(G). Every odd number below 2^m is d of a threshold graph on m vertices.
   Hence O contains every odd number whose binary expansion starts with 11.
3. **(Cones)** K_1 v H is odd whenever H is odd: with no universal vertex by THM-4996 B2.3, and with one by adjacent
   twins. Its value is 2^|H| + d(H).
4. **(Prime reduction)** The conjecture is equivalent to "every prime p >= 512 is in O".
   - Suppose an odd k other than 1, 5, 37 has a prime factor p other than 5 and 37. Then p is in O (below 512 by the
     census, and above by hypothesis), so k is in O by item 1.
   - Otherwise k = 5^a 37^b.
     - If a >= 2, use 25 (in O).
     - If a = b = 1, use 185 (in O, census).
     - If b >= 2, use 1369 = 37^2. It is in O: pad K_2 u G_115 to 10 vertices, where G_115 is any graph with
       d = 115, so the padded graph has d = 345. Cone over it: 2^10 + 345 = 1369.
5. **(What remains)** By items 2 and 3, a prime p fails to be in O through these routes only if every stage of the
   cone recursion p -> p - 2^(L(p) - 1) has a leading 10 and the chain avoids O. In particular this happens if the
   chain passes through 1, 5 or 37. The first open cases are of the forms 2^m + 1 (Fermat primes, such as 65537),
   2^m + 5 and 2^m + 37 (such as 1061 = 2^10 + 37, which is prime).

## Route

Use the THM-5055 hub moves with a persistent pair of adjacent twins (x, y):
- cones, E_2-joins, padding and PH(G, u, l) with u not in {x, y} keep x and y adjacent twins, hence keep the graph
  odd;
- the value formulas are those of THM-5055;
- so the offset method (R and PH(3, 2) carrying 2^(n-1) + c) should supply the families 2^n + {1, 5, 37} from small
  odd seeds.

An order-restricted exploration (`oddclos.py`; cones, joins with K_2 and E_2, and products by 3 only, with exhaustive
bases on <= 9 vertices) misses at order 10-13 the values 517, 1025, 1061, 2053, 4097, 4133 and 4613. Some of these are
not missing as values: 517 = 11 * 47 is in O. The sharp-order version is therefore more irregular than for even graphs,
and is not conjectured here.
