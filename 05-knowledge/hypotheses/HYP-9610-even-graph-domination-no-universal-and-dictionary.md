---
id: HYP-9610
title: "(A) No-universal-vertex form of THM-5040: for n >= 5, the dominating-set counts of even graphs on n vertices with no universal vertex are exactly the odd numbers up to M0(n) = 2^n - n - 1 (n even) or 2^n - n - 4 (n odd), except 3, 7, 13 and M0(n) - 2. (B) The binary/ternary dictionary of exceptional sets (tournaments {7, 21}, base-3 packings {10, 31}, even graphs {3, 7, 13}, odd graphs {1, 5, 37}, semicomplete digraphs and all graphs none), typed ANALOGY"
status: >
  (A) OPEN for general n.
  - PROVED for all n: the bound M0(n) and its attainment; for even n, the exclusion of M0(n) - 2.
  - FINITE-EXACT for 5 <= n <= 9 (census).
  - Containment verified at n = 10 and 11: the THM-5040 closure, plus a random search realising all 38 values that
    closure missed (e0search.c, independent automorphism code).
  (B) ANALOGY: no structural map is claimed.
session: opus-2026-11-09-S89
date: 2026-11-09
related: [THM-5040, THM-5055, THM-5056, THM-5041, THM-5042, THM-4996, THM-4950, THM-5071, THM-5070, HYP-9575]
script: 04-computation/experiments/egcp_20261109/
---

# HYP-9610 - Even graphs without a universal vertex; the binary/ternary dictionary

## (A) Values without a universal vertex

Let E0(n) be the set of d(G) over even graphs G on exactly n vertices with no universal vertex. These are the
co-connected or padded realisers, i.e. the possible cone inputs.

**Conjecture.** For n >= 5, E0(n) = {odd k <= M0(n)} minus {3, 7, 13, M0(n) - 2}, where M0(n) = 2^n - n - 1 (n even)
and 2^n - n - 4 (n odd).

**Proved for all n.** Write H = comp(G). H has no isolated vertex, so D(H) >= n + 1; THM-5040 Lemma 0 then gives
d <= 2^n - n - 1.
- **Odd n: D(H) is at least n + 4.**
  - D = n + 1 needs a perfect matching, which is impossible on an odd number of vertices.
  - D = n + 2 needs P_3 plus a matching, which is not co-even.
  - D = n + 3 needs P_4, C_4 or 2P_3 plus a matching. Each has even order.
  - Hence d <= 2^n - n - 4. Equality holds for comp(P_5 u sK_2), which is even since P_5 is co-even (5 = 1 mod 4).
- **Even n.** M0(n) is attained by the cocktail-party graph. M0(n) - 2 corresponds to D = n + 3, which forces P_4,
  C_4 or 2P_3 plus a matching; none is co-even (THM-5040 Lemma 3).

**Evidence.**
- n <= 9: exact from the census (THM-4996). The closure of THM-5040 (cover.py) misses only about 2n values of E0(n)
  per order, all near 2^(n-1) + 2^k.
- n = 10 and 11: random search (e0search.c) realised every one of those misses: 18 at n = 10 and 20 at n = 11. Each
  found graph was checked to be even, to have no universal vertex, and to have the stated d
  (e0search_n10.out, e0search_n11.out).
- **Open.** The exclusion of M0(n) - 2 for odd n >= 11, which corresponds to D = n + 6. Containment for n >= 12.

## (B) Binary/ternary dictionary (typed ANALOGY)

| Object | Count at order n | Parity theorem | Exceptional set | Source |
| --- | --- | --- | --- | --- |
| Tournaments (binary triangle: lines of the metagraph) | 2^C(n,2) labelled | H odd (Rédei) | {7, 21}; complete iff Chan-Pak | THM-830, Session234 |
| Semicomplete digraphs (ternary triangle; covering tiling pairs) | 3^C(n,2), times 2^n with loops | H = packings of defect cycles (mod 2) | none: every positive integer | THM-5071 (3), (4a) |
| One reciprocal pair, no missing pair | | | {1, 4, 7} | THM-5071 (4b) |
| Base-3 odd-cycle packings I(Omega, 3) | | = 1 (mod 3) | {10, 31}: the same alpha-sums 3, 10 as 7 = 1 + 2*3 and 21 = 1 + 2*10 | THM-4950, THM-5070 |
| Even graphs (equinumerous with tournaments, RPGFD) | A334335 | d odd (BCS) | {3, 7, 13}; one extra vertex on the top bands | THM-5040 |
| Odd graphs | | d odd | {1, 5, 37}; no extra vertex | THM-5056 (first), THM-5042 |
| All graphs | | d odd | none, on ceil(log2(k + 1)) vertices | Brouwer-Csorba-Schrijver 2009, Prop. 1.2; THM-4996 (1) |

**Reading.**
- Each binary object has a small exceptional set, and enlarging its alphabet removes it:
  - tournaments to semicomplete digraphs (cells from 2 states to 3);
  - even graphs to all graphs (dropping the parity constraint).
- In both even/odd pairs the exceptions are the counts of the few smallest configurations that a twin or symmetry
  obstruction forbids:
  - adjacent twins make K_2, K_3 and K_4 - e odd;
  - non-adjacent twins make P_3 and the empty graphs even.
- **The owner's lenses.**
  - 2 and 3 enter as cell alphabets (binary tiles, ternary semicomplete tiles) and as packing bases (H = I(Omega, 2),
    H_3 = I(Omega, 3)).
  - In Collatz they are the multiplier/divisor pair (THM-4532 and the hexagonal geometry of session140).
  - **Typed ANALOGY:** no map between these settings is claimed, and none is known (THM-4531, THM-4996 (4)).
