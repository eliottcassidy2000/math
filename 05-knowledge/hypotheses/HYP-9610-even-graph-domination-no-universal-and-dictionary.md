---
id: HYP-9610
title: "(A) No-universal-vertex form of THM-5040: for n >= 5, the dominating-set counts of even graphs on n vertices with no universal vertex are exactly the odd numbers up to M0(n) = 2^n - n - 1 (n even) or 2^n - n - 4 (n odd), except 3, 7, 13 and M0(n) - 2. (B) The binary/ternary dictionary of exceptional sets (tournaments {7, 21}, base-3 packings {10, 31}, even graphs {3, 7, 13}, odd graphs {1, 5, 37}, semicomplete digraphs and all graphs none), typed ANALOGY"
status: >
  (A) The exclusion side is PROVED for every n. Containment is FINITE-EXACT for n <= 13 and OPEN for n >= 14.
  - Exclusions: the bound M0(n), its attainment, and the exclusion of M0(n) - 2. Audit2 proved the odd-n exclusion
    (Note 1, excess additivity).
  - Containment: the census for n <= 9; audit2's direct search for n = 10-13, with every witness re-verified
    (a2_e0full_n10..13).
  - S89's own n = 10, 11 search had an array overflow (found[16] used with 18 and 20 targets), so six "realised"
    lines had no graph. The claim it supported is true by audit2's witnesses (MISTAKE-747).
  (B) ANALOGY: no structural map is claimed.
session: opus-2026-11-09-S89
date: 2026-11-09
related: [THM-5040, THM-5055, THM-5056, THM-5041, THM-5042, THM-4996, THM-4950, THM-4940, THM-1370, THM-5071, THM-5070, HYP-9530, HYP-9575]
script: 04-computation/experiments/egcp_20261109/
---

# HYP-9610 - Even graphs without a universal vertex; the binary/ternary dictionary

Audit: [AUDIT2](../results/egcp_20261109/audit2/AUDIT2.md).

## (A) Values without a universal vertex

Let E0(n) be the set of d(G) over even graphs G on exactly n vertices with no universal vertex. These are exactly the
inputs of an even cone (B2.3).

**Conjecture.** For n >= 5, E0(n) = {odd k <= M0(n)} minus {3, 7, 13, M0(n) - 2}, where M0(n) = 2^n - n - 1 (n even)
and 2^n - n - 4 (n odd).

**Exclusion side, PROVED for every n.**
- **Setting.** Write H = comp(G). H has no isolated vertex, so D(H) = 1 + n + e(H), where e(H) counts the sets of
  size >= 2 inside some open neighbourhood. G is even iff H is co-even (B2.5).
- **Facts about e (audit2, Note 1).**
  - e is additive over components (Brown-Nowakowski), and e(K_2) = 0.
  - The connected graphs with e <= 5 form 12 classes. Only three are co-even: K_3 and P_5 (e = 3), and P_6 (e = 4).
  - A disjoint union is co-even only if each component is.
- **Consequence.** e(H) is never 1, 2 or 5, so d is never 2^n - n - 2, 2^n - n - 3 or 2^n - n - 6.
  - Since d = 2^n - n - 1 - e is odd, n odd forces e odd, so e >= 3. This gives d <= 2^n - n - 4 = M0(n), attained by
    comp(P_5 u sK_2) and comp(K_3 u sK_2). M0(n) - 2 would need e = 5, which is impossible.
  - For n even, d <= 2^n - n - 1 = M0(n) (cocktail party, e = 0), and M0(n) - 2 would need e = 2, which is
    impossible.

**Containment, FINITE-EXACT for n <= 13.**
- n <= 9: the census (THM-4996; base_sets.json).
- n = 10-13: audit2's direct search (a2_e0hunt.c with a2_e0full.py) realised every conjectured value: 503, 1,013,
  2,038 and 4,084 values. Every witness was re-verified for d, evenness and the absence of a universal vertex.
- The THM-5040 closure (cover.py) alone misses about 2n values per order, all near 2^(n-1) + 2^k. These misses are
  limits of the closure, not gaps.

**Open.** Containment for n >= 14.

## (B) Binary/ternary dictionary (typed ANALOGY)

| Object | Count at order n | Parity theorem | Exceptional set | Source |
| --- | --- | --- | --- | --- |
| Tournaments (binary triangle: lines of the metagraph) | 2^C(n,2) labelled; A000568 unlabelled | H odd (Rédei) | 7 and 21 excluded at every order; complete iff Chan-Pak | THM-1370; THM-830, Session234 |
| Semicomplete digraphs (ternary triangle; covering tiling pairs) | 3^C(n,2), times 2^n with loops | H = packings of defect cycles (mod 2) | none: every positive integer | THM-5071 (3), (4a) |
| One reciprocal pair, no missing pair | | | 1, 4, 7 never attained; complete below 198 | THM-5071 (4b), (4c) |
| Base-3 odd-cycle packings I(Omega, 3) | | = 1 (mod 3) | 10 and 31 excluded; conjecturally complete (checked to 10^6) | THM-4940; HYP-9530 (A); THM-5070 (4) |
| Even graphs (equinumerous with tournaments, RPGFD) | A000568 unlabelled (connected: A334335) | d odd (BCS) | {3, 7, 13}; one extra vertex on the top bands | THM-5040 (THM-5055) |
| Odd graphs | | d odd | {1, 5, 37}; no extra vertex | THM-5056 (first), THM-5042 |
| All graphs | | d odd | none, on ceil(log2(k + 1)) vertices | Brouwer-Csorba-Schrijver 2009, Prop. 1.2; THM-4996 (1) |

**Reading.**
- Each binary object has a small exceptional set, and enlarging its alphabet removes it:
  - tournaments to semicomplete digraphs (cells from 2 states to 3);
  - even graphs to all graphs (dropping the parity constraint).
- **The even/odd split of d.** The small exceptions 3, 7, 13 and 1, 5 are the counts of the few smallest graphs whose
  twins fix the parity: adjacent twins make K_2, K_3 and K_4 - e odd; non-adjacent twins make P_3 and the empty graphs
  even.
- **Not explained by twins.**
  - 37: its four connected realisers have 6 vertices and automorphism groups of orders 4, 1, 4 and 2. One is
    asymmetric (audit2 Note 6).
  - {7, 21} for tournaments: these come from strong-component floors (THM-1370; THM-4996's table).
- **The owner's lenses.**
  - 2 and 3 enter as cell alphabets (binary tiles, ternary semicomplete tiles) and as packing bases (H = I(Omega, 2),
    H_3 = I(Omega, 3)).
  - In Collatz they are the multiplier/divisor pair (THM-4532 and the hexagonal geometry of session140).
  - **Typed ANALOGY:** no map between these settings is claimed, and none is known (THM-4531).
