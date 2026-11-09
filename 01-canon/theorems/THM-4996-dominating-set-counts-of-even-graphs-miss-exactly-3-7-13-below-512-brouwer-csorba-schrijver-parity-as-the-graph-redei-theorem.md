---
id: THM-4996
title: "Dominating-set counts of even graphs (von Brömssen's notion: no automorphism reverses an odd number of edges; equinumerous with tournaments by RPGFD 2023). (1) PROVED: the number d(G) of dominating sets (always odd, Brouwer-Csorba-Schrijver 2009) takes every odd value over all graphs, on ceil(log2(k + 1)) vertices (threshold graphs). (2) PROVED: a parity calculus for even graphs: adjacent twins are odd; unions, cones and joins are even under stated conditions; trees are odd iff some automorphism flips the central edge. (3) FINITE-EXACT: an exhaustive census of all connected graphs with d <= 511 (474,101 classes) shows that on even graphs d takes every odd value below 512 except 3, 7 and 13, and on odd graphs every odd value except 1, 5 and 37. The only realisers of 3, 7, 13 are K_2, K_3, K_4 - e (with isolated vertices), and 21 = d(K_1 v (P_3 u K_1)) is even. (4) The OCF statistic I(G, 2) misses 5, 7, 13, 15, 21, 23, 31 on even graphs below 32"
status: >
  PROVED + FINITE-EXACT + COMPUTATIONAL EVIDENCE. Independent audit P1 PASS WITH CORRECTIONS. The census was reproduced
  with independent canonical forms and full automorphism enumeration, all 256 per-value rows matching. Corrections
  applied (MISTAKE-746): citations (Brouwer-Csorba-Schrijver; von Brömssen; Alikhani-Peng) and wording.
  - (1) PROVED, elementary: the cone recursion d(K_1 v G) = 2^|G| + d(G). It is the x = 1 case of the
    domination-polynomial join theorem of Alikhani-Peng. No prior art was found for completeness or conciseness of d
    on the odd numbers. Correction (MISTAKE-747, opus-2026-11-09-S89, after its audit): both completeness and the
    conciseness ceil(log2(k + 1)) are A. E. Brouwer, P. Csorba, A. Schrijver, The number of dominating sets of a
    finite graph is odd, preprint 2009, Proposition 1.2 (same pad/cone induction as Lemma B1); restated as
    Alikhani-Fatehi-Klavzar arXiv:1512.07514, Prop. 2.1(iii). (Alikhani, Graphs Combin. 29 (2013) concerns D(G, -1).)
  - (2) PROVED: PROOF.md Part B, Lemma B2; checked by audit P1 by brute force, including all 436 trees up to 11
    vertices.
  - (3) FINITE-EXACT: dom_even.c.
    - Exhaustiveness rests on strict monotonicity under adding a joined vertex and on connected prefix orders
      (Lemma B3).
    - Every d is recomputed. The class counts agree with A001349 (n <= 9), the Euler transform with A000088, and the
      Euler transform of the connected-even counts with A000568 through n = 9: a count-by-count re-check of RPGFD's
      theorem (von Brömssen already verified the counts to n = 10; OEIS A334335).
    - Audit P1 reproduced everything independently (dom_census.c).
  - (4) FINITE-EXACT: stat_even.c with I(G, 2), exhaustive below 32.
  - Coverage beyond 512 is COMPUTATIONAL EVIDENCE only (dom_cover.c). Explicit even constructions reach every odd
    value <= 2^22 except 3, 7, 13 and 12,326 values clustered near 2^m - 2^k (including 2^m), which are not shown
    unattainable.
  KNOWN and credited:
  - A. E. Brouwer, P. Csorba, A. Schrijver, The number of dominating sets of a finite graph is odd, preprint (2009)
    (as cited by Heinrich-Tittmann, arXiv:1701.03453);
  - P. von Brömssen (né Andersson): even and odd graphs and the equinumerosity conjecture (OEIS A334335);
  - G. Royle, C. Praeger, S. Glasby, S. Freedman, A. Devillers, Tournaments and even graphs are equinumerous,
    J. Algebraic Combin. 57 (2023) 515-524, arXiv:2204.01947 (proof of the conjecture);
  - S. Alikhani, Y.-H. Peng, arXiv:0905.2251, Theorem 2 (the join formula of the domination polynomial); Arocha-Llano
    2000 (the domination polynomial);
  - THM-4531, THM-4721, THM-4722, HYP-9340 (repository even-graph lane);
  - S. H. Chan, I. Pak, Theoret. Comput. Sci. 1015 (2024) 114776, arXiv:2308.10214: Section 5.1 (complete, almost
    complete), Definition 5.7 (concise), Proposition 5.14 and Conjecture 6.2;
  - Session232 (THM-4987-4990): the automorphism-orbit side of 7/21. This theorem is a different, intrinsic
    statistic.
  NOT CLAIMED:
  - the exceptional set above 511 (HYP-9575 (B) conjectured {3, 7, 13}; now PROVED, with sharp orders, by THM-5040
    and independently by THM-5055);
  - any natural correspondence between d on even graphs and H on tournaments. The shared 7 is an analogy.
session: mac-mini-2026-11-09-paleyeven
date: 2026-11-09
related: [THM-4995, HYP-9575, THM-4531, THM-4721, THM-4722, THM-1370, THM-4986, THM-4987]
note: 05-knowledge/results/paleyeven_20261109/README.md
proof: 05-knowledge/results/paleyeven_20261109/PROOF.md
audit: 05-knowledge/results/paleyeven_20261109/AUDIT_P1.md
scripts:
  - 04-computation/experiments/paleyeven_20261109/dom_even.c
  - 04-computation/experiments/paleyeven_20261109/stat_even.c
  - 04-computation/experiments/paleyeven_20261109/dom_cover.c
outputs:
  - dom_K511.err, dom_K511_summary.txt (dom_K511.out is regenerable: ./dom_even 511 16), dom_K255.err
  - i2_K31.out, dom_cover_K511_22_basic.out, dom_cover_K511_22_genjoin.out
---

# THM-4996 - dominating-set counts of even graphs

## Statement

d(G) is the number of dominating sets of G; it is always odd (Brouwer-Csorba-Schrijver). A graph is **even** (von
Brömssen) if every automorphism reverses an even number of edges relative to a fixed reference orientation, and
**odd** otherwise. Even graphs and tournaments are equinumerous (RPGFD 2023).

1. **(PROVED)** Every odd k is d(G) for a threshold graph G on ceil(log2(k + 1)) vertices. So d is complete on the
   odd numbers and concise in Chan-Pak's sense.
2. **(PROVED)** Parity calculus:
   - swapping adjacent twins is odd, and swapping non-adjacent twins is even;
   - G is even iff each component is;
   - K_1 v G is even iff G is, when G has >= 2 vertices and no universal vertex;
   - A v G is even iff A and G are, when both complements are connected and A is not isomorphic to G;
   - sgn_G(g) sgn_{complement(G)}(g) = sgn(g);
   - a tree is odd iff it is bicentral with an automorphism exchanging the centres.
3. **(FINITE-EXACT)** Census of every connected graph with d <= 511:
   - order counts 1, 1, 2, 6, 21, 112, 853, 11,117, 261,080, 199,053, 1,855 (n = 1..11), and none at 12;
   - connected even counts 1, 0, 1, 2, 8, 43, 398, 6,413, 184,596 (n <= 9; A334335).

   The set {d(G) : G even} below 512 is all odd numbers except 3, 7, 13.
   - Connected even graphs miss exactly 3, 7, 13 and 25 = d(P_3 u P_3).
   - Up to isolated vertices, the graphs with d in {3, 7, 13} are exactly K_2, K_3 and K_4 - e, and each has adjacent
     twins.
   - 21 = d(K_1 v (P_3 u K_1)), which is even (|Aut| = 2, generated by a swap of non-adjacent twins).
   - Odd graphs (some component odd) take every odd value below 512 except 1, 5 and 37. The four connected realisers
     of the prime 37 are all even.
4. **(FINITE-EXACT)** I(G, 2) = sum over independent sets of 2^|S| (so H(T) = I(Omega(T), 2) by OCF). Over all
   graphs it takes every odd value. On even graphs below 32 it takes exactly 1, 3, 9, 11, 17, 19, 25, 27, 29.

## Reading

Tournaments and even graphs are equinumerous, with no natural bijection. Their intrinsic odd counts compare as
follows.

| | Rédei's H on tournaments | d on even graphs |
|---|---|---|
| multiplicative over | ordinal sums | disjoint unions |
| exceptional set (conjectured) | {7, 21} (Chan-Pak Conj. 6.2) | {3, 7, 13} (HYP-9575 (B)) |
| exceptional set (known exactly) | below 80,453 (klein's order-11 census) | below 512 |
| exclusion mechanism | floors of strong components | adjacent twins (odd transpositions) of the only realisers |

## Proof

PROOF.md Part B (Lemmas B1-B3, Theorem B4, Propositions B5-B6).
