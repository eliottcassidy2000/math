---
id: THM-4942
title: "First creators at order 11 and the generation of every tournament from them: with the known quadratic attachment form (HYP-2879, THM-4097, THM-4114) run exhaustively over all order-10 classes, every value new at order 11 is created by a one-vertex extension of an order-10 first creator, there is no hypo analogue at order 11, and every first creator of order <= 10 is fertile; every class of order 3-11 is a first creator, a one-vertex extension of a first creator, or obtained from strictly smaller tournaments by an ordinal sum, a substitution, or an ordinal sum with at most two arcs reversed. Through order 10 this is because the connectivity floor f_3(n) exceeds max H(n - 1) (with THM-4811 (B): for n <= 13, lambda >= 3 iff all scores lie in [3, n - 4]). At order 11 the floor inequality fails (f_3(11) = 11,475 in S(10)), but the description survives: the seven non-first-creator classes with lambda >= 3 are substitutions, the floor being C3[R5, R5, 1]. Strong fibers are finite antichains; among all odd h in S(11) (every odd h <= 80,451 except 7 and 21) exactly eleven are born at their Busch bound, 63 the largest"
status: >
  PROVED (elementary):
  - all coefficients of the attachment form are at most max H(n + 1) (the bound on M_ab is the audit's, m1);
  - strong fibers (odd h >= 3) are antichains under induced embedding, of order at most B(h) (THM-1900 (B), the
    insertion lemma, and Busch / Camion for finiteness);
  - the corollary: if f_3(n) > max H(n - 1), every non-first-creator of order n is at most two arc reversals from an
    ordinal sum, and for n <= 13 every tournament with all scores in [3, n - 4] is a first creator.
  FINITE-EXACT:
  - order 11 (pre-registered, af0112b596): R1-R5 all CONFIRMED.
    - Each of the 36,248 values new at order 11 extends an order-10 first creator, so 63 remains the only orphan
      value.
    - There are 0 hypo pairs, and no order-10 first creator is sterile.
    - Every order-10 first creator has at least 564 of 1,024 patterns fertile.
    - The only-first-creator share is 0.515.
    - The extensions reproduce klein's S(11) exactly.
    - The independent audit reran this exhaustively and checked all 36,248 witnesses by direct DP.
  - class-level genealogy at order 11, settled by the independent audit (not pre-registered): every order-11 first
    creator has a first-creator deletion. There are 0 orphans among the 506,912,096 first-creator extensions of
    the 1,085,824 order-10 non-first-creators. So through order 11 every first creator except the two creators of
    63 lies on a Moon chain of first creators from the 3-cycle;
  - fertility at orders 3-10: all first creators fertile; at least half the patterns fertile except the two creators
    of 63;
  - the generation census at orders 3-10. Pre-registered R6 and R9 were CONFIRMED and R7 and R8 FAILED: the pendant
    refinement misses 344 order-10 classes, all with lambda <= 2. The description holds through order 11: at order 11
    the 7 non-first-creator classes with lambda >= 3 are all substitutions;
  - the connectivity floors f_k(n) for n <= 10 (distinct from HYP-9472's two-cycle floor f_2), and lambda >= 3
    implies first creator at orders <= 10;
  - order 11 (pre-registered, c201ec7b34): R10 and R11 CONFIRMED. f_3(11) = 11,475. The seven 3-arc-strong values
    below 15,745 include four in S(10): 11,475, 15,213, 15,273 and 15,561;
  - the eleven values born at their Busch bound. This is exact on all of S(11), so for every odd h <= 80,451; it is
    conditional on HYP-9481 (a) for the 3,635 values of [80,453, 95,093] outside S(11).
  ANALOGY (typed): the dictionary between the 7/21 proof objects and fundamental constructions, in section 5 of
  PROOF.md.
  KNOWN and credited:
  - HYP-2879, THM-4097 (sections 2-3) and THM-4114 (section 6): the attachment form itself (start, end and slot counts,
    the subset-DP slot matrix, the cut-field identity, submodular directed-cut quadratic); THM-4118 (ear lattice). The
    pre-audit text presented it as new (MISTAKE-737);
  - THM-4811 (B): for n <= 13, lambda >= 3 iff min degree >= 3 (lambda = delta^0 whenever n <= 4 delta^0 + 1), and
    the class counts by lambda through order 10;
  - THM-1370, THM-4744, Moon and Busch (the 7/21 omission and the floor);
  - THM-029/343, THM-115 and Grinberg-Stanley (the OCF route);
  - THM-1900 (B) (neutral insertions);
  - THM-012b and THM-4094 (the insertion decomposition);
  - THM-030 and THM-050 (the split products behind W);
  - Session189 (degree-two responses);
  - klein's census, THM-4889 (S(11), n(h), W(h), B(h), HYP-9481);
  - THM-4869, HYP-9459 and HYP-9472 (cycle floors under connectivity);
  - THM-4922 and HYP-9512 (kind-pasteur S16);
  - Session224;
  - Sessions 225-229 (graft and two-path-cover transfer);
  - Session230 (concurrent, same prompt; the THM-4097 pointer and the minimal-basis remark) and opus-2026-11-08-S86
    (reserved the same prompt's IDs at 19:20:14).
  NOT CLAIMED: anything beyond order 11; Chan-Pak (OPEN).
  AUDIT: independent, PASS WITH CORRECTIONS (4 majors including the THM-4097 priority, 18 minors; MISTAKE-737). Every count, floor, witness and
  outcome was reproduced with independent code at orders 3-11.
session: kind-pasteur-2026-11-08-S17
date: 2026-11-08
source: 05-knowledge/results/first_creator_generation_20261108/PROOF.md
related:
  - HYP-9532
  - THM-4922
  - HYP-9512
  - THM-1370
  - THM-1900
  - THM-029
  - THM-115
  - THM-4744
  - THM-4811
  - THM-030
  - THM-050
  - THM-4889
  - THM-4869
  - HYP-9481
  - THM-4097
  - THM-4114
  - HYP-2879
  - 05-knowledge/results/session230_creator_grammar_20261109/README.md
  - 05-knowledge/results/path_spectrum_two_angles_metagraph_20261106/README.md
  - 05-knowledge/results/session229_combined_constructions_20261109/README.md
script: 05-knowledge/results/first_creator_generation_20261108/ext.c; 05-knowledge/results/first_creator_generation_20261108/gen.c; 05-knowledge/results/first_creator_generation_20261108/lam.c; 05-knowledge/results/first_creator_generation_20261108/lam11.c; 05-knowledge/results/first_creator_generation_20261108/audit/aud_o11.c
---

# THM-4942 - first creators at order 11 and the generation of every tournament from them

Proofs and data: [PROOF.md](../../05-knowledge/results/first_creator_generation_20261108/PROOF.md).
Audit: [AUDIT.md](../../05-knowledge/results/first_creator_generation_20261108/audit/AUDIT.md).

A **first creator** is a class of order n with H not in S(n - 1). Arc-strong connectivity lambda(T) is the fewest arc
reversals that make T non-strong. f_k(n) is the least H over n-tournaments with lambda >= k; it is not HYP-9472's
two-cycle floor.

1. **(KNOWN: HYP-2879, THM-4097 (4)-(5), THM-4114 section 6) The attachment form.** For a new vertex v with
   x_a = [a -> v],

       H(T + v_x) = H(T) + sum_a x_a L_a - sum_{a<b} x_a x_b M_ab.

   Here M_ab = W[a][b] + W[b][a] and L_a = f[V][a] - g[V][a] + sum_b W[a][b]. W[a][b] counts the pairs (path ending
   at a, path starting at b) that partition V. M is unsigned, unlike the transfer matrix of THM-030.
2. **(FINITE-EXACT, pre-registered) Order 11.**
   - Every value new at order 11 is created by a one-vertex extension of an order-10 first creator.
   - There is no hypo analogue.
   - Every order-10 first creator is fertile, with at least half its patterns fertile.
   - (audit) Every order-11 first creator has a first-creator deletion: class-level genealogy through order 11.
3. **(FINITE-EXACT) Generation, orders 3-11.** Every class is one of the following:
   - a first creator;
   - a one-vertex extension of a first creator;
   - obtained from strictly smaller tournaments by an ordinal sum, a substitution, or an ordinal sum with at most two
     arcs reversed (at most one through order 8).
4. **(KNOWN, PROVED and FINITE-EXACT) Connectivity floors.**
   - For n <= 13, lambda >= 3 is equivalent to all scores lying in [3, n - 4] (THM-4811 (B)).
   - For n <= 10, f_3(n) > max H(n - 1): the values are 171 > 45, 621 > 189, 2,055 > 661 and 5,753 > 3,357. So every
     such tournament is a first creator.
   - At order 11, f_3(11) = 11,475 < 15,745 = max H(10). Exactly 7 classes of 3-arc-strong 11-tournaments are not
     first creators. All 7 are substitutions, and 6 also extend a first creator.
   - The floor is C3[R5, R5, 1]: the 3-cycle, the regular first creator R5 of 15 twice, and a singleton. So item 3
     still holds at order 11.
5. **(PROVED) Strong fibers** of every odd value h >= 3 are finite antichains under induced embedding, of order at
   most B(h). Through order 11 the first creators of h are the members of order n(h) (THM-4922 (1)).
6. **(FINITE-EXACT) Values born at their Busch bound.** Among all h in S(11), so every odd h <= 80,451, they are
   exactly 3, 5, 9, 11, 13, 17, 19, 23, 35, 39 and 63.
