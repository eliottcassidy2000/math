---
id: THM-4942
title: "First creators at order 11 and the generation of every tournament from them: H(T + v) is an explicit quadratic form in the attachment pattern (split-path coefficients); every value new at order 11 is created by a one-vertex extension of an order-10 first creator, there is no hypo analogue at order 11, and every first creator of order <= 10 is fertile; through order 10 every tournament is a first creator, a first creator plus a vertex, an ordinal sum, a substitution, or an ordinal sum re-connected by at most two arcs, because the connectivity floor f_3(n) exceeds max H(n - 1) (for n <= 13, lambda >= 3 iff all scores lie in [3, n - 4]); this fails at order 11 (3-arc-strong 11-tournaments with H = 15,213 in S(10)); strong fibers are finite antichains; exactly eleven odd values up to 95,095 are born at their Busch bound, 63 the largest"
status: >
  PROVED (elementary):
  - the attachment form (section 1 of PROOF.md), with W, L, M given by split Hamiltonian paths;
  - for n <= 13, arc-strong connectivity lambda >= 3 if and only if every vertex has in- and out-degree >= 3;
  - strong fibers are antichains under induced embedding (from THM-1900 (B) and the insertion lemma);
  - the corollary: if f_3(n) > max H(n - 1) with n <= 13, every non-first-creator of order n is at most two arc
    reversals from an ordinal sum.
  FINITE-EXACT:
  - order 11 (pre-registered, af0112b596): R1-R5 all CONFIRMED. Each of the 36,248 values new at order 11 extends an
    order-10 first creator, so 63 remains the only orphan value. There are 0 hypo pairs and no sterile order-10 first
    creator; every order-10 first creator has at least 564 of 1,024 patterns fertile; the only-first-creator share
    is 0.515. The extensions reproduce klein's S(11) exactly.
  - fertility at orders 3-10 (all first creators fertile; at least half the patterns fertile except the two creators
    of 63);
  - the generation census at orders 3-10, where pre-registered R6 and R9 were CONFIRMED and R7 and R8 FAILED: the
    pendant refinement misses 344 order-10 classes, all with lambda <= 2;
  - the connectivity floors f_k(n) for n <= 10, and lambda >= 3 implies first creator at orders <= 10;
  - order 11 (pre-registered, c201ec7b34): R10 and R11 CONFIRMED. f_3(11) = 11,475. The seven 3-arc-strong values
    below 15,745 include four in S(10): 11,475, 15,213, 15,273 and 15,561.
  - the eleven values born at their Busch bound through 95,095.
  ANALOGY (typed): the dictionary between the 7/21 proof objects and fundamental constructions, in section 5 of
  PROOF.md.
  KNOWN and credited:
  - THM-1370, Moon and Busch (the 7/21 omission, the floor);
  - THM-029/343 and THM-115 (the OCF route);
  - THM-1900 (B) (neutral insertions);
  - THM-012b and THM-4094 (the insertion decomposition);
  - Session189 (degree-two responses);
  - klein's census (S(11), n(h), W(h), B(h), HYP-9481);
  - THM-4922 and HYP-9512 (kind-pasteur S16);
  - Session224;
  - Sessions 225-229 (graft and two-path-cover transfer).
  NOT CLAIMED: class-level genealogy at order 11; anything beyond order 11; Chan-Pak (OPEN).
  AUDIT: see audit/AUDIT.md.
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
  - HYP-9481
  - 05-knowledge/results/path_spectrum_two_angles_metagraph_20261106/README.md
  - 05-knowledge/results/session229_combined_constructions_20261109/README.md
script: 05-knowledge/results/first_creator_generation_20261108/ext.c; 05-knowledge/results/first_creator_generation_20261108/gen.c; 05-knowledge/results/first_creator_generation_20261108/lam.c; 05-knowledge/results/first_creator_generation_20261108/lam11.c
---

# THM-4942 - first creators at order 11 and the generation of every tournament from them

Proofs and data: [PROOF.md](../../05-knowledge/results/first_creator_generation_20261108/PROOF.md).
Audit: [AUDIT.md](../../05-knowledge/results/first_creator_generation_20261108/audit/AUDIT.md).

A **first creator** is a class of order n with H not in S(n - 1). Arc-strong connectivity lambda(T) is the fewest arc
reversals that make T non-strong. f_k(n) is the least H over n-tournaments with lambda >= k.

1. **(PROVED) The attachment form.** For a new vertex v with x_a = [a -> v],

       H(T + v_x) = H(T) + sum_a x_a L_a - sum_{a<b} x_a x_b M_ab.

   Here M_ab = W[a][b] + W[b][a] and L_a = f[V][a] - g[V][a] + sum_b W[a][b]. W[a][b] counts the pairs (path ending
   at a, path starting at b) that partition V.
2. **(FINITE-EXACT, pre-registered) Order 11.**
   - Every value new at order 11 is created by a one-vertex extension of an order-10 first creator.
   - There is no hypo analogue.
   - Every order-10 first creator is fertile, with at least half its patterns fertile.
3. **(FINITE-EXACT) Generation, orders 3-10.** Every tournament is one of the following:
   - a first creator;
   - a first creator plus a vertex;
   - an ordinal sum;
   - a substitution;
   - an ordinal sum re-connected by at most two arcs (at most one through order 8).
4. **(PROVED and FINITE-EXACT) Connectivity floors.**
   - For n <= 13, lambda >= 3 is equivalent to all scores lying in [3, n - 4].
   - For n <= 10, f_3(n) > max H(n - 1): the values are 171 > 45, 621 > 189, 2,055 > 661 and 5,753 > 3,357. So every
     such tournament is a first creator.
   - At order 11, f_3(11) = 11,475 < 15,745 = max H(10). There are 3-arc-strong 11-tournaments that are not first
     creators.
5. **(PROVED) Strong fibers** of every value are finite antichains under induced embedding, of order at most B(h).
6. **(FINITE-EXACT) Values born at their Busch bound** through 95,095: exactly 3, 5, 9, 11, 13, 17, 19, 23, 35, 39 and
   63.
