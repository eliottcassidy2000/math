---
id: THM-5040
title: "Sharp even-graph Chan-Pak (HYP-9575 (B) PROVED). For every n >= 1, the dominating-set counts of even graphs (von Brömssen) on at most n vertices are exactly the odd numbers up to M(n) = 2^n - n (n odd) or 2^n - n - 1 (n even), except 3, 7, 13 and, for odd n >= 5, X(n) = 2^n - n - 2. So every odd k other than 3, 7, 13 is d(G) for an even graph G. The least order of such a G is the least n with k <= M(n), plus one when k = X(n); it never exceeds ceil(log2(k + 1)) + 1. Tools: the complement formula d(G) = 2^n - D(comp G), E_2-multi-joins, a leaf stable recursion d -> 3*2^(n-1) + d, and a jump family d = 2^(m-1) - D(H) that fills the top band"
status: >
  PROVED, with a computer-assisted base. Independent audit PASS WITH CORRECTIONS (all applied; MISTAKE-747); see
  05-knowledge/results/egcp_20261109/audit/AUDIT.md.
  - Base: the exhaustive census of all graphs on <= 9 vertices (THM-4996 (3), audited P1), plus 268 base graphs on
    <= 9 vertices (248 minimal-order realisers and 20 leaf graphs), each re-verified by brute force
    (verify_base.py).
  - The inductive step is general and uses no computation beyond order 9.
  - Mechanical checks:
    - every branch of the case analysis for 10 <= n <= 24 (proof_check.py);
    - explicit witnesses built exactly as the proof prescribes, brute-forced for 10 <= n <= 12 (construct.py);
    - random tests of the four new lemmas (lemma_check.py);
    - a closure computation giving E(n) = O(n) for all n <= 25 (cover.py).
  - The audit re-implemented Section 4 independently and verified all 7,934 values that need exactly n vertices for
    10 <= n <= 14. It ran an independent exhaustive census for n <= 7, and checked Lemmas 7 and 8 exhaustively at
    small orders.
  - KNOWN and credited:
    - completeness and conciseness over all graphs are Brouwer-Csorba-Schrijver 2009, Proposition 1.2;
    - Lemma 0 is Heinrich-Tittmann and Brown-Nowakowski.
  - New: the even-graph statement (E(n), the exceptional set {3, 7, 13}, X(n) and the +1 cost) and its proof.
session: opus-2026-11-09-S89
date: 2026-11-09
related: [HYP-9575, THM-4996, THM-4531, THM-4721, HYP-9610, THM-5041, THM-5042]
script: 04-computation/experiments/egcp_20261109/
---

# THM-5040 - Sharp even-graph Chan-Pak

Packet: [note](../../05-knowledge/results/egcp_20261109/README.md). The scripts are in
`04-computation/experiments/egcp_20261109/`.

## Statement

**Notation.**
- d(G) is the number of dominating sets of G. It is always odd (Brouwer-Csorba-Schrijver).
- G is **even** if no automorphism reverses an odd number of edges under the reference orientation u -> v for u < v
  (von Brömssen; THM-4996, THM-4531).
- E(n) = {d(G) : G even with at most n vertices}. Padding with isolated vertices keeps evenness and d, so E(n) is also
  the set of values on exactly n vertices.
- M(n) = 2^n - n for n odd, and 2^n - n - 1 for n even.
- X(n) = 2^n - n - 2 for odd n >= 5, and X(n) is empty otherwise.
- O(n) = {odd k : k <= M(n)} minus {3, 7, 13} minus X(n).

**Theorem.** E(n) = O(n) for every n >= 1.

**Corollaries.**
1. **HYP-9575 (B) holds.** {d(G) : G even} is exactly the set of odd numbers other than 3, 7 and 13.
2. **Least order.** Let V(k) be the least order of an even graph with d = k. Then V(k) is the least n with
   k <= M(n), plus one if k = X(n). Note that O(n) is contained in O(n + 1).
3. **Cost of evenness.** ceil(log2(k + 1)) <= V(k) <= ceil(log2(k + 1)) + 1. All graphs attain the lower bound
   (Brouwer-Csorba-Schrijver 2009, Proposition 1.2; re-proved as THM-4996 (1)), so evenness costs at most one vertex.
   - The extra vertex is needed exactly when k lies in a **top band** T(m) = (M(m - 1), 2^(m - 1)) or k = X(m - 1),
     where m = ceil(log2(k + 1)) + 1. For m >= 6 there are 1 + floor((m - 2)/2) such values per order; for
     m = 2, 3, 4, 5 there are 0, 0, 0, 1.
   - Examples: for m >= 4, 2^m - 1 needs m + 1 vertices; 25 = X(5) and 119 = X(7) need 6 and 8.
4. **The unreached values of THM-4996 Proposition B5 are covered.** Its closure missed 12,326 values below 2^22,
   clustered near 2^m - 2^k. Two mechanisms account for such clusters, and Theorem T and Lemma 10 below supply both.
   This explains the clusters; it is not an exact classification of that closure's misses.
   - **Top bands.** T(m) is realisable only on m vertices, and only by graphs with no isolated vertex that are not
     joins (in particular not cones): a join on m vertices has d >= 2^(m-1) + 1.
     - Many top-band realisers are disjoint unions. For example P_3 u G_7 (G_7 the census graph with d = 101) gives
       505 in T(10) (audit, cor4_unions.out).
     - The B5 closure cones only over padded graphs, so it never cones over a top-band realiser.
   - **3, 7, 13 inputs.** Cone or join inputs that fall on 3, 7 or 13 give 2^m + {3, 7, 13} and
     3*2^(m - 2) + {1, 5, 11}.
   - The proof's own closure (cover.py) reaches every value of O(n) for n <= 25.

## 1. Tools

**Lemma 0 (complement formula).** This lemma is KNOWN; the credits follow the proof sketch. A set S fails to dominate
G iff some vertex v outside S has no neighbour in S, that is, S lies in the open neighbourhood of v in the complement
comp(G).
- Hence d(G) = 2^n - D(comp G), where D(H) = |union over v of 2^(N_H(v))| counts the sets contained in some open
  neighbourhood (the empty set included).
- **Unions.** For nonempty H_1, H_2, D(H_1 u H_2) = D(H_1) + D(H_2) - 1, because a nonempty set inside a
  neighbourhood lies in one component.
- **C_4-free graphs.** A set of size >= 2 lies in at most one neighbourhood. So
  D(H) = 1 + #(non-isolated vertices) + sum over v of (2^deg(v) - 1 - deg(v)).

Credits:
- D(H) = N(H, 1) for the neighbourhood polynomial N. The complement formula is the case x = 1 of Heinrich-Tittmann's
  identity D(G, x) + N(comp G, x) = (1 + x)^n (arXiv:1701.03453). Brouwer-Csorba-Schrijver 2009 already note that the
  non-dominating sets form the neighbourhood complex of the complement.
- The union rule and the C_4-free formula are Brown-Nowakowski, Australas. J. Combin. 42 (2008), Lemma 3 and
  Theorem 1. Their formula assumes no isolated vertices; the form above counts them correctly.

**Parity facts (THM-4996 Lemma B2).**
- (B2.1) Adjacent twins make a graph odd. So an even graph has at most one universal vertex.
- (B2.2) Padding keeps evenness.
- (B2.3) If |G| >= 2 and G has no universal vertex, then K_1 v G is even iff G is.
- (B2.5) sgn_G(g) sgn_comp(G)(g) = sgn(g) for every automorphism g.

**Lemma 5' (multi-join sign rule).** Let J = A_1 v ... v A_r with every A_i co-connected (comp(A_i) connected). Then J
is even iff every A_i is even and no two isomorphic A_i have odd order.

*Proof.* Aut(J) = Aut(comp J) is generated by two kinds of map: automorphisms of a single A_i, and swaps of two
isomorphic A_i, A_j along a fixed isomorphism. Label the vertices block by block, labelling isomorphic blocks
compatibly so that the fixed isomorphisms preserve order. The sign does not depend on the labelling.
- An automorphism of A_i alone keeps every cross edge's orientation, so its sign is sgn_(A_i).
- A swap of blocks of order s:
  - preserves the orientation of the edges inside the two blocks;
  - reverses all s^2 edges between them;
  - for every block lying between them in the labelling, reverses an even number 2 s |A_l| of edges.
  So its sign is (-1)^s.
- sgn_J is a homomorphism, so checking these generators suffices. QED

**Lemma 5 (E_2-multi-join).** For an even graph G on g >= 1 vertices and t >= 1, the graph J = E_2^(v t) v G is even
on n = g + 2t vertices, and d(J) = 2^n - 2^g - 2t + d(G). The value formula needs g >= 1; Section 4 uses g >= 3.

*Proof.*
- **Value.** comp(J) = t K_2 u comp(G). Lemma 0 gives D(comp J) = 1 + 2t + (D(comp G) - 1) = 2t + 2^g - d(G).
- **Parity.** The co-components of J are t copies of E_2 (even order) and the co-components of G. The latter satisfy
  Lemma 5' because G is even. QED

**Lemma 6 (cocktail-party graphs).**
- For n even, CP = comp(perfect matching) = E_2^(v n/2) has d = 2^n - n - 1 = M(n) and is even (Lemma 5').
- For n odd, K_1 v CP_(n-1) has d = 2^n - n = M(n) and is even (B2.3).

**Lemma 7 (leaf stable recursion SR).** Let G be even on n' >= 3 vertices, with a leaf q and no universal vertex. Form
SR(G) by adding a hub a adjacent to V(G) - q and a pendant leaf u at a. Then:
- d(SR G) = 3*2^(n' - 1) + d(G);
- SR(G) is even, has a leaf (u) and has no universal vertex.

*Proof.*
- **Value.**
  - Sets containing a: u is free, and q must be dominated by q or its neighbour, which gives 2(2^n' - 2^(n'-2)).
  - Sets without a: they contain u, and their trace on V(G) dominates G, which gives d(G).
- **a is fixed.** a has degree n'. Every other vertex has degree at most (n' - 2) + 1, since G has no universal
  vertex. So a is the unique vertex of maximum degree. Its unique non-neighbour q is also fixed.
- **Automorphisms.** The leaves at a (u and the isolated vertices of G) are non-adjacent twins, so permuting them is
  even. Modulo these permutations, Aut(SR G) = Aut(G)_q.
- **Sign.** Labelling a first, the edges at a never reverse. So the sign of each element is sgn_G = +1. QED

**Lemma 8 (jump family).** Let H have rho >= 2 vertices, let z be a vertex of H with h = deg_H(z) <= rho - 1, and
put k = rho - h >= 1. Form J(H, z): take comp(H) on V(H), add a hub a adjacent to V(H) - z, and add a set U of k
leaves at a. Then m = |J| = 2 rho + 1 - h, and:
- (i) d(J) = 2^(m - 1) - D(H);
- (ii) J has a leaf and no universal vertex;
- (iii) J is even iff sgn_H(g) = sgn(g) for every automorphism g of H fixing z.

*Proof.*
- **(i)** Count the sets containing a and those not containing a.
  - With a: the k leaves are free, and z must be dominated inside comp(H). The failing traces lie inside N_H(z), so
    this gives 2^k (2^rho - 2^h).
  - Without a: all leaves are taken, and the trace dominates comp(H). This gives 2^rho - D(H).
  - The sum is 2^(k + rho) - D(H), since k + h = rho.
- **(ii)** a misses z, every vertex of V(H) misses the leaves, and the leaves have degree 1.
- **(iii)** First, a is the unique vertex of maximum degree, m - 2.
  - If k >= 2, then m - 2 >= rho + 1, while every other vertex has degree <= rho.
  - If k = 1, then z is universal in H. So no other vertex of H is isolated in H, and every vertex of V(H) - z has
    degree <= rho - 1 < rho = m - 2.
  - (Audit: checked exhaustively for rho <= 5.)

  The rest of the argument:
  - So a is fixed, and so is z, its unique non-neighbour.
  - The degree-1 neighbours of a (the leaves U, and the vertices of V(H) - z that are universal in H) are pairwise
    non-adjacent twins. Composing with a permutation of them, which is an even automorphism by B2.1, every
    automorphism maps U to U and V(H) to V(H).
  - It therefore restricts to an element of Aut(comp H)_z = Aut(H)_z, and conversely every such element extends.
  - Labelling a, then the leaves, then V(H), the sign is sgn_comp(H)(g) = sgn(g) sgn_H(g) (B2.5). QED

**Lemma 9 (frontier families).** Each family below is even, has a leaf and no universal vertex (Lemma 8), and
Aut(H)_z is trivial. Write P = P_(rho - 1) u K_1.
- **(F-a)** m = 2 rho, rho >= 3: H = P, z = an end of the path (h = 1). Then D(H) = 2 rho - 3, so
  d = 2^(m - 1) - (m - 3).
- **(F-b)** m = 2 rho, rho >= 5: H = the spider with legs of lengths (rho - 4, 1, 2), z = the end of the first leg
  (h = 1). Then D(H) = 2 rho + 1, so d = 2^(m - 1) - (m + 1) = X(m - 1).
- **(F-c)** m = 2 rho - 1, rho >= 5: H = P, z = the second vertex of the path (h = 2; not fixed by the flip). Then
  D(H) = 2 rho - 3, so d = 2^(m - 1) - (m - 2).

The D values follow from Lemma 0's C_4-free formula. The bound rho >= 5 in (F-c) is needed: at rho = 4 the flip fixes
no vertex but swaps positions 2 and 3, and that J is odd (audit negative control).

## 2. Exclusions: E(n) is contained in O(n)

**Lemma 1.** If G is even on n vertices, then d(G) <= 2^n - n.

*Proof.* G has at most one universal vertex (B2.1). The empty set and the n - 1 other singletons fail to dominate.

**Lemma 2.** If n is even, then d(G) <= 2^n - n - 1.

*Proof.* Suppose d = 2^n - n. Then there is exactly one universal vertex u, and every set of size >= 2 dominates.
- Write G = K_1 v G'. Then G' has no universal vertex, and no pair of vertices lies in a common neighbourhood of
  H' = comp(G').
- So H' is a perfect matching on n - 1 vertices, which forces n to be odd. QED

**Lemma 3.** If n is odd, then 2^n - n - 2 is not in E(n).

*Proof.* Graphs on fewer vertices give d <= M(n - 1) < X(n). So suppose d(G) = 2^n - n - 2 on exactly n vertices.
- **Universal vertex.** Then G = K_1 v G', the non-dominating sets of G are those of G', and G' has no universal
  vertex. So exactly two sets of size >= 2 lie in neighbourhoods of H' = comp(G'). Hence H' has maximum degree <= 2,
  and H' is one of P_4 u sK_2, C_4 u sK_2, 2P_3 u sK_2.
- **No universal vertex.** Then exactly one such set exists, and H = comp(G) = P_3 u sK_2.
- **None is co-even (B2.5).**
  - P_4: its flip reverses 3 edges but is an even permutation.
  - C_4: a rotation reverses 2 edges and is an odd permutation.
  - P_3: swapping its ends reverses 2 edges and is an odd permutation.
  So G is odd. QED

**Excluded small values.** 3, 7 and 13 are realised only by K_2, K_3 and K_4 - e, up to isolated vertices. All three
are odd. This is THM-4996 (3), exhaustive for d <= 511 at every order (Lemma B3).

Lemmas 1-3 and this exclusion give E(n) contained in O(n).

## 3. The top band and the offsets

Write Tj(m) for the j such that 2^(m - 1) - j needs exactly m vertices:
- for m odd, Tj(m) = {1, 3, ..., m - 2};
- for m even, Tj(m) = {1, 3, ..., m - 3} u {m + 1} (the X-value 2^(m - 1) - m - 1).

An **L-graph** is an even graph with a leaf and no universal vertex.

**Theorem T (top band).** For every m >= 8 and every j in Tj(m), there is an L-graph on m vertices with
d = 2^(m - 1) - j.

*Proof.* Induction on m.
- **Base, m = 8, 9.** These are census graphs, listed as the L lines of base_witnesses.txt and verified.
- **Step, m >= 10.**
  - **New j.** Tj(m) minus Tj(m - 2) is {m - 2} for m odd and {m - 3, m + 1} for m even. Lemma 9 (F-c, F-a, F-b)
    supplies these.
  - **Old j.** SR (Lemma 7) carries 2^(m - 3) - j to 3*2^(m - 3) + 2^(m - 3) - j = 2^(m - 1) - j, and keeps the
    L-property. QED

**Lemma 10 (offsets).** For every m >= 8 and every c in {1, 3, 5, 7, 11, 13}, there is an L-graph on m vertices with
d = 2^(m - 1) + c.

*Proof.* Census for m = 8, 9 (base_witnesses.txt). SR maps 2^(m - 3) + c to 2^(m - 1) + c. QED

## 4. The induction

**Claim_3(n).** Every odd x <= 2^n - n - 3 with x not in {3, 7, 13} is in E(n).

For n <= 9 this is the census. Let n >= 10 and assume Claim_3 below n. Take x odd with
2^(n - 1) - n - 2 < x <= 2^n - n - 3; smaller x are covered by Claim_3(n - 1).

- **Case 1: x < 2^(n - 1).** Write x = 2^(n - 1) - j with 1 <= j <= n + 1. Either j is in Tj(n) (Theorem T), or
  x = M(n - 1) (Lemma 6, padded).
- **Case 2a: x > 2^(n - 1) and D := 2^n - x <= 2^(n - 2) + 1.** Note that D >= n + 3.
  - Take the least g with g = n (mod 2) and 2^g + n - g - 1 >= D. Then g <= n - 2.
  - Put d_G = 2^g - D + n - g. It is odd and 1 <= d_G <= 2^g - g - 3 (because D >= n + 3), so Claim_3(g) gives G
    on at most g vertices unless d_G is in {3, 7, 13}. Pad G to exactly g vertices; then E_2^(v (n - g)/2) v G works
    (Lemma 5).
  - If d_G is in {3, 7, 13} and g + 2 <= n - 2, use g + 2 instead. The new value is d_G + 3*2^g - 2 > 13, and it is
    still <= 2^(g + 2) - (g + 2) - 3.
  - Otherwise g = n - 2 and x = 3*2^(n - 2) + {1, 5, 11}. This is a cone over an L-graph on n - 1 vertices with
    d = 2^(n - 2) + {1, 5, 11} (Lemma 10, B2.3).
- **Case 2b: x > 2^(n - 1) and D >= 2^(n - 2) + 3.** Write x = 2^(n - 1) + y with 1 <= y <= 2^(n - 2) - 3, and cone
  over a graph on n - 1 vertices without a universal vertex.
  - If y <= 2^(n - 2) - n - 1 and y is not in {3, 7, 13}, pad a graph from Claim_3(n - 2).
  - If y is in {3, 7, 13}, then x = 2^(n - 1) + y is directly an L-graph (Lemma 10).
  - Otherwise y = 2^(n - 2) - j with 3 <= j <= n. Either j is in Tj(n - 1) (Theorem T), or y = M(n - 2) (Lemma 6,
    padded).

So Claim_3 holds for every n.

**Conclusion.** The odd numbers in (2^n - n - 3, M(n)] are:
- M(n) alone for n even;
- X(n) and M(n) for n odd.

X(n) is not in O(n), so the only value of O(n) above Claim_3's range is M(n), and Lemma 6 realises it. Hence O(n) is
contained in E(n), and together with Section 2, E(n) = O(n). QED

## 5. Verification (04-computation/experiments/egcp_20261109/)

| Check | Script | Output | Result |
| --- | --- | --- | --- |
| Census base: E(n) = O(n) for n <= 9; 248 minimal-order W witnesses and 20 L witnesses (orders 8, 9) recounted and checked for evenness | verify_base.py | verify_base.out | PASS |
| Complement formula and Lemmas 5, 7, 8 on random graphs (1,300 tests; Lemma 8 includes the evenness criterion in both directions, 206 even cases) | lemma_check.py | lemma_check.out | 0 failures |
| Every branch of Section 4 for every x at every 10 <= n <= 24 | proof_check.py | proof_check_24.out | all cases verified |
| Explicit witnesses built by the proof's own recursion and brute-forced (d, evenness, exact order) for 10 <= n <= 12: all delicate values plus a random sample, and all T(n) graphs | construct.py | construct_12.out | 0 failures |
| Corroborating closure under padding, cones, E_2/K_1/E_3 multi-joins, two-hub graphs TH(G, t), the top band and SR | cover.py | cover_25.out | E(n) = O(n) for every n <= 25; also without the injected top band (audit: `--notop`); the exact E0 conjecture misses about 2n values per order, which E never needs |
| Audit: independent re-implementation of Section 4, verifying every value that needs exactly n vertices; census of all labelled graphs on <= 7 vertices; exhaustive Lemma 7 and 8 checks | audit/indep_check.py, indep_census.py, exhaustive_small.py | audit/*.out | 0 failures |

The two-hub graph TH(G, t) is used only in cover.py's corroborating closure:
- **Construction.** Two non-adjacent hubs a, b, both joined to all of V(G), with t pendant leaves at a.
- **Value.** d = 2^(n-1) + 2^f - 2^t + d(G), where f = |G| and n = f + t + 2.
- **Parity.** TH(G, t) is even iff G is even, except at (f, t) = (1, 1), where TH = P_4 is odd. cover.py never uses
  (1, 1).
- **Proof of parity.** For f >= 2 the leaves all hang at a, and a's unique non-neighbour is b, so both hubs are
  fixed. This is the audit's derivation (audit/twohub_probe.out, twohub_parity.out).

## 6. Prior art and credits

- **Completeness and conciseness over all graphs.** Every odd m with 0 < m < 2^n is the dominating-set count of a
  graph on n vertices. This is A. E. Brouwer, P. Csorba, A. Schrijver, *The number of dominating sets of a finite
  graph is odd*, preprint 2009, Proposition 1.2: the same paper that proves the parity, with the
  isolated-vertex/cone induction of THM-4996 Lemma B1. It is restated as Proposition 2.1(iii) of
  Alikhani-Fatehi-Klavžar, arXiv:1512.07514. S. Alikhani, Graphs Combin. 29 (2013) 1175-1181, concerns D(G, -1)
  (connected G with D(G, -1) = n), and is cited by Alikhani-Fatehi for connected realisers. The first version of this
  section credited Alikhani with completeness and called the conciseness new; both are corrected (audit Finding 1,
  MISTAKE-747).
- **Credited from THM-4996:** von Brömssen and Royle-Praeger-Glasby-Freedman-Devillers (even graphs, equinumerous
  with tournaments); Alikhani-Peng (join formula).
- **KNOWN tools.**
  - Lemma 0: Heinrich-Tittmann; Brown-Nowakowski; the BCS remark.
  - Lemmas 7 and 8 are instances of the vertex-addition formula of THM-4996 Lemma B3.
- **No prior art found** for dominating-set counts of even graphs, the exceptional set {3, 7, 13}, the sharp orders
  M(n), X(n), or the +1 cost of evenness.
  - S89 searched the web and grepped the repo for "dominating" and "even graph".
  - The audit searched BCS 2009, Heinrich-Tittmann, Alikhani 2013 and its citing papers, Brown-Nowakowski and RPGFD.
- **Concurrent independent proof.** mac-mini-2026-11-09-digraph721 announced an all-n existence proof of (B) at
  12:49 the same day (THM-5055 draft, audit pending; not in canon when this theorem was pushed). It used:
  - offset moves R = PH(pendant, 1), which is this theorem's SR, and PH(3, 2), a hub missing a vertex with
    |N[u]| = 3, plus two leaves;
  - a seed lemma for each offset.
  Both proofs share the skeleton: cones, padding, E_2-joins. They differ in their seeds (here: the jump family and
  the three sparse H of Lemma 9). The messages crossed (MSG-2471, MSG-1987 and mac-mini's 12:49 reply); the two
  proofs should cross-cite once THM-5055 lands. mac-mini also re-derived the exclusion of X(n) by the deficit route.

## 7. Remarks

- **The tournament analogy, in the repo's binary/ternary dictionary (HYP-9610).**
  - Rédei's theorem (H(T) odd) and Brouwer-Csorba-Schrijver (d(G) odd) are the two parity theorems on the two sides
    of the RPGFD equinumerosity.
  - The exceptional sets are {7, 21} for H on tournaments (Chan-Pak) and {3, 7, 13} for d on even graphs.
  - 7 is shared. Each exceptional value is the count of the smallest odd structure forced by an adjacent-twin
    obstruction.
  - This is a typed ANALOGY. No map carries H to d (THM-4996 (4); THM-4531).
- **Odd graphs.** Odd graphs miss exactly {1, 5, 37}, on the minimal ceil(log2(k + 1)) vertices. This is THM-5056
  (mac-mini, first) and THM-5042 (independent second proof), which answer the question left open in the first version
  of this remark.
- **Second proof.** THM-5055 (mac-mini) is a second, independent proof of this theorem.
- **No universal vertex.** For 5 <= n <= 9, the values without a universal vertex are exactly the odd numbers up to
  M0(n) minus {3, 7, 13, M0(n) - 2}, where M0(n) = 2^n - n - 1 (n even) and 2^n - n - 4 (n odd). This is
  FINITE-EXACT; its all-n form is not needed for the theorem.
