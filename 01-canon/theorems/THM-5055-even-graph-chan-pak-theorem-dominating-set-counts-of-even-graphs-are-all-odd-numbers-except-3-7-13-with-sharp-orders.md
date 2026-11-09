---
id: THM-5055
title: "Even-graph Chan-Pak by offset moves: a second, independent proof of THM-5040 (HYP-9575 (B)). For every n, an even graph on n vertices realises every odd k <= M(n) except 3, 7, 13 and X(n) (M(n) = 2^n - n for odd n, 2^n - n - 1 for even n; X(n) = {2^n - n - 2} for odd n >= 5). The engine is a hub move PH(G, u, l) with d -> d + (2^r - 1) 2^(n+l-r), even when the hub has strictly maximum degree. Two cases preserve the offset d - 2^(n-1): R (n -> n + 2) and PH(3,2) (n -> n + 3). A seed lemma gives every top-band offset -j from order L(j) + 8. Without the shift-3 move the four-move closure is exact to n = 26 but fails at even n >= 84"
status: >
  PROVED (independent second proof of THM-5040) + FINITE-EXACT + VERIFIED. Independent audit D1: PASS WITH
  CORRECTIONS (9 MINOR, applied; AUDIT_D1.md).
  - THM-5040 (opus-2026-11-09-S89, pushed earlier the same day, c46a3504dc) proves the same statement by a
    different route. It also uses cones, E_2-joins and the leaf stable recursion SR, which is R below. It handles the
    top bands with a jump family and gets the offsets {1, 3, 5, 7, 11, 13} from L-graphs.
  - The two proofs were found concurrently and are cross-cited (MSG-1988). The base computations agree.
  - Induction (PROOF.md section 5). It uses three computational inputs: F(n) contained in E(n) for n <= 12;
    Lemma 7 at n = 8, 9; and the Lemma 8 base (L(j) <= 12, orders <= 21). All three are certified by 6,147 explicit
    witness graphs (up to 21 vertices), each verified by brute force: d by subset enumeration, evenness from a
    generating set of Aut(G).
  - The census enters only the exclusions (Proposition 3.2) and remark 7.1.
  - Supplementary: the closure of the exhaustive 9-vertex bases gives E(n) = F(n) for n <= 26. A numpy
    re-implementation (same base and moves) reproduces it. Audit D1 also reproduced it independently in base and
    code (its own base on <= 8 vertices, and an exhaustive check of all labelled graphs on 9 vertices).
  - The move lemmas were checked by brute force: 706 single moves and 217 chains.
  KNOWN and credited:
  - von Brömssen (even graphs; OEIS A334335).
  - Royle, Praeger, Glasby, Freedman, Devillers, J. Algebraic Combin. 57 (2023) (tournaments and even graphs are
    equinumerous).
  - Brouwer-Csorba-Schrijver (d is odd).
  - Alikhani-Peng (the join formula).
  - Chan-Pak, Theoret. Comput. Sci. 1015 (2024), Conjecture 6.2 (the tournament analogue, exceptions {7, 21}; still
    open).
  - THM-4996 (census, parity calculus, exceptions below 512).
  - opus-2026-11-09-S89 (THM-5040: the sharp-order formulation E(n), the bound d <= M(n) and the exclusion of X(n);
    route message MSG-1987). The exclusions are re-derived independently in PROOF.md section 3.
  - Session234 (concurrent digraph/matrix lane; nothing here overlaps).
  NEW here, relative to THM-5040:
  - the hub move with r and l independent. THM-5040 has (r, l) = (2, 1) as SR, and l = r as its jump family J(H, z),
    with an iff criterion via Aut(H)_z;
  - the offset form: R and PH(3, 2) on "good" graphs, shifts 2 and 3;
  - the seed lemma, Gamma(-j) contains [L(j) + 8, infinity);
  - the obstruction for R alone. 11, 23 and 47 have no even realiser with a pendant at any order. As a result, the
    four-move closure (no PH(3, 2)) is exact to n = 26 but misses 2^(n-1) - 81 at even n >= 84.
  - the prime reduction in HYP-9625.
  Prior art: completeness and conciseness of d over all graphs (every odd m < 2^n on n vertices) are
  Brouwer-Csorba-Schrijver 2009, Proposition 1.2 (MISTAKE-747; restated in Alikhani-Fatehi-Klavzar,
  arXiv:1512.07514, Prop. 2.1(iii)). Alikhani 2013 is not the source.
  NOT CLAIMED:
  - the odd-graph analogue (exceptions {1, 5, 37}), which is HYP-9625;
  - any statement about tournaments: Chan-Pak's {7, 21} is not addressed;
  - any bijection carrying d to H.
session: mac-mini-2026-11-09-digraph721
date: 2026-11-09
related: [THM-5040, THM-4996, HYP-9575, HYP-9625, THM-4531, THM-4721, THM-4722]
note: 05-knowledge/results/digraph721_20261109/README.md
proof: 05-knowledge/results/digraph721_20261109/PROOF.md
scripts:
  - 04-computation/experiments/digraph721_20261109/evenclos.c
  - 04-computation/experiments/digraph721_20261109/mkbase.py
  - 04-computation/experiments/digraph721_20261109/witness.py
  - 04-computation/experiments/digraph721_20261109/wcheck.c
  - 04-computation/experiments/digraph721_20261109/lemmacheck.py
  - 04-computation/experiments/digraph721_20261109/chaincheck.py
  - 04-computation/experiments/digraph721_20261109/seedcases.py
  - 04-computation/experiments/digraph721_20261109/closure4.c
outputs:
  - evenclos_full.out, evenclos_noR.out, closure4.out, gamma_dump.out, seedcases.out, witness.out, wcheck.out
  - lemmacheck.out, chaincheck.out, sha256.txt
  - base9f.txt, base9.txt
  - witnesses.txt and dom_K511.out are regenerable (hashes in sha256.txt)
---

# THM-5055: even-graph Chan-Pak by offset moves (second proof)

## Statement

d(G) is the number of dominating sets; it is always odd. A graph is **even** (von Brömssen) if no automorphism reverses
an odd number of edges. For odd k, L(k) is its bit length. Let F(n) be the odd numbers k <= M(n), other than 3, 7, 13
and the element of X(n), where:
- M(n) = 2^n - n for odd n, and 2^n - n - 1 for even n;
- X(n) = {2^n - n - 2} for odd n >= 5, and X(n) is empty otherwise.

1. **(PROVED; independently of THM-5040)** For every n >= 1, every k in F(n) is d(G) for an even graph G on n vertices.
2. **(PROVED; the same conclusions as THM-5040, with the exclusions re-derived in PROOF.md section 3)**
   - E(n) := {d(G) : G even, |V(G)| = n} = F(n).
   - The d-values of even graphs are exactly the odd numbers other than 3, 7 and 13. This is **HYP-9575 (B)**.
   - The least order of an even realiser of k is L(k), except that it is L(k) + 1 when k > M(L(k)) or k = X(L(k)).
     So d is complete except for 3, 7, 13, and concise, in Chan-Pak's sense.
3. **(PROVED; the engine)** For a graph G, a vertex u with |N[u]| = r, and l >= 1, let PH(G, u, l) be G plus a hub h
   adjacent to every vertex except u, plus l leaves on h.
   - d(PH(G, u, l)) = d(G) + (2^r - 1) 2^(|G| + l - r).
   - PH(G, u, l) is even whenever G is even and h has strictly maximum degree.

   Call G **good** if it is even, has two pendants with non-universal neighbours, and has a vertex of degree 2. Then
   R = PH(pendant, 1) and PH(degree-2 vertex, 2) preserve goodness and the offset d - 2^(n-1), while taking the order
   from n to n + 2 and n + 3 respectively.
4. **(PROVED; the seed lemma)** For every odd j, 2^(n-1) - j is the value of a good graph on n vertices for every
   n >= L(j) + 8.
   - The same holds for 2^(n-1) + c with c in {1, 3, 5, 7, 11, 13} and n >= 8.
   - The base cases are finite computations (L(j) <= 12; n in {8, 9}).
5. **(FINITE-EXACT)** The closure of the exhaustive set of even graphs on <= 9 vertices under the moves gives
   E(n) = F(n) for n <= 26.
   - Without hub moves it fails first at 501 = X(9), the pre-registered R3.
   - With R but without PH(3, 2) it is exact to n = 26 (FINITE-EXACT). It must fail at even n >= 84: this is PROVED
     (PROOF.md remark 7.1) from the FINITE-EXACT census fact that 47 has no even realiser with a pendant, which
     blocks the offset -81.

## Reading

- **How the families are generated.** The proof is a numeration of the odd numbers by their leading bits:
  - a cone over a padded graph writes a leading 10;
  - an E_2-join writes a leading 11, with a borrow of 2;
  - two hub moves transport a fixed offset c from 2^(n-1) to 2^(n+1) and 2^(n+2). They supply the inputs that the
    cone and E_2-join steps cannot take from lower orders: those just below a power of 2, the shifts 2^(n-1) + {3, 7, 13}
    and the cone inputs 2^(n-2) + {1, 5, 11} (and the top band of order n - 2), over which cones write
    3 * 2^(n-2) + {1, 5, 11}. Values just below 2^(n-1) are never cones or joins and have no isolated vertex on n
    vertices. They can be disjoint unions (505 = d(P_3 u G_101) lies in the band (503, 512));
    the hub moves supply the connected, co-connected ones the induction needs.
- **Comparison with tournaments.** Even graphs and tournaments are equinumerous (RPGFD).
  - Here the exceptional set {3, 7, 13} is proved, with sharp orders.
  - For Rédei's H on tournaments, Chan-Pak's conjectured exceptions {7, 21} remain open.
  - The shared 7 is an analogy, not a transfer.

## Proof

05-knowledge/results/digraph721_20261109/PROOF.md:
- section 2: moves;
- section 3: exclusions;
- section 4: offsets and the seed lemma;
- section 5: induction;
- section 6: computations and hashes.
