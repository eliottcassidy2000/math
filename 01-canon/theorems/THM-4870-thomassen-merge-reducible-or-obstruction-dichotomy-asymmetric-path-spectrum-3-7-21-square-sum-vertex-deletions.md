---
id: THM-4870
title: "The Thomassen merge. (A) Reducible-or-obstruction: a uniquely Hamiltonian tournament of order n has r non-critical vertices (each a reduction T -> T - v within the class) and k cut arcs (each a one-arc cut), with r + k = n + 1, so max(r, k) >= (n + 1)/2; over the F(2n - 6) classes, E[r] = ((5 + √5) n + 12 - 8√5)/10 + O(φ^(-4n)), i.e. the fractions φ/√5 and (5 - √5)/10, an exact analogue of the linearly-many-reducible-configurations-or-obstructing-cycles theorem for planar triangulations (Inoue-Kawarabayashi-Miyashita-Mohar-Thomassen-Thorup 2026). (B) The trivial group and {7, 21}: H(T) = 3 forces Aut(T) ⊇ Z_3, and for every odd H <= 1125, H is the Hamiltonian-path count of an asymmetric tournament iff H is not 3, 7 or 21. (C) Square-sum graphs: no Q_N is hypohamiltonian or hypotraceable; the vertex-deleted subgraphs of the non-traceable orders, of Q_23-Q_31, and the non-critical sets of the uniquely Hamiltonian Q_32 and Q_33 (6 and 9 vertices, pairwise non-adjacent, never cycle-consecutive). (D) Uniquely Hamiltonian graphs: no two consecutive non-critical cycle vertices for n <= 9, but 1,300 labelled examples at n = 10 (HYP-9462 (B) refuted). Also: a tournament automorphism group is a planar (Maschke) group iff it is cyclic."
status: >
  PROVED (elementary; kind-pasteur-2026-11-04-S10): (A) from THM-4868 (D) and (H) and singularity analysis of the
  rational generating function; (B) the forcing of Z_3; (C1) from THM-4505 and THM-4851 (a); the Maschke remark (odd
  order, Maschke 1896, CITED).
  FINITE-EXACT: (A) the exact mean checked for n <= 70; (B) for H <= 1125 (complete census below 225, sufficient
  construction on [225, 1125], pre-registered as HYP-9462 (A1) in ffedbba714 and CONFIRMED); (C2)-(C3); (D)
  exhaustive for n <= 10 on labelled graphs containing a fixed Hamiltonian cycle.
  CITED: arXiv:2603.24880 (the 2026 four-colour theorem); Thomassen 1974/1976 and MathWorld for the Thomassen graphs;
  Royle-Praeger-Glasby-Freedman-Devillers 2023 (equinumerosity); THM-4747 (|Aut| | H); THM-4812; THM-1370; THM-4810.
  INTERPRETATION / NUMEROLOGY: the four-colour analogy, the fives, and the coincidences of the note's section 7. The
  {6,4,8,9,12} + 105 co-occurrence in the LRC handoff file is a coincidence.
  REFUTED (recorded here as (D)): HYP-9462 (B).
  AUDIT: see the note, section 10.
session: kind-pasteur-2026-11-04-S10
source: 05-knowledge/results/thomassen_merge_hypohamiltonian_groups_20261104.md
related:
  - THM-4868
  - THM-4867
  - THM-4866
  - THM-1370
  - THM-4747
  - THM-4810
  - THM-4812
  - THM-4743
  - THM-4851
  - THM-4505
  - THM-4720
  - THM-4721
  - HYP-9462
  - HYP-9460
  - HYP-9458
  - HYP-9340
verification:
  - 04-computation/experiments/asym_spectrum_autspec_20261104.c
  - 04-computation/experiments/asym_spectrum_closure_20261104.py
  - 04-computation/experiments/uh_hpath_values_20261104.c
  - 04-computation/experiments/uh_graph_noncritical_20261104.c
  - 04-computation/experiments/uh_graph_downset_20261104.c
  - 04-computation/experiments/square_sum_vertex_deleted_20261104.c
  - 04-computation/experiments/thomassen_merge_20261104_explore.out
  - 04-computation/experiments/thomassen_merge_20261104_confirm.out
---

# THM-4870 - the Thomassen merge

## (A) Reducible-or-obstruction for uniquely Hamiltonian tournaments

Let T be uniquely Hamiltonian of order n >= 4, with r non-critical vertices and k cut arcs on its cycle.
- Every non-critical v gives T - v uniquely Hamiltonian, a reduction (THM-4868 Lemma 1).
- Every cut arc gives a one-arc cut, an obstruction (THM-4867).
- r + k = n + 1 (THM-4868 (D)), so max(r, k) >= (n + 1)/2.
- Both extremes occur: A_n has (r, k) = (2, n - 1), and 2^(n-5) classes have (n - 1, 2).

Over all classes of order n,

    E[r] = ((5 + √5) n + 12 - 8√5)/10 + O(φ^(-4n)).

*Proof.* E[r] = [x^n] ∂_t G(x, 1) / [x^n] G(x, 1) for the generating function of THM-4868 (H).
- G(x, 1) = x^4 / (1 - 3x + x^2) has a simple pole at φ^(-2).
- ∂_t G(x, 1) = x^4 (x^3 - 2x^2 - 3x + 2) / (1 - 3x + x^2)^2 has a double pole there.
- Partial fractions give the linear mean, with slope -ρ'(1)/ρ(1) = (5 + √5)/10. Here ρ(t) is the dominant root of
  1 - (1 + 2t) x + t x^2.
- The constant was checked exactly against the rows for n <= 70.

This is the exact counterpart of the 2026 triangulation theorem: linearly many reductions, or linearly many nested
obstructions.

## (B) The trivial group, Z_3 and {7, 21}

- If H(T) = 3, then T has a strong component with H = 3, which is C_3, so Aut(T) ⊇ Z_3.
- **FINITE-EXACT.** For odd H <= 1125, an asymmetric tournament with H Hamiltonian paths exists iff H is not 3, 7 or 21.
  - Below 225 this uses every strong class of order <= 10 with its automorphism group. Strong tournaments of order
    >= 11 have H >= 225.
  - On [225, 1125] it uses products of asymmetric strong classes of order <= 10 and of uniquely Hamiltonian
    tournaments, which are asymmetric (THM-4868 (E)).
  - The only strong value <= 1125 among orders <= 10 that needs a symmetry is 3.

Whether the statement holds for all H is HYP-9462 (A), which strengthens the open Chan-Pak completeness.

Tournament automorphism groups have odd order (no C_2), and among Maschke's planar groups the odd ones are cyclic. So a
tournament's automorphism group is planar iff it is cyclic. The first strong non-cyclic examples occur at order 7:
Z_3 × Z_3 (H = 99) and F_21 (P_7).

## (C) Square-sum graphs Q_N

- **(C1)** No Q_N is hypohamiltonian or hypotraceable.
  - For N <= 30 the minimum degree is 1 (THM-4505).
  - Q_31 has minimum degree 2 and is traceable.
  - For N >= 32, Q_N is Hamiltonian (THM-4851).
- **(C2)** For the non-traceable N = 18, 19, 20, 21, 22, 24, the vertices v with Q_N - v traceable number 1, 1, 1, 1,
  2, 4.
  - For Q_23, Q_25 and Q_26, Q_N - v is Hamiltonian only for v = 18.
  - For Q_31, exactly for v in {2, 4, 13, 15, 26, 28, 29}.
- **(C3)** The uniquely Hamiltonian Q_32 and Q_33 have the non-critical sets {2, 13, 15, 26, 29, 30} and
  {2, 13, 15, 18, 26, 29, 30, 32, 33}. Each is independent in the graph, and no two of its vertices are consecutive on
  the cycle.

## (D) Uniquely Hamiltonian graphs (the graph side of THM-4868 (C))

For labelled graphs containing the fixed cycle 0..n-1 with no other Hamiltonian cycle (exhaustive down-set census):

| n | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---:|---:|---:|---:|---:|---:|---:|
| such graphs | 3 | 11 | 69 | 561 | 6,903 | 113,107 | 2,689,657 |
| with two cycle-consecutive vertices both non-critical | 0 | 0 | 0 | 0 | 0 | 0 | 1,300 |

- At n = 10 the 1,300 are exactly the graphs with 6 > n/2 non-critical vertices.
- A sparsest witness is the 10-cycle plus the chords {1,6}, {2,8}, {3,7}, {4,9}.

So the tournament run of THM-4868 (C) has no graph analogue in either direction, and HYP-9462 (B) is refuted. QED.
