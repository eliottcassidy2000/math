# Partition-OCF Correspondence: Position Topology and Independent Sets

**Author:** opus-2026-03-06-S11b (continued^6)
**Status:** STRUCTURAL INSIGHT, computationally verified at n=5,7,9

---

## The Correspondence

### Statement

The surviving position patterns in the W-polynomial coefficient w_{n-1-2k}
are indexed by partitions of k into positive integers. Each partition
(j_1, j_2, ..., j_m) corresponds to an independent set type in the conflict
graph Omega(T): a collection of m vertex-disjoint directed odd cycles of
lengths (2j_1+1, 2j_2+1, ..., 2j_m+1).

### The Map

**Position side:** A run structure (2j_1, 2j_2, ..., 2j_m) in e_{2k} consists
of m blocks of even lengths 2j_1, ..., 2j_m. Each block of length 2j involves
2j+1 vertices forming an odd-sized subtournament.

**OCF side:** An independent set of Omega(T) consists of vertex-disjoint
directed odd cycles. A cycle of length 2j+1 lives on 2j+1 vertices.

**The bijection:** Block of length 2j <-> directed (2j+1)-cycle.

### Why It Works

1. **Even Block Vanishing:** Blocks of ODD length contribute zero because
   even-sized subtournaments have c_0 = 0 (only odd powers survive for even n).
   Therefore, only even-length blocks survive.

2. **Singleton Cancellation:** Isolated positions cancel by vertex-swap
   symmetry. So only contiguous blocks contribute.

3. **Factorization:** Disjoint blocks involve disjoint vertex sets, so
   their contributions multiply independently.

4. **OCF on subtournaments:** Each block of length 2j contributes c_0(S)
   for a (2j+1)-vertex subtournament S. Using OCF recursively,
   c_0(S) encodes the odd-cycle structure of S.

### The Partition Count

Level k (coefficient w_{n-1-2k}) has p(k) surviving pattern types,
where p(k) is the partition function:

| k | p(k) | Pattern types |
|---|------|---------------|
| 0 | 1    | () [universal] |
| 1 | 1    | (2,) [3-cycle] |
| 2 | 2    | (4,) [5-cycle], (2,2) [two 3-cycles] |
| 3 | 3    | (6,) [7-cycle], (4,2) [5+3 cycle pair], (2,2,2) [three 3-cycles] |
| 4 | 5    | (8,), (6,2), (4,4), (4,2,2), (2,2,2,2) |
| 5 | 7    | seven pattern types |

### Consequences

1. **The OCF sum H = sum w_{2k}/4^k recovers I(Omega,2) = sum alpha_m * 2^m**
   because each partition at level k contributes to a specific alpha_m
   (where m = number of parts in the partition).

2. **Perpendicularity of invariants:** Different partitions at the same level
   capture orthogonal aspects of tournament structure because they correspond
   to structurally independent position topologies.

3. **Recursive proof structure:** The (2j,)-type block uses OCF at size 2j+1
   (smaller n), while multi-part patterns are purely algebraic (or use OCF
   at even smaller sizes). This gives an inductive scheme.

---

## Connection to Hypercube Geometry

A tournament T is a point in {0,1}^{C(n,2)} (the edge hypercube).
The complement map T -> T^op is the antipodal map.
The W-polynomial, being invariant under complement (even powers only),
lives on the quotient hypercube / Z_2.

Each coefficient w_{2k} is a polynomial of degree 2k in the centered
edge variables s_e = A_e - 1/2. The partition decomposition shows that
this degree-2k polynomial splits into p(k) orthogonal components,
each capturing a distinct type of cycle interaction.

This is a spectral decomposition of the tournament's position in the
quotient hypercube, filtered by cycle complexity.

---

## Open Questions

1. Can the partition correspondence be used to prove OCF inductively?
   (The obstacle is c_0 = w_0, the "full block" contribution.)

2. Does the correspondence extend to even n? (At even n, only odd
   powers survive, giving a dual picture.)

3. Is there a generating function identity that encodes all levels
   simultaneously? (Something like: W(r) = sum over independent sets S
   of prod_{cycles C in S} g(C, r), for some function g?)
