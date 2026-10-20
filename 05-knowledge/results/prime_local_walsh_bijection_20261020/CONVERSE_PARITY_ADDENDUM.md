# Converse-compatible parity blocks of the Sylow-action bijection

**Status: PROVED; FINITE-EXACT for every class through n=6.** For a selected odd prime p, the full Sylow-action bijection can additionally send every self-converse tournament to a graph with an even number of edges, and send each pair of non-self-converse tournaments to two graphs of opposite edge-count parity. Both matched coefficients remain p-units.

## 1. Integral blocks after averaging

Let M[T,X]=Phi_X(T). Global reversal gives

    M[T^op,X]=(-1)^|E(X)| M[T,X].                       (1)

Choose the smaller canonical mask from each non-self-converse pair {T,T^op}. Replace its two rows by their half-sum and half-difference. For a self-converse type retain its row. This change of basis is invertible over R=Z[1/2]. By (1), the half-sum row equals the original row on even-edge columns and is zero on odd-edge columns; the half-difference row is zero on even-edge columns and equals the original row on odd-edge columns. All resulting entries are still integers.

Thus after reordering the matrix is block diagonal:

- E: rows are self-converse types and one representative per nonself pair; columns are even-edge sign-even graphs;
- O: rows are one representative per nonself pair; columns are odd-edge sign-even graphs.

Since the whole transformed matrix is invertible over R, each block is square and invertible over R. In particular

    #even-edge sign-even graphs = #self-converse types + #nonself pairs,
    #odd-edge sign-even graphs = #nonself pairs.

This recovers the inherited converse trace identity; it is not a new total-count theorem.

## 2. Weighted orthogonality survives

For n>=2 set c=2^(binomial(n,2)-1). Distinct original rows of M are orthogonal with graph-automorphism weights, and T and T^op have the same automorphism order. Therefore

    E diag(|Aut X|) E^t = c diag(a_r),

where a_r=2|Aut T| for a self-converse row and a_r=|Aut T| for a pair row. Likewise

    O diag(|Aut X|) O^t = c diag(|Aut T|).

All extra factors are powers of two. Applying the weighted-pairing theorem to each block gives p-unit matchings preserving the valuation of the actual tournament automorphism order. For n=0,1 the unique object handles the statement directly.

Every selected block coefficient is an original Phi_X(T). The signed p-group localization lemma therefore applies without change: the output's Sylow p action is conjugate to that of T. (For a nonself pair both tournaments have the same automorphism group as a subgroup of S_n.)

## 3. The full bijection

Run the canonical cofactor recursion separately on E and O. A self-converse row is sent to its E-matched graph. In each nonself pair, send the smaller canonical tournament representative to the E-matched graph and the other to the O-matched graph. This is a permutation of all types, because the block matchings are permutations of their disjoint target sets.

It satisfies simultaneously:

1. p does not divide the matched signed coefficient;
2. the Sylow p-subgroup actions are S_n-conjugate;
3. self-converse types have even-edge outputs;
4. the outputs of each nonself converse pair have opposite edge parity.

No intrinsic involution on even-graph types is asserted beforehand. The bijection **defines** an involution by transport of converse: it fixes the selected self-converse images and exchanges the two selected images of each nonself pair. Which even-edge graphs are fixed depends on the block matching. This is not graph complementation.

## 4. Refined signed identities

For each Sylow p-action conjugacy class C, let SC_C count self-converse tournament types of that action class. Let E_C^+ and E_C^- count sign-even graph types of that action class with respectively even and odd edge count. The matched blocks prove

    SC_C = E_C^+ - E_C^- .                              (2)

The corresponding counts for nonself converse pairs are E_C^-. Summing (2) against any function of the Sylow p-action class gives a refined trace identity. For example, for every k>=0,

    #{self-converse T: v_p(|Aut T|)=k}
      = sum over even graph types X with v_p(|Aut X|)=k
          of (-1)^|E(X)|.                              (3)

The same holds after filtering by any chosen family of p-subgroup-presence predicates, since these are functions of the Sylow action. This refinement is stronger than the unweighted self-converse count.

## 5. Finite controls

`walsh_converse.py` reads the frozen complete atlas `walsh_support.json`. It derives converse by actual vertex relabellings, constructs both integer blocks, checks their determinants are powers of two, verifies weighted orthogonality exactly, and constructs the blockwise matchings for p=3,5,7. It verifies p-unit support, valuation equality, self-converse parity, pair-opposite parity, and the signed valuation distribution. Output is `walsh_converse.json` and `walsh_converse.out`.

    python3 -B walsh_converse.py

The complete counts (self-converse types, nonself pairs, even-edge graphs, odd-edge graphs) are:

    n=3: (2,0,2,0)
    n=4: (2,1,3,1)
    n=5: (8,2,10,2)
    n=6: (12,22,34,22).

The script checks all 0<=n<=6, including empty odd blocks. Its numeric controls establish the stated valuation and parity checks; the Sylow-action refinement follows from the proved orbit-localization lemma, unless separately checked by the subgroup-control script.
