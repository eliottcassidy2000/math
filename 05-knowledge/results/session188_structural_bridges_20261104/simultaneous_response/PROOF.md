# Pairwise-disjoint reversals cannot cancel in the full named response

**PROVED:** any nonempty family of pairwise-disjoint nontrivial modules is separated from its simultaneous converse by an actual two-port count difference equal to 2 modulo 4. **FINITE-EXACT:** all simultaneous original-module families through five old names, and specified complete six-vertex controls. **ADDITIONAL PROVED EXTENSION:** `ANCHORED.md` covers every nested pair and, more generally, every laminar family with an exposed anchor. **OPEN in this packet:** fully covered laminar families and unrestricted overlapping families. No exact full-response collision is claimed.

## 1. Operation and observer

All names stay fixed. F_T(a,c) counts Hamilton paths after adding x->y, with a_u=1 meaning u->x and c_u=1 meaning u->y. Every comparison retains the complete named Boolean polynomial of bidegree at most (2,2), including all cubic and quartic terms.

A simultaneous family here consists of distinct nontrivial modules of the **original** tournament. An old edge is reversed precisely when its two endpoints belong to an odd number of chosen modules. For overlapping modules this means simultaneous edge XOR; it does not assume that an original module remains a module after an earlier reversal. Pairwise-disjoint families have no edge overlap, so their operations also form a legal sequence of module reversals.

## 2. Scalar invariance of a module converse

For a nonempty tournament A, q_k(A) is the number of ordered spanning covers by k nonempty directed paths. Reversing every path, while preserving the order of the k path slots, gives a bijection

    q_k(A)=q_k(A^op).

The occupied-gap substitution formula is

    H(Q[A])=sum_k c_k(Q,z) q_k(A),

where the nonnegative integer c_k counts exterior skeletons with k module runs. Thus reversing any module preserves H exactly, in every fixed surrounding tournament. The same applies after ports are added whenever that set remains a module in the actual extension. This is an exact scalar statement; it does not assert equality of arbitrary port responses.

## 3. Disjoint noncancellation theorem

Let M_1,...,M_s be pairwise-disjoint modules of T, each with at least two names, and s>=1. Let U reverse every M_i internally.

**Theorem.** There is a physical two-port query such that

    F_T(a,c)-F_U(a,c) = 2 modulo 4.

**Proof.** The incoming [E3-E4 parity theorem](https://github.com/eliottcassidy2000/math/blob/76b3b1a77d124100f03026420f6b821533cdec15/05-knowledge/results/two_port_parity_separation_20261104/PROOF.md) supplies an intrinsic query on M_1 with converse contrast 2 modulo 4. Give both ports M_1's original exterior incidences. Incoming E4 proves that this lifted query separates T from the tournament T_1 obtained by reversing only M_1, by 2 modulo 4.

For i>1, every old name outside M_i still sees M_i uniformly in T_1 and after any other disjoint reversal. Arcs between the two disjoint modules M_1 and M_i are uniform. Therefore each new port also sees every vertex of M_i uniformly under the chosen lift. Thus M_i remains a module in each actual port extension.

By Section 2, reversing M_i preserves the extended tournament's Hamilton count exactly. Do this for every i>1. It follows that F_(T_1)(a,c)=F_U(a,c) as integers. Combining with the E4 residue gives the stated contrast. QED.

The proof works for arbitrarily many modules, arbitrary internal sizes, and arbitrary strong or nonstrong ambient tournaments. It includes the incoming one-module theorem as s=1. It is a separator for this comparison class, not a reconstruction theorem for unrelated tournaments.

## 4. Exact bounded census

`search.py` enumerates every labeled tournament on one through five old names, every original nontrivial module including the whole set, and every nonempty family of distinct such modules. It compares **entire integer vectors**, not hashes or selected invariants. A complete response bank is retained for each order.

The finite universe contains 1,099 tournaments and 130,942 source/family presentations:

- 3,830 single-module families;
- 824 pairwise-disjoint families with at least two members;
- 12,228 non-disjoint laminar families;
- 114,060 families containing a crossing pair.

No nonempty family has an identity edge-XOR in this finite universe. None has an exact named full-response collision. At order five alone there are 5,120 exactly-two nested families, 2,040 exactly-two crossing families and 800 exactly-two disjoint families; these are presentations, not isomorphism classes.

This enumeration is a simultaneous-operation control with exact family semantics and complete retained vectors. It is **not** presented as a new general reconstruction bound or an unbounded proof for arbitrary nesting or overlap.

## 5. Six-vertex laminar hostile: uniform block probes lose the parity separator

Use the strong tournament C3[TT2,TT2,TT2] on pair modules {0,1}, {2,3}, {4,5}:

    T=(14,12,56,48,35,3).

Reverse all three pair interiors and the whole six-name set. The result preserves the three internal pair orders and converses the quotient:

    U=(50,48,11,3,44,12).

Both have H=45. Their exact full difference has 36 nonzero quadratic, 132 nonzero cubic and 132 nonzero quartic coefficients, and no constant or linear terms. All 484 named coordinates of bidegree at most (2,2) of each full response and their difference are retained.

If both attachment masks are constant on each pair, all 64 resulting queries agree modulo four. They are not equal as integers: 24 differ, including

    F_T(3,12)=509,       F_U(3,12)=545.

A nonuniform query restores a modular separator:

    F_T(1,5)-F_U(1,5)=2 modulo 4.

The receipt gives the exact counts and all 4,096 named queries on this pair, independently evaluated by Hamilton subset DP and compared with both full polynomial expansions. This example refutes transporting the incoming mod-four quotient separator to arbitrary block-uniform probes without proof. It is **not** an exact full-response collision and does not refute a general laminar noncancellation conjecture.

Two other six-vertex controls reverse the disjoint pairs {0,1},{2,3}, and the nested pair {0,1} and the whole set. Their differences retain 108 and 84 nonzero quartics respectively. A crossing-interval control on five names changes H from 1 to 11. Every full response vector and each complete interpolation battery are retained.

## 6. Boundaries and reproducibility

The disjoint theorem has been proved above. The separate `ANCHORED.md` proves every nested pair using exact singleton-port averaging, without a homogeneous-complement assumption. Exactly-two crossing modules have a separate positive-cover proof in companion work. None of these arguments is silently extended to fully covered laminar families or arbitrary crossing families. In a nested reversal, an outer operation changes the selected inner module's exterior incidences, so the same lifted ports need not form a module on both sides.

Run from this directory, retaining the sibling `response/verify.py` dependency:

    python search.py > receipt.json
    python -O search.py > receipt.optimized.json

The program writes `full_response_n1.json` through `full_response_n5.json`, containing every complete named vector including zeros. Their exact SHA256 values are pinned in the receipt. All gates use explicit exceptions. Normal and optimized results must agree byte for byte. The all-order conclusion is the invariant-profile argument and incoming E4, not the finite search.
