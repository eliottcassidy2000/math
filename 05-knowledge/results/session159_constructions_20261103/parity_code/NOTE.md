# Parity-code coefficients for crossing tournaments

Status: proved from explicit cycle classification and the previously established odd-cycle packing identity. This gives exact strong-tournament invariants and a coding-theoretic bridge, not a full Hamilton-path spectrum theorem. Novelty relative to the literature is not claimed. Finite enumerations below are controls, not proofs.

## General disjoint-endpoint theorem

Start with the transitive tournament on ordered vertices 0,...,n−1, and reverse exactly the k edges b_i→a_i where a_i<b_i and the 2k endpoints are all distinct. Let J_T(z)=∑_j p_j z^j, where p_j is the number of vertex-disjoint collections of j directed odd cycles. For each non-endpoint v let I(v)={i:a_i<v<b_i}. Then

    p_k = 2^(−k) ∑_(ε∈{−1,1}^k) (∏_i ε_i)
                         ∏_(v outside endpoints) (1+∑_(i∈I(v)) ε_i).   (1)

All p_j with j>k vanish.

Proof. Every directed cycle uses at least one reversed edge, since a directed cycle cannot use solely increasing edges. Disjoint cycles use distinct reversed edges, so at most k occur. In a packing of k cycles each cycle uses exactly one reversed edge and every reversed edge is used. A directed cycle with unique reversed edge b_i→a_i is precisely the increasing path from a_i through an odd number of internal vertices to b_i followed by the reversed edge. For a k-packing all 2k endpoints are occupied, and each remaining vertex is unused or belongs to one eligible interval; each interval must receive an odd number of vertices. These assignments are in bijection with k-packings.

The product in (1) is the generating polynomial for those assignments with the indeterminate for interval i specialized to ε_i. Averaging against ∏ ε_i is the ordinary Fourier projection onto odd exponent in every coordinate. Equivalently, use k binary incidence words, impose one odd-parity check per word, forbid a vertex from appearing in two words, and set entries outside their allowed intervals to zero. Formula (1) evaluates this precisely specified constrained parity code; it does not assert that the forbidden-symbol constraint itself is linear.

If every b_i−a_i≥2 and the intervals [a_i,b_i] cover the vertex set and have connected intersection graph, the tournament is strong. Each interval contains the directed cycle a_i→a_i+1→...→b_i→a_i because no reversed edge has consecutive endpoints. The interval cycles share vertices along a connected intersection graph and cover all vertices, so their union is strongly connected.

## A strong k-parameter family

For k≥1 and m≥0, use the vertex set 0,...,2k+m−1 and reverse the k edges

    (k+m+i) → i,  0≤i<k.

Exclude the degenerate case k=1,m=0, in which there is no directed cycle. All intervals overlap and cover, and all lengths are at least 2, so T is strong. The m non-endpoints lie in every interval. Formula (1) simplifies to

    p_k(k,m) = 2^(−k) ∑_(j=0)^k (−1)^j binom(k,j)(k+1−2j)^m.   (2)

This counts length-m words over {0,1,...,k} in which each nonzero symbol appears an odd number of times. In particular, p_k(k,k)=k!, since each nonzero symbol must appear exactly once. Thus this explicit strong tournament on 3k vertices has k! distinct packings of k disjoint directed odd cycles. This exhibits arbitrarily high odd-cycle packing degree in strong tournaments; it does not itself compute their full H.

## Full packing polynomial for the S158 crossing gadget

Let A,B,C≥1 and a=0,c=A,b=A+B,d=A+B+C. Reverse b→a and d→c. Put δ_A=[A=1], δ_C=[C=1] and δ=δ_Aδ_C. Every directed odd-cycle packing has at most two members, so

    J_T(z)=1+Nz+Pz²,
    P = [2^(A+C−2)3^(B−1) − δ_C 2^(A−1)
         − δ_A 2^(C−1) + δ(−1)^(B−1)] / 4.                 (3)

Proof of (3). The A−1 left-outer free vertices can be unused or assigned to the first edge; the C−1 right-outer ones to the second; and the B−1 common free vertices can be unused or assigned to either edge. Consequently P is the odd/odd coefficient sum of

    (1+x)^(A−1) (1+y)^(C−1) (1+x+y)^(B−1).

Applying the four signs (x,y)=(±1,±1) gives (3), including all boundary cases with 0^0 interpreted as 1. Each cycle is uniquely determined by its vertex support and its edge, and the two reversed edges distinguish the two cycles.

The S158 path count and the odd-cycle identity H(T)=J_T(2) give the further exact coefficient

    N = 2^(A+B−2) + 2^(B+C−2)
        + [2^(A+C−2)3^(B−1) − δ(−1)^(B−1)]/2.                (4)

The apparent half is always integral. For (A,B,C)=(1,2,1), J=1+6z and H=13; for (1,3,1), J=1+12z+2z² and H=33; for (2,2,2), J=1+14z+3z² and H=41. All these tournaments are strong because their two reversed-edge intervals overlap and cover the vertices.

## Sources and verification

The crossing Hamilton-path formula is from S158 [Square-sum intervals and two-edge interference](https://github.com/eliottcassidy2000/math/blob/6bb391f902df4a6fd34b947a7923f64a0d5e6057/05-knowledge/results/session158_creative_20261103/square_overlap/NOTE.md). The identity H=J(2) is the odd-cycle packing identity used in [S157](https://github.com/eliottcassidy2000/math/blob/3d6628b1804cf80457da16aa52489f016dfc99c4/05-knowledge/results/session157_gaps_20261103/cycles/README.md). This note introduces the explicit parity-filter calculations above and does not need a general theorem about interlace polynomials or delta-matroids.

Run `python3 check.py` and `python3 -O check.py`. A direct odd-cycle permutation enumerator computes N and P, and an independent endpoint-subset DP computes H, for all 56 positive triples of total order at most 9; both verify (3), the S158 formula, and H=1+2N+4P. Run `python3 top_check.py` and `python3 -O top_check.py` for nine independent general-overlap controls, including three-edge connected, nested, crossing and disconnected patterns. The cycle-packing enumerator checks (1) directly. Outputs under normal and optimized Python are byte-identical.
