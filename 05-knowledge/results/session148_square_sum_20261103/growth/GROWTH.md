# Exact growth of square-sum graphs

Let G_n have vertices 1,...,n and an edge {i,j} precisely when i != j and i+j is a square. Let d(n) be the number of edges introduced when vertex n is added. The floor formula d(n)=floor(sqrt(2n−1))−floor(sqrt(n)) is inherited from [Session 147's core audit](https://github.com/eliottcassidy2000/math/blob/50f4e1fa/05-knowledge/results/session147_square_sum_20261103/core_audit/AUDIT.md), along with the G_15–G_18 finite facts. The explicit threshold blocks, cancellation classification, sharpened error term, and path-cover identity below are deductions from these definitions, independently tested here. No general Hamiltonicity claim follows.

## 1. Exact degree thresholds, in blocks

For integers n,k >= 1 put a=floor(sqrt(n)). Then

    d(n) >= k  iff  max(a², floor((a+k)²/2)+1) <= n <= (a+1)²−1.

Indeed, the first candidate square strictly greater than n is (a+1)², and the kth is (a+k)². It is strictly less than 2n if and only if n >= floor((a+k)²/2)+1. This gives the entire block, including its two endpoints.

The block is nonempty exactly when

    (a+k)² <= 2(a+1)²−3,
    equivalently a²+(4−2k)a−k²−1 >= 0,
    equivalently a >= a_k := k−2+ceil(sqrt(2(k−1)²+3)).

Here a is a nonnegative integer, k >= 1, and the quadratic's other root is negative. Therefore the first such vertex is exactly

    n_k = max(a_k², floor((a_k+k)²/2)+1).

For k=1,2,3,4,5,6 these start at 3,13,33,61,99,163. The threshold describes first occurrence, not persistence: d(3)=1 but d(4)=0. In fact every n >= 5 has d(n)>=1, a separately established fact. The list of candidate neighbors is exactly {(a+j)²−n : 1 <= j <= d(n)}. These are distinct vertices in [1,n−1], so no diagonal or double counting occurs.

## 2. Birth/death recurrence and Pell cancellation

For n >= 2,

    d(n)−d(n−1)
      = 1_{2n−2 or 2n−1 is a square} − 1_{n is a square}.

The two possible birth squares cannot both occur, since consecutive integers cannot both be positive squares. The new upper endpoint of the interval of squares is 2n−1; the lower endpoint excludes n itself. At a square n=r², a birth and death cancel precisely in one of the following mutually exclusive cases:

* The new square has even root 2t and r²−2t²=1, with r>=2 and t>=1. Its solutions are r+t sqrt(2)=(3+2 sqrt(2))^j for j>=1. First cancellation vertices: n=9,289,9801,... .
* The new square has odd root s and s²−2r²=−1, with s>=3 and r>=2. Its solutions are s+r sqrt(2)=(1+sqrt(2))^(2j+1) for j>=1. First cancellation vertices: n=25,841,28561,... .

This uses the standard complete solution of the Pell equations. For example d(8)=d(9)=d(10)=1: at n=9 the newly admitted square 16 and the excluded square 9 cancel. Here 8=9−1 and 10=9+1 are vertices; this arithmetic observation does not make both endpoints of the G_15 Hamilton path. In fact its endpoints are 8 and 9.

## 3. Exact cumulative edge count and sharper asymptotic

Put K=floor(sqrt(N)) and M=floor(sqrt(2N−1)). Then

    |E(G_N)| = (M−K)N−K + K(K+1)(2K+1)/6
               − M(M+1)(2M+1)/12 + ceil(M/2)/2.

This expression is always integral. To derive it, sum d(n). A square m² contributes to floor(sqrt(2n−1)) once n >= floor(m²/2)+1, and to floor(sqrt(n)) once n >= m². Thus

    |E(G_N)| = sum_{m=1}^M (N−floor(m²/2))
               − sum_{m=1}^K (N−m²+1).

Use sum of squares and sum_{m=1}^M floor(m²/2) = [sum m²−ceil(M/2)]/2. In particular, expanding the resulting polynomials with K=sqrt(N)+O(1) and M=sqrt(2N)+O(1), and retaining the quadratic terms, gives

    |E(G_N)| = (2/3)(sqrt(2)−1) N^(3/2) + O(sqrt(N)).

For this cancellation the O(1) errors in K and M cannot simply be discarded termwise. Set K=sqrt(N)−u and M=sqrt(2N)−v with u=O(1) and v=O(1); the linear-in-u and linear-in-v terms in each cubic polynomial cancel. Both summands have a −N/2 term, which cancels too. This strengthens the previously stated O(N) remainder.

## 4. Exact Hamilton-path recurrence, with the missing case included

Let H_n be the number of ordered Hamilton paths in G_n, so reverse orientations count separately when n>1. Let w=n+1. Define

    A_n(w) = sum_{P=(v_1,...,v_n) Hamilton in G_n}
             (1_{w~v_1} + 1_{w~v_n}).

Let C_n(w) count all ordered pairs (P,Q) of nonempty directed simple paths in G_n whose vertex sets partition [n], with the last vertex of P and the first vertex of Q both adjacent to w. Then exactly

    H_(n+1) = A_n(w) + C_n(w).

Proof: delete w from an ordered Hamilton path. If w was first or last, the result is one Hamilton path with an eligible endpoint. If w was internal, the result is exactly one ordered pair (P,Q) of the stated type. Conversely, insert w at that endpoint or between P and Q. These constructions are bijections and their cases are disjoint. The ends of P and Q need not be adjacent in G_n, so counting only insertions into existing Hamilton paths omits valid paths. C_n depends on two-path covers, which is a substantive state; the identity does not establish that H_n>0.

For independent finite checks, H_14=0, H_15=H_16=H_17=2, H_18=0. In fact G_14 has leaves 8, 9 and 10, so H_14=0; but C_14(15)=2 and H_15=2. The ordered G_15 path shown below splits at 15 into paths (8,1) and (10,6,3,13,12,4,5,11,14,2,7,9). Their adjoining vertices 1 and 10 are not adjacent. This concrete example proves endpoint-only extension can miss all Hamilton paths at the next stage. Up to reversal, the G_15 path is

    8,1,15,10,6,3,13,12,4,5,11,14,2,7,9.

The graph has 15 edges and the path uses 14. Its sole unused undirected edge is {1,3}. G_18 has three degree-one vertices, so no Hamilton path. More generally, if a Hamilton path exists, it uses n−1 of E undirected edges, leaving E−(n−1) unused. A chosen orientation of every graph edge has E arcs, while the symmetric digraph has 2E arcs. An oriented Hamilton path uses n−1 arcs; its remaining arc counts are respectively E−(n−1) and 2E−(n−1). Along any unbounded sequence of n for which a Hamilton path exists, the used fraction is

    (n−1)/|E(G_n)| ~ 3/[2(sqrt(2)−1)sqrt(n)],

so the unused fraction tends to 1. This conditional statement does not assume Hamilton paths exist for all sufficiently large n.

Moreover G_14 is a connected tree (its 13 edges can be enumerated directly); every later vertex has at least one earlier neighbor, so G_n is connected for n>=14. Thus E−(n−1) is also the cycle-space dimension over F_2 for any n>=14, even when G_n has no Hamilton path.

A tournament is different again: it has exactly one directed arc on every pair, hence n(n−1)/2 arcs, and is not an orientation of sparse G_n unless G_n is complete. Reverses of used path edges count among unused arcs in the symmetric digraph.

## 5. Verification and scope

Run `python check_growth.py` and `python -O check_growth.py`. The verifier explicitly checks d(n) by enumerating actual candidate neighbors, not by assuming the floor formula; checks the edge formula against cumulative enumerated degrees and degree blocks for n<=10000; checks first occurrences for k<=100; checks the birth/death and Pell identities for n<=10000; and computes the two-path-cover recurrence from directed simple-path dynamic programming through n=18. All failures use explicit exceptions rather than `assert`, so optimization does not disable them. These tests corroborate but do not replace the proofs. Finite checks do not establish Hamilton paths for all sufficiently large n.
