# Hamilton-path cuts force size-dependent amplification of a homogeneous module

**Status: PROVED incidence inequality and its context corollaries; FINITE-EXACT implementation controls.** No literature-priority claim is made. The unrestricted claim that every positive odd Hamilton-path count except 7 and 21 occurs remains OPEN. This separate note does not modify the frozen Session142 binary-language packet.

## 1. The incidence bound

Let T be a tournament on n>=1 vertices, H(T) its number of directed Hamiltonian vertex words, and p_k(T) its number of covers by k unordered, nonempty, vertex-disjoint directed paths. The paths cover all vertices. Their internal orders are retained; only the order of the k paths is forgotten. Put p_k=0 for k>n.

**Theorem.** For 1<=k<=n,

    k! p_k(T) >= H(T) binom(n-1,k-1).                    (1)

Equivalently p_k(T)>=ceil(H(T) binom(n-1,k-1)/k!).

**Proof.** Cut a Hamiltonian path at k-1 of its n-1 gaps. There are exactly H(T) binom(n-1,k-1) such incidences. Forget the order of the resulting k nonempty runs. A fixed unordered path cover has at most k! preimages: choose an order of its k paths; their concatenation and the cut positions are then uniquely determined. Some orders fail the directed-edge conditions, which only reduces the number of preimages. Summing the fibre bound proves (1). No claim that every path cover concatenates is needed. QED.

In particular

    2 p_2(T) >= (n-1)H(T).                              (2)

Since H(T) is odd, when n is even the right side is odd while the left side is even, giving

    2 p_2(T) >= (n-1)H(T)+1.                            (3)

For n=1, use p_2=0 separately. The parity fact is the classical tournament Hamilton-path parity theorem; no new parity proof is asserted.

## 2. The inherited homogeneous run kernel

A one-port context is a tournament Q with distinguished vertex q. The tournament Q[T] replaces q by T; every exterior vertex has the same orientation to every vertex of T as it had to q. Let W_k(Q,q) count directed quotient words with q appearing k times, each exterior vertex appearing once, and no adjacent appearances of q. Such words record the maximal T-runs of a Hamilton path. All W_k are nonnegative integers, W_1=H(Q), and W_k=0 when k is too large.

The inherited maximal-run substitution kernel is

    H(Q[T]) = sum_(k>=1) k! W_k(Q,q) p_k(T).             (4)

For each quotient word, assign the k unordered internal paths to its k ordered q-runs in k! ways. This is a bijection, not a formal realization surrogate. Formula (4) is credited to THM-3121 and the Session136 gluing packet, whose copies and hashes are pinned in SOURCE_PINS.json.

Combining (1) and (4) gives the useful whole-kernel bound

    H(Q[T])/H(T) >= sum_(k=1)^n W_k(Q,q) binom(n-1,k-1). (5)

## 3. Active versus inactive ports

If q is a singleton strongly connected component of Q, every Hamilton path traverses the strong components in dominance order and visits the module only once. Hence

    H(Q[T]) = H(Q-q) H(T),                              (6)

where H(empty)=1. The multiplier H(Q-q) is odd. This is the multiplicative, not the amplifying, case.

Otherwise q lies on a directed triangle q->u->v->q. To see the triangle, take a directed cycle through q and follow it from an out-neighbor of q until its first in-neighbor of q; the transition gives such a triangle. Replace q by two transitive twins q1->q2. Start with the path q1,u,v,q2 and insert all other quotient vertices one at a time. The tournament insertion lemma preserves the old vertex subsequence: a new vertex can be placed before the first vertex, after the last, or at an edge where its relations switch from incoming to outgoing. Consequently the two twins stay separated. Collapsing their labels yields a word counted by W_2. Thus

    W_2(Q,q)>=1.

Also H(Q)>=3: a directed triangle has three Hamilton paths, and successive tournament path insertion extends each while preserving its restriction to the original triangle, so these extensions stay distinct.

The singleton-SCC condition and the cyclic condition exhaust the cases. A vertex in a nontrivial strong component lies on a directed cycle; conversely any directed cycle lies in one strong component.

**Corollary (uniform active-port amplification).** For every nonempty T of order n and every context whose port is cyclic,

    H(Q[T]) >= 3H(T)+2p_2(T) >= (n+2)H(T).              (7)

If n is even, the stronger bound is

    H(Q[T]) >= (n+2)H(T)+1.                             (8)

The context may depend on T and may have arbitrarily many vertices. No finite-template hypothesis is used. The directed triangle context realizes the first expression exactly:

    H(C3[T,K1,K1]) = 3H(T)+2p_2(T).                    (9)

For T=C3, n=3,H=3,p_2=3, (7) is equality: H(C3[C3,1,1])=15. For T the two-vertex transitive tournament, (8) is equality: H(C3[L2,1,1])=5. No claim that equality holds at all orders is made.

## 4. Consequence for construction strategies

Suppose cyclic homogeneous-module moves are required to have a bounded count multiplier H(Q[T])/H(T)<=B. Then (7) forces n<=B-2. Thus no such construction acts on a family of source tournaments with unbounded order, even if the context is chosen separately for every source.

In particular a fixed affine radix H_new=b H+d, with b fixed and d in a bounded finite digit set, has bounded multiplier for positive H and cannot be implemented by cyclic homogeneous one-port substitution along an unbounded-order family. Unbounded Hamilton-path counts necessarily require unbounded order, since H(T)<=n!.

Inactive ports remain exactly the odd-multiplication case (6); this note does not turn an arbitrary product into an affine digit. For the specific binary target 2H+-1 and H>1, neither case works: the cyclic case already gives at least3H, while in the inactive case divisibility by H would force H to divide1.

This is **not** an obstruction to the nonhomogeneous S0/R1 source operations in the separate binary-language note. Their new vertex loses to the old parent A but beats its old marked vertex x, so it does not see the whole old tournament A+x homogeneously. Multi-port surgery, changing old edges, and other actual constructions are also outside (7).

The earlier Session136 active-port floor 3H+2p_2>=2^n+1 remains valid. The incidence bound adds dependence on the source H, giving the new size-times-count control (7). Neither bound proves a new missing integer, a cofinite spectrum, or all-integer realization.

## 5. Verification

`verify.py` independently computes induced-subset Hamilton counts, unordered path-cover counts by least-vertex set-partition recursion, and actual C3 substitutions. Its default finite universe is every labelled tournament of orders1..5, all k, plus every pointed context of orders1..4 for the active-port/word-coefficient checks. The selector `--max-order N` permits a smaller bounded replay (1<=N<=5). Literal permutation counts give an independent Hamilton-count control through order4. Normal, -O, and relocated default receipts agree. These are finite controls of the implementation; (1), (4), and the insertion proof establish the universal claims.
