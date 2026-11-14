# One Hamilton count can read a tournament's whole cover polynomial

Let B be a nonempty tournament of known order n. Write P_B(t) for the number of ordered t-path covers (empty slots allowed), and H(B)=P_B(1). For a disjoint packing of directed odd cycles of lengths at least three, regard each cycle and each uncovered vertex as a block. The ordered-cover version of the odd-cycle formula is
[
 P_B(t)=\sum_{\mathcal C}2^{|\mathcal C|}t^{n-\sum_{C\in\mathcal C}(|C|-1)}.
]
Indeed, apply the Hamilton path odd-cycle formula separately in each color and assign every resulting block to one color. In particular, P_B(0)=0, P_B is monic and has only exponents of parity n. Its coefficients are nonnegative. Also P_B(t)≤n!t^n for positive integer t: there are t^n assignments of vertices to slots, and each slot contributes at most its size factorial.

## A controllable physical readout

For k≥1 let R_k(B)=C3[1,TT_k,B], in forward cyclic order. Let A be the set of B vertices and X the set of TT_k vertices on the cycle containing the singleton. If the singleton is not in a cycle, all cycles lie in B. Otherwise A and X are nonempty, |A|+|X| is even, the weight of that cycle is 2, and its two legs are Hamiltonian paths on A and X; the latter leg is unique. This cycle visits both other blocks just once because it cannot cross between them other than by going through its singleton to return. Deleting the cycle leaves a one-way join TT_{k−|X|}→B−A, so the remaining cover counts multiply. In enforcing parity by 1+(-1)^(|A|+|X|), the signed A-sum assigns all odd blocks the weight t−1 instead of t+1. The exclusion of the empty A removes P_B(t), and the exclusion of empty X removes t^k. This proves
[
 P_{R_k(B)}(t)=t^{k+1}P_B(t)+t\{[(t+1)^k-t^k][P_B(t+1)-P_B(t)]
                         +[(t-1)^k-t^k][P_B(t-1)-P_B(t)]\}.
]
For k=1 this is the background researcher's unit-wrap identity
(L P(t)=t^2P(t)+t[P(t+1)-P(t-1)]).
For scalar count, all k give
[
 H(R_k(B))=2h+(2^k-1)(Q-h),\qquad (h,Q)=(P_B(1),P_B(2)).
]
The latter gives an exact inverse for this construction only: for a specified source, (N−2h)/(Q−h)+1 must be an integral power of 2 with positive exponent. More generally N−2h must be factored with a Mersenne factor and the remaining factor must be the physical cover difference of a source with that same h. The factor test does not supply such a source.

For a literal same-order control take the row masks (22,24,10,1,12) and (6,12,24,17,3). Both have H=15, but their cover polynomials are t^5+8t^3+6t and t^5+10t^3+4t. Their wide wraps have the counts 93·2^k−63 and 105·2^k−75 respectively. Thus a count collision can hide two channels that respond differently at *every* height. Order five is the first possible same-order split of this type: through order four there is at most one nonleading monomial of the required parity, so H fixes its coefficient. The latter tournament is also the five-vertex exterior used in the 89/91 domino in the original session.

## A single absolute count, not a sequence of measurements

**Theorem.** An explicit physical context depending only on n, with O(n³ log(n+1)) added vertices, maps B to a tournament whose single Hamilton count determines P_B. Conversely, equal cover polynomials give equal readout counts. The construction does not assume nonvanishing of the coefficients of unpadded iterated wraps.

A transitive ordinal pad of p vertices acts on the polynomial by D_p f(t)=t^p f(t). A path cannot turn back from the later block of an ordinal sum; in each color the two paths concatenate uniquely, including empty paths. The ordinary unit wrap acts by L. We build an outer functional F, starting with F(f)=f(1). Suppose the current functional is ∑_{t=1}^s c_t f(t) with every c_t nonzero. Replace it by F D_p L. For each target index j=1,...,s+1 the new coefficient is
[
 j^2 j^p c_j+(j-1)^{p+1}c_{j-1}-(j+1)^{p+1}c_{j+1},
]
omitting terms off the old support (and f(0)=0). There is at least one nonzero term at every index. At j=s+1 there is only one, so it cannot vanish. For each of the other s indices the expression is a nontrivial sum of at most three distinct positive-base exponentials in p. Such a sum has at most two real zeros: divide by its smallest exponential and apply Rolle's theorem inductively to the derivative. Therefore some p∈{0,...,2s} makes all new coefficients nonzero. Search this set deterministically. Stop when s=d=max(1,floor((n−1)/2)).

Now append an innermost D_k. The resulting physical scalar is
[
 S=\sum_{t=1}^d c_t t^k P_B(t).
]
Take B_0=n! d^n∑|c_t| and k=d·bitlength(2B_0). For m≥2 and m≤d, (m/(m−1))^d≥2, so the total contribution from indices below m has absolute value at most B_0(m−1)^k<m^k/2. Divide S by d^k and round to the nearest integer to recover c_d P_B(d), even if c_d is negative. Subtract its full contribution and repeat, obtaining all d evaluations. Interpolate P_B(t) from P_B(t)−t^n after factoring its minimal allowed power of t (t for odd n, t² for even n), using the distinct nodes 1²,...,d². The cases n=1,2 are already determined by n. The pad sum before the last pad is O(d²); the coefficient mass grows per step by at most (s²+2s)s^p, so log B_0=O(n² log(n+1)). This gives the size bound.

One can make the **gap** of a 67/69 port-reader pair encode the same polynomial: feed the compiled tournament D as module B to the published source C3[1,1,B] and use the slope-72 reader pair. Its gap is exactly 2H(D). Alternatively one last unit wrap, not followed by a pad, ensures a strongly connected *scalar-readout* tournament: simply start the outer functional with evaluation after that wrap, already with nonzero coefficients (1,1), and continue as above; for d=1 keep the two channels and decode both. These are structural distinguishability theorems, not sharp upper bounds on the size of a witness for an arbitrary target. Counting the compiled graph is not proposed as an efficient enumeration algorithm.

## What the six-reader window does and does not imply

For the published slopes 72,73,74 at each offset 67h and 69h, the six distinct counts are consecutive odds exactly when h=1 and Q=4. Every split of a Hamilton path at one of its n+1 gaps injects into the ordered two-path covers, so Q≥(n+1)h≥2h. Sorting the six counts therefore produces alternating spacings 2h and Q−2h (ties if Q=2h). In the consecutive case n≤3. A three-vertex tournament has either H=3 (a cycle) or P(t)=t³ (transitive), so it cannot meet (h,Q)=(1,4). The singleton has Q=2. Thus precisely the transitive pair gives this mechanism's 355,357,...,365. Any unbounded continuation of this particular six-function interleaving needs something other than growing this one module: change the slope set, offsets, or physical readout.

The remaining Chan–Pak obligation is an inverse construction: given N, select a realizable source profile and a reader with the required absolute output, not just a decoder for a profile that already exists. Neither the large padding nor the six-odd window proves that inverse.
