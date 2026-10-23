# Square clocks, two-colored twins, and exact path reconstruction

## Carrier and statement

Let Q_N be the simple graph on {1,...,N} with an edge ab exactly when a+b is a square. A path is an ordered vertex sequence; its **square clock** is the ordered sequence of square roots of adjacent sums. This is not an unordered edge inventory, nor a clock modulo a fixed modulus.

The results below hold more generally for adjacent-sum words of permutations of {1,...,N}. For square clocks, squaring each entry gives that word.

**Theorem.** An ordered adjacent-sum word has at most one Hamiltonian realization when N is odd and at most two when N is even. If two exist, their difference is governed by a divisor d of N/2 and the exact two-coloring

A_d = union_{j=0}^{N/(2d)-1} {2dj+1,...,2dj+d},
B_d = {1,...,N} minus A_d.

The first path alternates A_d and B_d; its twin replaces each A_d vertex v by v+d and each B_d vertex v by v-d. Conversely, any Hamiltonian path alternating this partition has a distinct twin with the same adjacent-sum word. No Hamiltonian realizations are lost by this classification.

### Proof: reconstruction and multiplicity

Write the prescribed sums as s_1,...,s_{N-1}. Define c_0=0 and c_i=s_i-c_{i-1}. Every realization necessarily has

x_i=(-1)^i t+c_i, 0<=i<N,

where t=x_0. Conversely this formula has the prescribed sums.

For odd N, summing the vertices gives t+sum c_i=N(N+1)/2, so at most one t is possible. For even N, summing squared vertices gives

N t^2+2t sum_i (-1)^i c_i+sum_i c_i^2=N(N+1)(2N+1)/6.

This is a nonzero quadratic, so at most two real, hence at most two integer, starts exist. These moment equations are necessary; membership, positivity and distinctness still must be checked.

If two starts t and t+d exist with d>0, vertices at even positions form a set A and shift by +d, while vertices at odd positions form B and shift by -d. Both A,B and A+d,B-d partition {1,...,N}. Write A(z),B(z) for their generating polynomials and P(z)=z+...+z^N. Then

A+B=P, z^d A+z^(-d)B=P.

Subtract after multiplying the second identity by z^d to obtain

(z^(2d)-1)A=(z^d-1)P,

hence (1+z^d)A=P. A primitive 2d-th root of unity is not 1 and is a root of P, so 2d divides N. Polynomial division now gives exactly the displayed alternating blocks A_d. This argument also proves that the partition is forced, not fitted to the example.

Conversely the map tau_d swapping the two d-long halves of every 2d-block is a permutation. Across A_d,B_d its two shifts cancel, so every adjacent sum is preserved. It produces a distinct path. This proves all assertions.

## The square-sum threshold is exactly 36

The attached exact certificate exhausts even N=2,...,34 and all d with 2d dividing N. In the bipartite subgraph of Q_N retaining just A_d--B_d edges, every such candidate has an isolated vertex, more than two degree-one vertices, or is disconnected. These are rigorous Hamiltonian-path obstructions. Odd N cannot have twins by the theorem. At N=1 the sole path is unique.

At N=36,d=6 there is a Hamiltonian path

3,33,16,9,27,22,14,11,25,24,1,35,29,20,5,31,18,7,2,23,26,10,6,19,30,34,15,21,4,12,13,36,28,8,17,32.

Its twin under tau_6 is

9,27,22,3,33,16,20,5,31,18,7,29,35,14,11,25,24,1,8,17,32,4,12,13,36,28,21,15,10,6,19,30,34,2,23,26.

Both use every vertex once and have identical ordered square clocks. The checker verifies every edge and each lower-order certificate. Thus 36 is the smallest square-sum order admitting this particular information-loss phenomenon. This is NOT the first order admitting a Hamiltonian path, and is not a claim that every order >=36 has twins.

## Arbitrary words: interval minus forbidden points

For a word of k sums and a finite host [1,N], the same affine recursion describes all candidate vertex sequences. A walk in a simple graph must also exclude consecutive equal vertices; the simple-path collision test below enforces that and all other repetitions. Each even i gives 1-c_i<=t<=N-c_i; each odd i gives c_i-N<=t<=c_i-1. Intersect these integer intervals to get [L,U].

For vertex-distinct paths, two indices of the same parity collide for every t iff their c_i agree. If this occurs there is no simple realization. Otherwise each opposite-parity pair excludes at most one integer start, obtained from

((−1)^i−(−1)^j)t=c_j−c_i.

Let E be the finite set of integral excluded starts in [L,U]. The exact start enumerator is sum_{t=L}^U z^t minus sum_{t in E}z^t, with zero polynomial if L>U or a permanent collision occurs. Repeated excluded values must be counted only once. This is an exact finite reconstruction formula, not a moment approximation.

A marked vertex at any declared position fixes t outright. An unindexed mark need not. Multiple paths sharing a declared port yield equations of the form epsilon_a t_a−epsilon_b t_b=c, with epsilon in {+1,−1}; bounds and collision exclusions must remain attached. Around an even reflection word closure requires its alternating translation to vanish; around an odd word closure fixes twice the initial vertex. These are signed affine constraints, not a license to discard port positions.

## What a periodic coloring cannot do

For every modulus m>=1, vertices a=m^2 t^2 and b=3m^2 t^2 are distinct, have the same residue modulo m, and satisfy a+b=(2mt)^2. Thus any coloring depending only on n mod m has monochromatic edges at arbitrarily large sizes. The same holds for eventually periodic colorings, choosing t large enough.

This does not make periodic colorings useless: they can encode allegiance, prune choices or distinguish some edge types. It rules out their being proper colorings of the entire infinite square-sum graph. The exact block coloring A_d,B_d above is likewise a proper coloring only of the restricted cross-block graph, not of every square-sum edge. Its mathematical role is to classify twins of the ordered square clock.

## Relation to the session

The square clock is a concrete example of a projection with an exact, sharply bounded fibre. Its lost datum is one start when N is even; the twin displacement imposes a structured two-coloring. This makes the user's even/odd and coloring analogy testable without identifying it with tournament parity, Bott periodicity, or Berggren ancestry. Common marked ports restore the affine starts; separately chosen local starts can be incompatible.

No Goldbach, Collatz, all-size tournament-wiring or universal square-coloring claim follows. All lower-order minimality is certified by the finite divisor list, while multiplicity and block classification are all-size proofs.
