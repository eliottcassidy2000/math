# Strong tournaments with an exact partition-function transition

Let (A_R(z)) count packs of vertex-disjoint directed odd cycles in
(R), weighted by (z) per cycle. The inherited odd-cycle formula of
[Grinberg–Stanley](https://arxiv.org/abs/2307.05569) gives (H(R)=A_R(2)). The repo already connects this
polynomial to the hard-core model on the odd-cycle conflict graph
(`07-reflections/real-rootedness-of-independence-polynomial.md`);
THM-025 already refutes universal real-rootedness at order nine.
The result here is an explicit all-order *strong* family with a
computable zero distribution and pointwise pressure discontinuity.

## A general wrap identity

For any nonempty core (C), put
$$
K_C(z)=\sum_F 2^{|C|-|V(F)|}(2+z)^{|F|},
$$
where (|F|) is the number of packed odd cycles. Take (r) ordinal
copies of (C), followed by a transitive (m)-vertex tail, (m≥1);
call this whole block (B). Form the strongly connected cyclic
wrap (W=C_3[b,a,B]) with two singleton portals. Then
$$
 A_W(z)=A_C(z)^r+2^{m-1}zK_C(z)^r.
$$

A packed cycle meeting the portals must use both and is uniquely a
directed Hamiltonian path on some odd-cardinality subset (S) of
(B), closed through the portals. There is at most one such cycle.
Its contribution before summing is (zH(B[S])A_{B-S}(z)). For each
choice of (S) in the ordinal core, exactly (2^{m-1}) selections
from the transitive tail give odd total size. Tail choices affect
neither the path count nor the cycle packs. Expand (H) on (S) by
the inherited odd-cycle formula. In the resulting union pack, each
cycle has two placements: in (S), with weight 2; or outside (S),
with weight (z). Each uncovered core vertex chooses either side
independently. Both (A_C) and (K_C) multiply under ordinal sum
because no directed cycle crosses an ordinal cut. This proves the
identity, including the separate term for packs missing the portals.

### The jump occurs for every fixed core

For a fixed nonempty $C$ and $m$, and real $z>0$,
$K_C(z)>A_C(z)$ term by term (strict already for the empty pack).
Thus the limiting normalized log partition polynomial is
$\log K_C(z)/|C|$ for $z>0$, and its value at $z=0$ is zero.
The discontinuity is at least $\log 2$, since
$K_C(0)\ge 2^{|C|}$. At the boundary scale
$z=x/[2^{m-1}K_C(0)^r]$, the polynomials converge locally uniformly
in complex $x$ to $1+x$: the normalized analytic factors
$A_C(z)$ and $K_C(z)/K_C(0)$ are $1+O(z)$, uniformly on each such
compact $x$-set, so their $r$-th powers tend to 1. In particular,
Rouché's theorem locates a negative real zero asymptotic to
$-1/[2^{m-1}K_C(0)^r]$ (the unique zero in a small scaled disk
centered at $-1$, hence real by conjugation).

## A two-term strong family

For (C=C_3), (A_C(z)=1+z) and (K_C(z)=10+z). With fixed
(c=2^{m-1}) the formula becomes
$$
 P_r(z)=(1+z)^r+cz(10+z)^r,\qquad
 H(W)=3^r(1+2^{m+2r}).
$$
Its degree is (r+1). There are exactly two real roots when (r)
is odd and exactly three when (r) is even, all simple and negative.
No root is nonnegative, (-10), or (-1). On the intervals
(−∞,−10) and (−1,0) the left side of
(((z+1)/(z+10))^r=-cz) is strictly increasing whereas its right
side is strictly decreasing, with opposite endpoint inequalities.
On (−10,−1) odd (r) is impossible by sign. For even (r), set
(y=-z); the log derivative with respect to (z) of
(|(z+1)/(z+10)|^r/(-z)) is
(-9r/[(y-1)(10-y)]+1/y<0), because
((y-1)(10-y)/y ≤ 11−2√10<9). This ratio decreases from
infinity to zero. Thus this strong family has unboundedly many
nonreal zeros, not just the ordinal copies of a fixed counterexample.

## Zero distribution and pressure

Set (w=(z+1)/(z+10)). The transformed monic polynomial is
$$
 w^{r+1}-w^r-10cw+c=0.
$$
No root is (w=1); (z=-10) was not a root, so all roots correspond.
For any fixed radius in (1/10,1), Rouché's theorem applied to
(c(1-10w)) selects exactly one inner root for all sufficiently large
(r). This exceptional root tends to (1/10). No roots have
(|w|≥1+ε) for sufficiently large (r), since the
equation would require
$$
 |w|^r = c|10w-1|/|1-w| \le c(10+9/\varepsilon).
$$
Thus the other roots uniformly approach the unit circle. For every
fixed (k≥1) and (r>k), Newton's identities give the exact power
sum (p_k=1): the coefficient below the leading term is (-1)
and the next (k-1) coefficients are zero. Therefore the normalized
positive Fourier moments of the root angles tend to zero, by replacing
each root within the annulus by its radial projection and ignoring
the single exceptional root. The negative moments are their conjugates.
Trigonometric approximation now gives the uniform angle distribution.

The inverse Möbius transformation takes (w=e^(iθ)) to
$$
 z=-11/2+(9/2)i\cot(\theta/2).
$$
Consequently the empirical fugacity-zero measure converges weakly to
the Cauchy law on (Re z=-11/2), with vertical scale
(9/2). To handle the pole at (w=1), first omit an arc around 1:
its limiting mass equals its normalized angular length, which tends
to zero as the arc shrinks. The isolated inner root gives a separate
negative zero tending exponentially to zero. Indeed, once (|w|<1/5),
(10w-1=O(5^{-r})), hence (z=O(5^{-r})), and substitution into
(z=-w^r/c) sharpens this to (z∼-10^{-r}/c).

For real (z>0), the second term dominates and the pressure
(log P_r(z)/(3r+m+2)) tends to ((1/3) log(10+z)). At
(z=0) it is exactly zero for all (r). Thus the pointwise limit
has an upward jump at zero, despite the regularity of every finite
polynomial there. More quantitatively, for $z_r=e^{-\alpha r}$
with $\alpha>0$, the pressure tends to
$\max(0,\log 10-\alpha)/3$. This follows because the two terms
have exponential rates 0 and $\log 10-\alpha$, respectively.
At the finer scale $z=x/(c10^r)$, $P_r(z)\to1+x$ uniformly on
compact sets of $x$, by applying $(1+O(10^{-r}))^r\to1$
to both terms. A *single* packed crossing cycle can have exponentially
many realizations, so the origin is not an ordinary uniform
small-fugacity regime. The mechanism is an explicit illustration of the
[Yang–Lee partition-function zero approach](https://doi.org/10.1103/PhysRev.87.404), but uses
a nonlocal wrap and says nothing about a fixed-degree lattice model
or a mass gap. At (z=2) the exact scalar count above is only a thin
family, not a coverage theorem for Chan–Pak.
