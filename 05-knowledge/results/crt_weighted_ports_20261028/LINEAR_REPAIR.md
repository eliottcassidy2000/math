# Linear marking repair: concurrent agreement and an exact implementation

**PROVED scoped; same-agent alternate implementation.** The same bound arrived
concurrently in Session94's [linear protocol](../session94_recursive_ports_crt_clocks_20261029/loop_protocol/LINEAR_LOOP_PROTOCOL.md),
read at incoming commit733df0d442. This note credits that overlap and records
our recursive-factor implementation. There is no priority claim.

Let W_k have vertices F3^k, with tournament orientation determined at the
first differing digit by0→1→2→0. Let alpha negate digits, A_k=Aut(W_k), and
C_k=C_(A_k)(alpha). Inherited Session93 proves

    A_k=A_(k-1)^3 ⋊ C3,       C_k ≅ C_(k-1) × A_(k-1).

The paired embedding is psi(b)=(1,b,beta b beta), beta=alpha_(k-1).
Every accumulated based marking error belongs to C_k. To repair c, realize
c^-1 by the same construction below.

## The bound and proof

Every c in C_k is the holonomy of at most

    c_k=(3^k-2k-1)/4

triangles, hence at most3c_k legal square-root transport edges. For k=0,1
the group is trivial. With N=3^k, this is O(N), improving the previous
O(N log N) bound in [FULL_HOLONOMY](../lucas_prime_reflection_towers_20261028/FULL_HOLONOMY.md).

All permutations compose right to left. Write g(t,s)=sqrt(ts), the unique
odd-order square root, and h(u,t,alpha)=g(u,alpha)^-1 g(u,t)g(t,alpha).
Session93 proves that for a canonical zero-prefix child rotation v,

    x=root rotation, y=(v,1,1), t=x²alpha, u=y²alpha
    imply h(u,t,alpha)=psi(v).

This is a literal triangle of anti-involutions, not a formal generator with
uncharged execution. The calculation uses (y²x)^3=(v²,v²,v²), so its square
root is its fifth power.

For any child prefix p, choose a child automorphism a_p mapping the all-zero
prefix of the same length to p, with unchanged suffix. Then
psi(a_p) centralizes alpha and conjugates the canonical triangle to one
whose holonomy is psi(v_p). The base mark stays alpha. Reversing the triangle
gives its inverse at the same cost. The conjugator is used to construct
endpoints; it is not an uncharged path being secretly executed.

Every A_j element has one local C3 coordinate at each of the
(3^j-1)/2 internal tree nodes. Ordered recursive decomposition realizes it
with at most that many local rotations, using exponents+1 or-1. Decompose
C_k into the factors A_0,...,A_(k-1), conjugate the paired triangles as above,
and extend zero-child triangles by alpha on the other blocks. Products of
marks cancel there, so their transports are identity outside that child.
The number of triangles is at most

    sum_(j=0)^(k-1) (3^j-1)/2 = (3^k-2k-1)/4.

The decomposition reconstructs the full permutation, so every named port
is retained, including simultaneous ports. QED.

| k | vertices N | free ternary correction coordinates | edge upper bound |
|---:|---:|---:|---:|
|2|9|1|3|
|3|27|5|15|
|4|81|18|54|
|5|243|58|174|

At27 vertices, |C_3|=3^5=243. Session93's **minimum number of group generators**
is instead3, and its **minimum number of individual vertex anchors** is4.
These count different resources. The five free trits explain the243 possible
errors; they do not become three trits merely because three generators suffice.

## Verification and scope

`linear_repair.py` checks all3 corrections at9 vertices and all243 at27,
then24 seeded corrections each at81 and243. Every canonical triangle type
through243 is checked against all tournament arcs; actual roots are computed
cyclewise, returned marks are involutions, and the resulting full permutation
equals the requested correction. There are769762 always-active checks.

The all-size proof uses the recursive decomposition, not extrapolation from
these controls. This counts square-root **re-marking** edges, not native
Hamilton-path complement moves. It gives no black-core connectivity theorem.
Full permutation endpoints may cost O(N) to store and process per edge;
compact prefix descriptions need their own implementation cost model.
Original source weights are pushed along the actual permutations, never
reset to a uniform law. Physical arithmetic clocks are separate from this
explicit marking clock.
