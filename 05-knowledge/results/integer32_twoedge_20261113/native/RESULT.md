# Native metagraph as a marked affine cut system

Status: **PROVED IN NOTES after independent audit.**
The claims below concern the native nonpath-complement metagraph of
tournaments, with the dart and converse conventions of
`incidence_closure_20261113/closure/check.py`. They give an exact presentation,
not a new scalar A000568 recurrence or an all-order connectivity proof.

## 1. Coordinates and legal transport

Let n>=3, E=E(K_n), and x in F2^E encode reversal of the reference arc
i->j for i<j. For each pair of incident edges e={v,u}, f={v,w}, let
beta(v;u,w) be one precisely when the two arcs form a directed path.
These are the two-slot wedge observations recovered by Session320.
Write delta for the edge coboundary on the line graph L(K_n), so
`(delta x)_{ef}=x_e+x_f`. Then

    beta = beta_reference + delta x.

The image and its local star/triangle checks are proved in the companion
`reconstruction/PROOF.md`. Since L(K_n) is connected, its kernel is
{0,1}: beta remembers exactly a labelled tournament up to converse.

If p is a directed Hamilton path of the chosen orientation T, with
undirected edge set P, the native move reverses every nonpath arc. Thus

    x' = x + 1 + chi_P,
    beta' = beta + delta chi_P.                         (1)

The disappearance of the global 1 in (1) is specifically a converse
quotient, not permission to discard a marked path's orientation.

Every such move changes exactly

    sum_v d_P(v)(n-1-d_P(v)) = 2(n-2)^2                 (2)

wedge bits. At the two endpoints d_P=1, and at the n-2 interior vertices
d_P=2. A wedge is counted once at its unique common vertex. This proves
(2), including n=3. In particular the labelled beta-state graph has no
native loops for n>=3. Loops can appear after vertex relabelling.

## 2. Exact adjacency and the exceptional order four

For two realizable fields beta,beta', solve

    delta s = beta + beta'.

There are exactly two solutions s and s+1. The states are natively
adjacent iff at least one support S of these two masks is the edge set
of a spanning undirected path whose edges consistently orient into a
directed Hamilton path in a representative T of beta. The property is
unchanged on replacing T by its converse.

Necessity is (1). For sufficiency, use the oriented path in the native
move; (1) gives exactly beta'. Thus the guard is an equivalence. To test
the underlying path, require connectivity as well as two degree-one
vertices and n-2 degree-two vertices. The latter degrees alone admit a
triangle plus a disjoint edge at n=5. Even a genuine spanning path must
still pass the orientation guard.

Two different path supports P,Q give the same jump iff P and Q are
complements: delta(chi_P+chi_Q)=0. As both have n-1 edges, this requires
`2(n-1)=binom(n,2)`, hence n=4. At order four this ambiguity really
occurs: direct the edges of

    0 -> 1 -> 2 -> 3,       2 -> 0 -> 3 -> 1.

These complementary paths orient all six edges consistently as a
tournament, and both have the same beta jump. For n!=4, the jump
determines the path support uniquely; for n=4 there are at most two.
This is uniqueness at fixed labels, before quotienting target states by
permutation. Many different labelled jumps may enter one unlabelled
target class.

## 3. The extra phase that preserves native dart weights

The correct marked carrier is

    D_n = {(beta,p): p is an ordered Hamilton path in either of the
                      two converse representatives reconstructed by beta}.

For n>=2 the ordered p selects exactly one of those representatives.
Consequently `(beta,p) <-> (T,p)` is an equivariant bijection with the
ordinary labelled tournament-plus-directed-path carrier. The action of
S_n on D_n is free: a permutation fixing an ordered list of all vertices
fixes every vertex. Relabel p as (0,...,n-1); its consecutive arcs are
fixed and all other arcs are free. Therefore

    |D_n / S_n| = 2^((n-1)(n-2)/2).

This is exactly the repository's fixed-spine word set. Equation (1),
with p retained, defines its native involution N. The commuting
reflection is rho(beta,p)=(beta,reverse(p)); it selects the converse
orientation. On fixed-spine words it is order reversal plus arc
reversal. Fixed points here are on D_n/S_n; rho has no fixed ordered
list on D_n itself. Its fixed words are the blue darts; its other orbits are the
pairs of black darts. Do not replace D_n by unordered paths before
forming the dart census: this would quotient by rho and halve every
nonfixed pair.

The number of reflection-fixed words is
`2^floor((n-1)^2/4)`. Indeed reflection pairs the free edge coordinates
{i,j} and {n-1-j,n-1-i}; the number of their orbits is that exponent.
Hence an early reflection quotient would count

    (2^((n-1)(n-2)/2) + 2^floor((n-1)^2/4))/2

objects. At n=4 it gives 6 instead of 8 darts; at n=5 it gives 40
instead of 64. This is a source-weight error, not a new metagraph.

The node projection sends beta to its S_n orbit, namely a merged
tournament class. Counting the ordered marked carrier gives the
inherited integer matrices A and blue B; black K=A-B, with loops
counted by their two native darts. If T has H directed Hamilton paths,
automorphism order a, and its converse fibre has m=1 or 2 ordinary
classes, the merged source degree is mH/a. These weights follow from
the free ordered marking, and must survive every subsequent quotient.

## 4. Scope and relation to the other lanes

The representation supplies all nodes, a realizability guard, legal
transport, multiplicities and converse phase for every n. It does not
compress their enumeration to earlier A000568 values. The minimum
distance 2n-4 of the full wedge code differs from the native-jump
weight 2(n-2)^2: the latter restricts the cut support to Hamilton paths.
No connectivity result for the entire cut code implies connectivity
of the native graph or its positive-black-degree core, because only
the guarded cuts in section2 are available.

The F64 norm example in the census sidecar forgets scalar phase under
one operator and recovers it under a second. Here forgetting ordered
path phase immediately changes the integer transport weights. This is
a shared quotient-closure question, not an identification of the two
state spaces or their transition operators.

Independent finite checks are in `check.py`: complete ordinary and
merged native/blue matrices for n=3..6, compared entrywise to the
previous immutable census; every Hamilton path of every ordinary
representative through6; exhaustive candidate target beta states
through5; the n=4 double-support example; and the degree-only and
orientation-only hostile controls. Normal and optimized execution
must agree. The all-order statements above require their proof audit
in addition to these finite checks.
