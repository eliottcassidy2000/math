# Signed operation spectra and paired oddness

**PROVED scoped deductions; classical and repository dependencies identified below.
Same-agent audit only.** Local labels S1-S2 and P1-P3 belong to this file.

## S1. The relative signed support spectrum

Let D be a fixed oriented simple graph on V. Let A be any finite subgroup of
Aut(D), and let Omega be its free unordered pairs (nonedges of D), of size d.
Choose a reference direction on each pair. Completions of D form the cube
Q = F2^Omega, on which g acts by x -> P_g x + b_g. This is a signed coordinate
permutation action, not generally a linear permutation action.

For each cycle C of P_g, put ell(C)=|C| and
epsilon(C)=(-1)^(sum of b_g on C). Define

    F_D,A(z) = (1/|A|) sum_g product_C (1 + epsilon(C) z^ell(C))
             = sum_k f_k z^k.

A support S subset Omega is **relatively sign-even** if the orientation sign on S
is trivial on its stabilizer A_S. This uses A_S, not the full automorphism group
of the abstract graph S. Then:

1. f_k is the number of relatively sign-even A-orbits of k-element supports.
2. F(1)=|Q/A|.
3. F(-1) is the number of Q/A classes fixed by kappa(x)=x+1.
4. The single-free-edge-flip operator B, defined by
   (Bf)(x)=sum_e f(x+e), on A-invariant real functions has eigenvalue d-2k
   with multiplicity f_k.

**Proof.** The Walsh characters chi_S(x)=(-1)^(sum_{e in S} x_e) form a basis.
The action permutes their lines, with the stabilizer acting on each line by
the orientation sign. One invariant vector survives from a support orbit
exactly when that sign is trivial. On a g-stable support, every coordinate
cycle is either wholly present or absent, giving the product in the formula
by averaging the degree-wise trace. Both kappa and B commute with A.
They act on chi_S by (-1)^|S| and d-2|S|, respectively. The invariant function
space also has the orbit-indicator basis, so its dimension and the trace of
kappa give (2) and (3). This proves all four statements. QED.

Changing reference directions conjugates the action by a translation and
changes Walsh basis signs, leaving F invariant. Thus this is a gauge-independent
operation invariant. The weighted quotient matrix for B counts how many
individual coordinate flips lead to each orbit; loops and multiplicities count.
It is not the simple adjacency matrix obtained by forgetting those weights.

For D empty and A=S_n, relative sign-evenness is precisely the sign-even graph
definition used by [Royle et al.](https://arxiv.org/abs/2204.01947).
The Walsh mechanism is inherited from THM-4721
[twisted Burnside / components / Walsh basis](../../../01-canon/theorems/THM-4721-twisted-burnside-on-odd-fixed-point-sets-sylow-refinement-of-even-graphs-and-tournaments-components-walsh-basis.md).
The converse trace for ordinary tournaments is also inherited from THM-4531
[twisted Mallows-Sloane principle](../../../01-canon/theorems/THM-4531-twisted-mallows-sloane-principle-and-no-natural-bijection-for-switching-classes.md).
The present extension keeps a general oriented carrier and its flip operator
explicit. It is an application of standard signed permutation representations,
not a new general Fourier theorem.

**Complement boundary.** For nonempty D, kappa reverses only free pairs and keeps
D fixed. It is not global tournament converse. Any transport to native path
complement, graph complement, or global converse needs a separate commuting
diagram. At D empty it is global converse.

**Why E=V-1 enters.** If |V|=N and |E(D)|=N-1, then
d=binom(N-1,2), the dimension of tournaments on an (N-1)-set X.
That dimension equality supplies no action isomorphism.
A coordinate-permuting affine (monomial) A-equivariant map would preserve B
and kappa, so would force equality of F. A general affine map may mix flip
directions, and a weak A-set bijection need preserve neither operator.

The exact monomial criterion for full A=Aut(D) is already
THM-4835 [Euler criterion](../../../01-canon/theorems/THM-4835-euler-characteristic-criterion-for-marked-excess-minus-one-carriers-in-tournaments.md):
chi(G^K)=1 for every K<=A, equivalently V is isomorphic to E(D) plus one point
as an A-set. Trees of every orientation and functional orientations satisfy it.
This includes directed cycles plus an isolate: disconnected carriers are essential.
Equal F is only a necessary test here; no converse classification by spectra
is claimed.

## S2. An all-size flip-spectrum obstruction inside a weak trijection

Let B_k, k>=2, consist of two internally disjoint directed paths of length k
from u to w, and one isolated vertex. It has 2k+1 vertices, 2k arcs, and
d=k(2k-1) free pairs. Its full automorphism group A=C2 swaps the paths.
The swap has three fixed vertices and no fixed carrier edges.
On free coordinates it has three positive singleton cycles, k-1 negative
singleton cycles, and only cycles of length two otherwise. Therefore

    f_1(B_k) = (d+3-(k-1))/2 = k(k-1)+2.

An A-set X of size 2k has r transpositions and f=2k-2r fixed vertices.
If r>=1, its tournament action is free: a tournament cannot be fixed by a
nonidentity involution, since the arc on a swapped vertex pair would reverse.
Thus both completion spaces have 2^(d-1) A-orbits, and both are free A-sets
of size 2^d. They are weakly A-equivariantly bijective.

But the target's first-degree invariant dimension is

    f_1(Tour(X)) = (d+binom(f,2)-r)/2
                 = k(k-1)+(k-r)^2.

Equality with B_k would require (k-r)^2=2, impossible. For r=0 the action
is trivial and the dimension is d, also different. Consequently **no choice
of X can preserve the weighted single-edge-flip spectrum**, even though
every nontrivial C2 action on X gives a weak trijection. This is an explicit
operator obstruction refining the inherited non-monomial bigon result.

For the natural X=E(B_k), r=k. On the source, no kappa-fixed orbit exists:
the free pair {u,w} is fixed positively by both group elements and cannot
satisfy g x=kappa x. On the target, the swap has k negative singleton pair
coordinates and k(k-1) positive two-cycles. The anti-fixed equation is
consistent and has 2^(k+k(k-1))=2^(k^2) solutions. Division by the free C2
action yields **2^(k^2-1)** complement-fixed target classes.

This last mismatch concerns the natural edge-set target; if r<k, the target
also has a fixed free pair and its complement trace can be zero. Do not
generalize the complement obstruction to all target actions. The degree-one
obstruction does cover all of them.

At k=2 the source polynomial is 1+4z+8z^2+10z^3+7z^4+2z^5.
Its weighted quotient matrix has 32 states. Its spectral moments and complement
trace were checked independently of the signed-cycle polynomial. Positive
controls include a two-arm outward tree and directed cycles plus an isolate,
using their full automorphism groups, not just trivial stars.

## P1. Joint signs in the paired cube

An aligned consecutive pair is (T,W,s), with independent tournaments T,W
on [n] and a root-adjacency mask s in F2^n. The acting group is simultaneous
S_n. This represents T together with the rooted (n+1)-tournament obtained by
adjoining the root to W. No condition W=T is imposed.

A Walsh support is (I,J,K), where I,J are sets of core pairs and K is a
set of vertices. Its sign on its **joint** stabilizer is

    epsilon_I(g) epsilon_J(g) = epsilon_(I symmetric-difference J)(g).

Hence the class count P_n is the number of joint support orbits for which
this product is trivial. Separate positivity of both factors is sufficient
but not necessary. For n=2, taking I=J to be the single edge and K empty
gives a stabilizer transposition negative on each factor and positive on
their product. This is the least cancellation witness.

This is the same character mechanism as S1 in a two-edge-layer cube.
The concurrent Session90
[decorated Fourier theorem](../session90_even_graphs_odd_descent_20261029/even/DECORATED_EVEN_GRAPH_DUALITY.md)
uses a different basis: retain R=T XOR W and s as unsigned decorations,
and Fourier-transform only T. Our full Fourier transform in both T,W
becomes chi_(I XOR J)(T) chi_J(R). A further ordinary Fourier transform in
R relates the two bases. It is a linear change of basis, not an asserted
bijection of their decoration objects.

Our counts for demanding separate positivity on the joint stabilizer are
6,72,2452 for n=2,3,4. Session90's different invalid rule, demanding global
sign-evenness of its X without the joint decoration, gives 6,60,1092.
These filters are not the same; neither replaces the correct count
8,88,2752. The common surviving information is the actual joint stabilizer.

## P2. A free Klein four action on odd self-converse classes

For n>=2, define a(T,W,s)=(T^op,W,s) and
b(T,W,s)=(T,W,s+1). They commute and descend to aligned-pair classes.

**Claim.** The group <a,b> is free on these classes.

**Proof.** If g carries a pair to any nonidentity a^i b^j image, it fixes W,
so g belongs to the automorphism group of W. That group has odd order:
otherwise Cauchy's theorem gives an involution, which cannot fix a tournament.
Let m be the odd order of g. If i=1, iterating the equality m times gives
T=T^op, impossible because n>=2. If j=1, it gives s=s+1, also impossible.
One of i,j is nonzero. QED.

The inherited [Eplett / Session89 class correspondence](../session89_consecutive_tournament_pairs_20261028/encoding/ODD_SELF_CONVERSE_LIFT.md)
identifies these classes with self-converse tournament classes of order 2n+1.
Transporting a and b gives a free Klein action there, and proves that their
number is divisible by four for every odd order at least five.
The n=1 exception has two classes, so the range is sharp for this action.

This is a constructive symmetry explanation, not a new scalar divisibility
record: the existing [LOW_BITS.md](../session89_consecutive_tournament_pairs_20261028/counting/LOW_BITS.md)
has stronger all-order 2-adic information. An explicit self-converse
implementation needs a supplied anti-involution or an algorithm to find one.
Class independence follows from the inherited conjugacy/uniqueness theorem.
Named ports still require the centralizer and coherent transports; the
concurrent [nonflat transport examples](../session90_even_graphs_odd_descent_20261029/carriers/TRANSPORT_AND_HOLONOMY.md)
must not be erased by this class-level statement.

In the full Walsh basis, a acts by (-1)^|I| and b by (-1)^|K|.
Freeness therefore also implies that the four parity sectors have the same
dimension P_n/4. Complete computations give 2,22,688 in each sector for
n=2,3,4.

## P3. All three flips and an odd/even core boundary

Adjoin c(T,W,s)=(T,W^op,s), obtaining C2^3. Every nonidentity element
except ac and abc leaves at least one tournament unchanged; the preceding
odd-order argument makes its fixed set empty. The two exceptional elements
cannot both fix a class because their product b is free. Thus stabilizers
have size at most two and all orbits have size four or eight.

If n is odd, abc has no fixed class. Indeed an isomorphism g reversing T
and W has even order 2m with m odd, since g^2 is a tournament automorphism.
Then g^m is an anti-involution on odd n, so it has exactly one fixed vertex
(two fixed vertices would leave their arc unreversed). Its unique fixed
vertex is fixed by g, which commutes with g^m. A mask cannot be complemented
at that vertex.

If n is even, every vertex cycle of such g has even length. Otherwise the
odd power g^m would have fixed vertices, necessarily at least two at even
order, again impossible. For each qualifying g, both equations on s,
g s=s and g s=s+1, then have exactly 2^(number of vertex cycles) solutions.
The twisted Burnside trace formula gives |Fix(ac)|=|Fix(abc)|.

With flip bits ordered T,W,mask, the actual fixed-class vectors are

    n=2: (8,0,0,4,0,0,0,4)
    n=3: (88,0,0,32,0,0,0,0)
    n=4: (2752,0,0,128,0,0,0,128).

If U=|Fix(ac)| and V=|Fix(abc)|, the numbers of four- and eight-element
orbits are (U+V)/4 and (P_n-U-V)/8. At n=3 these are eight and seven,
giving fifteen total. This identifies exactly what the 7 and 15 count;
it supplies no addition law on those orbits and no Hamming-code identification.

## What this gives a closed-carrier search

S1 supplies cheap necessary tests before constructing an expensive coordinate
codec: compare f_1, then F(-1), then higher coefficients. S2 shows why equal
Burnside counts cannot certify an operation-preserving reduction. P1 explains
why separate port or factor tests discard recoverable information; a joint
stabilizer can cancel their signs. P2-P3 distinguish freeness on labelled data,
on classes, and on a chosen marking.

These tests do not by themselves minimize a deterministic attachment state
machine. The current
[minimal singleton/complement/fresh-re-marking carrier, M1-M3](../minimal_blue_selectors_20261028/THEORY.md)
has a weaker, explicitly delimited observer than arbitrary named-path tests.
The incoming minimum is two endpoint-rooted types plus a blue interface with
4n-8 states for n>=5, exactly32272 total behavioral cells at n=7. Its
seven-vertex collision and the 81/27 equal-triples hostile remain
mandatory controls. No formula above reconstructs an individual tournament
from F, nor turns a Fourier basis vector into a tournament.
