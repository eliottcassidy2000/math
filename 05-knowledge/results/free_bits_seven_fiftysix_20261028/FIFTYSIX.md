# A port-preserving 56-object correspondence

**PROVED, elementary; exact exhaustive controls.** No literature-priority claim.
The existing [boundary-support packet](../boundary_support_20261020_seven.md)
proves the matching split105=7+42+56. The construction below supplies a further
explicit map for its56-piece. Classical Fano/GL(3,2) geometry is not new.

Let V=F2^3 and G=GL(V), keeping the zero address and every port label.
Consider three G-sets:

1. oriented noncollinear triples of nonzero vectors: ordered bases(a,b,c)
   modulo cyclic rotation, but not reversal;
2. elements T of G of order three;
3. perfect matchings of the eight points V with four distinct XOR directions.

There are natural G-equivariant bijections between them.

## The forward map and its inverse

Send the oriented basis(a,b,c) to the linear map

    T(a)=b, T(b)=c, T(c)=a.

Its fixed line is span(a+b+c), and its invariant complementary plane is
K={0,a+b,b+c,c+a}. Define the matching

    {0,a+b+c}, {a,a+b}, {b,b+c}, {c,c+a}.                   (1)

Its four directions are{a,b,c,a+b+c}, a quadrangle of the Fano plane and a
weight-four word of the Hamming code. Formula(1) is unchanged by cyclic
rotation and changes under reversal. It is equivariant because it uses only
addition and retained labels.

Conversely, an order-three T over F2 has minimal polynomial dividing
(t+1)(t^2+t+1), whose factors are distinct. A nonidentity action in dimension
three must have a one-dimensional fixed space and a two-dimensional irreducible
space K=im(T-I). Let a be its unique nonzero fixed vector. Pair

    x with a+T(x), for x in K.                            (2)

The two cosets K and a+K are disjoint, so this is a matching. Its directions
are a+(I+T)x, all four points of a+K since I+T is invertible on K.

Finally start with a four-direction matching M. Its direction set D has XOR
zero. Four distinct nonzero vectors in F2^3 with XOR zero form an affine plane
avoiding0; put K=D+D. Every edge crosses K and its other coset. Write its
pairing as f:K->V\K and a=f(0). Every permutation of a four-point affine plane
is affine: AGL(2,2) embeds in S4 and both groups have order24. Thus
f(x)=a+Lx with L in GL(K). Distinct directions say I+L is invertible, excluding
identity and the order-two elements of GL(2,2); L has order three. Extend by
T(a)=a. This recovers T uniquely.

For a fixed T, the set V\(K union{a}) has size3 and is one T-orbit.
Its vectors form a basis: their sum is a, so they are not a Fano line, and
three distinct nonzero dependent binary vectors would have sum0. Following T
therefore recovers exactly one oriented basis modulo cyclic rotation. Reversing
it gives T^{-1} and the other matching over the same underlying triangle.
This completes the inverse proof.

Counts: |G|=7*6*4=168, hence bases modulo rotation give56. Alternatively choose
K in7 ways, a outside K in4 ways and one of the two order-three actions on K:
7*4*2=56. The theorem does not identify these objects with56 unlabelled
six-vertex tournament classes; that family has no supplied compatible G-action.

## Arithmetic and Frobenius coordinates

The Hamming columns1,...,7 correspond to primes
7,2,5,19,89,59,29. A quadrangle D has XOR zero, so its four selected primes
have square sigma. Picking which of its four points is a+b+c and choosing one
of two orientations produces the eight matchings over D. Thus

    7 quadrangles * 4 marked points * 2 orientations = 56.

The bit distinguishing orientations is lost by remembering only D and its
marked point. No parity or coprimality rule in this prime bank recovers it.
A later test sensitive to the pairing must keep it or supply a certificate
selecting one orientation.

In F8=F2[alpha]/(alpha^3+alpha+1), Frobenius x->x² fixes1 and has trace-zero
plane{0,2,4,6} in polynomial-basis masks. Formula(2) gives

    {0,1}, {2,5}, {4,7}, {6,3},

with directions{1,3,5,7}. Every order-three element is a GL conjugate of this
Frobenius. Hence the56 matchings are the56 such conjugates with the original
binary port addresses retained. Relabelling the actual primes is not a
sigma-value-preserving arithmetic operation; only the code predicate transports.

## The nearby 56 that does not match equivariantly

There are also56 unordered triples of the eight affine points. Each completes
uniquely to an affine plane by adjoining x+y+z; equivalently these are the
56 pointed planes(14 planes*4 marked points). This is a real bijection and a
useful three-to-one reconstruction certificate.

But under G fixing0, those triples have three orbits:

    21 containing0; 7 Fano lines; 28 noncollinear triples.

Our matching/order-three/oriented-triangle G-set is transitive. Thus no
G-equivariant bijection identifies the pointed-plane56 with this56 under the
same affine-address action. Choosing an unstructured bijection would erase a
real symmetry distinction. Manivel's Proposition8 instead identifies oriented
Fano triangles with triples on P1(F7) under the exceptional PSL(3,2)=PSL(2,7)
action: that is a different action on eight points, not the affine action with
one fixed point. This explains how an apparent obstruction can be repaired by
changing the action explicitly, rather than concealing the change.

## Why this is relevant to boundary forcing

The matching atlas acts on an eight-port relation by the equality kernel
B_M(s)=product_{xy in M}[s_x=s_y]. Formula(1) now indexes the56 kernels by
oriented triangles. The map preserves which actual ports are joined, so it
preserves each contraction sum against a fixed labelled boundary relation.
It does not preserve that relation if the ports are silently permuted, nor
does it turn a cubic three-edge-colour relation into a six-vertex-colour unit-
distance gadget. Geometry and colour semantics must still be authenticated.

The verification enumerates all105 matchings, all168 invertible matrices,
all56 order-three matrices, all56 oriented triangles, and all9,408
object/symmetry pairs. No random or bounded-height inference is needed.
