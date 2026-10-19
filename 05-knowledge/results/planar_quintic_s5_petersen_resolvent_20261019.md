# The quintic source flow: S5, a Petersen resolvent, and an A5 reduction

2026-10-19. **PROVED.** The quintic field of the nonzero-time source flow
of `H=x⁶t²(1+x²t)` has geometric and arithmetic Galois group S5. Each
source coordinate generates that quintic field and is not expressible
by radicals over the source rational field. An explicit degree-ten
pair-sum resolvent realizes the Petersen graph, reconstructs K5, and
retains the full Galois action. After a quadratic discriminant extension,
the group is A5, abstractly the icosahedral rotation group.

No historical novelty claim is made. The graph and Platonic statements
describe algebraic permutation structures; they imply no spatial symmetry
of the source flow and no graph-planarity criterion for planar JC.

## 1. The source quintic and geometric irreducibility

Let K have characteristic zero and fix `λ∈K*`. Use the source parameters
from [packet 127](planar_monomial_carrier_flow_trichotomy_20261019.md):

    z=x³t,   E=x⁶t²(1+x²t),   F=K(x,t)=K(z,E).

The actual source image `Z=exp(λ{−,H})(z)` satisfies

    Q(V)=AV⁵+5z⁵V²−3Ez⁵,
    A=3E−5z²−15λz⁵,                                 (1)

and `Q(Z)=0`. Its coordinate images are

    X=Z³/(E−Z²),   T=(E−Z²)³/Z⁸,   Z=X³T.            (2)

The initial z,E are algebraically independent. We now enlarge only the
constant field `K(z)` to an algebraic closure C and keep E transcendental.

Regard (1) as linear in E:

    Q(V)=3(V⁵−z⁵)E+
          [5z⁵V²−5z²V⁵−15λz⁵V⁵].                  (3)

The two coefficients in `C[V]` are coprime. At a root `V=ζz` of the
first, where `ζ⁵=1`, the second is

    5z⁷(ζ²−1−3λz³)≠0.

Indeed ζ is algebraic over K, while z is transcendental over K and
`λ≠0`. Gauss's lemma first applies to the primitive degree-one polynomial
in E and then to the polynomial in V. Its V² coefficient `5z⁵` is a
unit in `C[E]`. Thus Q is irreducible over `C(E)` and, consequently,
over F. Its geometric Galois group is transitive on five roots.

## 2. A formal local transposition

Set `κ=1+3λz³∈C*`. At `E=0`,

    Q(V)=5z²V²(z³−κV³).                              (4)

The three nonzero roots in (4) are simple and therefore lift uniquely
to roots in `C[[E]]`. This formal implicit-function assertion follows
coefficientwise: each next coefficient is determined linearly by the
nonzero derivative at the initial root.

For the double root at zero, substitute `E=s²,V=sW` and divide by s²:

    (3s²−5z²κ)s³W⁵+5z⁵W²−3z⁵=0.                   (5)

At `s=0`, the two simple solutions are `W=±√(3/5)`. They lift uniquely
to `W_±∈C[[s]]`, yielding two roots `V_±=sW_±` with opposite nonzero
linear coefficients. The automorphism `s↦−s` fixes the three roots in
`C[[s²]]` and exchanges these two roots by uniqueness of lifting.

All five roots thus lie in `C((s))`, with `E=s²`. Their generated field
over `C(E)` is a splitting field, stable under `s↦−s`. Its restriction
is a single transposition in the geometric Galois group. This is a
formal local calculation; no analytic monodromy theorem is invoked.

## 3. The group is S5 over both constant fields

Let G be the transitive geometric subgroup of S5 just obtained. Make a
graph on the five root labels whose edges are the transpositions
conjugate in G to the transposition from section 2. This graph is
G-invariant and has an edge. Transitivity makes all its connected
components have the same size, which divides five. Thus it is connected.
The transpositions along the edges of a connected graph generate the
full symmetric group. They all belong to G, so `G=S5`.

The geometric splitting-field group injects by restriction into the
splitting-field group over `F=K(z,E)`: an automorphism fixing C and E
permutes the same five roots and preserves their F-generated field.
The latter group is itself a subgroup of S5. Therefore

    Gal(Q/C(E))=S5,   Gal(Q/F)=S5.                    (6)

In particular the result remains valid for arbitrary characteristic-zero
K and every fixed nonzero λ, rather than only for a generic time.

## 4. Each coordinate generates the quintic and is not radical

Packet 127 proves that both X and T in (2) are nonrational over F and
that `F(X,T)=F(Z)` has degree five. By the tower law and primality of
five, each proper intermediate possibility has degree one. Hence

    F(X)=F(T)=F(Z),   [F(X):F]=[F(T):F]=5.             (7)

The normal closure of each field in (7) is the splitting field of Q,
whose group is S5. Neither X nor T is expressible by radicals over F.
Here the standard radical-solvability criterion is applied to the normal
closure: if an element belongs to a finite radical tower, adjoining
the necessary roots of unity and conjugates embeds its normal closure
in a solvable Galois extension. S5 is not solvable. See
[Romyar Sharifi, *Abstract Algebra*, §10.3, Theorem 10.3.10 and Example 10.3.12](https://math.ucla.edu/~sharifi/notes/algebra-ch10.html)
for the criterion and nonsolvable quintic setting.

## 5. An explicit irreducible degree-ten resolvent

Normalize Q to

    f(V)=V⁵+pV²+q,
    p=5z⁵/A,   q=−3Ez⁵/A.

Let its five distinct roots be `v_1,…,v_5`. The ten pair sums
`s_ij=v_i+v_j`, `i<j`, are distinct. Equality of sums sharing an index
would identify two roots. If two disjoint pairs had equal sums, apply
the S5 transposition exchanging one of their indices with the unused
fifth index. Subtracting the resulting equalities again identifies
two distinct roots, a contradiction.

Their monic resolvent is

    R(Y)=∏_(i<j)(Y−v_i−v_j)
        =Y¹⁰+pY⁷−11qY⁵−p²Y⁴+7pqY²−p³Y−q².          (8)

One exact coefficient certificate is

    Res_V(f(V),f(Y−V))
       =(Y⁵+8pY²+32q) R(Y)².                        (9)

Indeed the resultant is the product over all ordered pairs `(i,j)`;
the diagonal factor is `∏_i(Y−2v_i)=32f(Y/2)`, while every off-diagonal
pair occurs twice. Expansion of the Sylvester determinant verifies
the displayed coefficients in (8); both candidate square roots are
distinguished by their leading coefficient, so the monic one is unique.

S5 is transitive on unordered pairs, and the ten sums are distinct.
Thus R is irreducible of degree ten over F. Its splitting field equals
the splitting field of Q, since for distinct i,j,k,

    v_i=(s_ij+s_ik−s_jk)/2.

Accordingly its Galois group is S5 in its faithful action on ten pairs.

## 6. Petersen and K5 are intrinsic to this resolvent

Use the ten resolvent roots as vertices, joining `s_ij` and `s_kl`
exactly when `{i,j}` and `{k,l}` are disjoint. The result is the Kneser
graph `KG(5,2)`, namely the Petersen graph. Its vertices can equivalently
be viewed as the ten edges of K5; adjacency records disjoint edges.
This relation is intrinsic to the resolvent's Galois action: the stabilizer
of one root has nontrivial orbits of sizes three (disjoint pairs) and six
(intersecting pairs), so adjacency is the unique valency-three orbital
relation.

For an explicit standard Petersen labeling, use labels modulo five and
write

    o_i={2i,2i+1},   h_i={2i+2,2i+4}.

The o_i form the outer five-cycle, each `o_i` is adjacent to `h_i`,
and the h_i have inner-star adjacency `h_i~h_(i±2)`.

The graph itself reconstructs the five original labels. An independent
set is a pairwise-intersecting family of two-element subsets. Any such
family of size four has a common point: if `{a,b}` and `{a,c}` occur
and another member omits a, it must be `{b,c}`, and there is no fourth
distinct pair intersecting all three. Therefore the maximum independent
sets are exactly the five four-element stars

    S_i={ {i,j}:j≠i }.

Every graph automorphism permutes these five stars. Their pairwise
intersections recover each graph vertex, `S_i∩S_j={{i,j}}`, so the
induced action is faithful. Conversely every permutation of the five
labels preserves disjointness. Hence

    Aut(Petersen)=S5.

The Galois action from (6) realizes this entire automorphism group.
The five stars and their pairwise intersections also reconstruct the
vertices and edges of K5. There are other maximal independent sets of
size three; the reconstruction uses the maximum sets of size four.

## 7. The quadratic A5 reduction and the Platonic connection

For the monic quintic, the discriminant is

    Δ=q(108p⁵+3125q³).                               (10)

Let `d=∏_(i<j)(v_i−v_j)`, so `d²=Δ`. Under the full S5 group, d
transforms by the sign character. It is not in F, and

    F(d)=N^(A5),   Gal(N/F(√Δ))=A5,

where N is the splitting field. This is a genuine quadratic extension;
as a separate check, (10) has valuation one at `E=0`.

A5 is isomorphic to the orientation-preserving rotation group of a
regular icosahedron, equivalently its dual dodecahedron. A primary
academic exposition explicitly exhibiting that correspondence is
[J. M. S. Simões-Pereira, *The alternating group of degree 5 and the icosahedron rotations: seeing not only proving*](https://arxiv.org/abs/1602.03732).
This supplies the precise Platonic group identification after the
discriminant base change. The original S5 is not being identified with
the full spatial symmetry group of the icosahedron.

## 8. An exact golden-ratio degeneration and the coefficient 11

There is also a precise secondary explanation of the coefficient 11 in
(8). Set `h=1/E` and specialize the normalized coefficients at `h=0`,
holding z and λ fixed. Their denominators are units at h=0, and

    p=5z⁵h/[3−(5z²+15λz⁵)h] → 0,
    q=−3z⁵/[3−(5z²+15λz⁵)h] → −z⁵.

Consequently the actual quintic and its pair resolvent specialize to

    f_∞(V)=V⁵−z⁵,
    R_∞(Y)=Y¹⁰+11z⁵Y⁵−z¹⁰.                       (11)

The quintic limit is squarefree. Its roots are `zζ^i`, with ζ a primitive
fifth root of unity and indices modulo five. Rotation `i↦i+1` partitions
the ten unordered pairs into two five-element orbits, of cyclic distance
one and two. These are the sides and diagonals of the root pentagon. Each
orbit has constant `(Y/z)⁵`, because rotating a pair multiplies its sum
by ζ.

For the precise values, put `r=ζ+ζ⁻¹` and `φ=1+r`, choosing the square
root of five so that `r=(√5−1)/2`. Then `φ²=φ+1` and cyclotomic reduction
gives

    (1+ζ)⁵=−8−5r=−φ⁵,
    (1+ζ²)⁵=−3+5r=φ⁻⁵.

Thus `v=Y⁵/z⁵` satisfies

    v²+11v−1=0,
    v=−φ⁵ or φ⁻⁵,
    R_∞(Y)=(Y⁵+φ⁵z⁵)(Y⁵−φ⁻⁵z⁵).                (12)

Indeed `φ⁵−φ⁻⁵=11`. If `ψ=1−φ=−φ⁻¹`, the numbers
`L_n=φ^n+ψ^n` satisfy `L_0=2,L_1=1,L_(n+2)=L_(n+1)+L_n`, so this 11 is
exactly the fifth Lucas value `L_5`. It occurs in an actual specialization
of the source-flow resolvent, rather than through a numerical analogy.

On the ten pairs, disjointness gives the same Petersen graph as in
section 6. Each of the two rotation orbits induces a five-cycle, and
the five remaining edges match the two cycles. This is the usual
two-pentagon Petersen structure, with one cycle drawn as a pentagram.

The rotation C5 here is a symmetry of the limiting root configuration,
not a claim about its specialized Galois group: `V⁵−z⁵` already has
the base-field root z, and multiplication by ζ does not fix that root.
The generic S5 theorem remains the primary result. A shared coefficient
11 alone supplies no correspondence with unrelated structures such as
Collatz dynamics or the Golay code.

## Verification and limits

The universal group proof uses geometric irreducibility, formal simple-root
lifting, and the connected conjugate-transposition graph. A short exact
companion checks the resultant and discriminant formulas, the finite
Petersen labeling, and the golden-ratio degeneration; those checks do not
substitute for the Galois argument.

These statements concern the generic source function field and its actual
nonzero-time algebraic branch. Numerical specializations can have smaller
Galois groups and are not classified. This result neither changes the
nonpolynomial status of the source flow nor resolves planar JC.
