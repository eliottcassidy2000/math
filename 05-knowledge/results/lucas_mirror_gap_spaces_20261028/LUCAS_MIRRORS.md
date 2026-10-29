# Eleven golden points, a biplane, and twelve tournament frames

2026-10-28. **PROVED elementary constructions**, with explicitly finite tables
checked by `check_bridges.py`; **CITED** identification of the classical biplane
group and Jordan's theorem. Same-agent audit; no external-priority claim.
The new work here is the assembled Lucas/metric/frame interface and its exact
failure boundary. Neither PSL(2,11) nor its pentagonal reflection is new.

## 1. Inherited mechanisms and the source object

Let W_n be the cyclic binary words of length n with no adjacent occupied
positions, including the empty word. Its toggle graph is the classical Lucas
cube Lambda_n: two words are adjacent iff they differ in one bit. See
[Mollard, *Maximal hypercubes in Fibonacci and Lucas cubes*](https://arxiv.org/abs/1201.1494).
The word positions remain named; rotation is not quotiented out.

The all-period arithmetic reader is inherited from
[THM-4592, eleven-squares Fibonacci torus](../../../01-canon/theorems/THM-4592-the-eleven-squares-fibonacci-torus-is-the-period-five-collatz-cycle-set.md),
as corrected and proved in
[Session74](../session74_eleven_mechanisms_20261027/session74-clocks/RESULT.md).
It bijects W_n with Z[phi]/(phi^n-1) at odd n. Even n has exactly one
three-word fibre, at zero. The concurrent
[golden-fibre dynamics](../session91_packing_golden_carriers_20261029/golden_fibres/GOLDEN_FIBRE_DYNAMICS.md)
proves the minimal local-toggle repair: distinguish those three words at zero.
We use the stronger inherited theorem, not its old finite-only precursor.

For n=5, W_5 has one empty word, five singletons, and five independent pairs.
There are 15 toggle edges: five from empty to singleton, ten from singleton
to pair. Orient an edge by adding a bit. The ten nonempty vertices form C10;
the empty word joins every other vertex of that cycle. Degrees are
5 once, 3 five times, and 2 five times. This specifies every port and arc.

## 2. A faithful golden unit-distance realization

Put zeta=exp(2 pi i/5), and set

    z(w) = sum_(i=0)^4 w_i zeta^i.

These eleven points are distinct, and their unit-distance graph is exactly
Lambda_5. Toggling a bit changes z by a unit vector. For the converse one may
check the 55 differences in Z[zeta], reducing by
Phi_5=z^4+z^3+z^2+z+1: `check_bridges.py` computes d*conjugate(d) exactly and
finds precisely the 15 toggle pairs. This is an exhaustive algebraic proof of
this finite claim, with no numerical tolerance. Distinctness also follows
because a degree-at-most-four polynomial with coefficients in {-1,0,1}
vanishes at zeta only if it is a multiple of Phi_5; two legal words cannot
differ by the all-one vector.

Singletons lie on radius one. A pair {i,i+2} has

    zeta^i+zeta^(i+2) = phi^(-1) zeta^(i+1),
    phi=(1+sqrt5)/2.

Thus the two nonzero radial layers are pentagons. The graph is bipartite;
its plane drawing can overlap and is not an eleven-square packing. Neither
the packing's contact graph nor its clearance inequalities have been imported.

There is a useful dimensional distinction. The cubical complex consisting
of the five coordinate squares [0,1]^{ {i,i+2} } has link Ind(C5)=C5 at zero.
Its intrinsic square cone has angle 5 pi/2. The displayed planar realization
maps the squares to overlapping rhombi; it does not preserve their right
angles. On the link, radial normalization winds twice around the circle:
the cyclic order 0,2,4,1,3,0 advances by 144 degrees at each step. A singleton
and the corresponding pair have the same direction and different radii.
Discarding radius loses their bit-weight layer even though all eleven points
were distinct before that normalization.

## 3. Arithmetic translation closure is Paley11

Evaluation zeta->4 is a ring homomorphism Z[zeta]->F11. On W_5 let

    ell(w)=sum w_i 4^i mod11,       H={1,3,4,5,9}.

The label table is:

| Bit mask | 0 | 1 | 2 | 4 | 5 | 8 | 9 | 10 | 16 | 18 | 20 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| ell | 0 | 1 | 4 | 5 | 6 | 9 | 10 | 2 | 3 | 7 | 8 |

It is bijective. Singletons read as H, pairs as -H. Each directed toggle has
difference in H. Its closure under all arithmetic translations x->x+t is
exactly the Paley tournament x->y iff y-x is in H: containment is immediate,
and the five empty-to-singleton arcs already supply every difference in H.
The undirected translation closure is K11.

The literal operation changed. Arbitrary arithmetic translation is not a
Euclidean translation of these eleven fixed points. In particular the closure
does not preserve their two-colouring or their unit-distance relation.

Two golden roles must also remain distinct. Under zeta->4 the geometric
phi=1+zeta+zeta^(-1) reads as 8, whereas the golden word reader above uses
phi->4. These are the two conjugate roots of t²-t-1 modulo11. There is no
characteristic-zero ring homomorphism sending phi to zeta. The shared finite
word labels intertwine rotation, not a replacement of all real arithmetic.
This agrees with the concurrent
[two-channel observer calculation](../session91_packing_golden_carriers_20261029/zeckendorf/NORM_CONTROLLED_OBSERVERS.md)
and its distinction between the roots4 and8.

## 4. The mirror is nonlinear on points but preserves a biplane

Reflect positions by i->-i. On the finite labels its exact permutation is

    R=(2 8)(3 4)(5 9)(6 10), fixing 0,1,7,
    R(0)=0,
    R(x)=1/x for x in H,       R(x)=5/x for x in -H,
    R(x)=3x^9-2x^4 on all of F11.

For a singleton the first formula is inversion. A pair is
x=6*4^i; its reflection has label6*4^(-i-2), so their product is5.
Euler's criterion x^5=+1 or -1 gives the polynomial formula.

R preserves all fifteen directed toggle arcs. It preserves25 and reverses30
of the55 Paley arcs, so it is neither a Paley automorphism nor a global
converse map. For example0->1 stays forward, but1->2 maps to1->8 backward.

Introduce eleven blocks B_b=b+H. Every block has five points and every pair
of distinct blocks meets in two points, as the difference table of H shows.
This is the classical 2-(11,5,2) biplane. Translation T:x->x+1, multiplication
M:x->4x, and R preserve its incidence relation. On block indices R acts as

    S=(1 5)(2 6)(3 9)(7 8), fixing 0,4,10.

The essential distinction is R != S. A symmetry can act differently on the
two sides of an incidence relation. The point-block identification b<->B_b
was extra structure, not part of the unmarked design.

The very same point reflection appears in
[Ezra Brown, *The Fabulous (11,5,2) Biplane*, pp.91–92](https://personal.math.vt.edu/brown/doc/1152.pdf).
Brown identifies its classical group with PSL(2,11). Our Lucas derivation
recovers that known reflection from local digit geometry.

## 5. The twelve-frame repair and an exact projective line

A frame is a matching f from points to block indices. It decodes arcs by

    x ->_f y iff y belongs to B_(f(x)).

The identity frame f_infinity(x)=x gives Paley11. An incidence automorphism
with point action g and block action g_B transports frames by

    (g.f)(g(x))=g_B(f(x)).

Let f_0=R.f_infinity. Explicitly, on inputs0,...,10,

    f_0=(0,5,7,4,9,3,10,8,6,1,2),
    f_b(x)=b+f_0(x-b) mod11,       b in F11.

These twelve matchings are distinct. Direct substitution gives their actions:

| Generator | On frame label in P1(F11) |
|---|---|
| T | b->b+1, infinity fixed |
| M | b->4b, infinity fixed |
| R | b->6/b, exchanging zero and infinity |

The finite identities involve132 entries per generator and are checked
exactly. They prove closure under every word in the generators.

This action is faithful. If g fixes every frame, fixing f_infinity implies
g_B=g. Fixing f_b then says g commutes with f_b. But f_b has unique fixed
point b, so g fixes every b and is identity.

The projective transformations generate PSL(2,11). Their determinants are
squares (R has determinant-6=5). Powers of T give every upper unipotent;
conjugation by R gives every lower unipotent. These generate SL2 over a
field, as row elimination shows. Quotienting its scalar centre gives order
11*(11²-1)/2=660. Hence the generated point group is PSL(2,11), and its
identity-frame stabilizer has order55 and is exactly <T,M>. An independent
plain-Python breadth-first group enumeration obtains the same660 and55.

This is an explicit equivalence between two actions: a nonlinear exceptional
eleven-point action becomes fractional-linear on the twelve matching frames.
The twelfth frame is not an extra square in the physical packing.

### Exactly one bit survives the Lucas observations

Among these twelve frames, precisely f_infinity and f_0 agree with all15
directed Lucas arcs. Each of the other ten agrees with exactly six.
This is an exhaustive finite statement with the entire frame universe in
`receipt.json`; it makes no claim about all possible tournament completions.

The single additional test1->2 selects f_infinity: it is true there and
false in f_0. This is minimal because two candidates survived beforehand.
Any of the30 differently oriented pairs works. The probe is not a native
unit contact: its squared geometric length is4-phi, not1. To realize this
observation by a unit-distance forcing gadget requires an additional theorem.

## 6. The next Lucas prime is a hostile, with a short exact certificate

For n=7, the29 words biject with F29 by sum w_i24^i;24²=24+1 and24^7=1.
Let R_7 be word reflection transported through this reader, and T:x->x+1.
R_7 fixes five points and swaps twelve pairs, so both generators are even.
Their group is transitive of prime degree29, hence primitive.

Set g_0=id and successively g_j=R_7 T^(a_j) g_(j-1), where

    (a_1,...,a_12)=(19,17,22,23,18,1,18,21,6,1,15,17).

The explicit29-entry permutation has cycle lengths3,1,25; its25th power is
(0 20 23). Jordan's theorem therefore gives <T,R_7>=A29. The needed classical
form is stated as Theorem1.1 in
[Jones, *Primitive permutation groups containing a cycle*](https://arxiv.org/pdf/1209.5169).
The checker verifies the small word and3-cycle directly, independently of
the exploratory SymPy order computation.

Thus the period-five compact biplane closure does not extend unchanged even
to the next odd Lucas prime. The strongest survivor is the faithful finite
word reader and its guarded moves; unrestricted translations plus mirrors
can enlarge the symmetry group dramatically.
