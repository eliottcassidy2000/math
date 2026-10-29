# Joint vertex/run clocks: the exact dyadic minimum and a ternary obstruction

2026-10-29. **PROVED scoped completion/extension results; FINITE-EXACT native
modules; CITED faithfully flat descent. Same-agent alternative-path audit.**

Let V=coker(M) for the actual210-by210 Paley7 three-run endpoint matrix.
Retain its named port classes and the group Gamma=Aut(Paley7) x S3 of
vertex automorphisms and run-slot permutations. A completion here is a
Gamma-equivariant injection into a finite abelian group with a perfect
alternating Q/Z-valued pairing and a Gamma action by isometries. A pairing
on V is not prescribed. Ordinary-quotient injectivity is NOT required.

The exact minimum dyadic index is **2^40**. At3 every such completion has
index at least **3^8**, but attainment of that ternary bound is OPEN.
Canonical hyperbolic doubling gives the ternary upper bound3^18.

Inheritance: [PRIMARY_CLOCK_COMPLETION P1--P5](PRIMARY_CLOCK_COMPLETION.md)
gives the C7 minima and the purity boundary. The81/27 pair remains the
hostile source. The corrected near miss is assuming that a separately
extendible clock has a jointly compatible extension. The niche ancestor is
[THM-2779, Bockstein/symplectic decoder and Heisenberg root-degree gate](../../../01-canon/theorems/THM-2779-bockstein-symplectic-decoder-frame-torsor-and-heisenberg-root-degree-gate.md).
That theorem distinguishes retained central phase from a permutation
shadow; it is not a theorem about the present tournament cokernel.

Concept board: named ports / cyclic factors / Frobenius / primary depth /
root attachments / joint clock relations. Anchor: full-clock completion.
Niche: an exact index-two extension recurrence. Wildcard: the obstruction
to lifting two commuting clocks simultaneously. The META-PATTERNS card
used is controlled forgetting with its retained action/kernel sidecar.

## N1. The full group has an actual C21/Frobenius presentation

Write T:v->v+1, F:v->2v, R for the run three-cycle and S for the run swap.
All four preserve the native M with its actual port labels. Define

    C=T R,       H=F S.

On the210 ordered triples, C has ten free orbits of length21 and

    C^21=H^6=1,       H C H^-1=C^2,
    T=C^15, R=C^7, F=H^4, S=H^3.                     (N1)

Thus Gamma=C21 semidirect C6, with the normalizer acting by squaring.
This is a presentation of the supplied action, not an analogy inferred
from equal group orders. Literal port permutations verify all relations.

At2 let eta^2+eta+2=0 with eta even, and eta_bar=-1-eta. Besides X-1 and
Phi3=X^2+X+1, the C21 factor chart has reciprocal pairs

    f+ = X^3-eta X^2+eta_bar X-1,
    f- = X^3-eta_bar X^2+eta X-1,
    g+ = X^6+eta X^5-X^4-X^2+eta_bar X+1,
    g- = X^6+eta_bar X^5-X^4-X^2+eta X+1.            (N2)

The f product is Phi7 and the g product is Phi21. For example g+ is the
norm, over a primitive cube root w, of
X^3-1-eta*w*X^2+eta_bar*w^2*X. Expanding the norm and using
eta*eta_bar=2 proves the sextic formula. Their product with the first
two factors is X^21-1. The Frobenius degrees3 and6 prove irreducibility
modulo2. All factors and their inverses are checked modulo64.

The actual O_f-module partitions are:

| Factor | Degree | Partition | Binary logarithm of source order |
| --- | ---: | --- | ---: |
| X-1 | 1 | (2,2,1,1,1,1) | 8 |
| Phi3 | 2 | (1,1,1) | 6 |
| f+ | 3 | (3) | 9 |
| f- | 3 | (4,1) | 15 |
| g+ | 6 | (4,1) | 30 |
| g- | 6 | (1,1) | 12 |

The sum is80, recovering the entire previously certified dyadic group.
The chart reduces thirty7-orbit blocks to ten21-orbit blocks. It keeps
their inverse maps and checks every native block; it does not discard
port alignment. The C21 minimum is index2^24 by P1: contributions0,0,6,18
from the four self-dual/paired sectors. The extra normalizer matters.

## N2. The full dyadic minimum is index2^40

The four orthogonal sectors have exact minimum logarithmic indices

    4, 0, 18, 18.                                   (N3)

They are orthogonal because invariant pairings connect only reciprocal
characters of the odd normal subgroup C21. The following arguments give
both their lower bounds and explicit attaining completions.

### The X-1 sector: a determinant sign at precision4

The remaining F action splits this sector into its fixed part A and a
nontrivial elementary part B. Exact source bases give

    A=(Z/4)^2,       S_A=[[3,2],[2,1]] mod4,
    B=(Z/2)^4.

The receipt supplies a perfect F,S-invariant alternating form on B.
The determinant of S_A is3 modulo4. An alternating form on A is cJ;
S invariance forces2c=0, so A has no invariant perfect form.

An ambient paired group of order64 cannot repair A. The full list of
paired abelian groups of that order is (Z/8)^2, (Z/4)^2+(Z/2)^2 and
(Z/2)^6. The last has no order-four elements. In the middle case any
embedded (Z/4)^2 has independent primitive order-four directions, so
restriction of a perfect ambient form to it is perfect, a contradiction.
In the first case the source must be the entire order-four subgroup2W.
An ambient isometry has determinant1 modulo8, so its action on2W has
determinant1 modulo4, again contradicting S_A. This list is exhaustive.

The source order is16 and paired ambient orders are squares. Therefore
the next possible order is256, attained by A+A* with its literal
contragredient clock. Its index16 is optimal. Adding the already paired
B proves the first contribution4 in N3. This is an orientation loss
invisible modulo2, not merely an odd source dimension.

### The Phi3 sector: no repair

The native six-dimensional binary module has an explicitly stored
perfect alternating form invariant under all four original clocks.
Thus its contribution is0. The result comes from the full simultaneous
linear system, not from the C21 self-duality test alone.

### The cubic pair: two lawful root attachments are necessary and sufficient

Descent under the order-three F gives the following actual source groups
and residual involutions:

    A=Z/8,        S_A=3;
    B=Z/2+Z/16,  S_B=[[1,1],[8,5]].                 (N4)

The receipt constructs fixed generators, their complete unramified bases
and inverse charts. Thus these are the marked source modules, not just
abstract groups with equal sizes. To pair the two reciprocal cubic
sectors, it suffices to place B inside a larger involution group D and
give a surjective equivariant character D -> A*. Since3^-1=3 modulo8,
the target character also has involution3.

For a character f:B->Z/8 put f(a)=4s and f(b)=t. Equivariance requires
4s+2t=0 modulo8, so every such character has even image. There is no
surjective character on B itself.

No index-two overgroup suffices either. Any such overgroup is obtained
by adjoining x with2x=v. Because the quotient has order2, its involution
has the form Sx=x+w. Necessarily

    2w=(S-1)v,       (S+1)w=0.                     (N5)

Modulo2B the only allowed root classes are v=0 and v=a: S fixes just
that line. For either class, f(v) is divisible by4. An extended character
with value z on x must satisfy2z=f(v), hence z is even as well.
This excludes EVERY index-two overgroup, independent of the root gauge.
The exact recurrence below separately finds all four marked presentations.

An index-four overgroup does suffice:

    D=Z/4+Z/32,   S_D=[[3,1],[8,5]],
    B -> D: (a,b) -> (2a,2b),
    phi(x,y)=6x+y mod8.                            (N6)

The map is injective, S_D^2=1, and phi is surjective with phi(S_D z)=3phi(z).
The two root steps first adjoin a half of a, then a half of b; the first
step changes which second root classes satisfy the clock constraint.
All old generators and their clock values remain recoverable.

Tensoring D with the unramified cubic ring, use the reciprocal trace
pairing between that module and its dual. The source B embeds on one
side and A embeds on the other by the dual of phi. This constructs the
full-clock completion. Its source order is2^24; the ambient order is
|D|^(2*3)=2^42, giving index2^18.

For necessity, every proposed paired ambient cubic sector has dual halves
of the same length. Galois descent under F turns its B-half into an
overgroup of the displayed B. The absent surjection at indices1 and2
rules out all smaller lengths. Descent here is the standard effective
module descent for a finite unramified Galois extension: its semilinear
data are equivalent to base modules. This is a CITED classical input,
checked in [Stacks, Proposition35.3.9](https://stacks.math.columbia.edu/tag/023F).
The application additionally retains the explicit residual involution.

### The sextic pair: the C21 minimum already extends

Here H itself is full Frobenius on the unramified degree-six ring. The
actual fixed bases descend the two sources to Z/16+Z/2 and (Z/2)^2,
with no residual clock. On the larger side retain the first source;
on its dual embed the second by(x,y)->(8x,y). Reciprocal trace pairing
is Frobenius invariant, so this extends all of Gamma. Its index2^18
is already forced by the C21 lower bound. The explicit fixed bases
also verify descent directly on these two supplied finite modules;
no averaging by6 is used in characteristic2.

Adding N3 proves the exact full dyadic index2^40 and ambient order2^120.
The result concerns arbitrary injective completions. Its root embeddings
are not pure; the minimum with injectivity at every ordinary quotient
level remains a separate problem.

## N3. An all-size recurrence for involution-compatible root attachments

For ANY finite abelian group A with involution S, every source-fixed
index-two extension with an extending involution has a presentation

    A'=(A+Zx)/(2x-v),       Sx=x+w,
    2w=(S-1)v,             (S+1)w=0.               (N7)

Choose v from A/2A and solve the two displayed equations for w in A.
Changing x to x+a changes v to v+2a and w to w+(S-1)a. Thus choosing
one root representative loses no source-fixed extension. Residual gauge
can identify different w values; the algorithm counts marked
presentations, NOT unlabelled isomorphism classes.

Each A' element is uniquely(a,t), t=0,1. Its sum has the carry tv:
when both binary digits are1, add v and reset the digit to0. Its clock
is(Sa+t*w,t). These formulas prove sufficiency of N7 and give an exact
recursive construction retaining the entire old group and clock.

For an equivariant character f:A->Z/2^r with target involution u,
the lift value z=f'(x) is characterized exactly by

    2z=f(v),       (u-1)z=f(w).                     (N8)

Iterating N7/N8 therefore preserves the character selector as well as
the root and carry data. In N4, all four first-step presentations fail;
the second level has an explicit success at the35th tested presentation.
The second level is not claimed to be an exhaustive census. The successful
source chain is(v,w)=(a,a+4b), then(b,x+2b), yielding N6.

## N4. At3 the full clocks force an index of at least3^8

The saved source maps now induce every native clock on
V_3=(Z/3)^2+(Z/9)^6+Z/81. The original seven-cycle is trivial here,
but the vertex Frobenius and run clocks are not. Their joint relations
are checked at each coordinate's actual exponent.
The saved right inverses are integer coordinate lifts of quotient
generators, not homomorphic sections into the free port lattice.

Solve for ALL invariant alternating forms b:V_3 x V_3->Q/Z. Store the
upper-triangle variables with denominators3^min(e_i,e_j); impose every
C^T J C=J modulo81. Reversible local elimination produces a complete
generating set, retaining all nonunit pivots. The common radical is

    intersection_b rad(b) = R,
    R is abstractly (Z/3)^3+Z/9+Z/27,  |R|=3^8.     (N9)

An explicit jointly invariant form has radical exactly R. Thus3^8 is
the sharp minimum radical ORDER, not just a bound from reducing modulo3.
The common annihilator is computed by the full image of every form map;
both independent integer Smith and local elimination verify its order
and type. The attaining form has the same two independent checks.
Its image filtration has logarithmic orders10,2,0 under successive
multiplication by3. Thus its nondegenerate quotient V_3/R is
(Z/3)^6+(Z/9)^2. This is the quotient of the displayed form, not a
chosen direct-sum decomposition of V_3.

In any perfect paired ambient W containing V, its restriction b has
rad(b)=V intersect V-perp. Since |V-perp|=|W|/|V|, every completion has
index at least |rad(b)|. Hence N9 gives the universal full-clock bound
3^8, stronger than the old C7-only minimum3^4. Whether3^8 is attainable
with the entire clock action remains OPEN. Hyperbolic doubling gives
3^18. The obstruction below explains why sharpness of a radical bound
is not itself a completion theorem.

All six possible one-dimensional Gamma actions on an additional Z/81
were also tested: F acts by1,28 or55, S by+1 or-1, and R,T trivially.
None of these direct-module additions admits a perfect invariant form.
Their complete mod3 dual-map spaces each have dimension3; all27 members
have rank at most4 in the ten-dimensional ambient group. This finite
exclusion is consistent with N9; it does not rule out larger nonsplit
extensions.

## N5. Individual clock lifts can fail jointly

Let p be odd and V=F_p^3 with basis x,y,z. The commuting source clocks
send x->x+z and y->y+z respectively. In a four-dimensional symplectic
completion with <x,y>=<z,t>=1, ALL possible extensions have

    A(t)=t-y+alpha*z,       B(t)=t+x+beta*z.

Both preserve the pairing and have order p, but

    AB(t)-BA(t)=2z.                                (N10)

Changing alpha,beta cannot remove that central discrepancy. Any invariant
restriction to V has radical containing z; in dimension4 it must have
rank2, so rescaling gives precisely this calculation. Thus no minimal
dimension4 completion retains both commuting clocks. An ambient paired
group of order p^4 must be elementary to contain the rank-three source,
so allowing higher exponents does not help. Hyperbolic doubling attains
dimension6, proving the exact minimum index p^3 rather than the unclocked p.

At p=2 the displayed commutator vanishes, and both clocks do extend to
dimension4. The checker tests ALL p^2 lift pairs at2,3,5,11,13. This
dyadic exception belongs to this specific symplectic parabolic problem;
it is not THM-2779's different D8 refinement. The connection to that
older work is the retained central phase, not an identification of its
source, physical ancestry, or LRC claims.

## Field comparison, audit, and next decisive test

The sextic chart now gives an actual F64 residue-field occurrence:
21 divides63, and Frobenius has period6 on its primitive21 roots.
F256 cannot contain a primitive21 root, since21 does not divide255;
its primitive21 factors instead have degree3, because ord_21(256)=3.
In general the roots split over F_(2^m) iff6 divides m. These are actual
clock-field constraints. They are distinct from the earlier determinant
response experiments on all F64/F256 state triples. Likewise21 here
is a port-clock order, not a rank threshold for elliptic curves.

Reproduction: run `primary_native_clocks.py`, `dyadic_joint_cycle.py`,
`involution_root_extensions.py`, and `clock_kernel_audit.py` in this
directory, with and without `-O`. The last independently enumerates
every kernel for all256 matrices over Z/4 and all6,561 matrices over Z/9
of size2-by2, checking every source vector. New source clocks satisfy
literal intertwiners on all210 named ports, their heterogeneous cyclic
orders, and the full group relations. The81/27 singleton-slice controls
are replayed. Normal/optimized receipts agree; exceptions stay active.

Next: construct or exclude a3-primary completion of index3^8, retaining
the full joint action, rather than only its optimal restriction form.
Its common radical type(1,1,1,2,3), paired quotient(1^6,2^2), and explicit
nine-coordinate matrices give a concrete finite starting point. Test
the joint lift of the commuting vertex/run order-three clocks. The
second question is the pure
full-clock minimum: neither N2 nor the root recurrence promises ordinary
quotient compatibility. Native complement and arithmetic Selmer maps
remain separate operations requiring their own source diagrams. Sha
finiteness, Collatz and LRC(14) remain open.
