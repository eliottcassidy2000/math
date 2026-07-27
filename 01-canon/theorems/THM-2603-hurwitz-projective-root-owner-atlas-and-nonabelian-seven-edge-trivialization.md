---
id: THM-2603
title: "Hurwitz projective root-owner atlas and nonabelian seven-edge trivialization"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE; INDEPENDENT HOSTILE AUDIT REQUESTED;
  not yet promoted into the proved dependency graph.  There is an
  explicit surjective C_2*C_3 representation on P^1(F_13) in which the
  product A of the free-factor generators has order seven and two seven-point
  orbits, while their commutator C has order thirteen and one thirteen-point
  orbit plus one fixed point.  The transition-dependent product of the seven
  owner-conjugates of C, in the stated descending order, is the identity even
  though C^7 is nontrivial.  In the affine root chart, exhaustive ordered
  norms for representatives of all three order-seven conjugacy classes close
  on five nonzero root exponents per orientation; the six closure sets cover
  all of F_13^*, and their sharp selector cover number is three.  Trace and
  orientation are extra connection data, not physical labels supplied here.
  PSL_2(F_13) has no element of order 91 and the
  C_7 and C_13 centralizers are exactly those cyclic groups, so this closure
  is genuinely nonabelian and is not an internal C_91.  The homogeneous
  carrier Omega=PSL_2(F_13)/<C> has 84 points and a canonical six-frame
  projection to P^1(F_13).  Its commuting left C_7 and right C_6 actions give
  84=12*7=14*6, while left C splits it as six fixed boundary frames plus six
  thirteen-cycles over the affine chart.  These are exact permutation
  structures, not a typed identification with an LRC carrier.  The anharmonic S_3
  action splits P^1(F_13) into orbit sizes 3,3,6,2.  The three Barning
  Pythagorean branches are infinite-order unimodular maps whose mod-two image
  is only {identity,swap}; ternary branching is not C_3 torsion, and a
  three-step exact hostile shows that they do not preserve THM-2056's fixed
  Euclidean Farey defect.  This is an abstract projective atlas only: no
  common Boolean LRC carrier, positive semantic edge, target-current
  incidence, row exclusion, or proof of LRC(14) is constructed.
source: psl2z-farey-v4-audit-2026-07-28
depends_on: []
related:
  - THM-2055-determinant-gate-normal-fan-and-tangent-sector-reduction
  - THM-2056-kelvin-polar-farey-defect-certificate
  - THM-2315-marked-target-gain-corolla-and-pairwise-composition-boundary
  - THM-2329-boundary-triple-rerooting-and-transverse-gain-obstruction
  - THM-2542-seven-chart-cech-holonomy-and-c91-arrival-obstruction
  - THM-2584-b-word-depth-five-absolute-deep-root-tensor
  - THM-2585-saturated-normalized-target-projector-and-bockstein-noncommutation
  - THM-2590-boolean-bockstein-and-theta-selector-incidence-spectrum
  - THM-2591-theta-zero-selector-cech-coboundary-and-c91-holonomy-no-go
  - THM-2593-charged-target-section-atlas-and-minimal-c91-holonomy-trivialization
  - THM-2601-linear-bockstein-sheet-coordinate-and-nonlinear-target-monodromy
script: 04-computation/lrc14_hurwitz_projective_root_owner_atlas_thm2603.py
output: 05-knowledge/results/lrc14_hurwitz_projective_root_owner_atlas_thm2603.out
script_sha256: ca6ea3b94811845395eaeb74c1352c13da41153d35c23dfa63ed70137faa8817
output_sha256: 6006c4f3161de03a5471eaedcc247ad14f1118be8b2039b778d3b5f7d4175d47
hash_basis: working-tree bytes (LF)
---

# THM-2603 -- one projective object carries the root and owner clocks

**RESERVED / PROVISIONAL PROOF CANDIDATE; INDEPENDENT HOSTILE AUDIT
REQUESTED.**  Nothing in this file is available as a proved dependency until
the exact companion and proof text pass independent hostile audit and this
banner is explicitly promoted.

The post-THM-2593 frontier has two coprime cyclic coordinates:

```text
root deck:       C_13;
owner clock:     C_7.
```

The proved `C_91` construction treats them as commuting coordinates on an
external product.  The present theorem gives a different exact carrier.  In
the natural fourteen-point projective line, a root cycle and an owner cycle
coexist but cannot commute.  Owner transport rotates the root subgroup itself,
and the seven rotated root steps close by a nonabelian relation.

Passing from the projective line to its six-frame bundle produces a second,
84-point carrier on which the two recurrent tower grammars coexist exactly:
twelve owner orbits of length seven, and six root orbits of length thirteen
plus six boundary frames.  This sharpens the cardinality coincidence without
supplying the missing physical labels.

Within this repository, this is the first recorded exact group-level model of
the transition-dependent mixed square demanded by THM-2591; no literature
priority claim is made.  It is not yet a physical realization of that square.

## 1. A Hurwitz quotient of the two modular free factors

Work in `PSL_2(F_13)`.  Matrices below have determinant one; equality in the
projective group identifies `M` with `-M`.  Reduce the two integer matrices

```text
       [ 0 -1]                [ 3 -7]
x  =  [ 1  0],        y  =   [ 1 -2]                         (1)
```

modulo thirteen.  Direct multiplication over the integers gives

```text
x^2=-I,                    y^3=-I.                           (2)
```

Thus the abstract presentation

```text
PSL_2(Z) = C_2 * C_3                                      (3)
```

sends its order-two and order-three generators to the projective classes of
`x` and `y`.  Define

```text
A=xy,                 C=[x,y]=xyx^(-1)y^(-1).               (4)
```

Their reductions are

```text
       [12 2]                 [5 9]
A  =  [ 3 6],         C  =   [9 6].                         (5)
```

Also

```text
          [12 3]
A^(-1)C = [ 2 6].                                          (6)
```

Exact powers give

```text
ord_PSL(A)=7,          ord_PSL(C)=13,
ord_PSL(A^(-1)C)=7.                                      (7)
```

Enumeration of all determinant-one matrices modulo `+-I` gives `1092`
projective elements.  Breadth-first multiplication by `x,y` reaches all
`1092`; multiplication by `A,C` independently reaches all `1092`.  Hence

```text
<x,y>=<A,C>=PSL_2(F_13),          |PSL_2(F_13)|=1092.       (8)
```

So this is a surjective finite representation of the modular free product,
not a numerical analogy between branching arities.

## 2. The fourteen-point projective carrier

Use the affine coordinate

```text
P^1(F_13)=F_13 union {infinity},
[z:1] <-> z,                 [1:0] <-> infinity.            (9)
```

The two exact cycle decompositions are

```text
A:
 (0 9 12 1 3 6 2)(4 10 7 8 5 11 infinity),                (10)

C:
 (0 8 infinity 2 10 4 1 7 11 12 3 9 6)(5).               (11)
```

Thus:

```text
A = two owner-sized seven-cycles;
C = one root-sized thirteen-cycle plus one projective cusp. (12)
```

This is the exact sense in which one object carries both recurring LRC
alphabets.  It is stronger than the equality `14=2*7`: the two seven-rails
and the thirteen-root deck are permutation structures on the same labelled
projective line.

The root action can be put in the familiar affine form.  Let

```text
       [1 -1]
P  =  [3 -2].                                               (13)
```

Then, projectively modulo thirteen,

```text
P C P^(-1) = [1 1;0 1],          z -> z+1,                 (14)

P A P^(-1) = [7 5;10 11].                                  (15)
```

The owner cycles in this root chart are

```text
(0 4 6 10 7 5 3)(1 8 infinity 2 9 12 11).                 (16)
```

In particular the owner action sends

```text
8 -> infinity -> 2.                                        (17)
```

Therefore the affine thirteen-point root chart is not owner-invariant.  Any
attempt to put the root deck and owner action on one projective carrier needs
the fourteenth point.

This is a precise possible explanation for an existing asymmetry.  THM-2585
has thirteen target-shift sections, while THM-2315's exact target-gain carrier
has fourteen projective directions.  After granting an affine comparison of
the kind explicitly granted in THM-2593, projective completion adds one
missing boundary axis.  Equation (17) proves that an owner action of type
(15) cannot be tested on the thirteen-section atlas alone.  No theorem here
constructs that missing section or identifies the two existing alphabets.

## 3. The nonabelian seven-edge norm

For `i in Z/7Z`, transport the root step into the `i`th owner chart:

```text
C_i=A^i C A^(-i).                                          (18)
```

The multiplication order is load-bearing.  The exact relation is

```text
C_6 C_5 C_4 C_3 C_2 C_1 C_0 = 1             in PSL_2(F_13). (19)
```

There are two proofs.

First, direct multiplication gives the identity.  Second, put

```text
D=A^(-1)C.
```

Then `A^7=D^7=1` projectively by (7), and expansion gives

```text
C_6 C_5 ... C_0
 =A^7(A^(-1)C)^7
 =A^7D^7
 =1.                                                       (20)
```

At the `SL_2` matrix level both seventh powers in (20) equal `-I`, so their
product is literally `I`.

This is not the fixed-chart relation `C^7=1`: since `C` has order thirteen,

```text
C^7 !=1.                                                   (21)
```

Nor may (19) be reversed.  The ascending product

```text
C_0 C_1 ... C_6                                           (22)
```

has projective order thirteen.  The exact multiplication order records which
root chart is incoming and which is outgoing.

Each `C_i` is parabolic with one fixed projective cusp.  In the order
`i=0,...,6`, those cusps are

```text
5, 11, infinity, 4, 10, 7, 8,                              (23)
```

one complete `A`-orbit.  This cusp motion visualizes the transition-dependent
law, but it is **not** the reason (19) closes.  The load-bearing mechanism is
the separate order-seven identity `ord(A^(-1)C)=7` in (7), which turns the
ordered product into `A^7(A^(-1)C)^7`.  Cusp transport alone does not imply a
product relation.  Equation (19) is an edge law, not a change of seven vertex
origins.

### Comparison with the proved Cech no-go

THM-2542 and THM-2591 use one fixed additive coefficient group `F_13` with
trivial owner action.  Seven identical translations then have additive
holonomy `7a!=0`; a vertex selector contributes only a telescoping
coboundary.  Equation (19) changes that coefficient system:

```text
fixed root subgroup C_13
   -> owner-transported conjugates A^i C_13 A^(-i).         (24)
```

Therefore (19) neither contradicts nor weakens the Cech theorem.  It supplies
the exact abstract shape of a genuinely mixed correction: the transition
must know the owner chart and rotate the root subgroup.  A physical use would
need to realize (24) on one common positive ancestry carrier.

### The complete trace-and-orientation norm atlas

The single relation (19) sits inside a larger exact family.  Work in the
affine root chart (14), put

```text
       [1 1]                 [ 0 1]
U  =  [0 1],        g_t  =  [-1 t],       t in {3,5,6},    (24a)
```

and, for `a in F_13^*`, define

```text
H_(t,k)(a)=g_t^k U^a g_t^(-k),
N_t^F(a)=H_(t,0)(a) H_(t,1)(a) ... H_(t,6)(a),
N_t^R(a)=H_(t,6)(a) H_(t,5)(a) ... H_(t,0)(a).            (24b)
```

Both products telescope symbolically:

```text
N_t^F(a)=(U^a g_t)^7 g_t^(-7),
N_t^R(a)=g_t^7(g_t^(-1)U^a)^7.                            (24b1)
```

Moreover,

```text
tr(U^a g_t)=t-a,              tr(g_t^(-1)U^a)=t+a.        (24b2)
```

The noncentral determinant-one matrices of projective order seven have
precisely the six traces

```text
{+-3,+-5,+-6}.                                             (24b3)
```

The two moving matrices in (24b2) are visibly noncentral.  Since `g_t^7`
is central, (24b1) proves the exact criterion

```text
N_t^F(a)=1 iff t-a in {+-3,+-5,+-6},
N_t^R(a)=1 iff t+a in {+-3,+-5,+-6}.                       (24b4)
```

Thus the closure atlas is a trace-intersection theorem.  The exhaustive
companion independently recomputes every product and then supplies the finer
conjugacy classification below.

Each `g_t` has projective order seven.  Projective trace is defined only up
to sign, so the invariant is trace square.  The three values are

```text
t=3: tr^2=9,       t=5: tr^2=12,       t=6: tr^2=10.      (24c)
```

Their conjugacy classes have size `156`, are pairwise disjoint, and exhaust
all `468` elements of order seven in `G`.  Thus (24c) is one representative
of each order-seven projective conjugacy class, not three arbitrary traces.

Exhaustion of all twelve nonzero exponents gives the identity loci

| trace | forward closure `N_t^F(a)=1` | reverse closure `N_t^R(a)=1` | union |
|---|---|---|---|
| `3` | `{6,8,9,10,11}` | `{2,3,4,5,7}` | `{2,3,4,5,6,7,8,9,10,11}` |
| `5` | `{2,8,10,11,12}` | `{1,2,3,5,11}` | `{1,2,3,5,8,10,11,12}` |
| `6` | `{1,3,9,11,12}` | `{1,2,4,10,12}` | `{1,2,3,4,9,10,11,12}` |

The orientation dependence is structural rather than a table accident:

```text
N_t^R(a)=N_t^F(-a)^(-1).                                  (24d)
```

Hence the reverse closure set is the negative of the forward closure set.
For each fixed `(t,orientation)`, the twelve products have the same exact
conjugacy-signature census:

```text
identity:                         5;
order 13, conjugate to U^2:       2;   tr^2=4;
order 3:                          2;   tr^2=1;
order 6:                          2;   tr^2=3;
order 2:                          1;   tr^2=0.             (24e)
```

There are two order-thirteen conjugacy classes, represented by `U` and
`U^2`; each has size `84`.  Every nonidentity order-thirteen output in (24e)
lies in the nonsquare class represented by `U^2`.  The order-two, order-three,
and order-six outputs lie in their respective unique conjugacy classes.
Thus trace changes *which* exponents close, but not the per-state conjugacy
census.

The original relation is present with its exact pair data.  For

```text
Q=[4 5;0 10],       A'=P A P^(-1),
```

direct multiplication gives

```text
Q A' Q^(-1)=g_5,             Q U Q^(-1)=U^3.              (24f)
```

Accordingly, (19) becomes the reverse trace-five closure `N_5^R(3)=1`.
This is simultaneous conjugacy of the pair `(A',U)`, including the root
rescaling.  Conjugating only the order-seven element while silently keeping
the root chart fixed would discard load-bearing relative-position data.

Across all six trace/orientation states, the closure loci cover every
`a in F_13^*`.  The cover has a sharp selector invoice.  The two trace-three
sets are disjoint and partition `{2,...,11}`, missing only `{1,12}`; either
trace-six orientation supplies both missing exponents.  Exhaustion of the
`2^6` state subsets gives

```text
minimum universal state-cover size = 3,
minimal covers = {3F,3R,6F} and {3F,3R,6R}.                (24g)
```

No two connection states suffice, and trace five is redundant for universal
exponent coverage.  This is stronger than the bare union statement, but it
is also a sharper physical warning.  To use (24g), an LRC construction must
produce a typed selector choosing trace class and orientation from the
row/root transition.  Neither the Boolean packet nor the root sheet currently
supplies those two coordinates.  Choosing a closing state after inspecting
`a` would add an unlicensed sidecar, not prove a fixed connection.  The exact
atlas therefore gives a three-state connection target and obstruction invoice;
it does not yet trivialize the physical holonomy.

## 4. Why this is not a hidden C91

The complete projective element-order census is

```text
order 1:       1 element;
order 2:      91 elements;
order 3:     182 elements;
order 6:     182 elements;
order 7:     468 elements;
order 13:    168 elements.                                 (25)
```

There is no element of order `91`.  Exact centralizer enumeration further
gives

```text
Cent(A)=<A>,              |Cent(A)|=7,
Cent(C)=<C>,              |Cent(C)|=13.                     (26)
```

No nonidentity element of `<A>` commutes with a nonidentity element of
`<C>`.  Since `C_7 x C_13` is cyclic of order `91`, (25)--(26) imply:

```text
PSL_2(F_13) contains no internal C_7 x C_13 carrier.        (27)
```

The `C_91` mapping torus of THM-2542/2593 is a correct external abelian
cover.  It is not a subgroup of this projective atlas.  The new closure (19)
uses noncommutation rather than hiding it.

## 5. The homogeneous 84-frame carrier

The fourteen-point action has a canonical frame refinement relative to the
chosen root subgroup.  Put

```text
G=PSL_2(F_13),       H=<C>,       B=N_G(H).                (27a)
```

The element `C` has the unique projective fixed point `5`.  Exact enumeration
gives

```text
B=N_G(H)=Stab_G(5),        |B|=78.                         (27b)
```

More explicitly, the matrix

```text
       [0 3]
R  =  [4 4]                                                   (27c)
```

has projective order six, fixes `5`, and satisfies

```text
R C R^(-1)=C^4,
{4^j:0<=j<6}={1,3,4,9,10,12} in F_13^*.                  (27d)
```

Those are precisely the six nonzero squares modulo thirteen.  Every exponent
in (27d) occurs for thirteen elements of `B`, and

```text
B=H semidirect <R> = C_13 semidirect C_6,
B/H=<R H> congruent C_6.                                  (27e)
```

Now take the left-coset carrier

```text
Omega=G/H,                    |Omega|=1092/13=84.          (27f)
```

There is a `G`-equivariant projection

```text
pi:Omega -> G/B congruent P^1(F_13),
pi(gH)=g(5).                                                (27g)
```

Each fibre has six points.  Right multiplication by `B/H` is well-defined,
free on every fibre, and commutes with every left `G`-action.  Thus (27g) is
a six-frame bundle, not merely a set of cardinality `84`.

The owner element `A` acts freely on `Omega`: if `A^j gH=gH` for
`1<=j<7`, then the order-seven element `g^(-1)A^j g` would lie in the
order-thirteen group `H`.  Hence its exact orbit decomposition is

```text
left A:          12 cycles of length 7;
right R:         14 cycles of length 6;
left A and right R commute.                                (27h)
```

In particular,

```text
84=12*7=14*6.                                              (27i)
```

The root action is sharper still.  A coset `gH` is fixed by left `C` exactly
when `g` normalizes `H`, so the fixed set is the boundary fibre `B/H`.  Every
other orbit has length thirteen.  Consequently

```text
left C:       6 fixed frames + 6 cycles of length 13,
84=6*(1+13).                                               (27j)
```

The chart matrix `P` from (13) sends the fixed cusp `5` to `infinity` and
conjugates `C` to `q->q+1`.  If

```text
U=<[1 1;0 1]>,
```

then conjugation by `P` identifies `G/H` with the requested carrier `G/U`.
Therefore the six fixed points in (27j) are
exactly the six frames above the projective boundary, while the six
thirteen-cycles are exactly the `78` frames above the affine chart.  The
right `C_6` rotates the six boundary frames in one cycle, rotates the six
affine `C_13` orbits in one cycle, and acts on the twelve `A`-orbits as two
six-cycles.  The latter `C_6`-set is abstractly isomorphic to `F_13^*` under
multiplication by its square subgroup; no particular isomorphism is selected.

This is the precise co-occurrence promised by the two tower grammars:

```text
owner grammar:        12 displacement classes * 7 owner cells;
root/frame grammar:   6 frames * (13 affine roots + 1 boundary root). (27k)
```

Equation (27k) does **not** identify the grammars physically.  THM-2584/2586
has a typed positive bank `(s,ell) in F_13^* x F_7`; (27h) supplies twelve
unlabelled `A`-orbits and choices of origins on them.  It has not mapped `s`
to an orbit, `ell` to the left `A` coordinate, or any cell mass to a coset.
Likewise `14*6` suggests projective target direction times nonzero owner
colour, but identifying `B/H` with `F_7^*` requires a chosen `C_6`
isomorphism and supplies no Boolean incidence.

There is a second exact near-match with THM-2590 and the verified candidate
THM-2601.  After choosing one frame over one affine point and an isomorphism
`B/H congruent F_7^*`, the big cell in (27j) may be indexed by
`(q,kappa) in F_13 x F_7^*`; left `C` translates `q` and right `R` rotates
the frame.  THM-2601's verified provisional calculation, however, finds that
the thirteen physical coefficient sheets `Y_q` admit a separating linear
scalar but that target successor is a degree-eleven permutation in that
scalar, and its candidate proof rules out every linear functional that would
make successor affine.  Therefore the comparison from the
Bockstein carrier to (27g) cannot be inferred from the shared `13*6`
indexing.  A lawful comparison must be nonlinear in the coefficient factor
(for example via THM-2601's inverse sheet polynomial), or must retain the
target/frame data as an explicit sidecar.  No existing equivariant physical
map is claimed.

## 6. The exact anharmonic S3 orbit split

THM-2329 identifies the physical rerooting operation on a zero-sum character
triple with the anharmonic `S_3` action on `P^1(F_13)`.  In the coordinate
(9), take

```text
u(z)=1/z,                    v(z)=-1/(z+1).                 (28)
```

After determinant-one rescaling, `u` and `v` have projective orders two and
three and generate a group of size six.  Its four orbits are

```text
boundary:       {0,-1,infinity}          size 3;
harmonic:       {1,6,11}                 size 3;
generic:        {2,4,5,7,8,10}           size 6;
equianharmonic: {3,9}                    size 2.            (29)
```

The point stabilizers have respective orders

```text
2, 2, 1, 3.                                                (30)
```

The last pair is characterized intrinsically by

```text
z^2+z+1=0.                                                 (31)
```

Consequently the eleven transverse directions outside THM-2329's boundary
orbit split exactly as

```text
11=3+6+2.                                                   (32)
```

This locates the precise role of the primes two and three: they are the only
nontrivial stabilizer orders of projective rerooting.  The size-six object is
the six-element rerooting group, not a six-vertex tournament.  In particular,
the proved physical `S_3` orbit of a trivial-leg triangle remains the
three-point boundary orbit.  Moving from that orbit to all fourteen points
requires new operations; the full Hurwitz action from Sections 1--3 is not
automatically a physical Fourier symmetry.

The coordinate convention in THM-2315 writes mixed gains as `[1:g]`, whereas
(9) writes `[z:1]`; inversion swaps these affine conventions.  Every set in
(29) is inversion-stable, so the orbit and stabilizer statements are
unchanged.

## 7. Correcting the binary/ternary tree frame

The Pythagorean triple tree has three children, but branching arity is not
group torsion.  On the primitive parameter column `(m,n)^T`, its three Barning
branches are

```text
       [ 0 1]           [0 1]           [1 0]
B_1 = [-1 2],    B_2 = [1 2],    B_3 = [2 1].              (33)
```

They induce the fractional maps

```text
r=m/n -> 1/(2-r),       1/(2+r),       r/(2r+1),           (34)
```

with images in the three intervals

```text
(1/2,1),                (1/3,1/2),     (0,1/3).            (35)
```

All three have infinite order.  Indeed, `B_1-I` and `B_3-I` are nonzero
nilpotent matrices, so their positive powers grow linearly.  The eigenvalues
of `B_2` are `1+-sqrt(2)`, so its powers cannot be periodic.

Modulo two, however,

```text
B_1 -> swap,            B_2 -> swap,       B_3 -> identity. (36)
```

Thus their image in `GL_2(F_2)=Aut(V_4)=S_3` has order at most two.  The three
branches are not the order-three free factor and do not produce an `S_3`
permutation of the three nonzero `V_4` matching labels.  The exact statement
is:

```text
ternary triple tree = three-generator unimodular semigroup;
C_3 free factor     = one torsion generator of the modular group.     (37)
```

Likewise, the binary Stern--Brocot tree is generated by two parabolic moves;
its two children are not themselves the elements of a `C_2` action.  The
faithful `C_2*C_3` realization relevant here is (1)--(8).

The `K_4` matching identity for Pythagorean triples remains valid as a
separate carrier: the three matchings are the three nonzero elements of a
`V_4`, and `S_4/V_4=S_3` permutes them.  Equation (36) proves that the Barning
branch semigroup does not, merely by being ternary, realize that full
permutation action.

## 8. What the Barning maps do preserve -- and the Farey hostile

The determinants of (33) are

```text
1,-1,1.                                                    (38)
```

Hence the Barning maps preserve primitivity and the absolute determinant-one
relation between Farey neighbors.  They are legitimate enumeration charts
for primitive rays.

They do **not** preserve THM-2056's fixed-basis Euclidean defect.  Use the
calibration from that theorem with owner covector `w=(1,0)`:

```text
F(d)=||d||^2-91 w.d.                                       (39)
```

The visible primitive orbit

```text
d_0=(1,10),
d_1=B_1 d_0=(10,19),
d_2=B_3 d_1=(10,39)                                       (40)
```

has exact defects

```text
F(d_0)=10,             F(d_1)=-449,        F(d_2)=711.     (41)
```

Thus the same unimodular tree gives

```text
certified -> uncertified -> certified                       (42)
```

for the fixed scalar defect.  This is a hostile to transporting a gate
certificate along the tree, not an LRC row counterexample.  It is exactly the
basis warning in THM-2055/2056: a nonorthogonal `GL_2(Z)` change preserves
primitive/Farey incidence but changes the Euclidean norm, owner polygon,
normal cone, and defect.  A lawful Farey descent must transport those
sidecars at every branch.

## 9. Interface to the live LRC objects

The exact maps and losses are:

| source | target reading | preserved | lost / still required |
|---|---|---|---|
| `P^1(F_13)` | THM-2315 target-gain line | all fourteen projective labels, boundary/transverse split, cross-ratio action | physical word-current incidence |
| `Omega=G/<C>` | THM-2584/2586 84-cell bank | commuting `12*7` left-owner and `14*6` right-frame decompositions | typed `(s,ell)` bijection, positivity, cell masses, role data |
| affine `78`-frame big cell | THM-2590/2601 `(q,kappa)` sheet bank | six `C_13` cycles, one for each chosen frame label | a chosen `C_6 congruent F_7^*`, nonlinear coefficient-to-projective map |
| `C=<[x,y]>` | root-deck model | a thirteen-cycle and affine translation chart | actual THM-2542 root label on one Boolean carrier |
| `A=<xy>` | owner-clock model | two seven-cycles | actual THM-2584 owner-cell/rail identification |
| conjugates `A^i C A^-i` | mixed transition cochain | exact seven-edge closure (19) | positivity, ancestry, word, theta, target role, Bockstein unit |
| ordered norms `(t,orientation,a)` | sharp three-state cover of a six-state connection atlas | every nonzero root exponent closes in some state; sharp selector cover number three | typed trace/orientation selector; one fixed physical connection |
| affine patch of `P^1` | THM-2585/2593 thirteen-section atlas | translation orbit after a chosen affine comparison | the fourteenth projective boundary section and owner covariance |
| anharmonic `S_3` | THM-2329 rerootings | exact `3+3+6+2` gain stratification | no transport out of the proved physical boundary orbit |

The two seven-cycles in (10)/(16) are a plausible projective model for the two
rails of THM-2584, but equality of cycle sizes is not a map.  Abstractly, a
`C_7`-equivariant bijection

```text
Psi:F_13^* x F_7 -> Omega,
Psi(s,ell+1)=A Psi(s,ell),                                 (43)
```

exists automatically and noncanonically: choose a bijection from the twelve
`s` labels to the twelve `A`-orbits, then choose one origin in each orbit.
Consequently existence of `Psi` is not a finite test and carries no physical
content by itself.  The first nontrivial test is whether the existing packet
operations induce `Psi` without those choices, make its transported right
`C_6` action agree with a typed sidecar, and realize it on one common positive
Boolean carrier through an equivariant label `Phi` satisfying

```text
Phi(owner successor)=A Phi,
Phi(root transition in owner chart i)=A^i C A^(-i) Phi,     (44)
```

while retaining the `b` word, `theta=0` arrival/deep/later diagonal, and a
unit target-section coefficient.  If (44) exists, (19) pays THM-2591's
cycle-level invoice by a genuine transition-dependent square.  If no such
label exists, the hostile must print the first failed predicate: support,
positivity, endpoint role, owner clock, root action, or missing projective
section.

A partial-cube relabeling alone cannot pass this test.  THM-2584's four-edge
support path is a tree and hence a partial cube, but every Theta-cut cochain
on a tree is a vertex potential.  Pulled around the owner cycle it contributes
zero holonomy, exactly as THM-2591 already proves for every root selector.

## 10. Scope and exact verification

The exact candidate object is

```text
(P^1(F_13); Hurwitz generators x,y; owner A; root C;
 moving-cusp conjugates C_i; full trace/orientation ordered-norm atlas;
 homogeneous frame bundle G/<C>;
 anharmonic S_3 orbit atlas).                              (45)
```

It proves no physical identification between:

```text
THM-2542 root deck,
THM-2584/2586 owner and arrival rails,
THM-2585 target-shift sections,
THM-2315 target-gain directions,
THM-2334 relation residues.                                 (46)
```

In particular, (19) is not permission to replace the additive Cech class by
a projective action.  The transition-dependent action must be constructed
from the LRC packet itself.  No scalar row is excluded and LRC(14) remains
open.

The dependency-free exact companion enumerates every determinant-one matrix
modulo `+-I`, reconstructs the generated subgroups, cycle actions, conjugate
products, the symbolic norm telescopes and moving-trace criterion, all six
trace/orientation tables and their sharp set-cover invoice, conjugacy-class
census, element-order census, centralizers, normalizer, all `84` cosets and
their three commuting orbit decompositions, and the anharmonic orbits, then
checks the integral Barning hostile.  Reproduce with

```bash
python3 04-computation/lrc14_hurwitz_projective_root_owner_atlas_thm2603.py
python3 -O 04-computation/lrc14_hurwitz_projective_root_owner_atlas_thm2603.py
```

Both executions must byte-match

```text
05-knowledge/results/lrc14_hurwitz_projective_root_owner_atlas_thm2603.out.
```

The script performs `260` exact assertions.  The group, centralizer, and
normalizer checks range over all `1092` projective matrices; the homogeneous
checks construct all `84` cosets of size thirteen and exhaust their left and
right actions; the base action checks range over all fourteen projective
points.  No probabilistic or floating-point step is used.

**QED.**
