# JC boundary relations are a mod-two presentation, not a Gram kernel

> **STATUS: PROOF-QUALITY CANDIDATE + VERIFIED-EXACT, NOT CANON.**  The
> presentation theorem and its linear-algebra consequences are proved below.
> The exact companion verifies every finite `F4` identity and the four stated
> lattice controls.  No theorem ID is reserved or promoted here.  Application
> to a quartic Keller resolvent remains conditional on producing the stated
> rational equivariant SNC completion; reflection completion and geometric
> realization remain open.

The open phrase “boundary relations” after
[THM-2708](../01-canon/theorems/THM-2708-c3-hermitian-gain-holonomy-discriminant-gate.md)
contains two genuine Kummer sources and one Gram false positive.  They become
separate only after replacing the boundary intersection matrix by the actual
boundary-class presentation.

The clean answer is:

```text
boundary relations       K/2K              = unit squareclasses;
boundary nonsaturation   (L_sat/L)[2]      = Picard two-torsion;
Gram-only excess         rad(im A)         = no Kummer class at all.
```

Thus the exact object is the two-term complex

```text
C = Z^{boundary components}  --delta-->  P = Pic(Xbar),       (1)
```

reduced modulo two.  Its first mod-two homology is the whole Kummer module.
The Hermitian matrix of THM-2708 is the Gram shadow of this complex, not its
presentation matrix.

## 1. Inheritance and live concept board

The closest proved mechanism is the equivariant unit/class-group Kummer row
in
[THM-2655](../01-canon/theorems/THM-2655-quartic-keller-resolvent-v4-quasietale-torsor-and-kummer-class-group-gate.md)
and
[THM-2685](../01-canon/theorems/THM-2685-equivariant-kummer-boundary-parity-completion-and-divisor-residue-gate.md).
The closest exact matrix is THM-2708's Hermitian `F4` block.  The canonical
hostile is the `D4` discriminant plane: it is a genuine abstract `C3` standard
plane but its order-two discriminant pairing is symplectic.  The corrected
near miss is “`det B=0` means a Kummer plane”; it means only that the boundary
Gram shadow has a standard zero mode.  The least-used sidecar is the map from
boundary components to the ambient Picard lattice itself.

The active board is:

1. the integral presentation `C -> P`;
2. the relation lattice `K=ker(delta)`;
3. the saturation defect `L_sat/L`;
4. the `F4` class matrix `A`;
5. the Hermitian Gram shadow `B=A^dagger G A`;
6. the mod-four discriminant pairing and a stable metabolizer.

The first four give the exact answer.  The fifth is a fast exclusion shadow.
The sixth is a stronger necessary test when only the Gram lattice is known.

## 2. Geometric setup

Let `Xbar` be a smooth projective rational complex surface, let `sigma` act
on `Xbar` with `sigma^3=1`, and let

```text
U = Xbar \ D                                                  (2)
```

where `D` is a reduced `C3`-stable simple-normal-crossings divisor.  Write
`I` for the irreducible components of `D` and define

```text
C = direct_sum_(i in I) Z[D_i],
P = Pic(Xbar),
delta:C -> P,               [D_i] |-> its Picard class.       (3)
```

Both groups are free abelian, `delta` is `C3`-equivariant, and the
intersection form on `P` is integral and unimodular.  Put

```text
K = ker(delta),
L = im(delta),
L_sat = P intersect (L tensor Q).                            (4)
```

No independence, full-rank, nondegeneracy, or tree hypothesis is imposed in
Sections 3--7.

The same proof works whenever `Pic(Xbar)` is torsion-free and the relevant
intersection pairing is perfect.  Rationality is retained because it is the
completion box used by THM-2703/2708 and makes every cohomological input
literal.

## 3. The exact boundary-presentation theorem

There is a canonical `C3`-equivariant identification

```text
H^1_et(U,mu_2)  isomorphic_to
ker(delta_2:C/2C -> P/2P).                                 (5)
```

Moreover there is a natural exact sequence

```text
0 -> K/2K -> ker(delta_2) -> (L_sat/L)[2] -> 0.            (6)
```

Under `(5)`, the left term is exactly the global-unit squareclass branch and
the right term is exactly the Picard two-torsion branch.

### 3.1 Cohomological proof of `(5)`

The Kummer row on `Xbar` gives

```text
0 -> C^*/C^{*2} -> H^1_et(Xbar,mu_2) -> Pic(Xbar)[2] -> 0.
                                                                    (7)
```

Every complex constant is a square and `P` is torsion-free, so the outside
terms vanish and

```text
H^1_et(Xbar,mu_2)=0.                                      (8)
```

The localization sequence for `U subset Xbar`, together with purity in
codimension one, begins

```text
0 -> H^1_et(U,mu_2)
  -> direct_sum_(i in I) F_2[D_i]
  -> H^2_et(Xbar,mu_2).                                   (9)
```

The last arrow sends a component to its cycle class.  The Kummer injection

```text
P/2P -> H^2_et(Xbar,mu_2)                                 (10)
```

identifies this arrow with `delta mod 2`; no Brauer-class ambiguity enters
its kernel.  Intersections among components have codimension two and enter
later degrees of the support spectral sequence, not the degree-one residue
row `(9)`.  Thus exactness of `(9)` gives `(5)`.  Every map is functorial for
`sigma`, so the identification is equivariant.

Concretely, `(5)` sends a Kummer character to its parity residues along the
deleted boundary components.  A parity vector occurs exactly when its
component class is zero modulo two in `P`.

### 3.2 Arrow-by-arrow proof of `(6)`

For `x mod 2C` in `ker(delta_2)`, choose a lift `x in C`.  There is a unique
`p in P` such that

```text
delta(x)=2p,                                               (11)
```

because `P` is torsion-free.  Define

```text
theta(x mod 2C)=p+L in P/L.                               (12)
```

This is well defined.  Replacing `x` by `x+2y` replaces `p` by
`p+delta(y)`, which gives the same class modulo `L`.  Equation `(11)` makes
`2(p+L)=0`, so `theta` lands in `(P/L)[2]`.

It is surjective: if `p+L` has order at most two, then `2p=delta(x)` for
some `x in C`, and that `x mod 2C` lies in `ker(delta_2)` and maps to `p+L`.

Its kernel is the image of `K/2K`.  If `theta(x)=0`, write
`p=delta(y)`.  Then

```text
delta(x-2y)=0,                                             (13)
```

so `x mod 2C` is the reduction of an element of `K`.  Conversely every
element of `K` maps to zero.  Finally the map `K/2K -> C/2C` is injective:
if `k in K` equals `2c` in `C`, then `2delta(c)=0`; torsion-freeness of `P`
gives `c in K`.

The torsion subgroup of `P/L` is exactly `L_sat/L`, since `P/L_sat` is
torsion-free.  Therefore

```text
(P/L)[2]=(L_sat/L)[2],                                    (14)
```

and `(12)`--`(14)` prove `(6)` at every arrow.

### 3.3 The Tor correction is not optional

Regard `(1)` as a chain complex in degrees one and zero.  Its integral
homology is

```text
H_1=K,                    H_0=P/L.                        (15)
```

The universal-coefficient sequence for reduction modulo two is exactly

```text
0 -> H_1 tensor F_2
  -> H_1([C -> P] tensor F_2)
  -> Tor_1^Z(H_0,F_2) -> 0,                              (16)
```

which becomes `(6)` because

```text
H_1([C -> P] tensor F_2)=ker(delta_2),
Tor_1^Z(P/L,F_2)=(P/L)[2].                               (17)
```

Thus a computation which reduces only an integral basis of `K` modulo two
misses the saturation/Picard branch.  A computation which uses only the
intersection Gram matrix can add a different false-positive branch.  Both
errors are visible in `(16)`.

### 3.4 Identification with units and Picard classes

The divisor localization row is

```text
C^* -> O(U)^* -> C --delta--> P -> Pic(U) -> 0.           (18)
```

Hence

```text
O(U)^*/C^* = K,                 Pic(U)=P/L.               (19)
```

Constants are squares, so

```text
O(U)^*/O(U)^{*2}=K/2K.                                   (20)
```

Substituting `(19)`--`(20)` into the Kummer sequence on `U` recovers `(6)`.
This proves that the two terms are not merely isomorphic vector spaces: they
are the actual unit and class-group branches demanded by THM-2655.

## 4. The exact `F4` gate

Let

```text
W=the nontrivial irreducible two-dimensional F_2[C3]-module,
F_4=End_(F_2[C3])(W).                                     (21)
```

Because `2` does not divide `3`, every mod-two `C3` module is semisimple.
If the boundary contains `s` free three-component orbits, then

```text
(C/2C)_std isomorphic_to F_4^s.                           (22)
```

Fixed boundary components contribute only trivial modules.  Choose an
`F4` basis of the standard isotypic part of `P/2P`, say

```text
(P/2P)_std isomorphic_to F_4^t.                           (23)
```

The standard restriction of `delta_2` is an `F4`-linear matrix

```text
A:F_4^s -> F_4^t.                                        (24)
```

Taking standard parts in `(5)` gives the exact identity

```text
H^1_et(U,mu_2)_std
  isomorphic_to underlying_F2 ker(A),                     (25)

multiplicity_W H^1_et(U,mu_2)=nullity_F4(A).             (26)
```

Combining `(6)` and semisimplicity also gives

```text
nullity_F4(A)
 = multiplicity_W(K/2K)
 + multiplicity_W((L_sat/L)[2]).                          (27)
```

This is the complete boundary-relations classification:

```text
nullity(A)>0 through K/2K          = unit-standard plane;
nullity(A)>0 through L_sat/L       = Picard-standard plane;
nullity(A)=0                       = no quartic Kummer carrier.       (28)
```

The splitting of `(6)` as a `C3` module is noncanonical, but the two
multiplicities and their sum are canonical.  For a full `S3` conclusion one
must additionally retain the reflection action; a `C3` standard summand may
be exchanged with another summand by the reflection.

## 5. Why the Hermitian Gram block can overcount

Let `Q` be the unimodular intersection pairing on `P`.  Modulo two it is
perfect and `C3`-invariant.  The trivial and standard isotypic parts of
`P/2P` are orthogonal: a cross-pairing would be a `C3` map from `W` to the
trivial dual and no such nonzero map exists.  Perfectness of the full form
therefore makes its restriction to the standard isotypic part perfect.

Multiplication by `sigma` on `(P/2P)_std` has adjoint multiplication by
`sigma^{-1}`.  Under the `F4` coordinate `(23)`, this is the conjugation

```text
bar(omega)=omega^2.                                       (29)
```

Consequently the standard restriction of `Q` is represented by a
**nonsingular Hermitian** matrix

```text
G=G^dagger in Mat_t(F_4).                                 (30)
```

The hypotheses needed for nonsingularity are precisely:

1. `P` has a perfect integral pairing preserved by `C3`;
2. reduction modulo two remains perfect (automatic for a unimodular `P`);
3. `2` does not divide `|C3|`, so the isotypic splitting is orthogonal.

If any of these is dropped, a singular `G` must be retained as another
radical sidecar and the next identity changes.

Let

```text
M=delta^* Q delta                                         (31)
```

be the boundary intersection matrix.  Functoriality of the standard block
gives

```text
B=A^dagger G A.                                           (32)

```

This is the same Hermitian Gram-block construction used by THM-2708.  Here
equation `(32)` is also formed when relations make the full integral boundary
matrix singular; that scope extension is proved in this reflection and is
not attributed retroactively to THM-2708.  Put `I=im(A)` and

```text
rad(I)=I intersect I^{perp_G}.                            (33)
```

Then `A` induces an exact sequence

```text
0 -> ker(A) -> ker(B) --A--> rad(I) -> 0.                 (34)
```

Indeed, `Bx=0` says

```text
<Ax,Ay>_G=0 for every y,                                  (35)
```

which is exactly `Ax in rad(I)`.  Every vector in `rad(I)` has a preimage
under `A`, and the kernel of the displayed map is `ker(A)`.  Therefore

```text
nullity_F4(B)
 = nullity_F4(A) + dim_F4 rad(im A).                      (36)
```

The consequences are sharp:

```text
det(B)=1          => no standard Kummer plane,
rad(im A)=0       => ker(B)=ker(A),
det(B)=0          => relation/saturation OR Gram-radical false positive. (37)
```

In particular, the negative direction of THM-2708 extends through arbitrary
boundary relations: `det(B)=1` still excludes the standard Kummer carrier.
It is the positive interpretation of a singular `B`, not the exclusion,
which needs the presentation sidecar `A`.

## 6. Schur complements retain the same debt

Suppose the columns of `A` are partitioned as

```text
A=[A_0 A_1]                                               (38)
```

and the pivot Gram block

```text
B_00=A_0^dagger G A_0                                    (39)
```

is nonsingular.  The `G`-orthogonal projection away from `im(A_0)` is

```text
Pi=1-A_0 B_00^{-1} A_0^dagger G,
A_1^perp=Pi A_1.                                         (40)
```

The Schur complement is

```text
S=B_11-B_10 B_00^{-1}B_01
 =(A_1^perp)^dagger G A_1^perp.                          (41)
```

Because `B_00` is nonsingular, `A_0` is injective and its image is
nondegenerate.  An actual relation

```text
A_0x_0+A_1x_1=0                                          (42)
```

is equivalent to

```text
A_1^perp x_1=0,
x_0=-B_00^{-1}B_01x_1.                                   (43)
```

But `(36)` applied to `A_1^perp` gives

```text
nullity(S)=nullity(A_1^perp)
            +dim rad(im A_1^perp).                       (44)
```

Thus Schur complementation is a lawful elimination operation, but it does
not turn a Gram zero mode into a relation.  Its lost coordinate is still the
radical of the projected image.

## 7. Exact controls

The companion

```text
04-computation/jc_boundary_relation_presentation_f4_scout_20260728.py
```

with stored transcript

```text
05-knowledge/results/jc_boundary_relation_presentation_f4_scout_20260728.out
```

and LF-normalized hashes

```text
script_sha256=b220a3acbbdadc1819d6d4114dad8f07d01a84f379483c1bac1d45245514f310
output_sha256=92e2fcdfea21b62c4b015c18081d65ef240e31f7bf1d5e13dca5791a9068d3d3
```

separates four cases.

### 7.1 Genuine unit relation

Take two free boundary orbits and map both to the same free orbit in `P`.
On standard parts,

```text
A=[1 1],             G=[1],
B=[1 1;1 1].                                            (45)
```

Then

```text
nullity(A)=1,         nullity(B)=1,       rad(im A)=0.   (46)
```

Integrally, `delta=[I_3 I_3]`.  Its relation lattice has rank three, its
cokernel is torsion-free, and the mod-two kernel has dimension three.  The
standard part is a genuine unit plane generated by the difference of the two
boundary orbits.

### 7.2 Genuine saturation class

Take `delta=2I_3`.  Then `K=0`, while

```text
L_sat/L=(Z/2)^3.                                         (47)
```

On standard parts `A=[0]`.  Its one-dimensional `F4` kernel is the Picard
standard plane; the accompanying trivial line is the fixed part of `(47)`.

### 7.3 Minimal isotropic-radical false positive

Let

```text
P=I_(1,6)=<1> direct_sum <-1>^6                         (48)
```

with `C3` fixing the positive vector and cycling each of two negative
triples.  Map one free boundary orbit by

```text
d_i |-> e_i+f_i.                                         (49)
```

The map is primitive and injective over both `Z` and `F2`, so

```text
K=0,                  L_sat/L=0,
ker(delta_2)_std=0                                      (50)
```

so any geometric realization would have zero standard Kummer module.  Yet
its Gram matrix is

```text
M=-2I_3.                                                 (51)
```

Equivalently,

```text
A=[1;1],             G=diag(1,1),       B=[0].           (52)
```

Thus `nullity(B)=1` while `nullity(A)=0`; the entire Gram kernel is
`rad(im A)`.  This is the smallest `F4` false positive.  The lattice is
compatible with a rational-surface Picard signature, but no claim is made
that the three displayed classes are realized by an SNC divisor in a
quartic resolvent completion.

### 7.4 Nonsingular exclusion

For `A=G=[1]`, one has `B=[1]` and all three nullities in `(36)` are zero.
This is the minimal negative control.

The script exhausts every `A` of size at most `2 by 2` against every
nonsingular Hermitian target form of size one or two over `F4`, checking
`(36)` case by case.  It also checks the Smith/Tor identity `(16)` on all
four integral controls.  Normal and optimized executions are required to
byte-match the stored output.

## 8. The mod-four isotropy and metabolizer refinement

This section is a proved lattice lemma under its displayed hypotheses.  Its
quartic use remains conditional on the geometric completion hypotheses.

Let `M` be any symmetric integral matrix and put

```text
E=ker(M mod 2).                                           (53)
```

For `x,y in E`, define

```text
lambda_M(x,y)=x^T M y / 2 mod 2.                          (54)
```

The numerator is even because `Mx` is even.  Changing a binary lift by
`2z` changes `(54)` by `z^TMy`, which is even because `My` is even.  Hence
`lambda_M` is a well-defined symmetric `C3`-invariant `F2` pairing on `E`.

Every genuine presentation kernel is totally isotropic for `(54)`.  If
`x,y in ker(delta_2)`, write

```text
delta x=2p,                  delta y=2q.                  (55)
```

Then

```text
x^TMy=(delta x).(delta y)=4p.q,
lambda_M(x,y)=0.                                         (56)
```

Thus a standard plane in `ker(B)` on which `lambda_M` is symplectic cannot
come from either boundary relations or saturation.

Now assume additionally that `delta` is injective, `L` has full rank in
`P`, and `M` is its Gram matrix.  The inclusions

```text
L subset P=P^* subset L^*                                (57)
```

embed

```text
H=P/L subset A_L=L^*/L.                                  (58)
```

The subgroup `H` is `C3`-stable and isotropic for the discriminant pairing:
representatives in `P` pair integrally.  Since `P` is unimodular,

```text
|A_L|=|det M|=[P:L]^2=|H|^2.                             (59)
```

Nondegeneracy of the discriminant pairing gives

```text
H=H^perp.                                                 (60)
```

So `H` is a stable metabolizer.  Under the standard Tor identification

```text
ker(M mod 2) -> A_L[2],
x |-> [x/2] in L^*/L.                                    (61)
```

Equivalently, after identifying `L^*/L` with `Z^r/MZ^r`, the same class has
cokernel coordinate `[Mx/2]`.  Keeping these two coordinate conventions
separate prevents an integral cokernel vector from being mistaken for the
dual-lattice representative itself.

The order-two discriminant pairing is

```text
b([x/2],[y/2])=x^TMy/4 mod Z,                             (62)
```

and doubling `(62)` is exactly `(54)`.  This proves the pairing
identification, not merely its matrix resemblance.

Consequently, in the full-rank independent completion box:

```text
Pic(U)[2] contains W
 => the corresponding W in ker(M mod 2) is lambda_M-isotropic;

Pic(U)[2] contains W
 => 4 divides |H|
 => 16 divides |det M|.                                  (63)
```

The `D4` lattice has determinant four and its unique standard plane has

```text
lambda_W=[0 1;1 0].                                      (64)
```

It therefore has neither the order nor the isotropy needed for a
`C3`-stable metabolizer carrying `W`.  This does not contradict the affine
`d^2=abc` carrier in THM-2655: that carrier is outside the full-rank
independent boundary box.  Rather, it identifies exactly which completion
hypothesis must fail there.

For a non-full-rank but nondegenerate `L`, the subgroup

```text
L_sat/L subset A_L                                       (65)
```

is still isotropic, but need not be a metabolizer.  For degenerate `L`, the
discriminant group and `(62)` are unavailable; the presentation theorem
`(5)`--`(28)` remains valid and is the primary object.

The exact artifact

```text
04-computation/c3_discriminant_metabolizer_mod4_thm2711.py
```

supports this layer.  On current `origin/main`,
[THM-2711](../01-canon/theorems/THM-2711-c3-stable-discriminant-metabolizer-and-mod-four-isotropy-gate.md)
is **PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED** for precisely
the full-rank independent metabolizer scope of `(57)`--`(64)`.  It does not
cover dependent, non-full-rank, or degenerate boundaries; Sections 3--7 are
the presentation/Tor extension into that remaining scope.

## 9. Proposed theorem statement and dependency graph

### Proposed statement, without reserving an ID

> Let `U=Xbar\D` be a `C3`-equivariant SNC open surface with `Xbar` smooth
> projective rational.  If `delta:Z^{Irr D}->Pic(Xbar)` is the boundary-class
> map, then `H^1_et(U,mu_2)=ker(delta mod2)` equivariantly, and there is a
> natural exact sequence
> `0->ker(delta)/2ker(delta)->ker(delta mod2)`
> `->(im(delta)_sat/im(delta))[2]->0`.
> On standard parts, if `A` is the induced `F4` presentation matrix, the
> quartic standard-plane multiplicity is `nullity_F4(A)`.  If `B` is the
> Hermitian boundary Gram block, then
> `nullity(B)=nullity(A)+dim rad(im A)`.  Hence `det B=1` excludes the
> standard Kummer carrier without boundary-independence assumptions, while
> a singular `B` is positive only after the presentation or saturation
> sidecar is restored.

The clean dependency graph is

```text
THM-2655: quartic survivor requires standard W in H^1(U,mu_2)
       |
THM-2685: W lies in units or Pic[2]
       |
boundary presentation delta:C->P
       |-- K/2K ---------------- unit branch
       |-- Tor_1(P/L,F2) ------- saturation/Pic branch
       `-- A over F4 ----------- exact W multiplicity
                    |
                    `-- B=A^dagger G A
                           |-- det B=1: exact exclusion
                           `-- singular: add rad(im A) debt
                                      |
                                      `-- mod-four isotropy/metabolizer
                                          is a stronger necessary shadow
```

THM-2695 becomes relevant only after `A` has a genuine kernel: it decides
whether a class-branch character lifts from `mu_2` to `mu_4`.  It cannot
repair a Gram-radical false positive because no `mu_2` character exists
there.

## 10. Connection contract and exact frontier change

```text
source object:
  a C3-equivariant SNC boundary and its component-class presentation;

target object:
  the standard part of H^1_et(U,mu_2);

map:
  reduce delta:Z^D->Pic(Xbar) modulo two and take its F4 standard kernel;

preserved predicate:
  exact multiplicity of the V4 character plane;

destroyed by the Hermitian Gram quotient:
  whether a zero mode maps to zero in Pic(Xbar) or merely to an isotropic
  radical vector;

needed sidecars:
  the integral relation lattice K, the saturation quotient L_sat/L, and the
  C2 reflection action for a full S3 conclusion;

cheapest decisive test:
  compute the Smith form and C3 action of the actual boundary-class map,
  not only the intersection matrix, for each surviving quartic resolvent
  completion.
```

The mathematical frontier is narrower than “boundary relations are open.”
They are completely classified once the class map `delta` is known.  The
remaining geometric problem is to extract that presentation from an actual
quartic resolvent completion.  A Gram computation alone cannot do so, and
equation `(36)` states the exact amount by which it can overcount.
