---
id: THM-2655
title: "Semiregular-kernel quasi-etale return and quartic V4 Kummer class-group gate"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Let a complex
  affine-space Keller map have Galois closure M/K,
  sheet stabilizer H, and a normal subgroup N of monodromy which acts
  semiregularly on the sheets.  Over the normalization R of the target in
  M^N, the normalization Z in M is a connected N-Galois cover unramified in
  codimension one; hence Z over R_reg is an etale N-torsor.  Thus a regular
  normal kernel, dual to THM-2633's forbidden regular block quotient, must
  return as codimension-at-least-two topology.  For quartic A4/S4 monodromy,
  N=V4 and Q=G/V4 is C3/S3.  The irreducible character plane
  W=Hom(V4,C2) injects Q-equivariantly into H1_et(R_reg,mu2).  Kummer theory
  forces its canonical image either wholly into unit squareclasses or
  injectively into Cl(R)[2].  Therefore a quartic survivor requires a
  Q-standard rank-two Kummer carrier on its full Galois resolvent
  normalization.  The even-sign-change S4 quotient realizes the class-group
  alternative sharply but has nonconstant Jacobian.  No A4/S4 exclusion,
  degree bound, JC(2), general JC, or DC(2) follows.
source: root-long-frontiers-2026-07-28-quartic-v4-kummer-gate
depends_on:
  - THM-2633-derangement-character-obstruction-and-d4-keller-exclusion
related:
  - THM-2455-quartic-swallowtail-scaffold-and-endpoint-corrections
  - THM-2595-modular-v4-affine-lift-dichotomy-and-six-vertex-tournament-no-go
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
  - THM-2643-degree-five-six-keller-stabilizer-and-regular-block-quotient-census
script: 04-computation/jacobian_s4_resolvent_quasietale_hostile.py
output: 05-knowledge/results/jacobian_s4_resolvent_quasietale_hostile.out
script_sha256: f2f72df60e5443cfe35177cbb561c7e473b63006b2f7254d080a14a3be949c8f
output_sha256: 01748e08163b3ff45c971ea7597a2a1fe1d7791fef99ef17bf98baee0893a97e
hash_basis: LF-normalized bytes
---

# THM-2655 -- a regular kernel must return on the quotient normalization

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2633 rules out a regular **quotient** of a Keller sheet action: the
normal closure of a point stabilizer is all of monodromy.  A regular normal
**kernel** behaves dually.  It is invisible to every inertia group which
fixes a sheet, so it cannot ramify divisorially after passage to its quotient
field.  It must reappear as a quasi-etale torsor on the quotient
normalization.  For the two live quartic groups this kernel is `V4` and the
quotient is the full Galois cubic-resolvent group.  The missing invariant is
therefore not another discriminant statistic: it is a rank-two Kummer module
on the regular locus of that normalization.

## 1. The semiregular-kernel return theorem

Let

```text
F:A^n_C -> A^n_C                                           (1)
```

be a polynomial Keller map of finite generic degree.  Write its target,
source, and Galois-closure function fields as

```text
K subset L subset M,
G=Gal(M/K),                 H=Gal(M/L).                    (2)
```

The geometric sheet set is `Omega=G/H`.  Let `N` be a normal subgroup of
`G` which acts semiregularly on `Omega`; equivalently,

```text
N intersect gHg^(-1)=1                 for every g in G.  (3)
```

Put `E=M^N`.  Let `R` be the normalization of the affine target in `E`, and
let `Z` be its normalization in `M` (equivalently, the normalization of the
target in `M`).  Then

```text
Z -> R is connected, finite, N-Galois, and unramified in codimension one.
                                                               (4)
```

In particular, with `U=R_reg`,

```text
Z x_R U -> U
```

is a connected finite-etale `N`-torsor, and therefore

```text
pi_1^et(U) ->> N.                                        (5)
```

### Proof

The field extension `M/E` is Galois with group `N`, so the finite
normalization map is connected and `N`-Galois.  Let a height-one prime of
`R` lie below a prime of `Z`.  Its contraction to the finite target is a
target divisor `D`.  Away from the Jelonek set, the original cover and its
Galois closure are finite etale, so the inertia is trivial.

Suppose instead that `D` is a Jelonek component.  THM-2633 proves that its
whole geometric inertia group `I_D` fixes a finite affine inverse sheet.
After choosing that sheet,

```text
I_D <= gHg^(-1)                                           (6)
```

for some `g`.  The inertia in `M/E` is the standard intersection

```text
I_(M/E)=I_D intersect N.                                 (7)
```

Equations (3), (6), and (7) make it trivial.  Thus (4) holds.  On the
regular locus, purity leaves no codimension-at-least-two branch locus, so
the restriction is etale.  It remains connected because it is a nonempty
open subset of the integral normalization `Z`.  This proves (5).

This is the kernel-side complement to THM-2633.  A regular block quotient is
forbidden downstairs; a semiregular normal kernel is forced upstairs into
the etale topology of a quotient normalization.

## 2. The quartic `V4` specialization

Assume the generic degree is four and

```text
G=A4 or S4.                                               (8)
```

Let `V` be the normal Klein four subgroup.  Every nonidentity element of
`V` is a double transposition and hence a derangement in the natural
four-sheet action.  Thus `V` is regular, condition (3) holds, and Section 1
applies with `N=V`.

Write

```text
Q=G/V = C3  for A4,
Q=G/V = S3  for S4.                                      (9)
```

The normalization `R` is the **full Galois resolvent normalization**: it has
degree three in the `A4` case and degree six in the `S4` case.  In the
`S4` case it is not the non-Galois degree-three root field of one classical
resolvent root.  This distinction is load-bearing.

By (5), the regular locus carries a connected `V`-torsor.  Dualizing its
surjective monodromy gives a canonical injection

```text
i:W=Hom(V,C2) -> H^1_et(R_reg,mu_2).                     (10)
```

The construction is `Q`-equivariant.  The two-dimensional `F2`-module `W`
is irreducible: `C3` cycles its three nonzero vectors, while
`S3=GL_2(F2)` acts as the full automorphism group.

## 3. The Kummer unit/class-group dichotomy

Suppose `n>=2`, let `A=Gamma(R,O_R)`, and put `U=R_reg`.  Normality gives

```text
Gamma(U,O_U)=A,                  Pic(U)=Cl(A).            (11)
```

The Kummer sequence over `C` is

```text
0 -> A^*/A^(*)2 -> H^1_et(U,mu_2) -> Cl(A)[2] -> 0.      (12)
```

Intersect the canonical submodule `i(W)` with the unit term in (12).  The
intersection is a `Q`-submodule of the irreducible plane `i(W)`, so it is
either zero or all of `i(W)`.  Exactly one of the following alternatives is
therefore forced:

```text
i(W) lies in A^*/A^(*)2;                                  (13a)

the projection i(W) -> Cl(A)[2] is injective.             (13b)
```

This concerns the canonical torsor image.  It does not say that the two
ambient modules can never contain other abstract copies of `W`.

Consequently a quartic candidate is excluded if neither its unit
squareclasses nor its two-torsion divisor classes contain the required
standard `Q`-module.  Useful special cases are

```text
A^*=C^*       ==> Cl(A)[2] contains W;
Cl(A)[2]=0    ==> A^*/A^(*)2 contains W.                  (14)
```

Thus a factorial resolvent normalization with fewer than two independent
unit squareclasses is impossible.  More generally one must compute the
`Q`-module, not merely the dimensions: the three nonzero classes have to be
cycled by `C3` or carry the full `S3` action.

## 4. A sharp non-Keller `S4` hostile

Let `Z=A^3_(x,y,z)`.  Let `V` act by the four even sign changes, and let
`H=S3` permute the coordinates.  Then

```text
G=V semidirect S3 is isomorphic to S4.                   (15)
```

Every nontrivial element of `V` fixes a codimension-two coordinate axis, so
the quotient by `V` is quasi-etale.  Its invariant ring is

```text
A=C[a,b,c,d]/(d^2-abc),
a=x^2, b=y^2, c=z^2, d=xyz.                              (16)
```

The quotient `Q=S3` permutes `a,b,c`.  The ring is positively graded, so
`A^*=C^*`.  Nagata localization at `d`, or the displayed toric divisor
lattice, gives

```text
Cl(A)=Z^3/<2e_a,2e_b,2e_c,e_a+e_b+e_c>
     =(Z/2)^2.                                           (17)
```

The three nonzero classes are permuted faithfully by `S3`.  Hence (13b) is
realized with exactly the required standard module.  This proves the
Kummer/class-group gate is sharp as a finite-cover statement.

Here the invariant-ring and class-group statements have short algebraic
proofs.  A monomial fixed by all even sign changes has either three even or
three odd exponents, so the invariant ring is generated by `a,b,c,d` with
the sole relation in (16).  It is normal as a finite-group invariant ring;
equivalently, its displayed singular locus has codimension two.  After
inverting `d`, it becomes the Laurent UFD

```text
C[a^(+-1),b^(+-1),d^(+-1)].                              (17a)
```

Nagata therefore generates the class group by
`P_a=(a,d),P_b=(b,d),P_c=(c,d)`.  The Laurent units give exactly the
relations

```text
div(a)=2P_a, div(b)=2P_b, div(c)=2P_c,
div(d)=P_a+P_b+P_c,                                     (17b)
```

whose Smith form is (17).  Finally `S3` invariants are the symmetric
polynomials in `a,b,c` with `abc=d^2`, hence `C[A,B,d]`.

Both outer quotients are smooth affine spaces.  With

```text
e1=x+y+z,       e2=xy+xz+yz,       e3=xyz,
A=a+b+c,        B=ab+ac+bc,        d=xyz,                (18)
```

the induced quartic map `Z/H -> Z/G` is

```text
(e1,e2,e3) |-> (e1^2-2e2, e2^2-2e1e3, e3).              (19)
```

Its graph quartic and squared-pair resolvent are

```text
T^4-2AT^2-8dT+(A^2-4B),                                 (20)

U^3-4AU^2+16BU-64d^2,                                   (21)
```

with the roots of (21) equal to `4a,4b,4c`.  However,

```text
Jac(F)=4(e1e2-e3),                                       (22)
```

so this is not a Keller map.

For normalization-sensitive discriminants, set

```text
Q(W)=W^3-AW^2+BW-d^2.                                   (23)
```

The exact pullback identity is

```text
Disc(Q)(F)=Disc_H (Jac(F)/4)^2.                          (24)
```

The scaled cubic (21) is `4^3 Q(U/4)`, so its discriminant is
`4096 Disc(Q)`.  Recording this factor prevents the raw squared-pair
discriminant from being confused with the normalized quotient cover.

Finally, the `S3`-stable trace-zero slice `a+b+c=0` of (16) is

```text
d^2=-ab(a+b),                                            (25)
```

the `D4` rational double point.  Its three nonzero two-torsion classes are
permuted by diagram triality.  The label `D4` has therefore returned in a
legitimate but different role: it is the minimal **resolvent singularity**
carrying the `V4` character plane, not a surviving quartic monodromy group.
This three-dimensional quotient and its surface slice are not a planar
quartic polynomial-map realization.

## 5. Exact frontier

The theorem replaces a vague resolvent transfer with a concrete decision
problem:

```text
compute A^*/A^(*)2 and Cl(A)[2] as Q-modules
for the full Galois resolvent normalization.              (26)
```

If both lack the standard plane, the quartic branch is closed.  If one
contains it, the carrier has been found but must still be reconciled with
the graph quartic, Jelonek geometry, and the affine source.  The equality of
quartic and cubic-resolvent discriminants alone supplies none of this
topology.

The theorem does not make the classical resolvent a Keller map, identify a
resolvent sheet with a quartic source sheet, or transfer the special
grade-three cusp anatomy.  It excludes neither `A4` nor `S4` by itself and
proves no Jacobian or Dixmier conjecture.

## 6. Exact companion

Run

```bash
python3 04-computation/jacobian_s4_resolvent_quasietale_hostile.py
python3 -O 04-computation/jacobian_s4_resolvent_quasietale_hostile.py
```

Both executions byte-match

```text
05-knowledge/results/jacobian_s4_resolvent_quasietale_hostile.out.
```

The companion performs `38` hard exact checks: the group/coset action,
semiregular fixed loci, quotient invariant rings, quartic and resolvent
identities, both discriminant normalizations, the nonconstant Jacobian,
normal hypersurface singular locus, class-group Smith form, and irreducible
`S3` action on its three nonzero classes.  The declared hashes use
LF-normalized bytes.

An independent hostile audit rederived the tower-inertia intersection,
semiregular-kernel and purity arguments, connectedness over the regular
locus, the equivariant Kummer injection, and the irreducible-submodule
dichotomy.  It also checked the invariant rings, exact `S4` sheet action,
normality, units, class-group lattice, triality action, quartic/resolvent
formulas, the `D4` slice, and the nonconstant Jacobian.  The audit caught the
`4^6` scaled-resolvent discriminant factor and the canonical-image wording
in (13); both are repaired above.  Normal and optimized executions
byte-match the stored transcript and current hashes.

QED.
