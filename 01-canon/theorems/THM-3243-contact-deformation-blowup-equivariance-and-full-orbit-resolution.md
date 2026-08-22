---
id: THM-3243
title: "Contact-deformation blowup equivariance and full-orbit resolution"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  An affine automorphism lifts over the blowup at one marked point exactly
  when it fixes that point.  Hence the Singer action on the F13 contact
  deformation plane lifts over the delayed-contact origin, but the full
  Heisenberg action does not.  Blowing up the reduced union of all 169
  rational contact classes restores full AGL2-equivariance.  Its 2,366
  rational exceptional flags form one AGL2 orbit and split under H13 into
  a 169-point vertical affine orbit and a regular 2,197-point nonvertical
  orbit.  The latter is not the 2,197-point THM-3240 address carrier: their
  orbit and stabilizer spectra differ exactly.
source: root/multiscale-newton-flag/low-child-flag-extension/2026-08-02
audit: >
  The assertion-independent exact companion pins THM-3234, THM-3240, and
  THM-3241; checks the Heisenberg point orbit and stabilizer uniformly at
  p=3,5,7,13; generates all 26,208 elements of GL2(F13) from explicit
  matrices; verifies the full-center ideal identities by exact bivariate
  polynomial arithmetic; checks all exceptional direction orbits and
  stabilizers; and replays the order-168 Singer, order-14 projective, and
  order-12 radial factors together with the 169-class contact chart.
  Normal and optimized runs byte-match the stored transcript and
  LF-normalized hashes below.  An independent hostile audit rederived the
  universal-property lifting criterion, full-centre ideal invariance,
  rational exceptional boundary, both Heisenberg orbit stabilizers, the
  AGL flag stabilizer, Singer order split, and the exact-address orbit-spectrum
  obstruction, and found no defect.
depends_on:
  - THM-3234-singer-owner-compactification-and-pointed-heisenberg-carrier-gate
  - THM-3240-exact-address-heisenberg-clutch-on-carrier-imbalance
  - THM-3241-finite-field-contact-Singer-realization-and-order-gate
related:
  - THM-3228-four-jet-heisenberg-minimal-faithful-permutation-carrier-gate
  - THM-3236-contact-spectrum-primitive-element-and-root-reconstruction-gate
script: 04-computation/contact_deformation_blowup_equivariance_thm3243.py
output: 05-knowledge/results/contact_deformation_blowup_equivariance_thm3243.out
script_sha256: 35b1e2de125ea6e90c98c5d9545a3512dc1a4db72e25faa43642dfb59fb2c4c9
output_sha256: 76c345c64f9f0ebf9f77c346e88c30f647a932839635fff616fe7035b18fe955
hash_basis: LF-normalized bytes
---

# THM-3243 -- contact-deformation blowup equivariance and full-orbit resolution

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3241 identifies the order-two contact-deformation slice over `F_13`
with an affine plane and makes its zero class the delayed-contact point.
THM-3234 puts both a Singer cycle and the standard affine Heisenberg action
on an abstract affine plane.  These two structures behave differently under
the most elementary resolution of the delayed-contact point.

The obstruction is not cardinality.  It is stabilization of the blowup
centre.  Resolving the whole centre orbit repairs equivariance, but produces
a flag carrier whose Heisenberg orbit type is provably different from the
equal-cardinality exact-address carrier of THM-3240.

## 1. The one-centre lifting criterion

Let `k` be a field, let

```text
V=A_k^2,
phi(v)=Mv+b,              M in GL_2(k), b in V(k),       (1)
```

and write `pi_a:Bl_a(V)->V` for the blowup at the reduced rational point
`a`.  Functoriality of blowup gives a unique isomorphism

```text
tilde(phi)_(0,b):Bl_0(V) -> Bl_b(V)                     (2)
```

such that

```text
pi_b tilde(phi)_(0,b)=phi pi_0.                         (3)
```

There is an automorphism `Phi` of `Bl_0(V)` satisfying

```text
pi_0 Phi=phi pi_0                                         (4)
```

if and only if

```text
b=0.                                                       (5)
```

Indeed, `(2)` follows because `phi` carries the centre `0` isomorphically
to `b`.  If `b=0`, this already gives `(4)`.  Conversely, a lift in `(4)`
would make the pullback of the ideal of `0` along `phi pi_0` invertible.
When `b!=0`, its zero scheme is the single point

```text
pi_0^(-1)(-M^(-1)b),                                     (6)
```

away from the exceptional divisor.  There `pi_0` is an isomorphism and the
point ideal in a smooth surface has height two, hence is not invertible.
The universal property of the blowup therefore forbids `(4)`.

Thus a translation does not fail to act altogether: it gives the cospan
`Bl_0(V)->Bl_b(V)`.  What fails is a self-lift over the same marked base.
This distinction is load-bearing.

## 2. The minimal invariant centre

Now take `k=F_p`.  The standard Heisenberg action is

```text
rho(x,y,z)(r,w)=(r+x,w+z-yr).                           (7)
```

It sends zero to `(x,z)`.  Hence

```text
Stab_(H_p)(0)={(0,y,0):y in F_p},
H_p.0=V(F_p).                                            (8)
```

In particular the minimal reduced zero-dimensional `H_p`-invariant centre
containing zero is the reduced union of all rational points.  The same is
true for `AGL_2(F_p)`.  Its ideal is

```text
I_Z=(X^p-X,Y^p-Y),
Z=V(I_Z)=coprod_(v in F_p^2){v}.                        (9)
```

This ideal is invariant scheme-theoretically, not merely setwise.  For any
affine map `phi(v)=Mv+b` defined over `F_p`, its `i`th coordinate satisfies

```text
phi_i(X,Y)^p-phi_i(X,Y)
 =sum_j M_(ij)(X_j^p-X_j).                              (10)
```

Thus `phi^*I_Z` is contained in `I_Z`; applying the same identity to
`phi^(-1)` gives equality.  Consequently

```text
Bl_Z(V)                                                  (11)
```

is canonically `AGL_2(F_p)`-equivariant, and hence `H_p`-equivariant.
No choice of a distinguished rational centre remains.

## 3. Exceptional flags and their exact orbit spectrum

Because `Z` is a disjoint union of `p^2` smooth rational points, the
exceptional divisor of `(11)` is a disjoint union of `p^2` projective lines.
Every `F_p`-point of the blowup lies over a point of `Z`, so

```text
Bl_Z(V)(F_p)
 ={(v,L):v in F_p^2, L in P(T_vV)(F_p)},
#Bl_Z(V)(F_p)=p^2(p+1).                                 (12)
```

The full affine group is transitive on these point-direction flags.  Its
flag stabilizer has order

```text
#AGL_2(F_p)/(p^2(p+1))=p(p-1)^2.                        (13)
```

The Heisenberg subgroup sees strictly more structure.  The derivative of
`(7)` is

```text
D rho(x,y,z)=[[1,0],[-y,1]].                            (14)
```

It fixes the vertical direction `L_infinity=span(0,1)`.  On the remaining
directions `L_s=span(1,s)` it acts by

```text
s |-> s-y.                                               (15)
```

It follows that the rational exceptional boundary is exactly the disjoint
union of two `H_p`-orbits:

```text
E_vertical={(v,L_infinity):v in V(F_p)},   size p^2,
E_affine={(v,L_s):v in V(F_p),s in F_p},   size p^3.     (16)
```

The first orbit has stabilizer of order `p` and is the standard affine
`H_p`-set `V`.  The second has trivial stabilizer and is a regular `H_p`
torsor.  Equivalently, after choosing one nonvertical base flag,

```text
Bl_Z(V)(F_p) = V disjoint_union H_p                    (16a)
```

as an `H_p`-set; the first summand is canonical from the invariant vertical
direction, while the torsor identification of the second needs that
base-flag choice.  Adding the Singer element of THM-3234 generates
`AGL_2(F_p)` and
merges these two strata into the single flag orbit `(12)`.

## 4. The Singer lift and the radial normal gauge

Let `sigma` be multiplication by a primitive element of `F_(p^2)`, viewed
as a matrix in `GL_2(F_p)`.  It fixes zero, so Section 1 gives a self-lift
to `Bl_0(V)`.  On the exceptional line

```text
E_0=P(T_0V)=P^1,                                         (17)
```

its projective action has exact order `p+1` and is regular on the `p+1`
rational directions.  Moreover

```text
sigma^(p+1)=lambda I,
ord(lambda)=p-1.                                         (18)
```

The scalar in `(18)` fixes `E_0` pointwise and acts on its normal line by
the order-`p-1` radial gauge.  Thus blowup separates the two factors of the
Singer order

```text
p^2-1=(p+1)(p-1)                                        (19)
```

geometrically: projective direction on the exceptional divisor and radial
normal scaling transverse to it.

## 5. The exact p=13 contact-deformation specialization

Use THM-3241 with

```text
S=x^2-2,
m=2,
D_2={S+S^2H mod S^3:deg H<2}.                           (20)
```

The map

```text
H |-> c_2=(S')H mod S                                   (21)
```

identifies `D_2` with `F_169`, hence with `F_13^2`.  The zero class is the
unique deformation whose contact is delayed beyond order two.  The helper
`H=1+10x` maps to `alpha=1+2x`, and its multiplication matrix is

```text
[[1,4],[2,1]],                                           (22)
```

of order `168`.

Therefore:

1. `Bl_0(D_2)` carries the Singer lift.  Its exceptional line has one
   `14`-cycle of rational directions, and the fourteenth power is the
   radial scalar gauge of order `12`.
2. A nontrivial Heisenberg translation carries `Bl_0(D_2)` to the blowup at
   another deformation class and has no self-lift over `Bl_0(D_2)`.
3. Blowing up all `169` rational deformation classes restores the full
   formal `AGL_2(F_13)` action.  The resulting rational flag carrier has

   ```text
   13^2(13+1)=2366=169+2197                             (23)
   ```

   points.  Under `H_13` the two summands are respectively the vertical
   affine orbit and the regular nonvertical orbit; under `AGL_2(F_13)` they
   form one orbit with stabilizer `13*12^2=1872`.

The `2,197` in `(23)` is **not** the `2,197`-element exact-address quotient
`G_delta` of THM-3240.  The exceptional nonvertical flags are one regular
`H_13` orbit with trivial stabilizer.  THM-3240's `G_delta` is thirteen
orbits of size `169`, each with a stabilizer of order `13`.  Hence no
`H_13`-equivariant bijection between them exists despite equal cardinality.
The vertical `169`-orbit, by contrast, is the standard affine carrier of
THM-3234.

## 6. What the resolution does and does not remember

The blowup construction supplies an exact algebraic-moduli answer to the
centre-drift problem:

```text
one marked centre  ->  only its affine stabilizer lifts;
full rational orbit -> the full affine group lifts.      (24)
```

It does not identify the contact deformation plane with an LRC owner,
endpoint, root, or exact-address carrier.  Such an identification needs a
chosen basis and derivative trivialization even algebraically, and no
lawful physical intertwiner is supplied.  In particular this theorem does
not strengthen THM-3240's finite address clutch into a physical current.

Nor is a point of `E_0` automatically the next contact jet.  It is the
projective tangent direction of a one-parameter approach to the delayed
class.  Reading it as a later coefficient requires a specified Rees/family
map and its compatibility with the contact filtration.

Finally, `2366` is the number of rational points of the full-orbit blowup.
The `170` fixed-head lower bound in THM-3234 belongs to a different finite
permutation problem and is not a blowup cardinality.  The appearances of
`169`, `2197`, and `2366` expose exact orbit strata, but without a physical
map they are structure, not an LRC ledger decrement.

## 7. Exact companion

Run from the repository root:

```bash
python3 04-computation/contact_deformation_blowup_equivariance_thm3243.py
python3 -O 04-computation/contact_deformation_blowup_equivariance_thm3243.py
```

The companion uses exact finite-field and polynomial arithmetic only.  It
pins the three proved dependencies, contains no optimization-sensitive
assertion, verifies the affine generators and full finite orbit spectra,
and checks the contact/Singer specialization independently of the prose.
Ordinary and optimized runs must byte-match the pinned transcript.

QED.
