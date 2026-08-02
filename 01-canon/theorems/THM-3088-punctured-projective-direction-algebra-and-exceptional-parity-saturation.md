---
id: THM-3088
title: "Punctured projective-direction algebra and exceptional parity saturation"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  For every finite field F_q, punctured one-dimensional subspaces are
  orthogonal Boolean idempotents and identify Fun(P1(F_q)) integrally with the
  F_q^*-invariant algebra on F_q^2 minus the origin.  This algebra exhausts
  the parity-even punctured algebra exactly for q=2,3.  The residual
  within-line coordinate for larger q is F_q^*/{+-1}; its exact integral
  contrast ranks and the projector-normalized q-torsion are computed.  This
  is a pointwise finite-field algebra theorem, not a physical tree, quartic,
  Keller, owner, or LRC intertwiner.
source: root-projective-direction-algebra-2026-08-01
depends_on: []
related:
  - THM-2768-modular-c2-c3-quotients-to-a4-s4-and-bass-serre-cycle-ranks
  - THM-3045-k4-edge-isotypic-binary-ternary-integral-clutch
  - THM-3076-finite-affine-plane-line-quotient-tomography-and-p2-three-view-law
  - THM-3083-exceptional-binary-point-ternary-direction-s4-tomography-clutch
script: 04-computation/punctured_projective_direction_algebra_thm3088.py
output: 05-knowledge/results/punctured_projective_direction_algebra_thm3088.out
script_sha256: 046442362be18d7a263f6273563a419fcfbc2e6cee2ed162e375b6a7239e0954
output_sha256: 95fbf0a9100ceaae48024d60c7c19712c8a5f06aa2fa9a02f8c67c5f59a3562e
hash_basis: LF-normalized bytes
---

# THM-3088 -- punctured projective-direction algebra and exceptional parity saturation

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
HOSTILE AUDIT.**

## 1. Inheritance and exact claim

[THM-3076, finite affine-plane line tomography](THM-3076-finite-affine-plane-line-quotient-tomography-and-p2-three-view-law.md)
separates the `p+1` projective directions in `F_p^2` for prime `p`.
[THM-3083, the exceptional binary/ternary clutch](THM-3083-exceptional-binary-point-ternary-direction-s4-tomography-clutch.md)
then found two apparently competing integral facts at `q=3`: the normalized
line projectors carry an `F_3^3` defect, while punctured-line idempotents give
a saturated algebra clutch.  The present theorem identifies the universal
object behind the latter fact and proves exactly why its parity saturation is
exceptional at `q=2,3`.

Let `F_q` be any finite field, put

```text
V=F_q^2,                 D=P^1(F_q),                       (1)
```

and let `k` be a commutative coefficient ring with `1!=0`.  For `L in D`,
define the punctured-line indicator

```text
e_L=1_(L\{0}).                                                (2)
```

Pointwise multiplication gives

```text
e_L e_M=delta_(L,M)e_L,       sum_(L in D)e_L=1_(V\{0}).  (3)
```

Consequently

```text
Psi_q:Fun(D,k) -> Fun(V,k),       delta_L |-> e_L          (4)
```

is a `GL_2(F_q)`-equivariant algebra isomorphism onto

```text
A_dir={f:f(0)=0 and f(ax)=f(x) for all a in F_q^*}.        (5)
```

Its unit is `1_(V\{0})`.  As a map of integral Boolean lattices, `(4)` is
saturated.

Let

```text
A_par={f:f(0)=0 and f(-x)=f(x)}.                           (6)
```

Then `A_dir` is a subalgebra of `A_par`, and

```text
A_dir=A_par                         iff q in {2,3}.         (7)
```

For larger fields the exact missing coordinate on each direction is the
finite torsor

```text
F_q^*/{+-1}.                                               (8)
```

Thus the two/three exception is not a numerical coincidence: projectivizing
forgets the full scalar group, while parity forgets only its order-at-most-two
subgroup.  These two forgetful operations agree exactly when the scalar group
has no further coordinate.

## 2. Punctured directions are the scalar-orbit quotient

Distinct one-dimensional subspaces of `V` meet only at zero.  Their punctured
parts therefore partition `V\{0}`, proving `(3)`.  Moreover, for any nonzero
`x in L`,

```text
L\{0}={ax:a in F_q^*}.                                     (9)
```

Hence a function on `V\{0}` is constant on every punctured line if and only
if it is invariant under scalar multiplication.  This proves that `(4)` has
image `(5)` and is an algebra isomorphism.  Since `GL_2(F_q)` commutes with
scalar multiplication and permutes the lines, it also proves equivariance.

Over `Z`, order the point coordinates line by line.  The columns `e_L` have
disjoint supports.  Choosing one point from every punctured line exposes an
identity `(q+1)`-minor.  Therefore their span is a primitive direct summand of
the Boolean point lattice, proving integral saturation without division by
`q`.

## 3. Exact parity quotient and the two exceptional fields

The parity subgroup of `F_q^*` is

```text
{1}                         in characteristic 2,
{+-1}                       in odd characteristic.         (10)
```

Thus a punctured line contains

```text
q-1                         parity orbits if char(F_q)=2,
(q-1)/2                     parity orbits if q is odd.      (11)
```

The direction algebra keeps only the constant function on this orbit set.
Its complement consists of the within-line contrasts.  Over any field of
characteristic zero, their total dimension is

```text
rank(A_par/A_dir)=
  (q+1)(q-2)                         if char(F_q)=2,
  (q+1)(q-3)/2                       if q is odd.            (12)
```

Both quantities vanish exactly in the stated cases:

```text
q=2: every punctured line is one point;
q=3: every punctured line is one antipodal pair.            (13)
```

This proves `(7)`.  It also gives sharp hostiles for every larger field:
choose two distinct parity orbits in one line and assign them different
values.  The resulting function lies in `A_par` but not in `A_dir`.

The first failures expose both recurring small-prime grammars without
identifying them.  At `q=4`, parity is trivial and each line retains the
three-point torsor

```text
F_4^* ~= C3.                                                (14)
```

At `q=5`, the residual quotient is

```text
F_5^*/{+-1} ~= C2.                                         (15)
```

So `C2` and `C3` reappear immediately as *internal scalar fibres* after the
exceptional saturation.  Equations `(14)--(15)` do not supply a common
`C2*C3` action, a Bass--Serre tree, or a physical modular-group carrier.

## 4. Projector normalization and the exact q-torsion

The usual centered line channel is represented by

```text
h_L=q 1_L-1_V.                                             (16)
```

Let

```text
A_D={a in Z^D:sum_L a_L=0}.                                (17)
```

For `a in A_D`, the origin and constant terms cancel, so

```text
Phi_q(a):=sum_L a_L h_L
         =q sum_L a_L e_L
         =q Psi_q(a).                                     (18)
```

The lattice `Psi_q(A_D)` is saturated of rank `q`.  Therefore

```text
Smith[Phi_q:A_D -> Psi_q(A_D)]=diag(q,...,q)  (q entries), (19)

Psi_q(A_D)/Phi_q(A_D) ~= (Z/qZ)^q.                         (20)
```

At `q=3`, `(20)` is exactly the `F_3^3` projector-normalization defect in
THM-3083.  It is not an obstruction to the punctured algebra map: `Psi_q`
itself remains integrally saturated for every `q`.

## 5. The exceptional binary/ternary handshake

The all-prime theorem isolates the two ingredients used simultaneously by
THM-3083:

```text
|F_2^2|=4=|P^1(F_3)|,
AGL_2(F_2)=S4=PGL_2(F_3),                                  (21)
```

while `(13)` makes the four ternary punctured directions exactly the four
parity atoms.  On the binary side the four-point permutation module has
three nonconstant Walsh channels; on the ternary side the four direction
atoms have a three-dimensional augmentation module.  The chosen exceptional
`S4` identification clutches those two augmentations.

What is universal is `(3)--(20)`.  What is exceptional and gauge-dependent is
the cross-field identification `(21)`.  This distinction prevents the
projector defect, the pointwise algebra clutch, and the modular `C2*C3`
quotient of [THM-2768](THM-2768-modular-c2-c3-quotients-to-a4-s4-and-bass-serre-cycle-ranks.md)
from being silently conflated.

## 6. Stopping boundaries

The theorem deliberately retains the following losses.

1. The multiplication is pointwise.  Nothing here intertwines convolution,
   Fourier multiplication, transfer operators, or a physical current.
2. `Psi_q` is unital only after the punctured ideal is given its internal unit
   `1_(V\{0})`; it is not a unital map to the full function algebra on `V`.
3. For `q>3`, the residual torsor `(8)` is real data, not an error term.  A
   parity quotient cannot reconstruct it.
4. The `q=4` and `q=5` fibres `(14)--(15)` are not by themselves the two free
   factors of `PSL_2(Z)`.  A common action, orientation, and Bass--Serre
   incidence object would still be required.
5. No canonical exceptional `S4` gauge, quartic root owner, resolvent phase,
   Feuerbach point, Farey flank, tree move, Keller graph order, LRC word,
   terminal current, or positive carrier follows.

## 7. Exact evidence

Run

```bash
python 04-computation/punctured_projective_direction_algebra_thm3088.py
python -O 04-computation/punctured_projective_direction_algebra_thm3088.py
```

Both executions byte-match the stored transcript.  The companion uses only
explicit `require` gates.  It constructs prime fields `F_2,F_3,F_5,F_7` and
quadratic fields `F_4,F_9`; enumerates every projective line, punctured scalar
orbit, parity orbit, orthogonal-idempotent product, and integral pivot; checks
the exceptional equality and both formulas in `(12)`; verifies `(18)` on an
augmentation basis; and exhausts all `6,48,180` elements of `GL_2(F_q)` for
`q=2,3,4`.

```text
PROVED IN THE CANDIDATE:
  universal punctured-direction orthogonal-idempotent algebra;
  exact scalar-invariant quotient and GL2 equivariance;
  integral saturation and projector-normalized q-torsion;
  parity saturation iff q=2 or3;
  exact larger-field internal-contrast ranks and C3/C2 first failures.

NOT PROVED:
  independent hostile audit or promotion;
  a canonical cross-field S4 gauge;
  PSL2(Z), Farey, partial-cube, or tree realization on a common carrier;
  convolutional, quartic-owner, Keller, GMC, or LRC consequence.         (22)
```

QED (candidate).
