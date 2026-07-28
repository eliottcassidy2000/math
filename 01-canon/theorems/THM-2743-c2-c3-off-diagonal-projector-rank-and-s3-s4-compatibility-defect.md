---
id: THM-2743
title: "C2-C3 off-diagonal projector rank and S3-S4 compatibility defect"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For a C3 action r
  and an involution s, the operator
  (I-Pi3)sPi3 is exactly the charged part created from the equal-arm sector.
  It vanishes for every S3-normalizing reflection.  On the eight affine
  C2*C3 lifts to the Klein-four torsor, it vanishes exactly on the four
  coboundary/S3 lifts; on the four nonzero-H1/S4 lifts it has rank one and
  squared Hilbert--Schmidt norm 8/9, while the projector commutator and S3
  relation defect both have rank two.  This is a finite marked-resolvent
  compatibility detector and an exact equal-arm mixing gate, not a Keller
  exclusion or an LRC endpoint-current construction.
source: a4-resolvent-next-gate-scout-2026-07-28
audit: thm2694-full-lift-fibre-scout-2026-07-28 (independent theorem, matrix, group-closure, scope, normal/-O/stored, and hash audit: ACCEPT)
depends_on:
  - THM-2595-modular-v4-affine-lift-dichotomy-and-six-vertex-tournament-no-go
related:
  - THM-2596-modular-free-factor-farey-gram-owner-cocycle
  - THM-2696-reflection-completed-s4-relative-different-and-coordinate-invariant-jacobian-gate
  - THM-2711-c3-stable-discriminant-metabolizer-and-mod-four-isotropy-gate
  - THM-2714-c3-metabolic-length-parity-and-unique-plane-four-adic-escape
  - THM-2716-c4-arm-transporter-groupoid-and-relative-degree-holotopy-boundary
  - THM-2721-semantic-inner-triangle-equal-following-amplitude-and-current-reanchoring-no-go
  - THM-2741-highest-odd-faber-response-pole-capacity-closure
script: 04-computation/c2_c3_projector_compatibility_thm2743.py
output: 05-knowledge/results/c2_c3_projector_compatibility_thm2743.out
script_sha256: bf76bbcb60b5129643b9526bc009043455070cfad93532462d6a9e69d9128fff
output_sha256: 09582cce4a676a594241fd8a7d34f322f855bca31a1fe6343fbd34c0cc86c9e5
hash_basis: LF-normalized bytes
---

# THM-2743 -- C2-C3 off-diagonal projector rank and the S3/S4 defect

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The modular free product has one binary and one ternary generator, but passing
to its smallest symmetric quotient imposes an extra compatibility relation.
This theorem isolates that lost coordinate.  The off-diagonal Reynolds block

```text
                 (I-Pi_3) s Pi_3                              (1)
```

measures whether the binary move turns a ternary equal-arm vector into a
charged ternary vector.  On the affine Klein-four torsor of THM-2595 it is
exactly the cohomological distinction between the split `S3` lift and the
transitive `S4` lift.

## 1. The general projector gate

Let `k` be a field of characteristic different from three, let `E` be a
finite-dimensional `k`-space, and suppose

```text
r^3=I,                         s^2=I.                       (2)
```

Put

```text
Pi=(I+r+r^2)/3,                Q=I-Pi,
M=Q s Pi.                                                   (3)
```

Then `Pi` is the Reynolds projection onto `Fix(r)`, and

```text
M=0  iff  s(Fix(r)) is contained in Fix(r).                (4)
```

In particular, if the two moves factor through the dihedral quotient

```text
s r s = r^(-1),                                             (5)
```

then `sPi=Pi s` and `M=0`.  If `v` is an equal-arm vector, meaning `Pi v=v`,
then

```text
Q(sv)=Mv.                                                   (6)
```

Thus `(1)` is exactly the charged component created by the binary move.

The converse to `(5)` is false for arbitrary representations: preservation
of one fixed subspace need not determine how `s` conjugates `r`.  If `E` has
an invariant positive-definite inner product and `s,r` are orthogonal, then
`M=0` is equivalent to `[s,Pi]=0`; even there it need not imply `(5)` without
faithfulness or additional structure.  The equivalence below is a special
finite affine-torsor theorem, not a formal consequence of `(3)`.

## 2. The eight affine lifts

Use the Klein four group `V=F_2^2` and the matrices from THM-2595,

```text
S = [0 1],                 T = [0 1],
    [1 0]                      [1 1].                       (7)
```

For `a,b in V`, define affine permutations

```text
sigma_a(v)=Sv+a,             tau_b(v)=Tv+b.                (8)
```

The order conditions `sigma_a^2=tau_b^3=1` leave exactly eight pairs
`(a,b)`.  Translation of the origin by `v in V` changes the pair by the
coboundary

```text
(a,b) -> (a+(I+S)v, b+(I+T)v).                             (9)
```

There are four coboundaries and four representatives of the unique nonzero
class in

```text
H^1(C2*C3,V) = F_2.                                        (10)
```

Let `E=Q[V]` with its delta-function orthonormal basis.  Write `Sigma,Tau`
for the permutation matrices of `(8)` and

```text
Pi_b=(I+Tau+Tau^2)/3,
M_(a,b)=(I-Pi_b) Sigma Pi_b.                               (11)
```

The following conditions are equivalent for every one of the eight affine
lifts:

```text
(a,b) is a coboundary;
<sigma_a,tau_b> has order 6 and orbit type 1+3;
sigma_a tau_b sigma_a = tau_b^(-1);
rank M_(a,b)=0;
rank [Sigma,Pi_b]=0.                                       (12)
```

For each of the four nonzero cohomology lifts,

```text
<sigma_a,tau_b> = S4,             orbit type 4,
rank M_(a,b) = 1,
||M_(a,b)||_HS^2 = 8/9,
rank [Sigma,Pi_b] = 2,
rank(Sigma Tau Sigma-Tau^2) = 2.                          (13)
```

Consequently the smallest affine realization of the free binary/ternary
grammar has a one-dimensional equal-arm-to-charged leakage precisely when it
does **not** descend to the `S3` quotient.

## 3. Proof

The identity

```text
Pi^2=Pi,                 image(Pi)=Fix(r)                  (14)
```

follows from `r^3=1`.  For `y in Fix(r)`, equation `(3)` gives `My=Qsy`,
which vanishes exactly when `sy` is again fixed.  This proves `(4)` and
`(6)`.  Under `(5)`,

```text
s Pi s=(I+r^(-1)+r^(-2))/3=Pi,                             (15)
```

proving the normalizer boundary.

THM-2595 supplies the eight cocycles, the four coboundaries, and the
`S3/S4` image dichotomy.  In the coboundary class, translating the origin
conjugates the lift to the linear `S3` action, so `(5)` holds.  In the
nonzero class all four lifts are translation-conjugate.  It is therefore
enough to calculate after relabelling the four torsor points so that

```text
tau=(1 2 3),                    sigma=(0 1).               (16)
```

The fixed space of `tau` is spanned by

```text
c=e0+e1+e2+e3,             w=3e0-(e1+e2+e3).              (17)
```

The constant vector is killed by `M`.  Direct orthogonal projection gives

```text
||M w||^2 / ||w||^2 = (32/3)/12 = 8/9.                    (18)
```

Thus `M` has rank one and squared Hilbert--Schmidt norm `8/9`.  The same
four-by-four matrices give the two rank-two defects in `(13)`.  Translation
conjugacy preserves all these ranks and the Hilbert--Schmidt norm.  This
proves `(12)--(13)`.

## 4. Quartic-resolvent interpretation

The four roots of a quartic form a `V4` torsor after choosing one root, while
the cubic resolvent records the quotient `S4/V4 isomorphic to S3`.  Equations
`(12)--(13)` therefore provide an exact finite marked-fibre detector:

```text
split affine lift / S3 quotient       <=> off-diagonal rank 0,
transitive affine lift / S4 extension <=> off-diagonal rank 1.             (19)
```

This is the precise part of the proposed quartic--resolvent transfer that
survives scrutiny.  It says where the missing `V4` extension class sits: not
in the cubic equal-arm sector itself, but in its binary off-diagonal leakage.
THM-2741's exact real/imaginary sixth-power response is a useful quadrature
analogy, but neither theorem supplies the other's geometric realization.

It does **not** construct a degree-four graph quartic for a Keller map,
identify affine Jelonek inertia with the abstract generators, impose the
grade-three discriminant/cuspidal laws on that quartic, or exclude the
transitive `S4` case.  THM-2696 still requires divisor, source-unit,
specialization, and polynomial-graph sidecars.  Hence `(19)` is not a proof
of JC(2), DC(2), or any new degree-four Keller exclusion.

## 5. Equal-arm LRC interpretation

THM-2721 produces a positive three-arm corolla whose two charged abstract
`C3` transforms vanish.  If a lawful physical involution acted on the same
amplitude space, equation `(6)` would say exactly when that equal-arm mass
could be converted into charged mass:

```text
nonzero charged output  <=>  (I-Pi_3)sPi_3 is nonzero on that vector.       (20)
```

The `S3` reflection is the sharp no-go: it normalizes the ternary action and
cannot create charge.  A genuinely free-product/S4-type compatibility
defect can create only one charged direction in the minimal torsor model.

No such common physical action is currently present in canon.  THM-2716's
binary transporter and THM-2721's ternary corolla live on different carriers,
and THM-2721 also proves that its three endpoints miss the inherited current.
The missing sidecar is therefore concrete:

```text
one common amplitude space;
one lawful physical involution on it;
a nonzero off-diagonal block on the retained equal-arm vector;
an endpoint-current/owner typing for the charged output.                   (21)
```

Without `(21)`, the theorem supplies neither a row exclusion nor an LRC(14)
conclusion.  The scalar ledger is unchanged.

## 6. Exact reproduction

Run

```bash
python 04-computation/c2_c3_projector_compatibility_thm2743.py
python -O 04-computation/c2_c3_projector_compatibility_thm2743.py
```

Both executions byte-match the stored `15`-line transcript
`05-knowledge/results/c2_c3_projector_compatibility_thm2743.out`.  The
companion uses exact rational arithmetic and explicit exceptions, with no
optimized-away assertions.  It independently enumerates all eight affine
cocycles and four coboundaries, constructs both permutation generators and
the Reynolds projector, verifies the isolated subgroup orders `2` and `3`
before enumerating each joint image, and checks every rank, norm, and
relation-defect value in `(12)--(13)`.

## 7. Boundary ledger

```text
PROVED HERE:              abstract projector gate (4)--(6);
                          exact eight-lift S3/S4 rank table (12)--(13);
                          one-dimensional nonzero-class leakage.

NOT PROVED:               a Keller quartic carries the marked affine action;
                          S4 monodromy is excluded;
                          the LRC binary and ternary carriers are identified;
                          charged leakage reaches an endpoint current;
                          JC(2), DC(2), or LRC(14).                         (22)
```

QED.
