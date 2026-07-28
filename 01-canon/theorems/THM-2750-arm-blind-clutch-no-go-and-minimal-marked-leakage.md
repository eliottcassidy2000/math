---
id: THM-2750
title: "Arm-blind clutch no-go and minimal marked leakage"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On a common C3
  arm carrier, every transition commuting
  with the arm rotation has zero invariant-to-charged block.  Hence the
  scalar root-normalized gain -1 of THM-2744 can reverse an already charged
  class but cannot charge THM-2721's equal corolla.  The minimal signed
  escape is an unequal arm gain; the minimal positive permutation escape
  adds one C3-fixed reference and moves it into an arm.  The bare four-point
  leakage detects both A4 and S4 fixed-point-moving involutions, in agreement
  with THM-2746/2743, and does not distinguish them.
source: clutch-c3-holotopy-2026-07-28
audit: root-zero-clutch-audit/arm-blind-c3-hostile-audit-2026-07-28 (independent centralizer/normalizer, sign/group, four-point involution, field/lattice, cross-theorem scope, and normal/-O/hash replay: ACCEPT; post-audit script change is docstring-only)
depends_on: []
related:
  - THM-2567-deep-coloured-duty-replica-cycle-and-augmentation-cancellation
  - THM-2657-odometer-carry-root-lift-nonsplit-extension-and-cech-cocycle
  - THM-2716-c4-arm-transporter-groupoid-and-relative-degree-holotopy-boundary
  - THM-2721-semantic-inner-triangle-equal-following-amplitude-and-current-reanchoring-no-go
  - THM-2743-c2-c3-off-diagonal-projector-rank-and-s3-s4-compatibility-defect
  - THM-2744-relative-present-unit-repair-and-root-zero-overlap-clutch
  - THM-2746-c3-quotient-affine-v4-lifts-and-a4-projector-defect
  - THM-2749-fully-marked-root-zero-clutch-and-target-character-profile
script: 04-computation/lrc14_arm_blind_clutch_no_go_thm2750.py
output: 05-knowledge/results/lrc14_arm_blind_clutch_no_go_thm2750.out
script_sha256: 0104ce30b5f17697a767e277ac5ba66e403c75eecfed04732625d676bc3be699
output_sha256: 2e3ae1472dc285d4b1cc359567b567c8a1d08e60d882e4aa777e075c9be7beb2
hash_basis: LF-normalized bytes
---

# THM-2750 -- arm-blind clutches preserve the equal-arm sector

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2721 gives a positive equal three-arm corolla.  THM-2744 gives a genuine
open physical chart overlap whose recomputed raw coefficient vectors agree,
but whose root-normalized classes differ by the sign `-1`.  THM-2743 shows
that a binary move can create charged ternary amplitude only through the
off-diagonal block

```text
                         (I-Pi_3)sPi_3.                  (1)
```

This candidate computes their abstract composition.  The sign is central,
so it cannot create charge.  It also proves the sharp first signed and
positive repairs and reconciles their four-state permutation census with the
separate proved `A4` and `S4` marked lifts of THM-2746/2743.

## 1. The arm-blind clutch theorem

Let `k` be a field in which `3` is invertible, let `W` be a `k`-vector
space, and put

```text
V=k[C3] tensor W,
r:e_0 -> e_1 -> e_2 -> e_0,
Pi=(I+r+r^2)/3 tensor I_W,
Q=I-Pi.                                                 (2)
```

The image of `Pi` is the equal-arm sector

```text
k(e_0+e_1+e_2) tensor W.                               (3)
```

### Theorem

If `C in End_k(V)` commutes with `r`, then

```text
                         Q C Pi=0.                       (4)
```

This applies in particular to every operator

```text
C=sum_(j=0)^2 r^j tensor T_j,                           (5)
```

every arm-blind auxiliary transition `I_3 tensor T`, every scalar Cech gain,
and every product or holonomy of such gains.  It also applies to every
invertible `C` normalizing `<r>`, because a normalizer permutes `r,r^2` and
hence fixes `Pi`.

### Proof

Commutation with `r` gives commutation with the polynomial `Pi`.  Therefore

```text
Q C Pi=Q Pi C=0.
```

For `(5)`, if `u=e_0+e_1+e_2`, then more explicitly

```text
C(u tensor w)=u tensor sum_j T_j(w),                    (6)
```

so the invariant sector maps into itself.  For a scalar group-ring element
`f(r)` one has the sharper augmentation identity

```text
f(r)Pi=f(1)Pi,                    Q f(r)Pi=0.             (7)
```

This proves `(4)--(7)`.  `QED.`

The Cech interpretation is exact.  Scalar transition functions can define a
nontrivial line local system on `(3)`, but their holonomy is still block
diagonal in the `C3` isotypic decomposition.  A nontrivial scalar cocycle is
not a connecting map from the invariant line to either charged character.

## 2. Consequence for the THM-2744 sign

THM-2744 proves, in the primitive `C7` coefficient algebra over `F_13`, that
the source and target root-normalized classes are

```text
11 and 2=-11.                                            (8)
```

Thus its normalized chart gain is the scalar `-1`.  If one hypothetically
placed that line on a common three-arm carrier, its action would be

```text
C=-I_arm tensor I_coefficient,
```

and `(4)` would give

```text
(I-Pi_3)C Pi_3=0.                                       (9)
```

All six nontrivial `C7` clock characters surviving THM-2744 may remain in
the auxiliary factor `W`; none becomes an external charged `C3` arm
character.  If a noncentral gate `s` had already created charge, then

```text
Q(-I)sPi=-QsPi.                                         (10)
```

The sign can orient or reverse an existing charged class, but cannot create
one.

This conclusion is stronger than the scalar example.  Any arm-independent
clock, root-gauge, odometer, or coefficient operator obeys `(4)`.  Every
scalar-gained permutation of the three arms is also flat, because all six
three-arm permutations fix the equal vector.

## 3. Coefficient-field and lattice boundary

The canonical objects are not yet coefficients in one common module.
THM-2721's exact positive rational amplitude is

```text
A=2089891528601250990/1792160394037
 =2089891528601250990/13^11.                            (11)
```

It has no naive reduction modulo `13`.  Clearing the displayed common
denominator gives numerator residue `5`, but this is a choice of integral
lattice and is not a canonical identification with the THM-2744 normalized
line.  In particular, positivity over `Q` and root-normalized unitality over
`F_13` must not be silently merged.

The abstract no-go is independent of this typing problem: after **any**
choice of common field and module on which the transition remains arm-blind,
equation `(4)` still forces zero leakage.  A canonical application would need
both the common lattice/carrier and an arm-dependent map.

## 4. Minimal signed escape: unequal arm gains

Equations `(12)--(14)` remain valid over any field in which `3` is
invertible.  The sign, positivity, Euclidean Hilbert--Schmidt, and finite
permutation censuses from this point through Section 6 are over `Q`; the
displayed charged Fourier calculation is separately over `F_13`.  This scope
is necessary: in characteristic two, `-1=1` and every sign pattern below
collapses to the constant pattern.

Let

```text
D_g=diag(g_0,g_1,g_2),       gbar=(g_0+g_1+g_2)/3.        (12)
```

On the three-arm carrier,

```text
Q D_g Pi=(g-gbar u)u^T/3,           u=(1,1,1).           (13)
```

It has rank zero exactly when `g_0=g_1=g_2`, and rank one otherwise.  On an
equal input,

```text
Q D_g Pi(Au)=A(g-gbar u).                                (14)
```

For the one-arm sign `g=(-1,1,1)`,

```text
Q D_g Pi=
 [-4/9 -4/9 -4/9]
 [ 2/9  2/9  2/9]
 [ 2/9  2/9  2/9],                                      (15)

Au -> A(-4/3,2/3,2/3).                                  (16)
```

The operator has rank one and squared Hilbert--Schmidt norm `8/9`.  Over
`F_13`, with primitive cube root `omega=3`, the normalized Fourier
coefficients of `(-1,1,1)` are

```text
(trivial,chi,chi^2)=(9,8,8),                             (17)
```

so both charged characters survive.

Among all eight sign patterns, precisely the two constant patterns have
rank zero.  The six nonconstant patterns all have rank one and squared norm
`8/9`.  Together with the arm cycle, a two-minus pattern generates `A4` of
order `12`; a one-minus pattern generates `C2 x A4` of order `24`, with
central `-I`, rather than `S4`.

Thus `(15)` is the smallest signed linear escape, but it is not a positive
permutation transition.  THM-2744 supplies a common sign on its whole
coefficient line, not the marked pattern `(-1,+1,+1)` on three semantic
arms.  Producing `(15)` requires precisely an arm-dependent chart selector
on one common carrier.

## 5. Minimal positive permutation escape

No permutation of exactly three arm coordinates charges their equal vector.
The smallest positive bijective escape therefore adds one `C3`-fixed
reference coordinate.  On four points put

```text
r=(1 2 3),                    s=(0 1),
Pi_4=(I+r+r^2)/3.                                      (18)
```

Then `<r,s>=S4`, and

```text
M_4=(I-Pi_4)sPi_4
 = [  0    0     0     0  ]
   [ 2/3 -2/9  -2/9  -2/9]
   [-1/3  1/9   1/9   1/9]
   [-1/3  1/9   1/9   1/9].                            (19)
```

This again has rank one and squared Hilbert--Schmidt norm `8/9`.  Its action
on a fixed reference plus an equal corolla is

```text
(x,A,A,A)
 -> (0,2(x-A)/3,(A-x)/3,(A-x)/3).                      (20)
```

Hence the exact missing sidecar is

```text
one C3-fixed reference amplitude x on the same carrier;
a proof that x!=A;
a lawful involution exchanging that reference with one named arm.       (21)
```

If a future common-carrier theorem typed the THM-2744 sign as the fixed
reference value `x=-A`, then `(20)` would give

```text
(0,-4A/3,2A/3,2A/3),                                    (22)
```

the same charged profile as `(16)`.  This is a sharp conditional design,
not a current application: THM-2744 proves a scalar transition on a
different coefficient carrier, not a fourth coordinate or a swap with one
THM-2721 arm.

## 6. The complete four-point involution boundary

There are ten involutions on the four-point carrier of Section 5.  Their
exact census is

```text
fixed point retained:
  identity, group C3,                         rank 0;
  three arm transpositions, group S3,         rank 0;

fixed point moved:
  three double transpositions, group A4,      rank 1, HS^2=8/9;
  three fixed-arm transpositions, group S4,   rank 1, HS^2=8/9.          (23)
```

Thus on the bare permutation module the off-diagonal block detects whether
the involution moves the unique `C3` fixed coordinate.  It does not
distinguish `A4` from `S4`.

This independently reconstructs the boundary already separated in proved
canon.  THM-2746 treats the internal `V4` double translation and obtains the
`A4` row.  THM-2743 treats the nonzero compatibility class of the marked
affine `S3` lift and obtains the `S4` row.  Their binary linear shadow and
cohomological source, not the common rank/norm, distinguish the branches.
The special `S3/S4` equivalence of THM-2743 must not be exported to arbitrary
four-point actions.

## 7. Common-carrier audit

The nearest proved sidecars do not supply `(15)` or `(21)`.

| sidecar | retained datum | first missing datum |
|---|---|---|
| THM-2744 | open physical root-chart overlap, equal recomputed raw vectors, scalar normalized gain `-1`, one semantic/full-target cylinder witness | no three-arm index, global involution, or endpoint current |
| THM-2721 | one common source cylinder and exact positive `(A,A,A)` | comparison `C3` is abstract; chronological current reanchoring is empty |
| THM-2716 | binary transporter, unbased `C2` torsor, and sign line | no functor to the THM-2721 ternary amplitude carrier |
| THM-2657 | nonsplit odometer kernel and scalar Cech cocycle | central character values are arm-blind; pairwise overlap is not a common three-arm simplex |
| THM-2567 | exact coloured-face/augmentation-cancellation warning | its deep `C13` colour and quotient-duty carrier are not the corolla/root carrier |
| THM-2746 | exact positive four-state `A4` leakage model | no common physical LRC binary/ternary torsor or endpoint current |
| THM-2749 | proved two-sided semantic and primitive-target refinement of the root overlap | its `C13` target colour is not an external `C3` arm selector |

Even the promoted THM-2749, with equal source/target raw vectors, scalar sign,
and all primitive `C13` target characters, does not by itself defeat `(4)`.
Those target characters live in `W`; the missing map is the off-diagonal
operator in the external three-arm factor.

## 8. Cheapest decisive physical tests

The theorem leaves two sharply different tests.

1. **Three-arm signed test.** Rebuild three root-overlap coefficient charts
   on one literal semantic/current cylinder and compute their normalized gain
   triple `g`.  If it is constant, `(4)` proves every charged output zero.  If
   it is nonconstant, `(13)--(14)` give the charged class immediately.  The
   cheapest target pattern is `(-1,+1,+1)`.
2. **Four-state positive test.** Construct one `C3`-fixed reference
   coefficient `x` on the same support as the three equal arms, prove `x!=A`,
   and build a lawful positive involution moving it into one named arm.
   Equation `(20)` then gives the charged coefficient.

Either test must still retain endpoint-current and owner typing after the
charged projection.  Otherwise it reaches only THM-2743/2746's abstract
finite-fibre boundary.

## 9. Exact companion and boundary ledger

Run

```bash
python 04-computation/lrc14_arm_blind_clutch_no_go_thm2750.py
python -O 04-computation/lrc14_arm_blind_clutch_no_go_thm2750.py
```

Both executions byte-match the stored transcript.  The dependency-free
companion uses exact rational arithmetic and explicit exceptions, with no
optimized-away assertions.  It checks `343` group-ring controls, a
non-scalar arm-blind auxiliary operator, all six three-arm permutations, all
eight sign patterns over `Q`, the charged sign Fourier values over `F_13`,
the standard four-point `S4` matrix and general input formula, all ten
four-point involutions including the `A4/S4` split, and the fact that central
`-1` only negates existing leakage.  Normal, optimized, and stored outputs
agree after LF normalization, and the declared hashes match.

```text
PROVED HERE:              arm-blind/normalizer clutch no-go;
                          exact unequal-arm rank-one formula;
                          minimal positive four-state formula;
                          complete A4/S4 involution census;
                          exact coefficient-field warning.

NOT CONSTRUCTED:          a common THM-2721/2744 carrier or lattice;
                          a nonconstant three-arm gain in the physical bank;
                          a fixed reference and lawful reference-arm swap;
                          endpoint current, owner transport, row exclusion,
                          or LRC(14) proof.                              (24)
```

An independent hostile audit rederived the centralizer and normalizer gate,
the unequal-diagonal rank formula, all sign/group and four-point involution
censuses, the characteristic and lattice boundaries, and the scope against
THM-2743/2746/2749.  Its normal and optimized runs byte-matched the stored
`28`-line transcript.  The final script hash differs only because the stale
"scratch evidence" docstring was changed to "canonical exact companion"
after that audit; executable code and output are unchanged.

The LRC residual ledger remains `165`.  `QED.`
