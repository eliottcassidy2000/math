# The guard/source-clock debt cone is target-null

**Status: PROVED-ELEMENTARY exact operator-span no-go.**  This note
formalizes the smallest signed complex carried by THM-2701's recurrent
guard/source-clock debt and computes its image in the proved THM-2350/2461
target plane.  The image is exactly zero, even before ancestry is
marginalized.  This rules out a chain map generated only by the currently
retained debt operations; it does not rule out arbitrary added maps or a
new paired blocker--graft operation on the same ancestry.

## 1. The minimal signed debt complex

Retain only the two debt coordinates in THM-2701's recurrent grammar:

```text
B = guard-danger, unshifted-source-safe;
A = guard-safe, unshifted-source-danger with shifted-source-safe.       (1)
```

The target-tooth schedule and the fixed rail label remain external
decorations.  Declaring those scheduled teeth to be semantic target
coordinates would assume the chain map we are trying to construct; THM-2701
proves that the fixed `PAT_QB` rail does not change categorical target when
the delayed tooth label changes.

Over `F_13`, let `C_0=<B,A>` and let `C_1` have the two chronological arrows
`e_(BA):B->A` and `e_(AB):A->B`.  In the displayed bases,

```text
partial = [ -1   1 ]
          [  1  -1 ].                                      (2)
```

Thus `rank(partial)=1`, `dim H_0=1`, and `dim H_1=1`.  The surviving
one-cycle is exactly the signed version of the recurrent debt SCC.  If one
instead forgets chronology, the reduced kernel of the map
`<A,B> -> <omit debt>` is the single line `<A-B>`; (2) is its smallest
two-arrow mapping-torus realization.

The role increment of `B->A` is

```text
-e_H+e_j,                                                 (3)
```

because guard danger becomes safe while unshifted source safety becomes
source danger.  The reverse edge has the opposite increment.  Keeping the
guard and source repairs independent only enlarges (3) to
`span(e_H,e_j)` and will not change the calculation below.

## 2. The endpoint bank is full, but the debt subspace is zero

Use the role order

```text
(H,j,k_a,k_b,a,c),                                        (4)
```

where `j` is the source owner, `a,c` are the two targets, and `k_a,k_b`
their graft units.  THM-2461 equation (35), equivalently THM-2350's two
balanced dipoles, gives the target-covector matrix

```text
Lambda = [ 0  0 -1  0  1  0 ]
         [ 0  0  0 -1  0  1 ].                           (5)
```

The full matrix has rank two: the existing endpoint-current covectors span
the entire target plane.  But the debt-edge incidence matrix is

```text
rho = [ -1  1 ]     H
      [  1 -1 ]     j
      [  0  0 ]     k_a
      [  0  0 ]     k_b
      [  0  0 ]     a
      [  0  0 ],    c                                    (6)
```

so

```text
Lambda rho = 0.                                           (7)
```

Equation (7) is the sharp linear obstruction.  There is no shortage of
target directions globally; the particular generated operator line is
contained in their common kernel.  Every signed combination of guard and
source-clock toggles remains target-null.

## 3. The common-ancestry repair also vanishes pointwise

THM-2623 supplies the honest same-base Boolean cospan

```text
Q^safe + Q^danger = Q^(omit H).                            (8)
```

Its sector label must be retained.  The signed difference
`Q^danger-Q^safe` is a legitimate debt coordinate, but (8) supplies no
target differential.

THM-2379 gives much more than (8): on one nonnegative ancestry weight it
constructs an off-diagonal complement-repair table

```text
K_f^+(r,s),                    K_f^+(r,0)=0.                (9)
```

For a role covector `lambda_f:F_13^2->F_13`, the only lawful conditional
target pullback currently proved is

```text
Ktilde_f(r,delta)=K_f^+(r,lambda_f(delta)).                (10)
```

For both debt roles,

```text
lambda_H=lambda_j=0.                                      (11)
```

Combining (9)--(11) gives, on the same ancestry and before any integration
over the target plane,

```text
Ktilde_H(r,delta)=Ktilde_j(r,delta)=0
                  for all r,delta.                         (12)
```

The exact referee checks (12) on a basis of the entire anchored table
space.  That space has dimension

```text
13*(13-1)=156.                                            (13)
```

Pullback through either zero covector has rank zero.  Pullback through a
nonzero target covector has rank `156`, because every nonzero affine line
has thirteen-point fibres and every shift `s` is reached.  Thus the loss in
(12) is exactly target neutrality, not an accidental cancellation in one
chosen table.

This also clarifies the scope of THM-2379's positive guard coefficient.  Its
positive mixed colour is real and lives on common ancestry, but it is an
auxiliary factor-repair colour.  Composing it with the actual endpoint
quotient gives zero.

## 4. The smallest formal and lawful escapes differ

Adding one target-active graft role is enough algebraically:

```text
lambda_(k_a)=(-1,0),          lambda_(k_b)=(0,-1).          (14)
```

Either column gives a rank-one formal pullback.  It is not yet a physical
target action.  THM-2350 proves that the minimum-support lawful directions
are the balanced dipoles

```text
D_a=e_a-e_(k_a),             D_c=e_c-e_(k_b).              (15)
```

Applying (5) to (15) gives `2I_2`; since two is a unit modulo thirteen,
the paired directions span the whole target plane.  The factor two only
records the displayed role-coordinate convention; the two target lines are
the invariant content.

THM-2379 translates the graft complement while its paired blocker remains
inside the fixed ancestry weight.  Therefore it supplies (14) formally but
not (15) physically.  THM-2563 repairs this on its old-head packet by using
the paired moving-endpoint factor

```text
u_1(a x-s/13) u_(L_a)(k_a x+s/13),                         (16)
```

and proves a nonzero paired colour.  THM-2569 even places that paired corner
on common ancestry with a later event.  Neither theorem identifies its old
head with the THM-2701 debt SCC, and endpoint completion still lacks the
left residue `eta.u` and suffers the known full-`X` cancellation.  Importing
their conclusion without a new ancestry map would therefore change the
source object.

## 5. Sharp no-go and cheapest decisive test

The currently generated chain is

```text
guard/source debt cycle
  -> signed Boolean debt class on common ancestry
  -> positive off-diagonal factor-repair table
  -> zero under the target quotient.                         (17)
```

Within the operator span generated by guard complementation, source-clock
translation, reflection, and their signed compositions, (7) and (12) prove
a sharp no-go: no nonzero target-active current can be the image of the debt
class.

The cheapest genuine extension is now very specific.  On one exact positive
debt cylinder, identify the actual pivot labels `a,k_a` and compute

```text
R_(r,s)
 = integral w_debt(x) d_1(c x-r/13)
     u_1(a x-s/13) u_(L_a)(k_a x+s/13) dx.                 (18)
```

The test has three gates:

1. `w_debt` and both factors in (18) must live on one physical ancestry;
2. some nonzero `s`-character of (18) must survive; and
3. the left endpoint residue `eta.u` must be retained, so the character is
   `eta.(u-v)`, not merely the moving residue `-eta.v`.

Only thirteen shifts and the existing rational interval walls are involved.
A positive result would be the first lawful chain-map candidate from the
debt cone.  A zero result would identify the paired-dipole/full-`X`
cancellation directly on the new recurrent carrier.

No abstract linear map from the debt homology to `F_13^2` is ruled out; such
a map can always be written down.  What is ruled out is a map generated by
the proved physical debt operations and endpoint covectors.  No scalar row
is excluded and LRC(14) remains open.

## Reproduction

```bash
python3 04-computation/lrc14_guard_source_debt_cone_target_span_probe_20260728.py
python3 -O 04-computation/lrc14_guard_source_debt_cone_target_span_probe_20260728.py
```

Both runs byte-match

```text
05-knowledge/results/lrc14_guard_source_debt_cone_target_span_probe_20260728.out
```

The companion uses explicit modular Gaussian elimination and exhausts the
`156` basis tables against all `169` target shifts.  All logical guards run
unchanged under optimized mode.
