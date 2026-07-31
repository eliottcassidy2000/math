# The missing even edge is a removable response pole

**Status:** PROVED as THM-2752 + VERIFIED-EXACT + INDEPENDENTLY
HOSTILE-AUDITED.

## Inheritance and concept board

The closest proved mechanism is THM-2723: the third response on a polynomial
split exact prefix is a rational primitive with exactly one source pole.  The
canonical hostile is its `q=0` example `P=z^4+x, Q=z`; one-pole capacity alone
does not empty a split chart.  The corrected near miss is THM-2741's
highest-odd face: when every odd seed vanishes, its transverse `G2` response
pole vanishes with it.  The least-used sidecar is not a new coefficient but
the homogeneous combination `R+(d/2)Phi`.

The working board was:

| object | retained predicate | lost coordinate / obstruction | cheapest test |
|---|---|---|---|
| all-even germ at `P_infty` | two flux equations | odd Newton face | exact local valuations |
| third response | `U R'=kappa` | top response cancels `Phi` | compute `R_m+dPhi_m/2` |
| exact-prefix lift | polynomial regularity at finite source poles | forgotten `C,E` residues | smooth-point resultant |
| projective component | every component meets `h=0` | reducibility/nilpotents | normalize only the generic-image component |
| old `N_1=N_2` support atlas | exact coefficient equations | nonsplit terminal typing | treat only as a secondary control |

## The structural calculation

For the five even degrees

```text
m=2,6,10,14,22,
```

the local orders of both first and second Faber fluxes at `P_infty` are

```text
k_m=1,2,3,4,6.
```

They satisfy the single slope-four invoice

```text
(22-m)+4k_m=24.                                        (1)
```

This explains why the odd proof disappears: all lower even teeth land on the
same Newton face, rather than selecting one transverse branch by its highest
index.

The response contains one more order.  Globally, with weights
`wt(h,d,q,s)=(1,2,3,4)`, the typed combination is

```text
D_m=R_m+(d/2)Phi_m,

D_2=0,
D_6=-q^3/4,
D_10=5q^3s/8,
D_14=7q^3(dq^2-5s^2)/32,
D_22=-33q^3(6d^2q^4-84dq^2s^2+3q^4s+70s^4)/1024.     (2)
```

Thus `D_m` has local order at least `k_m+1`.  For zero first flux,

```text
R_25+(d/2)F_23=sum_m a_m h^(22-m)D_m.                 (3)
```

The factor is `d`, not `h^2`; both have weight two, but only `d` is the exact
Faber cancellation.

There is an all-degree reason for the common `q^3`.  For every
`m=4ell+2`, expand
`(A(z^2)+qz^3)^(m/4)` with
`A(y)=1+2dy+(d^2-s)y^2`.  Parity leaves only odd powers of `q` in
`R_m+(d/2)Phi_m`.  If `A(y)^beta=sum b_ny^n`, with
`N=m/2` and `beta=(N-2)/2`, its linear-`q` coefficient is proportional to
`b_N+d b_(N-1)`.  The coefficient recurrence from
`A(A^beta)'=beta A'A^beta`, evaluated at `n=N-1`, is exactly
`N(b_N+d b_(N-1))=0`.  Thus the response combination is divisible by `q^3`
for every such `m`.  This explains the syzygy but not the degree-22
slope-four alignment needed below.

## Why every last-end branch has slope at least four

On the finite `d=1` index cover let `e=v(h)`, `u=v(q)`, and `v_0=v(s)`.
If `a=min(u,v_0)<4e`, every lower even term has order

```text
(22-m)e+k_m a
 =6a+(6-k_m)(4e-a)>6a,                                (4)
```

and `W h^24` is also higher.  The initial system is therefore the pair of
degree-six top faces.

Unequal valuations cannot cancel: `Psi_22` has a unique pure `q^6` lowest
term when `u<v_0`, and a unique pure `s^6` lowest term when `v_0<u`.  At equal
valuation the two homogeneous faces are

```text
q s(q^2-3s^2)(3q^2-s^2),
(q^2-s^2)(q^4-14q^2s^2+s^4),                          (5)
```

and have gcd one.  Hence

```text
min(v(q),v(s))>=4v(h).                                 (6)
```

The cases with one coordinate identically zero are included in the unequal
argument.  If both vanish identically, every `D_m` vanishes directly.  A
residual `mu_2` may exchange the two index-cover branches, but a finite cover
multiplies every valuation by the same ramification index; `(6)` and
regularity descend.

On `F_23=0`, equations `(1)`--`(3)` now give

```text
v(R_25)>=28v(h),
v(R_25/h^25)>=3v(h)>0.                                 (7)
```

The response is not merely pole-free at `P_infty`; it vanishes there.

## Why this closes every component

The coefficient-independent top intersection consists of `P_infty` and five
transverse coarse points.  The latter are governed by the squarefree quintic

```text
20141047808z^5-14386462720z^4+1089822720z^3
-21288960z^2-35910z+81,
```

and the top response is nonzero at each, so each is a response pole of order
`25`.  This top calculation was repeated without importing an odd-seed
theorem under false premises.

Normalize the reduced irreducible component containing the generic physical
image.  The source is reduced, so nilpotents are killed; reducibility elsewhere
is irrelevant.  THM-2723 gives a finite surjective map from the source `P1`
and exactly one source response pole.

If the component contains two smooth top points, their disjoint pole fibres
contradict that capacity.  If it contains exactly one, pole order `25` forces
the source pole to be finite.  The exact-prefix lift equations

```text
beta^2+d=h^2C,                  q beta-s=h^4E          (8)
```

then have regular right sides.  At the smooth top point they reduce to
`d=-omega^2`, `s=q omega`; substitution in the two top fluxes has resultant

```text
-76608 omega^6,
```

which is nonzero because the smooth point has `d!=0`.  Thus no smooth top
point can occur.

Every remaining infinity point lies over `P_infty`, where `(7)` makes the
response regular.  It is already regular on `h!=0`, so it is a global regular
function on a projective integral curve and is constant.  Pullback contradicts
`U R_source'=kappa!=0`.

## Secondary control deliberately not used

There is a plausible second route through THM-2411's identities
`Phi_Q=-qN_1/7496192` and `Psi_Q-W=N_2/1319329792`: split into `q=0` and
`q!=0`, then try to reuse the support-three/four/five `N_1=N_2=0` atlas.
This is useful as an audit target, but it is not a dependency of THM-2752.
Those support theorems are stated in the genuine-nonsplit branch and end by a
nonsplit terminal contradiction.  Retyping every terminal step with the
third response would require a separate line audit.  The local regularity
proof avoids that transfer entirely.

## Boundaries and reusable mechanism

A nonzero first flux inserts `-(lambda d/2)h^23` into `(3)` and destroys the
order-28 gain; THM-2725 handles it.  A nonzero odd seed destroys the even
cancellation and creates the transverse response pole handled by THM-2745.
Constant Faber coefficients and the polynomial exact prefix remain
load-bearing.

The reusable pattern is: when a supposedly decisive response vanishes on a
boundary face, compute its exact syzygy with the defining equation.  One
extra local order can turn a failed pole-count proof into a projective
regular-function contradiction.

THM-2725 + THM-2745 + THM-2752 empty the entire chosen-sheet
split degree-22 response family.  They do not derive this family for every
Keller pair, perform integral raising, close the broader split branch, or
prove `JC(2)` / `DC(2)`.
