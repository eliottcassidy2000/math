---
id: THM-4173
title: "Repeated-top row-A complete planar Jacobian exclusion"
status: >
  PROVED RELATIVE TO THM-3827/3992/3997/4007/4103/4120/4122/4130/
  4147/4157/4171 + VERIFIED-EXACT + INDEPENDENTLY BRIDGE-AUDITED. The
  entire exact-weight-nine repeated-top row A is excluded, with no residual
  discriminant or universal-fibre-separation hypothesis. Repeated projected
  roots either carry several reduced critical points or force a zero
  Keller Hessian; neither loses affine critical length.
source: codex-frontier-synthesis-creative-20260826au
depends_on:
  - THM-3827-generic-fibre-genus-floor-for-nonlinear-cubic-plane-atlases
  - THM-3992-reduced-two-three-cusp-jet-repair-and-first-node-residual
  - THM-3997-reduced-two-three-hasse-repair-and-zero-residual-no-go
  - THM-4007-live-two-three-third-normal-row-five-weight-floor
  - THM-4103-jc23-theta-boundary-ramification-and-degree-response
  - THM-4120-jc23-extremal-target-degree-twenty-one-response
  - THM-4122-planar-keller-asymptotic-width-and-resonant-shear-contraction
  - THM-4130-theta-only-extremal-seam-critical-monodromy-obstruction
  - THM-4147-generic-exact-weight-nine-planar-jacobian-monodromy-exclusion
  - THM-4157-repeated-top-edge-wall-planar-jacobian-exclusion
  - THM-4171-row-a-inner-resultant-wall-planar-jacobian-exclusion
related:
  - THM-4155-generic-y-only-delta-zero-weight-nine-planar-jacobian-exclusion
  - THM-4159-inner-resultant-wall-planar-jacobian-exclusion
  - THM-4161-y-only-double-top-root-planar-jacobian-exclusion
  - THM-4164-y-only-triple-top-root-planar-jacobian-exclusion
  - THM-4165-y-only-inner-top-triple-intersection-planar-jacobian-exclusion
  - THM-4176-complete-repeated-top-wall-planar-jacobian-exclusion
script: 04-computation/jc23_repeated_top_row_a_complete_exclusion_thm4173.py
output: 05-knowledge/results/jc23_repeated_top_row_a_complete_exclusion_thm4173.out
independent_audit_script: 04-computation/jc23_repeated_top_row_a_complete_exclusion_thm4173_independent_bridge_audit.py
independent_audit_output: 05-knowledge/results/jc23_repeated_top_row_a_complete_exclusion_thm4173_independent_bridge_audit.out
script_sha256: b7fc2c372e398aff9181f001a0c20bd794a112b1a426dac288148f2a3d0ce8e4
output_sha256: 4a43b0f333c2cc2d8dbc7373d3f624c382ac5e9be4bbe85d9fbd41f6dbb02ba2
independent_audit_script_sha256: 6d7c01a3396e802a538c3a5530dcbf203df883d44c4da9c926635ef3a4b18810
independent_audit_output_sha256: 95810cae502f9a17d083906d4316088b73df4f30ffe21fd2f9ba960384b43625
hash_basis: raw LF bytes
primary_audit: >
  PASS. A self-contained exact symbolic audit reconstructs the complete
  normalized source, proves the differential bridge, verifies the two
  universal Morse pairs and leading rows, and independently recomputes the
  full X-resultant at an exact off-wall row-A control. Normal, optimized,
  and hash-seeded outputs byte-match. It deliberately does not compute any
  discarded discriminant or Q_19(-1/6).
independent_audit: >
  ACCEPT. A clean-room generic-coefficient referee proves the bridge identity
  and contrasts two reduced points sharing one projected coordinate with a
  single doubled point carrying the same resultant divisor. It independently
  checks both response budgets. Normal, optimized, and hash-seeded outputs
  byte-match.
---

# THM-4173 -- repeated-top row-A complete planar Jacobian exclusion

**PROVED RELATIVE TO THM-3827/3992/3997/4007/4103/4120/4122/4130/
4147/4157/4171 + VERIFIED-EXACT + INDEPENDENTLY BRIDGE-AUDITED;
JC(2) REMAINS OPEN.**

Work over `C` in the live `b=d=0` reduced `(2,3)` seam at exact residual
weight nine. Put

```text
P=T+X^2T^2,                 Y=XTP,
K=2848/45-(7/6)Delta,
```

and use the complete normalized source

```text
G=-X^2T/2-3P+(8/3)P^2-(1376/135)P^3
  +K Y^2+Phi P^2Y+Delta P^4+Theta P Y^2
  +eta P^3Y+zeta Y^3.                                  (1)
```

Set

```text
C=Delta+Theta,                 D_A=4Theta K^2-27eta^2. (2)
```

## 1. Theorem and inheritance pass

> **Theorem.** The complete repeated-top row-A locus
>
> ```text
> zeta=-eta,                  eta*Delta*C != 0           (3)
> ```
>
> contains no nonautomorphic planar Keller pair.
>
> No hypothesis is imposed on `Phi`, `D_A`, either residual discriminant,
> or the value `Q_19(-1/6)`.

The closest proved mechanism is
[THM-4157](THM-4157-repeated-top-edge-wall-planar-jacobian-exclusion.md),
which gives the row-A source, packet, critical resultant on its open chamber,
and carrier response. The canonical hostile is a collision of several
critical points in one projected `T`-fibre. The corrected near miss was the
implication

```text
fewer distinct projected T-values
  ==> fewer affine critical points / smaller critical length.          (4)
```

It is false: a repeated resultant root may record several distinct reduced
points. The least-used sidecar is the Jacobian of the critical generators;
it distinguishes that situation from one nonreduced intersection. The latter
forces a zero Hessian and is impossible for a Keller realization.

The cases `D_A=0` and `D_A!=0` are exhaustive.
[THM-4171](THM-4171-row-a-inner-resultant-wall-planar-jacobian-exclusion.md)
proves the result on the entire `D_A=0` locus, including all four endpoint
strata. It remains to remove the projection-coordinate gates from `D_A!=0`.

## 2. Affine critical length without a discriminant

Assume henceforth

```text
eta*Delta*C*D_A != 0.                                  (5)
```

Put

```text
f=G_X/T,                         h=G_T.
```

The exact row-A resultant identity of
[THM-4147](THM-4147-generic-exact-weight-nine-planar-jacobian-monodromy-exclusion.md)
is

```text
Res_X(f,h)=T^42(6T+1)^2 Q_19(T),                        (6)

Q_19(0)=-12288 C^6,
[T^19]Q_19=-1458 C eta^4 D_A^2.                        (7)
```

Thus `Q_19` has exact degree nineteen and has no zero at `T=0`. Direct
differentiation of `(1)` gives

```text
deg_X(f,h)=(7,8),
LC_X(f)=LC_X(h)=8CT^7.                                 (8)
```

Consequently, over every `T!=0`, the projective `X`-closures of `f=0` and
`h=0` have no common point at `X=infinity`.

Let

```text
Z=V(f,h) intersect {T!=0}.
```

For every `t!=0`, the local resultant formula gives

```text
ord_(T=t) Res_X(f,h)
 =sum_(z in Z, T(z)=t) i_z(f,h),                       (9)
```

where `i_z(f,h)` is the local intersection length. Since `Q_19(0)!=0`, all
roots of `(6T+1)^2Q_19` lie on the `T!=0` chart. Summing `(9)` gives

```text
length(Z)=2+19=21.                                     (10)
```

This uses the resultant divisor with multiplicity. It requires neither
squarefreeness of `Q_19` nor disjointness from `T=-1/6`.

Now suppose `(1)` comes from a Keller realization. At an affine critical
point, the determinant-one Hessian congruence gives

```text
det Hess(G)!=0.                                        (11)
```

There is also the polynomial identity

```text
boxed: T det D(f,h)=det Hess(G)+f G_XT.                (12)
```

Indeed, `f_X=G_XX/T`, `f_T=(G_XT-f)/T`, `h_X=G_XT`, and
`h_T=G_TT`. At a point of `Z`, equations `(11)--(12)` imply

```text
det D(f,h)!=0.                                         (13)
```

Every local intersection in `(9)` is therefore simple. Thus `(10)` counts
twenty-one distinct affine critical points, even when several have the same
`T`-coordinate.

At `T=0` one has

```text
f(X,0)=-X,                 h(X,0)=-(X^2+6)/2.           (14)
```

The actual critical ideal is `(Tf,h)`, so it restores precisely

```text
T=0,                      X^2=-6,                       (15)
```

two distinct points with Hessian determinant `+6`. Hence the complete affine
critical scheme is reduced and has

```text
boxed: L=21+2=23.                                      (16)
```

The Hessian congruence also says that every critical point lies over one of
the two target nodes. If `r_0,r_1` are their inverse counts, then

```text
r_0+r_1=23.                                            (17)
```

If `Q_19(-1/6)=0`, the two known points
`T=-1/6, X^2=6` are simple with Hessian determinant `-6`. Any resultant
valuation beyond their two contributions is therefore contributed by
additional points in the same `T`-fibre; it is not lost length.

## 3. Complete packet and labelled carrier

Finiteness of the critical scheme and THM-3827 give geometric connectedness
of the generic source fibre. The row-A strict transform of THM-4157 is
unchanged throughout `(3)` and gives normalization genus at most ten and the
displayed boundary packet

```text
(7,7,4,2,2,2,1),                   defect=18.           (18)
```

THM-4103 and Riemann--Hurwitz over the elliptic target force genus at least
ten. Thus the genus is exactly ten, `(18)` is complete, and no affine or
unlisted boundary ramification remains.

The nonrational face is

```text
q-1/2=K W^2-eta W^3.                                  (19)
```

Since `eta!=0`, this is one prime separable cubic closed point over `C(q)`,
including when `K=0`. MISTAKE-509 is firewalled: its three geometric
conjugates cannot be assigned independently.

The finite-separable carrier theorem of THM-4147 gives the two exhaustive
responses:

```text
finite:  origin packet=(7,7,4,1),   (n,beta)=(19,3);
full:    origin packet=(7,7,4,2,2,2,1),  n=25.          (20)
```

This invocation retains the proved order of operations: first resolve,
normalize, and shrink to a proper smooth family with a finite relative map;
only then delete the origin, the complete cubic carrier, its full inverse
image, and the extra Hurwitz values. The resulting restriction is finite
etale. Disjoint parallel Milnor cores transport all `L` inverse points as
distinct fixed sheets, and connectedness supplies transitive generation.

## 4. The two contradictions

Let `X,Y` now denote the handle permutations. From `(17)` and fixed-sheet
transport,

```text
|supp X|+|supp Y|<=2n-23.                              (21)
```

For the finite response, `n=19` and the three carrier meridians are
transpositions. If at least one handle is nonidentity, their total generator
index is at most

```text
2n-L-1+beta=38-23-1+3=17<18=n-1.                      (22)
```

If both handles are identities, the total index is only `beta=3<18`.
Neither case can generate a transitive action on nineteen sheets.

For the full response, `n=25` and the handle permutations alone generate
transitively, so their supports cover all sheets. Equation `(21)` gives

```text
|supp X intersect supp Y|<=n-L=2.                     (23)
```

The commutator-overlap lemma of THM-4147 gives

```text
ind([X,Y])<=4.                                         (24)
```

But the origin meridian is `[X,Y]^-1`, and completeness of `(18)` gives it
index eighteen. Thus `18<=4`, a contradiction.

This excludes every point of `(5)`. THM-4171 excludes every point with
`D_A=0`. The two cases prove the theorem. **QED.**

## 5. Exact scope and replay

The theorem closes all of repeated-top row A under `(3)`. It deliberately
makes no assertion about:

1. `eta=0`;
2. `Delta=0`, where the common Newton polygon contracts;
3. `C=0`, namely rows B, C, and D of THM-4157, not asserted here but closed
   in full by THM-4176;
4. entry into the reduced seam or another reduced cell;
5. exact residual weight at least ten;
6. failure of the inherited finite-separable transport package;
7. `JC(2)` or `DC(2)`.

The source projection `Res_s(Acrit,Ccrit)=p^8R_19` is not used. Hence its
endpoints and discriminant are not hypotheses of this theorem. Neither audit
calculates `disc(Q_19)`, `disc(R_19)`, or `Q_19(-1/6)`.

Run the exact bridge audits with

```text
python3 -B 04-computation/jc23_repeated_top_row_a_complete_exclusion_thm4173.py
python3 -B -O 04-computation/jc23_repeated_top_row_a_complete_exclusion_thm4173.py
PYTHONHASHSEED=4173 python3 -B \
  04-computation/jc23_repeated_top_row_a_complete_exclusion_thm4173.py

python3 -B \
  04-computation/jc23_repeated_top_row_a_complete_exclusion_thm4173_independent_bridge_audit.py
python3 -B -O \
  04-computation/jc23_repeated_top_row_a_complete_exclusion_thm4173_independent_bridge_audit.py
PYTHONHASHSEED=4173 python3 -B \
  04-computation/jc23_repeated_top_row_a_complete_exclusion_thm4173_independent_bridge_audit.py
```

The independent audit uses the same projected divisor `T^2` for two reduced
points with one projection coordinate and for one doubled point; only the
second has zero critical Jacobian. This is the exact hostile separating the
old projection observer from the repaired local-intersection argument.
