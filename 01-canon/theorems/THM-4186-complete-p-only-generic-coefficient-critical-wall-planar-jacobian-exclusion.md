---
id: THM-4186
title: "Complete P-only generic-coefficient critical-wall planar Jacobian exclusion"
status: >
  PROVED RELATIVE TO THM-3827/3992/3997/4007/4103/4120/4122/4130/
  4147/4176/4183 + VERIFIED-EXACT + INDEPENDENTLY NORMALIZED-AUDITED. On
  the live exact-weight-nine b=d=0 reduced (2,3) seam, the complete P-only
  generic coefficient chamber zeta=0, eta*Delta*K*Theta!=0 contains no
  nonautomorphic planar Keller pair, with Phi arbitrary and without any
  critical-resultant discriminant, Q_20(-1/6), or projected-coordinate
  separation hypothesis. A new symbolic (s,p) resultant and Morse bridge
  keep affine critical length 24 even on every projected collision wall.
  Exact hostiles prove both the repeated-projection and universal-fibre walls
  are nonempty. THM-4183 separately closes Delta=0. The P-only coefficient
  walls Delta!=0 with K*Theta=0, mixed B walls, other cells, entry, M>=10,
  JC(2), and DC(2) remain OPEN.
source: jc-next-wall-20260826
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
  - THM-4176-complete-repeated-top-wall-planar-jacobian-exclusion
  - THM-4183-p-only-delta-zero-planar-jacobian-exclusion
related:
  - THM-4155-generic-y-only-delta-zero-weight-nine-planar-jacobian-exclusion
  - THM-4159-inner-resultant-wall-planar-jacobian-exclusion
  - THM-4161-y-only-double-top-root-planar-jacobian-exclusion
  - THM-4164-y-only-triple-top-root-planar-jacobian-exclusion
  - THM-4165-y-only-inner-top-triple-intersection-planar-jacobian-exclusion
script: 04-computation/jc23_p_only_generic_critical_wall_complete_exclusion_thm4186.py
output: 05-knowledge/results/jc23_p_only_generic_critical_wall_complete_exclusion_thm4186.out
independent_audit_script: 04-computation/jc23_p_only_generic_critical_wall_complete_exclusion_independent_audit_thm4186.py
independent_audit_output: 05-knowledge/results/jc23_p_only_generic_critical_wall_complete_exclusion_independent_audit_thm4186.out
script_sha256: 76b5bdb6c58482624f06ac4dded6e3f5cb0d4a23b33981e486b9e63970983d9f
output_sha256: 5ea5464934f7dc7652d9144fd87b7a797370cc8bb696863836e3c00875be5a83
independent_audit_script_sha256: 6103862647fc5e206d185539f43fd93a5e4746865a09977545ea393bef9de6c8
independent_audit_output_sha256: 0e0e08eca0fedcf4c415d7ffc9141716a08361082c35159eeb41a8532a4002f1
hash_basis: raw LF bytes
primary_audit: >
  PASS. A standalone symbolic (s,p) calculation reconstructs the complete
  source, polynomial critical pair, incompatible leading rows, exact
  Hessian-bridge ideal membership, and the full four-parameter resultant
  p^4 R_20 with both endpoint formulas. It independently resolves the
  p=0/t=0 loss ledger and proves that the two hostile normalized collisions
  separate into squarefree p-residual points. Normal, optimized, and fixed-
  hash executions byte-match the frozen output.
independent_audit: >
  ACCEPT. A separate normalized (X,T) implementation imports no primary
  polynomial. It reconstructs exact-M=9 source completeness, both leading
  rows, the normalized Hessian bridge, both universal Morse pairs, and three
  exact resultants. It finds gcd(Q_20,Q_20')=T-1 on the projected-collision
  hostile and Q_20(-1/6)=0 with fibre gcd (X-1)(X^2-6) on the universal-
  fibre hostile. Normal, optimized, and fixed-hash executions byte-match.
---

# THM-4186 -- complete P-only generic-coefficient critical-wall exclusion

**PROVED RELATIVE TO THM-3827/3992/3997/4007/4103/4120/4122/4130/
4147/4176/4183 + VERIFIED-EXACT + INDEPENDENTLY NORMALIZED-AUDITED; JC(2)
AND DC(2) REMAIN OPEN.**

## 1. Theorem and inheritance pass

Work over `C` in the live `b=d=0` reduced `(2,3)` seam at exact residual
weight nine. Put

```text
P=T+X^2T^2,                    Y=XTP,
K=2848/45-(7/6)Delta.                                    (1)
```

Use the complete P-only source

```text
G=-X^2T/2-3P+(8/3)P^2-(1376/135)P^3
  +K Y^2+Phi P^2Y+Delta P^4+Theta P Y^2+eta P^3Y.       (2)
```

> **Theorem.** For every coefficient tuple satisfying
>
> ```text
> zeta=0,                    eta*Delta*K*Theta!=0,       (3)
> ```
>
> polynomial `(2)` is not the normalized exact-weight-nine source of a
> nonautomorphic planar Keller pair in the inherited reduced seam.
>
> The coefficient `Phi` is unrestricted. No discriminant of a source or
> normalized critical eliminant, no condition at `T=-1/6`, and no
> projected-coordinate separation hypothesis is imposed.

THM-4147 proved `(3)` only on its critical-open locus. The closest proved
mechanism is THM-4176's Morse-resultant divisor lemma, but the primary proof
below derives its source-chart version and computes a new complete symbolic
resultant. The canonical hostile is a repeated eliminant root caused by two
distinct reduced source points sharing one projection coordinate. The
corrected near miss is MISTAKE-421: etaleness of a finite source scheme does
not make a selected coordinate primitive or its eliminant squarefree.
MISTAKE-445/450 separately forbid reading a vanished fixed top slot as an
actual leading coefficient. MISTAKE-509 requires the two geometric branches
of the quadratic boundary point to respond together.

The least-used sidecar is the independent source coordinate

```text
p=T+(XT)^2.                                               (4)
```

It separates the first exact normalized discriminant hostile below and
turns the second into one ordinary nonzero-`p` root beside the collapsed
universal pair. This distinguishes source reducedness from projection
squarefreeness rather than assuming the distinction away.

## 2. Exact source and the Y-boundary scope correction

Put

```text
s=XT,              p=T+s^2,              y=sp,
t=p-s^2=T.                                                (5)
```

THM-3992/3997 and THM-4007 give the complete normalized source

```text
G=-s^2/(2t)+H(p,y),

H=-3p+(8/3)p^2-(1376/135)p^3+K y^2+Phi p^2y
  +Delta p^4+Theta p y^2+eta p^3y+zeta y^3.             (6)
```

For `wt(p)=2,wt(y)=3`, enumeration through weight nine, with the forced
deletions `y,py`, gives

```text
p, p^2, p^3, y^2, p^2y, p^4, py^2, p^3y, y^3.          (7)
```

The only weight-nine terms are `p^3y,y^3`; consequently

```text
exact M=9  iff  (eta,zeta)!=(0,0).                       (8)
```

This corrects the apparent “`zeta=0` boundary of Y-only” residual. Y-only
means `eta=0,zeta!=0`. Setting `zeta=0` there forces `eta=zeta=0` and exits
to `M<=8`; it is not a live exact-`M=9` wall and supplies no new exclusion
theorem. If instead `zeta=0,eta!=0`, the replacing top coefficient is
`eta p^3y`, and the source is precisely the P-only row `(2)--(3)`.

The correction is a filtration statement only. It neither constructs a
transverse Keller deformation from the lower-weight fibre nor enlarges the
scope of the exact-`M<=8` canon.

## 3. Complete source critical divisor

In `(s,p)` coordinates write `(6)` on the P row as

```text
H=-3p+(8/3)p^2-(1376/135)p^3+K s^2p^2+Phi sp^3
  +Delta p^4+Theta s^2p^3+eta sp^4.                    (9)
```

Define the polynomial critical pair

```text
A=t^2G_s/p,                       C=2t^2G_p.            (10)
```

Exact differentiation proves that both are polynomials and that

```text
t^2G_s=pA,                         2t^2G_p=C.           (11)
```

Their `s`-degrees and leading rows are

```text
(deg_s A,deg_s C)=(5,6),
LC_s(A)=2p(K+Theta p),
LC_s(C)=2p(2K+3Theta p),
3LC_s(A)-LC_s(C)=2pK.                                  (12)
```

Thus on `p!=0,K!=0` the two leading rows cannot vanish together. No common
critical point is lost at `s=infinity`, including on every projected-root
collision locus.

The primary exact calculation over
`Q[Delta,Phi,Theta,eta]` gives the complete resultant

```text
Res_s(A,C)=p^4 R_20(p),

R_20(0)=-31104K^2
       =-(96/25)(105Delta-5696)^2,

[p^20]R_20=-65610 eta^6 Theta.                         (13)
```

Both endpoints are units on `(3)`. Hence `R_20` has actual degree twenty,
no residual root has `p=0`, and no fixed-slot/actual-degree ambiguity
survives. Direct specialization to `t=0`, equivalently `p=s^2`, gives

```text
A=-s,                              C=s^2.               (14)
```

Therefore no residual point with `p!=0` lies on the singular coordinate row
`t=0`. Every root of `R_20` comes from the valid chart `p*t!=0`.

## 4. Full-source reducedness, not eliminant squarefreeness

Let

```text
J_AC=det D_(s,p)(A,C),             H_G=det Hess_(s,p)(G).
```

Differentiating `(11)` gives the exact ideal-membership identity

```text
t(2t^4H_G-pJ_AC)

= A(-4Cs+4C_p p s+2C_s p-C_s t)
  +C(-4A_p p s-2A_s p).                                 (15)
```

On a common zero of `A,C` with `p*t!=0`, this becomes

```text
p J_AC=2t^4 H_G.                                        (16)
```

The coordinate map `(X,T)->(s,p)` has Jacobian `t`. At a critical point its
Hessians are therefore related by congruence. Under a hypothetical Keller
realization `G=E(A_0,C_0)`, the unit source Jacobian and the two Morse target
nodes give

```text
det Hess_(X,T)(G)=det Hess(E)!=0.                        (17)
```

Equations `(16)--(17)` make every open-chart zero of `(A,C)` reduced. Since
`(12)` excludes `s=infinity`, the valuation of the `p`-resultant at any
finite `p_0!=0` is the sum of the local intersection lengths above `p_0`.
The sum of those valuations is the actual degree twenty in `(13)`. Thus
there are exactly twenty distinct open-chart affine critical points, even
when several share one `p`-coordinate.

This is the necessary/sufficient distinction:

```text
Keller--Morse  =>  full critical scheme reduced;
full scheme reduced  !=>  a chosen eliminant squarefree.              (18)
```

No converse is used.

## 5. The four collapsed universal points and exact length

The factor `p^4` in `(13)` is a coordinate artifact, not four points that
may be counted from the source resultant. Return to `(X,T)` and put

```text
f=G_X/T,                              h=G_T.             (19)
```

At `T=0`, exact specialization gives

```text
f=-X,              h=-(X^2+6)/2,
X^2=-6,            G=0,              det Hess(G)=+6.    (20)
```

These are two reduced points. On `p=P=0` with `T!=0`, one has
`T=-1/X^2`, and substitution into `(19)` gives

```text
f=-(X^2-6)/X,                   h=-(X^2-6)/2.           (21)
```

Thus the only remaining `p=0` critical points are

```text
T=-1/6,             X^2=6,
G=1/2,              det Hess(G)=-6.                    (22)
```

There are no others. Combining `(13),(20),(22)` yields the exact affine
critical length

```text
L=20+2+2=24.                                           (23)
```

The normalized identity

```text
T det D(f,h)=det Hess(G)+f G_XT                         (24)
```

independently verifies the same reducedness mechanism on `T!=0`. In
particular, if a residual factor meets `(6T+1)`, its extra resultant
valuation counts additional distinct points in that fibre. It does not
erase either universal point or lower the length.

Every affine critical point of a hypothetical Keller realization lies over
one of the two target nodes: `grad G=D(A_0,C_0)^t grad E`, with the matrix
invertible. If `r_i` is the number above the node of value `i/2`, then

```text
r_0+r_1=24.                                             (25)
```

## 6. Complete packet and quadratic carrier

The coefficient assumptions `(3)` are exactly the units required by
THM-4147's P-row Newton calculation. The coefficient `Phi` lies on no
boundary face. The lower Newton polygon, Pick ledger, and complete packet are

```text
(0,1),(2,0),(4,2),(4,3),(3,4),(1,5),(0,5),
(2Area,B,g)=(29,11,10),
packet=(8,5,4,3,2,2,1),                defect=18.       (26)
```

Equation `(13)` makes the affine critical scheme finite for every tuple in
`(3)`. THM-3827's closed-polynomial factor theorem therefore supplies
geometric connectedness. Pick gives normalization genus at most ten, while
THM-4103 and Riemann--Hurwitz over the elliptic target use the displayed
defect to give genus at least ten. Hence the genus and packet in `(26)` are
complete on the whole coefficient chamber; no critical-open hypothesis is
needed.

The only nonrational place is

```text
q-1/2=K W^2.                                            (27)
```

Since `K!=0`, `(27)` is one separable irreducible quadratic closed point over
`C(q)`: `q-1/2` has odd valuation at `q=1/2`. THM-4120 gives
`E_q(C(q))={O}`. If the quadratic place has finite horizontal image with
residue field `M`, then

```text
C(q) subset M subset C(q)(W).                           (28)
```

Prime degree and the absence of a finite rational target point force the
full quadratic residue field. MISTAKE-509's forbidden singleton allocation
never appears: the two geometric branches either both map to the origin or
become two conjugate carrier points, each with a transposition meridian.
Thus the exhaustive responses remain

```text
full:       n=25;
finite:     n=21,                     beta=2.           (29)
```

THM-4147's proper-smooth/finite-etale construction applies because Sections
3--6 have now established all its source-side hypotheses: finite reduced
critical scheme under Keller, geometric connectedness, complete packet, and
one separable prime carrier. Parallel Milnor-core transport gives handle
permutations `X_0,X_1` satisfying

```text
#Fix(X_i)>=r_i,
|supp X_0|+|supp X_1|<=2n-L.                           (30)
```

## 7. Both strict monodromy contradictions

For the full response, `X_0,X_1` generate transitively and their supports
cover all `25` sheets. Hence

```text
|supp X_0 intersect supp X_1|<=25-24=1.                (31)
```

THM-4147's commutator-overlap lemma gives

```text
ind([X_0,X_1])<=2.                                     (32)
```

The origin meridian has index equal to the complete packet defect `18`, a
contradiction.

For the finite response, the handles and two carrier transpositions must
generate transitively on `21` sheets. If at least one handle is nonidentity,
their total merger capacity is at most

```text
2n-L-1+beta=42-24-1+2=19<20=n-1.                      (33)
```

If both handles are identities, the capacity is only `beta=2<20`.
Transitivity is impossible in either case. This proves the theorem. **QED.**

## 8. Two exact hostile walls and failure anatomy

The new scope is nonvacuous in two different ways.

### 8.1 A reduced source with a repeated projected root

Take

```text
Delta=1271/180,          Phi=1733/7560,
Theta=-206281/7560,      eta=-1733/7560,
K=11891/216.                                               (34)
```

Every theorem coefficient is nonzero. In normalized coordinates,

```text
gcd(Q_20,Q_20')=T-1,
gcd_X(f(X,1),h(X,1))=X(X-1).                            (35)
```

Both points are Morse. In the independent source chart they are

```text
(s,p)=(0,1), (1,2),                  t=p-s^2=1,         (36)
```

and the source residual `R_20(p)` is squarefree with the two distinct roots
`p=1,2`. Thus the first failed implication is exactly

```text
two reduced points share T  !=>  the T-eliminant is squarefree.         (37)
```

Their values are `G=-1871/540` and `G=33724/945`, so this coefficient tuple
is only a wall witness, not a Keller candidate: neither point lies over a
target node.

### 8.2 An extra point in the universal projected fibre

Take

```text
Delta=1,                 Phi=-1176023/2700,
Theta=32981/450,         eta=1,
K=5591/90.                                                (38)
```

Then

```text
Q_20(-1/6)=0,
gcd_X(f(X,-1/6),h(X,-1/6))=(X-1)(X^2-6).               (39)
```

The new point `X=1` is Morse, and `Q_20` itself is squarefree there. In the
source chart it becomes

```text
(s,p)=(-1/6,-5/36),                   t=-1/6,           (40)
```

which is separated from the two universal `p=0` points; the source
`R_20(p)` is again squarefree. Its value

```text
G(1,-1/6)=12468079/30233088                           (41)
```

is neither target-node value. The mechanism is an extra reduced point in an
already occupied projected fibre, not a collision with or loss of the
universal pair.

These hostiles explain why deleting the discriminant and `Q_20(-1/6)` was a
coordinate convenience in THM-4147 rather than a geometric necessity.

## 9. Exact scope and surviving walls

Combining this theorem with THM-4183 gives the exact proved P-only union

```text
zeta=0, eta!=0, and

  Delta=0
  or
  Delta*K*Theta!=0.                                    (42)
```

This is not the whole P-only coefficient locus. The remaining P walls are
precisely

```text
Delta!=0,                     K*Theta=0.                (43)
```

On `Theta=0` the P polygon blows down; on `K=0` the quadratic carrier
equation `(27)` collapses. Neither response ledger may be transported without
a new boundary normalization. The theorem also does not close the mixed B
chamber or its coefficient/critical walls, another reduced cell, entry into
this seam, residual weight at least ten, a transport setup outside the
inherited proper-smooth contract, `JC(2)`, or `DC(2)`.

The Y-only `zeta=0` phrase is superseded only as an exact-`M=9` routing
description, not as a new theorem about lower weight. The honest next tests
are the two geometrically different coefficient contractions in `(43)`.

## 10. Exact artifacts and replay

Primary symbolic source-chart certificate:

```text
04-computation/jc23_p_only_generic_critical_wall_complete_exclusion_thm4186.py
sha256 76b5bdb6c58482624f06ac4dded6e3f5cb0d4a23b33981e486b9e63970983d9f

05-knowledge/results/jc23_p_only_generic_critical_wall_complete_exclusion_thm4186.out
sha256 5ea5464934f7dc7652d9144fd87b7a797370cc8bb696863836e3c00875be5a83
```

Independent normalized-coordinate certificate:

```text
04-computation/jc23_p_only_generic_critical_wall_complete_exclusion_independent_audit_thm4186.py
sha256 6103862647fc5e206d185539f43fd93a5e4746865a09977545ea393bef9de6c8

05-knowledge/results/jc23_p_only_generic_critical_wall_complete_exclusion_independent_audit_thm4186.out
sha256 0e0e08eca0fedcf4c415d7ffc9141716a08361082c35159eeb41a8532a4002f1
```

Replay both paths in normal, optimized, and fixed-hash modes:

```bash
python3 -B \
  04-computation/jc23_p_only_generic_critical_wall_complete_exclusion_thm4186.py \
  | diff -u \
      05-knowledge/results/jc23_p_only_generic_critical_wall_complete_exclusion_thm4186.out -
python3 -B -O \
  04-computation/jc23_p_only_generic_critical_wall_complete_exclusion_thm4186.py \
  | diff -u \
      05-knowledge/results/jc23_p_only_generic_critical_wall_complete_exclusion_thm4186.out -
PYTHONHASHSEED=4186 python3 -B \
  04-computation/jc23_p_only_generic_critical_wall_complete_exclusion_thm4186.py \
  | diff -u \
      05-knowledge/results/jc23_p_only_generic_critical_wall_complete_exclusion_thm4186.out -

python3 -B \
  04-computation/jc23_p_only_generic_critical_wall_complete_exclusion_independent_audit_thm4186.py \
  | diff -u \
      05-knowledge/results/jc23_p_only_generic_critical_wall_complete_exclusion_independent_audit_thm4186.out -
python3 -B -O \
  04-computation/jc23_p_only_generic_critical_wall_complete_exclusion_independent_audit_thm4186.py \
  | diff -u \
      05-knowledge/results/jc23_p_only_generic_critical_wall_complete_exclusion_independent_audit_thm4186.out -
PYTHONHASHSEED=4186 python3 -B \
  04-computation/jc23_p_only_generic_critical_wall_complete_exclusion_independent_audit_thm4186.py \
  | diff -u \
      05-knowledge/results/jc23_p_only_generic_critical_wall_complete_exclusion_independent_audit_thm4186.out -
```

All three primary streams print `checks=48` and
`THM4186_SOURCE_CHART_EXACT_ACCEPT`. All three independent streams print
`checks=57` and `THM4186_NORMALIZED_INDEPENDENT_ACCEPT`.
