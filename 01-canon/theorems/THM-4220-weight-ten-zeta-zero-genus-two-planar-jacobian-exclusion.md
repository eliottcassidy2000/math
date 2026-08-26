---
id: THM-4220
title: "Weight-ten zeta-zero genus-two planar Jacobian exclusion"
status: >
  PROVED RELATIVE TO THM-3992/3997/4012/4045/4218 + VERIFIED-EXACT
  + INDEPENDENTLY AUDITED.
  In the inherited b=d=0 reduced (2,3) exact-M=10 seam, the zeta=0
  wall with upsilon*xi*(upsilon+xi) nonzero contains no nonautomorphic
  planar Keller pair when D_V=Theta^2-4Kxi is nonzero or when
  K=Theta=0. The split-conic collision D_V=0,K!=0, other cells, seam
  entry, M>=11, JC(2), and DC(2) remain OPEN.
source: codex-planar-jacobian-weight-ten-wall-session-20260826
depends_on:
  - THM-3992-reduced-two-three-cusp-jet-repair-and-first-node-residual
  - THM-3997-reduced-two-three-hasse-repair-and-zero-residual-no-go
  - THM-4012-weighted-leading-face-good-elliptic-factor-observer
  - THM-4045-live-two-three-max-seven-hidden-elliptic-tail-no-go
  - THM-4218-exact-weight-ten-hidden-elliptic-tail-degree-three-planar-jacobian-exclusion
related:
  - THM-4217-complete-mixed-off-antidiagonal-delta-zero-planar-jacobian-exclusion
external: >
  Tim Dokchitser, "Models of curves over DVRs," arXiv:1807.00025v2,
  Definitions 3.7, 3.9, and 3.12 and Theorem 3.14, only through the audited
  use in THM-4045. The imported theorem supplies the toric face/edge model
  and rational chains; it supplies none of the face arithmetic, genus-two
  simplicity, or planar-Jacobian conclusion below.
script: 04-computation/jc23_weight10_zeta_zero_genus2_exclusion_thm4220.py
output: 05-knowledge/results/jc23_weight10_zeta_zero_genus2_exclusion_thm4220.out
script_sha256: 18b6326978cc85aad3e66da0c5abd5e61e3b345991ccfcc108ca3a62c949b6cf
output_sha256: d858845148074b4e6742bf5824cc90c7a2ad986fe179b2f63177010ba2feef3b
independent_audit_script: 04-computation/jc23_weight10_zeta_zero_genus2_exclusion_thm4220_independent_audit.py
independent_audit_output: 05-knowledge/results/jc23_weight10_zeta_zero_genus2_exclusion_thm4220_independent_audit.out
independent_audit_script_sha256: 086f1c68595fa908b03e14a80d94cf927d29783ced64a315e6968d7296be4272
independent_audit_output_sha256: 0e9b59e310d6daa22fcc34cfb72752cab3c4cc3bab05054217f049fdd4ca547e
---

# THM-4220 -- weight-ten zeta-zero genus-two exclusion

**PROVED RELATIVE TO THM-3992/3997/4012/4045/4218 + VERIFIED-EXACT
+ INDEPENDENTLY AUDITED; JC(2) REMAINS OPEN.**

## 1. Statement and scope

Work over `C` in the inherited `b=d=0` reduced `(2,3)` seam. In the notation
of THM-4218 put

```text
upsilon=[p^5]H,       xi=[p^2y^2]H,       zeta=[y^3]H,
D_V=Theta^2-4Kxi.                                      (1)
```

> **Theorem.** On the complete exact-weight-ten locus
>
> ```text
> zeta=0,             upsilon*xi*(upsilon+xi)!=0,       (2)
> ```
>
> either additional condition
>
> ```text
> D_V!=0                         or       K=Theta=0      (3)
> ```
>
> excludes a nonautomorphic planar Keller pair.

Every lower coefficient is arbitrary subject only to `(2)--(3)` and the
inherited normal-form relation. Because `xi` and `upsilon` are nonzero, this
is still exact `M=10`; setting `zeta=0` is not a filtration exit. The
complement of `(3)` inside `(2)` is exactly

```text
D_V=0,                         K!=0.                    (4)
```

That split-conic collision remains **OPEN**. The theorem proves neither seam
entry nor `JC(2)`.

The inheritance pass is:

- closest mechanism: THM-4218's complete lower Newton model and simple
  genus-two main component;
- canonical hostile: a rational-looking singular face can acquire a hidden
  positive-genus tail (the warning in THM-4012);
- corrected near miss: deleting the cubic owner does not always leave the
  main polygon--the `y^2` row creates a new rational side face;
- least-used sidecar: the entire positive-genus inventory, rather than a
  carrier-response congruence.

The replacement obstruction is total `Hom`-vanishing: after specialization
there is no component capable of mapping nonconstantly to the good elliptic
target.

## 2. Complete source and the two possible lower planes

Use source coordinates `s=XT`, `p=T+s^2`, `y=sp`. The complete polynomial is

```text
H=-3p+(8/3)p^2-(1376/135)p^3+K y^2+Phi p^2y
  +Delta p^4+Theta py^2+eta p^3y
  +upsilon p^5+xi p^2y^2,

K=2848/45-(7/6)Delta.                                  (5)
```

There is no ellipsis. The source pencil is

```text
F_Q=(s^2-p)(1-QH)-Qs^2/2,              Q=q^-1.         (6)
```

On its torus `p=s^2` is impossible, and `t=p-s^2`, `X=s/t` recovers the
actual source function field.

A monomial `p^i y^j` gives height-one endpoints

```text
(j+2,i+j,1),                    (j,i+j+1,1).            (7)
```

For `(r,k)` the possible lower height functions are

```text
M: nu_M=(r+2k-2)/10,            V: nu_V=(r-2)/2.       (8)
```

The gap above `M` is `(10-(2i+3j))/10` for both endpoints. The gaps above
`V` are `1-j/2` and `2-j/2`. With `zeta=0` one has `j<=2`, so all gaps are
nonnegative; equality on `V` occurs exactly at the first endpoint of a
`y^2` term. Thus `M` is always present, `V` is present exactly when
`(K,Theta)!=(0,0)`, and if `K=Theta=0` the sole `xi` equality is only a line
with the base vertex. Unique outer owners make this invariant under all
lower-coefficient cancellations.

## 3. Exact polygons, packets, and carriers

Put

```text
A=(0,1), B=(2,0), U=(4,2), W=(4,3), D=(4,4), E=(0,6).
                                                                  (9)
```

Pick's theorem and the intrinsic edge-index formula give:

| stratum | global polygon | `2Area` | boundary | genus | infinity packet | sum |
|---|---|---:|---:|---:|---|---:|
| `K!=0,D_V!=0` | `A,B,U,D,E` | 34 | 12 | 12 | `(9,9,3,3,2,2,1)` | 29 |
| `K=0,Theta!=0` | `A,B,W,D,E` | 32 | 10 | 12 | `(9,9,5,3,1)` | 27 |
| `K=Theta=0` | `A,B,D,E` | 30 | 10 | 11 | `(9,9,3,3,1)` | 25 |

For an oriented edge with primitive inward normal `(a,b)` and support
constant `c`, the displayed index is `a+b-c`; the vertical `EA` edge is the
affine divisor `s=0` and is omitted. Every row satisfies
`sum(e_P-1)=2g-2`.

The moving outer restriction is, after a monomial change,

```text
K!=0:                 q-1/2=K W^2;
K=0,Theta!=0:         q-1/2=Theta W;
K=Theta=0:            q-1/2=xi W^2.                   (10)
```

The first and third rows have a prime separable quadratic carrier over
`C(q)`; the middle has a rational section. The proof does not use a carrier
response.

## 4. Face normalizations and genus ledger

The main face is unchanged from THM-4218:

```text
g_M=(S^2-P)(1-upsilon P^5-xi S^2P^4)=R*C.             (11)
```

The rational branch `R` and the other branch meet at

```text
P=S^2,                         S^10=(upsilon+xi)^-1.    (12)
```

Their gradient determinant is `-10S^9(upsilon+xi)`, so these are exactly ten
transverse intersections. On `C`, put `Y=SP^2`; its normalization is

```text
xi Y^2=1-upsilon P^5,                                 (13)
```

a smooth genus-two curve.

When the side face exists, its torus equation is

```text
g_V=S^2(1-S^2P^2(K+Theta P+xi P^2)).                  (14)
```

Removing the monomial and putting `T=SP`, `Y=T^-1` gives

```text
Y^2=K+Theta P+xi P^2.                                 (15)
```

Under `D_V!=0`, this is a smooth irreducible conic, hence rational. The
common internal edge `BD` has scheme

```text
1-xi Z^2.                                              (16)
```

Its two reduced roots give two labelled paths from `V` to `C`.

The other special edge schemes, up to units and monomials, are

```text
K!=0:
 AB 1-Z,  BU 1-KZ^2,  UD K+Theta Z+xi Z^2,
 DE (1-Z)(upsilon+xi Z),  EA 1-upsilon Z^5;

K=0,Theta!=0:
 AB 1-Z,  BW 1-Theta Z,  WD Theta+xi Z,
 DE (1-Z)(upsilon+xi Z),  EA 1-upsilon Z^5.           (17)
```

Together with `(16)`, every scheme is reduced under `(2)` and `D_V!=0`, and
no root meets a corner. In the collapsed `K=Theta=0` stratum, `BD` is an
outer edge and the four schemes `AB,BD,DE,EA` are reduced under `(2)`.

The complete normalized core graph is therefore

```text
D_V!=0:       vertices R,C,V;  10 R--C paths, 2 C--V paths;
K=Theta=0:    vertices R,C;    10 R--C paths.           (18)
```

The internal `BD` slope pair is exactly

```text
0 > -6.                                                (19)
```

Dokchitser Definition 3.12 inserts the five intermediate integral slopes
for each of the two roots of `(16)`: ten rational multiplicity-one curves in
all. Replacing an edge by a path does not change the graph cycle rank. Thus

```text
D_V!=0:       b_1=12-3+1=10,     g=2+10=12;
K=Theta=0:    b_1=10-2+1=9,      g=2+9=11.             (20)
```

These are exactly the Pick genera in Section 3. All other face, toric, and
resolution components are rational. In particular, `C` is the only
positive-genus component; the full lower model has no hidden tail.

## 5. Multiplicity-one regular model

Take the common base change

```text
Q=sigma^30.                                             (21)
```

The integral face heights are

```text
M: 3r+6k-6,                       V: 15r-30.           (22)
```

Their primitive three-dimensional normals have last coordinate one, so the
face components have multiplicity one. The edge denominators in `(16)--(19)`
are one; all inserted chains have multiplicity one.

At each point of `(12)`, use

```text
s=sigma^-3 S,             p=sigma^-6 P,
H_M=sigma^30 H(sigma^-6P,sigma^-9SP),
U=S^2-P,                  V_0=(1-H_M)/S^2.             (23)
```

The exact scaled source equation is

```text
U V_0=sigma^30/2.                                      (24)
```

Each crossing is therefore an `A_29` smoothing. Its regular resolution adds
only multiplicity-one rational curves.

The side chart is equally explicit:

```text
s=sigma^-15 S,          p=P,
H_V=sigma^30 H(P,sigma^-15SP).                         (25)
```

After multiplying `(6)` by `sigma^30`, one obtains exactly

```text
(S^2-sigma^30P)(1-H_V)-sigma^30S^2/2=0.               (26)
```

Its special torus fibre is `(14)`. The face smoothness `(15)`, reduced edge
schemes `(16)--(17)`, and local charts `(24),(26)` meet the complete
face/edge regularity gate inherited through THM-4045. The generic moving
edges `(10)` are separable over `C(q)`; the affine `s=0` restriction is
`H(p,0)-q`, separable because its derivative cannot have a constant root
whose image is the transcendental `q`. Therefore the generic toric
completion is smooth and, by the function-field inverse after `(6)`, is the
actual source completion.

## 6. No component can carry the cover degree

The genus-two curve `(13)` is exactly the component in THM-4218. Its
order-five automorphism acts on `dP/Y,PdP/Y` by `rho,rho^2`, so
`Q(zeta_5)` embeds in `End^0(J(C))`. THM-4218's endomorphism-algebra audit
proves `J(C)` simple. Hence

```text
Hom(J(C),E_0)=0                                        (27)
```

for every elliptic curve `E_0`, and every morphism from `C` to an elliptic
curve is constant.

Since `q=sigma^-30`, scale the inherited target by

```text
A=sigma^-10 X,                    C_target=sigma^-15 Y. (28)
```

The exact target model becomes

```text
Y^2=X^3+1-(3a^2/4)sigma^20X-(a^3/4)sigma^30,           (29)
```

with smooth special fibre `E_0:y^2=x^3+1`.

A hypothetical Keller pair gives a finite nonconstant generic morphism from
the smooth source completion to this target. Resolve its rational extension
over the regular source model. The branch `R`, side conic `V`, every toric or
resolution chain, and every exceptional component are rational and hence
map constantly to `E_0` by Riemann--Hurwitz. Equation `(27)` makes `C`
constant as well.

For a relative target line bundle `L` of degree one, proper flat degree
conservation now gives

```text
deg(phi_generic)=sum_i m_i deg(phi_i^*L)=0.            (30)
```

This contradicts generic finiteness and nonconstancy. Both strata in `(3)`
are excluded.

## 7. Exact failure boundary

On `(4)`, necessarily `K,Theta,xi` are nonzero and

```text
K+Theta P+xi P^2=xi(P+Theta/(2xi))^2.                  (31)
```

The side conic splits and the `UD` edge scheme has a double root. The
one-step regular model above is no longer valid. A strict-transform audit
must retain the first lower term normal to that edge and prove that no new
positive-genus component carries degree. This theorem makes no claim there.

The typed connection is

```text
source:       complete zeta-zero exact-M=10 lower model;
target:       degree of a specialized cover of a good elliptic curve;
map:          Q=sigma^30 regular specialization;
preserved:    multiplicities, all component genera, and generic degree;
destroyed:    rational-chain parametrizations and off-face coefficients;
sidecars:     genus-two simplicity and reduced edge schemes;
decisive test: every special component has map degree zero;
hostile:      D_V=0,K!=0 split conic / double outer root.            (32)
```

## 8. Verification

The primary certificate checks the complete support, analytic plane gaps, all three
polygons and packets, Pick/defect ledgers, carrier discriminants, face
factorizations, conic discriminant, special edges, the `0>-6` chain ledger,
the exact `A_29` and side charts, component and graph genera, base-change
heights, target scaling, degree ledger, and hostile split factorization.
The clean-room audit independently enumerates the lower-support and collision
deletions, reconstructs the regular-model slopes, and rechecks component
completeness and the degree-zero obstruction.

Replay with

```bash
python3 -B 04-computation/jc23_weight10_zeta_zero_genus2_exclusion_thm4220.py
python3 -B -O 04-computation/jc23_weight10_zeta_zero_genus2_exclusion_thm4220.py
PYTHONHASHSEED=4220 python3 -B 04-computation/jc23_weight10_zeta_zero_genus2_exclusion_thm4220.py

python3 -B 04-computation/jc23_weight10_zeta_zero_genus2_exclusion_thm4220_independent_audit.py
python3 -B -O 04-computation/jc23_weight10_zeta_zero_genus2_exclusion_thm4220_independent_audit.py
PYTHONHASHSEED=4220 python3 -B 04-computation/jc23_weight10_zeta_zero_genus2_exclusion_thm4220_independent_audit.py
```

Compare stdout with the frozen output. **QED.**
