---
id: THM-4053
title: "Live reduced JC(2) maximum-weight-eight trichotomy and Eisenstein survivor"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED on the live b=d=0 reduced
  (2,3) seam at exact total residual weight eight; JC(2) remains OPEN. The
  complete lower model has an exact trichotomy. Off its repeated-edge wall,
  the p^4-only stratum has either one j=1728 elliptic component or only
  rational components and is impossible. Off delta+theta=0, the two-term top
  stratum has only the Bolza abelian factor Jac~E_8000^2 and is impossible.
  The sole nonresonant survivor is p*y^2-only: its unique positive-genus
  component is the target-compatible j=0 curve, and the generic fibre degree
  must be an Eisenstein norm a^2-ab+b^2; its total ramification is 14. The
  THM-4130 later excludes the survivor; THM-4134/4138 empty the Delta_V wall;
  THM-4140/4141 empty the Delta_D wall; and THM-4143 empties the two-term
  wall. Thus this exact-M=8 trichotomy is now empty on the stated seam.
  M>=9, entry, other cells, JC(2), and DC(2) remain OPEN.
source: long-precise-frontiers / 2026-08-24
audit: >
  PASS. The primary certificate expands the complete support, enumerates all
  lower faces in six support strata, verifies both harmless expanded-support
  cancellations, the four integral face functions after Q=sigma^24, all face
  equations and normalizations, three edge discriminants, long-edge leading
  discriminants, eight transverse node determinants, exact A_23 charts, Pick
  and genus ledgers, and 500 Eisenstein-norm controls. A SymPy-free audit uses
  direct supporting-plane inequalities and exact polygon tilings rather than
  triple hull enumeration, plus independent coefficient dictionaries and a
  literal norm sieve. Both paths have the same semantic digest; normal and
  optimized streams byte-match their frozen outputs.
depends_on:
  - THM-3992-reduced-two-three-cusp-jet-repair-and-first-node-residual
  - THM-4007-live-two-three-third-normal-row-five-weight-floor
  - THM-4012-weighted-leading-face-good-elliptic-factor-observer
  - THM-4045-live-two-three-max-seven-hidden-elliptic-tail-no-go
related:
  - THM-4017-sharp-weight-eight-specialization-obstruction-and-newton-ledger
  - THM-4044-sixty-clock-hasse-alias-and-planar-jc-boundary-firewall
  - THM-4130-theta-only-extremal-seam-critical-monodromy-obstruction
  - THM-4134-delta-v-collision-wall-strict-transform-and-high-branch-exclusion
  - THM-4138-delta-v-horizontal-carrier-monodromy-exclusion
external: >
  Tim Dokchitser, "Models of curves over DVRs," arXiv:1807.00025v2,
  Definitions 3.7, 3.9, and 3.12 and Theorem 3.14, through the same
  face/edge-model validity gate used in THM-4045.
script: 04-computation/jc2_weight8_complete_lower_model_thm4053.py
output: 05-knowledge/results/jc2_weight8_complete_lower_model_thm4053.out
script_sha256: 5d5d117d76200a4857d9f12bd1500511fa82d3c28880101816c78688487e11d4
output_sha256: fcc0bf5e890781ecc0e0313844c4f760104f96fe4201be7cfe9f443111bc970e
independent_audit_script: 04-computation/jc2_weight8_complete_lower_model_thm4053_independent_audit.py
independent_audit_output: 05-knowledge/results/jc2_weight8_complete_lower_model_thm4053_independent_audit.out
independent_audit_script_sha256: b8162608f0ea64d6e6d1aa0bfa46cdae8efe8784c8fa775454d7f944cd8fe094
independent_audit_output_sha256: 5d8d828091901594d66bad1bdc8e897e367d7b100275a2ddbc71985d966d8c02
semantic_sha256: 10c489b9be9b5214581e0f6a52609235b0506aaa5d1c6ef431c5fcb36ba6988e
hash_basis: raw LF bytes
---

# THM-4053 -- the exact maximum-weight-eight trichotomy

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED on one live reduced seam;
JC(2) OPEN.** Work over an algebraically closed field of characteristic zero.
This theorem determines the complete lower model at the first weight not
excluded by THM-4045. It kills two top-support strata away from exact
collision walls and turns the remaining stratum into a degree sieve.

## 1. Exact polynomial and forced row

Use THM-3992's normalized reduced `(2,3)` cell

```text
s=xt,               p=s^2+t,               y=sp,
G=gamma*s^2/t+H(p,y),                         gamma!=0.       (1)
```

Restrict to THM-4007's live seam `b=d=0` and suppose the **entire** polynomial
`H` has exact maximum weight eight for `wt(p)=2,wt(y)=3`. Then, with no
ellipsis,

```text
H=lambda*p+alpha*p^2+epsilon*p^3+kappa*y^2
  +phi*p^2*y+delta*p^4+theta*p*y^2,                    (2)
(delta,theta)!=(0,0).
```

Writing `A5=a^5`, the inherited rows are

```text
epsilon/gamma = 2752/(135 A5^3),
delta/gamma +(6/(7A5))(kappa/gamma)=-11392/(105 A5^4), (3)
```

or equivalently

```text
kappa=-5696 gamma/(45 A5^3)-(7A5/6)delta.              (4)
```

Two useful support walls are

```text
kappa=0       iff A5^4 delta/gamma=-11392/105,
kappa=epsilon iff A5^4 delta/gamma=-7936/63.           (5)
```

The second equality deletes one expanded support point but changes no lower
face.

## 2. The trichotomy

Put

```text
Delta_D=phi^2-4kappa delta,
Delta_V=phi^2-4epsilon theta.                          (6)
```

Every hypothetical Keller pair in this exact maximum-eight cell must lie in
one of the following unresolved sets:

```text
A. theta=0, delta*kappa!=0, Delta_D=0;                 (7a)
B. delta=0, theta!=0, Delta_V=0;                       (7b)
C. delta*theta!=0, delta+theta=0;                      (7c)
D. delta=0, theta!=0, Delta_V!=0, and the generic
   fibre degree n is an Eisenstein norm a^2-ab+b^2.    (7d)
```

Everything else at exact maximum weight eight is impossible. More precisely:

- If `theta=0,delta!=0`, then `kappa=0` is impossible, and
  `kappa!=0,Delta_D!=0` is impossible. Thus only `(7a)` remains.
- If `delta*theta!=0`, then `delta+theta!=0` is impossible. Thus only `(7c)`
  remains in the two-term top stratum.
- If `delta=0,theta!=0`, then `(7b)` is the collision wall and `(7d)` is the
  sole nonresonant survivor.

The theorem does not exclude `(7a)--(7d)`, and it does not prove JC(2).

## 3. Complete lower-face inventory

Put `Q=q^-1`. The generic source fibre is

```text
F_Q(S,P)=(S^2-P)(1-QH(P,SP))+gamma*Q*S^2=0.           (8)
```

The coefficient of `S^2` is the unit `1+gamma Q`, so it contributes one
height-zero point, not two lifted points. Expanding `(8)` gives the two base
points

```text
(2,0,0), (0,1,0)                                      (9)
```

and, at height one, the endpoints

```text
term          endpoints in (S exponent,P exponent)    shared coefficient
lambda p      (2,1),(0,2)
alpha p^2     (2,2),(0,3)
epsilon p^3   (2,3),(0,4)                             kappa-epsilon at (2,3)
kappa y^2     (4,2),(2,3)
phi p^2y      (3,3),(1,4)
delta p^4     (2,4),(0,5)                             theta-delta at (2,4)
theta py^2    (4,3),(2,4).                            (10)
```

Exact lower-hull enumeration gives four possible planes:

```text
L: z=(i+2j-2)/8,
D: z=(i+j-2)/4,
V: z=(j-1)/3,
R: z=i/2-1.                                           (11)
```

After removing monomial factors, the face equations are:

```text
top support              lower faces       normalized equations

delta!=0,theta=0         L,D               L:(S^2-P)(1-delta P^4)
                                            D:1-kappa S^2P^2-phi SP^3-delta P^4

delta=0,theta!=0         L,V,R             L:(S^2-P)(1-theta S^2P^3)
                                            V:-1+P^3(epsilon+phi S+theta S^2)
                                            R:1-S^2P^2(kappa+theta P)

delta*theta!=0           L,(R if kappa!=0) L:(S^2-P)(1-delta P^4-theta S^2P^3)
                                            R:1-S^2P^2(kappa+theta P).    (12)
```

For `delta!=0,theta=0,kappa=0,phi=0`, the `D` polygon contracts completely
and only `L` remains. For `kappa=0,phi!=0`, `D` remains rational. In the
two-term stratum `kappa=0` deletes `R` but not the main component.

The complete six-stratum genus ledger is

```text
stratum                         generic genus   graph rank   abelian dimension
delta-only, kappa!=0                  8             7              1
delta-only, kappa=0, phi!=0           7             7              0
delta-only, kappa=phi=0               4             4              0
theta-only                            8             7              1
both top terms, kappa!=0              9             7              2
both top terms, kappa=0               9             7              2.      (13)
```

The cancellations `kappa=epsilon` at `(2,3,1)` and `theta=delta` at
`(2,4,1)` preserve the displayed lower subdivisions. In particular,
`theta=delta` is not a resonance wall.

## 4. Normalized positive-genus components

### 4.1 The delta-only side

On `D`, put `T=SP` and

```text
W=2kappa T+phi P^2.
```

Then

```text
W^2=4kappa+(phi^2-4kappa delta)P^4
   =4kappa+Delta_D P^4.                               (14)
```

If `kappa Delta_D!=0`, this is a smooth genus-one quartic with `j=1728`.
It is the unique positive-genus component. If `kappa=0,phi!=0`, the `D`
equation is linear in `T` and rational. If `kappa=phi=0`, `D` disappears and
the entire genus is graph-theoretic. The repeated-edge case
`kappa delta!=0,Delta_D=0` is deliberately left open.

### 4.2 The theta-only survivor

On `V`, put

```text
X=P^-1,                 W=2theta S+phi.
```

Then

```text
W^2=Delta_V+4theta X^3.                               (15)
```

For `Delta_V!=0`, this is a smooth elliptic curve with `j=0`. It is the
unique positive-genus component. Over the algebraically closed residue field,
choose roots and set

```text
X=(Delta_V/(4theta))^(1/3) x,       W=sqrt(Delta_V)y;
```

equation `(15)` becomes

```text
E_0:y^2=x^3+1.                                          (16)
```

When `Delta_V=0`, the outer cubic edge is cuspidal; its stable replacement is
not determined here.

### 4.3 The two-term top

The nonrational factor of `L` is

```text
delta P^4+theta S^2P^3=1.                             (17)
```

THM-4012 identifies its normalization with the Bolza genus-two curve and
proves

```text
Jac(C_B) isogenous to E_8000^2,
End^0(E_8000)=Q(sqrt(-2)),
Hom(Jac(C_B),E_0)=0.                                  (18)
```

It is the only abelian component in this stratum. The right face, when
present, is rational.

## 5. Complete-model validity gate

The single base change

```text
Q=sigma^24                                             (19)
```

makes all four face functions integral:

```text
L=3i+6j-6,       D=6i+6j-12,
V=8j-8,          R=12i-24.                            (20)
```

Every face and edge denominator is one and every face multiplicity is one.
The nontrivial outer-edge polynomials are, up to units and reversal,

```text
delta-only:      kappa+phi T+delta T^2,
theta-only:      theta+phi T+epsilon T^2,
two-term:        (T-1)(delta T+theta).                 (21)
```

Their discriminants are exactly

```text
Delta_D,              Delta_V,              (delta+theta)^2.  (22)
```

These are the three collision walls in `(7)`. The equality `theta=delta`
gives discriminant `4delta^2`, so its support cancellation is harmless.
Every other edge is linear or the visibly squarefree quadratic
`1+gamma Q-Qkappa T^2`. The two long generic outer edges have leading
discriminants

```text
-256 delta^3 Q^3+O(Q^4)        when delta!=0,
-27 epsilon^2 Q^2+O(Q^3)       when delta=0,           (23)
```

and are therefore smooth over `k((Q))`.

The main factors meet where `P=S^2` and

```text
1-(delta+theta)S^8=0.                                 (24)
```

The determinant of their two gradients is

```text
-8(delta+theta)S^7.                                   (25)
```

Thus there are exactly eight transverse nodes in the delta-only and
theta-only strata, and in the two-term stratum off `delta+theta=0`. In the
main chart, write

```text
H_sigma=delta P^4+theta S^2P^3
       +sigma^3 phi SP^3
       +sigma^6(epsilon P^3+kappa S^2P^2)
       +sigma^12 alpha P^2+sigma^18 lambda P,
U=S^2-P,                    Z=(1-H_sigma)/S^2.         (26)
```

The exact completed local equation is

```text
UZ=-gamma sigma^24.                                   (27)
```

Each node is therefore an `A_23` smoothing; regular resolution adds only a
rational chain.

There is also a generic-completion gate. On `(8)`, `P=S^2` would force
`gamma Q S^2=0`, impossible on the torus. Hence

```text
t=P-S^2,                    x=S/t                     (28)
```

recovers the original source function field. The smooth toric completion is
the actual generic source curve, not an auxiliary curve. Dokchitser's
face/edge model theorem now covers every smooth face, edge, and boundary
locus, and `(27)` resolves the only remaining nodes. Consequently the genus
and positive-genus inventories in `(13)--(18)` are the complete regular
special fibres off the named walls; no hidden abelian component remains.

## 6. The two no-Hom contradictions

Since `q=sigma^-24`, put

```text
A=sigma^-8 X,                 C=sigma^-12 Y
```

in the target elliptic fibre. Its special fibre is the same good curve

```text
E_0:Y^2=X^3+1.                                         (29)
```

Resolve the rational extension of the finite nonconstant generic fibre map.
Every extra exceptional curve is rational. A rational component maps
constantly to `E_0`. In the delta-only stratum, the only possible abelian
component is either absent or has `j=1728`, hence has no nonzero Hom to the
`j=0` target. In the two-term stratum, `(18)` gives the same conclusion for
the Bolza component. Thus every special-fibre component maps constantly; by
connectedness their constants agree.

For a relatively ample target line bundle `L`, degree conservation gives

```text
deg(phi_generic^*L)=sum_i m_i deg(phi_i^*(L|E_0)).    (30)
```

The right side would be zero, while the finite generic map makes the left
side positive. This contradiction proves the exclusions in Section 2.

## 7. Eisenstein degree and ramification in the survivor

In the nonresonant theta-only stratum, `(16)` is the unique positive-genus
component and has multiplicity one. All rational components contribute zero
to `(30)`. If `n` is the generic fibre-map degree, then

```text
n=degree(phi_E:E_0->E_0).                             (31)
```

After translating the image of the origin, `phi_E` is an endomorphism of
`E_0`. In characteristic zero

```text
End(E_0)=Z[zeta_3],
degree(a+b zeta_3)=Norm(a+b zeta_3)=a^2-ab+b^2.       (32)
```

Therefore `n` must be an Eisenstein norm. Equivalently, every prime
`p=2 mod 3` occurs in `n` to even valuation. Since the generic source has
genus eight and the target genus one, `n>1` and Riemann--Hurwitz gives the
exact total ramification budget

```text
deg R_phi=2*8-2=14.                                   (33)
```

The elliptic specialization `(31)` is an isogeny and hence unramified. The
allocation of the fourteen generic ramification units among the rational,
node, and boundary charts is a new exact target, not a conclusion here.

## 8. Historical frontiers and later closure

At the time of this theorem, the maximum-eight problem was reduced to three
explicit collision walls:

1. stable reduction of the repeated quadratic edge `Delta_D=0` in `(7a)`;
2. stable reduction of the cuspidal cubic edge `Delta_V=0` in `(7b)`;
3. the eight-attachment collision `delta+theta=0` in `(7c)`.
THM-4130 later excludes the smooth survivor `(7d)`, after
THM-4103/4120/4126 fix its boundary response, degree, and nonproperness.
THM-4134 later stratifies `Delta_V=0` and excludes degrees `20,19`; THM-4138
then excludes its horizontal-BC degrees `16,15`, emptying that wall.
THM-4140/4141 empty `Delta_D=0`, and THM-4143 empties
`delta+theta=0`. Consequently every alternative `(7a)--(7d)` is now excluded
on this exact-`M=8` seam. Weight at least nine, entry, other reduced cells,
`JC(2)`, and `DC(2)` remain open.

**QED.**
