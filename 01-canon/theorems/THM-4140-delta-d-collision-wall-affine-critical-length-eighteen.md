---
id: THM-4140
title: "Delta-D collision-wall square top and affine critical length eighteen"
status: >
  PROVED RELATIVE TO THM-3992/4053 + VERIFIED-EXACT + INDEPENDENTLY
  AUDITED. On the live delta-only exact-M=8 wall Delta_D=0, the complete
  weight-eight response is a square, the strict critical resultant has
  degree fourteen, and every hypothetical Keller realization has exactly
  eighteen reduced affine critical points. This is a critical ledger, not
  an exclusion: stable reduction, the boundary packet, the finite response
  degree, the wall itself, JC(2), and DC(2) remain OPEN.
source: codex-planar-jacobian-cycle6-63-20260825
depends_on:
  - THM-3992-reduced-two-three-cusp-jet-repair-and-first-node-residual
  - THM-4053-jc2-live-max-eight-trichotomy-and-eisenstein-survivor
related:
  - THM-4130-theta-only-extremal-seam-critical-monodromy-obstruction
  - THM-4134-delta-v-collision-wall-strict-transform-and-high-branch-exclusion
  - THM-4138-delta-v-horizontal-carrier-monodromy-exclusion
script: 04-computation/jc23_delta_d_collision_wall_critical_length_thm4140.py
output: 05-knowledge/results/jc23_delta_d_collision_wall_critical_length_thm4140.out
independent_audit_script: 04-computation/jc23_delta_d_collision_wall_critical_length_thm4140_independent_audit.py
independent_audit_output: 05-knowledge/results/jc23_delta_d_collision_wall_critical_length_thm4140_independent_audit.out
script_sha256: 411e227de91fab8f4b642b816c8defd28d0c793f99ebae00661a17429f036bf7
output_sha256: 252d1e3669f578f6094075fb7ce800142209f887f648b6643eb2fcd070088833
independent_audit_script_sha256: 2641f7436814ec494051c91328feb909f162b64e6a8d72fc58dd63cf6c864740
independent_audit_output_sha256: 43de5d2f9e95c1976f45bd9d9cb74d01dc6227ac91152fd68b0436602a1e76ec
semantic_sha256: 1917ba89584d8537a7df2c686f5713a07362bca88ed129a93b1cfb1c95db21fa
hash_basis: raw LF bytes
primary_audit: >
  PASS. A standalone SymPy calculation parameterizes the complete live
  Delta_D wall, constructs the normalized (X,T) source polynomial, factors
  its critical resultant, restores the T=0 stratum, freezes the monic
  degree-fourteen coefficient ledger, and checks a squarefree rational
  control. Normal, optimized, and hash-seeded replays byte-match.
independent_audit: >
  ACCEPT. A separate (s,p) calculation imports no primary code, derives two
  different critical equations, eliminates s, identifies the fake p^2
  factor, restores the p=0 and t=0 strata, and reproduces the degree and
  length ledger with the same semantic hash.
---

# THM-4140 -- the `Delta_D` wall has critical length eighteen

**PROVED RELATIVE TO THM-3992/4053 + VERIFIED-EXACT + INDEPENDENTLY
AUDITED; THE WALL AND JC(2) REMAIN OPEN.** Work over `C` and retain the exact
reduced `(2,3)` cell and maximum-weight-eight hypotheses of THM-4053.

## 1. Theorem and inheritance

> **Theorem.** In the live delta-only collision wall
>
> ```text
> theta=0,        delta*kappa!=0,        Delta_D=phi^2-4kappa delta=0,
>                                                               (1)
> ```
>
> the complete normalized weight-eight row is a perfect square. The affine
> critical scheme of the exact source polynomial of any hypothetical Keller
> realization is reduced and has length
>
> ```text
> #Crit_aff(E(A,C))=18.                                 (2)
> ```
>
> Its critical values are the two target nodal values, so if `r_i` counts
> affine preimages of the two target nodes, then
>
> ```text
> r_0+r_1=18.                                           (3)
> ```

The closest proved mechanism is THM-4130's exact critical resultant on the
smooth theta-only seam. The canonical hostile is THM-4134: on a collision
wall a leading resultant row can vanish because roots escape a chart, so a
smooth-seam count cannot be specialized. The least-used sidecar is the
repeated quadratic edge itself. Parameterizing its root turns the entire top
response into one square before elimination.

## 2. Exact wall parameterization and square mechanism

THM-3992/4053 permit the same normalization as THM-4130. Put

```text
X=a^(5/2)x,       T=t/a^5,       P=T+X^2T^2,       Y=XTP,
Delta=a^17 delta, K=a^12 kappa, Phi=a^(29/2)phi.    (4)
```

After division by `a^3`, the forced row and collision equation are

```text
K=2848/45-(7/6)Delta,          Phi^2=4KDelta.          (5)
```

Because `KDelta!=0`, choose the repeated-root parameter `r` by

```text
K=Delta r^2,                  Phi=2Delta r.            (6)
```

Then `(5)` gives the complete rational chart

```text
Delta=5696/[15(6r^2+7)],
K=Delta r^2,                  Phi=2Delta r,
r!=0,                         6r^2+7!=0.               (7)
```

The last denominator condition is not an extra genericity assumption: if it
vanished, the forced nonzero right side in `(5)` would give a contradiction.
The whole maximum-weight-eight row factors as

```text
K Y^2+Phi P^2Y+Delta P^4 = Delta(P^2+rY)^2.            (8)
```

Thus the collision retains the same monomial support while gaining a square
multiplicity. Equation `(8)`, rather than the bare equality `Delta_D=0`, is
the structural reason to recompute the critical strict transform.

## 3. The `(X,T)` critical resultant

The complete normalized source polynomial is

```text
G=-X^2T/2-3P+(8/3)P^2-(1376/135)P^3
  +K Y^2+Phi P^2Y+Delta P^4.                           (9)
```

There are no omitted terms. Since `G_X` is divisible by `T`, put

```text
f=G_X/T,                         h=G_T.                 (10)
```

Their generic degrees in `X` are `7,8`, with common leading coefficient

```text
8Delta T^7.                                             (11)
```

Exact elimination over `Q(r)` gives

```text
Res_X(f,h)=T^42(6T+1)^2 Q_14(T;r),
deg_T Q_14=14,                                          (12)
```

where

```text
[T^14]Q_14=
  48947063421231612219454924509816611957374976 r^10
  /[8649755859375(6r^2+7)^9],                           (13)

Q_14(0)=
 -139887797298879228245180416
  /[3796875(6r^2+7)^6].                                (14)
```

Both are nonzero on `(7)`. The monic fifteen-coefficient row has digest

```text
1ef62ff920010fcc5eb64c2a8e09e90a03b6c1b4b0b49beb17f4c63ec509ec74. (15)
```

The factor `T^42` is a Sylvester degree-drop artefact, not a finite critical
multiplicity. Indeed,

```text
f(X,0)=-X,                    h(X,0)=-(X^2+6)/2,        (16)
```

so `f,h` have no common root there. For `T!=0`, `(7),(11)` exclude a common
root at `X=infinity`. Consequently `(12)` counts critical-scheme length
sixteen in the `T!=0` chart, with intersection multiplicity.

The universal pair accounts for the displayed `(6T+1)^2` divisor in `(12)`:

```text
T=-1/6,                       X^2=6,       G=1/2.       (17)
```

At special live values of `r`, `Q_14(-1/6;r)` may vanish as well. This is
additional projection multiplicity over the same `T`-coordinate, not a
failure of the universal pair or a change in the total length.

The actual critical ideal is `(Tf,h)`. Restoring `T=0` in that ideal gives

```text
T=0,                          X^2=-6,      G=0.         (18)
```

These are two further points. Equations `(12),(17),(18)` therefore give the
length `16+2=18` in `(2)`.

## 4. Why resultant length is point count

Let

```text
E(U,V)=V^2-U^3+(3/4)U+1/4.                             (19)
```

Its two critical points are `(-1/2,0),(1/2,0)`, with values `0,1/2` and
Hessian determinants `+6,-6`. For a hypothetical Keller pair `F=(A,C)`,

```text
grad(E o F)=DF^t grad(E),
Hess(E o F)=DF^t Hess(E)DF                              (20)
```

at a critical point. Since `det DF=1`, every affine critical point of `G`
is Morse and lies over one of the two target nodes. For `T!=0`, the Jacobian
of `(f,h)` is `T^-1 det Hess(G)`, hence each local intersection counted by
`(12)` has multiplicity one. Equations `(16)--(20)` prove both reducedness
and `(2),(3)`.

This implication is conditional only on the inherited Keller realization.
The displayed source polynomial at an arbitrary parameter `r` is not itself
claimed to come from a Keller pair.

## 5. Independent `(s,p)` projection

The independent audit sets

```text
s=XT,                  p=T+s^2,                  t=p-s^2. (21)
```

Then

```text
G=-s^2/(2t)-3p+(8/3)p^2-(1376/135)p^3
  +K s^2p^2+Phi sp^3+Delta p^4.                       (22)
```

For `pt!=0`, its critical equations are equivalently

```text
A=-s+t^2p(2Ks+Phi p),
B=-6+(32/3)p-(2752/45)p^2
    +6Ks^2p+7Phi sp^2+8Delta p^3.                     (23)
```

The exact identities

```text
t^2 G_s=pA,                   2t^2 G_p=t^2B-sA        (24)
```

audit that equivalence. Eliminating `s` gives

```text
Res_s(A,B)=p^2 R_14(p;r),       deg_p R_14=14,         (25)
```

with nonzero live-wall endpoints

```text
[p^14]R_14=
 -1539884872666062544522946019328 r^4
  /[512578125(6r^2+7)^6],

R_14(0)=-112127901696 r^4/[25(6r^2+7)^2].              (26)
```

The monic coefficient digest is

```text
e5d793f30b0db3def23bab3f1ed08ab247b98fb01ab91477099f757b4c623848. (27)
```

Here `p^2` is again fake: `A(s,0)=-s` and `B(s,0)=-6`. Before the division
used in `(23)`, the direct `p=0` gradient gives

```text
s^2=1/6,                   t=-1/6,                 G=1/2, (28)
```

which is the pair `(17)`. The birational chart omits `t=0`; reconstructing
the polynomial chart independently gives precisely `(18)`. The leading
`s`-coefficients `2Kp,6Kp` exclude a common root at infinity for `p!=0`.
Thus this different projection again gives `14+2+2=18`.

## 6. Controls, consequence, and boundary

At `r=1`, the primary `Q_14` is squarefree modulo `101`; independent checks
at `103,107,109` reproduce a `p`-valuation of two and squarefree residual
degree fourteen. The excluded specialization `r=0` has `K=Phi=0` and shows
why the live gate is load-bearing. The fake `p^2` factor and the collapsed
`t=0` chart are the two hostile projection controls.

Reproduce both exact routes with

```bash
python3 -B 04-computation/jc23_delta_d_collision_wall_critical_length_thm4140.py
python3 -B -O 04-computation/jc23_delta_d_collision_wall_critical_length_thm4140.py
python3 -B 04-computation/jc23_delta_d_collision_wall_critical_length_thm4140_independent_audit.py
python3 -B -O 04-computation/jc23_delta_d_collision_wall_critical_length_thm4140_independent_audit.py
```

Both routes emit the semantic digest
`1917ba89584d8537a7df2c686f5713a07362bca88ed129a93b1cfb1c95db21fa`.

This theorem supplies the exact nodal fixed-sheet budget for the next
`Delta_D` attack. It does **not** determine the wall's stable boundary
packet, finite horizontal ramification defect, generic cover degree, or
nonproperness transport. Until those are known, THM-4138's orbit-merger
lemma cannot be applied. The `Delta_D` wall, the `delta+theta=0` wall,
maximum residual weight at least nine, other cells, `JC(2)`, and `DC(2)` all
remain open. **QED.**
