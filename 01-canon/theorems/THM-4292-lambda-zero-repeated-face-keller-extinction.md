---
id: THM-4292
title: "Lambda-zero repeated-face Keller extinction"
status: >
  PROVED RELATIVE TO THM-4103/4230/4272/4291 + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED. On the exact-M=12 wall W=Lambda=0, U!=0, every
  exceptional divisorial component above the unique A23 contact has strictly
  positive Keller-differential order, including every repeated-root and
  subsequent Newton-Puiseux face. Hence every such component is
  Keller-constant in characteristic zero. THM-4294 combines this uniform tail
  extinction with exact central order nine to close the `W=Lambda=0` slice.
  The general coefficient wall, seam entry, and JC(2) remain open.
source: root/jc2-signal-20260830
depends_on:
  - THM-4103-jc23-theta-boundary-ramification-and-degree-response
  - THM-4230-exact-weight-twelve-cyclotomic-prym-planar-jacobian-squeeze
  - THM-4272-lambda-zero-a23-contact-and-e0-infinity-jet-obstruction
  - THM-4291-lambda-zero-genus-five-tail-degree-forty-two-equivariant-hostile
related:
  - THM-4289-a23-blowdown-observer-kahler-dualizing-quotient
  - THM-4290-exact-weight-twelve-deck-equivariant-visible-quotient-exclusion
  - THM-4294-lambda-zero-central-keller-extinction
primary_script: 04-computation/jc23_lambda_zero_repeated_face_keller_extinction_thm4292.py
primary_output: 05-knowledge/results/jc23_lambda_zero_repeated_face_keller_extinction_thm4292.out
primary_script_sha256: e334a32d20c1b0d8421622c94dbb6d38aa380730d663f416a5cf15481534ab19
primary_output_sha256: 0ed06f7ea39074ea79036e89091da5500dc7da54982a7e2d3acf734fc7ebe647
independent_audit_script: 04-computation/jc23_lambda_zero_repeated_face_keller_extinction_independent_audit_thm4292.py
independent_audit_output: 05-knowledge/results/jc23_lambda_zero_repeated_face_keller_extinction_independent_audit_thm4292.out
independent_audit_script_sha256: 16c81eab23c2c583df5f35a9bfb7d726396c4a91cebff09a9e207da46f423bd4
independent_audit_output_sha256: 30a2e6ec82ea9db6a0d0cc155d88912ec5508ccc8c81461e0abeb1a49a68b7fa
hash_basis: raw LF bytes
audit: >
  PASS. The primary SymPy path reconstructs the literal finite lower-row
  polynomial, derives the four moving-critical coefficients, verifies the
  maximal cancellation witness, and checks every collision scale and form
  order. A standard-library-only clean-room path independently reconstructs
  the critical series and valuation table. Normal, optimized, and
  fixed-hash-seed streams byte-match the frozen outputs.
---

# THM-4292 -- Lambda-zero repeated-face Keller extinction

**PROVED RELATIVE TO THM-4103/4230/4272/4291 + VERIFIED-EXACT +
INDEPENDENTLY AUDITED. THE REPEATED-TAIL GAP IS CLOSED LOCALLY. THM-4294
LATER CLOSES THE `W=Lambda=0` SLICE; THE GENERAL WALL, SEAM ENTRY, AND `JC(2)`
REMAIN OPEN.**

## 1. Statement

Work over `C` in THM-4230's reduced `(2,3)`, exact-weight-twelve family. Put

```text
W=0,                 Lambda=U+Z=0,                 U!=0,             (1)
```

so `Z=-U` and `D=4U^2!=0`. Retain **all** allowed lower rows. Let `Q` be the
unique length-twelve infinity contact of THM-4272.

> **Theorem.** After any finite base change, let a normal or regular model of
> the exact source carry the extended Keller graph. Every noncentral
> divisorial component centered above `Q` has strictly positive vertical
> order for the pullback of the good-target invariant differential. Therefore
> the Keller map restricts constantly to every such component in
> characteristic zero.

This includes initial quadratic squares, translations through all multiple
discriminant roots, and subsequent Newton-Puiseux faces. In particular every
positive-genus exceptional tail at the `A_23` contact is Keller-constant.

The closest proved mechanism is THM-4291's generically non-square valuation
lemma. Its canonical hostile is the deck-invariant prepared quadratic in
Section 7: outside the literal finite lower-row support it produces a genuine
genus-two, degree-two equivariant tail with **zero** differential order. The
corrected near miss is therefore that `B=t^12*unit` and deck symmetry alone
do not control repeated roots. The decisive sidecar is the exact four-step
critical-value ladder forced by the finite `Hhat` support.

## 2. Literal equation and global localization

Use THM-4291's infinity variables

```text
b=1/S,             r=P/S^2,             q=r-1,             t=sigma*b. (2)
```

The exact source equation is

```text
F=(1-r)(b^12-Hhat(r,t))-sigma^12*b^12/2
 =q(Hhat(r,t)-b^12)-t^12/2.                                (3)
```

On `(1)`, the complete lower-row polynomial is

```text
Hhat=U(r^6-r^4)
 +t (alpha_11 r^5+beta_11 r^4)
 +t^2(upsilon_5 r^5+xi_10 r^4)
 +t^3(eta r^4+zeta_3 r^3)
 +t^4(Delta r^4+Theta r^3)
 +t^5 Phi r^3
 +t^6(-(1376/135)r^3+K r^2)
 +(8/3)t^8r^2-3t^10r,

K=2848/45-(7/6)Delta.                                      (4)
```

Because `U!=0`, this is Weierstrass-quadratic in `q` at the contact. More
intrinsically it has the exact form

```text
F=q ell+q^2 V(q,t)-t^12/2,
ell=Hhat(1,t)-b^12,                 V(0,0)=2U.              (5)
```

The raw `sigma=0` toric fibre is exactly `R union C`: `R` is rational, `C`
is smooth of genus seven, and both have multiplicity one. THM-4272 proves
that their complete intersection on `(1)` is the single fat point `12Q`,
analytically the `A_23` contact. The conditions `U*Z*D!=0` make `C`
nondegenerate at every other torus and boundary stratum. Thus the total
family is smooth away from `Q`, and finite base change preserves that
smoothness.

Consequently every new component introduced by normalization, regularization,
or semistable reduction is centered above `Q`; fixed toric-subdivision
components are rational. The toric closure has arithmetic genus eighteen,
but its wall normalization can lose genus at this same boundary point. No
smooth-genus budget is assumed here. It is enough to control all divisorial
valuations of `(5)` centered at `sigma=b=q=0`.

## 3. Newton alternatives before repetition

Let `v=v_E` be such a divisorial valuation and write

```text
s=v(sigma)>0,            beta=v(b)>0,            tau=s+beta,
lambda=v(ell).                                                (6)
```

The Newton alternatives of `(5)` are

```text
lambda<6tau:  v(q_1),v(q_2)=lambda,12tau-lambda;
lambda>6tau:  q=t^6y, initial 2Uy^2-1/2;
lambda=6tau:  q=t^6y, initial 2Uy^2+ell_bar*y-1/2.       (7)
```

The first two rows are separable. In the third row the quadratic is double
exactly when

```text
ell_bar^2+4U=0.                                             (8)
```

If the initial quadratic is not a square, THM-4291 applies directly. With
`gamma=v(q)` and

```text
N=min(2gamma,gamma+v(A),12(s+beta)),                         (9)
```

one has `v(F_q)=N-gamma` and

```text
v_E(omega_0)>=3s+5beta>0.                                  (10)
```

Only iterated refinements of `(9)` remain.

## 4. Where a multiple discriminant root can persist

Put

```text
h(t)=Hhat(1,t)
 =c_1t+c_2t^2+c_3t^3+c_4t^4+c_5t^5+ct^6
  +(8/3)t^8-3t^10,                                         (12)

c_1=alpha_11+beta_11,        c_2=upsilon_5+xi_10,
c_3=eta+zeta_3,              c_4=Delta+Theta,
c_5=Phi,                     c=7168/135-(7/6)Delta.         (13)
```

### 4.1. The range `beta<s`

Normally `b^12`, or the first `c_k t^k` with `k<6`, gives
`lambda<6tau`. Cancellation can occur only at

```text
12beta=k tau,                       k<6.                    (14)
```

After translating through that cancellation and normalizing the Puiseux
branch, the terminal transverse coordinate is `x=ell/t^6`, and the face is

```text
2Uy^2+xy-1/2=0,              Disc_y=x^2+4U.                 (15)
```

Both discriminant roots are nonzero and simple because `U!=0`. Thus `(14)`
can create simple ramification, but not another multiple discriminant root.
Normalization cancels the corresponding simple derivative zero, and the
positive bound `(10)` survives.

### 4.2. The balanced range `beta=s`

Terms `c_k t^k`, `k<6`, again fall under the preceding non-square cases.
Once `c_1=...=c_5=0`, put `X=b/sigma`. The face becomes

```text
2Uy^2+(c-X^6)y-1/2=0,
D(X)=(c-X^6)^2+4U.                                        (16)
```

At a nonzero root,

```text
D'(X)=-12X^5(c-X^6)                                      (17)
```

cannot vanish, since that would force `U=0`. The only possible multiple root
is `X=0`; it occurs exactly when

```text
c^2+4U=0,                  c!=0,
D(X)=X^6(X^6-2c).                                           (18)
```

Centering at `X=0` is precisely passage to `beta>s`.

## 5. The exact four-step ladder for `beta>s`

Here repetition forces

```text
c_1=...=c_5=0,                    U=-c^2/4,                (19)
```

and hence

```text
beta_11=-alpha_11,       xi_10=-upsilon_5,
zeta_3=-eta,             Theta=-Delta,             Phi=0. (20)
```

Put

```text
q=t^6y,             y_0=1/c,             epsilon=b^12/t^6,
v(epsilon)=d_b=6(beta-s).                                  (21)
```

Direct substitution into the literal polynomial `(4)`, not an arbitrary
prepared quadratic, gives

```text
F/t^12=
 -(c^2/2)(y-y_0)^2
 +alpha_11 t y^2
 +t^2(upsilon_5 y^2+(8/3)y)
 +eta t^3y^2
 +t^4(Delta y^2-3y)
 -epsilon y+O(t^6).                                       (22)
```

Every omitted higher-`q` term first enters at `t^6`, because `q=t^6y`.
There is no `t^5` term. Solving for the moving critical point gives the
recursive critical-value coefficients

```text
C_1=alpha_11/c^2,
C_2=upsilon_5/c^2+8/(3c)              after C_1=0,
C_3=eta/c^2                            after C_1=C_2=0,
C_4=(Delta+32/9-3c)/c^2                after C_1=C_2=C_3=0. (23)
```

The `32/9` term is load-bearing. Under `C_2=0`,

```text
upsilon_5=-8c/3,       P_2'(y_0)=-8/3,       P_0''=-c^2,   (24)
```

so completing the moving critical square contributes `32/(9c^2)`.

Let `j` be the first nonzero index in `(23)`. Compare

```text
d_j=j(s+beta),                       d_b=6(beta-s).          (25)
```

Equality occurs only at

```text
beta/s=(6+j)/(6-j)=7/5,2,3,5.                                (26)
```

At equality a toric residue coordinate gives the critical-value face

```text
X^j(C_j-(1/c)X^(6-j)).                                      (27)
```

Every nonzero root is simple. The root `X=0` passes to a larger value of
`beta/s`. If all four coefficients vanish, the absent `t^5` row and

```text
6(beta-s)<6(s+beta)                                         (28)
```

force the `b^12q` term to split the root before any `O(t^6)` correction.
This exhausts all successive repeated faces.

## 6. Keller-form orders and constancy

THM-4103 and the exact source/target scalings give, up to the single global
residue-orientation sign,

```text
omega_0=dx_E/(2y_E)=-sigma^9 b^10 db/F_q.                   (29)
```

On `(21)`, write `P=F/t^12`. Then `F_q=t^6P_y`. If the first critical
splitter has gap `d`, normalization gives `v(P_y)=d/2`, while
`v(d_rel b)>=beta`. Hence

```text
v_E(omega_0)>=3s+5beta-d/2.                                 (30)
```

The complete table is

```text
first splitter    d                 Keller-form lower bound
b^12 q            6(beta-s)         6s+2beta
C_1 t             s+beta            (5/2)s+(9/2)beta
C_2 t^2           2(s+beta)         2s+4beta
C_3 t^3           3(s+beta)         (3/2)s+(7/2)beta
C_4 t^4           4(s+beta)         s+3beta.                 (31)
```

Every entry is strictly positive. At a simple collision root the normalized
local form is

```text
X-X_0=unit*Z^2,                    P_y=unit*Z.                (32)
```

Thus `dX=unit*Z dZ`; the apparent extra denominator zero cancels rather than
lowering `(30)`. Ramification multiplies `s,beta,d` together and preserves
positivity. At `X=0` the next ratio case applies, and blowups of an already
smooth normalized point cannot lower positive vertical order.

On a regular model carrying the Keller extension, `(29)` is the pullback of
a genuine relative Kähler differential, not merely a raw dualizing residue.
Its restriction to every exceptional component is zero. A nonconstant map
between smooth characteristic-zero curves has nonzero differential, so every
such component map to `E_0` is constant. This proves the theorem.

## 7. Hostile controls and the exact boundary

The constant-term valuation and deck action alone are insufficient. Put

```text
x=b^12,                  t=sigma*b,
P=q^2+2t^6q+t^12(1+t^8+xt^2).                               (33)
```

This deck-invariant prepared quadratic has `B=t^12*unit`, but

```text
Disc_q(P)=-4t^14(t^6+x)
         =-4sigma^14b^20(b^6+sigma^6).                       (34)
```

Under `b=sigma X` and `2q+A=sigma^20X^10Z`, it produces

```text
Z^2=-4(1+X^6),                 sigma^9b^10db/(2q+A)=dX/Z.    (35)
```

The tail has genus two and vertical order zero. After rescaling `Z`,

```text
(A_E,C_E)=(X^-2,Z/X^3)                                      (36)
```

is an equivariant degree-two map to `E_0:C_E^2=A_E^3+1`.
Multiplication by `-1/(2(1+t^8+xt^2))` even normalizes the constant term to
`-t^12/2`. This hostile is not realizable by `(4)`; equation `(22)` and
ladder `(23)` are exactly what exclude it.

Inside the exact family, the naive claim that the fixed `t^8q` term always
splits the root is also false. The allowed coefficients

```text
c=5152/405,                 Delta=4672/135,        Theta=-Delta,
upsilon_5=-8c/3,            xi_10=8c/3,            U=-c^2/4,
alpha_11=beta_11=eta=zeta_3=Phi=0                         (37)
```

make `C_1=C_2=C_3=C_4=0`. At `s=1,beta=6`, the alleged `t^8q` gap `14` and
the `t^4` critical gap `28` both cancel. The first splitter is `b^12q` at
gap `30`, yet `(31)` gives Keller-form order at least `18`.

For THM-4291's displayed specialization, by contrast,

```text
c=7168/135,                   C_2=8/(3c)!=0,                (38)
```

so its shorter two-gap calculation was valid. The longer ladder is needed
for the full lower-row theorem.

## 8. Scope and reproduction

This theorem closes the repeated-tail obstruction left by THM-4291 and proves
that no positive-genus exceptional tail can carry Keller degree. It does not
by itself determine the response degree on `Lambda=0`, descend the resolved
graph through the raw `A_23` contact, close the wall, prove exact-`M=12` seam
entry, or prove `JC(2)` or `DC(2)`. THM-4294 later closes this
`W=Lambda=0` slice without raw contact descent by proving the central map
constant and conserving degree.

Run

```bash
python3 -B 04-computation/jc23_lambda_zero_repeated_face_keller_extinction_thm4292.py
python3 -B -O 04-computation/jc23_lambda_zero_repeated_face_keller_extinction_thm4292.py
PYTHONHASHSEED=0 python3 -B 04-computation/jc23_lambda_zero_repeated_face_keller_extinction_independent_audit_thm4292.py
```

The first path uses exact SymPy algebra. The second is a clean-room
standard-library reconstruction with exact rational arithmetic. **QED.**
