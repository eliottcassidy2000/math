---
id: THM-3681
title: "Russell-cylinder Q_dagger universal sixth-stable debt"
status: >
  PROVED + VERIFIED-EXACT; PENDING INDEPENDENT HOSTILE AUDIT.  For the
  degree-eight Q_dagger fold, every arbitrary target two-form satisfies one
  universal retained order-six identity.  If J_0=1 and J_2=J_4=0, it forces
  Lambda(J_6)=-326763520/3^13.  Hence no target two-form, and therefore no
  target pair, has nonzero constant pullback.  Character dilation and formal
  conjugacy close the same Q_dagger compiler for every nonzero critical
  displacement H in t^2 C[t].  This closes one candidate family, not JC(2).
source: kps-s196 / complete retained order-six two-form cokernel scout, 2026-08-21
depends_on:
  - THM-3673-russell-cylinder-monomial-ramification-debt-dilation
  - THM-3675-russell-cylinder-critical-fold-formal-conjugacy-closure
  - THM-3677-russell-cylinder-degree-eight-fourth-debt-parabola
related:
  - THM-3680-russell-cylinder-qdagger-coupled-stable-lift
script: 04-computation/jc2_russell_cylinder_qdagger_sixth_stable_debt_thm3681.py
output: 05-knowledge/results/jc2_russell_cylinder_qdagger_sixth_stable_debt_thm3681.out
script_sha256: 98272d9d0b1f81bc7752e0b78f66832b21e28942987b3980a159197f79542618
output_sha256: 6e3c9869f6c7dd233ab9ebc6b2475654658d7747aa461de66522d6a10936ce8f
semantic_sha256: 204b836a68a19c64404a409d7a6298aa370a1ae9a8607b3100d288e2e384c982
hash_basis: raw LF bytes for files; canonical nonzero row/value ledger for semantic hash
---

# THM-3681 -- Q_dagger universal sixth-stable debt

**PROVED + VERIFIED-EXACT; PENDING INDEPENDENT HOSTILE AUDIT.**  The
zero-fourth-debt candidate survives the actual target-ring equations through
`J_2`, but it meets a universal arbitrary-two-form obstruction two even
layers later.  The obstruction is local at the retained ordinary triple and
does not use `J_1`, `J_3`, `J_5`, closedness, decomposability, or a selected
target representative.

All rings, formal germs, and two-forms are over `C`.  The exact certificate is
rational.

## 1. Fixed fold and notation

Use the Russell-cylinder compiler

```text
D=1+x^2q,
b=(D-1)(D+2)^2,       c=xD(D+2),       e=q(D+3),
y=c/3,                z=e+3,                         (1)
```

and the THM-3677 polynomial

```text
Q_dagger=
 (22868x^8-89583x^6+2916x^5+123684x^4-5832x^3
   -63530x^2+2916x-2187)/2916.                       (2)
```

It has the retained values and slopes

```text
Q(-1,0,1)=(-3,-3/4,-3),
Q'(-1,0,1)=(-9/2,1,9/2),                             (3)
```

and zero universal second and fourth retained debts.  First use the quadratic
fold

```text
q=Q_dagger(x)+t^2,                  w=t.              (4)
```

For an arbitrary target two-form `Omega`, write

```text
E^*Omega=j(x,t) dx^dt,
j(x,t)=sum_(n>=0)t^n J_n(x).                         (5)
```

For `a in {-1,0,1}` let

```text
T_(a,d)(P)=P^(d)(a)/d!,
Lambda(P)=(5P(-1)-18P(0)+13P(1))/18.                 (6)
```

Thus `T_(a,d)` reads the degree-`d` Taylor coefficient, not the raw
derivative.

## 2. The exact sixth-stable identity

Every arbitrary target two-form obeys

```text
0=Lambda(J_6)+R_4(J_4)+R_2(J_2)+R_0(J_0),            (7)
```

where

```text
R_4(P)=
 -(4024/3159)T_(-1,0)(P)+(4024/3159)T_(0,0)(P)
 -(5/27)T_(-1,1)(P)+(13/27)T_(1,1)(P),               (8)

R_2(P)=
 -(116061520/9979281)T_(-1,0)(P)
 +(116061520/9979281)T_(0,0)(P)
 +(26159/28431)T_(-1,1)(P)+(64/81)T_(0,1)(P)
 +(845/2187)T_(1,1)(P)
 +(10/81)T_(-1,2)(P)+(26/81)T_(1,2)(P),              (9)

R_0(P)=
 +(145460987072/3502727631)T_(-1,0)(P)
 +(572438466368/3502727631)T_(0,0)(P)
 +(749179514/89813529)T_(-1,1)(P)
 +(131584/85293)T_(0,1)(P)-(283114/531441)T_(1,1)(P)
 -(56348/85293)T_(-1,2)(P)+(3380/6561)T_(1,2)(P)
 -(20/243)T_(-1,3)(P)+(52/243)T_(1,3)(P).            (10)
```

The alternating lower blocks are not decorative.  The `R_4` row is the
negative of the fourth-debt identity's `J_2` block at `Q_dagger`, and the
first- and second-derivative part of `R_2` is the corresponding negative
`J_0` block.  This provides an internal recurrence control on the newly
computed row.

If

```text
J_0=1,                         J_2=J_4=0,              (11)
```

then all positive Taylor coefficients in `(8)--(10)` vanish and

```text
R_0(1)
 =145460987072/3502727631+572438466368/3502727631
 =326763520/1594323
 =326763520/3^13.                                    (12)
```

Equations `(7)` and `(12)` force

```text
Lambda(J_6)=-326763520/3^13 !=0.                      (13)
```

A constant source density would have `J_0` constant and every `J_n=0` for
`n>0`, contradicting `(13)`.  Therefore no arbitrary target two-form has
nonzero constant pullback through `(4)`.  In particular no form
`dF^# wedge dG^#` from a target pair does, so the `Q_dagger` critical-fold
construction cannot be a planar Keller map.

## 3. Why 252 monomials are complete

At each of the three retained branches, the common target point is

```text
(b,c,e,w)=(0,0,-3,0).
```

The functions `(y,z,w)` are regular formal parameters there.  Write an
arbitrary formal target two-form as

```text
Omega=A(y,z,w) dy^dz+B(y,z,w) dy^dw+C(y,z,w) dz^dw.   (14)
```

The pullbacks of `y,z,w` vanish to source order at least one at every branch.
Consequently a coefficient monomial of total target degree above six cannot
affect any source Taylor coefficient of total degree at most six.  It is
therefore enough, and necessary, to test

```text
3 * #{y^a z^b w^c:a+b+c<=6}
 =3*binom(9,3)=252                                     (15)
```

two-form monomials.

The complete retained row window consists of

```text
J_0 through x-Taylor degree 3:     12 rows,
J_2 through x-Taylor degree 2:      9 rows,
J_4 through x-Taylor degree 1:      6 rows,
J_6 values:                          3 rows,            (16)
```

for `30` rows total.  Exact rational elimination gives

```text
rank=26,                    left-cokernel dimension=4. (17)
```

Exactly one cokernel line is nonzero on the constant `J_0` vector.  Normalize
it so that its `J_6` block is `(5/18,-1,13/18)`; its complete nonzero ledger
is exactly `(7)--(10)`.  Thus the finite matrix proves the identity for every
formal two-form jet, and hence for every global regular target two-form.

The semantic SHA-256 in the front matter binds every nonzero
`(stable degree,branch,source degree):coefficient` entry of this normalized
row.

## 4. Every critical vertical displacement is closed

The character-decimation mechanism of THM-3673 is independent of the choice
of collision polynomial.  For

```text
q=Q_dagger(x)+alpha t^k,       w=t,
alpha!=0,                      k>=2,                   (18)
```

it transports the quadratic coefficient blocks

```text
(J_0,J_2,J_4,J_6) -> (J_0,J_k,J_(2k),J_(3k)).         (19)
```

Rescaling the quadratic source variable gives weights

```text
R_0 -> alpha^3 R_0,
R_2 -> alpha^2 R_2,
R_4 -> alpha R_4,
Lambda(J_6) -> Lambda(J_(3k)).                        (20)
```

Hence a nonzero constant source density already contradicts

```text
0=alpha^3 kappa (326763520/3^13).                     (21)
```

Now let `0!=H in t^2 C[t]`.  THM-3675 formally monomializes `H` and transforms
a hypothetical constant density into the `x`-independent series
`kappa psi'(u)`.  Every positive stable coefficient is therefore constant in
`x`.  The value coefficients of each positive block in `(7)` sum to zero:

```text
5/18-1+13/18=0,
-4024/3159+4024/3159=0,
-116061520/9979281+116061520/9979281=0.               (22)
```

All derivative terms also vanish on constants.  Thus formal conjugacy kills
the positive blocks but leaves the nonzero constant debt `(21)`.  Therefore

```text
for every 0!=H in t^2 C[t],
no arbitrary target two-form has nonzero constant pullback.              (23)
```

## 5. Hostile controls and scope

The complete order-four predecessor window has

```text
18 rows, 105 two-form monomials, rank 15,
left-cokernel dimension 3,
zero constant-active cokernel rows.                    (24)
```

This exactly reproduces `Q_dagger`'s zero fourth debt and is a hostile control
against manufacturing a false obstruction by truncation.  The companion also
mutates one nonzero coefficient of `(7)` and detects failure on the complete
monomial universe.

THM-3680 remains a genuine survival theorem: actual target-ring coefficients
do solve `J_0=1,J_1=J_2=0`.  The present result proves that no continuation of
that branch can have constant Jacobian; it does not erase the lower-order
construction or assert where a chosen lift first fails.

The theorem closes only the fixed compiler `(1)`, `Q=Q_dagger`, and critical
vertical displacements `(23)`.  It does not cover `H'(0)!=0`, another point on
the THM-3677 parabola, higher-degree zero-debt Hermite families, nonordinary
tangent collisions, another compiler, or arbitrary planar polynomial maps.
It constructs no counterexample and proves neither direction of the planar
Jacobian conjecture, which remains **OPEN**.

## 6. Reproduction

```bash
python -B 04-computation/jc2_russell_cylinder_qdagger_sixth_stable_debt_thm3681.py
python -O -B 04-computation/jc2_russell_cylinder_qdagger_sixth_stable_debt_thm3681.py
```

Normal and optimized transcripts agree with the pinned LF output.  The
companion reconstructs all `252` columns, the two exact cokernels, the unique
constant-active row, every one of the `252` identity evaluations, the
predecessor and mutation controls, and an assertion-free AST at `287` active
gates.  **QED pending independent hostile audit.**
