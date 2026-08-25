---
id: THM-3683
title: "Russell-cylinder sixth-debt quartic on the zero-fourth parabola"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On the complete
  THM-3677 zero-fourth-debt parabola, the retained order-six arbitrary-two-form
  cokernel over Q(r) has dimension four and one active quotient class.  Its
  constant debt is -256/3^13 times an explicit irreducible quartic F(r).
  Every point with F(r) nonzero is closed for every nonzero critical vertical
  displacement.  At the four roots of F, the specialized rank remains 26 and
  the entire four-dimensional retained cokernel is constant-inactive; these
  are the only surviving degree-at-most-eight Hermite-plane candidates through
  order six.  Two roots are real and two nonreal.  No actual lift or Keller
  pair is asserted, and JC(2) remains open.
source: kps-s197 / symbolic QQ(r) retained sixth-debt computation, 2026-08-21
audit: >
  PASS -- a separate differentiation-before-truncation implementation at
  r=2/3 recovered the 30-by-252 rank 26 matrix, the exact 23-entry active
  row, the embedded three-dimensional order-four cokernel, and every
  constant-response statement.  A second dense mod-137 implementation
  checked rank 26 and zero constant debt at all four distinct reductions of
  the quartic roots.  PARI independently confirmed irreducibility,
  discriminant, real-root count, isolating intervals, and numerical roots.
  Normal, optimized, and stored transcripts are byte-identical, and the two
  semantic hashes were independently reconstructed.  The conjugate
  degree-seven THM-3651 specialization supplies an additional exact
  certificate/debt scale control.
depends_on:
  - THM-3673-russell-cylinder-monomial-ramification-debt-dilation
  - THM-3675-russell-cylinder-critical-fold-formal-conjugacy-closure
  - THM-3677-russell-cylinder-degree-eight-fourth-debt-parabola
  - THM-3681-russell-cylinder-qdagger-universal-sixth-stable-debt
related:
  - THM-3651-russell-cylinder-degree-seven-double-zero-sixth-order-closure
  - THM-4043-exceptional-quartic-shifted-stable-identities-and-j6-lift
  - THM-4046-exceptional-quartic-j7-lift-and-j8-obstruction
script: 04-computation/jc2_russell_cylinder_sixth_debt_quartic_thm3683.py
output: 05-knowledge/results/jc2_russell_cylinder_sixth_debt_quartic_thm3683.out
script_sha256: b5e4132c322b3a01883688be9e8c993c5927a38f191c547eefcc84af432d9eb3
output_sha256: 3dfc849c88ca713c1a200fe05d26676624f9f83d9432e54a1c1c17058b320025
relation_sha256: 2fa2923fb48a36623a58b7fe4d8a9c52af20ec781b0fcbb35f65af897aebd8cd
quartic_sha256: 07c0716d557f56e5db7e30e8a1ea28322a34026b1266e3939e2446937d7483a5
hash_basis: raw LF bytes for files; canonical 23-entry symbolic row and quartic text for semantic hashes
---

# THM-3683 -- only four algebraic folds survive the sixth-debt screen

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  The fixed
`Q_dagger` obstruction of THM-3681 is one point of a symbolic family.  Over
the whole zero-fourth-debt parabola, the next debt is a quartic.  Its four
zeros are the exact next candidate set.

All rings, local germs, and two-forms are over `C` unless a smaller field is
displayed.

## 1. The zero-fourth-debt family

Use THM-3677's notation

```text
P=x^2(x^2-1)^2,
R_1=P(1-x^2),                    R_2=P(4-9x),

Q_r=Q_6+p(r)R_1+rR_2,

p(r)=520r^2/9-1688r/81-5717/729.                    (1)
```

Every `Q_r` has the retained values and slopes

```text
Q_r(-1,0,1)=(-3,-3/4,-3),
Q_r'(-1,0,1)=(-9/2,1,9/2),                          (2)
```

and both earlier debts vanish:

```text
D_2(Q_r)=0,                       D_4(Q_r)=0.         (3)
```

The parameter `r=0` gives `Q_dagger`.

Use the same Russell-cylinder compiler

```text
D=1+x^2q,
b=(D-1)(D+2)^2,       c=xD(D+2),       e=q(D+3),
y=c/3,                z=e+3,

q=Q_r(x)+t^2,         w=t.                           (4)
```

For a pulled-back arbitrary target two-form write

```text
j(x,t)=sum_(n>=0)t^n J_n(x),

T_(a,d)(S)=S^(d)(a)/d!,
Lambda(S)=5T_(-1,0)(S)/18-T_(0,0)(S)+13T_(1,0)(S)/18.
                                                               (5)
```

## 2. The symbolic sixth-order identity

Define the following three linear functionals.  First,

```text
R_4(S)=
   8(2340r-503)/3159 T_(-1,0)(S)
  -8(2340r-503)/3159 T_(0,0)(S)
  -5/27 T_(-1,1)(S)+13/27 T_(1,1)(S).               (6)
```

Put

```text
A_2=315394560r^3-310946688r^2-107434431r+14507690.
```

Then

```text
R_2(S)=
  -8A_2/9979281 T_(-1,0)(S)
  +8A_2/9979281 T_(0,0)(S)
  -(177840r-26159)/28431 T_(-1,1)(S)
  -64(9r-1)/81 T_(0,1)(S)
  +13(144r+65)/2187 T_(1,1)(S)
  +10/81 T_(-1,2)(S)+26/81 T_(1,2)(S).              (7)
```

Finally put

```text
A_-=2755286876160r^4-1919428213248r^3
    -1346654271456r^2+156577328193r-9091311692,

A_0=196806205440r^4+816178042368r^3
    -347643535824r^2-119357786847r+35777404148,

B_-=7569469440r^3-11632937472r^2-2408051880r+374589757.
```

The constant block is

```text
R_0(S)=
  -16A_-/3502727631 T_(-1,0)(S)
  +16A_0/3502727631 T_(0,0)(S)
  +2B_-/89813529 T_(-1,1)(S)
  -256(58968r^2+1233r-514)/85293 T_(0,1)(S)
  +26(331776r^2+182880r-10889)/531441 T_(1,1)(S)
  +4(121680r-14087)/85293 T_(-1,2)(S)
  +52(144r+65)/6561 T_(1,2)(S)
  -20/243 T_(-1,3)(S)+52/243 T_(1,3)(S).            (8)
```

Every arbitrary target two-form satisfies the exact identity

```text
0=Lambda(J_6)+R_4(J_4)+R_2(J_2)+R_0(J_0).          (9)
```

All coefficients in `(6)--(8)` are polynomials in `r` divided only by fixed
integers.  Thus `(9)` specializes at every complex value of `r`; it is not
merely a generic rational-function statement.

## 3. Completeness of the retained window

At each retained branch, `y,z,w` are regular local parameters.  The complete
coefficient six-jet universe for an arbitrary two-form therefore consists of

```text
3 binom(9,3)=252                                    (10)
```

monomial columns.  Retain Taylor rows

```text
J_0 through source degree 3:       12 rows,
J_2 through source degree 2:        9 rows,
J_4 through source degree 1:        6 rows,
J_6 values:                         3 rows.          (11)
```

Exact elimination over `Q(r)` gives

```text
30 by 252 matrix,       rank 26,       cokernel dimension 4. (12)
```

The inherited order-four window is `18 by 105`, has rank `15`, and its
three-dimensional constant-inactive cokernel embeds in `(12)`.  Modulo that
subspace, `(9)` is the unique active class.  The displayed representative has
exactly `23` nonzero entries.

The positive value blocks satisfy

```text
sum values in Lambda=0,
sum values in R_4=0,
sum values in R_2=0.                                (13)
```

These equalities will be load-bearing under formal conjugacy.

## 4. The quartic debt and complete degree-eight closure

For `J_0=1` and `J_2=J_4=0`, equation `(9)` becomes

```text
Lambda(J_6)=-D_6(r),

D_6(r)=R_0(1)=-256F(r)/1594323,                     (14)

F(r)=72783360r^4-77822208r^3-28419741r^2
     +7849770r-1276420.                              (15)
```

At `r=0`,

```text
D_6(0)=326763520/3^13,                              (16)
```

recovering THM-3681 with the same sign.

For `q=Q_r(x)+alpha t^k`, character dilation sends the stable blocks

```text
(J_0,J_2,J_4,J_6)->(J_0,J_k,J_(2k),J_(3k))
```

with weights `(alpha^3,alpha^2,alpha,1)`.  For a general nonzero
`H in t^2 C[t]`, reuse the `Q`-independent formal-conjugacy calculation of
THM-3675: a hypothetical constant density becomes `x`-independent, all
derivative terms vanish, and `(13)` kills every positive stable value block.
The remaining invoice is

```text
0=alpha^3 kappa D_6(r).                              (17)
```

Therefore

```text
F(r)!=0
  => no arbitrary target two-form has nonzero constant pullback
     for any 0!=H in t^2 C[t].                       (18)
```

Combining `(18)` with THM-3677 closes the entire degree-at-most-eight
principal Hermite plane except the four roots of `F`: points off the parabola
die at order four, and points on it with `F(r)!=0` die at order six.

## 5. The four exceptional folds are genuine sixth-window survivors

The polynomial `F` is irreducible over `Q`.  It has exactly two real roots,
isolated by

```text
-1/2<r_-<-2/5,                 13/10<r_+<4/3,       (19)
```

and one nonreal conjugate pair.  Numerically,

```text
r_-=-0.461526763807094...,
r_+= 1.311961237837735...,
r_C= 0.109398147600064... +/-0.130365061901059... i.
```

Its discriminant is

```text
-8003379749137430490967268590424928811403550720<0. (20)
```

More is true than the vanishing of the displayed debt.  Work over the
irreducible quartic field `Q(alpha)`, `F(alpha)=0`.  The specialized complete
matrix still has

```text
rank 26,       cokernel dimension 4,                (21)
```

and every cokernel row is constant-inactive.  Since every complex root is an
embedding of this field, `(21)` holds at all four roots.  Thus no retained
order-six arbitrary-two-form identity supplies a constant debt there.

This is a survival statement only for the complete retained six-jet window.
It does not construct target-ring coefficients solving `J_0=1`, does not
show that `J_1,J_2,...` can be killed, and does not rule out a later debt.
The sharp next probe is the actual lower-lift system at `Q_alpha`, followed,
if feasible, by orders eight and ten.

## 6. Verification and scope

The companion builds all `252` columns directly over `Q(r)`, verifies `(9)`
on each column, reconstructs both symbolic cokernels and their embedding,
checks the `23` coefficients and value-block sums, proves irreducibility and
root counts by exact polynomial arithmetic, and independently recomputes the
rank and all constant responses in the quartic field.

The independent hostile audit rebuilt the matrix from full pulled-back
polynomials at `r=2/3`, differentiating before truncation, and recovered the
same rank, active row, predecessor embedding, and value-block cancellations.
Over `F_137`, where `F` has four distinct roots, a separate dense jet path
gave rank `26` and zero constant response at every root.  PARI independently
confirmed the exact quartic arithmetic and root data.  As a further algebraic
control, the THM-3651 degree-seven subline `p=0`, `r=-beta/9` gives

```text
F(-beta/9)=2187(21544632beta-97639283)/33800,
<ell_+,12>=-(3645/11)D_6,                            (22)
```

with the minus identity obtained by conjugation.  This specialization is a
cross-control, not a dependency in the proof of the symbolic family.

Reproduce with

```bash
python -B 04-computation/jc2_russell_cylinder_sixth_debt_quartic_thm3683.py
python -O -B 04-computation/jc2_russell_cylinder_sixth_debt_quartic_thm3683.py
```

The theorem concerns only this Russell compiler, the principal retained
slope packet, degree at most eight, and critical displacements
`H'(0)=0`.  Mixed folds, other tangent charts, higher-degree Hermite families,
and arbitrary planar maps remain open.  No counterexample is constructed and
the planar Jacobian conjecture remains **OPEN**.  Downstream THM-4043 shifts
the identity proved here and obtains a stagewise actual lift through `J_6` at
the four roots.  THM-4046 extends this retained ladder at order eight,
reaches `J_7`, and proves a nonzero `J_8` debt; neither theorem produces a
global pair.  **QED.**
