---
id: THM-3619
title: "Russell-cylinder even-fold all-order collision-jet closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every even
  polynomial Q with Q(0)=-3/4 and
  Q(1)=Q(-1)=-3, the THM-3612 fold admits no pair of regular target
  functions whose source Jacobian is a nonzero constant.  The all-order
  comparison-curve invoice forces the Taylor germ
  -3/4-9/(4x^2) at x=1, which no polynomial with the stated values can
  have.  Exact source-order 4,6,8,10 matrices remain independent finite
  controls.  No JC(2) counterexample is claimed.
source: root / audit_thm3612 higher-jet continuation, 2026-08-21
audit: >
  PASS -- two independent hostile derivations checked the local chart, exact
  moving triple, parity and error orders, vertical coefficient extraction,
  all-order induction, polynomial contradiction, and finite controls; normal,
  optimized, and stored 310-gate transcripts are byte-identical.
depends_on:
  - THM-3612-russell-cylinder-even-fold-nongraph-collision-jet-rigidity
related:
  - THM-3561-rational-keller-danielewski-polynomial-completion
  - THM-3605-russell-cylinder-graph-slice-puncture-no-filling
  - THM-3610-russell-cylinder-full-linear-projection-collision-rigidity
  - THM-3623-russell-cylinder-even-general-vertical-fold-all-order-closure
script: 04-computation/jc2_russell_cylinder_even_fold_higher_jet_staircase_thm3619.py
output: 05-knowledge/results/jc2_russell_cylinder_even_fold_higher_jet_staircase_thm3619.out
script_sha256: 5430d5078ca017ae0aab39c95c4c1beb35520e5ad1823077c767e52dab2ffc18
output_sha256: 5332bfe8d30afb3401bb5fe818fd2fb8178ea30691e4caa73a6a07ebc0955a66
hash_basis: raw LF bytes
---

# THM-3619 -- Russell-cylinder even-fold all-order collision-jet closure

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This theorem
closes the doubly tuned boundary left by
THM-3612.  The proof is local at the retained triple collision, but its
quantifier over target pairs is global and unrestricted: the two outputs may
be arbitrary regular functions on the Russell target cylinder.

All rings, formal germs, derivatives, and differential forms are over `C`.

## 0. Statement and scope

Retain the THM-3612 polynomial fold

```text
E_Q(x,t)=(x,Q(x)+t^2,t),
Q even,                 Q(0)=-3/4,       Q(1)=Q(-1)=-3.      (1)
```

For arbitrary regular functions `F,G` on the exponent-one Russell target
cylinder, write

```text
J_Q(F,G)=Jac_(x,t)(F o Psi o E_Q,G o Psi o E_Q).            (2)
```

Then

```text
                         J_Q(F,G) notin C*.                  (3)
```

Thus **every even polynomial fold `(1)` is closed** as a possible source of a
nonzero constant planar Jacobian, for every regular target pair.  This does
not close non-even folds or implicit source planes outside `(1)`.  The
nonzero critical vertical replacements `t^2 -> H(t)` with `H in t^2 C[t]`
are proved separately in THM-3623; vertical terms with `H'(0)!=0` remain
outside both theorems.

Write

```text
q_n=Q^(n)(1).                                                (4)
```

The mechanism is an all-order necessary recurrence.  THM-3612 first forces
`q_1=9/2`; thereafter a comparison of the three collision branches forces

```text
q_n=-(n+1)q_(n-1),                         n>=2.             (5)
```

The resulting formal germ has a pole at `x=0` and cannot be polynomial.

## 1. A local Darboux chart and the universal pullback formula

Before the Russell-cylinder isomorphism of THM-3605, use the THM-3561
coordinates

```text
D=1+x^2q,
a=q/D^2,                  c=xD(D+2),
b=ac^2=(D-1)(D+2)^2,      e=a(b+4)=q(D+3).                 (6)
```

At the common target point `(b,c,e,w)=(0,0,-3,0)`, the function `b+4` is a
unit.  Hence

```text
a=e/(b+4),                 c=c,                 w=w          (7)
```

are regular local target parameters on `Y_2 x A1`; transporting them through
the polynomial isomorphism of THM-3605 gives an equally valid chart on the
Russell target cylinder.  The exact source determinant is

```text
Jac_(x,q)(a,c)=-3.                                        (8)
```

For the actual two-form `dF wedge dG`--or, more generally, for an arbitrary
formal target two-form--write

```text
Omega=P da wedge dc+K da wedge dw+R dc wedge dw.           (9)
```

On the fold `q=Q(x)+t^2,w=t`, let total `x` derivatives at fixed `t` be
denoted by `a_x,c_x`.  Equations `(8),(9)` give the exact coefficient

```text
E_Q^*Omega=j(x,t) dx wedge dt,
j=-6tP+a_x K+c_x R.                                      (10)
```

Indeed, `q_t=2t`, so
`Jac_(x,t)(a,c)=Jac_(x,q)(a,c)q_t=-6t`.

At the middle branch `x=0,t=0`, evenness and `(1)` give

```text
a=-3/4,             c=0,             a_x=0,             c_x=3. (11)
```

If `(2)` were a nonzero constant, scale one output so that it is `12`.
Then `(10),(11)` force

```text
R(y_0)=4.                                                 (12)
```

No decomposability or closedness condition on `(9)` is used below.  Allowing
arbitrary `P,K,R` is a strict enlargement of the forms `dF wedge dG`, so the
eventual obstruction applies in particular to every regular target pair.

## 2. The exact moving comparison triple

Put

```text
s=t^2,       g=(1-4s/3)^(-1/2),       X=g-1,
Q_infinity(x)=-3/4-9/(4x^2).                              (13)
```

The following three points in `(x,q,w)` share an **exact** target in the
local chart `(a,c,w)`:

```text
gamma_0=(0,-3/4+s,t),
gamma_-=(-g,-3/g^2,t),          gamma_+=(g,-3/g^2,t),

(a,c,w)(gamma_0)=(a,c,w)(gamma_-)
                 =(a,c,w)(gamma_+)=(-3/4+s,0,t).          (14)
```

For the middle point `D=1`; for the two side points `D=-2`.  Also

```text
Q_infinity(g)+s=-3/g^2.                                  (15)
```

Suppose, for some `n>=2`, that the side jets of `Q` and `Q_infinity` agree
through order `n-1` at `x=1`.  Then

```text
E(t)=Q(g)-Q_infinity(g)=O((g-1)^n)=O(s^n)=O(t^(2n)).      (16)
```

The actual side-fold points at `x=+-g` therefore differ from the comparison
points `(14)` only by `E(t)` in their `q` coordinate.  Taylor expansion in
that coordinate changes the `K,R` terms of `(10)` by `O(t^(2n))`; the `P`
term has the extra factor `t` and changes by `O(t^(2n+1))`.  In the side
derivatives, the actual slopes `Q'(+-g)` are retained rather than replaced by
the slope of `Q_infinity`.

Let `j_-,j_0,j_+` be the source germs of `(10)` in the local coordinates
`xi=x+1,x,x-1`, respectively, and define

```text
C_Q(t)=j_-(-X,t)+j_+(X,t)-2j_0(0,t).                     (17)
```

At the exact common target `(14)`, parity gives

```text
a(-x,q)=a(x,q),          c(-x,q)=-c(x,q),
a_x(-g)=-a_x(g),         c_x(-g)=c_x(g).                 (18)
```

Thus the `P` terms in `(17)` cancel among the three branches, the side `K`
terms cancel each other, and the middle `K` term vanishes.  Direct
differentiation of `c=xD(D+2)` at the side comparison value
`x=g,q=-3/g^2`, while retaining the actual slope `Q'(g)`, gives

```text
c_x=12-2g^3Q'(g);                   (side)
c_x=3.                              (middle)              (19)
```

Consequently, with `a_0=-3/4+s`,

```text
C_Q(t)=R(a_0,0,t)[18-4g^3Q'(g)]+O(t^(2n))
      =-4R(a_0,0,t)g^3(Q'-Q_infinity')(g)+O(t^(2n)),     (20)
```

because `Q_infinity'(x)=9/(2x^3)`.  This identity is the geometric source of
every row of the staircase.

## 3. Coefficient extraction and the all-order invoice

Let

```text
q_n^infinity=Q_infinity^(n)(1)
            =(-1)^(n-1) 9(n+1)!/4.                       (21)
```

These jets obey

```text
q_n^infinity=-(n+1)q_(n-1)^infinity.                     (22)
```

Under the order-`n-1` matching hypothesis of Section 2, put

```text
delta_n=q_n-q_n^infinity=q_n+(n+1)q_(n-1).               (23)
```

Then

```text
(Q'-Q_infinity')(1+X)
       =delta_n X^(n-1)/(n-1)!+O(X^n),
X=(2/3)s+O(s^2).                                         (24)
```

Using `(12),(20),(24)` yields the exact leading coefficient

```text
[t^(2n-2)]C_Q(t)
 =-c_n delta_n,
c_n=2^(n+3)/(3^(n-1)(n-1)!),                 n>=2.       (25)
```

For completeness, `(25)` is also the intrinsic vertical source-jet invoice.
For `N=2n-2`, define

```text
Delta_N(j)=[t^N]j_- -2[t^N]j_0+[t^N]j_+.                (26)
```

Suppose all source coefficients of total degree below `N` equal those of the
constant `12`.  A shifted monomial `xi^k t^ell`, with `k>=1`, can contribute
to order `t^N` in `(17)` only if `ell+2k<=N`.  If a higher coefficient of
`X(t)^k` supplies the remaining order, then

```text
k+ell <= N-k < N.                                        (27)
```

Its source coefficient is therefore already zero.  Only `k=0` survives, so

```text
[t^N]C_Q(t)=Delta_N(j).
```

Combining this with `(25)` proves the all-order quotient identity

```text
Delta_(2n-2)
 =-2^(n+3)/(3^(n-1)(n-1)!)
    (q_n+(n+1)q_(n-1)),                       n>=2,       (28)
```

modulo the constant-J source rows of lower total degree.

At `n=2`, `(28)` is

```text
Delta_2=-(32/3)(q_2+3q_1),                               (29)
```

which is exactly the corrected THM-3612 invariant
`-16(2q_2+27)/3` after `q_1=9/2`.

## 4. Induction and the polynomial contradiction

Assume for contradiction that `(2)` is the normalized constant `12`.
THM-3612 first forces

```text
q_1=9/2=q_1^infinity.                                    (30)
```

Apply `(28)` with `n=2`.  Since a constant germ has every `Delta_N=0` and
`c_n!=0`, it forces `q_2=q_2^infinity`.  Inductively, if the jets match
through `n-1`, the same identity forces

```text
q_n=q_n^infinity=-(n+1)q_(n-1).                          (31)
```

Thus the Taylor series of `Q` at `1` is exactly the Taylor series of
`Q_infinity` there.  Equivalently, the polynomial

```text
H(x)=xQ'(x)+2(Q(x)+3/4)                                  (32)
```

has zero Taylor series at `x=1`, because `Q_infinity` satisfies `H=0`.
Hence `H` is the zero polynomial.  If `a_d x^d` is the highest nonconstant
term of `Q`, its coefficient in `H` is `(d+2)a_d`, impossible in
characteristic zero.  Therefore the only polynomial solution of `(32)` is
`Q=-3/4`, contradicting `Q(1)=-3`.  This proves `(3)`.

The same contradiction is visible directly from `(21)`: every
`q_n^infinity` is nonzero, whereas the derivatives of a polynomial vanish
above its degree.

## 5. Independent finite matrices and the first four higher invoices

The comparison proof above is all-order.  Independently, exact unrestricted
target-two-form matrices in the polynomial smooth target parameters

```text
(C,Y,Z)=(C,Y,S+3/4)                                      (33)
```

verify the first four higher instances.  With the common constant normalized
to `12`, they are

```text
q_1=9/2, q_2=-27/2:
Delta_4 =-(32/9)(q_3-54);                                (34)

q_3=54:
Delta_6 =-(64/81)(q_4+270);                              (35)

q_4=-270:
Delta_8 =-(32/243)(q_5-1620);                            (36)

q_5=1620:
Delta_10=-(64/3645)(q_6+11340).                          (37)
```

For an invoice at source order `N`, the lower map uses every target
coefficient of degree at most `N-1` and every source row of total degree at
most `N-1`.  Its exact ranks are

| invoice | lower map | lower rank | degree-`N` symbol rank | full rank |
|---|---:|---:|---:|---:|
| `Delta_4` | `30 x 60` | `26` | `14` | `40` |
| `Delta_6` | `63 x 168` | `57` | `20` | `77` |
| `Delta_8` | `108 x 360` | `100` | `26` | `126` |
| `Delta_10` | `165 x 660` | `155` | `32` | `187` |

The homogeneous symbol has rank `3(N+1)-1`; its unique missing vertical
direction is `(26)`.  In each finite matrix, an explicit sparse left quotient
row expresses `(26)` in the lower rowspace.  Applying that row to the
constant target column gives, respectively,

```text
1536,              -6144,              51200/3,
-360448/9.                                                (38)
```

These are exact finite controls of `(34)--(37)`, not numerical evidence for
the all-order argument.

## 6. Hostile polynomial controls

The THM-3612 second-order survivor is

```text
Q_dag=-3/4-27x^2/2+18x^4-27x^6/4.                       (39)
```

It matches `q_1,q_2` but has `q_3=-378`; hence `(34)` gives
`Delta_4=1536`.  The following even polynomials successively match the
formal side jets farther:

```text
Q_3=-3/4-27x^2/2+9x^4+81x^6/4-27x^8+9x^10,             (40)

Q_4=-3/4-27x^2/2-45x^4/4+405x^6/4-297x^8/2
    +90x^10-81x^12/4,                                   (41)

Q_5=-3/4-27x^2/2-45x^4+270x^6-486x^8
    +855x^10/2-189x^12+135x^14/4,                       (42)

Q_6=-3/4-27x^2/2-189x^4/2+567x^6-2457x^8/2
    +2835x^10/2-1863x^12/2+1323x^14/4-99x^16/2.         (43)
```

Their first mismatches are

| control | matched formal side jets | first mismatching jet | invoice |
|---|---|---:|---:|
| `Q_dag` | through `q_2` | `q_3=-378` | `Delta_4=1536` |
| `Q_3` | through `q_3` | `q_4=7506` | `Delta_6=-6144` |
| `Q_4` | through `q_4` | `q_5=-127980` | `Delta_8=51200/3` |
| `Q_5` | through `q_5` | `q_6=2269620` | `Delta_10=-360448/9` |
| `Q_6` | through `q_6` | `q_7=-43454880` | `Delta_12=2293760/27` |

The last row uses the proved all-order invoice `(28)`, not an order-twelve
matrix.  Thus `Q_6` is no longer an open survivor: like every other
polynomial in `(1)`, it is closed at a finite order.  None of these controls
is asserted to carry a Darboux pair.

They are generated by

```text
Q_3=Q_dag+9x^4(x^2-1)^3,
Q_4=Q_3-(81/4)x^4(x^2-1)^4,
Q_5=Q_4+(135/4)x^4(x^2-1)^5,
Q_6=Q_5-(99/2)x^4(x^2-1)^6.                             (44)
```

The factor `x^4(x^2-1)^m` preserves the central value and all side jets below
order `m`; its order-`m` derivative at `1` is `2^m m!`.  The paired probes
`x^2(x^2-1)^m` change central jets while producing the same new side jet.
They independently check that `(34)--(37)` depend on the displayed side jet,
not on an omitted central jet.

## 7. Boundary and exact reproduction

The theorem proves exactly the all-regular-target-pair exclusion `(3)` for
the even quadratic-stable-coordinate polynomial folds `(1)`.  The proof does
not classify:

- non-even `Q`;
- vertical replacements with `H'(0)!=0`; or
- implicit source planes outside this fold family.

Those exits remain **OPEN**.  Nonzero replacements `H in t^2 C[t]` are now
closed by the separate all-order extension THM-3623.  The former tuned
polynomial locus is not open: the recurrence here closes it.

The deterministic companion verifies the compiler identities, the local
`b+4` chart and Jacobian, the exact common comparison target, the side
`c_x` formula, parity, the coefficient and source-order gates for
`2<=n<=12`, the rational ODE, every polynomial control `(39)--(44)`, and the
four independent unrestricted-two-form matrices `(34)--(38)`.

Reproduce with

```bash
python3 04-computation/jc2_russell_cylinder_even_fold_higher_jet_staircase_thm3619.py
python3 -O 04-computation/jc2_russell_cylinder_even_fold_higher_jet_staircase_thm3619.py
```

Both streams must be byte-identical to
`05-knowledge/results/jc2_russell_cylinder_even_fold_higher_jet_staircase_thm3619.out`.
