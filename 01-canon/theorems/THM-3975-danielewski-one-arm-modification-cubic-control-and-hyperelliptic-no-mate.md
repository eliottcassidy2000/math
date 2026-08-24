---
id: THM-3975
title: "Danielewski one-arm gradings, cubic controls, and hyperelliptic no-mates"
status: >
  PROVISIONAL VERIFIED-EXACT ALL-HEIGHT FINITE-TOWER EXTENSION / DELTA
  HOSTILE AUDIT REQUIRED, atop a PROVED + VERIFIED-EXACT + INDEPENDENTLY
  HOSTILE-AUDITED core. For every n>=2 the one-arm determinantal completion
  B_n has an elementary two-color DPD grading and the marked homogeneous LND
  x partial_t has kernel k[x] and plinth x^(n+1)k[x]. At n=2,3 the pair
  (p,x+y) is finite free of degree three with normal basis (1,x,z-1) and
  exact monogenic index p. The provisional extension proves that for every
  ell>=floor(n/2), (p,y+x^ell) is finite of exact degree 2ell+1. Uniformly
  in n, p has no rational constant-Jacobian mate, so none of these finite
  controls is Keller. No planar Jacobian counterexample is claimed.
source: jc-extra-debt-local / post-THM-3973 one-arm structural supplement, 2026-08-24
audit: >
  PASS (jc-zero-debt-lift, 2026-08-24). The audit independently reconstructed
  the two-chart normalization and affine modification, both DPD denominator
  colors, marked LND kernel and exact plinth valuation, the n=2 and n=3
  multiplication tables, finite-free normal bases, monogenic index and trace
  discriminants, the three generic P=0 addresses, and the cancellation
  div(dP wedge dC)=E with no boundary ramification term. It also rederived
  the conic residues at n=2 and the holomorphic hyperelliptic differential
  obstruction at every n>=3. Normal and optimized runs byte-match the frozen
  CHECKS=246 output; both raw hashes and the semantic hash agree. This audit
  predates and does not certify the provisional all-height finite-control
  extension or its additional exact gates.
depends_on:
  - THM-3757-pell-chebyshev-three-charge-hyperelliptic-obstruction-tower
related:
  - THM-3572-squarefree-danielewski-affine-modification-and-two-bracket-collapse
  - THM-3973-exact-volume-simple-cubic-determinantal-affine-plane-completion
  - THM-3974-height-tower-few-weight-darboux-support-obstruction
script: 04-computation/jc2_danielewski_one_arm_cubic_control_thm3975.py
output: 05-knowledge/results/jc2_danielewski_one_arm_cubic_control_thm3975.out
script_sha256: 712065d6bdb0b3aa28d46b78f91d6c76566531090e04ae97c0b92f994b4241ef
output_sha256: 9c8a05afce0541fffe6b62fc48b5dcc1dc50103ff719ebbf73c335412dbd361a
semantic_sha256: c4341f7091dadd755b23509b15656029fc90861822080545bd57aa8dd7213289
hash_basis: raw LF bytes
---

# THM-3975 -- two colors, two cubic heights, and one hyperelliptic debt

**PROVISIONAL VERIFIED-EXACT ALL-HEIGHT FINITE-TOWER EXTENSION / DELTA
HOSTILE AUDIT REQUIRED, atop a PROVED + VERIFIED-EXACT + INDEPENDENTLY
HOSTILE-AUDITED core.** Work over an algebraically closed field `k` of
characteristic zero. For an integer `n>=2`, put inside `k[x,t]`

```text
z=1+x^n t,                 p=zt,
y=x^(n-1)zt^2,             B_n=k[x,z,p,y],
X_n=Spec(B_n).                                             (1)
```

THM-3973 develops the completion and its support obstructions. This
supplement isolates four additional mechanisms and proves them
self-containedly.

1. The grading is an elementary two-color DPD presentation. If the base
   coordinate is `z`, then

   ```text
   D_+=0,
   D_-=-[z=0]/(n+1)-[z=1]/n.                             (2)
   ```

   The homogeneous locally nilpotent derivation `delta=x partial_t` has

   ```text
   ker(delta)=k[x],
   delta(B_n) intersect ker(delta)=x^(n+1)k[x].           (3)
   ```

   The second equality is the plinth of this **marked** action; no
   uniqueness assertion about all additive-group actions is made.

2. At the two first heights `n=2,3`, the natural target

   ```text
   P=p,                    C=x+y                           (4)
   ```

   makes `B_n` finite free of degree three over `R=k[P,C]`, with basis
   `{1,x,w}`, `w=z-1`. The full normal-order discriminants are given in
   Section 5. The smaller monogenic order `R[x]` has determinant index
   `P`, so its cubic polynomial discriminant has the extra factor `P^2`.
   Over generic `P=0`, three distinct normal addresses survive; the
   double root seen by `x` alone is an index collision, not ramification.

   More generally, for every `n>=2` and every
   `ell>=floor(n/2)`, the pair

   ```text
   P=p,                    C_ell=y+x^ell                 (4a)
   ```

   is finite of exact degree `2ell+1`. These maps are never Keller by
   item 3.

3. For every `n>=2`, there is no `Q in k(x,t)` for which

   ```text
   J_(x,t)(p,Q) in k(p)^*.                                (5)
   ```

   In particular neither `p` nor any nonconstant polynomial `f(p)` has a
   rational constant-Jacobian mate. This is the generic-fibre differential
   mechanism of THM-3757, now in the binomial family

   ```text
   W^2=1+4p x^n.                                         (6)
   ```

The two cubic controls have an explicit nonempty relative ramification
divisor meeting the affine-plane open. The all-height controls are ruled out
more uniformly by the no-mate theorem for their first coordinate `p`. Thus
the result supplies sharp positive finiteness controls and a uniform
negative mate theorem, not a counterexample to the planar Jacobian
conjecture.

## 1. Determinantal geometry and the one-arm boundary

The three identities

```text
x^n p=z(z-1),
xy=p(z-1),
zy=x^(n-1)p^2                                             (7)
```

are the maximal minors of

```text
[[z,p,x],[x^(n-1)p,y,z-1]].                              (8)
```

They give the exact presentation

```text
B_n = k[x,z,p,y]/(the three relations (7)).               (9)
```

Equivalently, begin with the smooth Danielewski surface

```text
Y_n=Spec k[x,z,p]/(x^n p-z(z-1)).                         (9a)
```

Then

```text
B_n=Gamma(Y_n,O)[p(z-1)/x].                               (9b)
```

Thus `X_n -> Y_n` is the affine modification along `x=0` with center
`V(x,p(z-1))`. Of the two lines over `x=0`, it deletes the `z=0` arm away
from `p=0`, replaces that center point by the line `D`, and retains the
`z=1` arm. This is the precise sense in which `(1)` is a one-arm
modification.

There is no hidden saturation. Indeed, on `D(z)` the quotient is

```text
k[x,z,z^-1,p]/(x^n p-z(z-1)),                            (10)
```

with `y=x^(n-1)p^2/z`. On `D(z-1)`, eliminate
`p=xy/(z-1)` to get

```text
k[x,z,y,(z-1)^-1]/(x^(n+1)y-z(z-1)^2).                  (11)
```

Both rings are domains and both presentations agree with the subring in
`k[x,t]`. Since `z` and `z-1` generate the unit ideal, `(10)--(11)` prove
`(9)` globally.

Both charts are smooth. For `(10)`, simultaneous vanishing of the
derivatives would require `x=0` and `2z-1=0`, while the equation at `x=0`
requires `z(z-1)=0`. For `(11)`, the derivative in `z` is
`-(z-1)(3z-1)`; at a possible point with `x=0`, the equation and the
inversion of `z-1` force `z=0`, where that derivative is nonzero. Hence
`X_n` is smooth and normal.

The source plane is exactly

```text
U=D_Xn(x) union D_Xn(z) isomorphic to A2_(x,t).           (12)
```

Indeed, `t=(z-1)/x^n` on the first chart and `t=p/z` on the
second, while `x` and `z=1+x^n t` have no common source zero. Its
complement is

```text
D=V(x,z,p) isomorphic to A1_y.                            (13)
```

Since `U^*=k^*` and `Cl(U)=0`, the localization sequence for the sole
boundary prime gives

```text
B_n^*=k^*,                   Cl(X_n)=Z[D].                (14)
```

For completeness, there can be no nonzero relation on `[D]`: a rational
function whose divisor is supported on `D` restricts to a unit on
`U=A2`, hence is scalar and has zero divisor.

The source form has the uniform, height-independent divisor

```text
eta=dx wedge dt,                 div_Xn(eta)=D.           (15)
```

To see the order at `D`, differentiate `z=1+x^nt` and use the local
equation `(11)`:

```text
eta=x^-n dx wedge dz
   =x/(unit) dx wedge dy             near generic D.     (16)
```

It also has a global primitive

```text
beta=-dz/((n-1)x^(n-1)),              d beta=eta.         (17)
```

Regularity is not a formal cancellation in the fraction field. From the
first relation in `(7)`,

```text
(2z-1)dz=n x^(n-1)p dx+x^n dp.                           (18)
```

At every point of `V(x)`, one has `z=0` or `z=1`, so `2z-1` is a unit.
Equation `(18)` therefore proves that `(17)` is regular across both
components over `x=0`; away from them, `x` is a unit. Finally

```text
d[-x^(1-n)dz/(n-1)]=x^-n dx wedge dz=eta.                (19)
```

This structural ledger will be used below to distinguish a Jacobian pole
from genuine relative ramification.

## 2. The exact two-color grading

Give the fraction field the grading

```text
wt(x)=1,             wt(z)=0,             wt(t)=-n.      (20)
```

Then

```text
p=x^-n z(z-1),
y=x^(-n-1)z(z-1)^2.                                    (21)
```

For every `m>=0`, the nonnegative piece is

```text
(B_n)_m=x^m k[z].                                        (22)
```

For `m>=1`, the negative piece is exactly

```text
(B_n)_(-m)
 =x^-m z^ceil(m/(n+1))(z-1)^ceil(m/n) k[z].              (23)
```

This formula is elementary. A generator monomial of weight `-m` can be
written

```text
x^c p^a y^b f(z),
c=na+(n+1)b-m>=0.                                        (24)
```

After extracting `x^-m`, its coefficient is

```text
z^(a+b)(z-1)^(a+2b) f(z).                                (25)
```

Among the pairs `a,b>=0` satisfying the inequality in `(24)`,

```text
min(a+b)=ceil(m/(n+1)),
min(a+2b)=ceil(m/n).                                     (26)
```

For the first minimum, a fixed value `a+b=l` has weight at most
`(n+1)l`, and equality is attained with `a=0`. For the second, a fixed
value `a+2b=l` has weight at most `nl`, and equality is attained with
`b=0`. Since `z` and `z-1` are coprime and `k[z]` is a PID, the ideal
generated by all coefficients `(25)` is exactly the product in `(23)`.

Equations `(22)--(23)` are equivalently the elementary DPD row `(2)`.
They make the two colors asymmetric but complementary:

```text
z=0 has denominator n+1,          z=1 has denominator n. (27)
```

No external classification of hyperbolic actions is used here; `(2)` is
simply a concise encoding of the proved homogeneous pieces.

## 3. The marked additive action and its exact plinth

The derivation

```text
delta=x partial_t                                             (28)
```

preserves `B_n`, because

```text
delta(x)=0,
delta(z)=x^(n+1),
delta(p)=x(2z-1),
delta(y)=(z-1)(3z-1).                                    (29)
```

It is locally nilpotent as the restriction of `x partial_t` on
`k[x,t]`, homogeneous of degree `n+1`, and

```text
ker(delta)=B_n intersect k[x]=k[x].                       (30)
```

Its plinth ideal is exactly the second row in `(3)`. Suppose
`delta(g)=f(x)` for `g in B_n`. Inside `k[x,t]`, integration in `t`
gives

```text
g=(f(x)/x)t+h(x).                                        (31)
```

Polynomiality first forces `x|f`. Along `D`, equation `(1)` gives

```text
ord_D(x)=1,                  ord_D(t)=-n.                 (32)
```

The regular term `h(x)` cannot cancel a negative valuation of the first
term in `(31)`. Thus `x^n | f/x`, or `x^(n+1)|f`. Conversely,

```text
delta(b(x)z)=b(x)x^(n+1)                                  (33)
```

for every `b(x) in k[x]`. This proves `(3)`.

The exponent `n+1` is therefore an exact invariant of the marked action.
We do not infer from it that the abstract surfaces `X_n` are pairwise
nonisomorphic; that would require control of all locally nilpotent
derivations, not only `(28)`.

## 4. Why heights two and three are finite cubic

Set `P=p` and `C=x+y`. The relations `(7)` imply the universal identity

```text
Py+xy^2=x^(n-1)P^3.                                      (34)
```

Indeed, the left side is `y(P+xy)=yPz=P(zy)`. Substituting
`y=C-x` gives

```text
F_n(X)=X(C-X)^2+P(C-X)-X^(n-1)P^3.                      (35)
```

At the first two heights these are the monic cubics

```text
F_2=X^3-2CX^2+(C^2-P-P^3)X+PC,                          (36)
F_3=X^3-(2C+P^3)X^2+(C^2-P)X+PC.                        (37)
```

Both are irreducible in `k[P,C][X]`. A reducible monic cubic over this
normal UFD would have a root in `k[P,C]`, and specialization at `P=1`
would give a polynomial root of

```text
G_2=X(C-X)^2+C-2X,
G_3=X(C-X)^2+C-X-X^2.                                   (38)
```

A root `r(C)` of degree at least two is impossible: the `r^3` term has
strictly greater degree than every other term. If `r=b` is constant, the
coefficient of `C^2` first forces `b=0`, after which the coefficient of
`C` is nonzero. Write the remaining possibility as `r=aC+b`. Its `C^3`
coefficient is

```text
a(a-1)^2.                                                (39)
```

The case `a=0` is the constant case. At `a=1`, the first row of `(38)`
becomes

```text
(b^2-1)C+(b^3-2b),                                      (40)
```

which cannot vanish, while the second row retains coefficient `-1` on
`C^2`. This proves irreducibility.

Put `w=z-1` and let `R=k[P,C]`. For `n=2`, multiplication on the proposed
`R`-basis `{1,x,w}` is

```text
x^2=Cx-Pw,
xw=Cw+C-(1+P^2)x,
w^2=PCx-(1+P^2)w.                                       (41)
```

For `n=3`, it is

```text
x^2=Cx-Pw,
xw=(C+P^3)w+C-(1+CP^2)x,
w^2=-P^2C+(PC^2+P^2+CP^4)x
    -(1+2CP^2+P^5)w.                                    (42)
```

These identities follow directly from `(7)` and `y=C-x`; conversely they
verify all three relations `(7)`. Initially let `R_0` denote the image of
the abstract polynomial ring `k[P,C]` in `B_n`. The `R_0`-span of
`{1,x,w}` is a subalgebra containing `x,z=1+w,P,y=C-x`, hence it is all of
`B_n`. In particular `B_n` is integral over `R_0`. Since `B_n` is a
two-dimensional domain, `dim(R_0)=2`; the prime kernel of
`k[P,C] -> R_0` therefore has height zero and is zero. Thus `P,C` are
algebraically independent and `R_0=R=k[P,C]` as asserted.

Over `k(P,C)`, one has

```text
w=(Cx-x^2)/P.                                            (43)
```

Irreducibility of `(36)--(37)` makes `{1,x,x^2}` linearly independent,
so `(43)` makes `{1,x,w}` linearly independent as well. Consequently

```text
B_n=R direct-sum Rx direct-sum Rw,       n=2,3.          (44)
```

The map `(P,C):X_n -> A2` is finite flat of degree three. Since `B_n` is
normal, it is the full integral closure of `R` in its cubic function
field.

For `ell=1`, the same target stops being monic after height three: at
`n=4` the leading coefficient in `(35)` is `1-P^3`, and for `n>=5` the
term `-P^3X^(n-1)` dominates. The correct all-height repair is to move the
positive-weight summand rather than abandon `P`.

### 4.1. An all-height finite tower of exact odd degree

Fix an integer

```text
ell>=floor(n/2),                    C_ell=y+x^ell.       (44a)
```

Substitution of `y=C_ell-x^ell` in `(34)` gives

```text
F_(n,ell)(X)
 =X(C_ell-X^ell)^2+P(C_ell-X^ell)-P^3X^(n-1).           (44b)
```

Because

```text
2ell+1>n-1,                                               (44c)
```

this is monic of degree `2ell+1` in `X`. Hence `x` is integral over
`k[P,C_ell]`, then `y=C_ell-x^ell` is integral, and finally `z` is
integral by

```text
z^2-z-x^nP=0.                                             (44d)
```

Thus `B_n` is finite over the image of `k[P,C_ell]`. The same dimension
argument used above shows that `P,C_ell` are algebraically independent.

The degree is not merely bounded by `2ell+1`; it is exact. Work over
`K=k(P)` and put `W=2z-1`. On the smooth generic fibre,

```text
W^2=1+4Px^n,
y=P(W-1)/(2x).                                            (44e)
```

Conversely,

```text
W=1+2x(C_ell-x^ell)/P,                                   (44e')
```

so `x` generates the generic function field over `K(C_ell)`.

After extending scalars algebraically, the two points over `x=0` have
`W=+1` and `W=-1`. At the plus point, `y` is regular; at the minus point,
`y` has one simple pole. If `n` is even, there are two infinity points,
`x` has pole order one at each, and `y` has pole order `n/2-1`. The
inequality `ell>=n/2` makes `x^ell` dominate, so `C_ell` has pole order
`ell` at each infinity point. If `n` is odd, there is one infinity point,
where

```text
ord_infinity(x)=-2,             ord_infinity(y)=-(n-2). (44f)
```

Now `2ell>=n-1`, so `x^ell` again dominates and contributes pole order
`2ell`. In both parities,

```text
deg div_poles(C_ell)=1+2ell.                             (44g)
```

The degree of a nonconstant function on a smooth projective curve is the
degree of its polar divisor. Therefore

```text
[k(P)(x,W):k(P)(C_ell)]=2ell+1.                          (44h)
```

This proves that `(44b)` is the generic minimal polynomial, that

```text
(P,C_ell):X_n -> A2
```

is finite of exact degree `2ell+1`, and that the threshold in `(44a)` is
the clean monicity-dominance boundary for this construction. At `n=2,3` and
`ell=1`, it recovers the two finite cubics above. Section 7 proves more:
because the first coordinate is `P=p`, no member of the all-height finite
tower can be Keller.

## 5. Index debt, normal discriminants, and the three addresses

The change from the normal basis `{1,x,w}` to the monogenic rows
`{1,x,x^2}` uses

```text
x^2=Cx-Pw                                                    (45)
```

and has determinant `-P`. In fact

```text
B_n/R[x] isomorphic to R/(P),              n=2,3,        (46)
```

so `(P)` is the exact index ideal. Let `Delta_n` be the trace
discriminant of the normal basis in `(44)`. The multiplication tables give

```text
Delta_2=
 4C^4P-8C^2P^4+20C^2P^2+C^2
 +4P^7+12P^5+12P^3+4P,                                  (47)

Delta_3=
 4C^5P+C^4P^4+22C^3P^2+22C^2P^5+C^2
 +4CP^8+22CP^3+P^6+4P.                                  (48)
```

Accordingly,

```text
Disc_X(F_n)=P^2 Delta_n,                 n=2,3,          (49)
```

where the left side is the polynomial discriminant, equivalently the
trace discriminant of `{1,x,x^2}`.

The geometric meaning of the factor `P^2` is exact. On `P=0`, both tables
specialize to

```text
x^2=Cx,              xw=Cw+C-x,              w^2=-w.    (50)
```

For generic `C!=0`, the three reduced points are

```text
(x,w)=(0,-1),             (C,-1),             (C,0).    (51)
```

Moreover,

```text
Delta_2(0,C)=Delta_3(0,C)=C^2!=0.                        (52)
```

Thus the normal finite algebra is etale at all three generic addresses.
The monogenic polynomial specializes to `X(X-C)^2` and forgets the
distinction between the last two points in `(51)`. The apparent double
branch is exactly the index debt `(46)`, not a ramification divisor of the
normal surface.

## 6. The genuine ramification left by the cubic controls

For `n=2,3`, direct differentiation gives

```text
x J_(x,t)(P,C)=-R_n,
R_n=x(2z-1)+y((n-2)z+1).                                (53)
```

There are two reduced components over `x=0`:

```text
D=V(x,z,p),
L_1=V(x,z-1,y),
div_Xn(x)=D+L_1.                                        (54)
```

The numerator does not contain `D`, since

```text
R_n restricted to D = y.                                (55)
```

It contains `L_1` exactly once. On the source plane,

```text
R_n=x[2z-1+x^(n-2)zt^2((n-2)z+1)].                      (56)
```

The bracket in `(56)` restricts generically to `1+t^2` on `L_1` for
`n=2`, and to `1` for `n=3`. Define the effective residual divisor `E_n`
by

```text
div_Xn(R_n)=L_1+E_n.                                    (57)
```

Equations `(53)--(57)` and `(15)` now give the precise cancellation

```text
div_Xn(J_(x,t)(P,C))=E_n-D,
div_Xn(dP wedge dC)=E_n.                                 (58)
```

In particular there is **no** ramification component along the added
boundary `D`: the simple pole of the source-coordinate Jacobian cancels
the simple zero of `dx wedge dt`. This uniform order-one cancellation is
easy to lose if one extrapolates from the weight `n` of `t` rather than
working in the smooth boundary chart.

The residual divisor is nonempty and meets `U=A2`. Indeed, on `U`

```text
J_2=-(t^3x^2+t^2+2tx^2+1),
J_3=-(t^4x^7+3t^3x^4+2t^2x+2tx^3+1),                   (59)
```

and each is a nonconstant polynomial over the algebraically closed field
`k`. Therefore `(P,C)` is finite but not etale and cannot be a Keller
pair. From `(58)` one also recovers the canonical class identity

```text
[E_n]=[D] in Cl(X_n).                                    (60)
```

No irreducibility claim for `E_3` is needed here. At `n=2`, THM-3973
additionally proves irreducibility of the residual curve.

## 7. The hyperelliptic generic fibre forbids every mate of p

Let `K=k(P)` with `P` transcendental, and restrict to the generic fibre
`p=P`. Put

```text
W=2z-1=1+2x^n t.                                        (61)
```

The equation `P=t(1+x^nt)` gives

```text
K(X_(n,P))=K(x,W),               W^2=1+4P x^n.          (62)
```

The polynomial on the right is squarefree, so `(62)` is the actual smooth
generic function field. Now suppose more generally that

```text
J_(x,t)(p,Q)=h(p),          Q in k(x,t),
h(p) in k(p)^*.                                             (63)
```

On the generic fibre, `p_t=W`. Since `dp=0`, the chain rule gives

```text
dQ=-h(P) dx/W.                                           (64)
```

Thus `(63)` would make the nonzero differential `dx/W` rationally exact
on the smooth projective model of `(62)`.

At `n=2`, after extending scalars by `sqrt(P)`, the conic has two points
at infinity. With local parameter `s=1/x`, one has

```text
W/x -> +/-2sqrt(P),
dx/W = -/+ [1/(2sqrt(P))] ds/s + regular terms.          (65)
```

Both residues are nonzero. A rational exact differential has zero residue
at every place, contradicting `(64)`.

For `n>=3`, the form `dx/W` is a nonzero holomorphic differential on the
smooth projective hyperelliptic curve `(62)`. At a finite branch point,
`W` is a parameter and `dx/W` is regular. At infinity its order is
nonnegative: for odd `n=2g+1` the order is `2g-2`, and for even
`n=2g+2` it is `g-1`. A rational function whose differential is
holomorphic has no poles in characteristic zero, hence is constant; its
differential is zero. This again contradicts `(64)`.

This proves the stronger statement `(5)`. If `f in k[s]` is nonconstant
and a rational `Q` satisfied

```text
J_(x,t)(f(p),Q)=c in k^*,                                (66)
```

then

```text
J_(x,t)(p,Q)=c/f'(p) in k(p)^*,                          (67)
```

contradicting `(5)`. Hence every nonconstant polynomial reparametrization
of the distinguished coordinate `p` has the same no-mate obstruction.

There is nevertheless a cheap near-slice:

```text
J_(x,t)(p,x(1-2z))=1+2(n+2)x^n p.                        (68)
```

It solves the constant term but leaves a supported multiple of `p`.
Equations `(62)--(65)` explain why no sequence of rational corrections can
remove that remainder globally: the missing primitive is logarithmic at
height two and holomorphic of positive genus thereafter.

## 8. Preservation, loss, and scope

The four mechanisms record different information about the same object.

```text
two-color grading:
  preserves  exact allowed zeros at z=0 and z=1 by weight;
  loses      addition between distinct weights;

finite cubic control:
  preserves  the full normal order and its three generic addresses;
  loses      etaleness along E_n;

all-height finite tower:
  preserves  the distinguished first coordinate p and exact map degree;
  loses      Keller etaleness by the generic p-fibre differential;

generic p-fibre:
  preserves  every rational candidate mate of p;
  loses      possible Darboux pairs whose first coordinate is unrelated
             to p.                                       (69)
```

The theorem therefore does **not** assert that no Darboux pair exists on
`B_n`, that all locally nilpotent derivations are conjugate to `(28)`, or
that the surfaces `X_n` are pairwise nonisomorphic. The broader support
and bounded-search hostiles belong to THM-3973 and THM-3974. The new sharp
content here is the marked plinth, the second finite-free cubic height,
the normal-versus-monogenic address separation at both heights, the exact
odd-degree finite tower `(44a)--(44h)`, and the all-height rational no-mate
theorem for `p`.

## 9. Exact verification contract

The companion script checks `414` exact identities and hostile gates. It
verifies:

1. all three minors, both chart equations, the primitive, LND formulas,
   near-slice, and the two ceiling exponents at heights `2<=n<=7`;
2. the two cubic equations and their specialized linear-root gates;
3. both full multiplication tables, commutation, all three determinantal
   relations, and cubic annihilation;
4. the trace discriminants `(47)--(49)`, the three addresses `(51)`, and
   the Jacobian numerator `(53)`;
5. for sampled `2<=n<=7` and four consecutive exponents beginning at
   `floor(n/2)`, the all-height elimination `(44b)`, monic degree
   `2ell+1`, generic-fibre coordinate `(44e)`, parity pole arithmetic, and
   nonconstant Jacobian.

The normal and optimized runs agree byte for byte. The semantic hash
commits to the uniform proofs beyond the sampled height range; the samples
are regression controls, not finite evidence for the all-height claims.

Reproduce with

```bash
python3 04-computation/jc2_danielewski_one_arm_cubic_control_thm3975.py
python3 -O 04-computation/jc2_danielewski_one_arm_cubic_control_thm3975.py
```
