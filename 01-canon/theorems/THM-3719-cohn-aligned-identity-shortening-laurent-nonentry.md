---
id: THM-3719
title: "Cohn complete identity-shortening Laurent nonentry"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Over every
  characteristic-zero field, no
  alternating two-factor left decoration of Cohn's matrix followed by zero
  or one arbitrary elementary right factor is a Jacobian matrix.  Two scalar
  determinant PDEs reduce by y-degree/x-degree descent to forbidden lowest
  Laurent layers in t=xy coordinates.  Together with THM-3709, this closes
  the complete alternating two-left/two-right Cohn cell for arbitrary right
  parameters, including identity factors.  It does not treat longer reduced
  words, other non-elementary classes, or JC(2).
source: root + kps-s191 / 2026-08-22
audit: >
  PASS.  An independent derivation checked all four determinant PDEs, both
  degree/valuation descents, the two forbidden Laurent layers, the
  quarter-turn signs, and completeness with THM-3709.  A hostile script audit
  replaced a tautological duality gate by both explicit substitutions.
  Normal, optimized, and frozen transcripts then agreed exactly.
depends_on: []
related:
  - THM-3653-cohn-factorial-repair-and-weighted-rectangle-holonomy
  - THM-3655-cohn-alternating-two-source-factor-row-obstruction
  - THM-3709-cohn-alternating-two-by-two-elementary-decoration-nonentry
script: 04-computation/jc2_cohn_aligned_identity_shortening_laurent_thm3719.py
output: 05-knowledge/results/jc2_cohn_aligned_identity_shortening_laurent_thm3719.out
script_sha256: dc55a042cebd9c5ff00ff275933c7f731c9dfc67afe771135a2b641f7b359c6c
output_sha256: 8476ea7f2b4fc0a40c0d3ddf70cef2e59429cc361cba4cbf0af3f8d651ed5a87
semantic_sha256: f971dee92c30b91f204c21b57171c8357d1fae9c545a54ec200cdeb3a2a05e94
hash_basis: raw LF bytes
---

# THM-3719 -- every identity-shortened two-left Cohn word is nonintegrable

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  THM-3709 closes
the alternating two-by-two decoration whenever both right parameters are
nonzero.  This theorem closes the entire boundary where at least one right
factor is the identity.

Work over a characteristic-zero field `k`.  Put

```text
A=1+xy,       B=x^2,       G=-y^2,       D=1-xy,
C=[A B;G D],                                      det C=1,
E_+(h)=[1 h;0 1],                 E_-(h)=[1 0;h 1].    (1)
```

For arbitrary `f,g,u,w in k[x,y]`, none of

```text
E_+(f)E_-(g) C E_-(u),          E_-(g)E_+(f) C E_-(u),
E_+(f)E_-(g) C E_+(w),          E_-(g)E_+(f) C E_+(w)  (2)
```

is a Jacobian matrix.  Taking `u=0` or `w=0` includes the completely absent
right word.

## 1. Expose one row and integrate it

If a row `(P,R)` is closed, characteristic zero supplies a polynomial
potential `Q` with `(P,R)=(Q_x,Q_y)`.  The two alternating left words expose
respectively

```text
bottom row rho_2+g rho_1,          top row rho_1+f rho_2. (3)
```

Because every matrix in `(2)` has determinant one, an exposed closed row and
its complementary row must satisfy a scalar determinant PDE.  We show that
neither of the two PDE types has a polynomial solution.  The other two cells
then follow by a quarter-turn duality.

## 2. The aligned `E_-` cell

Let `N=C E_-(u)`.  Its top row is

```text
rho_1=(1+xy+x^2u,x^2).                                  (4)
```

If the bottom row in the first matrix of `(2)` were `dQ`, then
`det(rho_1,dQ)=1` would give

```text
(1+xy+x^2u)Q_y-x^2Q_x=1.                               (5)
```

We prove `(5)` impossible.  Put `d=deg_y Q`.  The case `d=0` would say
`-x^2Q_x=1`, so `d>=1`.  Write `T=1+xy+x^2u`.  Its `y`-degree is at least one:
the coefficient of `y` contains `x(1+xu_1)` and cannot vanish polynomially.
If `deg_y T>=2`, the top term of `TQ_y` has degree strictly above every term
of `x^2Q_x`.  Hence

```text
deg_y T=1,                 u=u_0(x)+u_1(x)y.            (6)
```

Let `q_d(x)` be the leading `y`-coefficient of `Q`.  The coefficient of
`y^d` in `(5)` is

```text
x q_d'=d(1+xu_1)q_d.                                    (7)
```

The `x`-adic valuation of `(7)` is `ord_x(q_d)=d`.  Thus
`q_d=x^dR`, with `R(0)!=0`, and cancellation in `(7)` gives

```text
R'=d u_1R.                                               (8)
```

A nonzero polynomial cannot have a polynomial nonzero logarithmic derivative:
degrees in `(8)` force `u_1=0`, then `R` is constant.  In particular
`u=u(x)`.

Now set `t=xy` and write

```text
Q_tilde(x,t)=Q(x,t/x) in k[x,x^(-1),t].                 (9)
```

After division by `x`, equation `(5)` becomes

```text
[(1+x^2u(x))partial_t-x partial_x]Q_tilde=x^(-1).       (10)
```

Let `r` be the least Laurent exponent of `x` in `Q_tilde`, and let `F(t)` be
its nonzero coefficient.  The `x^2u partial_t` term raises the exponent by at
least two, while the remaining operator preserves it.  If `r<-1`, the lowest
layer of `(10)` would say

```text
F'-rF=0,                                                (11)
```

which has no nonzero polynomial solution for negative `r`.  If `r>-1`, the
right side of `(10)` cannot occur.  Thus `r=-1`, where `(10)` forces

```text
F'+F=1,                    hence F=1.                   (12)
```

But an honest monomial `x^i y^j` contributes to Laurent exponent `-1` only
when `j=i+1>=1`; its contribution is `t^j`.  Therefore the actual `F(t)` is
divisible by `t`, contradicting `(12)`.  The first cell in `(2)` is impossible.

## 3. The cross `E_+` cell

Let `N=C E_+(w)`.  Its top row is

```text
rho_1=(A,S),                  S=x^2+Aw.                 (13)
```

If the exposed bottom row in the third matrix of `(2)` were `dQ`, then

```text
A Q_y-S Q_x=1.                                         (14)
```

Put `d=deg_xQ>=1` and `s=deg_xS`.  Directly `s>=2`.  If
`deg_xw>=2`, then `s=deg_xw+1>2`; comparison of the top `x`-degrees in `(14)`
is impossible.  Hence `s=2` and

```text
w=w_0(y)+xw_1(y),             [x^2]S=1+yw_1.           (15)
```

For the leading `x`-coefficient `q_d(y)` of `Q`, the coefficient of
`x^(d+1)` in `(14)` gives

```text
y q_d'=d(1+yw_1)q_d.                                    (16)
```

The same valuation and logarithmic-derivative argument as `(7)--(8)` gives

```text
q_d=c y^d,                    w=w(y).                   (17)
```

Use coordinates `(y,t=xy)` and put `Q_tilde(y,t)=Q(t/y,y)`.  Equation `(14)`
becomes

```text
[(1+t)partial_y+
 (t/y-(1+t)y w(y))partial_t]Q_tilde=1.                 (18)
```

If `r` is the least Laurent exponent of `y` and `F(t)` its coefficient, the
`w`-term raises `r`, while the remaining terms produce at exponent `r-1`

```text
tF'+r(1+t)F.                                            (19)
```

For `r<0`, the homogeneous equation in `(19)` has no nonzero polynomial
solution: its `rtF` term has degree one above all others.  At `r=0`, `(19)`
forces `F` constant; subtract that harmless constant from `Q` and repeat.
Thus a normalized solution has `r>=1`.  Values `r>1` cannot produce the
constant right side of `(18)`, so `r=1`, where

```text
tF'+(1+t)F=1.                                           (20)
```

For `H=tF`, this is `H'+H=1`, whose only polynomial solution is `H=1`.
That contradicts `t|H`.  The third cell in `(2)` is impossible.

## 4. Quarter-turn duality and completeness

Set

```text
(X,Y)=(y,-x),              P(X,Y)=-Q(-Y,X).            (21)
```

The exposed-top scalar equation for `C E_+(w)` is

```text
(1-xy-y^2w)Q_x+y^2Q_y=1.                               (22)
```

Under `(21)`, it becomes `(5)` with
`u(X,Y)=-w(-Y,X)`.  Likewise the exposed-top equation for `C E_-(u)` becomes
the cross equation `(14)` with `w(X,Y)=-u(-Y,X)`.  Sections 2--3 therefore
exclude the second and fourth cells of `(2)` as well.

Every identity-shortened alternating two-left/two-right decoration is one of
the four cells `(2)`.  THM-3709 treats the complementary case in which both
right parameters are nonzero.  Their union proves that for arbitrary
`u,w,f,g in k[x,y]`, neither left order applied to either alternating
two-factor right order can be a Jacobian matrix.

## 5. Reproduction and boundary

Run

```bash
python3 -B 04-computation/jc2_cohn_aligned_identity_shortening_laurent_thm3719.py
python3 -B -O 04-computation/jc2_cohn_aligned_identity_shortening_laurent_thm3719.py
```

The assertion-free companion verifies all four row determinants and scalar
PDEs, both Laurent conjugacies on independent parameter/degree grids, the
lowest-layer kernels through exponent `-12`, and 98 nonlinear hostile row
probes.  The universal proof is the degree/valuation/Laurent argument above.

This closes a constructive double-coset cell, not every elementary word, every
non-elementary `SL_2/E_2` core, a Keller pair, or JC(2).  The first live Cohn
grammar now requires a longer reduced word on at least one side.  **QED.**
