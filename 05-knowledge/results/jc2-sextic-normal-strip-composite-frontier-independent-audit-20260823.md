# Independent audit of the sextic composite normal-strip frontier

**Status: FINITE-EXACT INDEPENDENT AUDIT — PASS WITH CERTIFICATE
STRENGTHENING.  No formula correction.  Not a sextic normal-strip theorem;
`JC(2)` remains open.**

This audits the frozen 76-check primary packet:

- `04-computation/jc2_sextic_normal_strip_composite_frontier_scout_20260823.py`,
  SHA-256 `0f23f081efe9359699cfbbe34f1169832ad1d86b042a66710856a1c9131b6e02`;
- its frozen output, SHA-256
  `0ee27143276ee5f0883dc956545d12704a7e6d00c5899e94a559dbaa24a0099a`;
- `05-knowledge/results/jc2-sextic-normal-strip-composite-frontier-scout-20260823.md`,
  SHA-256 `1428254b29753009db8980f894fe4816c5eaa7ea9d483296ca6f478aebd5455b`.

The independent companion does not import or execute the primary script.  It
uses a different integration and elimination route:

- high coefficient rows are solved by the radial Poincare homotopy on closed
  coefficient-space one-forms;
- the `(6,4)` quotient is proved by explicit two-sided ideal identities and
  three Buchberger reductions;
- the `(6,5)` target arm is eliminated first, after which a three-variable
  grevlex quotient and multiplication matrices replace the primary
  four-variable lex/unit-ideal calculations.

## 1. Top sieve and complete formal packets

Over the characteristic-zero UFD `k[s]`, after the constant target direction
is normalized, the top equation is

```text
6wq'-jw'q=0.
```

For every irreducible prime, `6 ord(q)=j ord(w)`.  Removing the gcd gives

```text
j=1: (ord_h(w),ord_h(q))=(6,1),
j=2:                            (3,1),
j=3:                            (2,1),
j=4:                            (3,2),
j=5:                            (6,5).
```

The first three packets are killed exactly by subtracting respectively a
scalar multiple of `C^6,C^3,C^2`.  The `j=0` row is impossible: if `C=b(s)`,
the `x^5` coefficient is `6wb'`, and `b'=0` makes the whole bracket zero.
Thus the top gcd/shear sieve is complete, subject to the stated nonzero-top
and characteristic-zero `k[s]` hypotheses.

For each residual packet, the independent audit sets the next unknown source
coefficient to zero, divides its row by `jQ`, checks that the resulting
coefficient one-form is closed, and recovers its potential by

```text
Phi(x)=integral_(tau=0)^1 sum_i x_i omega_i(tau*x) d tau.
```

This reconstructs all five lower source coefficients, including every
integration constant, and matches the primary formulas term by term.  It
then recovers both conserved polynomials in each packet by the same homotopy:

- `(6,4)`: the `x^2` and `x^1` rows are exact; the constant row is
  `P db-N da`;
- `(6,5)`: the `x^3` and `x^2` rows are exact; the `x^1` row has a nonzero
  `B,N` curl and is genuinely nonclosed; the constant row is `P db-V da`.

This validates “complete depressed rational packet” as a formal
differential-algebra statement.  It does not supply descent or polynomial
integrality of the depressed coefficients.

## 2. Independent `(6,4)` Artin and jet certificates

After eliminating the target arm and setting `u=X+2,v=Y`, the three face
generators are equivalent to

```text
g1=6v^2-u^3,
g2=u^2(2u+3v),
g3=u^4.
```

The audit proves both ideal inclusions explicitly.  If `a,p,q` denote the
source-arm and two conserved generators after target-arm elimination, then

```text
g1=a,
g2=-2a-3p,
g3=-3q-4ua-24p,
p=-(2g1+g2)/3,
q=(-g3-4ug1+16g1+8g2)/3.
```

For lexicographic order `v>u`, the leading monomials are
`v^2,vu^2,u^4`.  All three S-pairs have displayed reductions to zero, so the
standard monomials are exactly

```text
1,u,u^2,u^3,v,uv.
```

The quotient therefore has length six.  Its geometric support is only
`u=v=0`, equivalently `(X,Y,Z)=(-2,0,1)`, so it is nonreduced.  The independent
expansion also verifies the primary bracket certificate

```text
F0=(-v+2u/3-4/3)g1-(2/3)g2+(2/3)g3.
```

The new primitive-jet extension also passes.  It is reconstructed directly
from the full source arm `A(g)` and the first conserved polynomial, rather
than from the primary displayed jet expressions.  In the chart where the
depression shift has a simple pole and `t=g^{-1}` is a uniformizer, regularity
forces successively

```text
D=0, v1=0, F=0, u1=0, v2=0, I=0, v3=0,
```

hence `D=F=I=0`, `u=O(t^2)`, and `v=O(t^4)`.  This conclusion uses the
characteristic-zero domain property when a square or cube coefficient
vanishes.  It does not cover a higher-order pole, later coefficients
`u2,v4,c4`, the even constants, the second conserved value, the constant
bracket, or Kummer descent.

## 3. Independent `(6,5)` fourteen-point and unit certificates

Eliminate the monic target-arm variable first:

```text
W=-1-X-Y-Z.
```

The audit computes a grevlex basis of the resulting ideal in `QQ[X,Y,Z]`.
Its initial monomials are

```text
Y^2 Z^2, X Z^3, Y Z^3, Z^4, X^3,
X^2Y, XY^2, Y^3, X^2Z, XYZ.
```

Exactly fourteen monomials avoid this initial ideal, proving that the
quotient dimension is fourteen without relying on the primary terminal
degree.  On this quotient, multiplication by `X` has a squarefree
degree-fourteen characteristic polynomial.  Its primitive coefficient list
has SHA-256

`0e3c6807f6d25f289fd3e1116a59d4aa01749e838efd6e6a901a051df7eab4ef`,

independently matching the primary lex eliminant.  Squarefreeness in a
fourteen-dimensional quotient proves fourteen distinct geometric points over
the algebraic closure; it does not claim that the points are rational.

Finally, the exact multiplication matrices of `F1` and `F0` both have
nonzero determinant.  Hence each element is a unit in the apparent quotient:

```text
(C0,F3,F2,A0,F1)=(1),
(C0,F3,F2,A0,F0)=(1).
```

The joint determinant-packet SHA-256 is
`24d2c3f3b4cdaffd86867adf708a16b1fc407ab824056a8d750baf97bb21059a`.
Thus `F1` kills all fourteen branches, and `F0` is nonzero on all fourteen.

## 4. Scope and certificate audit

Every mathematical qualifier in the primary report is compatible with the
computation:

| primary claim | audit | exact boundary |
|---|---|---|
| twelve universal buckets | pass | degree at most six in each coordinate |
| constant target `SL2` direction | pass | characteristic-zero `k[s]`; nonzero top direction |
| gcd/Kummer table and three shears | pass | UFD factorization; `R,Q` nonzero constants |
| complete `(6,4)/(6,5)` packets | pass | formal depressed rational coefficient field only |
| `(6,4)` length-six cusp | pass | principal balanced face after target-arm elimination |
| primitive `(6,4)` jet | pass | simple pole of the depression shift only |
| `(6,5)` fourteen reduced points | pass | geometric points over an algebraic closure |
| `F1` and `F0` unit conclusions | pass | principal balanced face only |
| return to original polynomials | open as stated | simultaneous arm regularity, descent, integrality |
| sextic theorem / `JC(2)` | open as stated | no global conclusion |

There is no sign, factor, support, length, or unit-ideal correction.
There are two certificate/documentation nits:

1. The primary `(6,5)` gate calls its four-row lex basis “triangular” after
   checking only the number of rows.  Row count and terminal degree alone do
   not prove quotient length.  The actual basis does have leading monomials
   `W,Z,Y,X^14`, and the independent grevlex/multiplication calculation above
   supplies the missing implication.  The claim is true; its primary gate is
   underpowered.
2. The report first uses `g` for `gcd(6,j)` and later uses `g` without an
   explicit redefinition for the moving depression shift.  The face formulas
   require a statement such as `x=y+gamma(s)` and then
   `X=B/(Q gamma^2),...`.  The displayed source-arm polynomials are obtained
   modulo the target-arm equation.  This is notation/scope, not an algebraic
   error.

## 5. Exhaustive ledger of channels not closed by the scout

Write `gamma` for the depression shift, `d=-v(gamma)>0`, and
`m_f=-v(f)` for pole orders at a DVR.  In `(6,4)` the four target-arm term
orders are

```text
L_Q=4d,
L_B=m_B+2d,
L_N=m_N+d,
L_b=m_b.
```

The principal weak cone is `m_B<=2d,m_N<=3d,m_b<=4d`; zero normalized
coordinates are included in its face algebra.  Outside that cone, target-arm
regularity forbids a unique leader.  The exact nonprincipal leader set is one
of the following four:

```text
{B,N}, {B,b}, {N,b}, {B,N,b}.
```

For a named set `S`, all `L_i` with `i in S` are equal to a common
`L>4d`, and every other lower-coefficient order is strictly smaller.  This
enumerates all nonprincipal dominant ties for a pole shift.

In `(6,5)` the term orders are

```text
L_Q=5d,
L_B=m_B+3d,
L_N=m_N+2d,
L_V=m_V+d,
L_b=m_b.
```

The principal weak cone has bounds `2d,3d,4d,5d`.  The eleven exact
nonprincipal leader sets are

```text
pairs:   {B,N}, {B,V}, {B,b}, {N,V}, {N,b}, {V,b};
triples: {B,N,V}, {B,N,b}, {B,V,b}, {N,V,b};
quadruple: {B,N,V,b}.
```

Again their common order is greater than `5d`, with all omitted leaders
strictly below it.  The primary principal-face calculation closes none of
these fifteen nonprincipal cones.

The remaining channels, not represented merely by those leader subsets, are:

- `gamma` a DVR unit, regular with positive valuation, or identically zero;
- a pole of order `d>1`, whose intermediate uniformizer jets are invisible
  to `t=gamma^{-1}`;
- strict drops inside the principal weak cone.  A zero face coordinate is
  included algebraically, but the next nonzero valuation/strict transform is
  not;
- every proper subset of the lower coefficients having poles, and every
  coefficient being regular, a unit, zero at the prime, or identically zero;
- each conserved carrier being a pole, a unit, regular with a finite zero, or
  identically zero, including cancellation of its displayed leading form and
  the next surviving jet;
- for `(6,4)`, odd Kummer valuation (ramified quadratic extension), even
  valuation with nonsquare unit (unramified), split/already-rational square,
  constant extension, and the involution/parity descent of every coefficient;
- finite primes away from the scale divisor, finite zeros of coefficients,
  constant scale, the polynomial degree fan at infinity, and the all-constant
  channel;
- simultaneous regularity of both original arms after undoing depression,
  followed by descent from the rational or Kummer coefficient field to the
  original polynomial ring.

These omissions are exactly why the correct verdict is a verified frontier
packet rather than a sextic theorem.

## 6. Reproduction

```bash
python3 -B 04-computation/jc2_sextic_normal_strip_composite_frontier_independent_audit_20260823.py
python3 -B -O 04-computation/jc2_sextic_normal_strip_composite_frontier_independent_audit_20260823.py
```

Both executions must byte-match
`05-knowledge/results/jc2_sextic_normal_strip_composite_frontier_independent_audit_20260823.out`.
