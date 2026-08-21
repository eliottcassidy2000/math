---
id: THM-3598
title: "Danielewski rational-exact polar graph family and classification"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  On every exponent-n Danielewski surface with n>=2 there is an explicit
  infinite family of polynomial target functions with rational
  constant-bracket mates.  Choosing the transverse factor a=Sigma makes the
  target function critical on every arm, and the mate has the matching arm
  poles.  In every exponent, positive graph-linear functions h(b)+A(b)c
  with rational mates are classified by A^(n-1)/h' in C[h] of degree at
  most n-2; exponent three collapses to h'=A^2/beta.  None of the
  nonconstant-h exact graphs has a regular mate on a nonconstant squarefree
  target.  Complementarily, arbitrary polynomial positive-channel repairs above
  A/c^s have no rational mate when s divides n-1 and A contains a squarefree
  multiarm divisor.
  This is a polar near-counterexample, not a Darboux pair or a JC(2)
  counterexample.
source: root / planar-Jacobian mixed-coordinate construction-hostile session, 2026-08-21
audit: >
  An independent hostile audit rederived the all-exponent residue classifier,
  explicit primitive, exponent-three parity collapse, constant-h boundary,
  arm/residual primitive-torsor obstruction, shallow-negative field trace,
  and reciprocal divisor-degree lemma.  Normal and optimized companions are
  byte-identical to the stored 599-gate output.
depends_on: []
related:
  - THM-3589-danielewski-central-arm-every-line-and-kummer-trace-darboux-gates
  - THM-3595-danielewski-separated-transverse-time-form-rational-nonentry
  - THM-3596-a13-invoice-paid-mixed-coordinate-toric-time-form-nonentry
script: 04-computation/jc2_danielewski_rational_exact_polar_graph_thm3598.py
output: 05-knowledge/results/jc2_danielewski_rational_exact_polar_graph_thm3598.out
script_sha256: 95b8275c950d32ce862558cc3e17cf9f277495c93bbb3d817a85d58df0762a52
output_sha256: 22ad72a4d963548805480d8a3c0b2a08ae3a1eca75d3ef38d3f731a78fe04799
hash_basis: raw LF bytes
---

# THM-3598 -- Danielewski rational-exact polar graph family and classification

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This theorem is
a positive hostile to the idea that rational
time-form exactness is rare or nearly sufficient.  It gives many exact
rational mates, then identifies precisely why they cannot be regular.  It
constructs no polynomial Darboux pair and proves no case of `JC(2)`.

## 1. An all-exponent rational-exact family

Work over `C`.  Let `n>=2`, let `Sigma in C[b]` be nonzero, and put

```text
Y_(n,Sigma): c^n e=Sigma(b),                            (1)

{b,c}=c^n,       {c,e}=-Sigma'(b),
{b,e}=-n c^(n-1)e.                                      (2)
```

Choose any nonzero `a in C[b]`, choose a polynomial `H` with

```text
H'=a^(n-1),                                              (3)
```

and put `x=a(b)c`.  For every `G in C[T]`, define

```text
P=1/[(n-1)x^(n-1)],             Q=H(b)+G(x).            (4)
```

Then `P in C(Y)`, `Q in C[Y]`, and

```text
{P,Q}=1.                                                 (5)
```

Indeed, on the `c!=0` chart the bracket is

```text
{R,S}=c^n(R_bS_c-R_cS_b).                               (6)
```

Since `G(x)` commutes with `P(x)`, while

```text
{x,H}=-aH'c^n=-a^n c^n=-x^n,
dP/dx=-x^-n,                                             (7)
```

their product is one.  Thus the Hamiltonian time form of `Q` has the exact
rational primitive `(4)`.

### 1.1 The all-arm specialization pays with poles

Now assume `Sigma` is squarefree and nonconstant, and choose

```text
a=Sigma,          H'=Sigma^(n-1).                       (8)
```

Then

```text
x=Sigma c=c^(n+1)e,                                     (9)
```

so `Q=H+G(c^(n+1)e)` is a polynomial function.  This alternate presentation
records its contact with the arm divisor; intrinsically it is still the
positive-weight, all-arm-critical graph `H+G(Sigma c)`, not a third
independent channel.
On every arm

```text
D_beta={b=beta,c=0},             Sigma(beta)=0,         (10)
```

one has

```text
dH=Sigma^(n-1)db=0,
dx=d(c^(n+1)e)=0,
dQ=0.                                                    (11)
```

Consequently no **regular** `R in C[Y]` can satisfy `{R,Q}=1`: evaluating
the bracket on `D_beta` gives zero.  The mate `(4)` escapes exactly by having
a pole along `x=0`.  This is the arm-pole debt predicted by the central-arm
geometry, now visible without any support count.

For the A13 target `n=3`, `Sigma=b(b^2+1)`, a cheapest explicit row is

```text
H=b^7/7+2b^5/5+b^3/3,
x=Sigma c=c^4e,
P=1/(2x^2),                     Q=H+x+x^2.              (12)
```

It is a nonlinear polynomial target function with an exact rational mate,
but it is critical on all three arms.

## 2. Complete positive graph-linear classification in every exponent

The family above is not an isolated trick.  On `Y_(n,Sigma)`, put `k=n-1`
and let

```text
Q=h(b)+A(b)c,                  h nonconstant, A nonzero. (13)
```

Define `D_h=(1/h')d/db` and `R=A^k/h'`.  Then `Q` has a rational
constant-bracket mate if and only if

```text
R in C[h] and deg_h R<=k-1,       equivalently D_h^k R=0. (14)
```

Writing `R^[ell]` for ordinary differentiation of the one-variable
polynomial `R`, an exact mate is

```text
P=sum_(ell=0)^(deg R)
    (-1)^ell R^[ell](Q)/[ell!(k-ell)(Ac)^(k-ell)].      (15)
```

### 2.1 The residue classifier

Let `w` be transcendental.  The generic fibre of `(13)` has rational
function field `C(w)(b)`, because

```text
c=(w-h)/A,       eta=A^k db/(w-h)^(k+1).               (16)
```

At a generic simple root `b_0` of `h(b)=w`, the residue is

```text
res_(b_0) eta=(-1)^(k+1)(D_h^k R)(b_0)/k!.             (17)
```

A rational differential on `P1` is exact exactly when all residues vanish.
The roots `b_0` sweep a Zariski-dense set as `w` varies, so `(16)` is exact
if and only if `D_h^kR=0`.  The constants of `D_h` in `C(b)` are `C`.
Repeatedly subtracting the leading constant multiple of
`h^(k-1)/(k-1)!` therefore proves

```text
ker(D_h^k)=span_C{1,h,...,h^(k-1)}.                    (18)
```

This is `(14)`.  Conversely substitute `t=h` in `(16)` and integrate the
finite Taylor expansion of `R(w-(w-t))`; the result is exactly `(15)`.
There are no other finite poles, and the residue at infinity is the negative
sum of the finite residues.  This proves necessity and sufficiency.

For `n=3`, `k=2`.  Write `R=alpha h+beta`.  The alternative `alpha!=0`
would make

```text
A^2=(alpha h+beta)h'
```

have odd degree `2deg h-1`, impossible for the square on the left.  Thus
`R=beta in C*`, and the classification collapses to

```text
h'=A^2/beta,                  P=beta/[2(Ac)^2].         (19)
```

For `n>=4`, nonconstant `R` genuinely occurs.  For instance, at `n=4`,
`h=b^2`, `A=2^(1/3)b`, and `R=h` satisfy `(14)`.

If `h` is constant, the conclusion is different and elementary: for
`F'=A^k`, the function `F(b)/(Ac)^n` is a rational mate.  Thus the
nonconstant hypothesis in `(13)` is essential.

### 2.2 None of the nonconstant-h exact positive graphs has a regular mate

Now assume `Sigma` is squarefree and nonconstant, and suppose `(14)` holds.
Although `Q` has the rational mate `(15)`, it has no mate in `C[Y]`.

At an arm `D_beta`, if `A(beta)=0`, then `dQ` has no `c`- or `e`-component,
so every regular bracket with `Q` vanishes there.  Suppose instead that
`A(beta)!=0`.  Put `lambda=h(beta)` and `x=Ac`.  At the generic point of
`D_beta`, `x` is a uniformizer and

```text
Q-lambda=x+O(x^n).                                      (20)
```

If `R` has order `m<=k-1` at `lambda`, the leading order of `(15)` is
`x^(m-k)`.  Its coefficient is nonzero because

```text
sum_(ell=0)^m (-1)^ell binom(m,ell)/(k-ell)
 =(-1)^m m!(k-m-1)!/k! !=0.                            (21)
```

Thus the displayed mate has an arm pole.  Every other mate differs from it
by `r(Q)`, `r in C(Q)`: on the rational generic fibre, the constants of the
Hamiltonian derivation are exactly `C(Q)`.  If `r` is regular at `lambda`, it
cannot cancel the arm pole.  If `r` has a pole there, it creates a pole on
the residual divisor

```text
h(b)-lambda+A(b)c=0,             c!=0,                 (22)
```

where `(15)` is regular.  Hence the primitive torsor contains no global
regular function.  This is a vertical-saturation obstruction, not a finite
weight count.

### 2.3 Arbitrary positive repairs cannot rescue a shallow negative channel

The sign of the transverse weight is load-bearing.  Suppose `Sigma` is
squarefree with at least two distinct roots.  Let `1<=s<=n-1` divide `n-1`,
write `n-1=qs`, take arbitrary `F in C[b,c]`, and let `A` be a nonzero
multiple of `Sigma`.  Then

```text
Q_-=F(b,c)+A(b)/c^s
   =F(b,c)+(A/Sigma)c^(n-s)e                           (23)
```

is regular on `Y_(n,Sigma)`.  Nevertheless there is no
`P in C(Y_(n,Sigma))` with `{P,Q_-}=1`.

Put `K_0=C(w)(b)` and

```text
T(c)=c^s(F-w)+A.                                        (24)
```

This is the separable minimal polynomial of `c` for the rational map `Q_-`.
At a root `c_i`, `T'(c_i)=c_i^s(Q_-)_c(c_i)`.  The forced time form is
`db/[c^n(Q_-)_c]`, so partial fractions give the exact trace identity

```text
Tr_(K_0(c)/K_0)(eta)
 =sum_i db/[c_i^(n-s)T'(c_i)]
 =-[c^(n-s-1)](1/T(c)) db.                              (25)
```

The last coefficient is a polynomial in `w` of degree `q-1`, with leading
term `-w^(q-1)/A^q`.  Exactness upstairs would imply exactness of the trace;
differentiating `q-1` times in `w` would make `db/A^q` exact.

That is impossible when `A` has at least two distinct roots.  Indeed, if
`S'=1/A^q` and the distinct roots of `A` have multiplicities
`m_1,...,m_r`, then (apart from the immediate simple-pole logarithmic case)
`S` has finite poles of orders `qm_i-1`, while `S-S(infinity)` has a zero of
order `q deg(A)-1` at infinity.  Degree zero of a principal divisor would
require

```text
sum_i(qm_i-1)=q deg(A)-r >= q deg(A)-1,
```

so `r<=1`.  The root-count boundary is sharp: for one-root `A`, `db/A^q`
is rationally exact whenever `q deg(A)>1`; the remaining equality case is
logarithmic.  In exponent three, `(25)` specializes to

```text
s=1: (F(b,0)-w)db/A^2,          s=2: -db/A,            (26)
```

closing arbitrary nonnegative-weight repairs above either shallow negative
channel.  The first negative depth not covered on the A13 target is `-3`,
the honest `e` channel.

## 3. Assumption challenges and sharp boundaries

1. Rational exactness is not close to polynomial regularity.  In `(12)` the
   time form is exact, yet the target function is critical along every arm.
2. The exponent condition is sharp for this mechanism.  At `n=1`, the
   primitive of `-x^-1 dx` is logarithmic, not rational.
3. The square collapse `(19)` is special to exponent three.  Higher
   exponents are still completely classified by `(14)`, but admit
   nonconstant `R` and additional degree patterns.
4. Choosing `a` without all arm factors may move or remove some critical
   arms, but `(4)` still has poles along `ac=0`.  Section 2.2 proves the
   no-regular conclusion for graph-linear targets; arbitrary nonlinear `G`
   remains outside that torsor argument.
5. Section 2.3 allows arbitrary nonnegative powers of `c`, but only one
   shallow negative channel `c^-s` with `s|(n-1)`.  Deeper negative tails,
   several negative weights, and nonlinear dependence not of this form lie
   outside its trace obstruction.

```text
source       positive exact graph and shallow negative-channel repair
target       generic Q-fibre and its rational time form
map          primitive torsor for the positive graph; field trace for negative
preserved    rational exactness and the highest coefficient in fibre value w
destroyed    regularity or sheetwise pole locations before trace
sidecar      arm/residual vertical divisors and reciprocal-primitive degree
cheap test   compare the arm pole with its residual fibre; extract top w term
```

The correct construction lesson is therefore not merely “make the time form
exact.”  One must make it exact with a primitive regular simultaneously on
every arm, or arrange a genuine cross-arm pole cancellation unavailable in
the one-channel family.

## 4. Exact verification contract

The companion checks:

- `(5)--(7)` for `2<=n<=8`, several nontrivial `a`, and nonlinear `G`;
- the all-arm identity `(9)`, the critical-arm jet `(11)`, and row `(12)`;
- the exponent-three instance of residue identity `(17)` through an
  independent jet expansion;
- the all-exponent criterion `(14)`, primitive `(15)`, and nonconstant-`R`
  witnesses for `4<=n<=9`;
- exact-square exponent-three families `(19)` and active non-square hostiles;
- the degree-parity exclusion of `alpha!=0` on a broad integer box;
- the constant-`h` branch, arm leading coefficient `(21)`, and residual-fibre
  saturation;
- trace coefficient `(25)` for shallow negative channels dividing `n-1`;
- the multi-root reciprocal-primitive obstruction and its one-root boundary;
- the `n=1` logarithmic boundary.

The identities are universal; finite rows are exact controls, not
extrapolation.

Reproduce with

```bash
python3 04-computation/jc2_danielewski_rational_exact_polar_graph_thm3598.py
python3 -O 04-computation/jc2_danielewski_rational_exact_polar_graph_thm3598.py
```
