# Flat reciprocal index localization and the polynomial Keller firewall

**Source audit + PROVED local mechanism + FINITE-NUMERICAL global probes +
PROVED all-`Q` quadratic-`u` conductor-square no-go.**  This note does **not** prove
the announced global Caratheodory/Loewner counterexample, and it does not
construct a planar Keller counterexample.  It extracts what can be proved from
the displayed formula, identifies the still-unchecked global assertions, and
turns the mechanism into exact degree, divisor, symmetry, and search gates for
JC2.

Companion:

```text
python3 04-computation/caratheodory_flat_reciprocal_keller_transfer_audit.py
python3 -O 04-computation/caratheodory_flat_reciprocal_keller_transfer_audit.py
```

Both executions are required to reproduce
`05-knowledge/results/caratheodory_flat_reciprocal_keller_transfer_audit.out`
byte-for-byte.

## 1. Source correction: the requested arXiv identifier is unrelated

The exact public announcement is

<https://x.com/__alpoge__/status/2089971359921156203>

by Levent Alpoge, posted August 19, 2026.  The post credits John-Paul Smith and
Claude, says that a PDF writeup had circulated privately, and then records the
formula publicly.  The smallest relevant quote scope is that it claims an
umbilic of **“index `1+k/2`”**, and for `g_2` **“exactly one umbilic point.”**

The requested identifier
<https://arxiv.org/abs/2608.19068> is instead Simon Brendle and Pei-Ken Hung,
*A metric on `S^2 x S^2` with positive sectional curvature*.  Its primary PDF
contains no Alpoge--Smith Caratheodory/Loewner construction.  Exact-title,
exact-formula, author, and arXiv searches found no public PDF of the announced
construction as of August 21, 2026.  Thus this is an audit of the **primary
public formula**, not of the unavailable private writeup.  The distinction is
load-bearing.

The displayed data are

```text
f(x+iy) = -cos(2x)/4 + 3cos(2y)/10 - cos(4y)/32 + sin(x)sin(y),

g_k(z) = |z|^2 exp(-|z|^(-1/4) exp(-|z|^2))
         f((100/bar(z))^(k/2))/(1+|z|^2) + 10^10.       (1)
```

Here `k` is a positive integer.  Because `f(-w)=f(w)`, the sign monodromy of
the half-power for odd `k` disappears on the punctured plane.  Only `g_2` is
asserted in the announcement to be a smooth spherical support function with a
unique umbilic.

## 2. Exact seed mechanism: a uniformly nonzero spin-2 field

Write

```text
D=f_xx-f_yy,
C=2f_xy.
```

Direct differentiation gives the **PROVED** identities

```text
D=cos(2x)+(6/5)cos(2y)-(1/2)cos(4y),
C=2cos(x)cos(y).                                        (2)
```

Put `a=cos^2 x`, `b=cos^2 y`.  Then

```text
D=2a-4b^2+(32/5)b-27/10,
D^2+C^2=D^2+4ab=:M(a,b),       0<=a,b<=1.              (3)
```

There is an elementary exact minimization.  For fixed `b`, the unconstrained
minimizer in `a` is

```text
a_*(b)=(40b^2-74b+27)/20.
```

Consequently the constrained minimizer is `a=1` on `[0,1/10]`, `a=a_*` on
`[1/10,1/2]`, and `a=0` on `[1/2,1]`.  On the last piece,

```text
M(0,b)-49/2500
 =4(5b-4)^2(100b^2-160b+71)/625,
100b^2-160b+71=100(b-4/5)^2+7.                         (4)
```

On the middle piece,

```text
M(a_*,b)=b(5b-3)(8b-9)/5,
d/db M(a_*,b)=3(4b-1)(10b-9)/5,                        (5)
```

whose candidate values at `b=1/10,1/4,1/2` are
`41/100,49/80,1/4`.  On the first piece, `4ab=4b` already gives the desired
bound for `b>=49/10000`; below that cutoff `D(1,b)` is increasing but remains
less than `D(1,49/10000)<0`, whose square is larger than `49/2500`.
Equality in `(4)` is attained at `a=0,b=4/5`.  Therefore

```text
min_(R^2) (D^2+C^2)=49/2500.                            (6)
```

With the convention

```text
f_(bar w bar w)=(D+iC)/4,
```

equation `(6)` says

```text
|f_(bar w bar w)| >= 7/200                              (7)
```

everywhere.  This is the seed's core: it is a bounded periodic spin-2 field
with a quantitative nonvanishing margin.  Since it is nonzero on the entire
simply connected `w`-plane, its phase has winding zero on every closed plane
loop.  No numerical sampling enters `(2)--(7)`.

Two cheaper nonvanishing checks also expose why the coefficients were chosen.
If `C=0` and `cos x=0`, then

```text
D=-(40b^2-64b+27)/10,
```

whose quadratic has negative discriminant.  If `cos y=0`, then
`D=(20a-27)/10`, also nonzero for `0<=a<=1`.  Thus the simultaneous zero set
is empty before one even asks for the sharp margin.

## 3. The flat reciprocal envelope and exact local index multiplication

Let

```text
m=k/2,        c=100^m,
rho(r)=r^2 exp(-r^(-1/4)exp(-r^2))/(1+r^2),
w=c bar(z)^(-m).                                        (8)
```

For the round metric in stereographic `z`-coordinates, the relevant
trace-free spherical Hessian component is

```text
Q(h)=h_zz + 2bar(z)h_z/(1+|z|^2).                      (9)
```

An umbilic is a zero of `Q`.  If a positive-eigenline is oriented locally,
its unoriented line index around a zero is

```text
line-index = -wind(Q)/2.                               (10)
```

The leading term obtained by differentiating the reciprocal argument twice is

```text
Q(g_k)=rho(r)m^2c^2 z^(-2m-2) f_(bar w bar w)(w)+R_k.  (11)
```

All derivatives of `f` are bounded.  Moreover

```text
rho_z/rho=O(r^(-5/4)),       rho_zz/rho=O(r^(-5/2)).   (12)
```

The terms in `R_k` come from one reciprocal derivative, envelope derivatives,
and the spherical connection.  Relative to the magnitude
`rho r^(-2m-2)` in `(11)`, they are uniformly

```text
O(r^m)+O(r^(m-1/4))+O(r^(2m-1/2))+O(r^(m+2)).          (13)
```

Because `m>=1/2`, `(13)` tends to zero.  The exact seed margin `(7)` therefore
makes `(11)` a nonvanishing homotopy to its leading term on every sufficiently
small circle.

The positive scalar factors in `(11)` contribute no phase.  The factor
`z^(-2m-2)=z^(-k-2)` contributes winding `-(k+2)`.  The seed factor contributes
zero winding.  For odd `k`, the path in the `w`-plane ends at `-w(0)` rather
than `w(0)`, but `f_(bar w bar w)` is even.  More intrinsically, its global
logarithm is even because it is a nonzero even map on the simply connected
plane.  Hence its phase still closes with winding zero.  We obtain the
**PROVED local statement**

```text
wind Q(g_k)=-(k+2),
index_0(g_k)=1+k/2.                                    (14)
```

This proof is local near zero.  It does not show that `Q(g_2)` is nonzero at
every other point of the sphere.

### Why the envelope is smooth and why polynomiality destroys it

At zero, `exp(-r^(-1/4)+o(1))` beats every inverse power of `r`.  Every
derivative of the reciprocal oscillation grows only as some finite inverse
power, so their product is `C^infty`-flat.  This makes all Taylor coefficients
at zero vanish while the punctured neighborhoods retain arbitrary prescribed
phase winding.

At infinity, use `z=1/zeta`.  For `k=2`,

```text
r^2/(1+r^2)=1/(1+|zeta|^2),
r^(-1/4)e^(-r^2)=|zeta|^(1/4)e^(-1/|zeta|^2),
100/bar(z)=100bar(zeta).                               (15)
```

The middle expression is again smooth-flat.  Thus `(15)` proves the claimed
smooth extension of `g_2` across infinity.  It is important not to extend this
argument carelessly to odd `k`: sign invariance resolves monodromy on the
punctured plane, but a half-power at infinity need not be smooth merely because
the outer function is even.

The polynomial firewall is exact.  Along a positive ray the seed includes

```text
f(w)=-cos(2w)/4+43/160.
```

For every finite `L` and `m>0`, repeated differentiation of
`r^L cos(cr^(-m))` has a leading term of size

```text
const * r^(L-n(m+1))
```

after `n` derivatives.  Choosing `n>L/(m+1)` and suitable phase subsequences
shows that it is not `C^infty` at zero.  Hence no nonzero finite-order
real-analytic radial multiplier, and in particular no polynomial envelope,
can replace `rho` in this reciprocal construction.  The flat function is not
cosmetic damping; it is exactly the device that hides infinite oscillation
from every finite jet.

## 4. What the finite spherical computation does and does not establish

For `k=2`, the companion differentiates `(1)` symbolically and samples

```text
4Q=(g_xx-g_yy)-2i g_xy
   +4bar(z)(g_x-i g_y)/(1+|z|^2).                      (16)
```

At `131072=2^17` equally spaced points on each radius

```text
0.3, 0.5, 1, 2, 5, 10,
```

it finds winding `-4` and strictly positive sampled minima.  Thus the sampled
line index is `2`, independently reproducing the reported computation.  This
is **FINITE-NUMERICAL**: it certifies neither values between mesh points nor
the entire punctured sphere.

The spherical connection term is load-bearing.  If one mistakenly uses only
the Euclidean spinor

```text
(g_xx-g_yy)-2i g_xy,
```

the sampled winding is `-4` on radii through `2` but `-6` on radii `5` and
`10`.  This hostile replay explains why a numerically plausible calculation
with the wrong Hessian convention can give the wrong global index ledger.

The script also checks a positive control `exp(-4i theta)` and rejects a
zero-crossing hostile control `cos(theta)`.  Ordinary and optimized Python
replays use the same explicit `require` gates.

The announcement's two global claims remain at these honest statuses:

- `Q(g_2)` has no other spherical zero: **ANNOUNCED + FINITE-NUMERICAL probes,
  not proved here**.
- `10^10` makes the radius tensor positive: **ANNOUNCED, not independently
  certified here**.

Those are precisely the two missing pieces needed to promote the formula to a
fully audited smooth Caratheodory counterexample.

## 5. Support-function convexity: what is general and what is specific

For a smooth support function `h` on `S^2`, the radius tensor is

```text
R_h = nabla^2_(S^2) h + h g_(S^2).                     (17)
```

Its trace-free part is exactly the trace-free Hessian, so umbilics are its
zeros.  Adding a constant `C` gives the **PROVED general identity**

```text
R_(h+C)=R_h+Cg,                                        (18)
```

while leaving the trace-free part and every umbilic index unchanged.  By
compactness, every smooth `h` admits some sufficiently large `C` for which
`R_(h+C)` is positive definite.  This is the support-function convexification
mechanism.

Identity `(18)` proves existence of a large shift once smoothness is known.  It
does not by itself prove that the displayed numerical value `10^10` exceeds
the global negative eigenvalue bound of the oscillatory perturbation.  That
specific estimate belongs to the unavailable writeup or requires a separate
interval proof.

There is no useful Keller analogue of `(18)`.  In the trace gauge

```text
F=(lambda x+H_y, lambda y-H_x),
det DF=lambda^2+det Hess(H),                            (19)
```

adding a constant to `H` does literally nothing to `F`.  Adding constants to
the two target components is only target translation and cannot repair
properness or injectivity.  Adding a large scalar identity part changes the
constant-Hessian equation in `(19)`.  Thus the `10^10` trick should not be
imported into a Keller ansatz.

## 6. Polynomial index and divisor ledgers

The honest algebraic surrogate for reciprocal index multiplication is a
Laurent-support calculation.

Let `q(x,y)` be a complex polynomial of total degree at most `D`, nonzero on
the circle of radius `R`.  With `zeta=e^(i theta)`, substitute

```text
x=R(zeta+zeta^(-1))/2,
y=R(zeta-zeta^(-1))/(2i).                              (20)
```

Then `q` restricts to a Laurent polynomial supported in exponents
`[-D,D]`.  Multiplying by `zeta^D` and applying the argument principle gives

```text
wind(q)=N-D,       0<=N<=2D,
so |wind(q)|<=D.                                       (21)
```

More sharply, if the actual Laurent support lies in `[ell,u]`, then

```text
ell <= wind(q) <= u.                                   (22)
```

For a polynomial Hessian carrier `H` of degree `d`, either spinor

```text
q_+=H_xx-H_yy+2iH_xy,
q_-=H_xx-H_yy-2iH_xy                                  (23)
```

has degree at most `d-2`.  Therefore any isolated real line-field index on a
zero-free circle satisfies the **PROVED degree gate**

```text
|index|<= (d-2)/2.                                     (24)
```

Imitating the smooth value `1+k/2` costs at least

```text
d>=k+4.                                                (25)
```

This does not exclude a polynomial imitation; it prices it visibly in the
Newton polygon.  There is no flat Taylor-invisible index reservoir.

The constant-Hessian equation adds an exact factor ledger:

```text
q_+q_-=(Delta H)^2-4 det Hess(H).                      (26)
```

If `det Hess(H)=mu` is constant, then over a field containing `sqrt(mu)`,

```text
q_+q_-=(Delta H-2sqrt(mu))(Delta H+2sqrt(mu)).          (27)
```

Any proposed high-index spinor zero must therefore be compatible with the
two Laplacian level divisors and the Newton support of `(27)`.  Over the reals,
if `mu<0`, `(26)` is `|q_+|^2=(Delta H)^2-4mu>0`, so there are no real umbilic
zeros at all and every circle winding is zero.  If `mu>0`, zeros can occur only
on `Delta H=+/-2sqrt(mu)`.  These are search filters, not a proof of complex
JC2: a complex Keller pair need not have a useful real slice.

### The dephasing warning

The user's quantum-walk/resistor dictionary supplies a precise caution here.
Dephasing replaces hopping amplitude by its squared modulus; analogously,
`q_+q_-` is `|q_+|^2` on a real slice.  This preserves the zero set and modulus
margin but destroys the phase winding and its sign.  A search that records
only `(26)` has “classicalized” the spinor.  To use the index gate it must keep
the phase/Laurent-exponent sidecar `(20)--(22)`.  The product identity is a
powerful necessary filter, never an index certificate by itself.

## 7. Smooth one-point localization versus the algebraic `(35,1)` extremal

THM-3586 proves, in the reduced degree cell `(72,108)`, that the common
degree-`36` top base has at least two projective roots.  This uses Shastri's
**CITED** one-place-at-infinity criterion.  At the first cap `38`, the extremal
permitted base is

```text
K=t^35(t+lambda u),        lambda!=0.                  (28)
```

Thus the multiplicity partition `(35,1)` is the maximal algebraic
concentration possible for that degree: one may put `35` units at one infinity
root, but polynomiality and the global injectivity criterion force at least
one unit at a second root.

The smooth construction localizes arbitrarily large index at one point by
putting all reciprocal oscillation behind an envelope with zero Taylor series.
Equation `(28)` is the algebraic boundary analogue: the finite divisor cannot
disappear from the jet/leading-form ledger, and the most concentrated legal
degree-`36` divisor still splits into two roots.  This is a mechanism analogy,
not an equality between umbilic index and root multiplicity.  Its honest search
lesson is:

```text
smooth:     one point + arbitrary hidden index + flat sidecar;
polynomial: finite visible degree + at least two infinity roots + Newton ledger;
cap 38:     the sharp visible concentration is (35,1).                (29)
```

Any reciprocal-inspired cap-`38` ansatz should therefore start with `(28)`,
not a forbidden one-root power, and should budget its desired phase/index in
the subleading Newton faces rather than imagine a flat repair.

## 8. The involution gate: the half-power trick cannot remain equivariant

The smooth construction uses an involution essentially: going around an odd
half-power changes `w` to `-w`, and the even seed makes the two branches agree.
A direct algebraic imitation is dangerous in JC2 because involutive symmetry
is a solved stratum.

For a reflection-intertwining Keller map, choose coordinates so that

```text
sigma(x,y)=(x,-y),       tau(P,Q)=(P,-Q),
F sigma=tau F.                                             (30)
```

Then `P` is even and `Q` odd in `y`.  On the fixed line `y=0`,

```text
Q=0, P_y=Q_x=0,
det JF=P_x Q_y=constant !=0.                              (31)
```

Both one-variable polynomial factors in `(31)` are constant, so `P(x,0)` is
affine and the restriction of `F` to the fixed line is injective.
Gwozdziewicz's **CITED** injectivity-on-one-line theorem then makes `F` an
automorphism.
Thus a JC2 counterexample cannot intertwine source and target reflections.

For the fixed-point-only central involution, Miyanishi's **CITED** theorem is
stronger: an etale endomorphism of `A^2_C` commuting with an effective finite
group of even order is an automorphism.  In particular, an odd map satisfying
`F(-x,-y)=-F(x,y)` is an automorphism.  The source is
M. Miyanishi, [*Equivariant Jacobian Conjecture in dimension
two*](https://arxiv.org/abs/2110.06709), Theorem in the introduction.

Hence the even-seed branch resolution may inspire coordinates, but any actual
Keller search using it must break the corresponding global equivariance.  The
involution can be a local bookkeeping device, not a symmetry of the final map.

## 9. A concrete JC2 squeeze: the complete quadratic-`u`, conductor-square cell

Return to the THM-3586 normalization

```text
X=t^2-1,       Y=t(t^2-1),
A=a+Yc,        B=b+Yd,
E=1,           O=0.                                     (32)
```

Preserve the minimal nodal boundary and both branch jets:

```text
a0=u^2-1,                  b0=u^3-u,
a1=2uh,                    b1=(3u^2-1)h,
c0=1/2+2uk,                d0=3u/4+(3u^2-1)k.          (33)
```

Search the first genuine conductor-square layer

```text
a=a0+Xa1+X^2p,             b=b0+Xb1+X^2q,
c=c0+Xr,                   d=d0+Xs,                    (34)
```

in the entire rational universe

```text
h,k,z,w in Q[u],       deg(h),deg(k),deg(z),deg(w)<=2. (35)
```

There is **no numerator/denominator height cap** in `(35)`.  Solving the first
`X`-equations exactly gives gauges `z,w`:

```text
F=-(3u^2+1)(h^2+k^2)-(3/2)uk+k'/4-3/16,
G=-2(3u^2+1)hk-(3/2)uh+h'/4+3/8,

p=-F+2uz,                 q=-(3/2)uF+(3u^2-1)z,
r=-G+2uw,                 s=-(3/2)uG+(3u^2-1)w.       (36)
```

The conductor period is computed at the initial quadratic stage through the
closed THM-3586 coefficient pairing, before the Keller elimination or height
search.  It is retained as an independent global necessary check, not
inserted into the PDE substitutions.  Its defect has `u`-degree at most seven;
for example, its leading coefficient is checked exactly as

```text
[u^7](Phi'-2)=(512/105)(-2h2^2w2+4h2k2z2+3k2^3-2k2^2w2).
```

Write `h=h0+h1u+h2u^2` and similarly.  Three definite leading forms of the
Keller odd equation give

```text
[X^3u^12]O = -54(h2^4+6h2^2k2^2+k2^4),
[X^3u^8]O  = -54(h1^4+6h1^2k1^2+k1^4) after h2=k2=0,
[X^3u^6]O  = -24(z2^2+w2^2) after h1=k1=h2=k2=0.      (37)
```

Over `Q` (indeed any ordered field), these force successively

```text
h2=k2=0,             h1=k1=0,             z2=w2=0.
```

Thus the quadratic universe reduces exactly to the previously audited affine
cell.  Its remaining coefficient ladder is

```text
[X^4u^4]O = -30(3h0k0+w1)^2,
[X^3u^4]O = -6(3h0^2+3k0^2+2z1)^2,
[X^4u^2]O = -(15/8)(3h0+4w0)^2,
[X^3u^2]O = -(3/2)(3k0+4z0)^2.                        (38)
```

The remaining squares force

```text
w1=-3h0k0,
z1=-(3/2)(h0^2+k0^2),
w0=-3h0/4,
z0=-3k0/4.                                             (39)
```

The entire cell collapses to

```text
p=(16h0^2+16k0^2+3)/16,
q=3((32h0^2+32k0^2+3)u+8k0)/32,
r=(16h0k0-3)/8,
s=3(32h0k0u+4h0-3u)/16.                               (40)
```

Now the Keller equation has the parameter-independent residual

```text
[X^2u^0](E-1)=-87/32.                                  (41)
```

This proves the **all-`Q` quadratic-`u` conductor-square no-go**.  It is
stronger than any small-height enumeration in this cell.  It is not an
all-degree conductor-square theorem.  The definite quartics and sums of
squares use the order on `Q`; none of their zero implications is asserted
over `C`.

The period was computed and recorded, but not imposed in the PDE elimination,
at the initial stage.  As an independent check, the companion computes it a
second way by direct pullback and termwise exact `t`-moments.  The two
expressions agree identically.  On the forced collapsed cell they give

```text
Phi=2(16h0^3+12h0^2u-16h0k0^2+12h0k0u+23h0
      +96k0^3u^2+32k0^3+84k0^2u+24k0)/105,

Phi'-2=2(4h0^2+4h0k0+64k0^3u+28k0^2-35)/35.          (42)
```

Over `Q`, `(42)` first forces `k0=0`, then `h0^2=35/4`, which has no rational
solution.  Even after adjoining `sqrt(35)`, the PDE residual `(41)` remains.
The generic formulas `(36)` are a positive control: they annihilate all
boundary and first-`X` equations identically.  The hostile perturbation
`p -> p+1` is rejected immediately because it changes the first odd equation
by `-4X(3u^2-1)`.

The exact escape routes are now sharply limited:

1. keep conductor order two but allow at least one gauge to have `u`-degree
   at least three;
2. add a genuine `I^3` or higher conductor layer;
3. leave the minimal nodal boundary/root-line jet cell;
4. over `C`, inspect the nonreal branches of the quartics and sums of squares
   in `(37)` rather than importing the ordered-field collapse.

A rational height search should not spend a single cycle on `(35)`.  It should
symbolically eliminate the next degree layer, impose `(42)` at the start, and
enumerate heights only on residual positive-dimensional branches.

## 10. Transfer matrix and next decisive tests

| source mechanism | Keller target | preserved predicate | information destroyed | needed sidecar | cheapest decisive test |
|---|---|---|---|---|---|
| flat reciprocal envelope | reciprocal/Newton ansatz | local phase multiplication | polynomiality and finite jet | pole-order ledger | restrict to a ray; finite-order envelope fails `C^infty` |
| nonzero seed spinor | `q_+` of a Hessian carrier | zero-free phase carrier | constant-Hessian constraint | factorization `(27)` | reduce `q_+q_-` modulo `det Hess H-mu` |
| index multiplication | boundary Newton face | winding | hidden flat reserve | Laurent exponent interval | apply `(22)` before coefficient search |
| support shift | trace-gauge Keller map | trace-free Hessian | positivity repair | none available | note constant `H` is invisible in `(19)` |
| even seed / half-power | involutive descent | branch agreement | symmetry freedom | explicit broken-symmetry term | fixed-line test `(31)` / Miyanishi gate |
| one-point smooth localization | degree-36 top divisor | concentration | infinite jet invisibility | second infinity root | start at `K=t^35(t+lambda u)` |
| dephasing `amplitude -> modulus^2` | `q_+ -> q_+q_-` | zero set and size | phase/index sign | phase or Laurent sidecar | compare `(21)` with product-only filter |

The most promising bounded next computation is not a blind rational box.  It
is a symbolic degree-`3` gauge elimination in `(34)--(36)`, with the period
coefficients imposed simultaneously, followed by a small-height enumeration
only if an exact branch survives.  In parallel, any trace-gauge Hessian search
should print:

```text
deg H, Newton support(q_+), circle exponent interval,
factorization of (Delta H)^2-4mu,
real spinor zeros/indices (auxiliary only),
and every imposed/broken involution.                    (43)
```

These fields make the smooth construction useful without pretending that its
flatness, convexity shift, or global one-point localization survives in the
polynomial Keller category.

## Canonical connections

- `01-canon/theorems/THM-3586-nodal-cylinder-cap38-width-period-and-second-conductor-keller-gates.md`
  -- **PROVED/CITED/FINITE-EXACT** cap, period, conductor, and `(35,1)` root
  gates; residual Keller PDE open.
- `01-canon/theorems/THM-3587-consecutive-keller-fibre-factor-toric-and-coefficient-span-gates.md`
  -- **PROVED router** to `(19)`; constant-nonzero binary Hessian problem open.
- Gwozdziewicz and Shastri are used only in their **CITED** forms already
  scoped in THM-3586.
- Miyanishi, arXiv:2110.06709, is used only for the **CITED equivariant JC2**
  theorem for effective finite groups of even order.
