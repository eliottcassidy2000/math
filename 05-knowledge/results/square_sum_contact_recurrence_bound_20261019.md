# Elliptic contact search has a logarithmic anchor-dependent index bound

**Status: PROVED**, with exact algebraic controls and separately labelled
high-precision numerical checks. An independent adversarial audit passed
the proof and companion checker, including 49,622 additional controls. The predecessor proves that no constant
translate budget works for every anchor. This note supplies the complementary
upper bound: the prescribed signed search has a successful index

\[
 |q|=O\!\left(1+\log\frac1{f_0}\right)
\]

for every admissible primitive anchor of spread `f_0`. The constant is
absolute; either mixed-sign port is allowed. The successful index may have
either sign, and its spread direction is not prescribed in this global
bound. Separate explicit bounds prescribe increasing or decreasing spread,
avoid a finite forbidden set, or impose a local spread window.

This is an index bound for the existing orbit `Q+qP`, not a bound on the
smallest possible neighbor's hypotenuse, not a polynomial-time assertion,
and not a matching logarithmic lower bound. No novelty claim is made.

## 1. Controller, real components, and period convention

Use the controller of
[universal signed Pythagorean contacts](square_sum_ppt_bridges_20261019.md)
and the orbit considered in
[the no-uniform-budget theorem](square_sum_contact_step_bound_20261019.md).
For an anchor parameter `t_0>1`, let

\[
 f(t)=\frac{4t(t^2-1)}{(t^2+1)^2},\qquad f_0=f(t_0)\in(0,1/2).
\]

Fix `epsilon=+1` or `-1`, set `K=1+epsilon*f_0=k^2`, and use

\[
 E_k:Y^2=(X-K)(X^2-1),\qquad P=(0,k).
 \tag{1}
\]

The anchor `Q` is the image of
`(u_0,z_0)=(epsilon*t_0,t_0^2+1)` on

\[
 z^2=K(u^2+1)^2-4u(u^2-1).
 \tag{2}
\]

The predecessor proves that `P` is non-torsion. Write

\[
 e_1=\max(K,1),\quad e_2=\min(K,1),\quad e_3=-1,
 \qquad D=e_1+1=\begin{cases}2+f_0,&\epsilon=+1,\\2,&\epsilon=-1.\end{cases}
 \tag{3}
\]

There are two real components. The identity component is the unbounded
branch `X>=e_1` together with the point at infinity. The bounded oval has
`-1<=X<=e_2`. Since `0` is strictly inside the latter interval, `P` is on
the nonidentity component, and `A=2P` is on the identity component.

Normalize with the invariant differential

\[
 \omega=\frac{dX}{2Y}.
\]

Its positive real period is

\[
 \begin{split}
 \Omega
 &=\int_{e_1}^{\infty}
       \frac{dX}{\sqrt{(X-e_1)(X-e_2)(X+1)}}\\
 &=\frac{2}{\sqrt D}\,\mathsf K\!\left(\sqrt{1-\frac{f_0}{D}}\right),
 \end{split}                                                    \tag{4}
\]

where `mathsf K(kappa)=integral_0^(pi/2)
(1-kappa^2 sin^2 theta)^(-1/2) dtheta` uses the **modulus**, not the
parameter. The factor `2` in (4) is essential: the integral with
`dX/(2Y)` along one half of the identity component is only half a period.
The identity component is `R/(Omega Z)` in this coordinate, and every
other component is a torsor for that circle. Translation by `A` is an
irrational rotation: a rational rotation number would make a nonzero
multiple of `A`, and hence of `P`, the identity.

The real Lie-group description is standard; see
[Elkies's real-group discussion](https://people.math.harvard.edu/~elkies/M223.24/index.html).
The remaining recurrence estimate is proved below and needs no quantitative
equidistribution theorem or Diophantine estimate for the rotation number.

## 2. A signed circular-pigeonhole lemma

**Lemma 1.** Let `A` be a non-torsion point on a real elliptic identity
component of period `Omega`. Let an oriented, open arc beginning at `Q`
have invariant-differential length `ell`, with `0<ell<=Omega/2`.
There exist an integer `n` and a sign `sigma` such that

\[
 1\le n\le\left\lfloor\frac{\Omega}{\ell}\right\rfloor,
 \qquad Q+\sigma nA
\]

lies strictly inside that arc and is distinct from `Q`.

**Proof.** Put `N=floor(Omega/ell)`. The `N+1` distinct points
`0,A,...,NA` cut the circle into `N+1` gaps, one of length at most
`Omega/(N+1)<ell`. The difference of its endpoints, in the circular order,
is `nA` or `-nA` for some `1<=n<=N`, including when the gap crosses zero.
Choose its sign to agree with the specified orientation. Its distance is
positive because `A` is non-torsion. Translate by `Q`. QED.

Taking `A=2P` gives the bound

\[
 \boxed{\quad 0<|q|\le2\left\lfloor\frac{\Omega}{\ell}\right\rfloor,
       \qquad q\text{ even}.\quad}                         \tag{5}
\]

Both signs of `q` matter. A tiny positive irrational rotation can miss a
negative-side arc for arbitrarily many positive iterates. Lemma 1 reverses
the displacement when necessary; it does not assert a comparable bound
for `q>0` only, or for hitting an arbitrary arc not attached to `Q`.

## 3. Certified arcs, monotonicity, and distinct primitive neighbors

The explicit quartic-to-cubic map in the predecessor satisfies

\[
 \omega=-\frac{du}{z}.                                  \tag{6}
\]

This follows by differentiating its formula for `X` and using (2).
On the positive-`z` anchor branch, put `u=epsilon*t`. Then

\[
 z=(t^2+1)\sqrt{1+\epsilon(f_0-f(t))},\qquad
 |\omega|=\frac{dt}{(t^2+1)\sqrt{1+\epsilon(f_0-f(t))}}.
 \tag{7}
\]

The radicand is between `1/2` and `3/2` throughout the admissible region,
so this branch is smooth and `t` is a genuine local coordinate.
Let `theta=arctan(t)`. Since `f(t)=-sin(4theta)`, the two admissible
intervals are

\[
 \frac\pi4<\theta<\frac{7\pi}{24},\qquad
 \frac{11\pi}{24}<\theta<\frac\pi2.                       \tag{8}
\]

On each interval the spread is strictly monotone, and its endpoints have
spreads `0` and `1/2`. With `s=f(t)`, equation (7) therefore becomes

\[
 |\omega|=\frac{ds}{4\sqrt{1-s^2}\sqrt{1+\epsilon(f_0-s)}}.
 \tag{9}
\]

In particular, the exact lengths to the decreasing and increasing
spread boundaries are

\[
 L_- =\frac14\int_0^{f_0}
 \frac{ds}{\sqrt{1-s^2}\sqrt{1+\epsilon(f_0-s)}},
 \qquad
 L_+ =\frac14\int_{f_0}^{1/2}
 \frac{ds}{\sqrt{1-s^2}\sqrt{1+\epsilon(f_0-s)}}.           \tag{10}
\]

They depend on the anchor only through `f_0` and `epsilon`, not on which
swapped-leg parameter represents its petal. Both arcs are shorter than
`Omega/2`: their total length is at most `sqrt(2)*pi/24`, whereas
`Omega>=pi/sqrt(D)>=pi/sqrt(5/2)`.

Every interior point of either arc has `t>1`, `0<f(t)<1/2`, and spread
strictly different from `f_0`. The birational inverse is in its ordinary
chart there: the predecessor's `X=0` exception has
`|u|=1/(2K)<1`, and `X=K` cannot occur at finite nonzero `u` for `k!=1`.
Thus a rational translate in the arc yields a rational parameter and

\[
 \left(\frac{z}{t^2+1}\right)^2=1+\epsilon f_0-\epsilon f(t).
 \tag{11}
\]

Primitive normalization is essential: for reduced `t=n/m`, divide
`(n^2-m^2,2nm,n^2+m^2)` by `2` if `n,m` are both odd, and by `1` otherwise,
then reorder the legs. Spread is unchanged. Its strict change ensures a
new primitive petal because spread uniquely determines the unordered
normalized positive leg pair.

**Theorem 2. Directional and local bounds.** There is a larger-spread
opposite-port neighbor at an even signed index

\[
 0<|q|\le2\left\lfloor\frac\Omega{L_+}\right\rfloor,
\]

and a smaller-spread one with `L_+` replaced by `L_-`.
More generally, replace either endpoint in (10) by any permitted target
spread in that direction, and use the resulting shorter arc length.
Every interior point then stays within that spread window.

For a finite forbidden set, take the endpoint to be the nearest forbidden
spread in that direction or the admissibility boundary, whichever comes
first; ignore forbidden spreads equal to `f_0`. The open arc avoids all
forbidden petals. A strict requested window is enforced by making its
edge another candidate endpoint. The interval is nonempty because the
forbidden set is finite. Lemma 1 then proves the stated finite bound,
rather than only eventual termination.

## 4. A global logarithmic index bound

Each angular branch in (8) has width `pi/24`. At least one of the two
angular distances from the anchor to its endpoints is at least `pi/48`.
Equation (7), and the radicand bound `<=3/2`, give

\[
 \max(L_-,L_+)\ge\frac\pi{24\sqrt6}.                    \tag{12}
\]

**Theorem 3. Global signed-index bound.** For either anchor port there is
a distinct admissible primitive opposite-port neighbor in the orbit
`Q+qP`, with

\[
 \boxed{\quad
 0<|q|\le
 2\left\lfloor\frac{24\sqrt6}{\pi}\,\Omega\right\rfloor.
 \quad}                                                \tag{13}
\]

Here `q` may have either sign. Its spread direction is chosen by the
longer admissible arc and is not specified in advance.

There is also a wholly elementary upper bound. Write `r=sqrt(f_0/D)`.
The tangent substitution gives

\[
 \mathsf K(\sqrt{1-r^2})
 =\int_0^\infty\frac{dv}{\sqrt{(1+v^2)(1+r^2v^2)}}.
\]

Split this integral at `1` and `1/r`. The first piece is at most
`asinh(1)`, and the middle piece at most `log(1/r)`. In the last piece,
put `w=rv`; its integrand is at most `1/(w sqrt(1+w^2))`, whose integral
from `1` to infinity is `asinh(1)`. Hence

\[
 \Omega\le\frac{4\operatorname{arsinh}(1)+\log(D/f_0)}{\sqrt D}
 \le\frac{4\operatorname{arsinh}(1)+\log((5/2)/f_0)}{\sqrt2}.
 \tag{14}
\]

Combining (13) and (14), one successful even signed index satisfies

\[
 \boxed{\quad
 |q|\le
 2\left\lfloor\frac{24\sqrt3}{\pi}
 \left[4\operatorname{arsinh}(1)+\log\frac5{2f_0}\right]
 \right\rfloor.
 \quad}                                                 \tag{15}
\]

This proves the claimed logarithmic anchor-dependent bound globally,
including anchors arbitrarily close to either admissibility wall.
For a primitive triple `(a,b,c)`, `1/f_0=c^2/(2ab)<=c^2/2`, so this also
gives `|q|=O(log c)`. It is only an input-hypotenuse bound on the search
index; rational-coordinate sizes and the output hypotenuse are separate.

### An exact rational stopping budget

No unverified numerical elliptic integration is needed to implement a
certified finite search. Let `J` be the least nonnegative integer with
`2^J f_0>=D`, and define

\[
 W=\frac57(4+J).                                        \tag{16}
\]

The coarser estimate from the same three-piece split is
`Omega< (4+log(D/f_0))/sqrt(D) < W`, using
`log 2<1` and `sqrt(D)>=sqrt2>7/5`.
The arc integrand in (9) is strictly greater than `1/5`, since
`4 sqrt(1-s^2) sqrt(1+epsilon(f_0-s))<4 sqrt(3/2)<5`.
An allowed spread interval of width at least a rational `Delta>0` therefore
has length `ell>Delta/5`. For a real requested window, choose any positive
rational lower bound on its width. Thus the entirely rational budget

\[
 0<|q|\le2\left\lfloor\frac{5W}{\Delta}\right\rfloor     \tag{17}
\]

is valid for that direction/window. One of the two full spread distances
is at least `1/4`, giving the particularly simple global certificate

\[
 \boxed{\quad
 0<|q|\le2\left\lfloor\frac{100(4+J)}7\right\rfloor.
 \quad}                                                 \tag{18}
\]

Every quantity in (16)--(18) is obtained by exact rational/integer
arithmetic. Enumerating the even signed indices within this budget and
using exact contact tests must succeed. For sharper certified bounds,
use an upper enclosure for `Omega` and a positive lower enclosure for the
chosen arc length in (5). There is no need to decide whether an exact
real ratio is an integer.

## 5. What the cusp calculation does and does not prove

The standard complete-integral expansion gives

\[
 \Omega=\frac1{\sqrt2}\log\frac{32}{f_0}
       +O\!\left(f_0\log\frac1{f_0}\right)
 \qquad(f_0\longrightarrow0^+),                         \tag{19}
\]

for either port. See
[NIST DLMF, equation 19.12.1](https://dlmf.nist.gov/19.12#E1).
The period convention in (4) fixes the constant `32`.

The increasing-spread arc stays macroscopic:

\[
 L_+\longrightarrow C_\epsilon
 =\frac14\int_0^{1/2}\frac{ds}{\sqrt{1-s^2}\sqrt{1-\epsilon s}}>0,
\]

with explicit constants

\[
 C_+=\frac{\operatorname{artanh}(\sqrt3/2)
                   -\operatorname{artanh}(1/\sqrt2)}{2\sqrt2}
       =0.1540023096562408\ldots,
\]
\[
 C_-=\frac{\operatorname{artanh}(1/\sqrt2)
                   -\operatorname{artanh}(1/2)}{2\sqrt2}
       =0.1174035702670919\ldots.                         \tag{20}
\]

In fact `L_+=C_epsilon+O(f_0)`. Thus the directional upper budget for an
increasing-spread neighbor is

\[
 \frac{\sqrt2}{C_\epsilon}\log\frac{32}{f_0}+O(1),        \tag{21}
\]

whose leading coefficients are approximately `9.18307` and `12.04575`
for the plus and minus ports, respectively. This is a bound, not an
asymptotic formula for the first successful index.

By contrast, `L_-=f_0/4+O(f_0^2)`. The direct decreasing-spread recurrence
budget is consequently

\[
 \frac{4\sqrt2}{f_0}\log\frac{32}{f_0}
     +O\!\left(\log\frac1{f_0}\right).                   \tag{22}
\]

Its inverse-spread factor is a feature of this local-arc guarantee,
not a lower bound on actual decreasing-spread search times. Shrinking a
symmetric neighborhood and then using its shorter half would needlessly
introduce that factor into the undirected/global problem.

For its plus-port orbit, the predecessor establishes divergence of the
minimal successful signed index as `t_0->1+` or `t_0->infinity`, but gives
no lower rate. Its fixed-`q`
continuity argument cannot be applied uniformly to `q` growing with the
anchor. A matching `Omega(log(1/f_0))` lower bound would require a uniform
estimate for the moving orbit (including odd translates and all accepted
charts), or an equivalent quantitative control of its departure from the
singular Pell/Berggren orbit. No such estimate is asserted here.

## 6. Reproducible controls and scope

The companion script has three independent kinds of checks:

1. Exact symbolic identities: differential pullback, angular spread,
   period substitution, and the formula for `2P`.
2. Exact rational orbit witnesses for selected anchors and both ports,
   with cubic/quartic/contact equations, changed spread, primitive
   normalization, and the exact rational stopping budget all checked.
3. High-precision numerical controls of period integrals, arc integrals,
   the elementary upper bounds, and cusp limits. These floating-point
   values are diagnostic and are not used as proof certificates.

The two executions are byte-identical and retain **999 exact rational checks,
8 exact symbolic checks, and 198 high-precision numerical checks**. The
18 selected exact witnesses cover both ports, both positive-parameter
branches, and the odd/odd normalization divisor. They are witnesses in the
positive-`z`, even-index anchor suborbit, not claims of globally first
success under the predecessor's broader acceptance rule.

For orientation, these are diagnostic real-integral values (the rational
budgets in the last column are exact):

| Spread | Port | Period | Larger arc | Arc-based budget | Rational budget |
|---|---:|---:|---:|---:|---:|
| `1/1000` | `+` | 7.334426 | 0.153643 | 94 | 428 |
| `1/1000` | `−` | 7.335906 | 0.117201 | 124 | 428 |
| `2/5` | `+` | 3.021569 | 0.094361 | 64 | 200 |
| `2/5` | `−` | 3.192170 | 0.115741 | 54 | 200 |

These real-curve spread probes need not themselves be rational PPT
anchors; the separate exact-orbit controls use actual rational PPTs.
The displayed arc-based integer budgets are numerical diagnostics; only
the rational column is advertised as certified by this computation.

Run, with SymPy and mpmath installed:

```bash
python3 04-computation/experiments/square_sum_contact_recurrence_bound_20261019.py
python3 -O 04-computation/experiments/square_sum_contact_recurrence_bound_20261019.py
```

The exact rational budget (18), rather than the numerical values, proves
termination of the prescribed finite search. The construction still uses
only the even suborbit and the raw `t>1` anchor branch; allowing all integer
indices and the additional positive-spread chart can only improve its
first success. The bound therefore also applies to the broader acceptance
rule of the predecessor.
