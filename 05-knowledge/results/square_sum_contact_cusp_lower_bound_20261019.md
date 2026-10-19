# The prescribed contact search has logarithmic cusp complexity

**Status: PROVED.** For the prescribed plus-port curve and section, the
least successful signed index has order `log(1/f₀)` as the anchor approaches
either cusp. This supplies the matching lower rate missing from
[the logarithmic recurrence bound](square_sum_contact_recurrence_bound_20261019.md).
It covers every integer translate, including odd indices and both accepted
positive-spread charts. No novelty claim is made.

## 1. Statement and conventions

Let `t>1` be a rational admissible anchor and write

\[
f(u)=\frac{4u(u^2-1)}{(u^2+1)^2},\qquad
\delta=f(t)\in(0,1/2),\qquad K=1+\delta=k^2.
\]

Use exactly the prescribed plus-port objects

\[
E_k:Y^2=(X-K)(X^2-1),\qquad P=(0,k),
\]
\[
Q=\left(\frac{t^2+4t+1}{t^2+1},
\frac{8t(t+1)}{(t^2+1)^2}\right).
\]

Let `m(t)` be the least `|q|` for which `Q+qP` yields a different admissible
spread. Both charts `u>1` and `−1<u<0` are accepted; in the latter use
`t′=−1/u`. The preceding universal-contact and recurrence results prove
that `m(t)` is finite. Put

\[
\Omega=\int_K^\infty\frac{dX}{\sqrt{(X-K)(X^2-1)}},
\quad
\alpha(K)=\frac12\int_{K+1/(4K)}^\infty
\frac{dX}{\sqrt{(X-K)(X^2-1)}},
\tag{1}
\]
\[
C_0=\frac\pi4+\frac{\pi\sqrt2}{12}.
\tag{2}
\]

The period uses `ω=dX/(2Y)`, as in the recurrence note.

**Theorem.** If `0<δ≤1/1000`, every successful signed index satisfies

\[
\boxed{\quad |q|\ge \frac{\Omega-2C_0}{\alpha(K)}-1.\quad}
\tag{3}
\]

For an even successful index, the final `−1` can be omitted. Consequently

\[
\boxed{\quad
m(t)\ge \frac{1}{2\log(1+\sqrt2)}\log\frac1{f(t)}-O(1)
\quad(f(t)\longrightarrow0^+).\quad}
\tag{4}
\]

The constant in `O(1)` is absolute and works at both `t→1+` and `t→∞`.
Together with the preceding signed-index upper bound this proves

\[
m(t)=\Theta\!\left(\log\frac1{f(t)}\right)
\]

at both cusps. This does not assert an asymptotic equivalent for `m(t)`.

## 2. Two bounded admissibility clusters on each real component

The smooth completion of the real quartic is

\[
z^2=(u^2+1)^2(K-f(u)).
\]

Write `u=tan θ`, with `θ` modulo `π`. Its two real components are
parametrized by

\[
\frac{z}{u^2+1}=\sigma\sqrt{K+\sin4\theta},\qquad
\sigma\in\{+1,-1\}.
\]

The known differential identity `ω=−du/z` gives on either component

\[
|\omega|=\frac{d\theta}{\sqrt{K+\sin4\theta}}
         =\frac{d\theta}{\sqrt{K-f(u)}}.
\tag{5}
\]

The complete real curve, rather than just an affine inverse chart, is
being used here. Each component is a circle of invariant length `Ω`;
the nonidentity component is a torsor for the identity circle.

Enlarge the success set to all points with `0<f(u)<1/2`, even allowing
the anchor spread itself. This contains all accepted chart outputs.
The closed set

\[
H=\{\theta:f(\tan\theta)\le1/2\}
\]

consists of exactly two arcs:

\[
[-\pi/24,7\pi/24],\qquad [11\pi/24,19\pi/24]
\quad\pmod\pi.
\tag{6}
\]

Each arc has a negative-spread middle interval of angular width `π/4`
and two positive-spread end intervals, each of width `π/24`. On the
middle interval `K−f≥1`; on the ends `K−f≥1/2`. Thus each arc has
invariant length at most `C₀` from (2).

The shift `θ↦θ+π/2` preserves (5) and interchanges the arcs. It is
therefore exactly a half-period shift: corresponding points are separated
by `Ω/2`. In particular, if a base point `B` belongs to either arc, every
point of the enlarged success set has Abelian displacement from `B`
within `C₀` of `0` or of `Ω/2`, modulo `Ω`.

It follows that a displacement `±d` cannot be successful whenever

\[
C_0<d<\Omega/2-C_0.
\tag{7}
\]

This includes both accepted charts and both signs of `z`. Possible
exceptional affine inverse-chart points cause no gap: we have enlarged
the acceptance set on the complete curve.

## 3. The translating step is longer than a cluster

Doubling `P` gives exactly

\[
A=2P=\left(K+\frac1{4K},\frac{1-4K^2}{8Kk}\right).
\tag{8}
\]

Its `Y` coordinate is negative. Starting at infinity along that half of
the identity component, its Abelian coordinate is precisely `α(K)` in
(1), with `0<α(K)<Ω/2`. On either real component, translation by `A`
therefore adds `α(K)` modulo `Ω`.

At the cusp,

\[
\begin{aligned}
\alpha_0:=\alpha(1)
 &=\frac12\int_{5/4}^\infty\frac{dX}{(X-1)\sqrt{X+1}}\\
 &=\int_{3/2}^\infty\frac{dr}{r^2-2}
 =\sqrt2\log(1+\sqrt2).
\end{aligned}
\tag{9}
\]

In particular `α₀>C₀`. Here is a completely explicit range for the needed
strict inequality. Set `v=X−K` in (1), restrict the integral to `v≥1/4`,
and bound its two changing factors there. This gives

\[
\alpha(K)\ge
\frac{\alpha_0}{\sqrt{(1+4\delta)(1+4\delta/9)}}.
\tag{10}
\]

The elementary inequality `(1+x²)^(−1/2)>1−x²/2` for `0<x≤1` yields
`asinh(1)>5/6`, so `α₀>5√2/6>7/6`. If `δ≤1/1000`, then

\[
(1+4\delta)(1+4\delta/9)
\le\frac{565001}{562500}<\left(\frac{1003}{1000}\right)^2.
\]

Using `π<22/7` and `√2<10/7`,

\[
\boxed{\quad
\alpha(K)>\frac{3500}{3009}>
\frac{341}{294}>C_0.
\quad}
\tag{11}
\]

Every inequality is strict in the direction required by (7).

## 4. Even and odd translates

For even `q=2n`, take `B=Q`. Its spread is `δ<1/2`, so `B∈H`.
The case `n=0` is the excluded anchor. For every other `n`, the real
representative `d=|n|α(K)` is greater than `C₀`, by (11). If it were less
than `Ω/2−C₀`, (7) would forbid success. Therefore every successful even
index obeys

\[
|n|\alpha(K)\ge\Omega/2-C_0,
\quad\text{hence}\quad
|q|\ge\frac{\Omega-2C_0}{\alpha(K)}.
\tag{12}
\]

For odd indices, choose a base point on the other component that is known
to lie in `H` and is itself unsuccessful:

- If `1<t<4/3`, use `B=Q+P` and `s=+1`. The exact inverse parameter is
  \[
  u(B)+1=\frac{(t-1)(3t^2-1)}{t(t^2-4t+1)}<0.
  \]
  Hence `u(B)<−1` and `f(u(B))<0`.
- Otherwise, admissibility implies `t>√3`. Indeed, `f` is increasing
  on `[4/3,√3]` and `f(4/3)=336/625>1/2`. Use `B=Q−P` and `s=−1`.
  Its exact parameter is
  \[
  u(B)=\frac{t^2+4t+1}{t^3+2t^2+t-2}\in(0,1),
  \]
  since the denominator minus the numerator is `(t+1)(t²−3)>0`.
  Again its spread is negative.

Every odd index has the form `q=s+2n`, so `Q+qP=B+nA`.
The case `n=0` fails by the negative-spread calculation. For `n≠0`,
apply exactly (7) and (11) on this component to obtain

\[
|n|\alpha(K)\ge\Omega/2-C_0.
\]

Finally `|s+2n|≥2|n|−1`, proving (3). This argument handles all odd
indices directly. No fixed-index continuity assertion is being applied
to an increasing list of indices.

## 5. Cusp rate and precise scope

The shifted integral used in (10) has its lower endpoint bounded away
from zero as `K→1`. Differentiation under the integral, or an elementary
uniform comparison of its factors, therefore gives

\[
\alpha(K)=\alpha_0+O(\delta).
\]

The period expansion established with the same normalization in the
recurrence note is

\[
\Omega=\frac1{\sqrt2}\log\frac{32}{\delta}
       +O\!\left(\delta\log\frac1\delta\right).
\]

Substitute these into (3). The coefficient of `log(1/δ)` is

\[
\frac1{\sqrt2\alpha_0}
=\frac1{2\log(1+\sqrt2)},
\]

which proves (4). Alternatively, an elementary estimate
`Ω=(1/√2)log(1/δ)+O(1)` already suffices for the same lower statement.

The result is specific to this prescribed plus-port orbit. It is not a
lower bound for other starting sections, other translating points, the
minimum height of any neighbor, or bit complexity. The accompanying
upper bound makes the order of the signed translate index sharp; it does
not make this leading lower coefficient an asymptotic formula for the
first success.

The companion checker verifies the exact algebra and rational constant
inequalities, plus separately labelled high-precision period/step and
first-success diagnostics. No numerical quadrature is used as a proof
certificate.


## 6. Consequence along the outer ternary-tree spines

Starting from the root parameter `t=2`, the two outer Berggren spines have

\[
t_j^{(L)}=\frac{j+2}{j+1},\qquad t_j^{(R)}=2j+2.
\]

Their primitive hypotenuses are respectively

\[
c_j^{(L)}=2j^2+6j+5,\qquad c_j^{(R)}=4j^2+8j+5.
\]

In both cases the reduced numerator and denominator have opposite parity,
so no factor-two normalization occurs. Direct substitution gives
`f(t_j)=2/j+O(1/j²)`; both spines are admissible for `j≥3`. Hence

\[
\boxed{\quad m(t_j)=\Theta(\log j)=\Theta(\log c_j).\quad}
\]

More explicitly, the lower bound is
`m(t_j)≥[2 log(1+√2)]⁻¹ log j−O(1)` on either spine.
Thus this particular contact search requires logarithmically many elliptic
translations against linear Berggren depth. This is an index comparison,
not a bound on graph distance, rational-coordinate height, or bit cost.


## 7. Reproduce the controls

With SymPy and mpmath installed, run:

```bash
python3 square_sum_contact_cusp_lower_bound_20261019.py controls.json
python3 -O square_sum_contact_cusp_lower_bound_20261019.py controls.optimized.json
```

The checker retains 428 explicit checks, including 21 exact algebra/constant/spine controls. Normal and optimized Python runs produce identical JSON. The exact orbit probes use both signs and both accepted charts. Their
first-success absolute indices are `13, 8, 13, 19` for anchors `10000`,
`10⁶`, `10⁹`, and `1000001/1000000`, respectively. These are small targeted
regressions, not a further exhaustive census.
