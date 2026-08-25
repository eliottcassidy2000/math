---
id: THM-4028
title: "Sun 2-4-6-8 average-order criticality"
status: >
  PROVED + INDEPENDENTLY HOSTILE-AUDITED. For any fixed positive degrees
  d_1,...,d_r and fixed lower index bounds, the ordered tuple count below X is
  V_d X^sigma+O(X^(sigma-1/d_max)), where sigma=sum_i 1/d_i and V_d is the
  explicit anisotropic Dirichlet volume. For degrees 2,4,6,8, sigma=25/24,
  V=24.31102486226095..., and the summatory mean representation count is
  asymptotic to V X^(1/24). On intervals H=o(X), H>>X^(7/8), the local mean is
  asymptotic to J X^(1/24), J=(25/24)V=25.32398423152183.... No pointwise
  asymptotic or positivity follows; THM-4026 has count zero at its target.
source: root + independent anisotropic lattice-count audit, 2026-08-24
audit: >
  PASS. An independent proof audit checked the squeeze between two shifted
  pure-power sublevel sets, the unit-cube boundary error, all gamma arguments,
  the distinction V versus the shell derivative J, the mesoscopic threshold
  H>>X^(7/8), and the inability to difference the error pointwise.
depends_on: []
related:
  - THM-4026-sun-two-four-six-eight-binomial-counterexample
  - THM-4027-sun-two-four-six-eight-universal-modular-solubility
  - HYP-1953-additive-basis-normal-form-spectrum
  - MISTAKE-209
  - MISTAKE-219
---

# THM-4028 -- the 2-4-6-8 system is barely supercritical on average

**PROVED + INDEPENDENTLY HOSTILE-AUDITED.** This theorem separates tuple
volume, average multiplicity, and pointwise support.

## 1. General anisotropic count

Fix positive integers `d_1,...,d_r` and fixed nonnegative lower bounds `L_i`.
Let

\[
A_{\mathbf d}(X)=\#\left\{(t_1,\ldots,t_r):
t_i\in\mathbb Z,\ t_i\ge L_i,\
\sum_i{t_i\choose d_i}\le X\right\}.                  \tag{1}
\]

Put

\[
\sigma=\sum_i{1\over d_i},\qquad d_{\max}=\max_i d_i, \tag{2}
\]

and

\[
V_{\mathbf d}=
{\displaystyle\prod_i(d_i!)^{1/d_i}\Gamma(1+1/d_i)
 \over\displaystyle\Gamma(1+\sigma)}.                 \tag{3}
\]

Then

\[
A_{\mathbf d}(X)=V_{\mathbf d}X^\sigma
+O_{\mathbf d,\mathbf L}\!\left(X^{\sigma-1/d_{\max}}\right). \tag{4}
\]

Repeated zero values caused by indices below `d_i` are included literally in
(1). They live on lower-dimensional coordinate faces and affect only the
error term.

## 2. Proof

For `t>=0`,

\[
{(t-d+1)_+^d\over d!}\le {t\choose d}\le {t^d\over d!}. \tag{5}
\]

Thus (1) is squeezed between lattice counts for the pure-power body

\[
\mathcal R_X=\left\{u_i\ge0:\sum_i{u_i^{d_i}\over d_i!}\le X\right\}, \tag{6}
\]

with only fixed coordinate shifts and finitely thick coordinate faces.

Let `B(X)=#(R_X cap Z^r)`. Unit cubes based at counted lattice points cover
`R_X`, while every such cube is contained in

\[
\mathcal R_{X+C X^{1-1/d_{\max}}}}                  \tag{7}
\]

for a constant depending only on the degrees. Scaling
`u_i=X^(1/d_i)v_i` therefore gives

\[
B(X)=\operatorname{vol}(\mathcal R_1)X^\sigma
+O\!\left(X^{\sigma-1/d_{\max}}\right).                \tag{8}
\]

The fixed shifts in (5), the lower bounds, and the zero-coordinate fibres
contribute at most
`O(sum_i X^(sigma-1/d_i))`, which has the same upper order as the error in
(8). Finally set `s_i=u_i^(d_i)/d_i!` in the volume integral. The Dirichlet
simplex integral gives

\[
\operatorname{vol}(\mathcal R_1)
=\prod_i{(d_i!)^{1/d_i}\over d_i}
 {\prod_i\Gamma(1/d_i)\over\Gamma(1+\sigma)}
=V_{\mathbf d}.                                        \tag{9}
\]

This proves (4).

## 3. Specialization to Sun's degrees

Let `a(n)` be the number of ordered quadruples `w,x,y,z>=2` satisfying

\[
n={w\choose2}+{x\choose4}+{y\choose6}+{z\choose8}.       \tag{10}
\]

Here

\[
\sigma={1\over2}+{1\over4}+{1\over6}+{1\over8}
={25\over24},\qquad d_{\max}=8.                        \tag{11}
\]

Equations (3)--(4) become

\[
\sum_{n\le X}a(n)=VX^{25/24}+O(X^{11/12}),             \tag{12}
\]

where

\[
V={
(2!)^{1/2}(4!)^{1/4}(6!)^{1/6}(8!)^{1/8}
\Gamma(3/2)\Gamma(5/4)\Gamma(7/6)\Gamma(9/8)
\over \Gamma(49/24)}
=24.31102486226095\ldots.                               \tag{13}
\]

Consequently the cumulative mean is

\[
{1\over X}\sum_{n\le X}a(n)\sim V X^{1/24}.           \tag{14}
\]

The formal shell constant is different:

\[
J={25\over24}V
={\prod_{d\in\{2,4,6,8\}}(d!)^{1/d}\Gamma(1+1/d)
  \over\Gamma(25/24)}
=25.32398423152183\ldots.                               \tag{15}
\]

At the THM-4026 target, `V N^(1/24)=102.052225...`, while
`J N^(1/24)=106.304401...`. The former is the cumulative mean scale and the
latter is only a shell heuristic at a point; neither is `a(N)`, which is
exactly zero.

## 4. A valid mesoscopic consequence

If `H=o(X)` and

\[
{H\over X^{7/8}}\longrightarrow\infty,                 \tag{16}
\]

then subtracting (12) at `X+H` and `X` is legitimate at the interval scale:

\[
{1\over H}\sum_{X<n\le X+H}a(n)
\sim JX^{1/24}.                                         \tag{17}
\]

Indeed, the main increment is `J H X^(1/24)(1+o(1))`, while the inherited
error is `O(X^(11/12))`; their ratio is `O(X^(7/8)/H)`. At `H=1` the error is
vastly larger than the proposed main term, so pointwise differencing is
invalid.

## 5. The reciprocal-degree threshold

For a general degree tuple, (4) gives a sharp first gate:

- if `sigma<1`, the number of represented integers below `X` is
  `O(X^sigma)=o(X)`, so the representation set has density zero and cannot be
  cofinite;
- if `sigma=1`, its upper density is at most `min(1,V_d)`; in particular,
  `V_d<1` rules out cofiniteness; and
- if `sigma>1`, the mean multiplicity grows, but that is not sufficient for
  universality. Degrees `(2,4,6,8)` have surplus only `1/24`, are locally
  universal by THM-4027, and still miss the integer in THM-4026.

Thus reciprocal-degree volume is a necessary capacity coordinate, not an
exact coverage criterion. This is the proved correction to the tempting
“enough average representations force every target” heuristic.

## 6. Non-consequences and analytic frontier

Equation (12) does not imply a pointwise Hardy--Littlewood formula, a positive
singular series lower bound, density-one support, finitely many exceptions, or
even that zeros are rare. A small set of integers could carry much of the
ordered multiplicity. Local factors quantify residue bias but do not supply
independence or concentration.

The next analytic target is therefore a lower-tail theorem: control the
variance or an exceptional-set moment for `a(n)` after conditioning on the
true prime-power local factors. The THM-4026 zero is the required hostile
example for any proposed Poisson, independence, or mean-to-minimum transfer.
