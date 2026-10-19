# Sharp cubic output length in the prescribed plus-port cusp orbits

**Status: PROVED**, using the previously proved cusp-index lower bound and the
classical specialization theorem for canonical heights. Along either outer
Berggren spine, every successful point in the prescribed plus-port orbit has
primitive output hypotenuse with logarithmic size at least a constant times
`(log j)^3`. The least possible output size in that orbit, and the output size
at any least-absolute-index success, have exactly this order. No novelty claim
or explicit optimal constant is made.

## 1. Precise theorem

Put

\[
f(t)=\frac{4t(t^2-1)}{(t^2+1)^2},\quad
k(t)=\frac{t^2+2t-1}{t^2+1},\quad K(t)=k(t)^2=1+f(t),
\]
\[
E_t:Y^2=(X-K(t))(X^2-1),\quad P_t=(0,k(t)),
\]
\[
Q_t=\left(\frac{t^2+4t+1}{t^2+1},
            \frac{8t(t+1)}{(t^2+1)^2}\right).
\]

For either spine

\[
t_j^L=\frac{j+2}{j+1},\qquad t_j^R=2j+2,
\]

consider every integer `q` for which `R=Q_{t_j}+qP_{t_j}` yields a different
admissible spread through either accepted quartic chart `u>1` or `-1<u<0`.
Write `C'_j(q)` for the resulting primitive positive hypotenuse. On the first
chart set `v=u`; on the second set `v=-1/u`.

**Theorem.** There are absolute constants `c,c_1>0` and `J_0` such that, for both
spines, every `j>=J_0` and every successful integer `q` satisfy

\[
\boxed{\ \log C'_j(q)\ge c q^2\log j
                    \ge c_1(\log j)^3.\ }                 \tag{1}
\]

If `S_j` is the nonempty set of successful signed indices, then

\[
\boxed{\ \min_{q\in S_j}\log C'_j(q)=\Theta((\log j)^3).\ }\tag{2}
\]

If `q_j` is any success with least absolute index, then

\[
\boxed{\ \log C'_j(q_j)=\Theta((\log j)^3).\ }             \tag{3}
\]

The constants can be chosen simultaneously for the two spines. Binary output
length has the same order. Since the input hypotenuses are
`2j^2+6j+5` and `4j^2+8j+5`, their bit lengths have order `log j`; hence the
cubic output-length upper bound is sharp for this prescribed construction.
This is an output-size lower bound, not a matching lower bound for the
previous conservative running-time exponent.

Every successful index has the lower bound (1). We do **not** assert that all
successful indices have output length `O((log j)^3)`: the indices are unbounded.
Nor is this a lower bound for arbitrary rational neighbors, other starting
sections, other translating points, or the minus-port orbit.

## 2. The fixed generic section has positive canonical height

Use the normalization

\[
\widehat h_t(R)=\frac12\lim_{n\to\infty}
                    4^{-n}h(X([2^n]R)),
\qquad h(a/b)=\log\max(|a|,b),
\tag{4}
\]

where `a,b` are coprime and `b>0`, and `h(X(O))=h(infinity)=0`.
The same normalization over `Qbar(t)` replaces rational logarithmic height
by rational-function degree. Let `lambda=hat h_{Qbar(t)}(P)`.

For the generic cubic,

\[
\Delta=64(K^2-1)^2,\qquad c_4=16(K^2+3),\qquad
j_E=\frac{64(K^2+3)^3}{(K^2-1)^2}.                       \tag{5}
\]

This is nonconstant: it has a pole at `t=1`, where `K=1`. Thus the elliptic
curve over `Qbar(t)` is nonisotrivial.

Also `P` is not torsion. Its admissible smooth specialization at `t=5/4` has `k=49/41`.
Lemma 3 of [the universal-contact note](square_sum_ppt_bridges_20261019.md)
proves that `(0,k)` has infinite order for a reduced odd-over-odd positive
rational `k!=1`; its proof uses the unbounded negative 2-adic valuations of
successive doubles on an integral model. A generic torsion relation would
specialize to the same torsion relation on this smooth fiber. This is
impossible. Here `f(5/4)=720/1681<1/2`, so this test specialization also lies in
the established admissible range.

Consequently

\[
\boxed{\lambda>0.}                                      \tag{6}
\]

For an explicit primary-source justification of this standard function-field
positivity step, apply Theorem 2 of Silverman's
[*A Lehmer-Type Lower Bound for the Canonical Height on Elliptic Curves over
Function Fields*](https://arxiv.org/pdf/2402.14771), pp. 1–2, with algebraically
closed constant field `Qbar`, function field `F=Qbar(t)`, extension degree
`D=1`, and the non-torsion point `P`. Its nonisotriviality hypothesis is exactly
the nonconstancy in (5); its positive lower bound applies. An optimal or
numerical value of `lambda` is unnecessary here. Only positivity is imported;
it is unchanged by the customary factor-of-two normalization choices for
canonical height.

## 3. Specialization is uniform even when the real parameter approaches a cusp

Silverman's height-specialization theorem gives

\[
\lim_{\substack{t\in\mathbb Q,\ E_t\text{ smooth}\\h(t)\to\infty}}
 \frac{\widehat h_t(P_t)}{h(t)}=\lambda.                 \tag{7}
\]

The base is the fixed smooth projective curve `P^1`, with its degree-one
height `h(t)`; `P` is a fixed rational section. Remove the finitely many
singular fibers and coefficient poles to obtain the smooth family. The limit
ranges over all remaining high-height parameters, not merely parameters in a
fixed real compact subset. Therefore tending to a removed boundary point in
the real topology does not invalidate (7).

The original reference is Silverman, *Heights and the specialization map for
families of abelian varieties*, J. reine angew. Math. 342 (1983), 197–211,
[Theorem B](https://doi.org/10.1515/crll.1983.342.197). The same author states the
height limit explicitly in his
[2006 MSRI height lectures, slide 22 (PDF p. 27)](https://legacy.slmath.org/attachments/workshops/301/HtSurveyMSRIJan06.pdf#page=27).
A primary-source proof in a more general setting is Call–Silverman,
[*Canonical heights on varieties with morphisms*, Theorem 4.1, p. 183](https://www.numdam.org/item/CM_1993__89_2_163_0.pdf#page=22).

For clarity, that latter theorem also applies directly to the degree-four
rational map of `P^1` induced by duplication on this family, with the
constant section `X(P)=0` and the divisor `O(1)`. Its canonical height is
`2*hat h_t(P_t)` under (4); the same factor occurs generically and cancels.
The base-height normalization is degree one in both applications.

Thus, for `a=lambda/2>0` and all sufficiently large `h(t)`,

\[
\widehat h_t(P_t)\ge a h(t).                            \tag{8}
\]

On the two spines, respectively,

\[
h(t_j^L)=\log(j+2),\qquad h(t_j^R)=\log(2j+2).           \tag{9}
\]

In particular both heights tend to infinity although their real limits are
`1` and `infinity`. No statement of uniformity in a varying section `qP+Q`
is being inferred from (7): only the fixed section `P` uses specialization.

## 4. An elementary, uniform, one-sided height comparison

Only an upper bound for canonical height in terms of naive height is needed.
The duplication map is

\[
X(2R)=\frac{X^4+2X^2-8KX+4K^2+1}{4(X-K)(X^2-1)}.        \tag{10}
\]

Write `K=a_0/b_0` and `X=r/s` in reduced form. Its homogeneous numerator and
denominator, after clearing the common coefficient denominator, are

\[
F=b_0^2r^4+2b_0^2r^2s^2-8a_0b_0rs^3
                         +(4a_0^2+b_0^2)s^4,
\]
\[
G=4b_0s(b_0r-a_0s)(r^2-s^2).
\]

If `M=max(|a_0|,b_0)` and `H=max(|r|,|s|)`, then
`max(|F|,|G|)<=16 M^2 H^4`. These forms have no simultaneous zero when
`K^2!=1`; at the finite denominator roots `1,-1,K`, the numerator values are
`4(K-1)^2,4(K+1)^2,(K^2-1)^2`, and at infinity it is nonzero.
Reduction can only lower height. Thus for every point, including the
projective exceptional cases,

\[
h(X(2R))\le4h(X(R))+2h(K)+\log16.
\]

Iteration and (4) give

\[
\boxed{\ \widehat h_t(R)\le\tfrac12 h(X(R))+
                              \frac{2h(K(t))+\log16}{6}.\ }\tag{11}
\]

This proof never divides by a discriminant or assumes the real
coefficients stay away from the singular locus.
Since `h(k(t))<=2h(t)+log4`, we have `h(K(t))<=4h(t)+log16`.
Also `h(X(Q_t))<=2h(t)+log6`. Equation (11) therefore yields, for `h(t)>=1`,

\[
\widehat h_t(Q_t)\le\frac73h(t)+\tfrac12\log96\le5h(t). \tag{12}
\]

## 5. Large indices force large point height

The square root of the Neron–Tate height is a norm modulo torsion.
Its reverse triangle inequality, (8), and (12) imply

\[
\sqrt{\widehat h_t(Q_t+qP_t)}
\ge \max\{0,\ |q|\sqrt{a h(t)}-\sqrt{5h(t)}\}.
\tag{13}
\]

The positive part matters: one must not square a possibly negative
one-sided lower bound. For `|q|>=2sqrt(5/a)`, (13) gives

\[
\widehat h_t(Q_t+qP_t)\ge\frac a4 q^2h(t).              \tag{14}
\]

Using (11) in the opposite direction now gives

\[
h(X(Q_t+qP_t))\ge\frac a2q^2h(t)-O(h(t)+1),            \tag{15}
\]

uniformly in both `t` and `q` within this range.

## 6. Transfer to the primitive hypotenuse, in both accepted charts

A successful point corresponds to a finite nonzero quartic parameter `u`
and a rational `z`, on

\[
z^2=K(t)(u^2+1)^2-4u(u^2-1).
\tag{16}
\]

A polynomial of fixed degree has height at most its degree times the
argument height, plus its coefficient height and an absolute constant.
Here the coefficients have height `O(h(t)+1)`, so (16) gives

\[
h(z)\le2h(u)+O(h(t)+1).
\]

The forward map, valid for every successful `u!=0`, is

\[
X=\frac{k(t)z+K(t)+2u+K(t)u^2}{2u^2}.
\]

Elementary rational-height inequalities then give absolute constants
`A,B>0` with

\[
h(X)\le A h(u)+B(h(t)+1).                              \tag{17}
\]

This is the correct direction for a lower bound on output height; a bound
for the inverse map alone would not suffice. It also avoids any issue at
exceptional inverse-chart denominators. The sign of `z` is immaterial.

For the accepted positive parameter `v`, we have exactly `h(v)=h(u)` in
both charts: inversion and sign do not change rational height. Write
`v=n/m>1` with positive coprime integers `n,m`. Primitive normalization is

\[
C'=\frac{n^2+m^2}{d},\qquad
 d=2\ \text{if }n,m\text{ are both odd},\quad d=1\ \text{otherwise}.
\]

Consequently

\[
\boxed{\ |\log C'-2h(u)|\le\log2.\ }                    \tag{18}
\]

Thus neither fraction reduction, the alternate chart, the factor-two
primitive normalization, nor swapping the legs can erase a cubic lower
bound. Combining (15), (17), and (18) gives

\[
\log C'\ge c_2q^2h(t)-C_2(h(t)+1)                     \tag{19}
\]

for all sufficiently high-height parameters and all `|q|` above a fixed
constant.

## 7. Finish with the already proved index bounds

The [cusp lower-bound theorem](square_sum_contact_cusp_lower_bound_20261019.md)
proves, for every successful signed index in either accepted chart,

\[
|q|\ge\frac{1}{2\log(1+\sqrt2)}\log j-O(1)
\]

along both spines. In particular the hypotheses on `|q|` in (19) eventually
hold for **all** successes. The additive term in (19) is then absorbed into
the positive quadratic term, and (9) proves (1).

The inherited inverse-chart estimate applies to every success, including
any least-index success. Indeed a successful quartic parameter is finite
and nonzero, so its forward image is finite. Substitution in the forward
map and quartic equation shows that `X=K` would force
`(K^2-1)u^2=0`, which is impossible. Similarly `X=0` forces
`u=-1/(2K)`. Its normalized parameter is `v=2K` in `(2,3)`, and

\[
f(v)-\frac{24}{25}
=\frac{-4(v-3)(v-2)(2v+1)(3v+1)}{25(v^2+1)^2}\ge0.
\]

That spread is inadmissible. Thus none of the inverse-map denominators
excluded in the output-height note vanishes at a successful point.

The recurrence theorem supplies a successful index `|q|=O(log j)`. The
[output-height theorem](square_sum_contact_output_height_20261019.md) bounds
its output length by `O((log j)^3)`. This proves the matching upper bound in
(2). For any least-absolute-index success, the same index bound applies;
the point-height and inverse-map estimates in that note work in the second
chart as well because `h(-1/u)=h(u)`. This proves (3).

## 8. Scope and verification

The specialization limit is an established theorem, not an elliptic Lehmer
conjecture over number fields. Generic positivity is applied over a
function field and is established independently of the cusp analysis.
The argument does not require an effective specialization threshold, an
explicit canonical-height computation, or a numerical asymptotic fit.

The companion checker verifies the family identities, nonconstant
`j`-invariant, homogeneous duplication identity, exceptional-root values,
chart symmetry, and primitive-height normalization using exact arithmetic.
It is a finite algebraic regression check; the cited specialization theorem
and the uniform argument above supply the asymptotic proof. Run:

```bash
python3 square_sum_contact_sharp_output_height_20261019.py controls.json
python3 -O square_sum_contact_sharp_output_height_20261019.py controls.optimized.json
```
