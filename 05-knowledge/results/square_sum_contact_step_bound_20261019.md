# Elliptic contact search: no uniform finite translate budget

**Status: PROVED**, with exact bounded controls. This answers the next
construction question after [universal signed Pythagorean contacts](square_sum_ppt_bridges_20261019.md):
although every admissible primitive petal has infinitely many outgoing
`(+,-)` neighbors, the prescribed search through `Q+qP`, `q∈Z`, has no
uniform finite first-success bound. The failure persists when both real
positive-spread parameter charts are accepted. No novelty claim is made.

The geometric reason is explicit. At the degenerate parameter `k=1`, one
elliptic translation becomes the Berggren middle-branch map `u↦2+1/u`.
Its Pell iterates lie deep in the forbidden AP-positivity region. Thus
seeds near a parabolic cusp can defeat any fixed finite translate list.

## 1. Exact object and acceptance rule

Use the notation of the predecessor. For a rational anchor `t>1`, put

\[
f(t)=\frac{4t(t^2-1)}{(t^2+1)^2},\qquad 0<f(t)<\frac12,
\qquad k=\frac{t^2+2t-1}{t^2+1},\quad K=k^2.
\]

The prescribed curve, translating point, and anchor image are

\[
 E_k:Y^2=(X-K)(X^2-1),\qquad P=(0,k),
\]
\[
 Q=\left(\frac{t^2+4t+1}{t^2+1},
          \frac{8t(t+1)}{(t^2+1)^2}\right).
 \tag{1}
\]

For a translated point `R=Q+qP`, use the quartic inverse

\[
 u=\frac{kY+X-K}{X(X-K)},\qquad
 z=\frac{2u^2X-K-2u-Ku^2}{k}.
 \tag{2}
\]

Then

\[
 z^2=K(u^2+1)^2-4u(u^2-1),\qquad
 \left(\frac{z}{u^2+1}\right)^2=1+f(t)-f(u). \tag{3}
\]

There are two real positive-spread charts: `u>1` and `−1<u<0`.
The latter must not be discarded, because

\[
 f(-1/u)=f(u).
 \tag{4}
\]

After setting `t'=u` in the first chart and `t'=−1/u` in the second,
accept precisely `0<f(t')<1/2` and `f(t')≠f(t)`. This gives an outgoing
`(+,-)` neighbor on either side of the anchor spread. For reduced
`t'=n/m`, primitive normalization is

\[
 (n^2-m^2,2nm,n^2+m^2)/\delta,\quad
 \delta=1\text{ if }n+m\text{ is odd},\quad
 \delta=2\text{ if }n,m\text{ are both odd}. \tag{5}
\]

Reorder the legs to the odd/even convention afterward. The count restricted
to raw `u>1` is reported separately below; it is a different search.

## 2. Direct translation and the singular Berggren limit

For a quartic point `(u,z)`, put `a=K(u²−1)−u`. Direct substitution in
the elliptic group law gives

\[
 u(R+P)=\frac{a+kz}{u(u-2K)}=\frac{2Ku+1}{a-kz},
 \tag{6}
\]
\[
 u(R-P)=\frac{a-kz}{u(u-2K)}=\frac{2Ku+1}{a+kz}.
 \tag{7}
\]

The two charts agree because

\[
 a^2-k^2z^2=u(u-2K)(2Ku+1).
 \tag{8}
\]

These are identities on the ordinary charts; the alternate expression
handles the indicated cancellations. If `v=u(R+P)`, its new quartic
`z` coordinate is

\[
 z(R+P)=\frac{K(v^2-1)-v-uv(v-2K)}k; \tag{9}
\]

for `R−P`, the right side has the opposite sign. Consequently the limiting
quartic branch can be tracked, not just the parameter alone.

At `t→1+`, `k→1`, and the quartic degenerates into

\[
 z^2=(u^2-2u-1)^2.
\]

The anchor tends to `(u,z)=(1,2)` on the component
`z=−u²+2u+1`. On this component, (6) and (7) become

\[
 +P:G(u)=\frac1{u-2},\qquad
 -P:F(u)=2+\frac1u. \tag{10}
\]

This is exactly the middle Berggren branch `M` and its inverse in
[the AP classification](ap_square_pythagorean_20261019.md), Section 3.
The identification is at the singular limit. Nondegenerate elliptic
translation is not globally a Möbius transformation.

For completeness, the component maps to the nonsingular part of the nodal
cubic by

\[
 X=\frac{2u+1}{u^2},\qquad
 Y=\frac{(-u^2+2u+1)(u+1)}{u^3}. \tag{11}
\]

Starting at `u=1`, the limiting minus translates are

\[
 b_q=F^q(1)=3,7/3,17/7,41/17,\ldots,
 \qquad 2<b_q\le3. \tag{12}
\]

The plus translates are

\[
 a_1=-1,\qquad a_q=-1/b_{q-1}\quad(q\ge2). \tag{13}
\]

No `a_q` or `b_q` is a node or an omitted inverse-chart point; the only
rational zeros of `X` arise at `u=−1/2`, which these sequences never reach.
The rational recurrences have no zero denominators on these sequences.
Therefore each fixed translate has these limits as `t→1+`.

## 3. The unbounded-index theorem

**Theorem.** For every finite integer `B`, every sufficiently large integer
`m` gives an admissible primitive anchor with `t=(m+1)/m` for which no
translate `Q+qP`, `0<|q|≤B`, produces an admissible neighbor. This holds
even when negative-chart parameters are normalized by (4).

**Proof.** On `[2,3]`, `f≥24/25>1/2`: its derivative changes sign only at
`1+sqrt(2)` in that interval and `f(2)=f(3)=24/25`. Thus every `b_q`
has inadmissible spread. For `q≥2`, (4) and (13) give
`f(a_q)=f(b_(q−1))≥24/25`, so those translates also fail in a neighborhood
of the limit.

The remaining `q=+1` has limiting parameter `−1`, so its boundary sign
must be checked rather than inferred by continuity. Its exact parameter is

\[
 u(Q+P)=\frac{2t^3+t^2-2t+1}{t(t^2-4t+1)},
\]

and

\[
 u(Q+P)+1=
 \frac{(t-1)(3t^2-1)}{t(t^2-4t+1)}<0\quad(1<t<4/3).
\]

Hence `u(Q+P)<−1` and its spread is negative. For a fixed finite list of
`q`, intersect the finitely many neighborhoods just obtained. Every
sufficiently large `m` places `(m+1)/m` in that intersection. Such anchors
are admissible because `f((m+1)/m)→0+`. Consecutive parameters are coprime
and of opposite parity, so their PPT is primitive. Finally `q=0` is only
the anchor and is excluded by the distinct-spread rule. QED.

The predecessor proves eventual success for each fixed anchor, so the
least successful absolute index is well defined and unbounded. There is
no contradiction between these assertions.

An independent check at the opposite cusp `t→∞` gives, with Pell numbers
`p_0=0,p_1=1,p_(q+1)=2p_q+p_(q−1)`,

\[
 u(Q+qP)\longrightarrow p_{q+1}/p_q,\qquad
 u(Q-qP)\longrightarrow -p_{q-1}/p_q. \tag{14}
\]

All limiting spreads for positive `q`, except the second expression at
`q=1`, are

\[
 1-\frac1{(p_{q+1}^2+p_q^2)^2}\quad\text{or}\quad
 1-\frac1{(p_{q-1}^2+p_q^2)^2},
\]

and are at least `24/25`. The exceptional `Q−P` has negative spread for
`t>sqrt(3)`, from

\[
 u(Q-P)=\frac{t^2+4t+1}{t^3+2t^2+t-2}\in(0,1).
\]

Thus the opposite cusp independently defeats every fixed translate
budget. In particular, changing to the swapped-leg anchor parameter
`(t+1)/(t−1)` does not restore a uniform finite bound.

## 4. Complete finite controls and selected larger anchors

The companion standard-library script checks every one of the 528
admissible parity-oriented primitive PPTs with hypotenuse `c≤10000`.
All arithmetic uses `Fraction`. Each accepted neighbor is independently
checked on the cubic, the quartic, the contact equation, and the primitive
PPT normalization.

With both positive-spread charts accepted, the first absolute-index
histogram is

| First `|q|` | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Anchors | 65 | 98 | 282 | 6 | 14 | 21 | 9 | 3 | 25 | 5 |

The five anchors requiring `|q|=11` are
`(71,2520,2521)`, `(73,2664,2665)`, `(4899,140,4901)`,
`(5183,144,5185)`, and `(5475,148,5477)`.
Their first successful primitive neighbor hypotenuses have between 857
and 939 decimal digits, depending on the sign.

If raw `u>1` is required and negative charts are discarded, the maximum
is instead 22. Both `(73,2664,2665)` and `(5475,148,5477)` need that strict
budget. This larger number is an artifact of restricting the chart and
must not be advertised as the best all-chart bound.

Selected larger exact spine tests give:

| Anchor `t` | First `|q|` | First successful sign(s) | Neighbor hypotenuse digits |
|---|---:|---:|---:|
| 100 | 4 | +4 | 137 |
| 1000 | 4 | +4 | 202 |
| 10000 | 13 | +13 | 2734 |
| 10^6 | 8 | +8 | 1555 |
| 10^9 | 13 | +13 | 6113 |
| 10^12 | 17 | +17 | 13910 |
| 1000001/1000000 | 19 | −19 | 8925 |

Anchors `10^18`, `10^30`, `(10^12+1)/10^12`, and `(10^30+1)/10^30`
all fail every tested `0<|q|≤12`; no claim is made about their first success.
Six one-sided admissible rational approximants to the two cone walls,
with denominator grids `10^4,10^8,10^12`, all succeed at `|q|=4`.
No omitted inverse-chart point occurred anywhere in these finite probes.
No first-success witness in this chosen census needed the odd/odd divisor,
but (5) is implemented and checked separately on `t=9/7`, which normalizes
to the same primitive petal as `t=8`.

The script retains **46454 explicit checks**, including 50-step exact
Pell-limit controls. Ordinary and optimized Python runs produce identical
output. A separate symbolic checker independently verifies the translation
identities, small-multiple formulas, both limiting components, branch tracking,
and Pell algebra, with **233 exact checks**. Its algebra supplement also
proves that `|q|=1` always fails and that `q=+2` succeeds for every rational
anchor `24≤t≤32`; the latter is a local sufficient patch, not a uniform rule.

After placing the scripts in the experiments folder, reproduce with:

```bash
python3 04-computation/experiments/square_sum_contact_step_bound_20261019.py
python3 -O 04-computation/experiments/square_sum_contact_step_bound_20261019.py
python3 04-computation/experiments/square_sum_contact_step_algebra_20261019.py
```

The census uses only Python’s standard library; the symbolic checker requires
SymPy. The JSON census retains per-anchor first-success signs, chart choices,
normalization divisors, spread direction, and resulting hypotenuse digit counts.

## 5. Exact scope and stopping decision

This closes the proposed uniform small-multiplier route negatively. It
does not obstruct an anchor-dependent height or spread bound, another
starting point, another translating point, or a separate low-degree
construction. It does not supply a lower bound on the smallest possible
neighbor height, only on the index required by this particular orbit search.

The universal extension theorem and all its fan-gluing corollaries remain
intact. The finite `c≤10000` success budget is an implementation fact;
the no-uniform-budget conclusion follows from the degeneration proof.
A larger exhaustive census is unnecessary for deciding the stated question.
