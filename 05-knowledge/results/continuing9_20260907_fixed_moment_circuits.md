# Every circuit sign word survives fixing the first two root moments

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**
For every degree \(d\ge3\), sum \(S>0\), and square sum
\(S^2/d<Q<S^2\), the Newton-circuit image of polynomials with distinct positive
root parameters contains an explicit open neighborhood of the all-one vector.
Consequently every ternary circuit sign word occurs at these same two moments.
The quotients themselves are constrained: Section 5 gives the complete cubic
image and its sharp variance threshold for an unbounded upper quotient.

This concerns positive root parameters and their first two moments. It does
not impose the Laurent model's C/D interlacers, factorial coefficient carrier,
original phase, or response sign. Those remain additional predicates.

## 1. Inheritance, interpretation, and literature credit

The closest proved mechanism is the exact coefficient inversion in
[continuing8 Newton universality](continuing8_20260906_newton_universality.md).
That result fixes the first moment while increasing the initial Newton ratio
to force separated roots. Its explicit second-moment obstruction is the
degree-five anchor \(S=13,Q=59\), where \(R_1=338/275\) rules out some prescribed
amplitudes. The canonical hostile \((1,1,3,3,8)\) from
**THM-3004**, file
`01-canon/theorems/THM-3004-circuit-sign-change-cluster-law-and-classifier-refutation.md`,
already defeats a two-end sign classifier. The corrected near miss is the
norm/alternation implication repaired in
[continuing8 ballot repair](continuing8_20260906_ballot_repair.md).
The least-used sidecar here is the classical Gaussian multiplier sequence:
it supplies an interior point on *every* fixed-second-moment slice.

The concept board was: exact coefficient inversion; narrow clusters;
Gaussian multiplier/Jensen polynomials; fixed moments; discriminant boundaries;
and the original coefficient consumer. The cheap hostile was to demand only
strict Newton inequalities at fixed moments. It fails already in degree three,
and remains false in degree four (Section 6). The useful replacement separates
sign flexibility from amplitude constraints.

The root mechanism is classical, not a discovery claim. Iserles, Nørsett and
Saff, *On transformations and zeros of polynomials*, Rocky Mountain J. Math.
21 (1991), 331–357, explicitly recall the multiplier sequence
\(\{r^{k^2}\}\), \(0<r<1\), in Section 3, printed p. 342.
[Primary author-hosted paper](https://www.math.vanderbilt.edu/saffeb/texts/107.pdf).
Our multiplier \(q^{k(k-1)/2}\) differs by a positive geometric factor after
putting \(r=\sqrt q\). The corresponding entire generating function is the
deformed exponential; Wang and Zhang give its definition and earlier root
history in Section 1, printed p. 1, of
[Zeros of the deformed exponential function](https://arxiv.org/pdf/1709.04357).
The finite proof below is self-contained and also proves the separation
needed here. The new repository consumer is the fixed-moment circuit box,
not priority for this real-rooted family or for cubic discriminants.

Write
\[
 N(n)=\prod_{i=1}^d(n+r_i)=\sum_{k=0}^d e_k n^{d-k},\qquad
 h_k=\frac{e_k}{\binom dk},\quad
 R_k=\frac{h_k^2}{h_{k-1}h_{k+1}},\quad
 C_k=\frac{R_k}{R_{k-1}}.
\]
Here \(r_i>0\), \(1\le k\le d-1\) for \(R_k\), and \(2\le k\le d-1\)
for \(C_k\). The circuit sign is \(\operatorname{sgn}(C_k-1)\).
Scaling all \(r_i\) preserves both \(R\) and \(C\).
Since \(e_1=S\) and \(e_2=(S^2-Q)/2\),
\[
 R_1=\frac{S^2(d-1)}{d(S^2-Q)}.
 \tag{1}
\]
The range \(S^2/d<Q<S^2\) is necessary for distinct positive parameters;
the construction below proves its sufficiency in every degree \(d\ge2\).

## 2. A separated center on every fixed-moment slice

For \(0<q<1\) put
\[
 P_d(z;q)=\sum_{k=0}^d\binom dk q^{\binom k2}z^k.
 \tag{2}
\]
**Lemma.** The roots of \(P_d\) are simple and negative. If the positive roots
of \(F_d(x)=P_d(-x;q)\) are \(a_1<\cdots<a_d\), then
\[
 a_{i+1}>a_i/q.
 \tag{3}
\]

**Proof.** Pascal's identity gives
\[
 P_{d+1}(z)=P_d(z)+zP_d(qz),\qquad
 F_{d+1}(x)=F_d(x)-xF_d(qx).
 \tag{4}
\]
The assertion starts at \(F_1(x)=1-x\). Assume it in degree \(d\).
At \(x=a_i\), the point \(qa_i\) lies strictly between \(a_{i-1}\) and
\(a_i\), with \(a_0=0\), so \(F_{d+1}(a_i)\) has sign \((-1)^i\).
At \(x=a_i/q\), the second term of (4) vanishes and the first has sign
\((-1)^i\): for \(i<d\) this point is between \(a_i\) and \(a_{i+1}\),
and for \(i=d\) it is after the last root. Therefore the intermediate value
theorem and the leading coefficient provide one root in each of
\[
 (0,a_1),\quad (a_{i-1}/q,a_i)\ (2\le i\le d),\quad (a_d/q,\infty).
 \tag{5}
\]
These \(d+1\) disjoint intervals account for the entire degree, proving
simplicity and excluding additional nonreal roots. If these roots are
\(b_1<\cdots<b_{d+1}\), then \(b_i<a_i\) and \(b_{i+1}>a_i/q\),
which proves (3) at the next degree. This also covers the first and last
intervals. \(\square\)

Now let
\[
 q=\frac{d(S^2-Q)}{(d-1)S^2},\qquad \mu=S/d.
 \tag{6}
\]
Then \(0<q<1\) and
\[
 N_0(n)=n^dP_d(\mu/n;q)
 \tag{7}
\]
has \(d\) distinct negative roots. Its positive parameters are the reciprocals
of the positive roots of \(P_d(-x)\), multiplied by \(\mu\), in reverse order.
Its coefficients are \(e_k=\binom dk\mu^kq^{\binom k2}\), hence
\[
 \sum r_i=S,\quad \sum r_i^2=Q,\quad R_k=1/q,\quad C_k=1.
 \tag{8}
\]
At \(q=1\) the center becomes \((n+\mu)^d\), and the zero-variance slice has
only that repeated-root point. At \(q=0\) the positive degree-\(d\) construction
degenerates; \(Q=S^2\) is unattainable by \(d\ge2\) positive parameters.
For \(q>1\), the first Newton inequality fails.

## 3. Explicit coefficient chart and a whole-box radius

Normalize \(\mu=1\); scaling at the end recovers any \(S\).
For prescribed positive \(c_2,\ldots,c_{d-1}\), define
\[
 h_k=q^{\binom k2}
       \prod_{j=2}^{k-1}c_j^{-\binom{k-j+1}{2}},\qquad
 a_k=\binom dk h_k,\qquad P_c(z)=\sum_{k=0}^d a_kz^k.
 \tag{9}
\]
An empty product is one. Thus \(a_0=1,a_1=d,a_2=\binom d2q\).
Taking consecutive second differences of \(\log h_k\), or using the
recursion \(h_{k+1}=h_k^2/(R_kh_{k-1})\), proves
\[
 R_k=q^{-1}\prod_{j=2}^k c_j,\qquad C_k=c_k.
 \tag{10}
\]
In particular the first two moments remain exactly fixed, not approximately
fixed. These formulas give mutually inverse positive coordinate maps between
the tail coefficients \(a_3,\ldots,a_d\) and the circuits.

Choose positive rational points \(0<x_1<\cdots<x_d\) such that
\[
 (-1)^\ell P_d(-x_\ell;q)>0\quad(1\le\ell\le d).
 \tag{11}
\]
They exist by Section 2 and density of the rationals. For rational \(q\)
they can be found and certified by exact root isolation. Put
\[
 E_k=\binom k3,\quad
 M_\ell=|P_d(-x_\ell;q)|,\quad
 W_\ell=\sum_{k=3}^d E_k\binom dk q^{\binom k2}x_\ell^k,
\]
\[
 b=\min\left\{\frac12,\frac1{2E_d},
                   \min_\ell\frac{M_\ell}{4W_\ell}\right\}>0.
 \tag{12}
\]
Any \(0<\delta\le b\) is a valid radius.

**Theorem (fixed-moment box).** If \(|c_j-1|\le\delta\) for every \(j\), then
\(n^dP_c(\mu/n)\) has distinct negative roots, the prescribed \(S,Q\), and
exactly the prescribed \(C_j=c_j\).

**Proof.** The total exponent in the product at coefficient \(k\) is
\[
 \sum_{j=2}^{k-1}\binom{k-j+1}{2}=\binom k3=E_k.
\]
Writing \(a_k^0=\binom dkq^{\binom k2}\), we have
\[
 (1+\delta)^{-E_k}\le a_k/a_k^0\le(1-\delta)^{-E_k}.
 \tag{13}
\]
The larger of the two deviations from one is the upper one; for example,
\((1-\delta)^{-E}+(1+\delta)^{-E}\ge2\) by convexity when \(E\ge1\).
Bernoulli's inequality and \(E_k\delta\le1/2\) give
\[
 |a_k/a_k^0-1|\le(1-\delta)^{-E_k}-1
 \le\frac{E_k\delta}{1-E_k\delta}\le2E_k\delta.
 \tag{14}
\]
Coefficients \(k=0,1,2\) do not change. Therefore
\[
 |P_c(-x_\ell)-P_d(-x_\ell)|\le2\delta W_\ell\le M_\ell/2.
 \tag{15}
\]
All alternating signs (11) persist. Together with \(P_c(0)=1\), they give
\(d\) distinct roots on the negative axis; the degree excludes any additional
or multiple roots. Equations (9)–(10) prove the remaining claims. \(\square\)

For real \(S,Q\), (12) is an explicit positive real bound. For rational
\(S,Q\), all samples, \(b\), and a chosen smaller dyadic \(\delta\) can be
rational; rational target circuits then give rational coefficients.
There is no rational-root claim.

**Corollary.** For every ternary word
\(\sigma\in\{-1,0,1\}^{d-2}\), choose \(c_j=1+\delta\sigma_j\).
This realizes precisely that circuit sign word at the fixed first two moments,
with all root parameters distinct and positive.

The degree-five anchor \(S=13,Q=59\) has \(q=275/338\). The exact certificate
below gives \(\delta=1/2048\): its entire closed three-dimensional circuit
cube is feasible as positive root geometry. In particular all 27 ternary
words occur. This is a whole-box consequence of (14)–(15), not an inference
from the 27 displayed examples. The C/D interlacer predicate is not retained.

More generally, since simple negative-root polynomials form an open set,
the same chart shows that the fixed-moment circuit image is open at every
simple positive-root point. The separated center proves it contains the
all-one vector on every strictly feasible slice.

## 4. What this settles and what it does not

The source is the classical multiplier/Jensen family (2). The target is the
Newton-circuit map of a fixed-degree positive-root polynomial. The map is
reciprocal reversal, scaling, and the exact coefficient chart (9). It
preserves the first two moments, root positivity and simplicity, and every
requested small circuit perturbation. The discarded data are the actual
factorial row, the Laurent interlacers, the original phase equation and carry
coefficients. Those require separate sidecars; no noncancellation theorem
is transported here.

The earlier amplitude obstruction \(C_2=1/2\) at \(S=13,Q=59\) survives:
\(R_2=R_1C_2=169/275<1\). It does not obstruct a negative first circuit.
Thus fixing two root moments cannot restrict which circuit *sign words*
occur, although it restricts their distances from the tie hyperplanes.
The exact magnitude image for \(d\ge4\), and its intersection with the
actual C/D interlacer domain, remain open in this report.

## 5. Complete cubic image and the variance threshold

Scale \(S\) to three, put
\[
 q=\frac{3(S^2-Q)}{2S^2},\quad
 v=1-q=\frac{3Q/S^2-1}{2},\quad a=\sqrt v\in(0,1).
 \tag{16}
\]
Let \(p\) denote the product of the three *scaled* positive parameters.
Their polynomial is \(t^3-3t^2+3(1-v)t-p\), with discriminant
\[
 \Delta=27\{4v^3-[p-(1-3v)]^2\}.
 \tag{17}
\]
Put
\[
 L=(1+a)^2(1-2a),\qquad U=(1-a)^2(1+2a).
 \tag{18}
\]
The real-root condition is \(L\le p\le U\). When the roots are real and
\(p>0\), they are all positive: a putative negative pair \(-u,-w\) and
positive third root \(3+u+w\) would make their pair sum
\(uw-(u+w)(3+u+w)<0\), contradicting \(3(1-v)>0\).
A single negative root contradicts \(p>0\), and zero roots are excluded.
Consequently distinct positive roots occur **iff**
\[
 \max(0,L)<p<U.
 \tag{19}
\]
At \(U\) the roots are \(1-a,1-a,1+2a\). At \(L\) they are
\(1+a,1+a,1-2a\); this endpoint is positive exactly when \(a<1/2\).

Since \(C_2=q^3/p\), the complete strict-root image is
\[
 \begin{cases}
 (q^3/U,\ q^3/L),&0<a<1/2,\\
 (q^3/U,\ \infty),&1/2\le a<1.
 \end{cases}
 \tag{20}
\]
If repeated positive parameters are allowed, include the finite endpoints;
the infinite endpoint and a zero-product endpoint are never included.
At zero variance \(a=0\), only \(C_2=1\) occurs with equal parameters.
The upper amplitude is bounded exactly when \(Q<S^2/2\); this is the sharp
cubic threshold, including equality on the unbounded side.

Both signs and the tie occur at every nonzero feasible variance, as also
follows from Section 3. More explicitly,
\[
 C_{\min}=\frac{(1-a)(1+a)^3}{1+2a},\quad
 C_{\min}-1=-\frac{a^3(2+a)}{1+2a}<0,
\]
\[
 C_{\max}=\frac{(1-a)^3(1+a)}{1-2a},\quad
 C_{\max}-1=\frac{a^3(2-a)}{1-2a}>0\quad(a<1/2).
 \tag{21}
\]
The necessary Newton bound \(C_2\ge q\) is strictly weaker:
\(C_{\min}/q=(1+a)^2/(1+2a)>1\).

## 6. Exact hostile controls and reproduction

The cubic
\[
 n^3+3n^2+\frac94n+\frac{135}{256}
 \tag{22}
\]
has \(S=3,Q=9/2\), \(R=(4/3,16/15)\), and \(C_2=4/5\).
Every Newton ratio is strictly above one, but
\(\Delta=-25515/65536<0\), and \(4/5<27/32=C_{\min}\).
Its first failed implication is “all strict Newton inequalities imply
real-rootedness.” Multiplying (22) by \(n+1\) gives a degree-four hostile
at \(S=4,Q=11/2\), with
\[
 e=(1,4,21/4,711/256,135/256),\quad
 R=(8/7,784/711,18723/17920)>1.
\]
It has exactly two negative real roots and a nonreal conjugate pair.
The same moment slice nevertheless contains the full circuit box about
one furnished by \(d=4,q=7/8,\delta=1/1024\).

The [standalone source](continuing9_20260907_fixed_moment_circuits.py),
[exact certificate](continuing9_20260907_fixed_moment_circuits_certificate.json),
[normal output](continuing9_20260907_fixed_moment_circuits.out), and
[optimized output](continuing9_20260907_fixed_moment_circuits_optimized.out)
are frozen. Run:

~~~text
python continuing9_20260907_fixed_moment_circuits.py
python -O continuing9_20260907_fixed_moment_circuits.py
~~~

The finite universe is six boxes
\((d,q)=(3,3/4),(4,7/8),(5,275/338),(6,1/2),(3,1/4),(4,9/10)\),
with all 132 ternary words. The source also checks the recurrence and
derivative coefficient identities through degree ten, an independent
recursive coefficient reconstruction through degree eight, exact cubic
discriminants and repeated endpoints, both hostile polynomials, and the
scaled degree-five anchor. The whole-box/all-degree results rest on the
proofs above. Exact root isolators and sign samples provide finite controls;
no floating root calculation is accepted.

Both raw subprocess stdout streams are identical LF bytes, with explicit
stdout newline configuration; no output normalization was used.
They pass **1,190 always-active gates**.

- Source SHA256: a71ee4509bb46644f8e1fcf3a5f6931941363876540ea8565d1addf46e6e8102.
- Output SHA256: 5d51e2a2aa33a1dd55d4d3caaa1c64702567ca4787e78a8b5d1e8c4f0486abfe.
- Certificate SHA256: 291b62a5638f8a057ccf514851cb0b5d100a2a61a3da8eff9151b49040be11b9.
