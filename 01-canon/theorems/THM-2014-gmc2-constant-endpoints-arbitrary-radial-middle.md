---
id: THM-2014
title: "GMC(2) on the infinite radial slice aZ+b(ZW)+cW"
status: PROVED
source: codex-2026-07-21-S83
depends_on:
  - THM-1510  # circular-Gaussian charge/radial normal form
  - THM-1660  # constant-middle Hermite base case
related:
  - THM-1515  # broader {-1,0,1} claim; its published dominance proof is not used here
  - THM-1835  # Bessel/I0 identity; its abstract A2 does not supply this asymptotic
  - HYP-8765 # radial-channel return tower
---

# THM-2014 — constant charged endpoints and an arbitrary radial middle

Let (Z=(X+iY)/\sqrt2), (W=\overline Z), and (U=ZW), so

\[
\mathbb E[Z^rW^s]=\delta_{rs}r!,\qquad
L(f):=\mathbb E[f(U)]=\int_0^\infty f(u)e^{-u}\,du.
\]

For arbitrary (a,c\in\mathbb C) and (b\in\mathbb C[U]), put

\[
P=aZ+b(U)+cW.
\]

Then

\[
\mathbb E[P^m]=0\quad(m\ge1)
\qquad\Longleftrightarrow\qquad
b=0\ \text{ and }\ ac=0.
\]

Consequently every moment-null polynomial on this infinite-dimensional slice is
charge-one-sided, and the Gaussian Moments conclusion holds there: for every fixed
(Q\in\mathbb C[Z,W]),
(mathbb E[QP^m]=0) for all sufficiently large (m).

## 1. Exact charge reduction

Write \(\gamma=ac\). Charge balance forces the same number (k) of (aZ) and (cW)
factors. Hence, exactly,

\[
M_m:=\mathbb E[P^m]
=\sum_{0\le k\le m/2}
  \frac{m!}{k!^2(m-2k)!}\,\gamma^k
  L\!\left(U^k b(U)^{m-2k}\right).                 \tag{1}
\]

This is the coefficient form of

\[
\sum_{m\ge0}M_m\frac{t^m}{m!}
=L\!\left(e^{tb(U)}I_0(2t\sqrt{\gamma U})\right),  \tag{2}
\]

used only as a formal identity below.
Here (I_0(2t\sqrt{\gamma U})=\sum_{k\ge0}(\gamma U)^kt^{2k}/(k!)^2),
so no choice of square-root branch is involved.

## 2. A uniform factorial lemma

Let (d=\deg b\ge1), let \(\beta\ne0\) be its leading coefficient, and let
(b_+(u)) be obtained by taking absolute values of its coefficients. There is a
constant (C_b), depending only on (b), such that for every (n,k\ge0),

\[
\int_0^\infty u^k b_+(u)^n e^{-u}\,du
\le C_b |\beta|^n(dn+k)!.                           \tag{3}
\]

To see this, write

\[
b_+(u)=|\beta|u^d q_+(1/u),\qquad
q_+(x)=1+q_1x+\cdots+q_dx^d,\quad q_j\ge0,
\]

and normalize the integral by (N!=(dn+k)!). It becomes an expectation of
(q_+(1/V)^n) for (V\sim\Gamma(N+1,1)). Choose (0<\varepsilon<1).
On (V\ge\varepsilon N),
(q_+(1/V)^n\le\exp(Cn/V)\le\exp(C/(\varepsilon d))).
On (1\le V<\varepsilon N), (q_+(1/V)^n\le B^n), while the Gamma lower-tail
bound is

\[
\Pr(V<\varepsilon N)
\le \exp[-N(\varepsilon-1-\log\varepsilon)+O(1)].
\]

Choose \(\varepsilon\) so that
(d(\varepsilon-1-\log\varepsilon)>\log B\); then this dominates (B^n),
since (N\ge dn). Finitely many small values of (N) are absorbed into (C_b).
Finally, on (0<V<1), (q_+(1/V)^n\le B^nV^{-dn}); after normalization the
contribution is at most (B^n/((k+1)N!)), uniformly bounded. This proves (3),
including (n=0).

The same Gamma normalization gives the sharp radial asymptotic. If

\[
b(u)=\beta u^d+b_{d-1}u^{d-1}+O(u^{d-2}),
\]

then

\[
\frac{L(b^m)}{\beta^m(dm)!}
\longrightarrow
\exp\!\left(\frac{b_{d-1}}{d\beta}\right).          \tag{4}
\]

Indeed, with (V_m\sim\Gamma(dm+1,1)), the left side is
(\mathbb E[q(1/V_m)^m]), where
(q(x)=1+(b_{d-1}/\beta)x+O(x^2)). Gamma concentration gives
(V_m/m\to d), hence
(m\log q(1/V_m)\to b_{d-1}/(d\beta)). With the same fixed
\(\varepsilon\), the two lower regions used for (3) have normalized (L^1)
mass tending exponentially to zero, while the integrand is uniformly bounded
on (V_m\ge\varepsilon dm). This supplies uniform integrability. In particular,
the limit in (4) is nonzero.

## 3. Degrees at least two: entropy loses to the factorial gap

Assume (d\ge2). In the (k)-th summand of (1), set

\[
n=m-2k,\qquad N=dn+k=dm-(2d-1)k.
\]

Since (0\le k\le m/2), (N\ge m/2). By (3),

\[
\left|L(U^kb^{m-2k})\right|
\le C_b|\beta|^{m-2k}N!.
\]

Moreover,

\[
\frac{m!}{(m-2k)!}\le m^{2k},\qquad
\frac{N!}{(dm)!}\le (m/2)^{-(2d-1)k}.
\]

Therefore all charged corrections, uniformly including (k\asymp m), satisfy

\[
\frac{|M_m-L(b^m)|}{|\beta|^m(dm)!}
\le C_b\sum_{k\ge1}
 \frac{\left(2^{2d-1}|\gamma/\beta^2|\,m^{-(2d-3)}\right)^k}{(k!)^2}
=O(m^{-(2d-3)}).                                    \tag{5}
\]

Combining (4) and (5),

\[
\boxed{
M_m=\beta^m(dm)!
\left[
\exp\!\left(\frac{b_{d-1}}{d\beta}\right)+o(1)
\right].}                                           \tag{6}
\]

Thus (M_m\ne0) for every sufficiently large (m) whenever \(\deg b\ge2\).
This is a genuine uniform separation theorem: the (m^{2k}) multinomial entropy
is paid for by a factorial loss of (m^{(2d-1)k}), leaving the decisive
(m^{-(2d-3)k}).

## 4. Degrees zero and one

If (b=b_0), direct Gaussian integration gives

\[
F(t):=\sum_{m\ge0}M_m\frac{t^m}{m!}
=\exp(b_0t+\gamma t^2).
\]

Thus (F\equiv1) forces (b_0=\gamma=0).

If (b=b_0+b_1U), (2) and the Laplace transform of (I_0) give

\[
F(t)=
\frac{\exp\!\left(b_0t+\gamma t^2/(1-b_1t)\right)}{1-b_1t}.              \tag{7}
\]

If all positive moments vanish, then \(\log F=0\). Its first three coefficients are

\[
b_0+b_1,\qquad \frac{b_1^2}{2}+\gamma,\qquad
\frac{b_1^3}{3}+\gamma b_1.
\]

The last two imply \(-b_1^3/6=0\), so (b_1=\gamma=b_0=0).

Together with (6), this proves the nullcone statement in every radial degree.
If (b=0) and (ac=0), then (P=aZ) or (P=cW); charge additivity makes
(QP^m) have nonzero charge once (m) exceeds the finite charge range of (Q).
That proves the GMC conclusion. ∎

## Scope and consequence for the full problem

This theorem allows an **arbitrary** charge-zero radial polynomial, but the two charged
coefficients are constants. It does not prove the full \(\{-1,0,1\}\) stratum with
radial-polynomial charged coefficients (A(U)W+C(U)Z). For that generalization the
endpoint term is (H(U)=UA(U)C(U)), and the entropy boundary occurs when
(\deg H\) lies within two of (2\deg b). That boundary, and collisions among several
charge circuits, are the remaining radial-channel problem recorded in HYP-8765.

The proof deliberately does **not** use the invalid rule “the unique highest factorial
term cannot cancel.” What matters is the quantitative gap after multinomial entropy is
included; (5) is the missing uniform estimate.
