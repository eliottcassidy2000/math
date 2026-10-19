# Complete divergent clocks from finitely many Taylor coefficients

19 October 2026. **PROVED:** the local, full-passage theorems and corollaries proved here. **CITED:** the named classical integration, inverse-series, and analytic-dynamics inputs. **FINITE-EXACT:** the accompanying formal-algebra checks. Numerical orbit checks are supporting controls, not proofs or exact arithmetic certifications. No historical-priority claim is made.

## 0. Result and inheritance

Let an analytic increment have a positive even-order zero,

\[
 g(y)=a y^r+O(y^{r+1}),\qquad a>0,\quad r=2m\ge2.
\]

For sufficiently small fixed endpoints \(\ell<0<b\), start at \(y_0=\ell\) and iterate

\[
 y_{n+1}=y_n+\varepsilon+g(y_n),\qquad
 N_g(\varepsilon)=\min\{n\ge0:y_n\ge b\},\quad\varepsilon>0.
\]

The complete divergent part of the integer passage time is

\[
 \boxed{\quad N_g(\varepsilon)
   =\sum_{j=0}^{m-1} A_j(g)\,
       \varepsilon^{-(r-1-2j)/r}+O(1),\quad}
 \tag{1}
\]

where

\[
 A_j(g)=\frac{2\pi}{r}\csc\!\left(\frac{(2j+1)\pi}{r}\right)
   [y^{2j}]\left(\frac{g(y)}{y^r}\right)^{-(2j+1)/r}.
 \tag{2}
\]

The power uses the analytic branch positive at zero; \([y^k]\) denotes an ordinary Taylor coefficient, without a factorial. There is **no logarithmically divergent term** in this fixed-endpoint, two-sided passage. The coefficients depend only on the Taylor jet of \(g\) through degree \(2r-2\).

Consequently, two such germs with the same jet through degree \(2r-2\) have passage counts differing by \(O(1)\), even when one discrete germ embeds in an analytic flow and the other does not. More precisely, equality of the \(m\) coefficients in (2) is necessary and sufficient for bounded count difference; matching all jet coefficients is stronger than necessary.

These coefficients classify the specified **constant-additive parameter families in the specified coordinate**. They are not asserted to be coordinate-conjugacy invariants of the unperturbed germ. A nonlinear coordinate change also transforms the parameter unfolding; imposing a fresh constant-additive \(\varepsilon\) after that change generally defines a different family.

The closest established repository mechanism is the corrected real clock in
[`higher_order_bottleneck_clocks_20261019/HIGHER_ORDER_BOTTLENECKS.md`](../higher_order_bottleneck_clocks_20261019/HIGHER_ORDER_BOTTLENECKS.md), which treats pure even powers and \(y^4+\lambda y^6\). The present result retains the same summable discrete defect and makes the entire finite coefficient list explicit. The complementary obstruction is the distinction between a formal logarithm and an analytic first integral in
[`parabolic_cotangent_integrability_20261019/parabolic_cotangent_invariant_barrier.md`](../parabolic_cotangent_integrability_20261019/parabolic_cotangent_invariant_barrier.md). The relevant hostile tests are an asymmetric quartic, a perturbation in the last necessary jet degree, and a one-sided passage where a logarithm really does survive.

## 1. Exact quantifiers and domain control

**Theorem 1 (full-passage clock).** Suppose \(g\) is real analytic in a neighborhood of zero and has the expansion above. There is \(\delta>0\) such that for every fixed

\[
 -\delta\le\ell<0<b\le\delta
\]

there are constants \(\varepsilon_0>0\) and \(E<\infty\) for which the iterates through the first crossing of \(b\) are defined and

\[
 \left|N_g(\varepsilon)-
 \sum_{j=0}^{m-1}A_j(g)\varepsilon^{-(r-1-2j)/r}\right|\le E
 \qquad(0<\varepsilon\le\varepsilon_0).
 \tag{3}
\]

The constants may depend on \(g,\ell,b\), and \(r\), but not on \(\varepsilon\). The coefficients \(A_j\) do not depend on the fixed endpoints. No claim is made about an arbitrary orbit entering this interval, moving endpoints, an \(\varepsilon\)-dependent germ, or a uniform bound over all analytic germs with a given jet.

Choose \(\delta\) sufficiently small that \(g\) is analytic on a neighborhood of
\(J=[-2\delta,2\delta]\), that \(g(y)>0\) for \(y\ne0\) in \(J\), and that

\[
 \sup_J |g'|\le\eta<1,\qquad\eta\le\tfrac14.
 \tag{4}
\]

Write \(f_\varepsilon=\varepsilon+g\), and take
\(0<\varepsilon_0\le\delta/4\). Since \(g(0)=0\), at any pre-exit point
\(y\in[\ell,b)\subset[-\delta,\delta]\),

\[
 0<f_\varepsilon(y)\le\varepsilon_0+\eta\delta\le\delta/2.
\]

Thus every step used before stopping stays in \(J\), and

\[
 b\le y_N\le 3\delta/2.
 \tag{5}
\]

The passage is finite because every increment is at least \(\varepsilon\): if it had not exited, \(y_n\ge\ell+n\varepsilon\) would eventually exceed \(b\). In particular, an overshoot never calls for an unproved continuation of the germ.

There are finite constants \(B,C\), independent of \(\varepsilon\), such that on \(J\)

\[
 |f_\varepsilon''|\le B,\qquad
 (f_\varepsilon')^2\le C f_\varepsilon.
 \tag{6}
\]

For the second assertion,
\(g'(y)^2/g(y)\sim a r^2 |y|^{r-2}\) as \(y\to0\), so the quotient extends continuously and boundedly at zero. Since \(g\le f_\varepsilon\), (6) follows. This includes the quadratic case, where the limiting quotient is \(4a\).

## 2. A uniform discrete-versus-integral estimate

**Lemma 2 (summable corrected-clock defect).** Let \(f>0\) be \(C^2\) on an interval containing a whole step \([y,y+f(y)]\), and suppose throughout that step

\[
 |f'|\le\eta<1,\quad |f''|\le B,\quad (f')^2\le Cf.
\]

Set

\[
 \Psi(y)=\int_{0}^{y}\frac{ds}{f(s)}+\tfrac12\log f(y)
\]

whenever the primitive is defined on the interval in question (any basepoint in the interval gives the same increment). Then

\[
 |\Psi(y+f(y))-\Psi(y)-1|\le D f(y),\qquad
 D=\frac{5B}{12}+\frac{7C(1+\eta)}{12(1-\eta)}.
 \tag{7}
\]

**Proof.** Put \(F=f(y)\), \(A=f'(y)\), and
\(q(t)=f(y+tF)/F=1+s(t)\) for \(0\le t\le1\). The derivative bound gives

\[
 |s(t)|\le\eta t,\qquad1-\eta\le q(t)\le1+\eta.
\]

Taylor's theorem gives

\[
 |s(t)-tA|\le\tfrac12 BF t^2.
 \tag{8}
\]

Also \(q'(t)=f'(y+tF)\), so
\(|q'(t)|^2\le C f(y+tF)\le C(1+\eta)F\). Integrating from zero yields

\[
 |s(t)|^2\le C(1+\eta)F t^2.
 \tag{9}
\]

Use the exact identity \((1+s)^{-1}=1-s+s^2/(1+s)\). Equations (8)–(9) imply

\[
 \left|\int_0^1\frac{dt}{q(t)}-1+\frac A2\right|
 \le \frac{BF}{6}+\frac{C(1+\eta)F}{3(1-\eta)}.
\]

For real \(|s|\le\eta\),
\(|\log(1+s)-s|\le s^2/[2(1-\eta)]\). Therefore

\[
 \left|\tfrac12\log q(1)-\tfrac A2\right|
 \le \frac{BF}{4}+\frac{C(1+\eta)F}{4(1-\eta)}.
\]

Adding proves (7). The first-order \(A/2\) terms cancel. ∎

Apply the lemma to \(f_\varepsilon\). Since the increments telescope,

\[
 \sum_{n<N} f_\varepsilon(y_n)=y_N-\ell\le\tfrac52\delta.
\]

Hence

\[
 N=\int_\ell^{y_N}\frac{dy}{f_\varepsilon(y)}
   +\tfrac12\log\frac{f_\varepsilon(y_N)}{f_\varepsilon(\ell)}+O(1),
 \tag{10}
\]

with defect at most \(5D\delta/2\). The endpoint logarithm is uniformly bounded: \(\ell\ne0\), and \(y_N\in[b,3\delta/2]\) stays in a fixed compact interval away from zero. Both endpoint values of \(f_\varepsilon\) therefore stay between positive constants independent of \(\varepsilon\).

The overshoot integral is also bounded independently of \(\varepsilon\). It lies inside the final step, so the bound on \(q\) above gives

\[
 0\le\int_b^{y_N}\frac{dy}{f_\varepsilon(y)}
 \le\int_{y_{N-1}}^{y_N}\frac{dy}{f_\varepsilon(y)}
 \le\frac1{1-\eta}.
\]

We have proved the uniform reduction

\[
 \boxed{\quad N_g(\varepsilon)=
       \int_\ell^b\frac{dy}{\varepsilon+g(y)}+O(1).\quad}
 \tag{11}
\]

This is an estimate for the discrete recurrence, not an identification of that recurrence with an ODE or an exact analytic flow.

## 3. Exact flattening and all divergent terms

Let \(Q(y)=g(y)/y^r\), completed at zero by \(Q(0)=a\). The initial choice of \(\delta\) is also made sufficiently small that

\[
 u=U(y):=yQ(y)^{1/r}
 \tag{12}
\]

is a real-analytic, orientation-preserving diffeomorphism there. This is possible because \(U'(0)=a^{1/r}>0\). Let \(y=h(u)\) be its analytic inverse. By construction, **exactly** \(g(h(u))=u^r\); thus

\[
 \int_\ell^b\frac{dy}{\varepsilon+g(y)}
 =\int_{U(\ell)}^{U(b)}\frac{h'(u)}{\varepsilon+u^r}\,du.
 \tag{13}
\]

Choose a fixed \(d>0\) so small that \([-d,d]\) is contained in the transformed passage interval and in the domain of analyticity of \(h\). The substitution (12) is made inside the integral. It does not conjugate the discrete map to \(u\mapsto u+\varepsilon+u^r\): generally \(U(y+\varepsilon+g(y))\ne U(y)+\varepsilon+U(y)^r\).

The two outer portions of (13) are uniformly bounded because their distance from zero is positive. On the inner interval expand

\[
 h'(u)=\sum_{k=0}^{r-1}c_k u^k+u^rR(u),\qquad R\text{ bounded}.
 \tag{14}
\]

The last term contributes at most \(2d\sup|R|\), because \(r\) is even and
\(0\le u^r/(\varepsilon+u^r)\le1\). Every odd Taylor monomial integrates to zero on \([-d,d]\). In particular, the borderline degree \(r-1\), which could create a logarithm, is odd and cancels.

For \(k=2j\le r-2\), the full-line integral is convergent, and

\[
 \int_{\mathbb R}\frac{u^{2j}}{\varepsilon+u^r}\,du
 =\varepsilon^{-(r-1-2j)/r}
   \frac{2\pi}{r}\csc\!\left(\frac{(2j+1)\pi}{r}\right).
 \tag{15}
\]

This is the beta integral after scaling \(u=\varepsilon^{1/r}t\); its hypotheses are \(0<2j+1<r\). Truncating the full line to \([-d,d]\) costs at most

\[
 2\int_d^\infty u^{2j-r}\,du
   =\frac{2d^{2j-r+1}}{r-2j-1},
\]

a constant. Combining (11)–(15) proves (1), with \(A_j\) equal to the constant in (15) times \(c_{2j}\).

The no-log conclusion does **not** require symmetric original endpoints or an even germ. For example, the borderline monomial over the actual asymmetric transformed endpoints satisfies

\[
 \int_{U(\ell)}^{U(b)}\frac{u^{r-1}}{\varepsilon+u^r}\,du
 =\frac1r\log\frac{\varepsilon+U(b)^r}
                         {\varepsilon+U(\ell)^r}=O(1).
\]

The two logarithmic singularities cancel across zero. All other odd terms may likewise be canceled on the fixed symmetric inner interval, leaving only bounded outer portions.

## 4. Exact jet dependence, a complete criterion, and sharpness

Lagrange inversion applied to \(U(y)=yQ(y)^{1/r}\) gives

\[
 [u^n]h(u)=\frac1n[y^{n-1}]Q(y)^{-n/r}\quad(n\ge1).
\]

Differentiating coefficientwise gives the especially useful identity

\[
 \boxed{\quad c_k=[u^k]h'(u)=[y^k]Q(y)^{-(k+1)/r}.\quad}
 \tag{16}
\]

This proves (2). Computing \(c_{2j}\) uses only the coefficients of \(Q\) through degree \(2j\), equivalently the coefficients of \(g\) through degree \(r+2j\). The final possible divergent coefficient uses degree \(r+(r-2)=2r-2\). Terms of higher degree affect only the bounded remainder in the full passage.

**Corollary 3 (complete bounded-difference criterion for the additive families).** Let \(g_1,g_2\) have the same positive even order \(r\), with possibly different sufficiently small fixed negative and positive endpoints. Then

\[
 N_{g_1}(\varepsilon)-N_{g_2}(\varepsilon)=O(1)
 \quad\Longleftrightarrow\quad
 A_j(g_1)=A_j(g_2)\text{ for every }0\le j<m.
 \tag{17}
\]

The forward direction follows because any first unequal coefficient in the ordered, distinct divergent powers dominates all subsequent powers and the bounded remainder. The reverse direction is immediate from (1). In particular,

\[
 g_1-g_2=O(y^{2r-1})\quad\Longrightarrow\quad
 N_{g_1}-N_{g_2}=O(1).
 \tag{18}
\]

**Sharpness of the universal jet bound.** For \(r\ge4\), compare

\[
 g_0(y)=a y^r,\qquad g_\lambda(y)=a y^r+\lambda y^{2r-2},
\]

where \(\lambda\ne0\) is fixed and the neighborhood is small enough to preserve positivity. They agree through degree \(2r-3\). All clock coefficients before the last agree, but

\[
 c_{r-2}(g_\lambda)-c_{r-2}(g_0)
  =-\frac{r-1}{r}\lambda a^{-(2r-1)/r}.
\]

With \(K_r=(2\pi/r)\csc(\pi/r)\),

\[
 N_{g_\lambda}-N_{g_0}
 =-K_r\frac{r-1}{r}\lambda a^{-(2r-1)/r}
       \varepsilon^{-1/r}+O(1).
 \tag{19}
\]

Therefore the degree \(2r-2\) cannot be lowered in a theorem based solely on agreement of an arbitrary germ's full initial jet. For \(r=2\), the jet bound is degree two, and varying its leading coefficient changes the leading clock itself.

Matching the entire jet is not necessary. The quartic example in the next section exhibits genuinely different degree-six jets with the same divergent clock.

## 5. Quadratic, asymmetric quartic, and sextic formulas

For \(r=2\),

\[
 N_g(\varepsilon)=\frac{\pi}{\sqrt a}\varepsilon^{-1/2}+O(1).
 \tag{20}
\]

Thus arbitrary analytic cubic and higher terms leave the full-passage count unchanged up to a bounded error, provided \(a\) is unchanged.

For

\[
 g(y)=a y^4+b_5 y^5+b_6 y^6+O(y^7),
\]

formula (16) gives

\[
 c_0=a^{-1/4},\qquad
 c_2=a^{-3/4}\left(\frac{21b_5^2}{32a^2}-\frac{3b_6}{4a}\right).
\]

Hence

\[
 N_g(\varepsilon)=\frac\pi{\sqrt2}a^{-1/4}\varepsilon^{-3/4}
 +\frac\pi{\sqrt2}a^{-3/4}
   \left(\frac{21b_5^2}{32a^2}-\frac{3b_6}{4a}\right)
      \varepsilon^{-1/4}+O(1).
 \tag{21}
\]

At \(a=1,b_5=0,b_6=\lambda\), this recovers
\(-3\pi\sqrt2\lambda/8\) as the second coefficient. **Odd terms of the original germ cannot simply be discarded:** its \(y^5\) coefficient contributes quadratically to the divergent clock. The cancellation is of odd powers in \(h'(u)\), after exact flattening.

As a useful positive control, \(y^4\) and
\(y^4+y^5+(7/8)y^6\) have identical divergent clocks and therefore \(O(1)\)-close passage counts, despite different jets through degree six.

For a sextic, write

\[
 g(y)=a y^6(1+q_1y+q_2y^2+q_3y^3+q_4y^4+O(y^5)).
\]

Then

\[
 \begin{aligned}
 c_0&=a^{-1/6},\\
 c_2&=a^{-1/2}\left(\frac38q_1^2-\frac12q_2\right),\\
 c_4&=a^{-5/6}\left(
 \frac{21505}{31104}q_1^4-\frac{935}{432}q_1^2q_2
 +\frac{55}{36}q_1q_3+\frac{55}{72}q_2^2-\frac56q_4\right).
 \end{aligned}
\]

The corresponding clock is

\[
 N_g(\varepsilon)=\frac{2\pi}{3}c_0\varepsilon^{-5/6}
  +\frac\pi3c_2\varepsilon^{-1/2}
  +\frac{2\pi}{3}c_4\varepsilon^{-1/6}+O(1).
 \tag{22}
\]

## 6. Finite-jet clock accuracy does not decide analytic embedding

There is a concrete pair for every even \(r\ge2\). Fix \(a>0\) and put

\[
 F_{\rm emb}(y)=y\bigl(1-(r-1)a y^{r-1}\bigr)^{-1/(r-1)},\qquad
 F_{\rm poly}(y)=y+a y^r.
\]

The first germ is exactly the time-one map of the analytic vector field
\(a y^r\partial_y\). Its increment is

\[
 g_{\rm emb}(y)=F_{\rm emb}(y)-y
 =a y^r+\frac r2a^2y^{2r-1}+O(y^{3r-2}),
\]

whereas \(g_{\rm poly}=a y^r\). Thus (18) gives

\[
 N_{g_{\rm emb}}(\varepsilon)-N_{g_{\rm poly}}(\varepsilon)=O(1).
 \tag{23}
\]

At \(\varepsilon=0\), however, their normalized iterative logarithms have opposite convergence behavior. For \(F_{\rm emb}\), the logarithm is exactly \(a y^r\). For the nonlinear polynomial \(F_{\rm poly}\), it diverges. The latter is a **CITED** classical consequence recorded on page 3 of Aschenbrenner–Bergweiler: among globally meromorphic parabolic functions, only the parabolic Möbius germs have convergent iterative logarithm. The nonlinear polynomial here is not Möbius.

The cotangent-lift equivalence proved in the companion note then transfers this contrast: \(C_F(y,p)=(F(y),p/F'(y))\) has a nonconstant meromorphic first-integral germ exactly when the normalized iterative logarithm converges. The embedded example has the analytic invariant \(p a y^r\); the polynomial example has no nonconstant meromorphic first-integral germ.

This is an explicit information-loss statement. The finite coefficient vector determines every unbounded term of a full-passage clock, while omitting analytic embedding and cotangent-integrability information. The comparison concerns the \(\varepsilon=0\) germs; it does not claim that the additively perturbed map \(\varepsilon+F_{\rm emb}\) remains exactly a time-one flow of that same vector field.

## 7. Scaled integers and the hypotheses that matter

For an integer base \(B\ge2\), set \(\varepsilon=B^{-rn}\). Equation (1) becomes

\[
 N_g(B^{-rn})=\sum_{j=0}^{m-1} A_j(g) B^{(r-1-2j)n}+O(1).
 \tag{24}
\]

The complete divergent hierarchy has at most \(m\) levels. A quartic has levels \(B^{3n},B^n\); a sextic has \(B^{5n},B^{3n},B^n\). Omitting a nonzero coefficient loses an unbounded number of integer units. Retaining all of them gives bounded absolute uncertainty, but does not select a unique integer, determine residues, or imply anything about primes or arbitrary Collatz entry. Endpoint and discrete phase information remains in the bounded error.

The full-passage hypothesis is essential to the no-log conclusion. As a hostile control, take the pure quadratic increment, start **at zero**, and stop at a small fixed \(b>0\). The same summable-defect proof gives

\[
 N_{0\to b}(\varepsilon)
 =\frac{\pi}{2\sqrt\varepsilon}
   +\frac12\log(1/\varepsilon)+O(1).
 \tag{25}
\]

Indeed, its integral is \(\varepsilon^{-1/2}\arctan(b/\sqrt\varepsilon)
=\pi/(2\sqrt\varepsilon)+O(1)\), but the corrected-clock endpoint logarithm is
\(\tfrac12\log((\varepsilon+y_N^2)/\varepsilon)
=\tfrac12\log(1/\varepsilon)+O(1)\). This is why (1) must not be exported unchanged to a one-sided or moving-endpoint problem.

## 8. Attribution, dependencies, and reproducibility

The leading saddle-node bottleneck laws and near-parabolic clocks are established mathematics. Christian Kuehn's *Scaling of Saddle-Node Bifurcations: Degeneracies and Rapid Quantitative Changes*, Proposition 2.1, gives the higher even-order continuous-time exponent. The \(\pi\)-escape phenomenon in the Mandelbrot setting goes back to the literature discussed by Klebanoff and by Brockmoeller–Scherz–Srkalovic; the latter work develops a substantially broader parabolic-perturbation setting than the additive real family treated here. The current note provides a self-contained discrete error estimate and coefficient extraction for its explicit scope. It asserts no novelty or priority over that literature.

Primary/reference sources checked 19 October 2026:

- Christian Kuehn, Proposition 2.1: <https://arxiv.org/html/0807.1546>, original record <https://arxiv.org/abs/0807.1546>.
- Thies Brockmoeller, Oscar Scherz, Nedim Srkalovic, *Pi in the Mandelbrot set everywhere*: <https://arxiv.org/html/2505.07138v1>.
- NIST DLMF, Lagrange inversion, §1.10(vii): <https://dlmf.nist.gov/1.10#vii>.
- NIST DLMF, Euler beta integral, §5.12: <https://dlmf.nist.gov/5.12>.
- Matthias Aschenbrenner and Walter Bergweiler, *Julia's equation and differential transcendence*, page 3, paragraph before Theorem 1: <https://www.mat.univie.ac.at/~maschenbrenner/pdf/itlog-final.pdf>. The convergence/divergence assertion used here comes from that cited classical result, not from inferring divergence merely from Theorem 1's differential-transcendence statement.

Run `python3 verify_finite_jet_clocks.py` for exact coefficient identities and small numerical controls. Run with `python3 -O` as well: the verifier uses explicit checks, not removable Python assertions. The analytic proof establishes the quantified theorem; finitely many symbolic or numerical cases alone do not.
