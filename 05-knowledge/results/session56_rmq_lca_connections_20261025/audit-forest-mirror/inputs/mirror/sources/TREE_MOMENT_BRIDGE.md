# The tree moment measure, exact gamma bridges, and rational certificates

Session31, 2026-10-22. Status: **PROVED**, with exact rational finite checks and separately labelled high-precision diagnostics. No novelty claim is made for the Lambert-W Stieltjes representation, Nuttall–Bouwkamp integral, moment positivity, or classical Gauss/Padé mechanism. The main contribution of this note is a self-contained derivation and a careful, useful transport between the session's objects.

## 1. Main result and domains

Put

\[
 R(z)=-\frac{W_0(-z)}z=\sum_{m\ge0}a_mz^m,
 \qquad a_m=\frac{(m+1)^{m-1}}{m!},\qquad R(0)=1.
\]

Here \(W_0\) is the principal Lambert branch. The series has radius \(1/e\), whereas the analytic function is defined on \(\mathbb C\setminus[1/e,\infty)\). For \(0<v<\pi\), define

\[
 x(v)=e^{v\cot v}\frac{\sin v}{v},\qquad
 D(v)=v^2+(1-v\cot v)^2.
\]

Then \(x\) decreases bijectively from \(e\) to \(0\), and the probability measure

\[
 d\mu(x)=\rho(x)\,dx,\qquad \rho(x(v))=\frac v\pi
\tag{1}
\]

satisfies

\[
 \boxed{R(z)=\int_0^e\frac{d\mu(x)}{1-zx}},\qquad
 \boxed{M(s):=\int_0^e x^s\,d\mu(x)
       =\frac{(s+1)^{s-1}}{\Gamma(s+1)}}\quad(\Re s>-1).
\tag{2}
\]

Complex powers in \(M\)'s formula use the principal logarithm of \(s+1\); here \(\Re(s+1)>0\). In particular, \(M(m)=a_m\) and \(M(0)=M(1)=1\). This is a compactly supported Stieltjes moment sequence, or a Hausdorff moment sequence after scaling the support to \([0,1]\).

## 2. Boundary values, signs, and complete reconstruction

First,

\[
 \frac{x'(v)}{x(v)}=2\cot v-v\csc^2v-\frac1v
 =-\frac{D(v)}v<0.
\tag{3}
\]

This proves strict monotonicity without an assumed picture of the branch cut. Let

\[
 F(w)=\frac{R(1/w)}w=-W_0(-1/w),\qquad w\notin[0,e].
\]

For \(x=x(v)\), direct substitution gives

\[
 W_0(-1/x+i0)=-v\cot v+iv.
\]

Indeed, if \(u=-v\cot v\), then
\((u+iv)e^{u+iv}=-v\csc(v)e^{-v\cot v}=-1/x(v)\).
Since \(-1/(x+i0)=-1/x+i0\), the signs are

\[
 F(x+i0)=v\cot v-iv,\quad
 F(x-i0)=v\cot v+iv,\quad
 \rho(x)=-\frac1\pi\Im F(x+i0)=\frac v\pi>0.
\tag{4}
\]

A keyhole Cauchy integral around \([0,e]\) gives

\[
 F(w)=\frac1{2\pi i}\int_0^e
          \frac{F(x-i0)-F(x+i0)}{w-x}\,dx
      =\int_0^e\frac{\rho(x)}{w-x}\,dx.
\tag{5}
\]

There are no dropped endpoint contributions. At \(w=e\), \(F\) is bounded with a square-root correction; the small-circle contribution tends to zero. At \(w=0\), the principal-W large-argument expansion gives \(F(w)=O(\log(1/|w|))\); the small-circle contribution is \(O(\epsilon\log(1/\epsilon))\). At infinity, \(F(w)=1/w+O(1/w^2)\), so the outer contribution in the Cauchy reconstruction vanishes. These estimates also exclude hidden endpoint atoms. Comparing coefficients of \(1/w\) in (5) proves \(\int\rho=1\). Substituting \(w=1/z\) proves (2)'s first formula, with its removable value at zero.

For comparison, the jump in the original variable is positive on its upper bank. If \(t>1/e\) and \(x(v)=1/t\), then

\[
 R(t+i0)=\frac{v\cot v+iv}{t},\qquad
 R(t+i0)-R(t-i0)=\frac{2iv}{t}.
\tag{6}
\]

Confusing the \(w\)-plane and \(z\)-plane upper banks would reverse this sign.

### Endpoint behavior

As \(v\downarrow0\),

\[
 \log x(v)=1-\frac{v^2}{2}-\frac{v^4}{36}+O(v^6),
 \qquad
 \rho(x)\sim\frac1\pi\sqrt{\frac{2(e-x)}e}\quad(x\uparrow e).
\tag{7}
\]

If \(\epsilon=\pi-v\downarrow0\), then

\[
 x(v)=\frac{e\epsilon}{\pi}e^{-\pi/\epsilon}(1+O(\epsilon)),
 \qquad \rho(x(v))=1-\frac\epsilon\pi.
\tag{8}
\]

Thus \(\rho\) extends continuously with \(\rho(0)=1,\rho(e)=0\), is strictly decreasing in \(x\), and lies strictly between zero and one on the interior. In particular \(1-\rho(x)\sim1/\log(1/x)\) as \(x\downarrow0\). The finite endpoint value
\(R(1/e)=e\) follows by monotone convergence; the integrand is locally proportional to \((e-x)^{-1/2}\), not nonintegrable.

## 3. Fractional moments: a self-contained Mellin proof

The bounded density shows \(M(s)\) is holomorphic in \(\Re s>-1\), by dominated differentiation on compact subsets. For real \(-1<s<0\), Tonelli's theorem and the beta integral give

\[
 \int_0^\infty t^{-s-1}R(-t)\,dt
 =\Gamma(-s)\Gamma(1+s)M(s).
\tag{9}
\]

But \(R(-t)=W_0(t)/t\). Set \(t=u e^u\), for \(u>0\). The same integral becomes

\[
 \int_0^\infty u^{-s-1}(1+u)e^{-(s+1)u}\,du
 =\Gamma(-s)(s+1)^s
   +\Gamma(1-s)(s+1)^{s-1}
 =\Gamma(-s)(s+1)^{s-1}.
\tag{10}
\]

Cancellation in (9) proves (2) on this real interval. Both sides are holomorphic in \(\Re s>-1\), so the identity theorem proves the full stated domain. No unjustified extension from integer moments is used.

Integration by parts in (1), with both boundary terms zero, also yields

\[
 M(s)=\frac1{\pi(s+1)}\int_0^\pi x(v)^{s+1}\,dv,
 \qquad\Re s>-1.
\tag{11}
\]

Consequently the Nuttall–Bouwkamp integral is recovered in its continuous form:

\[
 \frac1\pi\int_0^\pi
       \left(e^{v\cot v}\frac{\sin v}v\right)^p dv
 =\frac{p^p}{\Gamma(p+1)},\qquad\Re p>0.
\tag{12}
\]

For real \(p\downarrow0\), the left side tends to one. This limit of (12) does not mean the negative moment \(M(-1)\) is finite.

### The divergent lower Mellin endpoint

Set \(p=s+1>0\). The exact formula is

\[
 M(p-1)=\frac{p^{p-1}}{\Gamma(1+p)}.
\]

Using the convergent local log-gamma expansion,

\[
 \log M(p-1)=-\log p+p(\log p+\gamma)
              -\frac{\zeta(2)}2p^2+O(p^3),
\]

and therefore

\[
 M(p-1)=\frac1p+\log p+\gamma
 +\frac p2\big((\log p+\gamma)^2-\zeta(2)\big)
 +O\big(p^2(1+|\log p|^3)\big).
\tag{13}
\]

In particular \(M(s)\sim1/(s+1)\) as \(s\downarrow-1\), and \(M(-1)=\infty\). A bounded density is entirely compatible with this divergence: \(\rho(0)=1\) and the integrand at \(s=-1\) is asymptotically \(1/x\).

## 4. Strict positivity and quantitative series control

For any real \(\beta>-1\), nonzero vector \(c=(c_0,\ldots,c_{n-1})\), and polynomial \(P(x)=\sum c_jx^j\),

\[
 \sum_{i,j=0}^{n-1}\overline{c_i}c_jM(\beta+i+j)
 =\int_0^e x^\beta|P(x)|^2\rho(x)\,dx>0.
\tag{14}
\]

Thus every shifted Hankel matrix is positive definite and every principal Hankel determinant is strictly positive. Strictness follows from positive density on a whole interval, not just from nonnegative moments. Cauchy–Schwarz gives strict log-convexity

\[
 a_m^2<a_{m-1}a_{m+1}\quad(m\ge1).
\]

More generally \(\log M(s)\) is strictly convex for real \(s>-1\): its second derivative is the variance of \(\log x\) under the tilted probability measure \(x^s\mu/M(s)\).

For \(b_m=e^{-m}a_m\), the scaled Hausdorff inequalities are strict:

\[
 (-1)^k\Delta^k b_m
 =\int_0^e (x/e)^m(1-x/e)^k\,d\mu(x)>0
 \quad(m,k\ge0).
\tag{15}
\]

The coefficients have the exact ratio

\[
 r_m=\frac{a_{m+1}}{a_m}
     =\left(1+\frac1{m+1}\right)^m\uparrow e.
\tag{16}
\]

There is also a sharp-scale elementary upper bound without relying on a numerical Stirling approximation. Equation (3) gives \((\log x)'\le-v\), whence \(x(v)\le e^{1-v^2/2}\). By (11),

\[
 0<a_m<\frac{e^{m+1}}{\sqrt{2\pi}(m+1)^{3/2}}.
\tag{17}
\]

The local expansion (7), Gaussian rescaling in (11), and the same Gaussian majorant prove
\(a_m\sim e^{m+1}/[\sqrt{2\pi}(m+1)^{3/2}]\). The series is consequently absolutely convergent on the entire closed disk \(|z|\le1/e\), although the function has its branch point there.

Let \(P_N(z)=\sum_{m=0}^N a_mz^m\). A finite algebraic identity, valid off the cut, is

\[
 R(z)-P_N(z)=z^{N+1}\int_0^e\frac{x^{N+1}}{1-zx}\,d\mu(x).
\tag{18}
\]

Hence, for \(|z|<1/e\),

\[
 |R(z)-P_N(z)|\le
 \frac{a_{N+1}|z|^{N+1}}{1-e|z|}.
\tag{19}
\]

For real \(0<z<1/e\), Jensen's inequality and the endpoint chord bound, applied to the tilted measure \(x^{N+1}\mu/a_{N+1}\), strengthen this to

\[
 \frac{a_{N+1}z^{N+1}}{1-r_{N+1}z}
 \le R(z)-P_N(z)
 \le a_{N+1}z^{N+1}
       \frac{1-(e-r_{N+1})z}{1-ez}.
\tag{20}
\]

Both inequalities are strict when \(z>0\). At the branch point a finite uniform bound remains:

\[
 0<e-P_N(1/e)<\frac{2e}{\sqrt{2\pi(N+1)}}.
\tag{21}
\]

This follows from (17) and the integral comparison for \(\sum_{n\ge N+2}n^{-3/2}\). For arbitrary real \(t>0\), (18) also proves

\[
 \operatorname{sign}(R(-t)-P_N(-t))=(-1)^{N+1},\qquad
 \frac{a_{N+1}t^{N+1}}{1+r_{N+1}t}
 \le |R(-t)-P_N(-t)|\le a_{N+1}t^{N+1}.
\tag{22}
\]

This finite remainder statement remains valid even when the infinite Taylor series at \(-t\) does not converge.

## 5. Exactly rational Gauss/Padé brackets

Let \(G_N(z)=\sum_{j=1}^N\omega_j/(1-zx_j)\) be the \(N\)-node Gauss quadrature for \(\mu\). Its nodes lie strictly in \((0,e)\), its weights are positive, and it integrates polynomials through degree \(2N-1\) exactly. Thus \(G_N=P_{N-1}/Q_N\) is the \([N-1/N]\) Padé approximant of \(R\). All its coefficients are rational because every \(a_m\) is rational.

For \(0<z<1/e\), the \(2N\)-th derivative with respect to \(x\) of \((1-zx)^{-1}\) is strictly positive on \([0,e]\). Hermite interpolation at the Gauss nodes has a nonnegative error, zero only at the nodes, so

\[
 G_N(z)<R(z).
\tag{23}
\]

To produce a rational upper bound without inserting the irrational endpoint, choose a rational \(b>e\) with \(bz<1\). The measure
\(d\nu_b=(b-x)d\mu\) is positive and has rational moments
\(c_j=b a_j-a_{j+1}\). Define

\[
 I_b(z)=\int_0^e\frac{(b-x)d\mu(x)}{1-zx}.
\]

Let \(G_{b,N}\) be its \(N\)-node Gauss approximant. The same argument gives \(G_{b,N}<I_b\). Since

\[
 1-zI_b(z)=(1-bz)R(z),
\]

we obtain the certified rational bracket

\[
 \boxed{G_N(z)<R(z)<
 U_{b,N}(z):=\frac{1-zG_{b,N}(z)}{1-bz}}.
\tag{24}
\]

The restriction \(bz<1\) is essential for this inequality direction. It is not enough merely that \(z<1/e\). The proof is classical Gaussian quadrature, not a new theorem about Padé approximants.

For exact computation, write \(Q_N(z)=1+\sum_{j=1}^Nq_jz^j\) and solve

\[
 a_k+\sum_{j=1}^Nq_j a_{k-j}=0\quad(k=N,\ldots,2N-1).
\]

Then \(P_{N-1}\)'s coefficient of \(z^k\) is \(\sum_{j=0}^kq_j a_{k-j}\). Apply the identical construction to \(c_j\). No numerical root finding, approximate nodes, or floating-point linear algebra is needed. Positive Hankel matrices prove nonsingularity.

The exact test uses \(b=2719/1000\), and certifies \(e<b\) from

\[
 e<\sum_{j=0}^{10}\frac1{j!}+\frac1{10\cdot10!}
   =\frac{4697191}{1728000}<\frac{2719}{1000}.
\tag{25}
\]

Thus (24) applies to \(z=1/n\) for every integer \(n\ge3\), including the hardest endpoint of that family, \(z=1/3\). Exact rational endpoints, outward-rounded decimal displays, and a comparison to the raw series are saved in `verification.json`. The numerical value of Lambert W is used only as a diagnostic cross-check, not in the certification.


### 5.1 The actual mirror input is log(3)/3

The rational point \(1/3\) above is only a convenient exact benchmark. For the mirror family the actual input is \(z_n=\log n/n\), and \(z_3=\log3/3\) is the largest such input for integers \(n\ge3\), since \(\log x/x\) decreases for \(x>e\). It lies much nearer the branch point and is a substantially harder numerical test.

There is no need to trust floating-point logarithms. With \(K=100\), put

\[
 L_K=2\sum_{j=0}^{K-1}\frac{(1/2)^{2j+1}}{2j+1},\qquad
 E_K=\frac{2(1/2)^{2K+1}}{(2K+1)(1-1/4)}.
\]

The positive atanh series and its geometric tail prove the strict rational enclosure

\[
 L_K<\log3<L_K+E_K.
\tag{25a}
\]

Write \(z_-=L_K/3\), \(z_+=(L_K+E_K)/3\). Exact arithmetic verifies
\(b z_+<1\) (its value is approximately \(0.99570893762953\)). Since \(R\) is increasing on \((0,1/e)\), (24) gives

\[
 G_N(z_-)<R(\log3/3)<U_{b,N}(z_+).
\tag{25b}
\]

Every quantity on the outside is rational, and the displayed decimal endpoints in `mirror_endpoint_verification.json` are rounded outward by integer arithmetic. The comparison raw-series bracket uses the same available moments \(a_0,\ldots,a_{2N}\): its lower endpoint is \(P_{2N-1}(z_-)\), and its upper endpoint is
\(P_{2N-1}(z_+)+a_{2N}z_+^{2N}/(1-Ez_+)\), where \(E=4697191/1728000>e\). Thus the comparison has a precise equal-moment meaning. The proof, not a Lambert-W numerical evaluation, certifies the resulting transcendental-input bracket.

The completed 64-node certificate is

\[
 2.478052666984991255845343341327
 <R(\log3/3)<
 2.478052699570190519812004143718.
\tag{25c}
\]

Its bracket width is less than \(3.258520\times10^{-8}\). The equal-moment degree-127 raw-series bracket has width approximately \(0.090553184\), so this particular certified width is reduced by a factor exceeding 2.77 million. These are bounds on the actual mirror input, distinct from the rational benchmark: at \(z=1/3\), 32 nodes give width below \(4.591\times10^{-18}\), versus approximately \(3.9964\times10^{-5}\) for the degree-63 raw-series bound.

## 6. An exact branch-kernel bridge, and a stronger full-object bridge

The MRB branch kernel in the companion audit is

\[
 A(m)=\frac2\pi m^{-1/m}\sin(\pi/m).
\]

For \(0<p<1\), Euler reflection gives the exact identity

\[
 \boxed{A(1/p)=\frac{2p^p}{\Gamma(p)\Gamma(1-p)}
              =\frac{2p^2 M(p-1)}{\Gamma(1-p)}}.
\tag{26}
\]

This is a valid connection through the complete fractional-moment function, not an inference from matching sine factors. Equivalently, if
\(J(p)=\pi^{-1}\int_0^\pi x(v)^p dv\), then
\(A(1/p)=2pJ(p)/\Gamma(1-p)\). At positive integer \(p\), the reciprocal-gamma expression is interpreted by continuation and vanishes, as the sine form does. Its complex version uses a consistent branch for \(p^p\).

A stronger bridge reconstructs the entire MRB analytic kernel. For \(\Re w>0\), put \(p=1/w\), so \(\Re p>0\). Formula (2) is zero-free there and gives

\[
 \boxed{w^{1/w}=\frac{w^2}
 {\Gamma(1/w)M(1/w-1)}}.
\tag{27}
\]

Indeed, \(\Gamma(p)M(p-1)=p^{p-2}\), while
\(\Log(1/w)=-\Log w\) on the right half-plane. Thus (27) reconstructs the full analytic function used by Abel–Plana, not merely one sample or one residue.

Consequently, with the summation convention of the companion MRB audit,

\[
 \boxed{C_{\rm MRB}=
 \sum_{k=1}^\infty(-1)^k
 \left(\frac{k^2}{\Gamma(1/k)M(1/k-1)}-1\right)}.
\tag{28}
\]

The moment argument \(1/k-1\) approaches the divergent boundary \(-1\); it is never legitimate to replace these moments by a finite limiting value. More explicitly, (13) gives

\[
 M(1/k-1)=k-\log k+\gamma
 +O((\log k)^2/k),
\]

and the term in parentheses in (28) is exactly

\[
 k^{1/k}-1=\frac{\log k}{k}
                 +O((\log k)^2/k^2).
\tag{29}
\]

It decreases to zero for integers \(k\ge3\), since the derivative of \(x^{1/x}\) has the sign of \(1-\log x\). Thus (28) converges by the alternating-series test but not absolutely. Pairing consecutive even and odd terms is absolutely convergent; the pair difference is \(O(\log k/k^2)\).

### What this bridge does and does not accomplish

The complete measure \(\mu\), or equivalently its determinate full integer moment sequence together with the resulting measure, determines all the negative fractional moments in (28). It therefore determines the full MRB kernel and its alternating functional. This is a useful full-object map, stronger than identifying the single CM/honeycomb value \(A(8)\).

It is nevertheless a **nonlinear re-expression**, not an evaluation of the alternating functional. Neither (26) nor (28) proves a modular transformation, a cancellation estimate at the nome scale, or an explanation of the particular defect \(A(8)-C_{\rm MRB}\). Using the explicit gamma formula in (28) immediately returns the original root series. A new defect theorem would need a genuinely useful way to estimate or transform this transported functional, rather than merely substitute the exact identity. The entire-perturbation obstruction in the companion audit remains intact: branch jumps and finitely many samples do not determine a general analytic profile, whereas (27) deliberately supplies the full profile.


### 6.1 A hostile test: finite quadrature destroys the MRB endpoint

The usefulness of finite Gauss quadrature for the mirror does **not** permit blindly using it in the full MRB reconstruction. Suppose
\(\mu_N=\sum_{j=1}^N\omega_j\delta_{x_j}\) is any fixed finite positive quadrature measure with every \(x_j>0\). Then

\[
 M_N(p-1)=\sum_{j=1}^N\omega_jx_j^{p-1}
 \longrightarrow B_N:=\sum_{j=1}^N\frac{\omega_j}{x_j}\in(0,\infty)
 \quad(p\downarrow0).
\]

Since \(\Gamma(1/k)\sim k\), the reconstructed profile satisfies

\[
 \frac{k^2}{\Gamma(1/k)M_N(1/k-1)}\sim\frac{k}{B_N}.
\tag{30}
\]

It grows linearly instead of tending to one. Its putative MRB summands do not even tend to zero, so the alternating series cannot converge. The number of nodes here is fixed before \(k\to\infty\); no assertion about a carefully controlled joint limit is being made. Adding a positive atom at zero is not a fix: it makes all the required negative moments infinite. The one-node Gauss rule, an atom of weight one at \(x=1\), already exhibits the failure explicitly.

### 6.2 A continuous endpoint model preserves an absolutely summable correction

A reference measure can instead preserve the singular endpoint. On \((0,1)\) define

\[
 \rho_0(x)=1-\frac1{1+\log(1/x)},
\]

and set it to zero elsewhere. This is a finite positive **nonnormalized** measure; it is not asserted to be the tree probability measure. Put

\[
 H(p)=\int_0^1 x^{p-1}\rho_0(x)\,dx
     =\frac1p-e^pE_1(p),\quad p>0,
 \qquad E_1(p)=\int_p^\infty\frac{e^{-u}}u\,du.
\tag{31}
\]

The equality follows by setting \(t=-\log x\) and integrating
\(e^{-pt}(1-1/(1+t))\). In particular \(H(p)>0\). The elementary expansion of \(E_1\) yields

\[
 H(p)=\frac1p+\log p+\gamma
       +p(\log p+\gamma-1)+O(p^2|\log p|).
\tag{32}
\]

Comparison with (13) proves

\[
 M(p-1)-H(p)
  =\frac p2\log^2p+O(p|\log p|).
\tag{33}
\]

For the endpoint-model kernel
\(K_0(k)=k^2/[\Gamma(1/k)H(1/k)]\), both denominator moments are asymptotic to \(k\). Subtracting the two reciprocals in (27) and using (33) therefore gives the more precise result

\[
 \boxed{k^{1/k}-K_0(k)
  =-\frac{(\log k)^2}{2k^2}
     +O\!\left(\frac{\log k}{k^2}\right)}.
\tag{34}
\]

The correction is absolutely summable. Thus the original MRB alternating series can be decomposed into the model alternating series
\(\sum(-1)^k(K_0(k)-1)\) plus an absolutely convergent correction. The model series converges because its difference from the already-convergent original series is absolutely summable. This is an analytic, structure-preserving repair of the finite-quadrature failure, not a numerical-acceleration claim and not an explanation of the CM/nome defect. The tests report only a few kernel errors and their scaled limiting coefficient; no superiority to the session's existing MRB algorithms is asserted.

## 7. Verification and primary-source trail

Run `python verify_tree_moments.py` and `python verify_mirror_endpoint.py` from this directory. The first writes `verification.json`; the second certifies the actual mirror input and tests the endpoint model, writing `mirror_endpoint_verification.json`. Both print short summaries. Exact rational checks cover Hankel determinants of sizes 1–8 with shifts 0–4, coefficient-ratio identities through index 100, and Gaussian brackets with 1, 2, 4, 8, 12, 16, 24, and 32 nodes. Separate 100-decimal mpmath diagnostics check the density transform, both branch signs, fractional complex moments, the gamma identity, and the full analytic-kernel bridge. Quadrature agreement is not represented as interval certification; the proof plus rational computation supplies the actual enclosures.

Primary literature, checked on 2026-10-22:

1. G. A. Kalugin, D. J. Jeffrey, R. M. Corless, and P. B. Borwein, *Stieltjes and other integral representations for functions of Lambert W*, Integral Transforms and Special Functions 23 (2012), 581–593, DOI [10.1080/10652469.2011.613830](https://doi.org/10.1080/10652469.2011.613830). [Authors' preprint](https://www.uwo.ca/apmaths/faculty/jeffrey/pdfs/StieltjesLamW.pdf). Theorem 2.1 proves the Stieltjes property and obtains the boundary-density representation for \(W(z)/z\). The reciprocal-variable compact measure above is a direct form of this established result.
2. G. A. Kalugin, D. J. Jeffrey, and R. M. Corless, *Bernstein, Pick, Poisson and related integral expressions for Lambert W*, Integral Transforms and Special Functions 23 (2012), 817–829, DOI [10.1080/10652469.2011.640327](https://doi.org/10.1080/10652469.2011.640327). [Authors' preprint](https://www.uwo.ca/apmaths/faculty/jeffrey/pdfs/BersteinPick.pdf). Theorem 2.2 supplies the elementary \(v\)-integral for \(W(z)/z\); Corollary 2.4 gives the integer Nuttall integral and discusses its earlier real-parameter origin.
3. C. J. Bouwkamp, *A Conjectured Definite Integral (A. H. Nuttall)*, SIAM Review 28(4) (1986), 568–569, DOI [10.1137/1028165](https://doi.org/10.1137/1028165). [Primary publisher entry](https://epubs.siam.org/doi/10.1137/1028165). This is the original solution attached to Nuttall's integral problem, not a new integral discovered in this session. [Original published solution](https://pure.tue.nl/ws/files/3257614/Metis219248.pdf), checked directly. It uses a deformed Hankel integral; the Mellin proof in this note is a separate derivation.
4. The earlier combined preprint has only three authors and a different title: [arXiv:1103.5640](https://arxiv.org/abs/1103.5640), *Stieltjes, Poisson and other integral representations for functions of Lambert W* (2011). It should not be conflated bibliographically with the four-author published Stieltjes paper.
5. NIST DLMF [5.5.3](https://dlmf.nist.gov/5.5.E3) records Euler's reflection identity; [5.7](https://dlmf.nist.gov/5.7) records the reciprocal-gamma and log-gamma expansions used at the Mellin boundary. These standard identities supply the gamma step, not a claim about the MRB defect.

6. NIST DLMF [5.9.2_5](https://dlmf.nist.gov/5.9.E2_5) records the complex half-plane integral directly, with phase Φ(t)=1−t cot t+log(t/sin t), attributed to Temme (1996, equation (3.38), page 71). Since x(t)=exp(1−Φ(t)), it is exactly (12), including the domain Re p>0. This is a direct authoritative provenance check for the complete fractional formula, not merely its integer specialization.

7. NIST DLMF [6.6.2](https://dlmf.nist.gov/6.6.E2) gives the convergent E₁ expansion used in (32). The endpoint-model conclusions are derived here by substituting that standard expansion and comparing it with (13).
