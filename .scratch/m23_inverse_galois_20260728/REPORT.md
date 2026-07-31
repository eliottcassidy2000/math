# M23 inverse-Galois audit

Status: **scratch research report, not canon**  
Date: 2026-07-28

## Bottom line

The phrase “find a polynomial whose Galois group is \(M_{23}\)” is
field-sensitive.

* Over \(\mathbf Q\), the problem is still open: no polynomial in
  \(\mathbf Q[x]\) with Galois group \(M_{23}\) is presently known.
* Over the rational function field \(\mathbf F _2(X)\), Abhyankar's very small
  answer is
  \[
       f_X(Y)=Y^{23}+XY^3+1,
       \qquad
       \operatorname{Gal}(f_X/\mathbf F _2(X))\simeq M_{23}.
  \]
* In characteristic zero, Elkies gives an explicit degree-23 polynomial
  \(P(x)-t\) over a quartic number field \(F(t)\) whose geometric and arithmetic
  Galois groups are \(M_{23}\).  Its factored form is recorded below.

This separates the solved geometric/group-realization problem from the missing
arithmetic descent to \(\mathbf Q\).  A search over integer coefficients is not
currently a routine finite computation; it is the unresolved \(M_{23}/\mathbf
Q\) inverse-Galois case.

## 1. Abhyankar's characteristic-two polynomial

Let \(X\) be transcendental over \(\mathbf F _2\).  Abhyankar proved
\[
  \operatorname{Gal}\bigl(Y^{23}+XY^3+1\,/\,\mathbf F _2(X)\bigr)\cong M_{23}.
\]
The sign variant \(Y^{23}+XY^3-1\) is identical in characteristic two.

Two cheap exact checks explain the geometry.  The derivative is
\[
  \partial_Y f=Y^{22}+XY^2=Y^2(Y^{20}+X).
\]
A common zero of \(f\) and \(\partial_Yf\) would have \(Y\ne0\) and
\(X=Y^{20}\), whence \(f=Y^{23}+Y^{23}+1=1\), a contradiction.  Equivalently,
on the cover
\[
  X=Y^{20}+Y^{-3},\qquad dX/dY=Y^{-4}\ne0.
\]
Thus the degree-23 cover is étale over the affine \(X\)-line; its essential
wild ramification is at infinity.

Abhyankar's proof does substantially more than a cycle census.  A
\(2\)-linearized polynomial controls the roots, and the binary relation code is
the punctured binary Golay code.  Its automorphism group supplies the
\(M_{23}\) upper bound, while transitivity and ramification force equality.
This is the mechanism; finite-field specializations alone do not prove it.

At \(X=1\), for example, exact factorization gives degrees \((3,5,15)\):
\[
\begin{aligned}
Y^{23}+Y^3+1={}&(Y^3+Y^2+1)\\
&\cdot(Y^5+Y^4+Y^3+Y+1)\\
&\cdot(Y^{15}+Y^{10}+Y^9+Y^8+Y^7+Y^5+Y^3+Y+1).
\end{aligned}
\]
The hostile near-miss \(Y^{23}+Y+X\), often recalled incorrectly as
Abhyankar's polynomial, has at \(X=1\) factor degrees \((2,8,13)\), a cycle
shape unavailable in \(M_{23}\).

Primary provenance:

* S. S. Abhyankar, “Mathieu group coverings in characteristic two,”
  *C. R. Acad. Sci. Paris Sér. I Math.* **316** (1993), 267–271, Theorem 2.
* N. Katz, A. Rojas-León, and P. H. Tiep,
  “Rigid Local Systems and Sporadic Simple Groups,” *Memoirs AMS* **308**
  (2025), explicitly cite Abhyankar's Theorem 2 for
  \(Y^{23}+tY^3-1\).
* I. Yie, “Mathieu group coverings and Golay codes,”
  *J. Korean Math. Soc.* **39** (2002), 289–317, develops the
  linearization/Golay-code carrier.

## 2. Elkies's characteristic-zero polynomial

Set
\[
 F=\mathbf Q(g),\qquad g^4+g^3+9g^2-10g+8=0.
\]
This field contains
\[
 \eta=(2g^3+4g^2+16g-7)/3,\qquad \eta^2=-23.
\]
Define
\[
\begin{aligned}
P_2={}&(8g^3+16g^2-20g+20)x^2\\
 &-(7g^3+17g^2-7g+76)x-13g^3+25g^2-107g+596,
\\
P_3={}&8(31g^3+405g^2-459g+333)x^3\\
 &+(941g^3+1303g^2-1853g+1772)x
   +85g^3-385g^2+395g-220,
\\
P_4={}&32(4g^3-69g^2+74g-49)x^4\\
 &+32(21g^3+53g^2-68g+58)x^3\\
 &-8(97g^3+95g^2-145g+148)x^2\\
 &+8(41g^3-89g^2-g+140)x
   -123g^3+391g^2-93g+3228.
\end{aligned}
\]
Then
\[
 P(x)=P_2(x)^2P_3(x)P_4(x)^4,\qquad
 \operatorname{Gal}(P(x)-t/F(t))\cong M_{23}.
\]
Elkies proves the geometric monodromy is \(M_{23}\).  The arithmetic group is
also \(M_{23}\): geometric monodromy is normal in arithmetic monodromy, while
\(N_{S_{23}}(M_{23})=M_{23}\).

There are coprime polynomials \(P_7,P_8\) of degrees 7 and 8 with
\[
  P=P_7P_8^2+\tau,
  \quad
  \tau=\frac{2^{38}3^{17}}{23^3}
  (47323g^3-1084897g^2+7751g-711002).
\]
The three branch-cycle types are
\[
  1^3\,2^2\,4^4,\qquad 1^7\,2^8,\qquad 23.
\]
Their defects \(14+8+22=44=2\cdot23-2\) satisfy Riemann–Hurwitz.
Furthermore
\[
 \operatorname{disc}_x(P(x)-t)=
  (\text{square in }F^\times)\,t^{14}(t-\tau)^8.
\]
The passport and square discriminant put monodromy in \(A_{23}\), but leave
exactly the hard \(M_{23}\)-versus-\(A_{23}\) boundary.  Elkies separates them
with the action on 5-subsets: reduction modulo a degree-one prime over
\(10^8+7\), an exact factor count on the 5-set quotient, and the Weil bound
contradict the irreducibility/genus consequence that \(A_{23}\) would impose.
The 5-set resolvent is therefore the decisive sidecar; the discriminant is not.

Primary source:

* N. D. Elkies, “The complex polynomials \(P(x)\) with
  \(\operatorname{Gal}(P(x)-t)\cong M_{23}\),”
  *Open Book Series* **1** (2013), 359–367,
  [doi:10.2140/obs.2013.1.359](https://doi.org/10.2140/obs.2013.1.359),
  [publisher PDF](https://msp.org/obs/2013/1-1/obs-v1-n1-p18-s.pdf).

## 3. Exact local audit and its boundary

Run:

```bash
python3 .scratch/m23_inverse_galois_20260728/verify_elkies_passport.py
python3 -O .scratch/m23_inverse_galois_20260728/verify_elkies_passport.py
```

The script uses exact arithmetic in \(\mathbf Q(g)\) and independently checks:

* irreducibility of the quartic field polynomial and \(\eta^2=-23\);
* \(\deg P=23\);
* the exact identity \(P_2^2P_3P_4^4=P_7P_8^2+\tau\);
* squarefreeness and pairwise coprimality of all five factors;
* the passport, Riemann–Hurwitz defect sum, and square leading discriminant.

It deliberately reports that these checks leave \(M_{23}\) versus \(A_{23}\).
Replaying the final group separator requires Elkies's finite-field 5-set count,
or an independent certified monodromy computation (interval-certified branch
continuation followed by exact permutation-group recognition).

## 4. Current frontier over number fields

Häfner's braid-orbit analysis states that \(M_{23}/\mathbf Q\) remains open and
constructs regular realizations over \(\mathbf Q(\sqrt{-7})\) and
\(\mathbf Q(\sqrt{-15})\), using the Nielsen classes
\((14A,2A,2A,2A)\) and \((15A,2A,2A,2A)\):

* F. Häfner, “Braid orbits and the Mathieu group \(M_{23}\) as Galois group,”
  [arXiv:2202.08222](https://arxiv.org/abs/2202.08222).

A newer constructive-inverse-Galois paper again records \(M_{23}\) as the
remaining rational exception among the cases it discusses, while noting
realization over every number field \(K\) in which \(-1\) is a sum of two
squares:

* R. van Bommel, E. Costa, N. Elkies, N. Keller, M. Schiavone, and
  J. Voight, “The constructive inverse Galois problem via Hilbert modular
  forms: realizing 17T7,”
  [arXiv:2411.07857](https://arxiv.org/abs/2411.07857), current version
  dated 2026-06-03.

So the useful reframe is: the group theory and several geometric realizations
are known; the missing object is a rational descent/rational point on the
appropriate Hurwitz component.  A regular \(\mathbf Q(t)\) realization would,
by Hilbert irreducibility, yield number-field specializations, but is stronger
than merely finding one polynomial over \(\mathbf Q\).

## 5. Connections and non-connections to the repo

* **Torsor/JC:** the Galois closure of Abhyankar's affine cover is a finite
  étale \(M_{23}\)-torsor, and the 23-sheet cover is its quotient by a point
  stabilizer \(M_{22}\).  The entire obstruction is carried at infinity by
  wild ramification.  This is a sharp guardrail for Jacobian-conjecture
  analogies: affine étaleness in characteristic \(p\) does not control the
  boundary, and this example must not be transported to characteristic zero
  without its wild-inertia sidecar.
* **Resolvents/JC:** the square discriminant and branch passport forget the
  \(M_{23}\)-versus-\(A_{23}\) distinction.  Elkies's 5-subset quotient restores
  precisely the missing orbit information.  This is a concrete instance of
  the repo's rule that a discriminant or quotient invariant is not the full
  torsor carrier.
* **Codes rather than cosmetic tournaments:** Abhyankar's binary linearization
  produces a genuine relation object—the punctured Golay code—whose exact
  automorphism group is \(M_{23}\).  It is a stronger carrier than an imposed
  tournament on 23 roots because it records the preserved incidence relations.
* **GMC(2):** there is no established mathematical bridge to
  Mathieu–Zhao/Gaussian-moment work; “Mathieu” is only a shared name.  The
  transferable method is exact relation extraction (linearized roots
  \(\rightarrow\) binary code \(\rightarrow\) automorphism group), not a theorem
  linking the two subjects.
* **Speculative niche:** the 11-dimensional \(\mathbf F _2\) Todd/Golay module
  suggests asking whether modular-Galois constructions can be made to land
  exactly in \(M_{23}\), followed by a descent audit.  This is a research
  direction, not an existing construction.

## 6. The three supplied 2026 arXiv links

The links in the bundled prompt are not sources for the \(M_{23}\)
inverse-Galois problem:

* [arXiv:2607.22818](https://arxiv.org/abs/2607.22818) is on fully quantum
  walks;
* [arXiv:2607.14316](https://arxiv.org/abs/2607.14316) is on idempotent Schur
  multipliers;
* [arXiv:2607.24528](https://arxiv.org/abs/2607.24528) is on Feige's
  conjecture.

They may belong to other lanes in the user's multi-problem bundle, but they
provide no provenance for an \(M_{23}\) polynomial.
