---
id: THM-3200
title: "Fixed LRC-channel cleared-overlap quasipolynomial and mass recurrence boundary"
status: >
  PROVED + VERIFIED-EXACT.  For one fixed primitive cap-two channel (P,Q)
  and one ordered label lane (e,f) in the integral reflected-arc model of
  THM-3171, the unreduced cleared overlap numerator along dilation g is
  eventually a degree-at-most-two quasipolynomial of period dividing
  |Qe-Pf| (period one when Qe=Pf).  Its OGF is rational and it has a
  constant-coefficient tail recurrence.  The normalized mass is
  residue-wise quadratic-over-quadratic and P-recursive, but is C-finite
  exactly when it is eventually constant on every residue class.  This
  compiles each fixed dilation ray; it is not a uniform primitive-channel
  cone theorem or LRC(14).
audit: >
  The exact companion hash-pins THM-3171's digital-tent and affine-stability
  engine, reconstructs every certified residue polynomial on four positive
  and hostile primitive-3:5 lanes, checks all intervening physical heads,
  verifies 4,796 exact formulas and 4,724 constant-recurrence instances
  through dilation 1,200, and proves exact periods 3, 1, 19, and 1.
  Ordinary and optimized Python transcripts are byte-identical.
source: root/frontier-synthesis/lrc-channel-sequences/2026-08-02
depends_on:
  - THM-3171-global-high-channel-cell90-floor-and-all-width-uniform-two-star-law
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
script: 04-computation/lrc_fixed_channel_cleared_overlap_quasipolynomial_thm3200.py
output: 05-knowledge/results/lrc_fixed_channel_cleared_overlap_quasipolynomial_thm3200.out
script_sha256: 79c856f898a58cda0a0fcef21a72935848607ee9e5ec7397fc0163c56c5de51d
output_sha256: 4c4ab0c4c80e3eb0a538b6fee58bcec9c250e4983df09747839e2fa36e41b2d4
engine_dependency: 04-computation/lrc_small_channel_symbolic_tail_thm3171.py
engine_dependency_sha256: d73273a4cf4b88bea2890e001166d96cb07dd9b61f3a248ff1538ec44579796a
hash_basis: LF-normalized bytes
---

# THM-3200 -- fixed LRC-channel cleared-overlap quasipolynomial and mass recurrence boundary

**PROVED + VERIFIED-EXACT.**

## 1. Statement and normalization

Fix the integral reflected-arc model used at cell \(90\) in THM-3171:

\[
 L=168,\qquad
 A_e(p)=\bigcup_{k=0}^{p-1}
 \left[\frac{R+168k-12}{168p-e},
       \frac{R+168k+12}{168p-e}\right],                    \tag{1}
\]

where \(R=90e\bmod168\).  Fix coprime integers

\[
 P<Q\le2P                                                       \tag{2}
\]

and one ordered label lane \((e,f)\).  Along its common-dilation ray put

\[
 p=gP,\quad q=gQ,\quad
 z_g=168Pg-e,\quad w_g=168Qg-f,\quad D_g=z_gw_g.             \tag{3}
\]

Let \(I_g=\mu(A_e(gP)\cap A_f(gQ))\).  Its natural cleared integer
numerator is

\[
 N_g=D_gI_g.                                                  \tag{4}
\]

This is the unreduced digital-tent numerator, not the numerator obtained
after cancelling \(\gcd(N_g,D_g)\).  Define

\[
 C=Qe-Pf,\qquad
 M=\begin{cases}|C|,&C\ne0,\\1,&C=0.\end{cases}              \tag{5}
\]

Then the following assertions hold.

1. There is \(g_0\) and, for each residue \(r\bmod M\), a polynomial
   \(A_r(h)\in\mathbb Q[h]\) of degree at most two such that

   \[
   N_{r+Mh}=A_r(h)\qquad(h\ge g_0).                          \tag{6}
   \]

   Thus \(N_g\) is eventually quasipolynomial of degree at most two and
   period dividing \(M\).

2. The ordinary generating function of the full numerator sequence is
   rational.  After cancellation its denominator divides

   \[
   (1-x^M)^3.                                                \tag{7}
   \]

3. With \(E\) denoting the forward shift, the tail is annihilated by

   \[
   (E^M-1)^3.                                                \tag{8}
   \]

   Equivalently, for every sufficiently large \(g\),

   \[
   N_{g+3M}-3N_{g+2M}+3N_{g+M}-N_g=0.                       \tag{9}
   \]

   The recurrence order is at most \(3M\), and can be strictly smaller.

4. On each residue class the normalized overlap has the form

   \[
   I_{r+Mh}=\frac{A_r(h)}{D_r(h)},\qquad \deg D_r=2.         \tag{10}
   \]

   Hence every residue subsequence, and therefore their finite
   interlacing, is P-recursive; its OGF is D-finite.

5. The mass sequence \(I_g\) is eventually C-finite—equivalently, has a
   rational OGF—if and only if it is eventually constant on every residue
   class modulo \(M\).  Equivalently,

   \[
   A_r(h)=c_rD_r(h)                                         \tag{11}
   \]

   for constants \(c_r\) and every residue \(r\).

6. For every fixed rational target \(a/b\), the exact cleared comparison

   \[
   F_g=bN_g-aD_g                                             \tag{12}
   \]

   is another eventual degree-at-most-two quasipolynomial of period
   dividing \(M\).  An infinite fixed-channel floor assertion therefore
   reduces to finitely many exact quadratic minima and a finite head.

In roots-of-unity coordinates, assertion 1 is equivalently the eventual
closed form

\[
 N_g=\sum_{\zeta^M=1}
 (a_\zeta g^2+b_\zeta g+c_\zeta)\zeta^g.                    \tag{13}
\]

The finite residue state has been Fourier-diagonalized; its eigenvalues are
roots of unity and its Jordan blocks have size at most three.

## 2. Exact tent and the fixed cross-coordinate

Let

\[
 R=90e\bmod168,\qquad S=90f\bmod168,\qquad
 B_g=\frac{Rw_g-Sz_g}{168}.                                 \tag{14}
\]

The numerator in (14) is divisible by \(168\), since
\(z_g\equiv-e\), \(w_g\equiv-f\), \(R\equiv90e\), and
\(S\equiv90f\pmod {168}\).

THM-3171 derives the exact two-tent formula

\[
 N_g=\sum_{k=0}^{gP-1}\sum_{\ell}
 \left(
 [12(z_g+w_g)-168|B_g+kw_g-\ell z_g|]_+
 -
 [12(w_g-z_g)-168|B_g+kw_g-\ell z_g|]_+
 \right),                                                   \tag{15}
\]

where only an interval meeting \(A_e(gP)\) contributes.  The cap-two
condition leaves at most one contributing \(\ell\) for each \(k\), but
retaining the displayed finite candidate sum is convenient.

Write

\[
 k=r+tP,\qquad \ell=s+tQ,\qquad 0\le r<P.                   \tag{16}
\]

The load-bearing cancellation is

\[
\begin{aligned}
 B_g+kw_g-\ell z_g
 &=B_g+rw_g-sz_g+t(Pw_g-Qz_g)\\
 &=\alpha_{r,s}g+\beta_{r,s}+tC,                            \tag{17}\\
 Pw_g-Qz_g&=Qe-Pf=C.
\end{aligned}
\]

Only finitely many \(s\) occur, independently of \(g\).  Indeed,
contribution to (15) places \(\ell\) within a fixed distance of
\((B_g+kw_g)/z_g\).  After subtracting \(tQ\), its centre is

\[
 \frac{B_g+rw_g+tC}{z_g},                                   \tag{18}
\]

which is uniformly bounded for \(0\le t<g\), because its numerator and
denominator are both \(O(g)\).

Fix \(r,s\) and one of the two tents.  The constraints
\(0\le k<gP\), \(0\le\ell<gQ\), the tent support, and the sign split in
the absolute value make its contributing \(t\)-range an intersection of
intervals whose endpoints are either

\[
 ag+b\quad\hbox{or}\quad
 \left\lfloor\frac{ag+b}{C}\right\rfloor,\
 \left\lceil\frac{ag+b}{C}\right\rceil.                    \tag{19}
\]

The first type is integral affine.  The second type occurs only when
\(C\ne0\).

## 3. Quasipolynomial proof

Suppose first that \(C\ne0\), and restrict \(g\) to one residue
\(g=r_0+Mh\).  Every expression of the second type in (19) is then an
integral affine function of \(h\), because \(M=|C|\).  A maximum or minimum
of finitely many affine functions eventually selects one branch.  The
absolute-value sign in (17) stabilizes by the same comparison.

On a stabilized branch, each summand of (15) is affine in \(g,t\), and
the contributing \(t\)-interval has affine integral endpoints in \(h\).
Its sum uses only

\[
 \#\{t:L(h)\le t\le U(h)\},\qquad
 \sum_{t=L(h)}^{U(h)}t,                                     \tag{20}
\]

so it is a polynomial of degree at most two in \(h\).  There are finitely
many \(r,s\) and two tents.  Adding the outer terms and subtracting the
inner terms proves (6).

If \(C=0\), expression (17) is independent of \(t\).  After finitely many
affine branch crossings, an affine value times the affine count (20) is
quadratic.  This proves (6) with \(M=1\).

For a degree-at-most-two polynomial \(q(h)\),

\[
 \sum_{h\ge0}q(h)y^h
\]

has denominator dividing \((1-y)^3\).  Splitting the full sequence by
residues and putting \(y=x^M\) proves (7).  Three \(M\)-step finite
differences kill every quadratic residue polynomial, which proves
(8)--(9).  Discrete Fourier inversion on \(\mathbb Z/M\mathbb Z\) gives
(13).

In the high-channel scope of THM-3171, \(I_g\) has a positive floor, so
\(N_g\) has degree exactly two.  The theorem retains the degree-at-most-two
form because it also applies to fixed lanes outside that positive-floor
hypothesis.

## 4. The normalized-mass boundary

Equation (10) follows immediately from (3) and (6).  It also gives the
explicit first-order polynomial-coefficient recurrence on one residue:

\[
 A_r(h)D_r(h+1)I_{r+M(h+1)}
 -
 A_r(h+1)D_r(h)I_{r+Mh}=0.                                \tag{21}
\]

Thus each residue sequence is P-recursive, and finite interlacing preserves
P-recursiveness.

It remains to prove the C-finite criterion.  Let
\(R(h)=A_r(h)/D_r(h)\).  Suppose its values satisfy an eventual
constant-coefficient recurrence.  That recurrence holds at infinitely many
integers, so

\[
 \sum_{i=0}^d c_iR(h+i)=0                                   \tag{22}
\]

is an identity of rational functions.  Normalize so \(c_0\ne0\).  If \(R\)
had a finite pole, choose one with maximal real part.  At that pole the
\(i=0\) term in (22) cannot be cancelled by any \(i>0\) term, since such
cancellation would require a pole of \(R\) with larger real part.  Hence
\(R\) has no finite poles and is a polynomial.  An overlap mass is bounded,
so this polynomial is constant.  This proves necessity in (11).

Conversely, if (11) holds on every residue, then \(I_g\) is eventually
periodic.  Its OGF is rational and it is C-finite.  This proves assertion 5.

Multiplication or division by a fixed nonzero rational preserves the
C-finite numerator structure.  Division by the moving quadratic \(D_g\)
does not, unless the cancellation (11) occurs.  Reducing the fraction
\(N_g/D_g\) introduces the nonlinear sequence \(\gcd(N_g,D_g)\).  After
removing any common polynomial factor, the remaining value-gcd divides a
fixed resultant and is periodic after a finite residue refinement; the
reduced numerator therefore remains eventual quasipolynomial, but the
clean period-dividing-\(M\) bound need not survive.  This is why (4) is the
canonical normalization.

Finally, (12) is a fixed rational linear combination of two
degree-at-most-two quasipolynomials, because \(D_g\) is quadratic with
period one.  Assertion 6 follows.

## 5. Four exact positive and hostile rays

The exact companion reconstructs the certified affine branches and then
checks every physical head.  All four formulas below hold for \(g\ge2\).

### 5.1 Sharp label-\(\{1,6\}\) orientation

Take \((P,Q;e,f)=(3,5;6,1)\).  Here \(C=27\), while the exact period is
three:

\[
 N_g=4284g^2-2520g+C_{g\bmod3},\qquad
 (C_0,C_1,C_2)=(0,-336,84).                                \tag{23}
\]

The minimal annihilator is

\[
 (E-1)^3(E^2+E+1),                                         \tag{24}
\]

of order five rather than the universal order \(81\).  Moreover,

\[
 I_g=\frac{N_g}{(504g-6)(840g-1)}
 =\frac{17}{1680}-\frac{8213}{1411200g}+O(g^{-2}).          \tag{25}
\]

The nonzero inverse-linear term is an explicit witness that the mass is not
C-finite.

### 5.2 Reverse label-\(\{1,6\}\) orientation

For \((P,Q;e,f)=(3,5;1,6)\), \(C=-13\), but all residue modes cancel:

\[
 N_g=12096g^2-24g,                                         \tag{26}
\]
\[
 I_g=\frac{N_g}{(504g-1)(840g-6)}
 =\frac1{35}+\frac1{4900g}+O(g^{-2}).                       \tag{27}
\]

The numerator has the order-three annihilator \((E-1)^3\), while the mass
again fails the C-finite criterion.

### 5.3 Genuine label-12 residue tail

For \((P,Q;e,f)=(3,5;12,1)\), \(C=57\), and the exact period is \(19\).
Writing \(r=g\bmod19\),

\[
 N_g=\frac{149184g^2+5664g+c_r}{19},                        \tag{28}
\]
\[
\begin{aligned}
(c_0,\ldots,c_{18})={}&
(0,-7788,-900,7440,14268,19584,23388,25680,26460,25728,\\
&14592,7872,12180,14976,16260,16032,14292,11040,6276).
\end{aligned}                                               \tag{29}
\]

The period vector is rational, nonconstant, and \(19\) is prime, so every
primitive nineteenth-root mode is present.  The minimal annihilator is

\[
 (E-1)^3\Phi_{19}(E),                                       \tag{30}
\]

of order \(21\).  The mass has

\[
 I_g=\frac{N_g}{(504g-12)(840g-1)}
 =\frac{37}{1995}+\frac{103}{88200g}+O(g^{-2}),             \tag{31}
\]

and is not C-finite.

### 5.4 Reverse label-12 cancellation

For \((P,Q;e,f)=(3,5;1,12)\), \(C=-31\), but again the exact period is one:

\[
 N_g=7644g^2-72g,                                           \tag{32}
\]
\[
 I_g=\frac{N_g}{(504g-1)(840g-12)}
 =\frac{13}{720}+\frac{1571}{12700800g}+O(g^{-2}).          \tag{33}
\]

Thus the safe period bound is often far from minimal, and orientation can
retain a genuine residue tail or collapse it completely.

## 6. What the compiler replaces

For one fixed channel, ordered lane, and integral reflected cell, the proof
is constructive:

1. enumerate the finite \((r,s)\) tent bank in (16);
2. certify affine branch stability on each residue modulo \(M\);
3. interpolate its quadratic from three exact values and check a fourth;
4. verify the finite physical heads; and
5. minimize each comparison quadratic (12) on the integers.

This replaces every unbounded dilation scan on that fixed ray.  Once the
formula is compiled, evaluation is \(O(1)\) in \(g\), rather than a fresh
Euclidean floor-moment evaluation.

The argument also applies to any fixed reflected cell whose arcs have the
same integral endpoint form as (1): the cell changes
\(\alpha_{r,s},\beta_{r,s}\) in (17), but not \(C=Qe-Pf\).  A finite split
of boundary arcs preserves the conclusion.  If the cleared endpoint data
have an extra fixed denominator \(d\), a safe period bound is \(dM\).

This does not replace a scan over all primitive \((P,Q)\).  The period
\(|Qe-Pf|\), the finite state bank, and the stabilization height can all
grow.  A uniform cone theorem still needs a finite primitive-channel
reduction or an analytic estimate outside a finite bank.  No arbitrary
channel closure, physical-survivor classification, or proof of \(LRC(14)\)
is claimed.

## 7. Exact replay

The companion hash-pins

~~~text
04-computation/lrc_small_channel_symbolic_tail_thm3171.py
SHA256 d73273a4cf4b88bea2890e001166d96cb07dd9b61f3a248ff1538ec44579796a
~~~

and uses its exact digital-tent evaluator and future-stability certificate.
It proves the residue tails, checks every intervening physical head, and
then directly verifies:

~~~text
formula checks:     4,796 through g=1,200
recurrence checks:  4,724 through g=1,200
exact periods:      3, 1, 19, 1
~~~

Run

~~~bash
PYTHONDONTWRITEBYTECODE=1 \
python3 04-computation/lrc_fixed_channel_cleared_overlap_quasipolynomial_thm3200.py
PYTHONDONTWRITEBYTECODE=1 \
python3 -O 04-computation/lrc_fixed_channel_cleared_overlap_quasipolynomial_thm3200.py
~~~

and compare both byte-for-byte with the declared output.  No floating-point
comparison participates in the audit.

**End of proof.**
