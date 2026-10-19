# Polynomial-bit construction of a signed Pythagorean contact

**Status: PROVED**, using the previously proved logarithmic signed-index bound
and an elementary division-polynomial height estimate. No canonical-height
comparison, minimal-model computation, factorization, numerical elliptic
integral, or novelty claim is used. The companion controls are finite-exact;
they do not replace the uniform proof.

## 1. The theorem and its precise scope

Let `T=(a,b,C)` be a primitive positive Pythagorean triple with
`C²>4ab`, and choose either port `epsilon in {+1,-1}`. Put

\[
 L=1+\lceil\log_2 C\rceil,\qquad f_0=\frac{2ab}{C^2},\qquad
 k=\frac{|a+\epsilon b|}{C},\quad K=k^2.
\]

**Theorem.** There is an absolute constant `c` and a deterministic exact
algorithm that returns a distinct admissible primitive triple
`U=(a',b',C')` and a rational number `r>0` satisfying

\[
 r^2=1+\epsilon f_0-\epsilon\frac{2a'b'}{C'^2},
 \qquad \log C'\le cL^3.                                  \tag{1}
\]

The algorithm uses `O(L)` elliptic group additions, `O(L³)`-bit integer
operands, and polynomial bit time. A deliberately conservative bound with
schoolbook arithmetic and an unamortized Euclidean-algorithm analysis is
`O(L^10)` bit operations and `O(L³)` working bits, excluding the read-only
input and returned output. The output and its rational contact certificate
also have `O(L³)` bits.

This is a one-neighbor, unrestricted-spread-direction theorem for **either
prescribed port**. It does not give polynomial-time downward-spread search,
an arbitrary narrow prescribed window, avoidance of an arbitrary supplied
forbidden set, or a polynomial bound on the numerical value of `C'`.
Equivalently, the numerical bound is `C' <= exp(O((log C)^3))`.
A constrained window of width `Delta` may introduce a `1/Delta` search
factor, exponential in its binary encoding length; the theorem therefore
does not assert polynomial complexity for constructing arbitrary long
distinct-petal fans.

The input-length estimate is genuine: the three input integers have total
binary length `Theta(L)`. All constants below are independent of the anchor,
its port, and proximity to either admissibility boundary.

## 2. Inherited exact termination certificate

Use the controller and its anchor point from
[universal signed contacts](square_sum_ppt_bridges_20261019.md):

\[
 E_k:Y^2=(X-K)(X^2-1),\qquad P=(0,k),\qquad
 t_0=\frac{C+a}{b}>1,\quad (u_0,z_0)=(\epsilon t_0,t_0^2+1).
\]

The image of `(u_0,z_0)` is `Q`. The predecessor proves `P` non-torsion.
Let `D=2+f_0` for the plus port and `D=2` for the minus port. Let `J` be the
least nonnegative integer such that `2^J f_0 >= D`, and set

\[
 N=\left\lfloor\frac{100(4+J)}7\right\rfloor.              \tag{2}
\]

The [rational recurrence budget](square_sum_contact_recurrence_bound_20261019.md),
equation (18), guarantees success among

\[
 Q+2iP,\quad Q-2iP,\qquad 1\le i\le N,                   \tag{3}
\]

using the ordinary inverse chart with `t=epsilon*u>1`, `z>0`,
`0<f(t)<1/2`, and `f(t)!=f_0`. Indeed its certified successful open arc is
contained in that chart. Since `f_0>=2/C²` and `D<=5/2`,
`J<=2 ceil(log_2 C)+1`, so `N=O(L)`. Computing `J,N` needs only integer
shifts, products, comparisons, and division by 7.

An index count alone does not establish bit complexity. The remainder
supplies the missing uniform height bound.

## 3. An integral short model with controlled weighted size

Write `A=|a+epsilon*b|`, so `0<A<2C` and `A!=C`. First set

\[
 x=C^2X,\qquad y=C^3Y.
\]

Then

\[
 y^2=x^3+a_2x^2+a_4x+a_6,
 \quad a_2=-A^2,\quad a_4=-C^4,\quad a_6=A^2C^4,
 \quad \widetilde P=(0,AC^2).
\]

Now set

\[
 \xi=9x+3a_2,\qquad \eta=27y,
\]

obtaining the integral short model

\[
 \eta^2=\xi^3+\alpha\xi+\beta,
 \quad \alpha=81a_4-27a_2^2,
 \quad \beta=54a_2^3-243a_2a_4+729a_6.                    \tag{4}
\]

Its image `P_s` has coordinates `xi_0=-3A²`, `eta_0=27AC²`.
For `R=6C` the explicit inequalities

\[
 |\xi_0|\le R^2,\qquad |2\eta_0|\le R^3,\qquad
 |\alpha|\le R^4,\qquad |\beta|\le R^6                 \tag{5}
\]

hold: bounds `12C²,108C³,513C⁴,5400C⁶`, respectively, suffice.
These are rational isomorphisms, not isogenies or minimal-model changes.
The original cubic has distinct roots `K,1,-1`; hence (4) is nonsingular.
The image `Q_s` has rational coordinates of logarithmic height `O(L)`,
because the input, the displayed fixed-degree forward map, and both changes
of coordinates all have height `O(L)`.

## 4. Elementary uniform division-polynomial bound

For a rational `v=p/q` in lowest terms with `q>0`, write
`h(v)=log max(|p|,q)` and `h(0)=0`. All logarithms in this section are natural.

Use the classical short-Weierstrass division polynomials `psi_n` and

\[
 \xi([n]S)=\xi(S)-
 \frac{\psi_{n-1}(S)\psi_{n+1}(S)}{\psi_n(S)^2}.        \tag{6}
\]

The definitions, parity normalization, and multiplication formula are in
Andrew Sutherland's [MIT 18.783 lecture 5, sections 5.5--5.6, Lemma 5.20 and
Theorem 5.21](https://ocw.mit.edu/courses/18-783-elliptic-curves-fall-2025/mit18_783_f25_lec05.pdf).
They apply here over `Q`, of characteristic zero, on the nonsingular short
model (4). In particular `psi_n(P_s)!=0` for `n>=1` because `P_s` has infinite
order. The coefficient-norm estimate below is proved here rather than cited
as a height-comparison theorem.

Put `G=xi³+alpha*xi+beta` and normalize

\[
 F_n=\begin{cases}\psi_n,&n\text{ odd},\\
                  \psi_n/(2\eta),&n\text{ even}.
       \end{cases}
\]

Then `F_n in Z[xi,alpha,beta]`. Its weighted degree, for weights `(2,4,6)`, is
`n²-1` for odd `n`, and `n²-4` for even `n`. This follows directly by induction
from the following polynomial recurrences and initial values:

\[
 F_1=F_2=1,\quad
 F_3=3\xi^4+6\alpha\xi^2+12\beta\xi-\alpha^2,
\]
\[
 F_4=2(\xi^6+5\alpha\xi^4+20\beta\xi^3
          -5\alpha^2\xi^2-4\alpha\beta\xi-\alpha^3-8\beta^2),
\]
\[
 F_{2m}=F_m(F_{m+2}F_{m-1}^2-F_{m-2}F_{m+1}^2),\quad m\ge3,
                                                                    \tag{7}
\]
\[
 F_{2m+1}=\begin{cases}
  16G^2F_{m+2}F_m^3-F_{m-1}F_{m+1}^3,&m\text{ even},\\
  F_{m+2}F_m^3-16G^2F_{m-1}F_{m+1}^3,&m\text{ odd},
 \end{cases}\quad m\ge2.                              \tag{8}
\]

Let `||F||_1` be the sum of absolute values of its integer coefficients.
We claim the absolute bound

\[
 \log\max(1,\|F_n\|_1)\le4n^2\log(n+1).               \tag{9}
\]

Here are details of the elementary induction. The initial norms are
`1,1,22,88`; (7)--(8) give the claim also for `5,6,7`. For `n>=8`, every
child index is at most `n/2+2`; the sum of its squared indices, counted with
multiplicity in a product, is at most `n²+6`. Since `||16G²||_1=144`, taking
logs in either recurrence costs at most `log 288` in addition to the sum
of the child log norms. Also

\[
 n^2\log(n+1)-(n^2+6)\log(n/2+3)\ge20/9\quad(n\ge8).   \tag{10}
\]

For completeness, `(n/2+3)/(n+1)<=7/9`, `log(9/7)>=2/9`, and
`log(n/2+3)<=n/4` for `n>=8`. Consequently the left side of (10) is at least
`2n²/9-3n/2>=20/9`. Multiplying by 4 absorbs `log 288<8<80/9`, proving (9).
This bound concerns actual reduced universal polynomials; no division by a
polynomial or by a discriminant is hidden in the norm argument.

Using (5), weighted homogeneity, and the extra factor `2eta_0` for even `n`,
we get the explicit evaluation bound

\[
 \boxed{\quad
 \max(1,|\psi_n(P_s)|)
 \le(6C)^{n^2-1}(n+1)^{4n^2}\qquad(n\ge1).
 \quad}                                                \tag{11}
\]

The evaluations are integers. Substitution in (6) therefore gives

\[
 h(\xi([n]P_s))=O\bigl(n^2[\log C+\log(n+1)]\bigr).     \tag{12}
\]

For `n=1` use (5) directly; for `n>=2`, bound the integer numerator
`xi_0 psi_n²-psi_{n-1}psi_{n+1}` and denominator `psi_n²` before reduction.
There is no small-denominator estimate to prove. Reducing a fraction can
only decrease its height.

From (4),

\[
 2h(\eta([n]P_s))
 \le3h(\xi([n]P_s))+\log(1+|\alpha|+|\beta|),            \tag{13}
\]

so the same order bounds the second coordinate. Negative `n` only changes
the sign of `eta`.

## 5. Translation, inverse map, and primitive output

A fixed-degree rational expression with coefficients of height `O(L)`
in arguments of height at most `H` has height `O(H+L)` wherever its
denominator is nonzero. This elementary fact follows by clearing the
argument denominators and bounding the finitely many numerator terms.

Apply this fact to the chord-and-tangent addition formulas with `Q_s`.
The exceptional addition cases are harmless: a point at infinity is handled
separately; equal x-coordinates mean equal or opposite points; opposite
points add to infinity, and equal points use the doubling formula, with a
zero y-coordinate again giving infinity. Exact equality decides every case.
Thus, for every nonzero integer `q` for which the translate is finite,

\[
 h(X(Q+qP))+h(Y(Q+qP))
 =O\bigl(q^2[L+\log(|q|+1)]\bigr).                     \tag{14}
\]

The changes of model only add `O(L)` to this estimate. The same estimate
holds for the coordinates of `[q]P`.

For the ordinary quartic chart, use

\[
 u=\frac{kY+X-K}{X(X-K)},\qquad
 z=\frac{2u^2X-K-2u-Ku^2}{k}.                           \tag{15}
\]

At a finite point with `X=0` or `X=K`, skip the candidate; skip infinity as
well. The successful open arcs used to prove (2)--(3) contain none of these
points, so this does not lose the termination guarantee. All remaining
formulas and tests are fixed-degree exact rational operations. In particular,
`k>0` and all possible zero denominators are tested before division.

For `|q|<=2N=O(L)`, equations (14)--(15) give `h(u),h(z)=O(L³)`.
After the exact acceptance tests in section 2, reduce `t=epsilon*u=n/m>1`.
Then `n,m` are positive, coprime, and have `O(L³)` bits. Put

\[
 d=\begin{cases}2,&n,m\text{ both odd},\\1,&\text{otherwise},\end{cases}
 \qquad (a',b',C')=\frac1d(n^2-m^2,2nm,n^2+m^2),          \tag{16}
\]

and reorder the two legs if desired. This is a primitive positive triple,
with `log C'<=2h(t)+log 2=O(L³)`. Its spread differs from the anchor's, so
it represents a distinct petal even after swapping legs. The positive
rational certificate `r=z/(u²+1)` has `O(L³)` bits and satisfies (1).

## 6. Exact algorithm and bit complexity

1. Validate positivity, the Pythagorean equation, coprimality, and `C²>4ab`.
   Construct `P,Q`, and the rational stopping budget `N`.
2. Compute `B=2P` and maintain `R=iB`, starting with `R=O`.
3. For `i=1,...,N`, replace `R` by `R+B`. Test the two points `Q+R` and
   `Q-R`, handling every group-law and inverse-chart exception as above.
4. Return the first accepted primitive triple and rational root. Recheck
   its primitive, admissibility, changed-spread, and contact identities.

The recurrence theorem proves that step 4 occurs by `i=N`. This uses at most
`1+3N=O(L)` group additions, and a constant number of rational operations
and comparisons per candidate. It needs no real-number arithmetic.

Reduce rational numbers after every arithmetic operation. Equations
(12)--(15) bound every stored coordinate by `O(L³)` bits; each intermediate
within a fixed group-law or inverse-map expression is a fixed product or
sum of such integers, so it also has `O(L³)` bits. This is the essential
control preventing exponential rational-expression growth.

For an explicit uncomplicated upper bound, a long multiplication or division
of `B`-bit integers uses `O(B²)` bit operations. Euclid takes `O(B)` divisions
(the larger argument shrinks by at least a factor two every two steps), so
a deliberately unamortized gcd bound is `O(B³)`. Therefore one reduced
rational operation costs `O(B³)`, and the full search costs
`O(L*(L³)^3)=O(L^10)` bits. Standard improved gcd/arithmetic implementations
sharpen this substantially, but no improvement is needed for the theorem.
The search retains only a constant number of points and rational values,
so its working storage is `O(L³)` bits.

All estimates remain uniform as the real normalized discriminant tends to
zero. The proof bounds explicit integer polynomial evaluations and never
inverts the discriminant. It also avoids the incorrect inference that a
small number of group additions alone ensures small intermediate integers.

## 7. Finite-exact controls

The companion `square_sum_contact_output_height_20261019.py` checks:

- the short-model scaling identity symbolically;
- weighted degrees and the integer form of (9) for initial division polynomials;
- division-polynomial multiplication against independent exact group additions;
- the explicit integer bound (11), on both ports and near both parameter cusps;
- exact searches, primitive normalization, and rational contact certificates;
- all special group-addition and inverse-chart cases.

The normal and `python3 -O` runs are byte-identical: **2,373 exact checks**
(13 symbolic, 556 integer, 1,804 rational) and 16 returned witnesses. Both
normalization divisors are checked among the input anchors; all 16 selected
outputs happen to have divisor 1. The largest tested input hypotenuse has
34 bits; its plus-port output has 193,450 bits at `q=-54`, within the exact
budget 600. These are finite observations, not estimates of an optimal
constant or a worst-case lower bound. Every check uses an explicit exception,
not a removable Python assertion.
