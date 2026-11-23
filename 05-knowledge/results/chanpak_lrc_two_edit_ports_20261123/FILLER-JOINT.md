# One physical selector for phase and relative height

**PROVED** for the selection and capacity arguments, using the indicated
**FINITE-EXACT** certificate for this bank. This extends the separate
[max-gap selector](FILLER-CIRCLES.md) and [CRT wheels](FILLER-CRT.md).
It is not an exact absolute inverse or a new Lonely Runner inequality.

## Simultaneous same-phase straddle

Use bank types \((0,2,6,11,13)\) (zero-based) with 240 copies of each.
Write \(x_i\) for the number of second alternatives of type \(i\). All
five pairs have nine vertices and match in \(D,Q\) within their type.
The ordinal product feeds the [eight-site frame](THIRD-SLOT-JACOBIAN.md),
so the order is \(9\cdot1200+10=10810\) and \(D,Q\) are fixed across
all choices. Its physical corner heights are
\[
 h=11D+42Q+20R,\quad 71D+258Q+116R,\quad
 23D+90Q+44R,\quad 7h.
\]
The normalized mixed coefficient \(d=-3D-6Q\) is also fixed, hence the
gain-map Jacobian is fixed. This does not say the gain offsets are fixed.

Let \(C=11D+42Q\), \(x_*=(0,0,211,196,103)\), and let \(R_*\) be the
product of third readings at \(x_*\). For any unit \(s\bmod77\) and any
real target \(X\) satisfying
\[
  1.0065 \le (X-C)/(20R_*) \le 1.0070,
\]
there are explicit choices \(x^-,x^+\), both with \(R\equiv s\bmod77\),
whose actual starting heights satisfy
\[
 h^-\le X<h^+,\qquad h^+-h^- < (164/10^9)(X-C).
\]
In particular the nearer of the two is within 0.082 ppm of the target's
height above \(C\). There is no claim that the filler switch relating
the two backgrounds is a single arc. Each background separately
has the same two literal frame edits that lead to its own sevenfold height.

To reproduce, write \(P(v)=\prod_i r_i^{v_i}\), where
\[
 (r_i)=(292461/292221,153207/150807,299109/298869,284037/283797,
 317643/317403).
\]
Read the following rows as signed changes of the copy counts:
\[
 \begin{pmatrix}
 0&0&2&-1&0\\
 1&0&-1&1&-1\\
 3&0&-1&-10&9\\
 5&0&-5&-1&1\\
 8&1&-14&-6&-8
 \end{pmatrix}.
\]
Their rational products \(g_i=P(v_i)\) are all \(>1\) and \(1\bmod77\).
Their step sizes \(g_i-1\) are approximately \(760.382,107.727,9.654,1.803,
0.163931\) ppm. The exact certificate shows:

1. Among \(a\in\{0,1,2,3,4\}^5\) with \(1\le P(a)\le1.0065\), every unit
   \(R_*P(a)\bmod77\) occurs.
2. With \(b=(9,7,11,5,11)\), \(g_0^{b_0+1}>1.0070\) and
   \(g_i^{b_i+1}>g_{i-1}\) for \(i>0\); \(g_4<1+164/10^9\).
3. The whole box \(x_*+a+\sum_i q_i v_i\) fits in \([0,240]^5\), even
   allowing \(q_4\le12\) and the other \(q_i\le b_i\).

Take an anchor of the required residue. Divide \((X-C)/(20R_*P(a))\)
successively by the largest powers of the \(g_i\) that leave the
quotient at least one. The inequalities bound the digits by \(b\) and
the final quotient by \([1,g_4)\). The lower count is obtained using the
resulting digits, and the upper by adding one copy of the last signed
move. Their difference is \((g_4-1)(h^--C)\). Neither move changes the
chosen phase. This gives a finite selector in a common absolute band
for the fixed layout, not an asymptotic density statement about that layout.

The replay also verifies \(D\equiv385,\ Q\equiv1156\bmod2310\) and
\(R\equiv21\bmod30\) for every copy choice. Exactly 32 unit phases
give all three lower/intermediate counts coprime to \(2\cdot3\cdot5\cdot7\cdot11\):
exclude \(s\in\{2,4\}\bmod7\) and \(s\in\{1,10\}\bmod11\).
The upper count is a multiple of seven by design.

## Sharp capacity of the fixed bank, without quotas or guards

Let \(R_1,\ldots,R_{28}\) be the 28 listed readings, \(\lambda_j=\log R_j\).
Allow arbitrary multiplicities of every alternative in the full bank.
The number of distinct starts of the eight-site frame in \([X,2X]\) is
\[
 \left(\frac{\log2}{27!\prod_{j=1}^{28}\lambda_j}+o(1)\right)(\log X)^{27}.
\]
The starting count is independent of the ordering of factors, since the
three ordered slot readings multiply. This statement uses the original
bank, with no newly added filler, no arbitrary ambient padding that
changes the readings, and no coprimality conditions.

Every \(R_j>1\), and \(20R\le h\le73R\). A direct upper bound
\(O((\log X)^{27})\) follows by fixing 27 multiplicities, each \(O(\log X)\),
and noting that the last one lies in an interval of bounded length.
Consequently there are intervals of length \(\Omega(X/(\log X)^{27})\)
in \([X,2X]\) with no start from the bank, even with variable quotas.

For the sharp asymptotic, the certificate gives 27 individually private
primes across the 28 scalar readings, not merely across the 14 ratios.
Their valuations kill every exponent in a multiplicative relation but the
one for \(66825>1\), which must then vanish too. The \(R(e)\) are therefore
distinct. Count first the exponents for \(20R(e)\in[X,2X]\), a strip of width
\(w=\log2\) at \(T=\log(X/20)\) in \(\sum\lambda_j e_j\). Fix \(e_3,\ldots,e_{28}\)
and put \(U=T-\sum_{j\ge3}\lambda_j e_j\). For \(U\ge0\), count \(e_2\) for
each \(0\le e_1\le U/\lambda_1\) by the number of integers in an interval of
length \(w/\lambda_2\); the remaining admissible \(e_1\) contribute \(O(1)\).
Irrationality of \(\lambda_1/\lambda_2\) gives the count
\(wU/(\lambda_1\lambda_2)+o(U+1)\) as \(U\to\infty\), uniformly in starting
phase (irrational-rotation equidistribution for intervals, by approximating
with finite trigonometric polynomials). Summing over the remaining coordinates
reduces the leading term to the positive linear-weight simplex Riemann sum
\[
 \sum_{e_3,\ldots,e_{28}\ge0}(T-\textstyle\sum_{j\ge3}\lambda_je_j)_+
 \sim \frac{T^{27}}{27!\prod_{j\ge3}\lambda_j}.
\]
For any fixed residual cutoff the bounded-\(U\) slices, including
\(-w\le U<0\), contribute \(O(T^{26})\); for the others uniformity
makes the total error arbitrarily small relative to \(T^{27}\).

Finally, every bank entry has \(D_j\le Q_j\) and \(Q_j^7<R_j^6\) by the
certificate. By multiplication, \(h=20R+O(R^{6/7})\), giving
\(h/(20R)=1+O(X^{-1/7})\) when \(h\) or \(20R\) is in \([X,2X]\).
For all products with heights in either interval, coefficients are
\(O(\log X)\), and the nonzero difference of their log-products
is \(\gg(\log X)^{-B}\) for fixed \(B\), by
[Matveev's fixed-log lower bound](https://doi.org/10.1070/IM2000v064n06ABEH000314).
This dominates the power-decaying corruptions of \(\log(20R)\). Hence
the actual starts are eventually distinct and retain the \(R\)-ordering within
this pool, and at most a bounded number can cross either dyadic endpoint.
The elementary upper bound and hole statement do not require Matveev.

## What kind of expansion could matter?

Within the same multiplication-only grammar, suppose the menu for
\([X,2X]\) is allowed to grow to \(m(X)\) distinct triples, each from a
nonempty filler. Each satisfies \(R_j\ge3\): assign each vertex to one
of the three ordered slots and choose one Hamilton path per occupied
slot. Thus there are at most \(L=\lfloor\log(X/10)/\log3\rfloor\) factors
in that dyadic range and at most \(\binom{m(X)+L}{L}\) different count vectors.
If \(m(X)=o(\log X)\), the binomial bound is \(X^{o(1)}\), which by
prime density is far too small to cover all primes there. More precisely
any all-prime menu must satisfy \(\liminf m(X)/\log X\ge c_*\),
where \(A=1/\log3\) and the unique \(c_*>0\) is approximately \(0.58419\) and solves
\((A+c_*)\log(A+c_*)-A\log A-c_*\log c_*=1\), by Stirling.
This is a necessary alphabet-size budget, not a sufficiency result.
A grammar with additional order-sensitive or nonmultiplicative joins is
outside this obstruction.

**Transfer and loss.** The LRC maximum-gap mechanism works here after
passing to logarithmic height and choosing a modular phase. The signed
moves lie in the kernel of the finite phase map and still move in the real
direction. Neither coordinate counts the density of the physical integers;
the thin-strip asymptotic supplies that separate sidecar and explains why
arbitrarily fine fixed relative precision does not yield all absolute targets.

Replay: python3 filler-joint-check.py. The script uses exact rational and
integer arithmetic, verifies both sides in all 60 phases at four rational
targets, and checks the inequalities that apply on the whole band.
