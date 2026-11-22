# The exact rank-loss set and the last-gap congruence selector

For the [shared-gap chain](SHARED-CONTINUANTS.md), a positive boundary
has state \((D,0,G)\), with \(D\) positive even, \(G\) positive odd and
\(D\ge2G\). The empty-gap map is \(Z(D,U,V)=(D+2U,V,D+U)\).
After \(j\) empty gaps, write its height and middle coordinate as
\[
 F_*=a_jD+b_jG,\qquad U_*=c_jD+d_jG.
\]
A next nonempty gap of profile \((h,Q)\) gives
\((D',0,G')=(Q(F_*+U_*),0,hF_*)\).

**Exactly two singular run lengths.** The determinant
\(\Delta_j=a_jd_j-b_jc_j\) obeys
\(\Delta_{j+3}=-\Delta_{j+2}-\Delta_{j+1}+\Delta_j\)
with initial values \(0,2,-3\), by the exterior-square transfer.
Let \(T_0=0,T_1=T_2=1\) and extend the tribonacci recurrence backwards.
The same initial values and recurrence identify
\(\Delta_j=T_{-j-4}\). The integer zero set of this tribonacci sequence
is \(\{0,-1,-4,-17\}\), by [Mignotte and Tzanakis, Example 2](https://doi.org/10.4064/aa-57-4-357-364),
with this indexing recorded by [Bilu et al.](https://arxiv.org/pdf/2210.16959).
Hence the only \(j\ge0\) with rank loss are **0 and 13**. Applying the
next gap does not change rank over the rationals because its map
\((F_*,U_*)\mapsto(Q(F_*+U_*),hF_*)\) has determinant \(-hQ\ne0\).

**Full rank does not imply preservation of content.** For a primitive
input pair and a prime \(\ell\nmid hQ\), put
\(w=Z^{-j}(1,0,-1)\) and write its coordinates as \(w_D,w_U,w_V\).
Then
\[
v_\ell(\gcd(D',G'))=
\min\{v_\ell(w_U),v_\ell(Dw_V-Gw_D)\},
\]
with \(v_\ell(0)=+\infty\).
The projection to \((F_*,U_*)\) has primitive kernel vector
\((1,0,-1)\) and is split over the integers. Since \(\det Z=1\),
vanishing of both outputs modulo \(\ell^a\) means that
\((D,0,G)\) is a multiple of \(w\) modulo \(\ell^a\); the multiplier
is a unit for primitive input. Equivalently the two displayed
expressions vanish modulo \(\ell^a\). For the converse, when \(w_U\)
vanishes one of \(w_D,w_V\) must be a unit, and the cross equation
determines the common multiplier. The cross product of the two
projection rows gives \(w_U=\Delta_j\). Thus a nonsingular run can
create content only at primes dividing the determinant or the
profile factors. For example \((2,0,1)\) after six empty gaps becomes
\((48,20,37)\); the next transitive width-one gap gives \((250,0,105)\),
gaining a factor 5 at full rank.

## One congruence replaces the blind search for the second register

For an absolute target \(N\) and a fixed last template \((j,h,Q)\), put
\[
 P=h a_j+Q(a_j+c_j),\qquad R=h b_j+Q(b_j+d_j).
\]
A possible previous positive boundary must satisfy \(N=PD+RG\).
In particular \(R\) is odd because \(h,b_j\) are odd and \(Q\) is even
(the orbit \(Z^je_3\) retains odd total height). With \(d=\gcd(P,R)\),
reject if \(d\nmid N\). Otherwise every candidate has
\[
 M=2P/d,\qquad G\equiv (R/d)^{-1}(N/d)\pmod M,\qquad
 1\le G\le\left\lfloor N/(2P+R)\right\rfloor.
\]
Indeed evenness of \(D\) gives the congruence modulo \(2P\), and
\(D\ge2G\) gives the bound. Conversely each member of this progression
in the interval gives an exact integral pair in the positive cone
by \(D=(N-RG)/P\); its earlier physical compatibility still needs
the backward test. If there is no preceding positive gap, test the
initial scalar \(h+Q\) separately after normalization. The list has at most one element if
\(2P(2P+R)>dN\). For a prime \(N\), necessarily \(d=1\), since
\(N\ge 2P+R>d\) for any candidate. This includes composite as well as
prime targets and identifies when a large last template is worth
checking with no free second-register loop.
