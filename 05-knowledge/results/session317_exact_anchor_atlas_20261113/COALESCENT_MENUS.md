# One exact strong count with arbitrarily many different continuation values

The two-slot invariant is too rigid to preserve on every step. The
opposite strategy is to hold the **output count** fixed while creating
many different two-slot responses for the next step. This can be
done at unbounded scale with literal tournaments and an unbounded
number of responses, not just the isolated coalescences in the
[readiness note](SLOT_INVARIANTS.md).

Let $J(T,U)=C_3[T,U,1]$ be the strong cyclic join
$T\to U\to x\to T$, for nonempty $T,U$. Write $P=Q_T,\ R=Q_U$.
Its exact ordered-slot polynomial is
$$
 Q_{J(T,U)}(z)=z\left[
 P(z)R(z)+
 \sum_{\epsilon=\pm1}
 (P(z+\epsilon)-P(z))(R(z+\epsilon)-R(z))
 \right].
$$
Indeed every mixed cycle contains $x$ and one nonempty directed path
from each block. Its total length is odd precisely when the two path
orders sum to an even number. The [packing expansion](https://arxiv.org/abs/2307.05569)
weights each selected odd cycle by $2z$; the two slot shifts sum
the even-parity choices by $Q_{T[S]}(-1)=(-1)^{|S|}H(T[S])$
and the slot convolution. Packings without a mixed cycle give
$zP(z)R(z)$, proving the formula. At $z=1$ it reads
$$
 H(J(T,U))=(Q_T(2)-H(T))(Q_U(2)-H(U))+2H(T)H(U).
$$

For $r\ge1$ and $0\le i\le r$ take $T_i$ to be the ordinal
sum of $i$ directed triangles and $2(r-i)$ singleton vertices.
Thus
$$
 H(T_i)=3^i,\quad Q_{T_i}(2)=3^i4^r,\quad
 Q_{T_i}(3)=3^i11^i9^{r-i}.
$$
The ratios at the two-slot response coincide, but those at three
slots vary with the split. For $L\ge r$ use the strong tournaments
$$
 Y_i=J(T_i,T_{r-i}\to TT_L),\qquad 0\le i\le r.
$$
They all have exactly the same count
$$
 H_{r,L}=3^r[(4^r-1)(4^r2^L-1)+2].
$$
Their next two-slot responses are **all different**. At $z=2$
in the wrap formula, every term is constant in $i$ except
$$
 -2\cdot3^r4^r
 \left(2^L11^i9^{r-i}+3^L11^{r-i}9^i\right).
$$
The bracket strictly decreases with $i$ when $L\ge r$. For
consecutive $i$, the positive contribution to its difference is
proportional to $2^L11^i9^{r-i-1}$ and the negative to
$3^L11^{r-i-1}9^i$. Their ratio is
$(2/3)^L(11/9)^{2i-r+1}<1$ for $0\le i<r$.
This proves strict distinction and hence nonisomorphism, without
searching classes. The $r=1$ instance gives equal count
$36\cdot2^L-3$ and a two-slot difference of
$48(3^L-2^L)$.

Thus even after imposing strongness, a target can support an
unboundedly large menu of exact continuation values. The mechanism
uses **equal $Q(2)/H$, different $Q(3)/H$** and complementary
source counts. This is not the fixed-full-polynomial neutral
switch of Session 316: here the whole point is that the next
readout differs. The anchors $H_{r,L}$ remain sparse, so it does
not provide menus at arbitrary target values.

There is an important prime-target guard. The polynomial $t^3+2t$
is divisible by 3 at every integer $t$, so $Q_{T_i}(t)$ has
integer-value divisibility $3^i$. The two-block wrap formula
preserves products of such divisors, giving $3^r\mid Q_{Y_i}(t)$
at **every** integer $t$. Further cyclic wraps or ordinal singleton
padding preserve that common divisor. Thus the menu does not
give a prime-target grammar by homogeneous repetition; to use it
there one would have to change the evaluation divisor through
selective contacts, as opposed to merely changing the scalar
$Q(2)$. The [contact capacity theorem](CONTACT_CAPACITY.md) is a
separate guard on how much such selective contacts can destroy.
