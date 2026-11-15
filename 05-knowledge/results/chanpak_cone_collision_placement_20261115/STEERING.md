# Target-first complementary carries: a smaller absolute bracket

Use the [physical mask-zero carrier](../chanpak_local_rows_span3_20261115/DYADIC_LOCAL.md)
$$
F(x,y)=100xy+89x+11\big(p(x)+p(y)+1+(x\bmod2)\big),\qquad x,y\ge1 .
$$
Here $p(2u+1)=p(u)$, $p(2u)=p(u)+u+(u\bmod2)$, and $p(0)=0$. In particular $0\le p(v)\le v$.

**All-scale result.** For each sufficiently large prescribed $N$, a deterministic decoder finds literal counts strictly on either side of $N$, at distance $O(N^{1/3})$. It uses a logarithmic number of digit or quadratic evaluations, and also gives **every exact solution in its target-selected fixed-sum half-slice** in the same complexity. This improves the earlier [all-target $N^{3/8}$ bracket](../chanpak_absolute_steering_20261115/TWO_TAPE_STEERING.md) for a physical two-tape carrier. It does not assert that the chosen slice contains an exact solution. In fact the union of exact hits for this particular center-selection algorithm is sparse.

## 1. Let the low bits move, keep their carry fixed

Take $q=2^r\ge2$, an integer $u$, and set $s=2qu+q-1$, $a=s/2$. In the half-slice $x+y=s,\ x>y\ge1$, write
$$
 x=q(u+j)+v,\qquad y=q(u-j)+q-1-v,\qquad 0\le v<q .
$$
Every integer $x>a$ in the slice has this representation; there is no jump of size $q$ when $x$ advances by one. Put $d=x-a>0$, $R=\operatorname{rev}_r(v)$, $U=u+j,\ V=u-j$. The low words are complementary. The inherited split identity therefore gives the exact equality
$$
p(x)+p(y)=p(U)+p(V)+p(v)+p(q-1-v)+u(q-1)+j(q-1-2R)+(U\bmod2).
$$
Indeed the reversed low words sum to $q-1$, and the opposite high quotients have the same parity. The sign-sensitive piece $j(q-1-2R)$ has absolute value at most $|j|q\le d+q/2$. The other unfrozen terms are bounded by $2u+q$. Thus the exact count differs from the smooth quadratic
$$
G(d)=100a^2+89a+11+11u(q-1)-100d^2+89d
$$
by at most
$$
 |F(x,y)-G(d)|\le 11(d+2u+2q+2).
$$
The earlier all-target proof froze both low addresses, trading a small digital error for a physical step $q$. Complementary carries allow **unit movements in the ordinary addresses**, with the differential high quotient paying only the displayed $O(d)$ error.

## 2. Choose the center from the absolute target and solve both envelopes

Let $b=\lfloor\sqrt{(N-11)/100}\rfloor+1$ and take $q$ dyadic of size $b^{1/3}$. Choose the least $a=qu+(q-1)/2$ with $a\ge b+q$, so $a-b$ lies between $q$ and $2q$. The extra $q$ ensures a useful lower derivative bound, not just an upper bracket.

At the starting half-integer $d=1/2$ the quadratic exceeds $N$ by $\Theta(bq)$; its error there is $O(b/q+q)=o(bq)$. At the far end $y=1$ the upper envelope is $O(b)=o(N)$. Each of
$$
 G_{\pm}(d)=G(d)\pm11(d+2u+2q+2)
$$
is a concave quadratic whose upper level set on this half-slice is an initial interval. All exact candidates must have $G_-(d)\le N\le G_+(d)$. The first crossing of the lower envelope and the last crossing of the upper envelope are on scale $\Theta(\sqrt{bq})=\Theta(b^{2/3})$, where each slope has absolute value $\Theta(b^{2/3})$. The envelope gap there is $O(\sqrt{bq}+b/q+q)=O(b^{2/3})$. Hence the two crossings are only $O(1)$ ordinary addresses apart. Bisection finds their integer boundary positions; checking those $O(1)$ points with the **literal digital formula** gives the complete exact fiber in the slice. The adjacent exterior points give strict lower and upper physical witnesses, with error $O(b^{2/3})=O(N^{1/3})$. In the replay the envelopes are multiplied by four so every comparison is integral.

## 3. Why a decoder is not a supply proof

On a dyadic $b$-interval of scale $A$, the selected $q$ takes only constantly many values, and only $O(A/q)=O(A^{2/3})$ centers can be selected. At each center exact hits have $d=O(A^{2/3})$. Consequently the entire set of prescribed targets at most $X$ that hit in **their own selected slice** has cardinality $O(X^{2/3})$, summing the $O(A^{4/3})$ count geometrically. It cannot supply every large target. This is a barrier for the center-selection rule, not for the carrier with two unrestricted tapes.

The constructive next step is a host change that supplies the $O(N^{1/3})$ residual without disturbing the selected complementary carry. The existing fixed-core log-width slide cannot be substituted into this conclusion without another amplification theorem. The paired improvement and exact interval inversion are all-order results; the [script](verify_steering.py) checks only the implementation of the split and the decoder, not the conjecture.
