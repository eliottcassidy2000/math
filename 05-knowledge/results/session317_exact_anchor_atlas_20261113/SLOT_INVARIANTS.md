# Readiness congruences, cyclic return barriers, and a changed-interface bypass

## Universal divisibility, before choosing a grammar

By the odd-cycle packing expansion $Q_T(t)=t^\epsilon P_T(t^2)$
with $\epsilon$ the order parity. For two tournaments of the same
parity, if $Q(t_i)$ agrees at distinct integers $t_i\ge2$, then
$P_+-P_-$ is divisible in $\mathbb Z[z]$ by
$\prod_i(z-t_i^2)$. Evaluating at 1 gives
$$
 \prod_i(t_i^2-1)\mid H_+-H_-.
$$
The difference of Hamilton counts is even as well. For different
parities, use $Q_+^2-Q_-^2$ to get the corresponding divisibility
of $H_+^2-H_-^2$. In particular $Q(2)\equiv(-1)^nH\bmod3$.
A $+2$ step preserving $Q(2)$ must change parity and must go from
$H\equiv2$ to $H+2\equiv1\bmod3$. A second such step cannot follow.
These are universal necessary conditions, not properties of one
attachment grammar.

A minimum same-parity gap is attained: the strong order-seven rows
$(112,81,35,7,12,26,60)$ and $(104,81,35,6,13,26,60)$ have
$(H,Q(2))=(153,1812),(159,1812)$. Flip arc $0$--$4$ first to get
$(159,1836)$, then $0$--$3$ to get the final state; the intermediate
is strong as well. Thus a literal $+6$ change followed by a
Hamilton-neutral repair restores the two-slot count once.
Appending $L$ universal sinks to each source and cyclically wrapping
produces strong exact pairs $(1812\cdot2^L+153,1812\cdot2^L+159)$.

## Repeated cyclic wrapping forces increasingly incompatible divisibility

Let $p(t)=Q_+(t)-Q_-(t)$, $d=p(1)$, and use the
[cyclic wrap operator](../session316_hinge_and_signed_repairs_20261113/CYCLIC_WRAP.md)
$\mathcal Lp(t)=t^2p(t)+t(p(t+1)-p(t-1))$.
If $\mathcal L^j p(2)=0$ for $0\le j\le r$, the values through
$t=r+2$ are forced to be $p(t)=c_t d$, where
$$
 c_1=1,\quad c_2=0,\quad
 c_{j+1}=c_{j-1}+(1-j^2)c_j/j.
$$
Indeed $(\mathcal L-I)^j p(2)=0$ by binomial expansion, and the
coefficient of the furthest value $p(j+2)$ is $2\cdot3\cdots(j+1)$,
nonzero. The branch at 0 terminates since $p(0)=0$.
The displayed sequence, extended with $c_0=0$, satisfies
$\mathcal Lc=c$ at the positive integers, so it gives the unique
forced values by triangularity. For $j\ge3$ its nonzero signs
alternate and magnitudes grow at least by factor $j-1$ per step.
No nonzero $d$ can obey this forever with a fixed polynomial.

At the first two returns this already gives a sharper integer filter.
For nonempty equal-parity sources write
$p(t)=t^\delta(t^2-4)R(t^2)$ with $\delta=1$ odd and $2$ even.
Then $R(1)=-d/3$, $R(9)=d/(5\cdot3^\delta)$ after one return,
and $R(16)=-8d/(36\cdot4^\delta)$ after two. Newton slopes
and second divided differences at integer nodes are integers for an
integer polynomial. The corresponding necessary multiples of $d$
are $(60,90)$ after one return and $(3780,7560)$ after two,
for odd and even order respectively. For the same actual order,
the leading monomial cancels and all other packing coefficients
are even, so apply the test to $p/2$. Including the basic parity
requirement, the same-order multiples become $(120,90)$ and
$(7560,15120)$; the even one-return value stays 90 because the
least multiple from $p/2$ in that case is 45 before imposing parity.
This rules out repeatable cyclic readiness for the gap-six pair.


## A literal changing-interface bypass at one exact target

A continuation need not retain all inputs separately if the next count only sees their combination. The cyclic wrapper reads the sum $H+Q(2)$ for its count. The strong order-seven tournaments with rows
$$
(98,92,17,5,9,30,60),\qquad (82,36,17,7,10,29,62)
$$
have responses $(H,Q(2),Q(3))=(99,1332,9153)$ and $(103,1328,8877)$. Their wraps both have exactly 1431 paths, while their new $Q(2)$ values are respectively 23436 and 22860. Thus two different source counts can return to **the same exact target** with distinct continuation resources, by allowing the input $Q(2)$ to compensate. This is a change of the preserved functional, not a violation of the rigidity of a fixed $Q(2)$ fiber.

The source gap four is minimal for this particular neutral arrival: $Q_T(2)$ is divisible by four for every tournament of order at least two, since complementing the two slot subsets pairs $2^{n-1}$ terms $H(S)H(T-S)$, each odd. Equal $H+Q(2)$ therefore requires a count difference divisible by four. For even order at least four it requires a multiple of 24: the odd-cycle formula gives $8\mid Q(2)$ (no-packing term $2^n$, a single odd cycle leaves at least one unused vertex, and multiple cycles contribute at least four powers of two). Modulo three, on even orders $Q(2)\equiv H$, so neutrality of $H+Q(2)$ also forces $3\mid\Delta H$. The odd seven-vertex pair attains the small-gap bound. It does not claim an infinite sequence of such arrivals.

The companion replay verifies all literal tables and wraps by exact subset dynamic programming, and independently computes the integer interpolation constants.
