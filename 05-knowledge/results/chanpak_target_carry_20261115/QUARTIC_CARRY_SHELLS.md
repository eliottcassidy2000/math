# Target-selected fine phase in quartic carry shells

A physical bracketing theorem, not exact universal placement.
The [previous two-tape carrier](../chanpak_absolute_steering_20261115/TWO_TAPE_STEERING.md)
has count F(x,y) = B + Qxy + Tx + Sy + 45[p(x)+p(y)], with
(B,Q,T,S) = (2855,5760,4034,3914). Its exact binary surcharge satisfies
p(0)=0, p(2b+1)=p(b), p(2b)=p(b)+b+(b mod 2), and 0≤p(b)≤b.

**Theorem.** For each fixed integer M≥1, powers of two q and integers
⌈q/200⌉≤u≤⌊q/100⌋ index explicit intervals of width Θ_M(q²) at height
Θ_M(q⁴), disjoint at fixed sufficiently large q. Every odd N in their
union has physical bracketing counts F₋≤N≤F₊, both congruent to N modulo
2M and within O_M(q)=O_M(N^(1/4)). This treats Θ_M(q³) explicit odd targets
per scale. It adds arbitrary fixed-modulus phase to the prior quartic
shell result at the same error exponent. It makes no assertion about
targets outside the intervals.

## The target chooses the phase while the complementary carry stays fixed

Write M=2^a m, with m odd. Fix g=gcd(45,M) and A=2^k≥max(2^(a+1),g).
Use

$$
K=2qu+q-1,\quad X=qu+v,\quad Y=qu+q-1-v,\quad
x=AX+d,\quad y=AY+d,\quad 0\le v<q.
$$

The low q-word of X+Y is all ones and splits without internal carries.
The appended common address d is chosen only as the target requires.
For c_d=A-1-rev_k(d) and δ_d=[d<A/2], the binary split gives exactly

$$
p(x)+p(y)=c_d K+2p(d)+\delta_d+p(X)+p(Y).
$$

Prescribe X≡0 mod m, so Y≡K mod m. All nonvariable terms of F modulo
2M can then be computed at X=0,Y=K: A kills the dyadic component, and
the X-residue accounts for the odd one. Call the resulting odd integer

$$
H_d=B+Qd(AK+d)+Td+S(AK+d)+45[c_d K+2p(d)+\delta_d].
$$

Choose d in [0,g-1] such that B+SAK+(T+S)d≡N mod g. This is possible
because g divides Q, and gcd(T+S,45)=1. Then 2g divides N-H_d, so

$$
45z\equiv (N-H_d)/2\pmod M
$$

has a solution. In each aligned v-block of length D=2^L, where L=L(M)
is the fixed length of the [complementary compiler](../chanpak_absolute_steering_20261115/PAIR_PROGRAMMING.md),
the high prefixes of X,Y sum to (K+1)/D-1, and their D-addresses
are complementary. The compiler can therefore program X mod m and
(p(X)+p(Y))/2 mod M together, assigning the latter z. This gives one
physical program in every full aligned block, with F≡N mod 2M.

## Absolute localization

The key split cancellation is

$$
p(X)+p(Y)=u(q-1)+2p(u)+(u\bmod 2)+p(v)+p(q-1-v).
$$

Remove the last two terms times 45 from F, calling the remaining
quadratic in v G_d(v). Then 0≤F-G_d(v)≤45(q-1); its smooth one-step
change on q/2≤v≤q-1 is O_M(q). In particular
G_d(q/2)-G_d(q-1)=QA²(q²-2q)/4-O_M(q). Uniformly in v,
|G_d(v)-G_0(v)|≤2QA²K+O_M(1): the product shifts by Q(AdK+d²),
the linear terms by (T+S)d, and c_d changes by at most A-1;
the other differences depend only on M. Since K≤q²/50+q,
every d shifts the model by at most QA²q²/25+O_M(q).

Take the explicit common target interval

$$
I_{q,u}=[G_0(q-1)+QA^2 q^2/16,\ G_0(q/2)-QA^2 q^2/16].
$$

For sufficiently large q, D divides q/2 and the 1/16 margins
dominate the possible shifts, the O_M(q) checksum, and the endpoint
loss of at most D−1 positions. A programmed witness in the first
block of [q/2,q) thus exceeds N, while one in the last lies below it.
Bisection on actual counts preserves these opposite signs without
assuming monotonicity; it ends in adjacent blocks. Their position
difference is less than 2D, and their count difference is O_M(q).
The algorithm invokes the compiler at O_M(log q) block positions.

The interval width is QA²q²/8-O_M(q). Increasing u by one increases
the center by Θ_M(q³): the product contributes
QA²[q²(2u+1)+q(q-1)], and the other differences are O_M(q), using
p(u),p(u+1)≤u+1. This proves disjointness at fixed q, and the stated
height and size.

## Exact inversion: the phase is not enough

For a fixed d and target N in this common interval, an exact witness
needs N-45(q-1)≤G_d(v)≤N. It must satisfy v-q/2≥q/8 for large q: the
upper guard and shift bound require a drop from G_d(q/2) of at least
(1/16-1/25)QA²q²-O_M(q), while the drop at q/8 is only
QA²q²/64+O_M(q). Here the quadratic is decreasing and each step
drops at least QA²q/4-O_M(1)>45q. There is at most one candidate v
on this half-shell, found by a quadratic cutoff followed by an
exact checksum test. Recognition spends no all-residue budget;
it is a partial inverse, not an existence proof.

There are only O_M(q) exact outputs in this alphabet per shell and
at most O_M(q²) at its whole indexed scale, even with all g common
phases. The bracketing theorem treats Θ_M(q³) target positions.
That gap in exponents shows why programming more residues alone
cannot make this low all-ones carry family universal.
