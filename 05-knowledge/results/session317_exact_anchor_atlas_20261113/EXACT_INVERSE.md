# An exact square-root search for a prescribed absolute count

The [terminal law](TERMINAL_ATLAS.md) supplies a target-aware inverse,
not just a forward enumerator. Fix a physical terminal profile
$(h,E,s,c)$ with $E>0$; write $A=h+E$. To decide whether a prescribed
integer $N$ equals $H(Pu)$ for some binary prefix, one can find **all**
such prefixes, up to leading zeroes, in
$O_u(\sqrt{N+1}\log^2(N+2))$ integer operations. This is a
single-terminal result, not a search across all tournaments.

The two-sided coordinate is literal. For $B=2^mX+Y$ with
$0\le Y<2^m$, let $R_m(Y)$ reverse exactly those $m$ bits and let
$\eta=\lfloor Y/2^{m-1}\rfloor$. As derived in the
[cut-one work](ONE_RETURN_STAR.md), the binary defect satisfies
$$
 K(B)=K(X)+K(Y)+X R_m(Y)+\eta(X\bmod2).
$$
Therefore at fixed $Y$ the sought count is
$$
 [c+AY-hK(Y)-s(Y\bmod2)]
 +[A2^m-hR_m(Y)]X-hK(X)-h\eta(X\bmod2).
$$
This is exact, including the one-bit seam term. Since $K(B)\le B$
by induction and $c\ge s$, any witness has $B\le L=\lfloor N/E\rfloor$.

For the monotonicity needed to invert this expression, write
$X=2^rZ+2^r-1$ with $r$ the number of trailing 1s and $Z$ even.
Then $K(X+1)=K(Z+1)$, while $K(X)\ge K(Z)$, and if $Z=2D$
the increment $K(Z+1)-K(Z)$ is $D+(D\bmod2)\le X+1$.
(The inequality $K(X)\ge K(Z)$ also follows directly by appending
$r$ ones.) The positive coefficient in brackets is at least
$E2^m+h$. Choose $m\ge1$ so that
$E2^m\ge h(\lfloor L/2^m\rfloor+2)$, which needs
$2^m=O_u(\sqrt{N+1})$. Then on
$0\le X\le\lfloor(L-Y)/2^m\rfloor$, the expression is
strictly increasing: the $K$ increment and the final parity
correction cost at most $h(X+2)$. For each low word $Y$,
binary search gives the exact high word $X$ or rejects.
Evaluating $K$ uses $O(\log(N+2))$ integer operations.
The same search finds the neighboring values in this truncated
prefix range if a target is absent. Requiring strong outputs
amounts to the terminal reachability guard and $B>0$; it does not
change the arithmetic search.

The older [golden denominator filter](../denominator_arithmetic_filter_20261004.md)
distinguishes reading an arithmetic phase from realizing its
integer. This is the tournament version: low bits choose the
reversed slope, high bits are recovered by a monotone solve, and
the exact terminal graph still has to be supplied. No residue or
continuous limit is mistaken for an absolute witness.

## The orthogonal inverse: select a prefix and solve for the terminal

Instead, fix $B\ge1$ and a proposed reservoir count $h$. Put
$D=N-h(B-K(B))$. The terminal equation is $D=EB+c-s\ell$.
Since $0\le s\le h$ and $h\le c\le E+s$, every possible end
response lies in the finite interval
$$
 \left\lceil\frac{D-h(1-\ell)}{B+1}\right\rceil
 \ \le E\le\
 \left\lfloor\frac{D-h(1-\ell)}{B}\right\rfloor .
$$
For the lower bound use $D\le E(B+1)+s(1-\ell)$; for the upper
use $D\ge EB+h-s\ell\ge EB+h(1-\ell)$. Its width is
$O(N/B^2+1)$, since $K(B)\le B$. In particular there are only a
bounded number of choices of $E$ per prefix in a square-root-sized
window $B\asymp\sqrt N$. At fixed $h$ each such $E$ has only a
polylogarithmic normalized terminal search, by the cut decoder in
the terminal note. Test the exact equation and the strongness
guard after decoding; one can also insist $W-b$ be transitive with
$b$ in the out-cut to stay in the two-apex word class.
This narrows the next construction problem to moving between
actual reservoir counts $h$, not adding more congruence filters.
