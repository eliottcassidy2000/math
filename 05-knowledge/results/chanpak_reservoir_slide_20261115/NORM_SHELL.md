# A shifted square-shell decoder for two growing repair banks

The [wrapped reservoir law](RESERVOIR_CHARGE.md) can turn two matched
intervals into an exact shifted product. This makes a precise
analogy with square-sum and primitive-triple constructions: the
division problem becomes a shell with a discriminant and a bounded
inverse, not just a collection of apparently related numbers.

Use two independent copies of the [unit slide](UNIT_SLIDE.md),
at the same height \(k\), and one marked singleton of neutral
padding on each side inside the banks. Let \(U\) be the upper
apex count at that height and \(S\) their matched marked start.
The added singletons do not change any member's factor. The
outer construction gives
\[
N=(U-2i)(U-2j)+C,\qquad
C=[2(45+S)]^2,\qquad 0\le i,j\le k.
\]
The two subtracted digits come from independent physical slide
positions, not arbitrary arithmetic factors.

## Inversion without factoring

Let \(P=N-C\). Suppose \(U>2k^2\), as holds for the unit slide
at every height (its guard cycles alone give \(U\ge 1+2^{12k}\)
at \(k\ge1\); at \(k=0\) the inequality is trivial). Put
\[
u=\left\lceil\frac{U^2-P}{2U}\right\rceil,\qquad
v=\frac{P-U^2+2Uu}{4}.
\]
Any working pair has \(u=i+j\) and \(v=ij\): from the product
formula the rounding error is exactly \(2ij/U\in[0,1)\).
Conversely test that \(v\) is an integer, that
\[
D=u^2-4v
\]
is a nonnegative square \(m^2\), that \(u\equiv m\pmod2\),
and that \(i=(u-m)/2,\ j=(u+m)/2\) lie in \([0,k]\).
These tests are necessary and sufficient and return the
literal two slides. No factoring of the large shifted target
and no search through the \(k+1\) variants are needed.

Equivalently the candidates lie on the exact equation
\[
N+(i-j)^2=(U-i-j)^2+[2(45+S)]^2.
\]
This is the transferable geometry: first choose the sheet
(the sum \(i+j\)), then the residual must be a square with an
actual bounded digit inverse. Even at the correct sheet and
parity, most residuals do not pass that square test.

For fixed \(u\) the output is a translate of the negatives of
squares of the allowable differences. Each difference has
the same parity as \(u\); two successive allowed nonnegative
differences \(m,m+2\) change the count by \(4m+4\).
Distinct sums \(u\) lie in disjoint bands separated at scale
\(U\), because the entire \(4ij\) term is at most \(4k^2<2U\).
Thus **exposing a repair interval as a factor changes its
geometry**. Two growing factor intervals do not give a
near-square growing interval of products. An external
variable charge moves these bands without adding independent
digits inside one fixed band.

The [replay](replay_norm_shell.py) checks the inverse against
the complete small factor windows over many heights and a
literal height-one target from their exact transfer profiles.
This solves target decoding **given** the height and the
charge; it does not prove there is a physical choice of both
for every target. The asymmetric tagged construction keeps
the same shifted product but avoids requiring a square charge.
