# Consecutive successor rows, exact root decisions, and retained norm data

**Status:** PROVED elementary results below; FINITE-EXACT declared experiments. The general nonexistence of integer rows at every power at least five is not proved here.

**Source and ownership:** the user's two consecutive-sum pyramids, 2026-10-14. Companion [checker](../../04-computation/experiments/pyramid_successor_rows_20261014b.py) and [saved JSON](pyramid_successor_rows_20261014b.json). This packet repairs the universal unit-width bracket in [HYP-2454](../hypotheses/HYP-2454-triangular-power-balance-towers.md) and [HYP-2457](../hypotheses/HYP-2457-faulhaber-anchor-expansion.md); their earlier finite computations and exact low-power identities remain valid.

## 1. Inheritance and the marked row

The closest inherited mechanism is the cancellation of all even midpoint moments in HYP-2454/2457. The corrected near miss is extrapolating a bracket checked only for powers 3 through 8 to every power. The least-used useful sidecar is the ordered endpoint address: it distinguishes equality of powers, equality of unsquared sums, and exact coincidence of raw blocks. The related Pell endpoint program in HYP-2456 remains separate.

The live objects are row width, centre, power, odd moments, integer lattice, and marked quadratic norm. The anchor is the two displayed towers; the niche is a complete exact decision procedure for each supplied higher-power row; the wildcard is a reversible rational norm map and its arithmetic failures.

Fix integers $n,p\ge1$. The left side has $n+1$ terms, the right side $n$ terms, with centre $c$ and positive start $a=c-n$:

\[
D_p(c,n)=\sum_{i=0}^{n}(c-i)^p-\sum_{i=1}^{n}(c+i)^p.
\tag{1}
\]

The source is $(n,p,c)$, not its common total alone. A balanced row requires $c>n$ and $D_p(c,n)=0$. The compiler retains all three coordinates and checks the actual endpoints.

## 2. The two exact towers and their successor rules

Write $u=n(n+1)$. Direct expansion gives

\[
D_1=c-u,\qquad D_2=c(c-2u).
\]

Consequently these are the unique positive rows for powers one and two:

| Power | Start $a$ | Centre $c$ | Common total |
|---|---:|---:|---:|
| 1 | $n^2$ | $n(n+1)$ | $n(n+1)(2n+1)/2$ |
| 2 | $2n^2+n$ | $2n(n+1)$ | $n(n+1)(2n+1)(12n^2+12n+1)/6$ |

The first row covers the entire square shell $[n^2,(n+1)^2-1]$, split after $n^2+n$. Thus the first three rows are exactly the user's $1+2=3$, $4+5+6=7+8$, and $9+10+11+12=13+14+15$. The square starts are $3,10,21,36,55,\ldots$.

For $p\in\{1,2\}$ the complete successor rule is

\[
a_p(n)=pn^2+(p-1)n,\qquad
a_p(n+1)-a_p(n)=2pn+2p-1.
\tag{2}
\]

It retains $n$ and increases the two side lengths by one. Its second difference is $2p$. The common totals have third difference 6 for $p=1$ and fifth difference 480 for $p=2$. These are polynomial recurrences, with no search for a new row.

The map $(n,c=u,p=1)\mapsto(n,2c,p=2)$ is an exact map between the two classified families. It changes the tested predicate. Forgetting the squares in the second family leaves ordinary side difference $u$, not zero; the sum of its two unsquared sides is four times the common total of the first family.

Only one entire raw block is shared between these two families:

\[
Q_L(3)=F_R(4)=[21,22,23,24].
\]

For the left-left and right-right comparisons, equality of lengths forces the same index, and the starts differ. The right-square/left-linear comparison forces linear index $n-1$ and incompatible starts. The remaining comparison forces linear index $n+1$; its start equation is $n^2-2n-3=0$, hence $n=3$. This is an all-index proof, independently checked through index 64. The smaller blocks $10,11,12$ and $13,14$ have equal sums of squares, 365, but their ordinary sums are 33 and 27.

## 3. A uniform replacement for the false unit-width bracket

Put $S_j(n)=\sum_{i=1}^n i^j$. The inherited odd-moment identity is

\[
D_p(c,n)=c^p-2\sum_{j\text{ odd}}\binom pj S_j(n)c^{p-j}.
\tag{3}
\]

For $c>0$, divide by $c^p$. The resulting function is strictly increasing from $-\infty$ to 1, because every subtracted coefficient is positive. Hence there is exactly one positive real root. For $p\ge3$ it satisfies the **proved uniform bound**

\[
pu<c<p(u+1).
\tag{4}
\]

At $c=pu$ the $j=1$ term equals one, while the $j=3$ term is positive, so the defect is negative. To prove the upper sign, use $2S_j\le u n^{j-1}$ and $\binom pj\le p^j$. At $c=p(u+1)$,

\[
2\sum_{j\text{ odd}}\binom pj\frac{S_j}{c^j}
\le \frac{pu/c}{1-(pn/c)^2}
=\frac{u(u+1)}{(u+1)^2-n^2}<1.
\]

The denominator exceeds the numerator by $u+1-n^2=n+1>0$. The geometric series is legitimate since $n/(u+1)<1$. In particular the unique root exceeds $n$, so it is the unique positive-start row candidate.

The former claim $pu<c<pu+1$ is false. Exact counterexamples are

\[
D_{15}(181,3)=-61355600116348100287252240367,
\quad D_{16}(33,1)=-2148849686515840130559.
\]

Their roots lie respectively in $(181,182)$ and $(33,34)$, beyond the proposed upper endpoints. The first failure in the declared rectangle $3\le p\le32,1\le n\le128$ is $(p,n)=(15,3)$; this finite ordering is not a global minimality theorem. The old formal correction $(p-1)(p-2)/(12p)$ already exceeds one at $p=15$. No uniform asymptotic remainder is needed for the exact counterexamples or (4).

### A complete decision for each supplied row

For fixed integer $n,p$, (1) is a **monic integer polynomial in $c$**. The rational-root theorem therefore makes every rational root an integer. Moreover $D_p(c,n)\equiv c^p\pmod2$, so an integer root must be even.

The checker computes moments by the exact telescoping recurrence

\[
S_0=n,\qquad
S_j=\frac{(n+1)^{j+1}-1-\sum_{r<j}\binom{j+1}{r}S_r}{j+1}.
\tag{5}
\]

It then binary-searches the integer interval in (4). It either supplies the exact balanced centre or adjacent integers with strictly negative/positive defects, proving the unique root irrational. The search needs $O(\log p)$ defect evaluations after the moments, with $O(p^2)$ arithmetic operations for (5); operand sizes still depend on $p$ and $\log n$. This is not a constant-bit-cost claim. The receipt recomputes the source and signs before acceptance.

For powers three and four the complete irrationality result is elementary. Their nonzero root equations are

\[
D_3=c^3-3uc^2-u^2/2,\qquad
D_4=c(c^3-4uc^2-2u^2).
\]

The signs at $3u,3u+1$ are $-u^2/2$ and $(17u^2+12u+2)/2$; those at $4u,4u+1$ are $-8u^3$ and $(4u+1)(14u^2+8u+1)$. Thus no rational positive row exists for $p=3,4$. The finite rectangle finds none for $5\le p\le32$ either, but does not establish all-power nonexistence.

## 4. Changing the split is a different exact problem

The one-extra-left condition is essential. Let the earlier block have $m=n+d$ terms, the following block $n$ terms, and positive integer start $a$. If $m\le n$, equality is impossible because later terms are larger. For $d>0$, the power-one condition is exactly

\[
a=\frac{n^2}{d}-\frac{d-1}{2},
\]

with positivity and integrality required. For squares set

\[
L=d^2-d-2n^2,\qquad C=2S_2(m-1)-S_2(m+n-1).
\]

The defect is $da^2+La+C$. Thus the compiler tests $\Delta=L^2-4dC$ for being a nonnegative square and tests

\[
a=\frac{-L+\sqrt\Delta}{2d}
\tag{6}
\]

for positive integrality. This is complete: for positive real $a$, the ratio of the earlier power sum to the later power sum strictly increases. Its logarithmic derivative is the difference of two weighted averages of $p/(a+i)$; every earlier index is smaller than every later index. The ratio tends to $m/n>1$, so there is at most one positive crossing, and any crossing is the plus root in (6).

For example,

\[
18^2+19^2+\cdots+34^2=35^2+36^2+\cdots+42^2
\]

has 17 versus 8 terms. It is outside the user's consecutive-width tower. Counts are retained rather than silently identifying these distinct families.

## 5. A reversible norm map, with its lattice sidecar

For a square row form the ordered vectors

\[
L=(a,a+1,\ldots,c),\qquad R=(c+1,\ldots,c+n,0).
\]

They have the same squared Euclidean norm. Their difference has $n$ entries $-(n+1)$ and last entry $2n(n+1)$. Therefore

\[
v=L-R,\quad v\cdot v=n(n+1)^2(4n+1),\quad
H_n=I-\frac{2vv^T}{v\cdot v}
\tag{7}
\]

is a rational orthogonal involution carrying $L$ to $R$ and back. This gives an explicit reusable operation, not merely equal scalar totals: the source is the ordered padded vector, the preserved predicate is squared norm, and $n$ plus the lattice and positivity conditions are retained sidecars.

At $n=1$,

\[
H_1=\begin{pmatrix}3/5&4/5\\4/5&-3/5\end{pmatrix},
\qquad H_1(3,4)=(5,0).
\]

But $H_1(1,0)=(3/5,4/5)$ is not integral and $H_1(0,1)=(4/5,-3/5)$ is not positive. Thus the map is not a general integer-positive recursion. The row $(3,4,5)$ is the unique consecutive positive Pythagorean triple: expanding $a^2+(a+1)^2=(a+2)^2$ forces $a=3$. It is primitive, whereas the primitive triple $(5,12,13)$ is not consecutive. No claim identifies all primitive triangles with the row tower.

The number 1729 supplies a different, cubic predicate:

\[
1^3+12^3=9^3+10^3=1729.
\]

The inherited [THM-463, two-cube divisor reconstruction](../../01-canon/theorems/THM-463-two-cube-representations-are-a-divisor-property-on-the-split-axis.md) retains $s=x+y\mid N$, reconstructs $xy=(s^3-N)/(3s)$, and checks the quadratic discriminant. The two positive pairs have $s=13,19$. The factorization $(x+y)(x^2-xy+y^2)$ uses an Eisenstein norm, not the Gaussian norm $x^2+y^2$. Indeed 1729 is not a sum of two integer squares: it is divisible by 7 but not 49, and $x^2+y^2\equiv0\pmod7$ forces both $x,y$ divisible by 7. Shared numerals therefore do not transport the quadratic predicate.

## 6. Verification, corrections, and scope

Run `python -B -X utf8 04-computation/experiments/pyramid_successor_rows_20261014b.py` and repeat with `-O`. The saved JSON is deterministic. The explicit universes are: moments through degree 16 at $0\le n\le32$; direct defects at $1\le n\le16,1\le p\le10$; all root decisions at $1\le n\le128,1\le p\le32$; successor rows through 64; all split counts at most 10 and starts at most 100, plus the two wider square examples; norm maps through $n=8$; 15 malformed controls. The giant input $n=10^{100},p=17$ checks that moment computation does not enumerate the row.

The all-power result is the unique real root, the exponent-width bracket, and a complete finite decision for each given source. No universal nonexistence at powers at least five, Carmichael implication, Collatz certificate, Jacobian theorem, or plane-colouring construction follows. The demonstrated historical error was extrapolating finite exponent checks to a uniform unit-width bound; the replacement keeps the exponent as part of the exact bound. During checker development, the successor increment was independently corrected to $2pn+2p-1$ and tested by literal row comparison.

Normal and optimized runs agree at **34,734 exact checks**. Saved LF-normalized JSON SHA256: `24e5d0f89e58b642568bf8841bebd8857eed82511dc89d09beb89adcb845b92a`.
