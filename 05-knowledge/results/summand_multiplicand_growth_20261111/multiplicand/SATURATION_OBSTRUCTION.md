# Independent multiplicands never fill their numeric frontier

Status: **PROVED** all-round statements; **FINITE-EXACT** initial-state
controls in `saturation_receipt.json`. This concerns synchronous two-parent
closure, not ordinary directed reachability.

Let $2\le a<b$ be multiplicatively independent: $a^u b^v=a^{u'}b^{v'}$
for nonnegative integer exponents implies $(u,v)=(u',v')$. Put

\[
S_0=\{a,b\},\qquad
S_{t+1}=S_t\cup\{xy:x,y\in S_t,\ x<y\},\qquad M_t=\max S_t.
\]

The full strict closure is

\[
C=\{a,b\}\cup\{a^u b^v:u,v\ge1\}.
\tag{1}
\]

This is the two-seed instance of the sibling
[mixed-support closure theorem](../summand/CLOSURE.md). Directly, no pure
higher power can be generated from one seed type; every mixed monomial is
obtained by starting with $ab$ and adjoining the remaining seed factors.
All these later factors differ from the running mixed product.

## Exact least missing element at every round

**Theorem.** For every $t\ge0$, the least element of $C$ absent from $S_t$ is

\[
\boxed{L_t=a^{t+1}b}.
\tag{2}
\]

For every $k\ge0$, the exact first round containing $a^k b$ is $k$.
Indeed its construction has exactly one $b$ leaf. At any product node on
the path from that leaf, the other child has no $b$ leaf and therefore must
be the lone pure seed $a$. Thus the mixed spine adds exactly one $a$ factor
per round. The successive products $b,ab,a^2b,\ldots$ attain that bound.
Multiplicative independence is what makes the leaf counts intrinsic.

Conversely a mixed monomial $a^u b^v$ with $u+v\le t+1$ is available by
round $u+v-1\le t$, by the sequential construction. If $u+v\ge t+2$, then

\[
a^u b^v=a^{u+v-1}b\,(b/a)^{v-1}\ge a^{t+1}b.
\]

Thus every closure element smaller than (2) is present, while (2) itself
requires round $t+1$. This proves the exact numeric first-hole law.

## Permanent failure of saturation below the live maximum

The first three nontrivial maxima are

\[
M_1=ab,\qquad M_2=ab^2,\qquad M_3=a^3b^3.
\tag{3}
\]

For the last identity, the two largest members of
$S_2=\{a,b,ab,a^2b,ab^2\}$ are $a^2b$ and $ab^2$.
Their product is the largest possible round-three product.

Rounds zero and one are saturated relative to $C$: $S_t=C\cap[1,M_t]$.
At round two this remains true exactly when $b<a^2$. Equality is excluded
by independence. If $b>a^2$, the first hole $a^3b$ already lies below
$M_2=ab^2$.

At round three,

\[
\frac{M_3}{L_3}=\frac{b^2}{a}>1.
\]

Moreover $ab$ remains live and $M_t>ab$ for $t\ge2$, so
$M_{t+1}\ge abM_t$, while $L_{t+1}=aL_t$. Therefore for every $t\ge3$,

\[
\boxed{L_t<M_t},\qquad
\boxed{\frac{L_t}{M_t}\le\frac{a}{b^{t-1}}\longrightarrow0}.
\tag{4}
\]

There is no eventual saturation of the numeric frontier. Indeed the
first unavailable closure element becomes an arbitrarily small fraction
of the largest available one.

For $\{a,b\}=\{2,3\}$ the complete first states are

    S0 = {2,3}
    S1 = {2,3,6}
    S2 = {2,3,6,12,18}
    S3 = {2,3,6,12,18,24,36,54,72,108,216}.

The first missing closure values are $6,12,24,48$ respectively. At round
three, 48 is missing below the live maximum 216. At the next round the
first hole is 96 while the maximum is 23328.

## A positive amplifier and the exact lost coordinate

Any distinct live pair $(x,y)$ supplies, in two further rounds,

\[
(x,y)\longmapsto(x^2y,xy^2).
\tag{5}
\]

First form $xy$, then multiply it separately by $x$ and $y$. The ratio is
preserved, while logarithmic sum grows by three. In log coordinates the
matrix is $\left(\begin{smallmatrix}2&1\\1&2\end{smallmatrix}\right)$,
with eigenvalues 3 and 1. Repeating (5) gives at round $2j$ the pair

\[
\left(a^{(3^j+1)/2}b^{(3^j-1)/2},
      a^{(3^j-1)/2}b^{(3^j+1)/2}\right).
\]

Large mixed products therefore have a balanced parallel construction,
while the one-$b$ boundary has an unavoidable serial spine. This is the
mechanism behind the gap; it is not merely nonuniform spacing of logarithms.

Multiplicative dependence can destroy the obstruction. With seeds
$\{2,4\}$, the putative one-$b$ target $a^3b=32$ also equals $ab^2$ and
is already live at round two, as $4\cdot8$. Leaf exponents are no longer
uniquely determined by the integer. The prime-power valuation chart then
reduces the dynamics to scalar additive closure, where the sibling eventual
saturation theorem applies. For independent seeds the two exponent
coordinates remain faithful, and (2) proves the opposite behavior.

## Complete two-seed saturation dichotomy

For any two distinct integers `2<=a<b`, eventual saturation relative to
their full strict product closure holds **if and only if a and b are
multiplicatively dependent**.

Only the dependent direction remains. Dependence means their prime-exponent
vectors are positively proportional. Divide the vector of a by the gcd
of its entries to obtain a primitive nonnegative integer vector v. The
vector of b is also an integer multiple of v: Bezout's identity for the
entries of v shows that any rational scalar making all its entries
integers must itself be an integer. Hence there is an integer c>=2 and
positive integers r<s with `a=c^r`, `b=c^s`. The map `k -> c^k` preserves
order, distinct parents, and each synchronous sum/product generation.
The additive saturation theorem for seeds `{r,s}` therefore gives exactly
the required product saturation. Necessity is the independent obstruction
already proved above. The common base c need not be prime; `{6,36}` is
one such dependent example. No assertion of necessity for arbitrary
larger, dependent seed sets is made here.

Run `python3 -B 05-knowledge/results/summand_multiplicand_growth_20261111/multiplicand/saturation_check.py`.
It builds only four rounds for four small independent seed pairs and the
dependent hostile $\{2,4\}$; all-round scope comes from the proof above.

## Extension to any finite independent seed set

Let a finite set of at least two integers >=2 be jointly multiplicatively
independent, and let a<b be its two smallest values. The exact least
missing value is still `a^(t+1)b` at every round t. Indeed independence
forces the same one-b branch for that target. Every mixed monomial of
total degree at most t+1 is built sequentially by round t, beginning
with its two largest distinct seed types. Every higher-degree mixed
monomial has value at least `a^(t+1)b`, since at least one factor is
at least b and all other factors are at least a. Seed values themselves
are already present, regardless of their size.

The two-seed subsystem supplies the round-three lower bound for the
maximum, and multiplying a live maximum by ab gives the same subsequent
bound as in (4). Thus additional independent seed directions do not
restore eventual saturation. The explicit formulas for the first three
maxima and the round-two criterion above remain two-seed statements.
The small checker covers the two-seed examples; this extension follows
from the displayed proof and is not described as an exhaustive census.
