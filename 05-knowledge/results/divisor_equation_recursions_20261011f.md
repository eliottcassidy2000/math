# Divisor records, a reciprocal-equation recursion, and an actual field escape

**Status: PROVED elementary identities and recursion interfaces; FINITE-EXACT bounded controls; CITED sequence identification.** A positive integer solution and its actual triangle rotation are exhibited. Neither a least positive solution nor the full elliptic rank is claimed. Escaping the earlier Moser field does not produce a six-colour equality tester or settle the plane problem.

## 1. Inheritance and the two lists

The closest proved divisor mechanism is DB1--DB3 in [arithmetic braids](arithmetic_braids_20260917_divisors.md): divisors form a prime-exponent box, and a chosen proper-divisor statistic is computed from that box. Its counts exclude endpoints in ways the ordinary divisor function does not. [Divisor layers and pair geometry](divisor_layers_pair_geometry_20261004.md) already separates box cardinality from its incidence geometry. The closest equation mechanism is [prime shells and elliptic curves](prime_shells_20260921_elliptic.md), but that treats the **different** equation with three symmetric fractions summing to four. [Creation numbers](creation_numbers_20260925.md) retains the original denominators and positivity after a cubic change of coordinates. Those guards are essential here too.

The live objects are divisor records, exponent-box incidence, a guarded elliptic chart, a positivity arc, and the exact angle field of a solution triangle. The canonical hostiles are an exponent profile that is not a record, a cleared cubic with a zero original denominator, a valid signed point outside the positive arc, and a unit rotation still defeated by a stronger field screen.

Write $\tau(n)$ for the number of **all** positive divisors, including $1$ and $n$. The strict records of $\tau$ are the highly composite numbers: [OEIS A002182](https://oeis.org/A002182). Their successive divisor counts are [OEIS A002183](https://oeis.org/A002183). The supplied lists have unequal lengths; their exact pairing is

| Record argument | 1 | 2 | 4 | 6 | 12 | 24 | 36 | 48 | 60 | 120 | 180 | 240 | 360 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Divisor count | 1 | 2 | 3 | 4 | 6 | 8 | 9 | 10 | 12 | 16 | 18 | 20 | 24 |

Thus the final count $24$ belongs to the next argument $360$, which was absent from the first list. These are record arguments and record values, not two interchangeable number sequences.

**Complete bounded recursion.** Every record has the form

$$n=2^{a_1}3^{a_2}\cdots p_r^{a_r},\qquad a_1\ge\cdots\ge a_r\ge1. \tag{1}$$

A skipped prime could be replaced by a smaller prime without changing $\tau$; an inverted exponent pair could be exchanged to lower $n$ without changing $\tau$. Either contradicts strict record status. Recursively append an exponent $1\le a_{i+1}\le a_i$ and reject products exceeding a supplied bound $X$. Sort all retained products and keep strict increases of $\prod_i(a_i+1)$. This enumerates every record through $X$: every nonrecord is dominated in divisor count by an earlier record, itself among the candidates. The shape condition alone is not sufficient: $72=2^3 3^2$ has $12$ divisors, already attained at $60$. Nor do successive records form a divisibility chain: $24\nmid36$ and $48\nmid60$.

**Incidence obstruction.** The divisor box of $n$ is $\prod_i\{0,\ldots,a_i\}$. Its cover graph is bipartite, with colour $\Omega(d)\bmod2$. Conversely, the graph joining every comparable pair contains a $K_4$ whenever $\Omega(n)\ge3$, by taking four points of a divisor chain. Four mutually unit-separated points cannot lie in the plane: the two possible third vertices of a unit equilateral triangle are at distance $\sqrt3$. Thus these two direct incidence choices give either a two-colour graph or a non-realizable unit graph. The record count alone does not select useful unit-edge incidence.

## 2. The exact equation and its guarded chart

The equation under discussion is

$$\frac A{B+C}+\frac B{A+C}=\frac C{A+B}. \tag{2}$$

Initially $A,B,C$ are rational and all three pair sums must be nonzero. Positivity and integrality are additional predicates. Put $s=A+B$ and

$$F=s^3+Cs^2-C^2s-C^3-AB(2s+C).$$

Direct expansion gives the denominator-sensitive identity

$$\frac A{B+C}+\frac B{A+C}-\frac C{A+B}
=\frac{F}{(A+B)(A+C)(B+C)}. \tag{3}$$

The tuple $(1,1,-1)$ satisfies $F=0$ but makes the original fractions undefined. It must be rejected. A valid solution also has $2s+C\ne0$: substituting $C=-2s$ into $F$ gives $3s^3\ne0$.

Define

$$x=\frac{4(s+2C)}{2s+C},\qquad y=\frac{12(A-B)}{2s+C}. \tag{4}$$

For every triple in this chart,

$$y^2-x^3-5x^2+32x=\frac{576F}{(2s+C)^3}. \tag{5}$$

Consequently the target is the nonsingular cubic

$$E:\quad y^2=x(x^2+5x-32). \tag{6}$$

Conversely, a finite rational point on (6) with $x\notin\{-4,0,8\}$ gives

$$A:B:C=(8-x+y):(8-x-y):(4x-8). \tag{7}$$

Exactly these excluded finite $x$ values give zero original pair sums. Infinity is also outside this affine chart. Clearing rational denominators, dividing by the integer gcd, and choosing $A+B>0$ gives a unique primitive signed representative. Equations (4) and (7) are inverse on these domains; the script rechecks (2), not merely (6).

**Positivity is a separate interval test.** The representative in (7) is positive exactly when

$$2<x<4\quad\text{on the real locus of (6)}. \tag{8}$$

Indeed $C/s=(2x-4)/(8-x)$ and $(A-B)/s=y/(8-x)$, while

$$ (8-x)^2-y^2=(4-x)(x+4)^2.$$

These prove sufficiency in (8). Conversely, (2) with positive entries implies $C< A+B$, since

$$AB=\frac{(s+C)^2(s-C)}{2s+C}>0,$$

and then (4) gives (8). Also $C>\max(A,B)$: if, for example, $C\le A$, the first fraction in (2) is already at least $C/(A+B)$ and the second is positive. Every positive solution is therefore a nondegenerate triangle.

## 3. A recursive positive construction, with the signed boundary retained

On (6), for distinct nonopposite points $(x,y),(u,v)$ let $m=(v-y)/(u-x)$; for doubling use $m=(3x^2+10x-32)/(2y)$. The resulting point is

$$x'=m^2-5-x-u,\qquad y'=m(x-x')-y. \tag{9}$$

Opposite points add to infinity; infinity is the identity. Substitution verifies preservation of (6). These formulas, followed by the guarded inverse (7), form an exact rational recursion. They do not preserve positivity automatically.

Starting with $P=(4,4)$ gives the boundary tuple $(1,0,1)$. Doubling gives $2P=(36,-228)$ and the valid but signed tuple $(32,-25,-17)$. At $8P$ the recursion gives the positive primitive solution

$$\begin{aligned}
A&=301361533449900458837600,\\
B&=49105016933436320224063,\\
C&=316629033253501281102807.
\end{aligned} \tag{10}$$

The code evaluates all three original fractions exactly. Among the tested multiples $P,\ldots,24P$, precisely $8P,15P,22P$ land in the positive interval. This finite statement does not classify all rational points or all positive solutions. The independent search with $A+B\le999$ has no positive solutions; it is not a global leastness result.

There is also an elementary proof that this recursion produces infinitely many distinct **signed** rational points, without a rank computation. At $3P$, $v_2(x)=-4$ and $v_2(y)=-6$. Whenever $v_2(x)=-2e$ with $e\ge2$, the cubic term uniquely has least valuation in (6), giving $v_2(y)=-3e$. The numerator of the doubling slope has valuation $-4e$ and its denominator valuation $1-3e$, so $v_2(m)=-e-1$. In (9), $m^2$ uniquely has least valuation, yielding

$$v_2\bigl(x(2^j\cdot3P)\bigr)=-4-2j.$$

Thus these points are all distinct. Nothing here claims that these particular doublings remain positive. Swapping $A,B$ is exactly $(x,y)\mapsto(x,-y)$; forgetting that ordering loses a marked side of the triangle.

## 4. Actual angle transport and the next field screen

For the positive triangle, Heron's identity and (2) give

$$16\mathcal A^2=(s^2-C^2)\bigl(C^2-(A-B)^2\bigr)
=(s^2-C^2)^2\frac{2s+3C}{2s+C}. \tag{11}$$

For (10), the reduced last ratio is

$$\frac{23671178496485557}{14590865133564723}.$$

The product of its numerator and denominator is $D\cdot69739671^2$, where

$$D=3\cdot1237\cdot10193\cdot11329\cdot165713
=71013535489456671. \tag{12}$$

The displayed factors are distinct primes (checked by complete trial division). The unit rotation through the actual triangle angle opposite $C$ is

$$z=\frac{-4540156681460417}{19131021815025140}
 +\frac{69739671}{19131021815025140}\sqrt{-D},\qquad z\bar z=1. \tag{13}$$

This is an exact metric transport: its real part is $(A^2+B^2-C^2)/(2AB)$ and its positive imaginary part is $2\mathcal A/(AB)$. It is not an unrelated rational parametrization of the unit circle.

Since $D$ is squarefree and differs from $3$ and $11$, the imaginary quadratic field $\mathbf Q(\sqrt{-D})$ is not one of the three quadratic subfields of $L_0=\mathbf Q(\sqrt{-3},\sqrt{-11})$ (the third is the real field $\mathbf Q(\sqrt{33})$). Thus (13) genuinely leaves the exact field screened in [unit equality candidates](unit_equality_candidates_20261011e.md). The construction preserves unit norm and the original triangle angle; it does not preserve or supply a graph's forcing relation.

**It does not escape every known colouring screen.** Here $D\equiv7\pmod8$, but at $11$ both $-3$ and $-D$ are nonsquares. The companion [field-escape recursion](unit_field_escape_recursion_20261011f.md) retains this marker: in $L_0(\sqrt{-D})$, conjugation acts as the unramified quadratic Frobenius together with the ramified $\sqrt{-11}$ inertia flip. The inherited residue-colouring mechanism, **THM-4558, [residue colourings of field planes](../../01-canon/theorems/THM-4558-residue-colourings-of-field-planes-the-heegner-tower-is-at-most-four-chromatic.md), Lemma R**, then gives the applicable five-colour screen. This local-field consequence is an inherited/companion result, not established by the two finite Legendre checks alone. Consequently merely adjoining (13) to the old field cannot furnish an actual unit-separated six-colour equality tester.

## 5. Transfer contracts and next usable obligations

| Source | Map and preserved predicate | Lost information / retained sidecar | Cheapest decisive test |
|---|---|---|---|
| Divisor record | Prime-exponent recursion; exact $\tau$ | A count forgets primes, divisibility incidence and metric coordinates | $72$ has the correct profile but is not a record |
| Exponent box | Cover graph, parity $\Omega\bmod2$ | Has only two-colour complexity; comparability changes the graph | Four-point chain gives impossible planar unit $K_4$ |
| Original rational equation | (4)--(7), preserving the guarded equality | Must retain nonzero pair sums and the positivity interval | $(1,1,-1)$; $2P$ signed output |
| Positive triangle | Actual rotation (13), preserving its angle and unit norm | One angle is not a unit graph or a forcing certificate | Rational Cayley parametrization has a different real coordinate |
| Old field plus new rotation | Algebraic closure in the joined field | Escaping one field does not escape another residue bound | The two nonsquare markers at $11$ |

The useful construction output is (10) together with its source-owned rotation (13), plus a proved recursion and exact filters that can generate further candidates. For the six-colour target, a next candidate must escape every applicable inherited field screen and then supply a finite graph with an actual equal-port colouring and a proof that all different-port colourings are absent. [Existence/count amplifiers](existence_count_amplifiers_20261011e.md) explains why many equal colourings or a small probability of inequality cannot replace that zero-support obligation. No such tester is produced here.

## 6. Reproduction and finite scope

Implementation: [divisor_equation_recursions_20261011f.py](../../04-computation/experiments/divisor_equation_recursions_20261011f.py). Saved output: [divisor_equation_recursions_20261011f.out](divisor_equation_recursions_20261011f.out).

```text
python -B -X utf8 04-computation/experiments/divisor_equation_recursions_20261011f.py
python -B -O -X utf8 04-computation/experiments/divisor_equation_recursions_20261011f.py
```

The checker compares the complete recursive record list through $100000$ with an independent divisor sieve; checks all cover edges and a four-chain when present for the first thirteen records (133 divisor vertices); checks the rational polynomial identity on the declared $[-4,4]^3$ box; compares sequential and binary elliptic addition through $24P$; verifies the source fractions, positivity and angle square classes; and includes fourteen exact-type/domain hostiles. There are **1,897 always-active checks**, unchanged by optimized mode. Normal, optimized and saved LF output agree at SHA256 `7cd96483d87e713c872d6dc606a13aa45499101c7fb34fe4481c2d17e9af4f5e`. The bounded empty search and the finite multiple census are explicitly separate from the all-height algebraic proofs.
