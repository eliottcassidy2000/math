# Congruent numbers are exactly the squarefree gap classes of positive AP square-sum triangles

**Status:** PROVED, with classical real elliptic-group structure explicitly cited for the bounded-multiplier refinement. Exact examples and a finite test universe are recorded separately.

**Date:** 2026-10-19.

This continues the [Pythagorean classification](ap_square_pythagorean_20261019.md). Distinct label triples here are compared up to multiplication, not affine translation: every nonconstant three-term AP has the same affine shape.

**Scope:** A square-sum triangle means a three-vertex clique in the graph where two distinct positive integers are adjacent when their sum is a square. It does **not** mean that these three integers are side lengths of a Euclidean triangle. In particular, `(482,3362,6242)` is a square-sum triangle but does not satisfy the Euclidean triangle inequality.

**Priority:** This is a new derivation in this research session, not a claim of literature novelty. The classical congruent-number / three-squares-in-AP correspondence is distinguished below from the extra positivity constraint handled here.

## 1. Main theorem

For a positive squarefree integer `N`, the following are equivalent:

1. `N` is a congruent number: a rational right triangle has area `N`.
2. There are positive integers `a,d`, with `d>0`, such that the three numbers
   `a, a+d, a+2d` have square pairwise sums and the squarefree part of `d` is `N`.
3. There are infinitely many such triples with greatest common divisor `2`, pairwise nonproportional as triples of integer labels, with squarefree gap class `N`.

More constructively, from **any** rational right triangle of area `N`, form its standard positive-branch point `P` on

`E_N: y² = x³ − N²x`.

At least one of `P,2P,3P` gives a rational right triangle whose larger-leg / smaller-leg ratio is strictly greater than `2+sqrt(3)`. Clearing denominators and taking the primitive integer right triangle `(A,B,C)` then gives the positive AP square-sum triangle

`2(C²−4AB), 2C², 2(C²+4AB)`.

The bound `3` is sharp for a universal assertion involving the first consecutive multiples: the area-30 triangle `(5,12,13)` fails at `P` and `2P` and succeeds at `3P`.

This is an exact reformulation and a bounded conversion from a supplied congruent-number witness. It does not decide which positive squarefree integers are congruent and does not solve the congruent-number problem, BSD, or LRC(14).

## 2. The elementary maps and the positivity coordinate

### 2.1 Necessity needs no elliptic curve

Suppose the pairwise sums of `(a,a+d,a+2d)` are squares:

`r²=2a+d, s²=2a+2d, t²=2a+3d`.

Then `(t+r,t−r,2s)` is an integer right triangle:

`(t+r)²+(t−r)²=2(t²+r²)=4s²`.

Its area is `(t²−r²)/2=d`. If `d=N h²`, divide all three sides by `h` to obtain rational area `N`. Thus every realized squarefree gap class is congruent.

### 2.2 Forward construction from a primitive right triangle

Let `A,B,C` be positive coprime integers with `A²+B²=C²`. Define

`U=2(C²−4AB), V=2C², W=2(C²+4AB)`.

Then

- `U+V=4(A−B)²`;
- `U+W=4C²`;
- `V+W=4(A+B)²`;
- `V−U=W−V=8AB=16(AB/2)`.

Consequently this is a positive AP square-sum triple exactly when

`C²>4AB`, equivalently `max(A/B,B/A)>2+sqrt(3)`.

The other inequality root is `2−sqrt(3)=1/(2+sqrt(3))`. The boundary cannot occur for positive rational legs, because it would force an irrational leg ratio.

Since a primitive Pythagorean triple has coprime legs of opposite parity and odd `C`, `gcd(C²,4AB)=1`; hence `gcd(U,V,W)=2`.

If `(A,B,C)` is obtained by rational scaling `lambda` from an area-`N` rational triangle, its area is `N lambda²`. This integer area and squarefreeness of `N` imply `lambda` is an integer: writing `lambda=p/q` in lowest terms would give `q² | N`. Its AP gap is therefore `N(4 lambda)²`.

### 2.3 Curve and triangle formulas

For a point `(x,y)` on `E_N` with `x>N` and `y != 0`, use positive sides

`a=|(x²−N²)/y|, b=|2Nx/y|, c=|(x²+N²)/y|`.

Direct substitution gives `a²+b²=c²` and `ab/2=N`. Conversely, a positive rational triangle `(a,b,c)` of area `N` gives

`x=N(a+c)/b, y=2Nx/b`.

Here `x>N`, and substitution recovers the original triangle with the designated leg ordering. Put `u=x/N>1`. The ordered leg ratio is

`r=a/b=(u−1/u)/2`.

The scale-free coordinate `u` is therefore sufficient to test the positivity chamber. Equivalently, the exact rational test is

`(u²+1)² > 8u(u²−1)`.

No floating-point comparison with `sqrt(3)` is needed in computation.

## 3. Elementary non-torsion proof using only doubling and v₂

This section avoids invoking the general congruent-curve torsion theorem, Dirichlet's theorem, Nagell–Lutz, or Mordell–Weil.

For a finite rational point `P=(x,y)` with `y != 0`, the usual chord-and-tangent law gives

`x(2P)=(x²+N²)²/(4y²)=(x²+N²)²/[4x(x²−N²)]`.

Write `v=v₂`, `k=v(x)`, and `e=v(N)`, where squarefreeness gives `e in {0,1}`. The valuations below are valid for rational numbers; an odd 2-adic unit, including a rational with odd denominator, has square congruent to `1 mod 8`.

- If `k<0`, the `x²` term has strictly smaller valuation than `N²`, so
  `v(x(2P))=4k−(2+3k)=k−2<0`.
- If `e=0` and `k>0`, the sum and difference `x² +/- N²` have valuation `0`, so `v(x(2P))=−2−k<0`.
- If `e=0` and `k=0`, then `v(x²+N²)=1` and `h=v(x²−N²)>=3`; thus `v(x(2P))=−h<0`.
- If `e=1` and `k=0`, both sum and difference have valuation `0`, so `v(x(2P))=−2`.
- If `e=1` and `k=1`, write `x=2u, N=2v₀`, with `u,v₀` odd units. Then `v(x²+N²)=3` and `v(x²−N²)=2+h`, where `h=v(u²−v₀²)>=3`. Hence `v(x(2P))=1−h<=−2`.
- If `e=1` and `k>=2`, both sum and difference have valuation `2`, so `v(x(2P))=−k<0`.

All cases yield `v₂(x(2P))<0`. Subsequently each doubling reduces this valuation by exactly `2`. Such a point cannot become a rational 2-torsion point, whose x-coordinate is `0,N`, or `−N`; and the coordinates of `2P,4P,8P,...` are distinct. Therefore `P` has infinite order.

This proves, in particular, the only rational torsion points of `E_N` are `O,(0,0),(N,0),(−N,0)`.

## 4. A low-dependency existence proof by compactness

This proof establishes sufficiency and infinitely many shapes without density or any estimate on the first successful multiplier.

The nonsingular projective real curve `E_N(R)` is a closed subset of the compact real projective plane. Its group law and inverse are continuous. For the infinite-order rational point `P`, the distinct points `P,2P,3P,...` have a convergent subsequence `n_j P -> Q` with strictly increasing positive indices.

Then

`R_j=(n_(j+1)−n_j)P -> Q−Q=O`.

Every `R_j` is nonzero by infinite order. A nonzero sequence tending to `O` cannot consist of only finitely many points, so it contains infinitely many distinct values. For this Weierstrass curve, approaching `O` through real affine points forces `x -> +infinity`: real affine x-coordinates lie in `[-N,0] union [N,infinity)`, and the portion with x bounded is compact in the affine chart and excludes `O`.

Consequently

`a/b=(x²−N²)/(2Nx) -> +infinity`.

We obtain infinitely many rational right triangles of area `N` satisfying `C²>4AB` after primitive reduction. Distinct unbounded leg ratios ensure infinitely many distinct shapes, and the construction in Section 2 completes the main equivalence.

There is no claim that the doubling subsequence alone is dense or that every multiple works. The compactness proof uses arbitrary positive integer multiples and differences.

## 5. Stronger conversion: one of the first three multiples works

### 5.1 Precisely cited group input

The identity component `E_N(R)^0` is a circle group. It consists of `O` and the two real affine branches with `x>=N`. There is a group coordinate `t in R/Z`, with `O` at `0`, the unique nontrivial 2-torsion point in this component `T=(N,0)` at `1/2`, and negation `t -> −t`. On `0<t<1/2`, choose orientation so `u(t)=x(t)/N` strictly decreases from infinity to `1`.

This standard real elliptic-group fact is cited, for example, in Noam Elkies's Harvard Math 223 lecture notes (2024), near the concluding discussion of real elliptic curves:

https://people.math.harvard.edu/~elkies/M223.24/index.html

One can construct this coordinate by integrating the real invariant differential around the component and dividing by its full real period. The only uses here are the group law, the half-period location of `T`, and the stated monotonicity. No equidistribution assertion is used.

### 5.2 Locate a guaranteed admissible arc

Let `R=2+sqrt(3)` and `beta=R+sqrt(R²+1)`, the unique `u>1` for which `r=(u−1/u)/2=R`. Let `delta in (0,1/2)` be the circle coordinate of the upper branch point with `u=beta`.

Translation by `T` gives

`u(P+T)=(u(P)+1)/(u(P)−1)`.

This replaces `r` by `1/r`. Thus the admissible locus is exactly the punctured union of the radius-`delta` arcs around `0` and `1/2` in `R/Z`. We now show the strict bound `delta>1/8`.

Normalized duplication is

`u(2P)=(u²+1)²/[4u(u²−1)]=(r+1/r)/2`.

At `u=beta`, this equals `(R+1/R)/2=2`. At a quarter-period `t=1/4`, `2P=T`. Solving `u(2P)=1` gives

`(u²−2u−1)²=0`, hence `u(1/4)=1+sqrt(2)`.

Because `beta>1+sqrt(2)`, monotonicity gives `delta<1/4`. Hence `0<2delta<1/2`, so the x-coordinate comparison at the doubled boundary is on the same monotone half-circle. Since `u(2delta)=2<1+sqrt(2)=u(1/4)`, we have `2delta>1/4`, or `delta>1/8`.

### 5.3 A four-point circular pigeonhole argument

For any `t in R/Z`, put `s=2t`. Among the four circle points

`0,s,2s,3s`,

a cyclic gap has length at most `1/4`. The endpoints differ by `q s` for an integer `q` with `1<=q<=3`, so

`||q s|| <= 1/4`.

It follows that `q t` has circle distance at most `1/8` from `{0,1/2}`. Since `1/8<delta`, `qP` is strictly inside the admissible locus unless it is exactly `O` or `T`. For rational `P` with nonzero y-coordinate this exception is impossible by Section 3.

Every rational right triangle supplies `P` in the identity component by `x>N`; therefore one of `P,2P,3P` succeeds. For an arbitrary rational non-torsion point on the other real component, apply the statement to `2P`, obtaining one of `2P,4P,6P`.

### 5.4 Infinitely many shapes without an unbounded search per seed

For each `n>=1`, apply the three-multiple construction to `nP`, obtaining some `q_n nP` with `q_n in {1,2,3}`. The positive integer indices `q_n n` are unbounded, and infinite order makes these infinitely many distinct points.

A fixed unordered positive right triangle of area `N` has at most four preimages in the identity component: choose one of two leg orderings and one of two signs of y. A fixed primitive integer shape uniquely fixes its positive area-`N` rational scaling. Thus the construction yields infinitely many primitive shapes, not merely infinitely many square dilations. Distinct unordered primitive right triangles give distinct AP triples: the middle vertex recovers `C²` and the gap recovers `AB`, which together recover the unordered legs. Because every output has gcd `2`, different outputs cannot be rational multiples of one another.

### 5.5 A uniform strict positivity margin

The same argument guarantees more than positivity. At circle distance `1/8` from `O` or `T`, duplication gives a quarter-period, and therefore the larger-to-smaller leg ratio is

`kappa=1+sqrt(2)+sqrt(2+2sqrt(2))`, approximately `4.6115817893`.

The selected rational point cannot have circle coordinate an odd multiple of `1/8`, since that would give order-8 torsion. Its distance to `{0,1/2}` is therefore strictly less than `1/8`, and its unordered leg ratio is strictly greater than `kappa`.

Consequently one can choose a successful multiplier among the first three whose output satisfies

`U/V > 3−2sqrt(2)` and `W/U < 5+4sqrt(2)`.

An exact test for this stronger margin, requiring no irrational arithmetic, is

`(C²+2AB)² < 2C⁴`.

This guarantee concerns a suitable multiplier among the first three, not necessarily the first multiplier passing mere positivity. No optimality over rational points is asserted here.

## 6. Sharpness and the fixed class N=5

### 6.1 Uniform bound 3 cannot be replaced by 2

On `E_30`, take `P=(45,225)`, corresponding to `(5,12,13)`.

- `P`: primitive triangle `(5,12,13)`; `13²<4*5*12`.
- `2P=(169/4,−1547/8)`: primitive triangle `(14161,40560,42961)`; `42961²<4*14161*40560`.
- `3P=(2790045/121,−4660327575/1331)`: primitive triangle `(90025452,34597057445,34597174573)`; the positivity inequality holds.

This is a hostile control for a first-two-multiples assertion, not an exception to the theorem.

### 6.2 Class 5 has infinitely many primitive realizations

The starting rational triangle is `(3/2,20/3,41/6)`, with area `5`, and

`P=(25/4,75/8)` on `E_5`.

It gives primitive `(9,40,41)` and the AP square-sum triple `(482,3362,6242)`, with gap

`2880=5*24²`.

The double is not admissible:

`2P=(1681/144,−62279/1728)`

has primitive triangle `(2307361,2420640,3344161)`, whose leg ratio is too close to `1`. Its formal AP triple has negative first term `−22315497060478`, so positivity must not be silently dropped.

The third multiple is

`3P=(127351225/2439844, 1430549626725/3811036328)`.

Its primitive right triangle is

`(124286848883560,642780541615689,654686219104361)`.

The corresponding positive AP square-sum triple is

`(218114746641670847960177053922,`
` 857228090970326756026018436642,`
` 1496341435298982664091859819362)`.

Its gap is

`639113344328655908065841382720 = 5*357522962711112²`.

Using the convention `m>n` in this paragraph, the Euclidean parameters are `(m,n)=(25470245,2439844)`, which are not consecutive. Thus this explicitly leaves the consecutive-parameter / repeated-left-branch spine containing `(9,40,41)`. No claim is made that all points of `E_5(Q)` are generated by this one point; a rank computation is unnecessary for the infinite-family result.

## 7. Exact test universe and controls

Reproduce using only Python's standard library:

`python3 04-computation/experiments/congruent_square_gap_classes_20261019.py`

The checked universe is all coprime opposite-parity Euclidean parameters `2<=m<=100`, `1<=n<m`, with both leg orderings. For each seed the script checks curve membership, rational area, primitive reduction, two successive valuation drops, exact pair-sum squares, AP gaps, gcd `2`, and at least one positive output among the first three multiples.

Result:

- 4,080 ordered seed points;
- 1,908 distinct squarefree congruent-number classes represented;
- first eligible multiplier counts: `1,486` at `1`; `1,690` at `2`; `904` at `3`;
- first multiplier satisfying the stronger `kappa` margin: `1,228` at `1`; `1,530` at `2`; `1,322` at `3`;
- eligible multiples of the displayed class-5 point through `24`: `1,3,4,5,7,8,11,12,15,16,19,20,23,24`;
- all 14 class-5 outputs have distinct primitive shapes;
- the explicit class-30 first-two failure and the class-5 inadmissible double pass as hostile controls.

The [source](../../04-computation/experiments/congruent_square_gap_classes_20261019.py) and [machine output](congruent_square_gap_classes_20261019.json) are included. This finite enumeration is independent evidence for the formulas and boundary bookkeeping, not the proof of the universal claim.

## 8. Connection ledger and limitations

- **Source:** a positive AP square-sum clique, or a rational right triangle of prescribed area.
- **Target:** a congruent-number squareclass and a point on the congruent elliptic curve.
- **Preserved:** right-triangle identity; area up to rational-square scaling; AP gap squareclass; after the admissible-point selection, strict positivity.
- **Lost by passing to squareclass:** primitive shape, size, leg ordering, point sign, and elliptic multiple/address.
- **Necessary sidecar:** the rational curve point (or primitive right triangle) and the strict chamber test `C²>4AB`.
- **Cheapest decisive hostile:** the class-5 double gives the same squareclass but a negative AP vertex; the class-30 seed refutes a uniform two-multiple bound.
- **What is not proved:** a criterion deciding congruence of arbitrary `N`; a height bound for the initial rational witness; efficient global enumeration by increasing AP gap; completeness of a single cyclic subgroup on a positive-rank curve; any general prize conjecture.

### Nonsquare-gap corollary

The common difference `d` of a positive AP square-sum triangle is never a perfect square. Indeed, Section 2.1 gives an integer right triangle of area `d`, whereas Fermat's square-area right-triangle theorem says no rational right triangle has nonzero rational-square area.

This classical input is proved as Theorem 2 in Keith Conrad, *The Congruent Number Problem*, Harvard College Mathematics Review 2.2 (2008), printed pages 60–61 (PDF pages 62–63):

https://people.math.harvard.edu/~knill/various/eulercuboid/HCMR_Fall2008.pdf#page=62

Only the noncongruence of class `1` is used in this corollary. No assertion about classes `2`, `3`, or prime classes is needed.

## 9. Literature context, with exact separation of roles

Alf van der Poorten, *Fermat's Four Squares Theorem*, arXiv:0712.3850 (2007), introduces the classical correspondence between three rational squares in AP and rational right-triangle area, including the `(9,40,41)` area-180 example:

https://arxiv.org/abs/0712.3850
https://arxiv.org/html/0712.3850v1

That correspondence explains the necessary gap squareclass. It does not by itself discharge the strict positive-vertex inequality in our forward construction. Sections 4 and 5 supply that additional argument.

For the classical torsion statement independently corroborated by our Section 3, see Larry Rolen, *Modular Forms Lecture 8: Torsion Points on the Congruent Number Elliptic Curves* (Vanderbilt, 2020):

https://math.vanderbilt.edu/rolenl/ModularFormsLecture8.pdf

For real elliptic circle groups, see the Elkies notes linked in Section 5.1.

A related immediate graph consequence of Fermat's four-squares theorem is that four distinct vertices in arithmetic progression cannot form a square-sum clique: their pairwise sums include five distinct squares in arithmetic progression. The present theorem therefore describes exactly which squarefree gaps occur at the maximal AP-clique size of three. This uses the classical four-squares obstruction; it is not a new proof of Fermat's theorem.
