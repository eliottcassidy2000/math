# Square-sum bridges between Pythagorean AP petals

**Status: PROVED**, with **FINITE-EXACT** examples and a complete signed-contact
census through primitive hypotenuse 10000. The elliptic change of coordinates,
non-torsion proof, real recurrence, primitive-normalization step, and gluing
argument received an independent adversarial audit. No novelty claim is made.

The new conclusion is a port-preserving extension theorem: every admissible
primitive Pythagorean petal has infinitely many admissible primitive neighbors
through either mixed-sign port. It follows that arbitrarily long distinct
petal chains can be assembled and then scaled to positive integer square-sum
fans whose consecutive outer pairs are arithmetic-progressions with their
common center. This is a theorem about these sparse, AP-paired subgraphs, not
Hamiltonicity of the entire interval graph `Q_N`.

## 1. Inheritance, target, and the coordinate that must be retained

The predecessor is
[AP square-sum triangles](ap_square_pythagorean_20261019.md), Theorem 1 and
Section 9. Its proved map sends an admissible primitive positive Pythagorean
triple `T=(a,b,c)` and integer scale `g>=1` to

\[
2g^2(c^2-4ab,\ c^2,\ c^2+4ab),\qquad c^2>4ab.
\]

Leg order does not affect a petal. Use the odd-leg/even-leg convention when
listing primitive triples, and regard swapped legs as the same petal.
The earlier common-center construction yields a friendship subgraph. Its
bare windmill has a cut vertex and has no Hamiltonian cycle with two or more
petals. That obstruction does not say the actual square-sum graph has no
additional edges between petals.

The useful object is the **signed contact graph**, retaining which endpoint
of each AP triangle participates in a contact. Unsigned connectedness loses
exactly the information needed to concatenate these triangles.

The live concept board is:

- AP triangles: exact primitive node plus square scale; common-center alignment.
- Signed contacts: endpoint signs and rational square predicates; chain gluing.
- Elliptic curves: rational group law and recurrence; universal new neighbors.
- Distinct differences: finite forbidden spreads; no interval-label coverage.

The anchor lane is the AP-paired gluing theorem. The niche lane is the complete
finite signed graph, including a hostile unsigned component. The wildcard
lane is the elliptic curve, which supplies the unbounded extension mechanism.

## 2. Signed contact graph and exact scale transfer

For an admissible primitive triple define its normalized spread

\[
 f_T=\frac{2ab}{c^2},\qquad 0<f_T<\frac12.
\]

For a finite selection of petals, choose a positive common multiple `H` of
their hypotenuses. Put

\[
 M=2H^2,\qquad
 v_T^\varepsilon=2H^2(1+2\varepsilon f_T)
 =2(H/c)^2(c^2+4\varepsilon ab),\quad\varepsilon\in\{-1,+1\}.
 \tag{1}
\]

All vertices are positive integers. Each triple `(v_T^-,M,v_T^+)` is an AP
square-sum triangle. Its edge-sum roots are

\[
 2(H/c)|a-b|,\quad 2H,\quad 2(H/c)(a+b).
 \tag{2}
\]

For distinct primitive petals `T,U`, put a contact between signed ports
`(T,epsilon)` and `(U,eta)` exactly when

\[
 1+\varepsilon f_T+\eta f_U\in\mathbb Q^2.
 \tag{3}
\]

This is an **if and only if** for a square-sum edge between their integer
ports, since their sum is `4H²` times (3). No additional denominator-clearing
factor is needed: an integer that is the square of a rational number is the
square of an integer. Equivalently, `H²` times (3) is already integral and
`H` times its rational square root is integral. One must not simply assume a
product of hypotenuses divides their least common multiple.

The spread determines the primitive petal uniquely up to swapping legs:
`(a+b)/c=sqrt(1+f_T)` and `|a-b|/c=sqrt(1-f_T)` determine the unordered
normalized leg pair. Thus distinct petals have distinct spreads. Formula
(1) then gives distinct lower ports, distinct upper ports, and no equality
between a lower port, the center, and an upper port.

### Theorem 1. Port-compatible chain transfer

Suppose distinct petals `T_1,...,T_k` have contacts between successive petals
such that the two contacts incident to every interior petal use opposite
ports. Then, at `H=lcm(c_1,...,c_k)`, their `2k+1` positive integer vertices
have a Hamiltonian cycle which traverses every petal's two outer endpoints
consecutively. In particular, a chain of mixed-sign contacts

\[
 (T_i,+)\sim(T_{i+1},-),\qquad 1\le i<k,                 \tag{4}
\]

has the explicit cycle

\[
 M,v_{T_1}^-,v_{T_1}^+,v_{T_2}^-,v_{T_2}^+,\ldots,
 v_{T_k}^-,v_{T_k}^+,M.                                \tag{5}
\]

**Proof.** Use each petal's internal edge, the specified contacts, and the two
end edges to the center. Equations (1)--(3) prove integrality and all square
sums; distinct spreads prove no vertex repetition. QED.

Conversely, within the class of Hamiltonian cycles which use every petal's
outer-pair edge, deleting the center and contracting those pair edges
produces precisely a port-compatible simple chain. This converse is not
claimed for cycles which do not respect all pair edges.

## 3. The exact elliptic controller

Write a rational Pythagorean parameter as `t=n/m>1` and set

\[
 f(t)=\frac{4t(t^2-1)}{(t^2+1)^2}.
 \tag{6}
\]

It equals the spread of the normalized triple
`(n²-m²,2mn,n²+m²)`. Admissibility is exactly `0<f(t)<1/2`.

Fix an admissible anchor `T=(a,b,c)` with spread `f_0`, fix a port sign
`epsilon`, and put

\[
 k=\frac{|a+\varepsilon b|}{c}>0,
 \qquad k^2=1+\varepsilon f_0.
 \tag{7}
\]

The opposite-port neighbor equation is

\[
 r^2=1+\varepsilon f_0-\varepsilon f(t).
 \tag{8}
\]

With `u=epsilon*t` and `z=r(u²+1)`, this becomes

\[
 \mathcal C_k:\quad
 z^2=k^2(u^2+1)^2-4u(u^2-1).                            \tag{9}
\]

The anchor itself supplies the rational point
`u_0=epsilon*t_0`, `z_0=t_0²+1`, where `f(t_0)=f_0`. It is a starting point,
not a permitted distinct neighbor.

### Lemma 2. Explicit birational map

For `k>0`, `k!=1`, the smooth completion of (9) is isomorphic over `Q` to

\[
 E_k:\quad Y^2=(X-k^2)(X^2-1).                          \tag{10}
\]

On the indicated affine charts the maps are

\[
 X=\frac{kz+k^2+2u+k^2u^2}{2u^2},\qquad
 Y=\frac{(X-k^2)(Xu-1)}k,                              \tag{11}
\]

\[
 u=\frac{kY+X-k^2}{X(X-k^2)},\qquad
 z=\frac{2u^2X-k^2-2u-k^2u^2}{k}.                      \tag{12}
\]

**Proof.** Substitute directly in both directions. Alternatively set
`L=(z+k+2u/k)/u²`; (9) reduces to

\[
 (L^2-k^2)u^2+(4-4L/k)u+(4/k^2-2kL-2k^2)=0.
\]

Taking `X=(kL+k²)/2` yields (10)--(12). The quartic discriminant is
`2^14(k^4-1)^2`, and the cubic has distinct roots `k²,1,-1`, so both smooth
completions are nonsingular. The rational inverse maps extend across the
finite omitted sets of these smooth projective curves. QED.

At the anchor, the formulas lie in ordinary inverse charts. In particular
`u_0!=0`; a finite quartic point with `X=0` has
`u=-1/(2k²)`, whose magnitude is less than 1 because `k²>1/2`, whereas
`|u_0|>1`. A finite point with `X=k²` in the displayed chart `u!=0` is impossible
when `k!=1`.

### Lemma 3. A universal non-torsion point

If `k=A/B` with positive odd integers `A,B` and `k!=1`, then the point
`P=(0,k)` on (10) has infinite order.

**Proof.** Scale `(X,Y)` by `(B²,B³)` to obtain the integral equation

\[
 \widetilde E:\quad y^2=(x-A^2)(x^2-B^4),\qquad
 \widetilde P=(0,AB^2).
 \tag{13}
\]

The tangent formula gives

\[
 x(2\widetilde P)=A^2+\frac{B^4}{4A^2},
 \qquad v_2(x(2\widetilde P))=-2.                     \tag{14}
\]

More generally, if a rational point `Q` on (13) satisfies
`v_2(x(Q))=-2s<0`, the leading term in the cubic strictly dominates the
other three terms, so `v_2(y(Q))=-3s`. The doubling slope

\[
 \ell=\frac{3x^2-2A^2x-B^4}{2y}
\]

has valuation `-s-1`: its numerator has a unique least-valuation term `3x²`.
In `x(2Q)=ell²+A²-2x`, the first term then has the unique least valuation
`-2s-2`. Hence

\[
 v_2(x(2^j\widetilde P))=-2j\quad(j\ge1).
 \tag{15}
\]

All these points exist and are distinct. Thus `P` is non-torsion. QED.

For every primitive anchor, (7) is already reduced with odd numerator and
odd denominator: a common odd prime dividing `a+epsilon*b` and `c` would
contradict primitive coprimality. Consequently Lemma 3 applies to either
port of every admissible primitive petal. No rank computation or inference
from the mere existence of a nonzero point is being used.

## 4. Universal neighbors and arbitrarily long AP-paired fans

### Lemma 4. Rational recurrence near a rational point

On a real elliptic curve containing a rational non-torsion point `P`, every
rational point `Q` is an accumulation point of other rational points from
either side in a real local coordinate.

**Proof.** The projective real curve is a compact topological group. Among
the distinct multiples `nP`, choose a convergent subsequence. Differences
of successive terms tend to the identity and are nonzero. A further
subsequence consists of distinct differences: otherwise a finite set of
nonidentity points would converge to the identity. Translating by `Q`
produces distinct rational points tending to `Q`. To get both sides, also
use the negatives of these small increments. In a real local group
coordinate, inversion is negation; translation is a local diffeomorphism.
After passage to a subsequence, the two choices approach `Q` from opposite
sides. QED.

Only the standard one-dimensional real elliptic Lie-group structure is imported;
see [Noam Elkies, Math 223, September 18 lecture and real-group discussion](https://people.math.harvard.edu/~elkies/M223.24/index.html).
The compact-group recurrence argument above is included in full. Stronger
claims of equidistribution are unnecessary.

### Theorem 5. Every admissible primitive petal has infinitely many opposite-port neighbors

For every admissible primitive petal `T`, either `epsilon in {-1,+1}`, and
every finite set of forbidden petals, infinitely many other admissible
primitive petals `U` outside that set satisfy

\[
 1+\varepsilon f_T-\varepsilon f_U\in\mathbb Q^2.
 \tag{16}
\]

Moreover, for every positive real `h`, infinitely many such neighbors have
`f_T-h<f_U<f_T`, and infinitely many have `f_T<f_U<f_T+h`.

**Proof.** Apply Lemma 3 to (10), then Lemma 4 at the image of the anchor
point of (9). Lemma 2 brings rational points arbitrarily near that anchor
back to the quartic. The inequalities `t>1` and `0<f(t)<1/2` are open and
strict, so a sufficiently small neighborhood consists of admissible
parameters. Equation (9) gives (16). Since `z_0!=0`, the implicit function
theorem makes `t` a real local coordinate on the quartic at the anchor. Also

\[
 f'(t)=-\frac{4(t^4-6t^2+1)}{(t^2+1)^3}.
\]

Its only zero for `t>1` is `1+sqrt(2)`, where `f(t)=1`, outside the
admissible region. Thus `f` is strictly monotone on a sufficiently small
anchor neighborhood. The two-sided form of Lemma 4 supplies infinitely many
points with spread on each side of `f_0`, as close as desired.

Reduce each resulting rational parameter `t=n/m>1` to lowest terms and
normalize

\[
 (n^2-m^2,2nm,n^2+m^2)/\delta,\qquad
 \delta=\begin{cases}1,&m+n\text{ odd},\\2,&m,n\text{ odd}.\end{cases}
 \tag{17}
\]

This is a primitive PPT, with the legs reordered to the chosen parity
convention if necessary. The rational neighborhood need not preserve
opposite parity; the division by 2 in (17) is essential.

For each fixed spread `s`, equation `f(t)=s` is a nonzero polynomial equation
of degree at most four. Thus it has finitely many parameters, and each
primitive petal accounts for finitely many rational points of (9). Removing
the finite fibers of the forbidden petals still leaves infinitely many
admissible primitive neighbors. QED.

### Corollary 5.1. Unbounded AP-paired fan extension

Every finite mixed-sign chain (4) extends to arbitrarily long chains of
distinct admissible primitive petals. In particular, for every `k>=1`
there are positive integer labels on a square-sum fan `K_1 join P_{2k}`
such that each consecutive outer pair, together with the center, forms an
AP square-sum triangle. The fan has the Hamiltonian cycle (5).

**Proof.** Apply Theorem 5 repeatedly with the current last petal's plus
port and all previous petals forbidden. Apply Theorem 1 only after the
finite chain has been chosen. QED.

The center and integer scale depend on the chosen finite chain and grow
with it. There cannot be infinitely many distinct positive integer petals
at a single fixed center. The infinite object in Theorem 5 is the rational
shape/contact graph, not a fixed finite integer graph.

This is also an effective, though not efficiency-bounded, construction:
enumerate integer multiples of the explicit point `P`, translate by the
anchor point, apply (12), and retain admissible parameters outside the
finite forbidden fibers. Lemma 4 proves that this search eventually finds
an extension. No bound on the first successful multiple is claimed.

### Corollary 5.2. Increasing gaps and the separate graceful rank labeling

For every `k>=2`, the chain in Corollary 5.1 can be chosen with

\[
 0<f_1<f_2<\cdots<f_k<2f_1.                            \tag{18}
\]

After the common integer scaling, its outer Hamiltonian path from the
center has strictly increasing absolute edge differences. The closing edge
has a different difference from every path edge. Rank-compressing the
vertices gives a graceful labeling of that path, explicitly

\[
 k,k-1,k+1,k-2,k+2,\ldots,0,2k.                         \tag{19}
\]

**Proof.** Fix `B=min(2f_1,1/2)`. At each step the two-sided strengthening
of Theorem 5 supplies a fresh mixed-sign neighbor with spread in
`(f_i,B)`. For any desired finite length this constructs (18).
Put `d_i=v_i^+-M=M-v_i^-=4H²f_i`. The path's edge-difference list is

\[
 d_1,2d_1,d_1+d_2,2d_2,\ldots,d_{k-1}+d_k,2d_k.        \tag{20}
\]

It is strictly increasing because the `d_i` are positive and strictly
increasing. The closing difference is `d_k`, and `d_1<d_k<2d_1`, so it
coincides with no path difference when `k>=2`.

There are `k` lower ports and `k` upper ports. Increasing spreads put the
center at rank `k`, lower port `i` at rank `k-i`, and upper port `i` at
rank `k+i`. Thus the path's rank sequence is (19), the reverse-complement of the
classical Walecki order, whose successive absolute differences are exactly
`1,2,...,2k`. QED.

For the classical construction and its complement symmetry, see
[Hicks, Ollis, and Schmitt, Example 4.2 and its preceding paragraph](https://community.middlebury.edu/~jschmitt/papers/HicksOllisSchmitt2018.pdf#page=12).
The earlier [graceful exchange packet](graceful_square_exchange_20261019.md)
also treats that order. Its graceful property is not a new existence result.

There are **two distinct labelings**. The original positive integer labels
preserve square sums and the simultaneous AP petals. The ranks `0,...,2k`
make the path graceful, but rank compression does not preserve square sums.
The cycle is not claimed graceful: its closing rank difference is `k`,
which repeats a path difference. No bounded-label theorem for the original
square-sum realization follows.

## 5. Exact small certificates

### Five vertices

Take `T=(21,220,221)` and `U=(1073,264,1105)`. Their `(+,-)` contact has
rational root `941/1105`. With `H=1105`, a cycle is

\[
2442050,1518050,3366050,175874,4708226,2442050.
\]

The five successive pair-sum roots are

\[
1990,2210,1882,2210,2674.
\]

The AP petals are `(1518050,2442050,3366050)` and
`(175874,2442050,4708226)`. In particular
`3366050+175874=1882²` is a genuine cross-petal bridge.

For an independent elliptic control, (11) applied to
`k=241/221`, `u=33/4`, `z=941/16`, then translated by `x=X-k²` and scaled
by `(221²,221³)`, gives

\[
 (x,y)=\left(\frac{345856}{1089},\frac{20503467712}{35937}\right)
\]

on `y²=x(x+9240)(x+106922)`.

### Seven vertices

The three petals

\[
 (1443,76,1445),\quad(1419,8260,8381),\quad(495,4888,4913)
\]

have contacts `(-,-)` and `(+,-)` in that order, using opposite ports at
the middle petal. Here `H=712385`, and the cycle is

\[
\begin{split}
1014984776450,&\ 1228222358306,\ 801747194594,\\
337514444450,&\ 1692455108450,\ 608014784450,\\
1421954768450,&\ 1014984776450.
\end{split}
\]

The seven successive pair-sum roots are
`1497734,1424770,1067362,1424770,1516730,1424770,1561070`.

## 6. Complete bounded census and hostile port control

The companion script enumerates every primitive parity-oriented positive
PPT with `c<=10000`, keeps exactly `c²>4ab`, and tests all four sign pairs
between every two distinct petals. It finds:

- 528 admissible primitive petals;
- 15 signed cross-petal contacts;
- 503 isolated vertices, six components of size two, three of size three,
  and one of size four;
- a maximum port-compatible simple chain of three petals in this bounded
  universe, witnessed by the seven-vertex example.

The size-four component has hub `(1073,264,1105)` and leaves
`(21,220,221)`, `(3477,236,3485)`, `(9393,776,9425)`. All three hub contacts
use its minus port. Consequently no two-edge path through that hub is
port-compatible. Unsigned connectedness does not supply the required chain.
The unbounded theorem proves that the apparent bounded shortage is a height
cutoff, not a genuine terminal obstruction.

Two independent predicates verify the contact list: integer-square testing
of the numerator over `(c_i*c_j)²`, and reduced-fraction numerator/denominator
square testing. The output retains all 15 contacts and their rational roots.
The program also checks both signs of the elliptic map and inverse on all
528 anchors, the first four valuation-controlled doublings for every sign,
and every integer edge in the five- and seven-vertex certificates. Separate
controls check the increasing-gap/rank formula for `k=2,...,50` and the
odd/odd rational-parameter normalization.

From the repository root, after placing the script in the experiments folder:

```bash
python3 04-computation/experiments/square_sum_ppt_bridges_20261019.py
python3 -O 04-computation/experiments/square_sum_ppt_bridges_20261019.py
```

Both executions give byte-identical output and retain all 17153 explicit
checks. The finite census supports the implementation and examples; the
unbounded results follow from the proofs, not extrapolation from that census.

## 7. Connection contract and remaining boundary

Source: signed, admissible primitive-PPT contact chains. Target: positive
integer AP-paired square-sum fans. Map: normalized spread, common-hypotenuse
LCM scale, and formulas (1)--(5). Preserved predicate: every required
square-sum edge, AP pairing, distinct vertex labels, and compatible port
traversal. Lost data: any height bound, density among interval labels, and
coverage of unused vertices. Required sidecar: endpoint signs; dropping them
fails on the explicit 1105-centered component.

The elliptic correspondence is a separate exact bridge from the quartic
contact equation to a rational group law. Its sidecars are primitive
normalization, strict positivity, and finite forbidden fibers. The valuation
argument supplies infinite order; real compactness supplies recurrence into
the correct positive cone. Neither ingredient alone would complete gluing.

Ordinary square-sum fan existence is easier by rational conic constructions
and common square scaling. Indeed, for existing rational labels `M,u`, choose
nonzero rational `s`, set `p=(s+(M-u)/s)/2`, `q=(s-(M-u)/s)/2`, and add
`v=p²-M`. Then `M+v=p²` and `u+v=q²`; taking `|s|` large ensures positivity,
and finitely many parameters can be avoided to keep labels distinct. A finite
construction can be cleared to integers by common square scaling.
The substantive additional constraint here is the
simultaneous AP pairing with a single center, expressed by the primitive-PPT
and port-preserving extension theorem. This note does not claim that bare
fan existence is new, nor that these sparse examples prove a Hamiltonian
path through every label of `Q_N`, the graceful-tree conjecture, or another
unrelated global problem.
