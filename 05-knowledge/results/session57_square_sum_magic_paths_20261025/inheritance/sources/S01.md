# AP square-sum triangles: a complete Pythagorean classification and Berggren construction

**Status: PROVED, with FINITE-EXACT reproducibility checks and an independent audit.**
The four-square arithmetic-progression obstruction is **CITED**, with the
specific consequences proved below. No novelty claim is made. The elementary
three-square/Pythagorean correspondence is classical; the point here is the
exact graph positivity cone, scale normalization, inherited tree operations,
minimum certificate, counting measure, and local-move interface.

Here “triangle” means three vertices spanning a `K3` in the square-sum graph
`Q_N`: its vertices are `1,...,N`, and distinct vertices are adjacent when
their sum is a square. The labels are not claimed to be lengths of a Euclidean
triangle. In particular, `482+3362<6242` is irrelevant to graph adjacency.

## 1. Inheritance and the concrete bridge

The immediate predecessor is
[period3_graph_square_bridge_20261019.md](period3_graph_square_bridge_20261019.md),
which gave the unique smallest AP square-sum triangle by largest label,
`(482,3362,6242)`, a finite minimum certificate, and an infinite quartic family.
The relevant repository mechanisms are:

- [THM-3334, Berggren parabolic spine](../../01-canon/theorems/THM-3334-berggren-parabolic-spine-gaussian-collision-torsor.md): consecutive Euclid parameters and the repeated outer branch;
- [THM-3357, three-branch Walsh/parent circuit](../../01-canon/theorems/THM-3357-berggren-three-branch-walsh-level-collapse-and-parent-circuit.md): the exact `L,M,R` convention;
- [THM-3756, odd-square ordinals and affine descent](../../01-canon/theorems/THM-3756-odd-square-ordinal-berggren-affine-descent.md): uniqueness of the canonical primitive tree and its scale caveat;
- [graceful_square_exchange_20261019.md](graceful_square_exchange_20261019.md): midpoint endpoint rotations, whose local difference-preservation condition is exactly an AP equation.

The motivating pair-sum roots are

\[
62^2=482+3362,\quad82^2=482+6242,\quad98^2=3362+6242.
\]

Dividing by two gives `31,41,49`; then

\[
40=(49+31)/2,\qquad9=(49-31)/2,\qquad40^2+9^2=41^2.
\]

This is a lossless change of coordinates, rather than a numerical coincidence.

## 2. Exact classification, including all scales

### Theorem 1

Every positive, strictly increasing AP square-sum triangle is uniquely

\[
(x,y,z)=2g^2(C^2-4AB,\ C^2,\ C^2+4AB),                 \tag{1}
\]

where `g` is a positive integer, `(A,B,C)` is a primitive positive
Pythagorean triple with **legs ordered by size**, `A>B`, and

\[
C^2>4AB\quad\Longleftrightarrow\quad A/B>2+\sqrt3.      \tag{2}
\]

Its increasing pair-sum roots are

\[
\sqrt{x+y}=2g(A-B),\quad\sqrt{x+z}=2gC,\quad
\sqrt{y+z}=2g(A+B).                                    \tag{3}
\]

Each pair of distinct vertices has gcd exactly `2g²`; each pair of
pair-sum roots has gcd exactly `2g`. Thus there are no gcd-one triangles of
this kind. “Base triangle” below means `g=1`, with the unavoidable vertex
gcd `2`.

**Proof.** Write `x+z=2y`, and let `p,q,r` be the positive increasing roots
of `x+y,x+z,y+z`. Then

\[
p^2+r^2=2q^2,\qquad q^2=2y.
\]

The second identity makes `q` even. Reduction of the first modulo four
forces `p,r` even too. Put `(u,v,w)=(p,q,r)/2`. Reduction modulo two shows
`u,w` have the same parity, so

\[
P=(w+u)/2,\quad Q=(w-u)/2,\qquad P^2+Q^2=v^2
\]

are positive integers with `P>Q`. Divide this triple by its gcd `g`, giving
`(A,B,C)`. In a primitive Pythagorean triple the legs have opposite parity
and `C` is odd. Solving the three edge-sum equations gives (1) and (3).
Every operation was uniquely reversible.

Conversely, (1) has common difference `8g²AB>0`; (2) gives positivity;
and expanding (3) proves all three graph edges. The quadratic inequality
`(A/B)²-4(A/B)+1>0`, with `A/B>1`, proves the equivalence in (2).

The odd integers `C²-4AB,C²,C²+4AB` are pairwise coprime. A common odd
prime would divide `C` and `AB`, contradicting the pairwise coprimality
of a primitive Pythagorean triple. For the first and last entries, use
that their difference is `8AB`. Likewise, `A-B,C,A+B` are pairwise
coprime and odd. This proves both exact gcd claims. QED.

The even leg of a primitive Pythagorean triple is divisible by four, and
one leg is divisible by three. Consequently

\[
z-y=y-x=8g^2AB\quad\hbox{is divisible by }96g^2.       \tag{4}
\]

All three vertex cores in (1) are odd, so the vertices also have the
identical two-adic valuation `1+2 v_2(g)`, which is odd. In particular
all labels are even.

Scaling a base triangle by an arbitrary integer does not preserve square
sums; the exhaustive scale law is square scaling `g²`. For example twice
`(482,3362,6242)` is not a square-sum triangle.

## 3. The precise Berggren cone and an effective word generator

Keep the repository's parity-oriented convention separate from the
size-oriented `(A,B)` above:

\[
(a,b,c)=(n^2-m^2,2mn,n^2+m^2),\qquad
0<m<n,\quad\gcd(m,n)=1,\quad m+n\text{ odd}.             \tag{5}
\]

Thus `a` is the odd leg and `b` the even leg; either can be larger. Formula
(1) is symmetric in the legs, so apply it with `AB=ab,C=c`.
THM-3357's children are

\[
L(m,n)=(n,2n-m),\quad M(m,n)=(n,2n+m),\quad
R(m,n)=(m,n+2m).                                       \tag{6}
\]

For `t=n/m>1`, their slope maps are

\[
L:t\mapsto2-1/t,\quad M:t\mapsto2+1/t,\quad R:t\mapsto t+2.
\]

Set

\[
\mu=\cot(5\pi/24)=2-\sqrt3+\sqrt6-\sqrt2
 =1.3032253728\ldots,
\]
\[
\lambda=\cot(\pi/24)=2+\sqrt3+\sqrt6+\sqrt2
 =7.5957541127\ldots.
\]

### Theorem 2

A primitive node gives a positive AP square-sum triangle if and only if

\[
1<t<\mu\quad\hbox{or}\quad t>\lambda.                  \tag{7}
\]

In the first cone, exactly its `L` child is admissible; in the second,
exactly its `R` child is admissible. No `M` child of any positive node is
admissible. The induced graph on admissible nodes is therefore a disjoint
union of infinite monochromatic rays, rather than a ternary subtree.
Every node has an admissible `L⁴` descendant and an admissible `R⁴`
descendant.

**Proof.** For `theta=arctan(m/n)` in `(0,pi/4)`,

\[
a/c=\cos2\theta,\qquad b/c=\sin2\theta,\qquad
4ab/c^2=2\sin4\theta.
\]

Positivity is `sin4theta<1/2`, giving
`theta` in `(0,pi/24)` or `(5pi/24,pi/4)`, exactly (7).
The middle slope is in `(2,3)`, outside both cones. If `t<mu`, then
`1<L(t)<t<mu`, whereas `R(t)` is between `3` and `mu+2` and fails.
If `t>lambda`, its `R` child remains there, while `L(t)<2` and is
larger than `mu`, hence fails. Finally,

\[
L^k(t)=1+\frac1{k+1/(t-1)},\qquad R^k(t)=t+2k,          \tag{8}
\]

so `L⁴(t)<5/4<mu` and `R⁴(t)>9>lambda`. QED.

There is no rational slope on either boundary: equality in (2) would
make the rational leg ratio equal to `2+sqrt3` or its reciprocal.
The strict inequalities therefore lose no integral equality case.

### Corollary 2.1: a three-or-four-step entry rule

Every admissible word ends in at least three consecutive `L` letters or
at least three consecutive `R` letters. Every word ending `LLLL` or
`RRRR` is admissible. In each maximal run of one outer letter, first
entry occurs at its third or fourth step, and admissibility persists
until the next change of letter.

Indeed, a maximal `L` run starts at the root slope `2`, an `M` slope in
`(2,3)`, or an `R` slope in `(3,infinity)`. Applying `L` twice leaves
slope at least `4/3>mu`. Applying it four times always succeeds.
A maximal `R` run starts at a slope in `(1,3)` and its second iterate
is below `7<lambda`, while its fourth always succeeds. At the third
step the exact prefix slope still matters:

\[
L^3(t)<\mu\iff t<1+\frac{\mu-1}{4-3\mu},\qquad
R^3(t)>\lambda\iff t>\lambda-6.                        \tag{9}
\]

The two thresholds are approximately `4.3570897074` and `1.5957541127`.
Thus no assertion of a finite-state language based solely on the last
three letters is intended. The Euclid-slope boundaries have degree four
over the rationals; the fact that the **leg-ratio** wall is quadratic does
not license a quadratic-continued-fraction or regular-language claim.

The root path is

\[
(3,4,5)\xrightarrow L(5,12,13)\xrightarrow L(7,24,25)
\xrightarrow L(9,40,41).
\]

The old family with consecutive parameters `(m,n)=(k,k+1)`, `k>=4`, is
exactly this `L` spine: root word `L^(k-1)`. The ternary tree plus (7)
now classifies every base triangle, not only that one infinite family.

## 4. An analytic proof of the unique 6242 minimum

The cone description replaces the predecessor's finite search with a short
unbounded argument. In the low cone, `t<mu<4/3`. Because `n>=m+1`, this
forces `m>=4`. If `m=4`, the only integer `n` with `4<n<4mu` is `5`,
which is primitive and admissible. It gives

\[
(a,b,c)=(9,40,41),\qquad(x,y,z)=(482,3362,6242).
\]

Every other low-cone choice has `m>=5,n>=6`, hence `c>=61` and
`z>2c²>=7442>6242`. In the high cone, `t>lambda>7` forces `n>=8,m>=1`,
hence `c>=65` and `z>8450>6242`. Scales `g>=2` only increase every label
by at least four. This proves the unique minimum by **largest vertex**.

It is not a minimum by smallest vertex: the node `R³` has
`(a,b,c)=(63,16,65)` and produces `(386,8450,16514)`.

## 5. A shared legal move for square-sum and graceful paths

For an AP square-sum triangle `x<y<z`, let `d=y-x=z-y`. The two paths

\[
(y,z,x)\quad\hbox{and}\quad(y,x,z)
\]

both use square-sum edges and have the same edge-difference multiset
`{d,2d}`. Their path centers are respectively `z` and `x`; this does not
contradict the predecessor's observation that the different order
`(x,y,z)` has repeated edge differences `d,d`.

More generally, suppose a simple square-sum path ends at `x`, contains
consecutive vertices `y,z` in that orientation, and its suffix is

\[
\ldots,y,z,\ldots,x.
\]

Delete edge `{y,z}`, add `{y,x}`, and reverse the suffix to obtain
`...,y,x,...,z`. This is a valid Posa endpoint rotation when `{y,x}`
was not already a path edge. The square-sum predicate is preserved because
`x+y` is a square. The edge-difference multiset is preserved because
`|y-z|=|y-x|=d`; all other undirected edges are unchanged.

Consequently the rotation preserves a graceful labeling **if the entire
host path already has one**, with its required vertex-label range and
all edge differences present exactly once. This supplies infinitely many
exact local switch templates. It does not construct a graceful Hamiltonian
path, insert unused vertices, or make the three actual labels a literal
graceful labeling on `{0,1,2}`. Endpoint placement, host edges, label range,
and the unused-edge condition remain load-bearing side data.

## 6. A classical obstruction to gluing, and a construction that survives

Let

\[
p=A-B,\quad q=C,\quad r=A+B,\quad
D=C^2-4AB,\quad E=C^2+4AB.
\]

The sequence

\[
D,\quad p^2,\quad q^2,\quad r^2,\quad E               \tag{10}
\]

is an AP of common difference `2AB`. Fermat's four-square theorem says that
four distinct rational squares cannot be in arithmetic progression; see
[Alf van der Poorten, *Fermat's Four Squares Theorem*, arXiv:0712.3850v1](https://arxiv.org/abs/0712.3850v1),
which gives a descent proof. Thus neither positive `D` nor `E` is a square.
After multiplication by `g²`, that conclusion still holds.

**Vertex-type obstruction.** In every positive AP square-sum triangle,
its middle vertex is twice a square and neither endpoint is twice a square.
Therefore an endpoint of one such triangle cannot coincide with the middle
of another, at any permitted square scale. This is an exact obstruction to
one proposed chain-of-triangles construction.

Gluing at a common middle survives. The smallest hypotenuse with two
admissible primitive representations is **FINITE-EXACT** `C=145`:

\[
(17,144,145)\mapsto(22466,42050,61634),
\]
\[
(143,24,145)\mapsto(14594,42050,69506).
\]

A stronger unbounded construction needs no fixed-hypotenuse multiplicity
theorem. Choose any `k` distinct admissible primitive triples
`(A_i,B_i,C_i)`, and set `H=lcm(C_1,...,C_k)`. Apply (1) with
`g_i=H/C_i`. Every resulting center is `2H²`. Their distances from the
center are

\[
\delta_i=8H^2 A_iB_i/C_i^2.
\]

These are distinct: `AB/C²=r/(1+r²)` is strictly increasing in the
smaller/larger leg ratio `r` in `(0,1)`, and equal primitive leg ratios
would give the same primitive triple. All lower endpoints lie below the
common center, and all upper endpoints above it. Thus all `2k+1` vertices
are distinct and the `k` triangles share only their center.

**Construction theorem.** Some square-sum graph contains a friendship
(windmill) graph of `k` AP triangles for every `k>=1`. This is a subgraph
claim, not an induced-subgraph claim. For example, selecting successive
members of the `L` spine supplies arbitrarily many input triples.

**Sharp host warning.** The bare windmill has no Hamiltonian path for
`k>=3`: deleting its center leaves `k` components, whereas deleting a
vertex of a Hamiltonian path leaves at most two. Extra square-sum cross-edges
or a different host construction are necessary for Hamiltonian use.

## 7. Counts: one-third by hypotenuse, much less by tree depth

Let `A(X)` count admissible primitive triples with `C<=X`. Let `P(N)`
count base AP square-sum triangles with largest label at most `N`, and
`T(N)` count all such triangles including every square scale.

### Theorem 3

\[
A(X)=\frac{X}{6\pi}+O(\sqrt X\log X),                 \tag{11}
\]

\[
P(N)=\kappa\sqrt N+O(N^{1/4}\log N),                  \tag{12}
\]

\[
T(N)=\frac\kappa2\sqrt N\log N+O(\sqrt N),             \tag{13}
\]

where

\[
\kappa=\frac1{\sqrt2\pi^2}
\int_0^{\pi/6}\frac{du}{\sqrt{1+2\sin u}}
=0.0309567443988252174\ldots.                           \tag{14}
\]

**Proof.** Coprime integer pairs of opposite parity have planar density
`4/pi²`. More precisely, for any fixed bounded region with piecewise smooth
boundary and dilation `R`, their count is

\[
\frac4{\pi^2}\operatorname{area}(\Omega)R^2+O(R\log R).
\]

To see the constant and error, opposite-parity points in a dilated region
have count `area/2+O(perimeter+1)`. Apply Mobius inversion over odd common
divisors; the coefficient is

\[
\frac12\sum_{d\text{ odd}}\frac{\mu(d)}{d^2}
=\frac12\frac{1/\zeta(2)}{1-1/4}=\frac4{\pi^2}.
\]

The boundary-error sum is `O(R log R)`; the tail of the convergent
area sum is `O(R)`. The Euclid sector is `0<theta<pi/4`; the two
admissible sectors have total width `pi/12`. Radius squared is `C`,
so their area below `C<=X` is `pi X/24`, proving (11). The full sector
instead gives `X/(2pi)+O(sqrt X log X)`, hence an asymptotic fraction `1/3`.

For the largest-vertex bound, put `m=r sin theta,n=r cos theta`. Then

\[
z=2r^4(1+2\sin4\theta).
\]

The admissible region below `z<=N` is a dilation by `N^(1/4)` of a fixed
piecewise smooth region. Polar integration and the symmetry of the two
sectors give precisely (14), proving (12).

Uniqueness of the scale decomposition gives

\[
T(N)=\sum_{g\le\sqrt{N/6242}}P(N/g^2).
\]

The main term is `kappa sqrt N` times this harmonic sum, yielding (13).
The error is bounded by a constant times

\[
N^{1/4}\sum_{g\le G}g^{-1/2}\log(2+N/g^2)=O(\sqrt N),
\]

with `G=floor(sqrt(N/6242))`; use
`log(2+N/g²)=O(1+log(G/g))` and the corresponding integral bounds.
This proves the stated remainder. QED.

This is a different measure from uniformly sampling a ternary level.
Corollary 2.1 proves for every depth `d>=4` that its admissible fraction
lies between

\[
\frac2{81}\quad\hbox{and}\quad\frac2{27},              \tag{15}
\]

which excludes `1/3`. The exact depth counts for `d=3,...,10` are
`2,5,15,44,132,396,1188,3566` out of `3^d` nodes. No limit for the
ternary-depth fraction is claimed. This is a concrete instance of the
repository's warning that changing the address/height changes the measure.

## 8. Congruent-number interface and the next bounded proof route

The common difference in (1) is

\[
8g^2AB=(4g)^2\,(AB/2).
\]

Thus its squarefree class is the congruent-number class of the associated
right-triangle area. For the first example,

\[
AB/2=180=5\cdot6^2,\qquad y-x=2880=5\cdot24^2.
\]

The same primary paper cited above explicitly relates `(9,40,41)` to
congruent number `5`. Dividing those side lengths by six gives the rational
right triangle `(3/2,20/3,41/6)` of area `5`. On the standard curve
`E_5: Y²=X³-25X`, it gives the rational point

\[
(X,Y)=(25/4,75/8).
\]

The familiar inverse formulas

\[
a=|X^2-n^2|/|Y|,\quad b=|2nX/Y|,\quad
c=|X^2+n^2|/|Y|
\]

recover the sides and area `n` directly by substitution. This is an
explicit elliptic-curve proof route for seeking additional triangles in
a fixed gap-squareclass. The present packet proves the forward invariant
and the displayed example; it does not rely on or claim a complete
real-component-density proof that every congruent-number class enters
the positivity cone. That would require carrying rational group-law,
non-torsion, and real-angle data, rather than merely recognizing a curve.

## 9. Connection contract and stopping boundary

The lossless object is `(primitive Berggren node, positive scale g)`.
The map (1) preserves exactly AP shape and all three square-sum edges;
its inverse is supplied in Theorem 1. Forgetting `g` forgets the graph-label
height. Forgetting the node but keeping only `C` forgets the two endpoint
labels, as the `C=145` example shows. Forgetting the slope in favor of only
word depth loses positivity and the third-step selector in (9).

Concrete transfers now available are:

1. an exhaustive generator and unique inverse, instead of a single family;
2. an analytic minimum proof without a finite search;
3. a one-edge local move preserving both square sums and difference labels;
4. an endpoint/middle gluing obstruction from the four-square theorem;
5. arbitrary common-center friendship subgraphs by square scaling;
6. explicit asymptotics and a rigorous distinction between two sampling measures.

None supplies a Hamiltonian path through every label of `Q_N`, proves the
graceful-tree conjecture, or closes LRC, Collatz, or an elliptic-curve rank
problem. Those would need an independently proved global extension or
coverage mechanism.

## 10. Reproduction and audit scope

From repository root:

```bash
python3 04-computation/experiments/ap_square_pythagorean_20261019.py
python3 -O 04-computation/experiments/ap_square_pythagorean_20261019.py
python3 04-computation/experiments/ap_square_pythagorean_20261019_audit.py
python3 -O 04-computation/experiments/ap_square_pythagorean_20261019_audit.py
```

The main companion checks all `159139` primitive triples with `C<=10^6`,
of which `53052` are admissible; all `88573` tree nodes through depth ten;
and an independent direct graph census of all `79` AP square-sum triangles
with largest vertex at most `10^6`. It checks both gcd laws, inverse
normalization, branch predicates, suffix controls, endpoint types, gap
divisibility, rotation differences, the `C=145` shared center, and windmill
constructions for `k=1,...,8`. There are `1071641` explicit checks, retained
under optimized Python. The max-label counts are:

| N | base triangles P(N) | all scales T(N) |
|---:|---:|---:|
| 10^6 | 31 | 79 |
| 10^8 | 309 | 1498 |
| 10^10 | 3093 | 22125 |
| 10^12 | 30952 | 292514 |

The independent audit separately checks `9841` tree nodes through depth
eight, `15919` primitive triples through `C<=100000` (of which `5305`
are admissible), and the same complete graph census by middle-square
enumeration. Both scripts give byte-identical output with and without
Python optimization.

These tables audit the formulas and conventions; the unbounded theorems
follow from the proofs. Firstness at `C=145` is finite-exact, while the
6242 minimum has the independent analytic proof in Section 4.
