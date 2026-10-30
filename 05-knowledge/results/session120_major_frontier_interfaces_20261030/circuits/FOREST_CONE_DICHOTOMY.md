# Forest-cone obstruction to minimal tight-prefix presentations

**Status: proved bounded-family rewrite and exact affine-state optimum.**
This extends the Session118 exclusive common-factor rewrite to every forest
interaction graph, including paths of unbounded length. It is not a theorem
that arbitrary GP/P2 sequences expose such a graph, nor a new unrestricted
Boolean-circuit lower bound. The algebraic ingredients are elementary.

## Model and boundary

Let F=(V,E) be a simple forest with no isolated vertices, e=|E|>0,
v=|V|, and maximum matching size nu. Let delta be 0 or 1. Consider an exclusive
cone computing

    Q = XOR_{ij in E} (x_i AND x_j) XOR (z if delta=1).

There is one distinct AND per edge, combined by a binary XOR tree, with the
optional independent input z as an additional leaf. Its size is 2e+delta-1.
All nonroot product outputs and nonroot XOR outputs have no outside consumers. The
root may have arbitrarily many consumers. The boundary input wires may have
arbitrary outside consumers for the ordinary-coordinate rewrite below.
The root is the only output of this cone; all removed gates are distinct.

## Theorem 1: matching-sensitive, fanout-safe rewrite

There is a fixed-input-coordinate replacement of size

    e + nu + delta - 1,

so the displayed saving is e-nu. The replacement retains all boundary wires
and all outside consumers. In particular, a gate-minimum DAG cannot contain
such an exclusive monomial forest cone unless F is a matching.

**Proof.** Choose a leaf x and its neighbor y. Write R=N(y) minus {x}. Then

    Q_F = y (x + XOR_{u in R} u) + Q_{F-{x,y}}.

Creating the parenthesized form costs deg(y)-1 XOR gates. Create one AND for
its product with y; recurse on F-{x,y}, ignoring isolated vertices. A maximum
matching can always be chosen to contain xy: if y is matched elsewhere, swap
that edge for xy; if y is unmatched, add xy. Consequently the recursion takes
exactly nu pivot steps. Every original edge is removed at exactly one step,
and the sum of degrees of pivot neighbors y is e. Thus there are e-nu
input-transform XORs, nu ANDs, and nu+delta-1 output XORs, proving the count.

The expressions depend only on retained predecessor boundary values; redirect
all uses of the old root to the new root. Exclusivity ensures that no deleted
value is still required. No direct-sum lower bound is being assumed. If F is
not a matching, its edges cannot all belong to a matching, so e>nu and this is
a strict reduction, contradicting gate minimality.

**Sharp subfamily.** For a disjoint union of nontrivial stars, nu equals the
number c of components and e=v-c. The replacement size is v+delta-1. When the boundary values are independent formal inputs, every
boundary input is essential, so the usual essential-input bound makes this
ordinary-coordinate cone implementation optimal. This lower-bound assertion
is not made for dependent internal boundary signals. For general forests the formula
is an upper bound, not a claim of exact original-coordinate complexity.

## Theorem 2: owned-input affine normalization

Strengthen the boundary assumption: every forest vertex input is a free,
unpaired parameter with no computational consumers outside this cone
(including affine equations). The optional z is also an independent free
parameter. Then semantic rewriting followed by legal M/D4 preparation gives

    Q = p_1 q_1 + ... + p_nu q_nu + (z if delta=1)

with exactly 2nu+delta-1 charged gates, no new affine equations remaining, and
v-2nu unused parameters. The affine embedding retains ALL parameters, including
unused ones, and remains invertible. The original domain dimension and safe
core dimension do not decrease. The exact changes are

    Delta s = -2(e-nu),  Delta I = -(v-2nu),
    Delta q = 0,        Delta a = 0.

For W=s+4I+2q, the additional loss is

    2(e-nu) + 4(v-2nu).

These are extra losses for this rewrite, not extra losses already implicit
in a GP/P2 row count.

**Proof and legal realization.** At a leaf pivot replace the incident star
by y times w, with w=x+u_1+...+u_k. The leaf x has no other remaining consumer.
Compute the sum as a chain. Move its first XOR to the affine system (M), then
exchange its output with x and delete the equation (D4). The newly introduced auxiliary block is [1] both before and after exchange;
any pre-existing disjoint affine system retains its invertible block. The fresh parameter
is now the sole input to the next chain XOR apart from u_2, so repeat. Each
chain step removes one charged equation and records the inverse substitution
x=w+u_j in the affine embedding. Each shear is invertible. After the chain,
w occurs only in its new product with y; recurse on the remaining forest.
Earlier product parameters are independent new coordinates, so they do not
retain hidden uses of later leaf parameters. Ownership is essential here:
with an outside consumer of old x, D4 is unavailable.

The number of removed transform equations is e-nu. Combined with Theorem 1's
e-nu ordinary savings this gives the stated s change. Exactly the 2nu pivot
coordinates (and optional z) are essential at the endpoint; every other
original dimension remains as an unused coordinate of the recorded embedding.
All parameters are unpaired throughout, so I loses v-2nu and a is preserved.
No dense affine transformation is being declared free in ordinary circuits.

**Exact optimum across affine coordinates.** The polar matrix of Q has rank
2nu: each invertible leaf shear splits off one nonsingular 2-by-2 block. A
function depending on k coordinates has polar rank at most k, so delta=0
requires at least 2nu essential coordinates in ANY invertible affine chart.
When delta=1, the derivative along the radical vector z is the constant 1.
If only 2nu coordinates were essential, full polar rank on these coordinates
would force every radical vector to have zero projection onto them, and hence
zero derivative, a contradiction. At least 2nu+1 coordinates are necessary.
The essential-input bound therefore gives 2nu+delta-1 gates in every affine
chart, attained by the displayed normal form. This counts computational gates
in the new chart, NOT the original-input cost of implementing that chart.

## Consequences and precise limits

The Session118 residual ab+yb+ac+z is the four-vertex path y-b-a-c with
(e,nu,v,delta)=(3,2,4,1). The formulas recover its 6 -> 5 ordinary-gate rewrite,
then 5 -> 4 affine-state normalization. The theorem also covers every longer
forest and every disconnected forest, with shared root fanout preserved.

For a gate-minimum state, an exclusive monomial quadratic cone therefore has
this bounded structural alternative: its interaction graph is a matching,
or it contains a cycle, or the claimed exclusive monomial-cone description
fails. This is a necessary condition only. A cycle does not imply hardness;
for example a four-cycle factors as (a+c)(b+d).

Initial circuit minimality need not survive restrictions. Thus the theorem
cannot automatically be applied to every later state by assuming it minimum.
Nor does an arbitrary quadratic DAG have an exclusive monomial expansion
without blowup. Marked-pair coordinates cannot be freely mixed by these shears.
For an explicit sharing boundary, let bx and by feed an XOR root but let bx
also have an outside user. Replacing the root by b(x+y) needs two new gates
and retaining bx needs one: three gates, equal to the original three. The
advertised strict one-gate saving is no longer justified by this rewrite.
Likewise, if a marked pair (x,y) excludes 11, the shear x'=x+z turns its
constraint into NOT(x'+z=1 AND y=1), generally not a marked-coordinate pair.
This is why the free/unpaired ownership hypothesis is substantive.

No universal amortized surplus, low-d disperser theorem, candidate coverage
certification, or leading-coefficient improvement is established here.

## Verification

`check_forest_cones.py` independently enumerates all simple forests without
isolated vertices on 2 through 6 labeled vertices. It brute-forces maximum
matchings, computes binary polar ranks, constructs a literal fixed-coordinate
replacement DAG, tests every assignment with delta=0 and delta=1, and checks
forward/inverse shears and exact gate inventories. Star forests additionally
check attainment of the essential-input lower bound. Explicit exceptions,
not Python assertions, implement every guard; normal and -O runs must agree.

## Cycle controls: the forest hypothesis cannot be dropped

Two owned, unpaired cyclic examples exhibit different failures of a naive
extension of the formula 2nu+delta-1 using the original graph's matching size.

* **Triangle:** q=xy+yz+zx has matching number one and polar rank two, but
  its radical contains (1,1,1), along which the derivative is the constant
  one. Under the invertible shear x'=x+z, y'=y+z, z'=z,
  q=x'y'+z'. Its optimum across affine charts is TWO gates, not
  2nu-1=1: the radical-linear flag has appeared although the original
  displayed polynomial had no linear monomial. Three essential coordinates
  are necessary by the rank/radical argument above.
* **Four-cycle:** q=x1x2+x2x3+x3x4+x4x1 equals
  (x1+x3)(x2+x4). Its matching number is two, but its polar rank is only two.
  Taking p=x1+x3, r=x2+x4 and retaining x3,x4 gives an invertible affine
  chart in which q=pr, with ONE gate, not 2nu-1=3. Its rank forces at
  least two essential coordinates, so this is the exact affine optimum.

Thus odd cycles can introduce a nonzero radical-linear restriction, and even
cycles can invalidate the matching-rank equality. These examples are specific
controls, not assertions that all odd or all even cycles behave identically.
The checker verifies both polar ranks, the triangle radical derivative, all
assignments of both identities, and invertibility of the displayed charts.
