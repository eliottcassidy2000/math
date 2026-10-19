# Every AP square-sum triangle extends to arbitrarily many common neighbors after finite square scaling

**Status: PROVED**, using the previously proved non-torsion mechanism and the
standard continuous real elliptic group law. The four-anchor obstruction is
**CITED** Fermat; the displayed certificates and bounded test universe are
**FINITE-EXACT**. No literature-priority claim is made.

**Date:** 2026-10-19.

## 1. Exact statement and quantifiers

Let `Q` be the simple graph on positive integers with an edge between distinct
labels `u,v` exactly when `u+v` is a perfect square. Rational realizations use
positive rational labels and rational squares instead.

### Theorem 1. Common-neighbor extension of an AP triangle

For every positive integer AP square-sum triangle

`A={a,a+d,a+2d}`, with `a,d>0`,

and every positive integer `k`, there are a positive integer `L` and distinct
positive integers `b_1,...,b_k`, disjoint from `L² A`, such that every `b_i`
is adjacent to all three vertices of `L² A`.

Consequently `Q` contains a copy of `K_3 join I_k` whose distinguished triangle
is the square-dilated original AP triangle. Here this is a **subgraph**, not
an assertion that all leaf-pair sums are nonsquares. Additional edges may be
discarded in taking a subgraph.

The quantifiers are `for every A, for every k, there exists L`. The one common
square factor is chosen **after a finite rational family** has been selected.
There is no assertion of one fixed integer dilation serving every `k`.

### Theorem 2. Sharp AP-anchor threshold

For arbitrary positive rational `a,d`, the following are equivalent:

1. There is a rational common neighbor, distinct from the anchors, of
   `a,a+d,a+2d`.
2. There are infinitely many positive rational common neighbors, unbounded
   above, and avoiding any prescribed finite set.
3. A rational right triangle has area `d`, equivalently the positive rational
   squareclass of `d` is a congruent-number class.

The anchors need not themselves form a triangle. In contrast, four distinct
rational anchors `a,a+d,a+2d,a+3d` with `d!=0` have **no rational common
neighbor at all**, even if positivity and graph distinctness are dropped.

This is the sharp threshold for **anchors in arithmetic progression**. It
is not a bound for arbitrary anchor sets.

## 2. Inheritance and the new move

The closest proved mechanism is
[Congruent square-gap classes](congruent_square_gap_classes_20261019.md),
Sections 2--4: an AP square-sum triangle gives a rational right triangle of
area equal to its AP gap; every nonzero-y rational point of a congruent curve
has infinite order; compactness supplies rational points approaching infinity.

The companion
[Pythagorean port-bridge theorem](square_sum_ppt_bridges_20261019.md) constructs
arbitrarily long AP-paired fans by moving through different primitive petals
and preserving endpoint signs. The present operation instead fixes **all
three rational anchors** and adjoins vertices adjacent to the entire triple.
It leads directly to `K_4` and `K_{3,3}`, without constructing a fan chain.

The live concept board is: AP anchor constraints; congruent elliptic curves;
finite common square scaling; complete bipartite subgraphs; and Fermat's
four-squares obstruction. The anchor lane is the extension theorem; the niche
lane is the complete fixed-integer hostile; the wildcard is the exact
three-versus-four AP-anchor boundary.

The essential retained sidecar is the denominator scale. The canonical hostile
below has **zero** external integer common neighbors before rescaling.

## 3. Proof of the rational extension criterion

### 3.1 Three square equations and a right triangle

A common-neighbor equation has the form

`x+a=p², x+a+d=q², x+a+2d=r²`.

Thus the three rational squares are in AP with difference `d`. Choose
nonnegative roots. Then `(r+p,r-p,2q)` is a positive rational right triangle,
since

`(r+p)²+(r-p)²=2(r²+p²)=4q²`,

and its area is `(r²-p²)/2=d`. This proves necessity. A putative `p=0` would
force the rational ratio `r/q` to equal `sqrt(2)`, so that degeneracy cannot
occur when `d>0`.

Conversely, a positive rational right triangle `(A,B,C)` of area `d` gives

`p=|A-B|/2, q=C/2, r=(A+B)/2`,

with `q²-p²=r²-q²=AB/2=d`. Therefore

`x=p²-a`                                                        (1)

solves all three common-neighbor equations. One such triangle need not make
`x` positive or distinct from the anchors; the next step supplies infinitely
many choices and handles both conditions.

### 3.2 Infinite-order input and real recurrence

Use the congruent curve

`E_d: Y²=X³-d²X`.

A positive area-`d` right triangle gives the rational point

`X=d(A+C)/B, Y=2dX/B`, with `X>d, Y!=0`.

This point has infinite order. To specify the dependency exactly, write
`d=N h²`, with `N` a positive squarefree integer and `h>0` rational. The
change `(X,Y)=(h² U,h³ V)` identifies `E_N` and `E_d`. Section 3 of the
[congruent-gap packet](congruent_square_gap_classes_20261019.md)
proves directly by doubling and 2-adic valuations that every rational point
on `E_N` with nonzero y-coordinate has infinite order: after the first
doubling its x-coordinate has negative 2-adic valuation, and each further
doubling lowers that valuation by two. No rank, Mordell--Weil generation, or
BSD input is used here.

For completeness, the needed recurrence is short. The nonsingular projective
real curve is compact, and its group law is continuous. A subsequence of the
distinct multiples `nP` converges. Differences of consecutive terms of that
subsequence give nonzero rational points `R_j -> O`. After taking another
subsequence these are distinct. As `R_j` approaches `O`, its affine
x-coordinate `X_j` tends to positive infinity: real affine points have
`X in [-d,0] union [d,infinity)`, and the bounded-X portion stays in a compact
affine set disjoint from `O`.

The standard real elliptic-group structure used here is described in
[Elkies, Math 223 (2024), the September 18 group-law and real-curve discussion](https://people.math.harvard.edu/~elkies/M223.24/index.html).
The recurrence argument itself is included above; no equidistribution or
claim about the doubling subsequence is needed.

### 3.3 The common-neighbor coordinate is unbounded

For `R=(X,Y)` with `X>d`, the area-`d` sides are

`A=|(X²-d²)/Y|, B=|2dX/Y|, C=|(X²+d²)/Y|`.

Since `Y²=X(X²-d²)`, as `X -> infinity` we have `A²~X` and `B²~4d²/X`.
Consequently `p²=(A-B)²/4 -> infinity`. Equation (1) supplies positive,
unbounded rational common neighbors. A finite forbidden set can be avoided,
including all three anchor values and all previously chosen neighbors.

An exact useful coordinate identity is

`x = X(2R) - (a+d)`,                                           (2)

because

`X(2R)=(X²+d²)²/(4Y²)=C²/4=q²`.

This also gives an effective unbounded search: enumerate positive multiples
of the supplied infinite-order point, evaluate (1) or (2), and keep fresh
positive values. Recurrence proves eventual success for every finite target
count. No bound on the first successful multiplier or output height is claimed.

Together, these arguments prove Theorem 2's three-way equivalence.

## 4. Passage back to one integer graph

For an integer AP square-sum triangle, write its edge sums as

`R²=2a+d, S²=2a+2d, T²=2a+3d`.

The integer right triangle `(T+R,T-R,2S)` has area `d`. Hence Section 3 applies.
Notice that this initial right triangle produces `x=a+d`, the middle anchor:
it solves the three square equations but is not a graph common neighbor in
the simple graph. Distinctness cannot be silently dropped.

Choose any finite list of fresh positive rational neighbors `x_1,...,x_k`,
with respective rational roots `p_i,q_i,r_i`. Let `L` be a common multiple
of all those root denominators in lowest terms. Then

`L² x_i=(L p_i)²-L² a`

is a positive integer, as are the three scaled anchors. The three cross-edge
roots are the integers `L p_i,L q_i,L r_i`; the core-edge roots are `LR,LS,LT`.
Multiplication by `L²` preserves all inequalities and distinctness. This
proves Theorem 1.

The same argument for finitely many rational anchors also works after
including their denominators in the chosen `L`; square scaling preserves
every required rational-square edge. The integer-core formulation above is
the one needed for the stated extension.

### Graph consequences and limits

- With one leaf, the scaled original triangle is contained in a `K_4`.
- With three leaves, discard the three core edges to obtain `K_{3,3}` with
  bipartition the three core vertices and the three leaves. Thus every AP
  square-sum triangle has a square dilation contained in a finite nonplanar
  square-sum graph.
- With `k` leaves, the guaranteed `K_3 join I_k` has chromatic number exactly
  four for `k>=1`; extra ambient edges are uncontrolled. Increasing `k` does
  not prove increasing chromatic number.
- No `K_5` is supplied. These are square-sum graphs, not planar unit-distance
  graphs: a `K_4` cannot be realized with all six distances equal in the
  Euclidean plane. There is no whole-graph adjacency-preserving transfer of
  these graphs to the planar unit-distance graph; this does not exclude
  transfers of selected smaller subgraphs.
- The construction does not provide a Hamiltonian cycle through every label
  of an interval graph `Q_N`, bounded-height realizations, or LRC(14).

## 5. Why four AP anchors are impossible

If a rational `x` were adjacent in the square-equation sense to all four
anchors, the four distinct numbers

`x+a, x+a+d, x+a+2d, x+a+3d`

would be rational squares in nonconstant AP. This is exactly excluded by
Fermat's four-squares theorem.

**Primary proof source:** Alf van der Poorten, *Fermat's Four Squares
Theorem*, arXiv:0712.3850 (2007), theorem headed “Fermat's four squares
theorem,” immediately after the introduction:

https://arxiv.org/html/0712.3850v1

The source proves the rational statement by clearing denominators and
descent. We use this classical theorem, not a new proof of it. Its separate
three-squares/area correspondence agrees with Section 3.1.

## 6. A complete unscaled hostile

For the AP core `(482,3362,6242)`, the gap is `2880` and the core-edge roots
are `62,82,98`. If an integer `x` satisfied all three common-neighbor
equations, their rational square roots would be integers. Put

`q²-p²=2880`, so `(q-p)(q+p)=2880`.

There are precisely 15 unordered positive factor pairs of `2880` with
matching parity. Testing the third square on all 15 gives only

`(x,p,q,r)=(3362,62,82,98)`.

Since `3362` is the middle anchor, the unscaled simple graph has **no external
integer common neighbor** of this core. An independent scan of all
`0<=p<=719` gives the same unique equation solution; the bound follows from
the matching-even-parity factor pair, whose smallest factor is at least two.

More generally, any fixed distinct integer anchors already have only
finitely many integer common neighbors, by factoring their fixed difference
as a difference of two squares. Thus a fixed `L` cannot serve arbitrarily
large `k`. The finite-family rescaling is logically necessary.

Quantitatively, for two positive integer anchors with difference `D>0`, a
positive integer common neighbor gives `0<p<q` and a strict factor pair
`q-p<q+p` of `D`, necessarily of matching parity. There are at most
`floor(tau(D)/2)` such pairs, where `tau` counts positive divisors. Therefore
the square-scaled AP core has at most `floor(tau(L²d)/2)` positive integer
common-neighbor equation solutions. Its middle vertex is always one such
solution and must be excluded, giving the slightly sharper graph bound

`number of external integer common neighbors <= floor(tau(L²d)/2)-1`.

In particular, an extension with `k` leaves requires `tau(L²d)>=2(k+1)`.
The strict factor inequality matters: allowing a zero first root for
arbitrary signed labels would add the equal-factor case when `D` is a square.

## 7. Explicit six-vertex certificate from class 5

On `E_5`, take `P=(25/4,75/8)`. Multiply curve coordinates by `(24²,24³)`
to get `E_2880` points. At the multiple `nP`, formula (2) is

`x_n=576 X(2nP)-3362`.

The first multiple returns the excluded `x_1=3362`. The second gives the
negative value `x_2=-1856688659521/3878673841`. The next three work.

For each row below the three common-neighbor roots are the listed numerators
divided by the common denominator `D_n`, and `x_n=(p_n)²-482`.

| n | D_n | numerator of p_n | numerator of q_n | numerator of r_n |
|---|---:|---:|---:|---:|
| 3 | 14896790112963 | 1036987385464258 | 1309372438208722 | 1534134780998498 |
| 4 | 223092910950900087178802 | 249563579992463717493803519 | 249850594047271558364480641 | 250137278774864229623059201 |
| 5 | 4168056500083002214944744728553146087995 | 230076377241990452361605372947651026178498 | 320887053228866028337428058295226484802002 | 391154525085689757012277699003111694342498 |

Choose

`L=lcm(D_3,D_4,D_5)`

` =13851986720466207228947843753473356093372173219780057654419750125010709636370`.

The six positive integer labels are exactly

`L²*(482,3362,6242,x_3,x_4,x_5)`.

This is a complete factored certificate; the companion JSON additionally
prints every integer label and all 12 required integer edge roots. All three
leaf-pair sums are nonsquares, so **this particular six-vertex example** is
even an induced `K_3 join I_3`. That extra finite fact is not included in the
universal theorem.

For a smaller `K_4` certificate, use `L=D_3` and only `x_3`. In both examples
the core is exactly one square dilation of the original specified core.
The four labels for this smaller certificate are

```text
106962719432781995724908175858
746076063761437903790749558578
1385189408090093811856590941298
968380118179215608261287314706
```

In this order, the roots of edge sums `(0,1),(0,2),(1,2),(0,3),(1,3),(2,3)` are

```text
923600987003706, 1221536789262966, 1459885431070374,
1036987385464258, 1309372438208722, 1534134780998498.
```

## 8. Reproduction and exact test scope

From the repository root, after installing the files at the indicated paths:

```bash
python3 04-computation/experiments/ap_triangle_common_neighbors_20261019.py
python3 -O 04-computation/experiments/ap_triangle_common_neighbors_20261019.py
```

The program uses only the standard library and keeps its explicit checks
under optimization. Normal and optimized runs produce byte-identical output.

The complete Euclidean-parameter universe `2<=m<=30`, `1<=n<m`, coprime
opposite parity, filtered by strict AP-core positivity, has 66 primitive
cores. For each, all first eight elliptic multiples are checked, giving 528
rational common-neighbor equation samples. The tests verify area, curve
membership, the duplication coordinate, both AP differences, all three
square equations, collision exclusion, and integer square scaling.

Within that bounded multiple window, selecting at most three fresh positive
neighbors gives counts `0:4, 1:5, 2:8, 3:49`. The zero cases are not
counterexamples to the unbounded theorem. They explicitly prevent mistaking
this finite window for a uniform first-eight conversion guarantee.

Further controls verify the 12-edge certificate, its three nonedges, the
complete unscaled hostile by two independent paths, and the criterion
without a triangle: anchors `(1,7,13)` have positive common neighbor
`1305201/19600`, with roots `(1151,1201,1249)/140`, although `1+7` is not a
square. The four-anchor sample `1/4,25/4,49/4,73/4` checks a concrete failed
fourth square; universal impossibility uses the cited theorem.

Total: **6,615 exact checks**. These finite tests verify the formulas and
certificates; the compactness argument, not enumeration, proves infinity.

## 9. Connection contract

- **Source:** three positive rational AP anchors with congruent gap class;
  in the main application an existing positive integer AP square-sum clique.
- **Target:** arbitrarily large common-neighbor families and, after finite
  scaling, `K_3 join I_k` square-sum subgraphs.
- **Map:** right-triangle/three-square conversion; elliptic recurrence;
  `x=X(2R)-(a+d)`; one common denominator square.
- **Preserved:** the distinguished AP-core shape and gap squareclass, strict
  positivity, all required edges, and distinct labels.
- **Lost:** fixed integer scale, output height, and control of leaf-pair edges.
- **Necessary sidecars:** the finite selected family, rational square roots,
  their common denominator, and exclusion of the middle-anchor collision.
- **Cheapest decisive hostiles:** the original core has zero external integer
  common neighbors; the second class-5 multiple is negative; four AP anchors
  are impossible regardless of scaling.

The new graph operation extracts a nonplanarity consequence from the same
congruent-number carrier that previously classified AP gap classes. It does
not solve the congruent-number problem or any of the unrelated global
conjectures mentioned above.
