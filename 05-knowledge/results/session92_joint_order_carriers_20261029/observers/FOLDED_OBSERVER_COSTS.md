# Folded tournament paths: observer cost and a native quadratic obstruction

**Status:** finite-exact, independently reconstructed and executable. All claims
below concern the specified 64-state carrier, not arbitrary tournament sizes.
The 16/35 targets were supplied by the user; their original latest packet was
not located by this lane. The 64/24 construction is inherited from
`session90_even_graphs_odd_descent_20261029/oddness/CENTRALIZER_ORBIT_PATH_DP.md`
and its checker. Source snapshot inspected: `b07ec51b204`.

## Construction and contract

Let W=C3[C3,C3,C3], vertex v=3i+j with i,j modulo3. Within a child use
j→j+1; between children use i→i+1. The anti-involution is τ(i,j)=(-i,-j),
with fixed root0. Ordered pair representatives give
`[(1,2),(3,6),(4,8),(5,7)]`.

A state (M,v) consists of a nonempty mask of used τ-pairs and a current vertex
belonging to M. Its successor (M∪{pair(w)},w) exists iff that pair is unused
and v→w. There are64 states and96 directed transitions. Level is |M|; final
acceptance means M=15 and v→0. All singleton-mask states start with mass1.
A successful half path determines the entire τ-fixed Hamilton path, so the
exact count is33.

A partition is admissible iff its blocks have constant level and acceptance,
and every state in a block has the same total transition weight into each
block. This is a forward equitable partition, with outgoing row sums, NOT a
claim that all intermediate state masses agree. To propagate arbitrary source
masses, initialize a block with the SUM of its members' source masses.

The root/τ centralizer C3 has24 state orbits. Its orbit partition is admissible.
Further equitable refinement from output colours yields

| observer menu | refinement history | final blocks |
|---|---|---:|
| level, acceptance | 5→10→15→16 |16|
| plus current vertex=3 |9→23→34→35|35|

All three 64/24/16 carriers compute33. The35 carrier also computes33 but
supports the extra named test. The new observer can split old symmetry orbits;
there is no contradiction in35 exceeding24.

**Minimality certificate.** Start with output equality P0 and repeatedly colour
x by its previous colour and its vector of outgoing weights into previous
blocks. Every admissible partition refines P0. If it refines Pt, its block
weights determine all Pt-block weights, so it refines Pt+1. Thus it refines
the stable result. Stability makes that result admissible. This proves the
coarsest admissible partition is unique and minimal in number of blocks.
The supplied checker performs every refinement and checks every final block.
The depth-three stabilization is exact here, not a universal depth bound.

## A quadratic obstruction in the actual tournament carrier

Write A f(x)=Σ_{x→y}f(y). Begin with indicator functions of all output colours.
The rational span of A^k of these indicators is a *linear future-observer
space*. The nilpotence A^4=0 makes k=0,...,3 sufficient.

For the base menu its dimension is12, while the partition algebra dimension
is16. Pointwise products of pairs of linear observables already span all16.
The12-dimensional space nonetheless has16 distinct evaluation signatures:
there is no pairwise blindness at this base menu, just linear dependence.
For example g=A1_{level3}; g² is outside that12-dimensional span.

For the vertex3 menu, the distinction becomes a concrete hostile:

* linear future-observer dimension21;
* distinct linear evaluation signatures34;
* minimal equitable partition35.

Take x=(8,7), y=(4,8). They agree on **every** A^k applied to **every** named3
output-colour indicator. Let f indicate `(level4, rejecting, current≠3)`.
The three successors of x are (9,1),(9,2),(12,8); those of y are
(5,1),(5,2),(6,6). Their values under A²f are respectively

```
             x successors       y successors
A²f             1,1,1              0,0,3
sum               3                  3
sum of squares    3                  9
```

Hence A((A²f)²)(x)=3 but A((A²f)²)(y)=9. Linear total counts miss the
concentration of continuation mass at a shared successor. The product
observes two independent downstream choices sharing their first successor:
it is a rooted branching query, not two unrelated marginal samples.

**Important limitation:** this linear space only permits powers of A followed
by one output indicator. It is not automatically the Hankel-minimal
representation for arbitrary words interleaving transitions and named tests.
Likewise16 and35 are minimal *state partitions*, not lower bounds for every
rational linear realization, nor faithful tournament reconstruction results.

## All named-current-vertex menus, finite-exact

For each subset J⊆{1,...,8}, add all predicates v=j for j∈J and run the same
minimality procedure. All256 menus are tested. The range of minimal partition
sizes is:

| Number of names | smallest | largest |
|---:|---:|---:|
|0|16|16|
|1|21|37|
|2|22|46|
|3|44|53|
|4|52|58|
|5|54|62|
|6|57|64|
|7|64|64|
|8|64|64|

Single names1,2 cost21; names3,4,5 cost35; names6,7,8 cost37.
Exactly six inclusion-minimal menus recover all64 states. Each has six names;
their bitmasks (bit j−1 means name j) are126,190,222,237,245,249.
Equivalently their omitted pairs are {1,8},{1,7},{1,6},{2,5},{2,4},{2,3}.
Any seven names suffice. This does NOT say six arbitrary bits are necessary:
it concerns this named-predicate vocabulary and this fixed carrier only.

The block-count cost c is not submodular:

`c({1})+c({3})=21+35=56 < c(∅)+c({1,3})=16+43=59`.

Nor is it supermodular: c({1})+c({2})=42 >16+22=38.
One cannot safely optimize sensor choice by assuming diminishing returns or
increasing returns. The whole finite menu table is cheap enough to enumerate.

## Separate stable observations can fail under conjunction

The checker includes a tiny weighted-DAG hostile with source states x,y and
four terminal states carrying coordinate labels00,01,10,11. x has successors
00,11; y has successors01,10. Either coordinate alone gives a stable three-block
partition and merges x,y. The common refinement of these two stable partitions
has five blocks but is NOT stable. Joint refinement has six, separating x,y.
Thus closing each marginal and intersecting once is not a general algorithm
for the jointly closed carrier. Closure must run again after the merge.

## Links to earlier exact carriers (not numerical analogies)

* **Zeckendorf/CRT:** `session91_packing_golden_carriers_20261029/zeckendorf/`
  proves arithmetic register sharing via gcd/lcm compatibility, with exact
  minimal complete-DFA counts. This supplies an algebraic *reason* that some
  joint state information is redundant. A raw product is an upper bound, not
  an exact minimum. Here the vertex-menu census shows both redundancy and
  synergy; the correctness criterion is weighted outgoing equivalence, not
  the deterministic-language Myhill–Nerode criterion.
* **Golden-fibre repair:** in that packet's `golden_fibres/`, a single named
  toggle plus the cyclic reader separates the three zero-fibre words, forcing
  a larger operation-closed carrier. The present named vertex3 test enlarges
  16 to35. Both are explicit operation-relative information costs, but their
  state spaces and operation algebras differ.
* **Joint run order:** the current user's Hamilton-path counterexamples have
  equal marginal cursor responses yet unequal assemblies. The native x/y
  calculation above gives a precise adjacent mechanism—linear sums forget
  shared-successor concentration. It is not a proof that those assembly
  witnesses are identical to this quadratic witness.

## Loss ledger and verification

Source64states → base16blocks preserves source-mass aggregation, level, final
flag and weighted transitions; loses named-current-vertex queries. Restoring
name3 gives35, and selected six-name menus give64. Centralizer compression24
retains more distinctions than base16 but cannot implement arbitrary named
marks. Linear12-space loses closure under multiplication; quadratic products
restore the whole base16 partition algebra.

`check_observer_algebra.py` uses exact integers and rational elimination,
stdlib only, with100 explicit active checks in both normal and `python3 -O`
runs. `normal.json` and `optimized.json` are byte-identical. Individual loop
comparisons inside those checks are not counted as additional checks.
The parent audit lane independently rebuilt the folded refinement and
rational ranks using a separate implementation.
