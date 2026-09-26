# 233, triangle claims, and height-dependent components of a marked Collatz path

**Status: PROVED, root-audited elementary reconstructions and fixed-word monotone cut theorem;
FINITE-EXACT independent controls stated below. No Collatz closure.**
2026-09-26, bounded graph/geometry lane. No canon ID is claimed here.

The useful reconstruction is an exact correspondence of **edge coordinates**,
followed by an exact threshold formula for strongly connected components of
a particular Collatz comparison tournament. Triangle inequality, a unique
line, a unique path, and absence of loops are different properties.

## 1. Inheritance and scope

* Closest proved graph mechanism: [THM-354, good-cut SCC count](../../01-canon/theorems/THM-354-good-cut-scc-count.md).
  In a tournament with a directed Hamiltonian path, component boundaries are
  precisely path cuts crossed by no backward arc. The interval fact below is
  inherited, not a new SCC theorem.
* Closest proved arithmetic mechanism: the affine gate and direction lemmas in
  [the procedural order-law search, sections 3 and 4.6](collatz_procgen_20260922_order_laws.md).
  Fixed-word plus-sheet descent sets are upward closed in their residue class.
  Our contribution is the explicit component-threshold synthesis of this
  mechanism with THM-354, including a cylinder-level 233 certificate.
* Corrected near miss: [the compression audit, section S6](collatz_mod6_20260922_compression_and_lean_audit.md)
  and [the scaffolding audit, theorem 5.1](collatz_mod6_20260917_scaffolding_audit.md)
  already show that fixing a Hamiltonian spine leaves `binom(N-1,2)` free
  chords. A layout or a marked path is not free information.
* Canonical hostile: the backward edge in the reduced two-vertex tournament
  lifts to a cyclic triangle. The reduced tournament is acyclic nonetheless.
  On the arithmetic side, the already recorded `231 -> 233` decay-clock
  crossing is the positive-carry hostile; its discovery belongs to the
  procedural search above, not this session.
* Least-used relevant sidecar: retain the **labelled cut thresholds**, not just
  the number of SCCs or the coefficient signs. This reveals exactly where a
  fixed word changes its comparison geometry as source height increases.

Concept board: metric betweenness; marked path/chord coordinates; SCC cuts;
positive affine carries; exact valuation cylinders; the two meanings of 233.
No outside literature theorem is imported: the elementary claims used here
are proved below or attributed to the precise repository result above.

## 2. Minimal hostile examples and strongest metric repair

A loop-free complete graph with one edge per pair can have edge lengths
1,1,3. Its edge data violate the triangle inequality 3<=1+1. Thus edge
uniqueness and absence of loops do not imply a metric inequality. A tournament
specifies directions between distinct vertices, and supplies no edge lengths
until those are defined separately.

1. Three equidistant points satisfy the triangle inequality. Their complete
   weighted graph has a cycle, and each pair has a unique shortest path (the
   direct edge). Thus triangle inequality, even together with unique shortest
   paths, does not imply absence of graph cycles or a line embedding.
2. A connected undirected graph has a unique **simple** path between every
   pair exactly when it is a tree: a cycle supplies two simple paths, whereas
   two different simple paths contain a cycle. The four-vertex star is the
   smallest branching tree, so this condition does not make the graph a line.
   A finite tree is a path exactly when every degree is at most two.
3. Four points with the shortest-path metric of a unit square have adjacent
   distances one and opposite distances two. Every triangle is degenerate,
   with side lengths `(1,1,2)`, but the metric does not embed in a line. If two
   opposite points are put at coordinates zero and two, each of the other two
   must be at coordinate one, contradicting their mutual distance two.
   This is minimal: every degenerate three-point metric embeds in a line.

The exact line statement is: a finite metric admits an isometric embedding in
the real line iff there is **one common ordering** `v_1,...,v_N` for which

```
d(v_i,v_k) = d(v_i,v_j) + d(v_j,v_k)       whenever i<j<k.
```

The forward implication follows from ordered real coordinates. Conversely,
put `v_1` at zero and `v_j` at the sum of the first `j-1` adjacent distances;
the equalities prove that every pairwise distance is preserved. Positive
adjacent gaps follow from the metric axiom. A single compatible betweenness
order is the missing information in the square example.

Even on a line, completing the metric to a weighted complete graph need not
give unique shortest combinatorial paths: three points at `0,1,2` have two
shortest routes from zero to two, one direct and one through one.

## 3. The exact smaller-tournament chart

Fix a marked directed spine `v_1 -> ... -> v_N`, with `N>=2`. For `1<=a<b<=N-1`,
send the reduced edge `{a,b}` to the original chord `{v_a,v_(b+1)}`. Declare

```
a -> b in R     iff     v_a -> v_(b+1) in T.
```

The inverse sends a nonadjacent pair `{v_i,v_j}`, `i<j-1`, to `{i,j-1}`.
Thus this is a bijection between all orientation arrays on the reduced
`K_(N-1)` and all tournaments containing the marked spine. It accounts for

```
binom(N,2) = (N-1) + binom(N-1,2)
```

edge coordinates. The chart has `2^binom(N-1,2)` possibilities. The marked
order is a sidecar: forgetting it gives variable fibres equal to the number
of Hamiltonian paths. The prior double count is
`sum_T H(T)=N! * 2^binom(N-1,2)`; it is not a new compression principle.

This is **not a vertex map or graph quotient**. In fact, a reduced triangle
on `a<b<c` maps to the four-vertex undirected path

```
v_(b+1) -- v_a -- v_(c+1) -- v_b.
```

Incidence and directed-cycle structure are not preserved. For `N=3`, the
single reduced edge `2->1` is acyclic but its lift is
`v_1->v_2->v_3->v_1`. The precise survivor is:

> The lift is acyclic iff every reduced edge points in increasing **label**
> order. Merely being an abstractly transitive reduced tournament is weaker.

Any backward original chord closes the corresponding spine interval into a
directed cycle; if there are no backward chords, the given ordering is a
transitive ordering. This proves the assertion without a counting argument.

For completeness, in a tournament the following are equivalent: acyclic;
transitive; no directed triangle; a unique directed Hamiltonian path.
A shortest directed cycle is a triangle, since either orientation of a chord
in a longer shortest cycle creates a shorter cycle. A triangle has three
Hamiltonian paths. Each such path extends to a Hamiltonian path on all
vertices by successive insertion of the remaining vertices: insert a new
vertex before the first existing vertex that it beats, or append it if there
is none. The three restrictions on the original triangle stay distinct.
Thus a unique Hamiltonian path excludes a triangle; the other directions
follow from a transitive order. Undirected cycles are a separate notion:
every tournament on at least three vertices has complete underlying graph.

### Which triangle tests actually suffice?

Let `a_ij=1` when `v_i->v_j`, for `i<j`. A triple `i<j<k` is cyclic precisely
when `a_ij=a_jk != a_ik`. Consecutive triples alone do not suffice: on four
vertices, orient everything forward except `v_4->v_1`. The consecutive
triples are transitive, but `(v_1,v_2,v_4)` is cyclic.

The triangular test family `(v_i,v_(i+1),v_j)` for all `j>=i+2` does suffice.
With `a_i,i+1=1`, absence of its cyclic triangle implies
`a_ij >= a_(i+1),j`. Starting from `a_(j-1),j=1` forces every `a_ij=1`.
This is a graph transitivity test with `binom(N-1,2)` labelled tests; it is
not an application of metric triangle inequality.

## 4. Cut thresholds for a fixed positive Collatz word

Let `w=(k_1,...,k_(N-1))`, with every `k_i>=1`. Its positive-sheet affine
nodes, initially treated as functions of a positive real variable `x`, are

```
m_0(x)=x,
m_j(x)=a_j x+b_j=(3^j x+C_j)/2^S_j,
S_j=sum_(ell=1)^j k_ell,
C_j=sum_(ell=0)^(j-1) 3^(j-1-ell) 2^S_ell.
```

Define the **marked chronological comparison tournament** `T_w(x)`:

* vertices are the positions `0,...,N-1`;
* adjacent positions always have the spine arc `i->i+1`;
* for nonadjacent `i<j`, orient `i->j` iff `m_i(x)>m_j(x)`.

Exclude the finite set of positive **nonadjacent** contacts `m_i(x)=m_j(x)`
for now; this avoids a chord tie convention. Vertices remain labelled by
position. This graph uses a declared chronological gauge on
adjacent pairs and actual height descent on the chords. The tournament using
height comparisons on **all** pairs would simply be transitive and generally
would not contain the chronological spine. Its type must not be confused
with `T_w(x)`.

The normalized intercepts increase strictly:

```
b_j/a_j = C_j/3^j = sum_(ell=0)^(j-1) 2^S_ell / 3^(ell+1).
```

All slopes are distinct, since a nonzero power of three is not a power of
two. For each nonadjacent pair `i<j`, define the extended threshold

```
theta_ij = infinity                          if a_j>a_i,
theta_ij = (b_j-b_i)/(a_i-a_j)               if a_j<a_i.
```

If `a_j>a_i`, increasing normalized intercepts give `b_j>b_i`, so the future
node is larger for every positive `x`: the chord is permanently backward.
If `a_j<a_i`, subtracting the affine functions shows that the chord is
forward exactly when `x>theta_ij`. Negative or zero thresholds mean it is
already forward for all positive `x`.

**PROVED monotone-cut theorem.** As `x` increases away from contacts, backward
chords can only disappear. The tournaments form a nested chain of chord
arrays with at most `binom(N-1,2)+1` different values. Their strongly
connected components can only split. For the cut immediately before vertex
`r`, `1<=r<=N-1`, set

```
Theta_r = max({0} union {theta_ij : i<r<=j, j>=i+2}).
```

Then that cut is a component boundary exactly when `x>Theta_r`.
Consequently there are at most `N-1` distinct positive component-splitting
heights, and at most `N` component partitions.

Proof: a backward chord `j->i` and the forward spine make the entire interval
`[i,j]` strongly connected. As in THM-354, components are consecutive spine
intervals and a cut separates them iff no backward chord crosses it. The
threshold rule makes this equivalent to `x>Theta_r`. Once open, a cut never
closes. An infinite threshold is a permanent obstruction at that cut. This
also proves the two event bounds. The asymptotic partition is determined by
the intervals of chords with `a_j>a_i`, since all contracting chords are
eventually forward. QED.

This theorem is a finite-word statement, not a claim about changing the
starting height along one orbit. The values of `x` giving the exact odd
valuation word form one odd residue class modulo `2^(S_(N-1)+1)`. This follows
inductively: after a realized prefix, increasing the source by its modulus
changes the final odd node by twice an odd multiple; the next exact valuation
selects one of the `2^k` further residue classes. Restricting the real theorem
to that class preserves it, but some real chambers may contain no realization.

## 5. The recovered 231-to-233 crossing, with a new cut certificate

The prior procedural search records the exact nodes

```
231,347,521,391,587,881,661,31,47,71,107,161,121,91,137,103,155,233
```

and word

```
(1,1,2,1,1,2,6,1,1,1,1,2,2,1,2,1,1).
```

Thus `S_17=27` and direct integer arithmetic gives

```
m_17(x)=(129140163 x+1441352971)/134217728,
theta_(0,17)=1441352971/5077565.
```

The coefficient contracts but the endpoint rises at `x=231`, since this
source lies below the contact. This is **not** a no-descent orbit segment:
it already visits 31 at position seven. The old crossing and the old dip
are inherited, not discoveries here.

The exact new component certificate is:

* `Theta_7=1441352971/5077565`, attained by the chord `(0,17)`;
* every other `Theta_r` is infinite;
* the same rational number is the largest positive contact of any pair of
  affine nodes, including adjacent ones.

These are finite rational inequalities checked explicitly by the companion
script, independently of the earlier procedural engine. They imply:

```
x=231:              one strong component, all 18 positions;
x=231+2^28 t,t>=1:  exactly two components, [0,...,6] and [7,...,17].
```

There is also a short direct proof of the component certificate. The backward
chords `(0,6)` and `(7,17)` are permanent because their block multipliers are
`729/256>1` and `59049/8192>1`. Together they cover every cut except cut seven
and make each of the two displayed intervals strongly connected. The slope
`a_0` is the least among positions zero through six, so `m_0(x)` is the least
height in that interval for every `x>0`. The slope `a_17` is the greatest
among positions seven through seventeen; the increasing normalized carries
therefore make `m_17(x)` the greatest height in that interval. Hence every
cross-interval comparison points forward exactly when `m_0(x)>m_17(x)`.
If the reverse strict inequality holds, chord `(0,17)` joins the intervals
into one strong component. This yields the claimed unique finite cut
threshold without checking the cross-interval contacts one at a time.

Every `t>=0` realizes the exact word: at each position the change from the
`t=0` orbit is `3^j 2^(28-S_j)t`, an even integer. Hence the intermediate
nodes remain odd and all prescribed divisions are exact. No pair collides
on this cylinder: the `t=0` nodes are distinct, and all `t>=1` sources exceed
the largest positive contact. The endpoint formula is particularly simple:

```
m_17(231+2^28 t)=233+258280326 t,
x-m_17(x)=-2+10155130 t.
```

So the larger same-word sources descend at this endpoint, and the SCC cut
opens exactly at the earlier dip. The component change is a concrete result
of retaining source height and carry alongside the word.

## 6. The other recovered 233 is a Hamiltonian-path count

[THM-316, staircase antipalindrome](../../01-canon/theorems/THM-316-staircase-antipalindrome.md)
defines a different exact object: the all-zero interleaved staircase on eight
vertices has **233 Hamiltonian paths**, with starting-vertex counts
`[80,29,39,24,21,18,11,11]`. A fresh subset dynamic program reconstructs this
from the literal pair/rank edge definition. The neighboring counts are
`1,5,29,233,2489` for `2,4,6,8,10` vertices, so this is not the Fibonacci
sequence merely because `233=F_13`.

Do not substitute the separate base-path staircase:
[THM-337, base-path staircase recurrence](../../01-canon/theorems/THM-337-base-path-staircase-recurrence.md)
is explicitly conjectured and its eight-vertex Hamiltonian count is 57.
Neither the count 233 nor a vertex value 233 supplies a map between these
objects. The actual transfer in section 4 is the cut/carry mechanism and
works for arbitrary fixed words, with 233 providing a useful inherited test.

## 7. Connection contract and stopping boundary

| Item | Exact content |
|---|---|
| Source | Positive-sheet affine odd-Collatz word, including its carries and real source height |
| Target | Tournament on chronological positions, with declared adjacent spine and height-oriented nonadjacent chords |
| Map | Evaluate affine nodes; compare each nonadjacent pair; apply the edge chart from section 3 if desired |
| Preserved predicate | Every nonadjacent block's actual height-descent sign; SCC cut structure under this observer |
| Destroyed information | Magnitudes, primes, exact residues if the word is omitted, and the actual signs of adjacent steps |
| Required sidecar | The word, exact realization class, height parameter, affine carries, and adjacent-step signs when discussing descent |
| Cheapest hostile test | Reduced backward K2 lifts to C3; fixed word of 231 gives different SCCs at 231 and 231+2^28 |

The observer's cycles are **comparison cycles**, using chords that are not
Collatz edges. They are not periodic dynamical orbits. An acyclic observer
requires every nonadjacent future node to be smaller than its earlier node,
which is far stronger than eventual descent. Its first spine arc even
forgets whether the first actual Collatz step rises or falls. Pure height
tournaments are always transitive and supply no such criterion either.

The strongest survivor is an exact, monotone, carry-sensitive component
filtration along one valuation cylinder. Obtaining a uniform bound over all
long words, or proving that every positive orbit eventually descends, remains
OPEN. The finite-word filtration supplies neither implication. No new open
conjecture is promoted merely from the numerical recurrence of 233.

## 8. Reproduction and independent controls

Run from the repository root:

```
python 04-computation/experiments/crossroads233_20260926_graph.py
```

Companion output: [crossroads233_20260926_graph.out](crossroads233_20260926_graph.out).
The script uses exact integers and fractions, no external packages, and
checks retained under `python -O`. Its explicit universe is:

* every fixed-spine tournament on `N=2,...,7`: respectively
  `1,2,8,64,1024,32768` orientation arrays, with SCCs computed independently
  by transitive reachability; the triangle test family is checked on all;
* uniqueness of Hamiltonian paths via independent subset DP through `N=6`;
* the two metric hostile controls and all 24 possible line orders of the
  four-cycle metric;
* all 1,364 words of lengths 1 through 5 with exponents 1 through 4, at one
  exact rational sample in every positive real chord chamber, totaling 8,324
  samples; direct height comparisons agree with the threshold predictions,
  and reachability SCCs only split;
* fresh staircase counts for 2 through 10 vertices and the anti-automorphism;
* the 233 fixed-word rational cut certificate and four direct integer
  realizations `231+2^28 t`, `t=0,1,2,3`.

The four realization samples are finite controls; the preceding affine and
parity argument proves the asserted entire cylinder. The exhaustive finite
censuses check the elementary proofs rather than replace their quantifiers.
