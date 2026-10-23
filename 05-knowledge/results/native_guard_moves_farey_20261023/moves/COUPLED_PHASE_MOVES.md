# Sharp coupled-move budgets on native pinned cyclic-cell families

Status: **PROVED** for the actual regular hosts and supplied path families
below. The weight-filter extension is **CONDITIONAL** on an authenticated
factorization of the guard through the stated integer weights. Exact native
controls are **FINITE-EXACT**. No arbitrary electrical guard is replaced by
an integer weight, and no theorem about every Hamiltonian path or every
native metagraph is claimed.

## 1. Inheritance, carrier and move budget

The closest proved mechanism is the actual regular-cell construction and
three-state compiler in
[`session43.../native-branch-width/THEOREM.md`](../../session43_observers_phase_guards_20261023/native-branch-width/THEOREM.md).
Its binary difference words require exponentially many affine-coset
branches, although the native path family has a small transfer description.
The relevant boundary is
[`native-complexity/THEOREM.md`](../../session43_observers_phase_guards_20261023/native-complexity/THEOREM.md):
general shared C3 templates are not necessarily tractable chain families.
Our earlier
[`compiler/README.md`](../../hamiltonian_observer_restrictions_20261023/compiler/README.md)
has exact Hamming connectivity inside a native automorphism orbit, but an
arbitrary additional guard can split that orbit.

Here the anchor is fixed-source connectivity under actual coupled moves.
The niche is the integer lift of a finite cyclic phase path. The wildcard
is removing primality from the local move theorem. The live concept board
is **actual phase / integer displacement / pins / simultaneous support /
guard factorization / converse**.

Let `q,t>=3` be odd integers. The outer host is the regular circulant
tournament on `t` cells. Each cell is the regular circulant tournament on
`q` owners, with arcs `u -> u+s` for `1<=s<=(q-1)/2`. On every outer arc
`i -> j`, orient a cross pair by

\[
 (i,u)\longrightarrow(j,v)
       \quad\Longleftrightarrow\quad v-u\in\{1,2\}\pmod q.
                                                               \tag{1}
\]

Reverse all remaining cross pairs. This is a genuine regular tournament of
order `qt`: each score is

\[
 \frac{q-1}{2}+2\frac{t-1}{2}+(q-2)\frac{t-1}{2}
             =\frac{qt-1}{2}.
\]

For actual phases `a_0,...,a_r`, where `r=t-1`, concatenate the cyclic-cut
paths

\[
 (i,a_i),(i,a_i+1),\ldots,(i,a_i+q-1)
\]

in cell order. The boundary to the next cell is an arc precisely when

\[
                 d_i=a_{i+1}-a_i\in\{0,1\}\pmod q.           \tag{2}
\]

Thus these actual Hamiltonian observers correspond to `a_0` and a binary
difference word `d` of length `r`. Retain both; the common-phase lift is not
silently quotiented out. Primality is unnecessary for (1), (2), or any move
theorem below. It was needed in the inherited affine-coset-width proof,
which is a separate statement.

Pin actual first/last owners, hence `a_0=alpha`, `a_r=beta`. An indexed owner
at clock `iq+u` determines the phase `a_i` by its actual label. Indexed arc
pins determine the corresponding one or two phases in the same way.
Multiple observations of a shared cell refer to the *same* phase coordinate.

A **budget-b coupled move** changes at most `b` actual cell phases
simultaneously, keeps all pins, and has an accepted Hamiltonian observer at
both endpoints. Its certificate lists the changed cell rotations and
checks every resulting native boundary arc. It is an atomic block move;
the statement does not assert a sequence of accepted one-cell intermediate
states, an elementary tile flip, or a pre-existing metagraph edge.

This budget is literal Hamming support in the actual phase tuple. Each
nonzero phase change moves all `q` owners in that cell. It is not the
number of coordinates in a different hierarchical parameterization, and
is not the maximum distance between arbitrary two paths in the fibre.

## 2. Integer displacement is the missing component label

Let

\[
                   k(d)=\sum_{i=0}^{r-1}d_i\in\{0,\ldots,r\}.
\]

The endpoint pins say only `k(d)=beta-alpha (mod q)`. They lose the integer
lift of this displacement. With fixed endpoints, a one-cell move can only
replace adjacent differences `01` by `10` or conversely: changing the
middle phase adds a nonzero residue to one difference and subtracts it
from the next, and these are the only two distinct binary pairs with the
same sum modulo odd `q>=3`. Consequently one-cell moves preserve `k`.

Conversely adjacent swaps connect every binary word of a fixed weight `k`.
They sort it to `1^k0^(r-k)` in at most `k(r-k)` swaps. The components of the
one-cell move graph are therefore **exactly the feasible integer weights**,
not merely a collection of finite examples.

There is a native repair. On `q` consecutive differences replace

\[
                      0^q\longleftrightarrow1^q.              \tag{3}
\]

For the forward change, add `1,2,...,q-1` modulo `q` to the `q-1` internal
cell phases. Both endpoint phases stay fixed, every new difference is
binary, and all native path arcs remain valid. This is one legal coupled
move of exact budget `q-1`.

Swaps and (3) connect every fixed-endpoint fibre: sort a word, change the
weight by multiples of `q` at the unique `1/0` seam of its canonical word,
then unsort to the target. A path from weights `k` to `l` has length at most

\[
          k(r-k)+|k-l|/q+l(r-l).
\]

Hence a uniform diameter bound for this specified primitive move graph is
`2 floor(r^2/4)+floor(r/q)`. This bound counts block events; their phase
support is charged separately as `1` or `q-1`.

## 3. The budget q-1 is necessary, even for nonlocal moves

For two accepted words `d,e` with the same endpoint phases, define the
integer cumulative-difference walk

\[
 z_0=0,\qquad z_j=\sum_{i<j}(e_i-d_i).
\]

Its increments belong to `{-1,0,1}`, and the phase at index `j` changes iff
`z_j` is nonzero modulo `q`. If the weights differ by `gq>0`, the walk goes
from zero to `gq`. It must visit every integer `1,...,gq-1` at distinct
times. Of those integers, exactly `g(q-1)` are not multiples of `q`.
Therefore any such move changes at least

\[
                     g(q-1)                                  \tag{4}
\]

actual internal phase coordinates. Negative weight differences give the
same bound by reversing the two words. This applies to arbitrary
noncontiguous simultaneous support, not only to our chosen local moves.

The bound is attained between the canonical words of weights `k` and `l`:
flip the `|k-l|=gq` differences at their seam in one block. The cumulative
walk increases or decreases monotonically, hitting a multiple of `q`
exactly at each completed turn. Thus the exact minimum phase distance
between two different weight layers is

\[
                   \frac{q-1}{q}|k-l|.                        \tag{5}
\]

**Sharp unfiltered theorem.** The minimum budget making a pinned-endpoint
fibre connected is `q-1` when it contains at least two feasible weights,
one when it has a single weight but multiple words, and zero for a
singleton. An empty fibre has no connectivity claim.

For `q=3`, the first nontrivial regular-host example is `t=5`, order15.
Equal first/last phases permit weights zero and three among four
differences. One-cell moves cannot join them, whereas two-cell moves do.
The local witness `000 <-> 111` changes phases `(0,0,0,0)` to `(0,1,2,0)`.

## 4. Authenticated weight guards: larger blocks bridge missing layers

This section assumes an additional guard has actually been shown to factor
through `k`. Let its accepted feasible weights be a nonempty ordered set
`K={k_1<...<k_s}`. This factorization accepts the entire word layer at each
accepted weight, not an arbitrary subset of it.

If `s>=2`, the exact minimum connectivity budget is

\[
             \boxed{B=(q-1)\max_{i<s}\frac{k_{i+1}-k_i}{q}}.  \tag{6}
\]

If there is one accepted weight, use the singleton/nontrivial-layer rule
in section 3. For sufficiency, sort within a layer and jump between
consecutive accepted canonical weights with the single larger block of
(5). No rejected intermediate weight is visited: this is a simultaneous
coupled move. For necessity, cut `K` at a largest gap. Any edge crossing
this cut changes weight by at least that gap, so (4) forces budget at least
`B`. Thus rejected layers can impose an exactly computable larger block
cost; increasing one-cell search time cannot repair the obstruction.

## 5. Shared internal pins and correlated guards

Additional pinned phases divide the difference word into intervals of
lengths `r_1,...,r_h`, with prescribed endpoint residue differences. Retain
the actual pinned values. On each interval let `k_j` be the integer
weight. The compatible family at weight vector `k=(k_1,...,k_h)` is the
Cartesian product of its fixed-weight word layers, connected by the same
one-cell swaps. Its size is

\[
                         \prod_j\binom{r_j}{k_j}.
\]

For two feasible vectors `k,l`, the exact minimum actual phase distance
between their layers is

\[
                c(k,l)=\frac{q-1}{q}\sum_j|k_j-l_j|.          \tag{7}
\]

The lower bound applies (4) inside each independently pinned interval;
the changed internal coordinate sets are disjoint. Simultaneous canonical
seam flips on the intervals attain the sum, while every shared pin stays
fixed.

Suppose an authenticated guard factors through this **entire weight
vector**, accepting a set `K`. For any budget `b>=1`, the actual guarded
move graph is connected iff the graph on accepted vectors with edges
`c(k,l)<=b` is connected. Necessity follows by projecting every move to its
weights. Sufficiency follows by within-layer swaps and the attaining
canonical jumps. Equivalently the minimum positive budget is the maximum
edge cost of a minimum-bottleneck spanning tree on `K`; a lone nontrivial
layer needs one, and a lone actual state needs zero.

An exact native order27 hostile makes the coupling necessary. Take `q=3`,
`t=9` and pin phases zero at cell indices `0,3,5,8`. Require equality of the
integer weights in the first and last intervals. The only accepted weight
vectors are `(0,0,0)` and `(3,0,3)`, and each layer is a singleton. Exactly
four internal phases change between them. Moving either interval alone
violates this guard; one budget-four coupled move succeeds. This guard is
an explicitly chosen integer condition, **not** the electrical mod-eight
guard.

## 6. Exact electrical guard, kernels and converse

The recovered original electrical phase has the form

\[
          \phi_x(T,P)=x_{[T]}-x_{[K_P(T)]}\pmod8,
\]

with a specified actual merged-metagraph potential `x` in the corresponding
Laplacian kernel. Its constant-on-automorphism-orbits result is supplied by
the separate guard packet; it is not a consequence of the number `q` or
the difference-word weight.

**An original electrical-guard example with actual fixed ports.** The
guard packet recovered a source of the original order-eight atlas: a source
singleton, then two directed C3 cells with cross differences `{1,2}`, then
a sink singleton. Pin the actual endpoint owners at clocks zero and seven.
The six legal contiguous-cell paths have phases `(a,b)` with
`b-a=0 or 1 (mod3)`. Using the original frozen electrical potential, the
complete table is

| Difference | Native merged source | Native merged target | Forward tile mask | Original phase mod8 |
|---|---:|---:|---:|---:|
| 0 | 3517 | 532 | 1765183 | 1 |
| 1 | 3517 | 1290 | 1799615 | 5 |

Difference two fails the actual boundary arc. Thus the unfiltered six-state
one-cell graph is connected, but either frozen phase level has three
isolated vertices under one-cell moves. The actual simultaneous move
`(a,b)->(a+1,b+1)` preserves the source, both indexed endpoint owners, the
target class and the original electrical phase. It connects all three
retained lifts in that level. Its exact minimum phase budget is two: every
distinct pair in a level differs in both coordinates.

This is a fully authenticated positive guarded move, not a synthetic
weight predicate. It stays within a common-phase automorphism orbit, so it
does not establish a move between different target classes. The source is
not regular, unlike the unbounded family of sections 1--5. The checker
independently rebuilds its actual arcs and uses the original owner map and
potential bytes, retaining their hashes and their **forward-bit** gauge.
Converse and the two actual endpoint pins are checked literally.

In our hosts, adding a common phase to every cell is an actual host
automorphism because all cross relations depend on phase differences.
Thus an isomorphism-invariant electrical phase can factor through the
complete difference word `d` after the common-phase lift is recorded.
It need not factor through `k(d)`, or through its weight vector after pins.
Sections 4 and 5 apply to this original electrical guard **only after that
stronger factorization is authenticated**. Otherwise the actual guard
must be evaluated on every proposed endpoint, and the accepted induced
move graph may still be disconnected. The checker does not pretend to
supply an order15-or-larger frozen electrical atlas.

Our nontrivial block moves generally change the native backward word, as
the order27 control explicitly verifies. They are not merely changes of
an invisible automorphism-kernel lift. The original tournament remains
fixed, every cell owner is retained, and the template clocks remain the
literal indices `iq+u` throughout. No unobserved common phase is selected.

Converse sends each retained actual path to its reversal in `T^op`, and
an owner clock `j` to `qt-1-j`; an arc clock also reverses its endpoints.
Conjugating our moves by this bijection preserves their phase-support
budget, accepted paths and transported pins. The factorization and guard
certificate must be transported too. This is a statement about actual
reversed paths, not about keeping the same printed clock or residue.

## 7. Connection ledger and reproduction

| Item | Retained structure |
|---|---|
| Source | Actual regular tournament, contiguous cyclic-cut observer family, indexed owner/arc pins |
| Target | Binary difference paths with integer weights, plus the actual common-phase lift |
| Map | `a -> (a_0, a_(i+1)-a_i)`; retain the unique binary lift of each legal difference |
| Preserved | Fixed native tournament, every owner and clock, accepted endpoint arcs, shared pins |
| Information lost by weights alone | Arrangement of differences, actual path, potential electrical phase |
| Needed guard sidecar | Original predicate evaluated on the actual path, or an authenticated weight-factorization table |
| Minimal hostile | Same endpoint residue but different integer weights; no one-cell path crosses them |
| Positive repair | Sharp coupled seam flips and the weighted quotient criterion (7) |

`check.py` constructs the actual regular hosts for `q=3,5,7` and
`r=2,4,6,8`, including every binary difference word and every actual
common-phase lift. It verifies degree balance, all literal Hamiltonian
junctions, primitive support, deterministic sorting, exact cross-layer
distance minima, and necessity/sufficiency for every nonempty weight
filter in each tested endpoint fibre. A composite `q=9,r=10` control checks
the residue-zero fibre on the actual order99 host. There are 13 host
families and 103 explicitly defined weight filters, plus the correlated
order27 pin hostile and actual converse control.

```text
python -B 05-knowledge/results/native_guard_moves_farey_20261023/moves/check.py
python -O -B 05-knowledge/results/native_guard_moves_farey_20261023/moves/check.py
```

The initial write run passed 35,898 exact checks. The receipt replay adds
one equality check. Adding the independently realized original order-eight
guarded example brings the write count to **35,949**. No floating
calculation or optimization-sensitive assertion is used. These
computations authenticate the stated small electrical example, but do not
authenticate a weight factorization on the larger regular hosts.
