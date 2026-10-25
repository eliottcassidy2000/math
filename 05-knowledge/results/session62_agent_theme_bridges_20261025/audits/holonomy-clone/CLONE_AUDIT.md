# Independent audit: exact clone cancellation and real exterior invariants

2026-10-25. **PASS / PROVED in the stated scope.** The exact root order,
legal fixed-XOR operations, exterior generalization, unbounded tournament
realization and synchronization hostile all survive independent audit.
No historical novelty or computational-complexity claim is certified.

## Inheritance and scope

The reviewed proof is the clone lane's `README.md`, with `check.py`.
The weighted inflation fiber/cut-orbit correspondence is inherited from
`05-knowledge/results/ports_regular_actions_hn_20261025/counting/CLONE_INVOLUTION.md`.
Its module-path input is **THM-4463**, file
`01-canon/theorems/THM-4463-homogeneous-tournament-contacts-and-pinned-scale-floor.md`.
These inputs were inspected. The new finite tests do not constitute a
fresh factorial enumeration of all weighted quotient graphs at order 33.

For a fixed tournament with canonical transitive modules, the gap set is
intrinsic, automorphisms preserve directed module orders, and cut subsets
parametrize interval decompositions. A weight-preserving quotient
isomorphism lifts through the unique transitive-block order isomorphisms.
Thus cut orbits, rather than labeled cuts or module counts alone, are the
correct inherited fiber carrier. Tournament automorphism groups have odd
order: an involution would reverse the edge between an exchanged pair,
and every finite group of even order contains an involution by inverse
pairing.

## Independent root-order derivation

Let an odd-order finite group `A` act on `G`, and let `P(t)` count its
subset orbits by size. Applying Burnside separately at each cardinality
gives a summand `product_D(1+t^|D|)` for each group element. Its cycles
have odd lengths. At `t=-1+x`, every factor begins with `|D|*x`.
Therefore the first nonzero coefficient of a minimum-order summand is
positive. Such summands cannot cancel, proving

```
d = ord_(-1) P = min_g c(g),
kappa = (1/|A|) sum_(c(g)=d) product_D |D| > 0.
```

This is independent of faithfulness. The empty-set convention gives
`P=1,d=0,kappa=1`. Since `P` has integer coefficients, both exact division
by `(1+t)^d` and its leading coefficient are integral, hence `2^d` divides
`P(1)`.

If `r` is the number of whole-group point orbits, each element has at
least one cycle per orbit. Equality `d=r` requires one *single* element
that is cyclic on every orbit simultaneously. Odd cycle lengths imply
`d=r=|G| mod 2`. For transitive nontrivial actions, the fixed-point average
is one, so some element is a derangement; its odd cycles have length at
least three and give `d<=|G|/3`. This argument does not supply simultaneous
derangements on distinct orbits.

The signed binomial moments in the actual clone fiber follow by expanding
`P(-1+x)`. Their conversion factor is `(-1)^(n+j)` and equals `+1` for
`j=d`, because `n=d mod 2`. Thus all moments below `d` vanish and the first
nonzero one is exactly positive `kappa`, not its negative.

## The operation sidecar and an audit strengthening

For a fixed subset `F`, the map `C -> C symmetric_difference F` descends
to the subset-orbit quotient **iff `F` is invariant**. Necessity does not
need odd group order: `F` and `gF` lie in the same input orbit, while their
images are empty and `F symmetric_difference gF`. The empty orbit contains
only the empty set, forcing `F=gF`. Sufficiency is equivariance. This
stronger descent statement was supplied during audit and incorporated in
the author's proof.

There are exactly `2^r` invariant masks, giving the full rank-`r` fixed-XOR
group. For odd `A` it acts freely on subset orbits. Indeed, if
`gC=C symmetric_difference F`, iterate through the odd order of `g`,
using `gF=F`, to conclude `F` is empty.

For orbit sizes `n_i` and intersection counts `a_i`, a Boolean-operation
orbit has polynomial `product_i(t^a_i+t^(n_i-a_i))`. The derivative of
one factor at `-1` is `(-1)^a_i(n_i-2a_i)`, nonzero because `n_i` is odd.
Each such Boolean orbit therefore has order exactly `r`. Extra order
`d>r` is cancellation **between** these operation orbits, not extra
fixed-XOR generators. Arbitrary nonlinear permutations of the orbit set
are outside this classification.

## Real exterior-algebra theorem

For a finite group acting on a real finite-dimensional vector space,
averaging the trace on each exterior power is projection onto invariants.
Expanding `det(I+t g)` as exterior traces proves the displayed exterior
Molien identity directly. An invariant inner product makes each matrix
orthogonal. At `t=-1`, each fixed eigenvalue contributes one zero; every
remaining real eigenvalue is `-1` and contributes `2`, while nonreal
conjugate pairs contribute `|1-lambda|^2>0`. Consequently

```
ord_(-1) sum_j dim((exterior^j V)^A)t^j = min_g dim V^g.
```

Odd group order is **not** needed for this real-representation theorem.
Splitting off the common invariant space gives a genuine exterior factor
`(1+t)^dim(V^A)` with nonnegative integral residual coefficients.

For an odd permutation group, every subset stabilizer preserves the
orientation of its wedge: all its cycles are odd. One invariant line
therefore occurs per subset orbit. This proves equality with the subset
polynomial as graded dimensions, not as a multiplication rule on subsets.

Both advertised hostiles are correct:

- A transposition on two points has subset polynomial `1+t+t^2` and root
  order zero; its real exterior invariant polynomial is `1+t`, order one.
- Scalar multiplication by a primitive cube root on complex `C^3` has
  minimum fixed dimension zero but exterior polynomial `1+t^3`, order one.
  Real conjugate-pair positivity cannot simply be deleted.

The standard identity is externally cited in the candidate; the argument
above makes the audit independent of that bibliographic attribution.

## Literal tournament realization

The core `R_k` has data vertices `F3^k` with lexicographic cyclic orientation
and ordered observer triangles `Y_i`, with `x_a -> y_(i,c)` iff `a_i=c`.
Its stated degree formulas are exact. The smallest observer degree minus
the data degree is `(3^(k-1)+3-2k)/2`, positive for every `k>=2`.
All parts are therefore intrinsic. Observer automorphisms are rotations;
their shifts force the entire data map to be the corresponding translation.
All translations preserve the arcs, proving the **full** automorphism group
is `C3^k`, not merely a subgroup.

Pair-freeness is also justified: differing data coordinates have an
observer separator; two observers in one triangle have a data separator;
different observer parts have degree difference at least three, whereas
homogeneous pairs have degree difference exactly one; a data/observer pair
is distinguished by a different observer triangle. The last step requires
`k>=2`.

Inflating all data vertices to ordered pairs introduces exactly those
homogeneous pairs. Any pair crossing base blocks retains a base separator.
The homogeneous-pair paths therefore have precisely the intended two-point
components and singleton observers. Automorphisms descend to size-preserving
core automorphisms and lift uniquely; every translation preserves sizes.
The gap action is regular `C3^k`. Nonidentity translations have only
three-cycles, yielding `r=1,d=3^(k-1)` and the stated positive coefficient.
The base `R_1` really has three unexpected pairs, so omitting `k>=2` would
invalidate this proof.

The four-form construction at `k=2` is equally sound. Its degree classes
are `8` for data and `16,13,10,7` for the observers. Two coordinate rotations
force data translation, and the other two rotations are then determined.
Inflating the twelve observer vertices yields the advertised 33-vertex
tournament. Every nonzero translation belongs to exactly one of the four
one-dimensional kernels. It has three fixed gap points in that orbit and
one three-cycle in each other orbit: `r=4,d=6,kappa=24`. This is a literal
tournament hostile to replacing simultaneous cyclicity by separate choices.

## Independent computations and precise limits

The standalone `check_independent.py` imports no candidate module. It
constructs subset orbits by removing whole orbits, divides the resulting
integer polynomial repeatedly by `1+t`, and compares the first nonzero
coefficient with the cycle prediction. It also partitions subset orbits
into invariant-toggle cubes and verifies their order is exactly `r`.

Its action universe consists of the empty and trivial three-point action;
regular `C3`, `C9`, and `C3^2`; four projections of `C3^2`; the nonabelian
order-21 affine group on seven points; `C15` acting on disjoint cycles of
lengths three and five; diagonal and independent actions on two triples;
and the even transposition hostile. Each action is explicitly constructed.
It verifies the fixed-XOR necessity witness for every mask of every action.

The real-representation controls use independent exact determinant
expansion and rational rank computation for a sign representation,
a real `C3` rotation, the real `D3` standard representation, and that
rotation plus a trivial line. The complex hostile is calculated by its
exterior characters.

A different tournament validation path exhausts observer rotations after
independent degree and separator checks, then tests every forced map on
every arc. It verifies every homogeneous pair in the core and inflated
tournaments for `k=2,3,4` (24, 63 and 174 inflated vertices), and the
33-vertex four-form example. This is a structure-assisted exhaustive
argument, not a claim to have used an unrelated generic automorphism
package. The candidate's separate generic adjacency-backtracking checker
was replayed in both modes, including all automorphisms of its 24-, 33-
and 63-vertex examples. Candidate receipts matched exactly.

No `2^27` or `2^81` subset census is claimed. The larger examples' root
orders follow from their proved gap actions and the arbitrary-action theorem.
The combined independent suite has **123,237 explicit gates**, with
byte-identical ordinary/optimized outputs. Its clone evidence is under
`odd_actions`, `literal_tournaments` and `exterior`.

No faster enumerator, black-path move, connectivity theorem, geometric
realization, or prize-problem resolution follows from these results.
