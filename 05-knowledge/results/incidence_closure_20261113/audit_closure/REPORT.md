# Independent review of closure and fixed-anchor minimization

**Audit verdict: PASS, with explicit scope.** The elementary coarsest-partition
proof, source normalization, fixed-anchor reduction, dominance and branch/prune
logic are valid for the declared universe. No mathematical or implementation
repair is requested. This report is independent review evidence; the producer's
RESERVED status is not changed here. The root owns any promotion.

The finite minimum sequence is `0,0,1,2,9` for orders `3,4,5,6,7`, in each of
four palettes: raw uncolored, source-normalized uncolored, raw blue/black and
source-normalized blue/black. The minimum is the number of **named physical
anchor columns whose direct measurements, together with full Q, are stable**.
A colored anchor supplies two numeric outputs at the cost of one named anchor.
This does not optimize arbitrary scalar encodings, bit counts at all orders,
linear Krylov dimension, or the number of anchor seeds before further refinement.

## 1. Frozen source boundary

All paths below are relative to the enclosing incidence-closure packet. The
following exact versions were read and checked:

| File | SHA256 |
|---|---|
| `closure/PROOF.md` | `960634c39d54650ca1bc5ff44b56a2de1b3ef66e9e7bda1cd5604f14f78d2a53` |
| `closure/analyze.py` | `8a06553cdd048e5345188311f9f026317cf7fc24d5f6dace64e94bf48d1a611c` |
| `closure/check.py` | `4638279a9cd5fc6c121a0dff84fbf500c029ac4915e0db78df47a4c7ff434747` |
| `closure/orbits.cpp` | `afbe0cee83c798cac153b1950e5cc0012d4339acf7f766560f1e8188cac68f12` |
| `closure/data.json` | `c0484a99712edff88e8158b5f6000097f32f88874150d059da5aefa20f3de8f7` |
| `closure/analysis.json` | `41093c9b0b11e055ceaafcd3b6e247da1077cbe978a7addd9661fe1300e65023` |
| `audit/independent.py` | `793dc9f3090387e826446d8a749ca80870ee5d6d3a148b5597162f791a2b4831` |

The first three files and the relevant sections of `audit/independent.py` were
read mathematically. The orbit helper is pinned as an input dependency, **not**
claimed as separately audited source by this lane. This review does not repeat
the complete order-seven labelled tournament census or independently exhaust
its eight-anchor search. It reviews that search's completeness and checks its
recorded failed budgets, witnesses and lower-bound packing. The root's different
cover solver is an additional independent verification, not presumed completed
by this report.

## 2. The scalar weights and collision universe are correctly typed

The merged matrix is an integer matrix of actual fixed-spine darts; each loop
has two diagonal darts. Its row sum is

    d(T) = |{[T],[T^op]}| H(T) / |Aut(T)|.

The ordinary matrix lacks the converse-orbit factor. This distinction is
retained in the data and in the separate ordinary audit. The normalized walk
is `P(i,j)=A(i,j)/d(i)`. Both color subkernels divide by this same total source
degree. Color-specific row denominators would be a different process and are
not used here.

In `analyze.py`, normalized comparisons use the exact cross product

    M(i,c) d(j) != M(j,c) d(i),

which is precisely inequality of source-normalized entries. The independent
certificate checker rebuilds actual row sums, verifies them against exported
degrees, and uses rational fractions rather than cross products. Every degree
is positive in the supplied tournament universe.

The optimization universe is the unordered distinct pairs in the same **full**
Q cell. The audit reconstructs this universe using ordered-cover vectors,
rather than the producer monomial-coefficient vectors. The change of basis
between full ordered-cover coefficients and full Q coefficients is invertible;
the two collision universes agree exactly. At seven there are 184 such pairs.
No H-only, score-only or sampled-Q substitution occurs.

## 3. Coarsest stable partition: proof audit

Let a partition be refined by retaining the current source cell and appending
all row sums into every current target cell for every specified operator.
The retention of the source cell is necessary and is present in both producer
and root implementations. Thus each iteration only refines.

Every strict iteration increases the number of cells, giving the stated bound
`|V|−|Pi_0|`. At termination, all required cell sums agree within each cell.
For minimality, suppose Theta is any stable refinement of Pi_0. If Theta refines
Pi_r, every Pi_r target cell is a union of Theta target cells. Each summand has
constant row sum on a Theta source cell, so the signature defining Pi_(r+1)
is constant there. Hence Theta refines Pi_(r+1). Induction proves Theta refines
the fixed partition. This is exactly the **unique coarsest stable refinement**,
not merely a refinement that happens to stabilize.

The intertwining `ML=LB` for cell-indicator lift L follows directly from the
stable row sums. Its products preserve cell-projected walk responses; for a
palette, the corresponding identity holds for every colored word. For P this
is every-starting-state strong lumpability, not only stationary aggregation.
The A/P normalization rule is also correct: constant positive degree in each
initial cell remains constant in every later cell, so multiplication of a
source signature by that cell's common reciprocal degree preserves equality
at every iteration.

The pointwise-algebra reformulation in section 7 is valid over the real or
rational fields used here. A finite unital subalgebra of functions is the whole
space of functions constant on its indistinguishability partition: polynomial
interpolation of separating functions recovers the cell indicators. Operator
invariance gives precisely the stable-row condition. A merely linear generated
space need not be a pointwise algebra, so no Krylov-rank conclusion is smuggled
into the partition counts.

As a cheap independent hostile/control family, the audit enumerates every simple
undirected labelled graph on 1..4 vertices, every initial partition, and raw or
two-color edge palettes. Across 2,010 cases it computes closure and independently
compares with **all** stable refinements. The coarsest universal property and the
strict-round bound hold. These cases support the implementation; the induction
above is the all-finite proof.

## 4. Exact anchor reduction and pruning

An anchor covers a collision pair exactly when its output differs on that pair.
A set of anchors separates every Q cell iff its masks cover the collision
universe. At the audited orders, the coarsest stable Q refinement is discrete,
so a Q-refining measurement partition is stable iff discrete. This supplies the
claimed equivalence between direct separation and stability; it would require
rechecking at a new order.

Each pruning operation in `exact_cover` is sound:

1. **Duplicate masks.** Equal coverage has equal cost; retain any representative.
2. **Dominance.** Replace any strict-subset mask by a containing mask. The covered
   set cannot shrink and the number of selected columns cannot increase, even if
   the replacement was already selected. A maximal superset exists in the finite
   family. Thus discarding all nonmaximal masks preserves the optimum.
3. **Packing.** Each chosen uncovered pair must be covered by a column in its
   candidate neighborhood. Pairwise disjoint neighborhoods force distinct columns.
   The greedy packing need not be maximum; its size is still a lower bound. The
   audit validates the supplied packing using **all original columns**, independently
   of the producer's dominance reduction.
4. **Maximum coverage.** If no column covers more than b remaining pairs, k columns
   cover at most kb, including overlaps. The test `k*b < |U|` safely rejects k.
5. **Branching.** Every cover must contain some column covering the selected
   uncovered pair. Trying all such candidates is exhaustive. Their order is only
   a heuristic. Removing the chosen column's coverage gives precisely the remaining
   obligation. Previously used columns cannot cover a currently uncovered pair,
   so no omitted selected-column state is needed in the memo key `(U,k)`.
6. **Budgets.** The greedy witness supplies a finite upper bound; the packing
   supplies a valid lower bound. Exhausting every intermediate budget proves the
   first successful budget minimal. Lower bound five at seven excludes budgets
   zero through four, and the stored failures explicitly exhaust five through eight.

The producer failed-budget-eight state counts are 45,300 (A), 45,300 (colored A),
44,905 (P), and 44,415 (colored P). These are memoized search-state counts, not
counts of all eight-element subsets. Exhaustiveness comes from the covering-pair
branch argument and safe pruning, not from their magnitude.

The audit loads the exact producer `exact_cover` function **only for a differential
algorithm test**; it does not pretend this is an independently implemented solver.
Its reference routine enumerates every subset of columns. They agree on all
4,944 distinct families of at most five nonzero four-bit masks, plus ten controls
with zeros, duplicates, dominance, reordered columns and noncontiguous bit indices.
The concrete order-three through order-six minima are independently checked by
literal column-subset enumeration as well. At seven this lane directly verifies
the nine-column upper witnesses and five-pair packing, and audits the exhaustive
recursion. The root's separate solver owns any additional independent lower-bound
computation at seven.

## 5. Root independent incidence and arbitrary-depth repair

The pinned `audit/independent.py` does not import producer census code. Its
`literal` routine builds ordinary relabeling orbits by applying every permutation
to the least remaining labelled orientation; deleting complete orbits justifies
its no-overlap check. Aut is then `n!/orbit_size`, converse is full bit complement,
and ordinary representatives use unrestricted minima.

For each ordinary representative, every Hamiltonian permutation is enumerated.
The native target flips exactly the nonpath bits (`top XOR pathmask`). Dividing
physical path counts by Aut is legitimate because an automorphism fixing an
ordered Hamiltonian path fixes every vertex; the action is free. Blue detection
constructs the actual path-reversal permutation and tests whether it carries the
tournament to its converse. These are different formulas from the fixed-spine
producer construction and compute the intended same incidence. Ordered covers
are counted by forcing breaks at every backward adjacency and choosing any of
the remaining breaks; the binomial formula in the root audit is correct.

The symmetric ordinary matrix, converse equivariance, doubled loops, and the
projection obtained by summing all ordinary entries in two converse fibres are
checked explicitly. Chiral signed-block reconstruction also has the correct
factor: the merged entry is twice the ordinary same/opposite sum, while their
difference is the odd block. This review accepts the logic of the root's literal
n=3..6 audit; it does not enlarge that independent literal universe to n=7.

For the separate reversible path `0—1—...—m`, start with cells `{0}` and
`{1,...,m}` and use the degree-normalized walk. In each round the next vertex
adjacent to an already isolated vertex is split off. Every other vertex in the
remaining interval still sends probability one to that interval, so it does not
split early. Hence exactly `m−1` strict rounds are required. This directly
justifies the root finite checks for m=2..14 and the arbitrary-depth example.

At the initial partition, actual two-step return probability to 0 is `1/2`.
The stationary coarse two-state walk returns with probability `1/(2m−1)`, since
the total degree mass of the nonzero cell is `2m−1` and its boundary flow is one.
Their strictly positive difference for m>=2 gives a second-step detection while
repair still takes `m−1` rounds. The n=7 tournament trace `156→271→272` is therefore
an instance of two repair rounds, not a universal repair bound. The sheaf plus++
theorem does not apply to this operation without a site and its restriction maps.

The reviewed root function hashes, useful if its main gains further controls,
are:

- `literal`: `fd7d88b40272b90960cb6be021abccfc18ee7c4e613e0aa50cc9cce25035a3bc`
- `refine`: `9eb421ffeda8500c1dc9712856f1759ad5210bf24727f21f1fcf2f03257d2b07`
- `initial`: `51076567717e5bc6c864bb9fe6cdf7acba43ca271ecd5f475e2072fd674a6c6d`

## 6. Reproduction and limits

From the isolated worktree root:

```sh
python3 05-knowledge/results/incidence_closure_20261113/audit_closure/check.py > 05-knowledge/results/incidence_closure_20261113/audit_closure/normal.json
python3 -O 05-knowledge/results/incidence_closure_20261113/audit_closure/check.py > 05-knowledge/results/incidence_closure_20261113/audit_closure/optimized.json
```

Both receipts are byte-identical and record 16,097 explicit checks. The checker
uses exceptions instead of removable assertions, writes no producer files, and
pins the producer source and data. Its SHA256 is
`67091a89429d1972c669b793cd9b53a9052285fc655edc8d6ceb30c69330b5d7`.
The matching receipt SHA256 is
`21824119e7d9d46be31f2ad89c197f2f514fac7dcc985c01e8ac73cdf1208569`.

No all-order tournament separation, order-eight/nine closure, graph-generation
recurrence, or spectrum conclusion is obtained. The proof review is general
for the finite partition algorithm; the displayed tournament minima are finite.
No producer file or status was changed and no commit was made by this lane.
