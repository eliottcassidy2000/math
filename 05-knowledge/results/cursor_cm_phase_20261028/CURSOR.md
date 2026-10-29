# Read-only cursor: complete faithfulness through order nine

**FINITE-EXACT**, with an all-size search reduction proved below. All-size
faithfulness itself remains **OPEN**. This is local result C1 of this file,
not a promotion of an earlier reserved theorem.

## The object and the observation

A Hamilton presentation is a tournament T on positions 0,...,n-1 with every
arc i->i+1. Native complement tau reverses all other arcs. The path order is
part of the state: two different bit rows are different presentations.

Let R_i(T) be the isomorphism class of T with its first vertex, last vertex,
and position i marked with three roles (roles may coincide). The inherited
[read-only cursor](../moving_selector_atoms_20261028/THEORY.md), C1, observes

    ((R_i(T), R_i(tau T)))_(i=0,...,n-1).

Its published carrier also retains cursor position and a native blue flag.
We test the displayed smaller record. Faithfulness of that record implies
faithfulness after those extra observations are added. We do not identify
independently chosen isomorphisms at different positions with one common map.

## Exact group reduction, at every size

Write G=Aut(T;0,n-1). If G is trivial, the source scan alone is faithful:
choose an endpoint-preserving isomorphism to a competing presentation.
For each position, uniqueness of that isomorphism forces it to fix that
position. It must be the identity. The same argument applies to tau T.

More generally, equality of source scans implies that, after choosing one
endpoint isomorphism, a competing presentation is T relabelled by a
permutation preserving each G-orbit setwise. Conversely every such
relabelled presentation has the same source scan. To see the converse,
each displaced vertex can be restored by some element of G; the restoring
element need not be the same for different vertices.

Therefore a complete collision search need only:

1. Retain T for which both endpoint groups are nontrivial.
2. Enumerate permutations inside the source endpoint-orbit blocks.
3. Discard unchanged presentations and those without the prescribed path.
4. For each survivor U, compute all endpoint-preserving isomorphisms
   tau T -> tau U. For each i, require at least one of these maps to fix i.

This is a simultaneous existence test at every marked position, not a
comparison of group orders or coarse orbit-size multisets.

The two endpoint groups can BOTH be nontrivial, and their vertex orbits
can overlap nontrivially. Already at seven there are examples with group
orders 3 and 5. Thus neither 'one side must be rigid' nor 'the two orbit
partitions always meet discretely' proves the all-size conjecture.

## Two complete finite paths

[cursor_probe.py](cursor_probe.py) enumerates all
2^((n-1)(n-2)/2) path presentations directly through eight. All endpoint
automorphisms and transporters use the inherited exact backtracking oracle,
without a probabilistic isomorphism hash.

| n | All path presentations | Both endpoint groups nontrivial | Distinct legal candidates tested | Collisions |
|---|---:|---:|---:|---:|
|1..6|1, 1, 2, 8, 64, 1024|0|0|0|
|7|32768|4|20|0|
|8|2097152|32|160|0|
|9|268435456|580|2256|0|

The last row uses [cursor_symmetry.py](cursor_symmetry.py). A tournament
automorphism group has odd order: an involution would swap some pair and
reverse its arc. By Cauchy's theorem, a nontrivial endpoint stabilizer has
an element of odd prime order. At nine its action is on seven interior
vertices, so that prime is 3, 5, or 7.

Enumerate every such permutation, including products of disjoint equal
prime cycles. Its action partitions unordered pairs, with the orientation
sign tracked. Prescribed path arcs either force a consistent orientation
of a pair orbit or reject the permutation. Enumerate every free pair-orbit
choice, and deduplicate the resulting tournaments. This covers every
possible nontrivial source endpoint group; all omitted sources were already
settled by the rigidity argument.

| n | Prime-order permutations | Sources with nontrivial endpoint group |
|---|---:|---:|
|7|44|262|
|8|224|7448|
|9|1574|386054|

The two methods agree on the dual-nontrivial counts and tested candidates
at seven and eight. They share the GI oracle; this is an alternate coverage
path, not an independent implementation of graph isomorphism.

## What is and is not settled

The record is faithful through nine, even without the extra blue flag.
That is an exhaustive result over the stated finite universe. The Cauchy
reduction is valid at all sizes, but its finite execution is not an all-size
proof. The next useful obstruction is a pair of endpoint-isomorphism cosets
whose pointwise fixing requirements all hold without a common identity map.

The [endpoint response theorem](ENDPOINTS.md) concerns a different menu:
one-stage unlabelled source/partner classes after intact-run assembly. Its
eight-vertex collision does not refute this cursor census: the cursor can
inspect internal positions that the endpoint response deliberately forgets.

Reproduction (from repository root):

```text
python -B 05-knowledge/results/cursor_cm_phase_20261028/cursor_probe.py --max-order 8
python -B 05-knowledge/results/cursor_cm_phase_20261028/cursor_symmetry.py
```

Receipts: [direct census](cursor_probe.json), [symmetry census](cursor_symmetry.json).
