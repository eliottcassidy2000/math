# Four-step tournament returns: a finite-local proof

**Status: PROVED IN NOTE / independently internally reviewed.** The finite
local lemmas have exact certificates in four_local.py and matching
normal/-O receipts. See [the independent review](AUDIT_FOUR.md). The
all-order step is justified below, not extrapolated from a small census.

**Claim.** Every sequence of at most four ordinary Fomin–Zelevinsky mutations
whose initial and final matrices are tournaments preserves the number H of
Hamiltonian paths. Hence any H-changing ambient return has length at least
five. This does not claim such a return exists at length five or any length.

## The local-to-global mechanism

Fix a word whose mutated vertex set is S. Mutation restricted to any
induced set containing S agrees with restriction of the full mutation.
In particular every initial, intermediate, or final entry (i,j) is determined
by the induced tournament on S∪{i,j}. This set has at most |S|+2 vertices.
Thus a claimed implication “final matrix is a tournament ⇒ intermediate
matrix is a tournament” can be certified for all orders by checking all
labelled tournaments up to |S|+2, with the word fixed. A counterexample
would have a bad intermediate entry, and restriction to its endpoints
plus S would be a counterexample in that finite universe.

We use this only for properties of entries. The H conclusion comes from
the explicit global operations below; local H equalities alone would not
prove global preservation.

Cancel consecutive identical mutations first. Every remaining four-letter
word, after renaming vertices in first-appearance order, has one of the
five shapes 0101, 0102, 0120, 0121, 0123. Shorter words were handled in
[PROOF.md, section 2](PROOF.md).

## The five cases

**0101.** The inherited rank-two pentagon says 01010 sends a tournament B
to τB, where τ swaps vertices 0 and 1. See
[Session291 cluster theorem 3](../../session291_cluster_ports_20261113/cluster/PROOF.md).
If D is the result of 0101, then μ_0(D)=τB and hence

    D=μ_0(τB)=τ(μ_1(B)).

Therefore D is a tournament iff vertex 1 is a source or sink of B. Its H
then agrees with B by native mutation and relabelling. This is an exact
all-order characterization, not just a finite test.

**0102.** Exhaustive direct substitution on every labelled tournament of
orders 3,4,5 proves: if the final matrix is a tournament, the matrix after
010 is a tournament. The set S has size three, so the entry-locality lemma
upgrades this finite implication to every order. The three-step return
preserves H by the previous theorem, and the last native step preserves H.

**0121.** The same complete orders-3-through-5 substitution proves: a
tournament output forces the first step to be native. Again entry locality
proves this for all orders. The remaining three-step return preserves H.

**0123.** On four vertices, a tournament output is possible only for the
two total orders 0→1→2→3 and its reverse. On five vertices, where vertex 4
is exterior to the word, only the total order 0→1→2→3→4 and its reverse
work. For an arbitrary larger input, restriction to the four named vertices
fixes one of these orientations; restriction together with each exterior
vertex places every exterior vertex after the four-block in the positive
case, or before it in the negative case. Thus the original tournament is

    TT_4→R or R→TT_4,

with the four mutations applied successively at sources or at sinks. All
four steps are native. No condition on internal edges of R is imposed.

**0120.** At order three the complete return table is:

| Initial edge mask | Final edge mask | H |
|---:|---:|---:|
|0|3|1|
|1|1|1|
|2|5|3|
|5|2|3|
|6|6|1|
|7|4|1|

The bit order is (01,02,12), with bit one denoting the forward edge.
At order four, a tournament output occurs exactly when vertices 0,1,2
form a directed triangle and the exterior vertex either dominates all
three or loses to all three. At order five, the same criterion holds and
the two exterior vertices must lie on the same side; their mutual edge is
unrestricted. The local certificates also show that the output reverses
all edges meeting the triangle and leaves the exterior mutual edge intact.

Consequently, for every order n≥4 the global input must be C_3→R or
R→C_3. The output moves that extreme triangle to the other side and
reverses its orientation, leaving R unchanged. Necessity follows by the
four-vertex restriction for each exterior vertex and the five-vertex
restriction for each exterior pair. Sufficiency follows from the same
local entry tests: they determine every output edge. Therefore both
counts equal 3H(R). Together with the explicit three-vertex table, this
settles the final word shape. All cases preserve H. ∎

## Bounded exploration beyond the proof

scout_returns.py independently enumerates genuine four-step excursions in
all labelled tournaments through order five. The counts by word shape are:

| n | 0101 | 0120 |
|---:|---:|---:|
|3|0|24|
|4|96|96|
|5|1920|480|

Every such excursion preserves H. It also explores all mutation words
through length seven from one representative of each of the 18 ordinary
types at orders 3,4,5. Across these seeds the balls contain 183,018 matrices
(sum, not a disjoint union); there is no cutoff on matrix weights and no
H-changing tournament return was found. Exact layer counts and maximum
entries are retained in scout_receipt.json. This is a bounded search,
not closure of the mutation classes or an all-length invariance theorem.

Stopping reason: the decisive short-word classification is complete, while
the next search window supplies no positive count-changing signal. The
actual construction still needs a physical operation beyond these native
steps, transitive-pair transfers, extreme-triangle transfers, and pentagon
relabelings.

Reproduce from the repository root:

```text
python -B 05-knowledge/results/triples_ports_sieve_20261112/ports/four_local.py
python -O -B 05-knowledge/results/triples_ports_sieve_20261112/ports/four_local.py four_local_receipt.optimized.json
python -B 05-knowledge/results/triples_ports_sieve_20261112/ports/scout_returns.py
python -O -B 05-knowledge/results/triples_ports_sieve_20261112/ports/scout_returns.py scout_receipt.optimized.json
```
