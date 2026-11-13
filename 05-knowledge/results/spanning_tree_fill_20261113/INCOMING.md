# Incoming work changes the allowed operation, not the truth boundary

**PROVED scoped comparisons / inherited results credited.** Initial commit
dfaf7c7e88dc; incoming798965f7a5 was read during this session. Neither
ordinary Collatz nor Chan–Pak is proved here.

## Three edge alphabets must stay distinct

The concurrent [fill-tree packet](../fill_tree_matrix_20261113/README.md)
studies the literal graph L: ordinary Collatz, even doubling, and the extra
odd n≡1 mod6 move to(n−1)/3. Its [forced-excursion proof](../fill_tree_matrix_20261113/collatz/NOTE.md)
shows that, for even K≥2, every successful L path from n=2^K−1 reaches
at least2·3^K−2. Universal L reachability was left OPEN.

The user's clarification in this session allows **every genuine reverse
Collatz edge too**. Our graph B contains L and additionally permits odd
doubling and genuine odd predecessors from even inputs≡4 mod6. Its
[complete descent](collatz/PROOF.md) has peak≤12n+4. For the same family
n=2^K−1,3 divides n, so the seven-edge macro sends n directly to n/3.
The step12n+4→4n+1 is a genuine reverse odd edge unavailable in L.
It is precisely an exit from a dyadic chain that invalidates applying
L's doubling-elimination lemma in B. There is no conflicting peak theorem.

Graph A, ordinary Collatz with only its genuine edges made reversible,
still has exactly the classical basin question. The extra odd-to-even
quotient in B is not a genuine inverse edge and remains an explicit
assumption. The model comparison therefore identifies which operation
closes this relaxed reachability problem.

## A simplex can be filled in different exact ways

The incoming [multiplicity proof](../fill_tree_matrix_20261113/multiplicity/PROOF.md)
weights insertion into gap i by g_i+1 and proves exact uniform distribution
over every fixed-size padding simplex, with d! labelled histories per cell.
Our [canonical forest](multiplicity/PROOF.md) instead chooses a unique parent
by deleting from the first nonempty gap. The former retains all physical
birth histories; the latter removes their quotient collisions. Both preserve
the same core word and Hamilton count. Neither creates a missing count.

The incoming rank-r difference formula recovers nonnegative primitive core
cells; forgetting r yields a negative coefficient already at H9. Our
permanent-port formula retains a different lost coordinate: whether each
named port lies in a core or in a singleton gap. Combined, the two exact
views identify both the cell dimension and its physical interface.

## The same target tree can lack a physical edge

The incoming all-order induced-count boundary, built from the earlier
November2 obstruction, excludes13→27 under any sequence of pure births.
Our [four-rule target tree](TARGET_TREE.md) also contains this edge. This
rules out an induced-only realization of that tree, including arbitrary
same-H13 source preparation followed by births. In contrast the literal
H3 preparation in [BIRTH_FLOW](multiplicity/BIRTH_FLOW.md) succeeds for
17/29 because it changes the old arcs before a cut that is physically
available. It is a scoped constructor, not a way around the universal
13→27 containment obstruction.

The concurrent tail-marked Matrix–Tree and square-zero path formulas agree
with our [matrix treatment](matrix/PROOF.md). These elementary applications
are independently recovered, not asserted as exclusive new results.
Our tail/head hostile, ordinary Boolean cofactor extraction, and exact
full birth-flow comparison expose complementary retained coordinates.

## Other incoming source mechanisms

At startup the [surplus-selector packet](../surplus_selector_20261113/README.md)
already supplied arbitrary23-prefix repair, compensated third cores and
actual neutral-context slope transport. Its exact selector is complete
for its declared language and finite head, not every odd target. The
[growing-t packet](../growing_t_returns_20261113/README.md) had integrated
the complete k17 unit-ready exclusion. These results reinforce the need
for changing interfaces: increasing the number of count-preserving
contexts does not produce a unit-ready source at an all-order ready hole.

During final integration,68ab145445 supplied the
[endpoint-context grammar](../endpoint_context_recursion_20261113/PORT_GRAMMAR.md):
changing the output port closes three fixed-port obstructions and gives a
renewable four-vertex context with target recurrence X_(k+2)=8X_(k+1)+X_k.
Its all-order four-coordinate compiler and our retained-singleton handoff
both require the physical output port to be part of the state. Neither
supplies a universal absolute-target selector or removes the ready hole.

The clock continuation through8961e87e88 adds a proved two-cherry I/O
factorization and declared finite order-twelve search results. Its wrap
and anchor patterns remain hypotheses. During integration, unsupported
global spectrum and all-clock composition wording was narrowed to its
actual scope and recorded in [MISTAKES](../../../01-canon/MISTAKES.md).
No numerical result or clock construction was changed.
