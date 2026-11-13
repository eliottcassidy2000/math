# Strong collisions for the entire first-moment slot hierarchy

Let Q_T(q) count ordered covers of T by q path slots (empty slots allowed).
For each slot keep, separately, every labelled undirected arc's total use
count and every labelled start and end total. Do this for all q. It still
does not identify a strong tournament up to converse, even after adding
the full labelled induced-subtournament slot-polynomial table.

## A converse replacement invisibility theorem

Suppose a module U in a tournament is vertex-transitive internally.
Replace U by its converse, without changing any outside incidence.
For each cover, first retain the exterior vertex orders in all slots and
the positions of the maximal runs in U. This is a quotient layout.
If there are j runs, fillings of those ordered nonempty runs are exactly
the F_j(U) path covers of U; no arcs between U and its exterior depend on
which vertices fill the run. Reversing each run gives a bijection to
fillings in the converse. Thus:

- The unmarked layout counts and exterior-to-exterior arc uses are unchanged.
- An internal unordered arc use in any slot is unchanged by this bijection.
- A contact to a particular vertex of U is unchanged: in any fixed
  nonempty run, each of the |U| vertices is its start in exactly
  F_j(U)/|U| fillings, and also its end in that many. This uses vertex
  transitivity, not a false claim that reversal preserves a named start.
  The same argument gives all per-slot endpoint marginals.
- Even the full labelled induced Q table is unchanged. Indeed on restricting
  to any subset, the intersection with U is still a module. For this claim
  transitivity is not needed: the run-reversal bijection gives the
  unmarked count for any intersection.

This also works for any sequence of such replacements. It does not
identify the two-end distribution of one run or the joint usage of
two adjacent arcs: reversing a run exchanges its start and end.

## A strong order-24 witness, and exponential fibres

Use the non-self-converse vertex-transitive 21-vertex Cayley tournament U
from [Session 319](../session319_arc_usage_and_signed_cycles_20261113/ARC_USAGE_SHADOW.md).
Take a directed triangle of blocks with sizes 21,1,2, putting U in the
first block and transitive tournaments in the other two. Reverse U alone.
Both outputs are strong. The three blocks are their maximal proper
modules: a module meeting two blocks must meet the third (it otherwise
has a distinguisher there), then any outside vertex would distinguish
vertices in the other two, so the module would be the full tournament.
Since all three sizes differ, an isomorphism fixes the block positions
and would identify U with its converse, impossible. A global converse
reverses the quotient triangle and no size-preserving permutation restores
it. This is a collision at order 24 for the entire first-moment hierarchy
above. No order minimality is claimed.

Replace the 21-block by an ordinal sequence of k copies of U. Each copy
is itself a module and can independently be reversed, invisibly. The
ordered SCC decomposition inside the large block preserves all k positions
under isomorphism, and the distinct-sized quotient blocks again rule out
global converse. At order 21k+3, one strong fibre therefore contains at
least 2^k distinct converse-merged classes, even with the induced Q table
and every per-slot first moment attached.

## A selective test that breaks the collision

Write the outer cycle as w→S→A→w, with A=(a→a') and S the large block.
For any pair u,v in S, simultaneous usage of the two undirected edges
{u,v},{v,a} in one Hamilton path is possible iff u→v. Indeed v→a
is fixed, so co-usage forces u→v→a. If the orientation is correct, begin
with u,v,a,a',w and follow any Hamilton path of S minus {u,v}, which
exists by Rédei. This is an exact two-contact read of the internal
direction, not a comparison of the scalar H values. It explains how the
higher observer sees the module chirality without any labelled arc
directions in its own input.
