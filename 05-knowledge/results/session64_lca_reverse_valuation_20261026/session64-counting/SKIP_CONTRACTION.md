# Direct skip-normalized ordinary tournament DAG

2026-10-26. **PROVED / IMPLEMENTED / FINITE-EXACT:** contract deterministic
unit-cost-free unaffordable-stage and early-terminal skip chains in the
ordinary `(ell,S,D)` DAG. `contracted_ordinary.py` constructs the resulting
smaller graph directly, without first materializing the full graph, and then
uses the frozen reverse-potential engine on it. This is an elementary
administrative compression already implicit in the recurrence's skip rule,
not a new defect-state equivalence, a new Burnside identity, or a claim to
historical priority.

## Exact contraction, not arbitrary stage erasure

A nonterminal state is *active* when ell<=S. If ell>S, the only possible
batch is k=0. The affordable D coordinates do not change when the stage
decreases from ell to ell-2: the stage has not yet passed an affordable odd
length. Thus b=b', u'=u, T=1 and J=0. An early S=0 state also has a unique
unit zero-cost suffix to the common terminal.

Define the canonical image

    c(ell,0,()) = terminal,
    c(ell,S,D)  = (min(ell,largest_odd_at_most(S)), S, D), S>0.

This is exactly the endpoint of the forced skip chain. At every active
vertex, retain every ordinary choice k and its original integer transfer
factor. Canonicalize only that edge's destination. The complete weighted
path sets before and after contraction are in a bijection: delete forced
unit edges, or uniquely reinsert them. Hence all exact counts, normalized
counts, and stage-compatible future choices are preserved.

The direct builder applies c immediately to each new destination, retaining
all affordable D coordinates. It never constructs the removed vertices. It
still keeps the active stage; this is not a quotient that forgets stage
indiscriminately. Since canonical target stage is strictly below source
stage, the graph remains acyclic.

This operation changes neither signed-core semantics nor literal clone
mergers: it is implemented only on the ordinary e=0 recurrence. All ordinary
tournament classes remain present, including classes containing clones.
Ports and owner identities were already absent in this counting carrier.

## Compatibility with reverse potentials and merged weights

The suffix minimum V is constant along a forced zero-cost chain. Consequently
it descends through c, every retained edge has the same reduced cost J', and
the baseline B is unchanged. The contracted graph's reverse pass and its
normalized modular transfer are therefore exactly the previous theorem on
a smaller finite graph.

Prefix minimum P and accumulated Y/Z need a more careful statement. Several
incoming histories can merge at the canonical destination, so a skipped
source vertex's partial aggregate need not equal its destination's full
aggregate. The code checks P and exact Y/Z only at surviving active vertices
and the common terminal. It checks V at every original vertex. At N=100,
440 removed vertices have prefix minimum strictly above that of their
canonical destination; this is a positive hostile to the false pointwise
P-descends claim.

## Cheap hostiles

* In the actual N=6 graph, `(3,1,(0,)) -> (1,1,(0,))` is a forced unit skip.
  Blindly dropping stage makes it a self-loop. The specified contraction is
  valid, because it deletes the forced edge and retains the active stage.
* Not every k=0 edge is contractible. In the actual N=8 graph the active
  vertex `(3,3,(0,2))` has two choices, and its k=0 edge has J=3 as the
  lower-bound gauge changes from b=1 to b'=0. Both that branch and its
  competing k=1 branch are retained.
* V descends along every forced chain; P and merged weights need not do so
  pointwise. Equality holds at the surviving full aggregate, not at each
  source contribution separately.

## Controls and measured structural reduction

    python contracted_ordinary.py > contracted-normal.json
    python -O contracted_ordinary.py > contracted-optimized.json
    cmp contracted-normal.json contracted-optimized.json

Both runs pass identically with explicit exception gates. All 54 ordinary
targets N=0..50,60,80,100 agree exactly with the frozen full-DAG engine;
432 normalized modular comparisons cover K=1..8. All 24,102 surviving
prefix coordinates agree in both exact Y and exact Z. Every omitted vertex
is explicitly verified to have one k=0 edge with T=1 and J=0; every retained
edge is checked against the full graph with its destination canonicalized.
The full engine separately supplies known-count, unmerged-partition, and
unrestricted full class-enumeration controls.

| N | Full input states | Direct contracted | Full edges | Direct contracted | K=2 updates, full→contracted | K=8 updates, full→contracted |
|---:|---:|---:|---:|---:|---:|---:|
|30|212|138|472|398|44→14|143→76|
|50|926|648|2732|2454|158→31|440→210|
|60|1720|1237|5499|5016|156→30|704→335|
|80|5119|3876|18109|16866|374→51|1408→650|
|100|13444|10547|50044|47147|568→64|2576→1161|

Input-state columns exclude the common terminal. Update columns use the
separate proved valuation-budget filter from `REVERSE_IMPLEMENTATION.md`.
Unlike that filter alone, skip contraction avoids constructing the removed
vertices and edges. It does not remove any active defect distinctions;
the existing minimality theorem for individual normalized continuations is
untouched. There is no measured wall-clock, peak-memory, or asymptotic
bit-complexity claim. Count arithmetic and graph bookkeeping are reported
separately, and the constructor still processes every retained edge.
