# Supplement — Paley11 bounded depth search and exact compiler lower bounds

**Status:** PROVED necessary clock bounds; FINITE-EXACT bounded search with
NO WITNESS; Paley11 two-complement forest reachability remains OPEN.
Frozen main packet and scaling addendum are unchanged.

## D1. A bounded reverse search did not attain a forest

The main shared-port packet proves that Paley11 cannot be one legal path
complement away from any forest-framed tournament. A cheap next question is
whether exactly TWO complements suffice. We started from literal P11,
x→y iff y−x is a nonzero square modulo11, retaining every vertex name.

For a tournament T, exact subset DP computes its minimum feedback-arc count
f(T) and samples optimizing rank orders through DP ties. For each chosen
rank, a second exact subset/endpoint DP maximizes the number of backward
rank edges used by a legal Hamilton path. That overlap makes the inherited
transported exception formula k+(n−1)−2s favorable. This criterion chooses
paths; it does not certify globally smallest feedback count after complement.

The deterministic seed is7307311. The declared maximum domain was:

-64 first optimal-rank draws, each with8 maximum-overlap Hamilton-path draws;
- for each distinct first labelled target,8 optimal-rank draws, each with4
  maximum-overlap Hamilton-path draws;
- for each distinct second target with f≤10,32 rank draws (16 optimal, then
  permitted rank-cost slack at most2, always capped so cost≤10).

We deduplicated exact labelled tournament rows, not isomorphism classes.
Every generated Hamilton path was checked for legality before complement.
The search stops on the first explicit forest, or when that declared finite
budget is exhausted. No extra batch was run after the negative result.

### Exact observed outcome

| Quantity | Value |
|---|---:|
| First rank draws |64|
| First path draws |512|
| Distinct labelled first targets |337|
| Smallest first-target f |15|
| Second rank draws |2,696|
| Second path draws |10,784|
| Distinct labelled second targets |8,517|
| Smallest second-target f |13|
| Final promising forest-rank draws |0|

A forest on11 vertices has at most10 edges, so none of the sampled second
targets could have ANY forest frame. The absence of a forest frame for these
8,517 exact targets is certified by their exact f values, but the domain does
NOT exhaust legal two-complement words. In particular maximum-overlap path
selection and sampled optimal-rank ties are filters. Two-complement existence
remains OPEN. The finite negative result is not a new all-word obstruction.

The receipt retains an explicit best pair of paths P,Q, the intermediate and
final literal tournament rows, an optimal final rank with13 exceptions, its
actual cycle debt, and the two-complement clock. This is a reproducible near
miss, not a forest certificate. Starting at P11 the reverse sequence is
`tau_Q(tau_P(P11))`; reversing it would use Q then P. Each path is legal at
its own named stage; a whole-path remarking separates the two complements.
No elementary cut-and-insert relocation count is inferred from that remarking.

## D2. A source-preserving necessary complement clock

Use the proved main-packet tariff

    |f(tau_P T)−f(T)|≤c(n),
    c(1)=c(2)=0, c(3)=1, c(n)=n−3 for n≥4.

A forest-framed source at fixed order n has f≤n−1. Legal path remarking and
full converse preserve f. Hence any word with d actual whole path complements,
interspersed only with such operations, that reaches a target T must satisfy

    d≥ceil(max(0,f(T)−(n−1))/c(n))                 (D2)

when c(n)>0. This is a necessary lower bound conditional on a construction,
not an existence theorem or a replacement for the user's physical move clock.
The actual specified primitive remains one whole legal path complement.

For the inherited family P7[C3^[k]], n=7·3^k, the main packet proves
f=n(n−1)/6. Therefore:

| Target order | Exact f | Minimum possible complement count |
|---:|---:|---:|
|7|7|at least1|
|21|70|at least3|
|63|651|at least10|

The order-seven bound is attained by the inherited explicit forest→Paley7
witness. No attainment is claimed at21 or63. For Paley11 exact finite DP
finds f20, and(D2) gives at leastceil((20−10)/8)=2. D1 does not settle whether
that lower bound is attained.

The root synthesis proves the analogous resource bound for monotone disjoint
J/star construction trees from singletons: joins add f, stars add at most1,
and every complement is charged at its actual intermediate size. There are
at most n−1 star nodes, and all intermediate sizes are at most the final n.
Thus the same displayed lower count applies when these are the ONLY allowed
primitives (with relocation and converse). It does not apply to arbitrary
substitution supplied as a free constructor, copying whole subobjects at no
cost, or uncharged changes of the named source state.

All source masses remain attached to their original named records. None of
the randomized tie choices rerandomizes the mathematical source law; they
are finite search choices only. Integer complement clocks and the actual
shared port/path names are retained.

## Reproduction

    python3 search_depth.py > search.stdout.json
    cp search_result.json receipt.json
    python3 -O search_depth.py > search.optimized.stdout.json
    cp search_result.json receipt.optimized.json
    cmp receipt.json receipt.optimized.json
    python3 verify_depth_receipt.py
    sha256sum -c SHA256SUMS

`search_depth.py` is self-contained, uses exact integer subset DP, and has no
removable asserts. `verify_depth_receipt.py` independently checks the saved
best literal paths, tournament rows, final ranking, triangle/forest boundary,
and necessary clock arithmetic; it does not claim to certify an unperformed
search. `DEPENDENCIES.json` pins the frozen source theorem report and the
base repository commit used by this lane.
