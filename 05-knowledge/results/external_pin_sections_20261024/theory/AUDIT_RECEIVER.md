# Independent audit of conditional receiver migration

Status: **VERIFIED proof and implementation review.** Reviewed
`RECEIVER_BANKS.md` and `check.py`. Normal and optimized frozen replays both
pass **1,040,572** checks across **42,554** complete ordered endpoint routes.
The generation count is one less because replay also compares the receipt.
No producer or receipt changes were required.

## Quantifiers and full records

The all-size statement is conditional on a complete receiver section at
EVERY permitted actual frame, grade in the declared band, and requested
charge. This is stronger than marginal charge coverage and is explicitly
retained. The helper graph contains full compatible records and ignores
only the total charge; its edges already preserve the helper population,
internal predicates and frozen boundary. Its connectedness is an input,
not an inference from projected charge supports.

The deficiency lemma counts distinct charge decompositions with one chosen
full witness per supported charge. Each decomposition determines exactly
one chosen pair, so subtracting the bound on incompatible witnessed pairs
is sound. The strict inequality and both equality hostiles are correct.
Boundary-dependent intervals require suitable authenticated splits; no
universal supply or native realization is inferred. The single-slot
difference-set criterion likewise requires the two actual options to glue:
their labels alone cannot establish that condition.

## Termination and exact restoration

For an endpoint with receiver grade below l, enough selected outside slots
exist because k is at least l. For a grade above u, enough outside vacancies
exist because k is at most u+e. The normalized grade
r0=max(l,k-e) lies in the section band. An endpoint below l necessarily has
r0=l, so the procedure never needs an unsupported upward normalization.
Above r0, an outside vacancy exists at every decrement. These facts justify
both entry branches and the terminating normalization bound u-l+1.

After normalization the outside sets have the same size and therefore join
by Johnson exchanges. Conditional sections pay every arrival at its actual
frame. The given helper path then connects the required helper endpoints.
The receiver-only reset restores the literal target-normalized record before
reversing its saved normalization. This reset is necessary even when both
skeleton distances vanish; the explicit final +1 in the event bound is valid.

When E is empty, every feasible grade is already in the band, so no outside
entry, decrement, or maximum over an empty outside cost set is invoked.
A singleton band permits entry directly to that grade and needs no
consecutive-grade changes. Empty receivers are also consistent when the
section hypothesis is satisfiable, for example for a trivial charge group.

## Costs, complements and native corollary

The sharp grade-only bound min(r+s,2h-r-s) is the maximum possible symmetric
difference. Summing the largest clamped number of nonnegative receiver
switching costs gives a valid weighted upper bound. The lower-entry and
upper-entry displayed formulas are safe, possibly weaker bounds; each uses
the appropriate grade branch. Adjacent changes, outside exchanges, helper
edges and receiver-only refresh retain their actual occurrence costs.
The theorem explicitly charges attached controls through the full frame-edge
cost rather than treating a shared owner as free.

Core complementation is a coordinate change with NEGATED core labels,
unchanged helper charge and a fixed all-core offset. It is never inserted as
an unguarded event. The finite-native three-cheap-erasure corollary separately
uses authenticated repair labels and private full records. Its interior
20/(49+c) estimate and outer direct-event20/85 estimate are compatible:
outside the complete band, the smaller free-core population has at most
seven selected or unselected records, bounding its changed core support by
14. The outer argument asserts connectivity only between feasible endpoints,
not charge coverage. No all-order native supply follows from this corollary.

## Independent boundary controls

In addition to both frozen replays, imported the unmodified checker and
exhausted every ordered feasible endpoint pair in three additional models:

| Model | Complete endpoint routes |
|---|---:|
| C2, receiver labels(0,1,0,1), outside(0,1), singleton band[2,2] |3,834|
| Trivial group, empty receiver, two outside slots, band[0,0] |54|
| Trivial group, two receiver slots, empty outside, band[1,1] |36|

These controls pass **82,801** checks over **3,924** routes. They are frozen
in `audit.py`, which also independently enumerates literal binary states,
their populations and charges, then checks every compiled arrival, event
type and occurrence cost: **94,565 independent checks / 10,702 events**.
Normal and optimized executions agree. The singleton model exercises both
low and high entry branches. An additional literal same-frame transition
between the two different grade-one receiver records in the empty-outside
model requires exactly one receiver-only refresh, confirming the endpoint-
restoration sidecar independently of the main models.

Frozen replay commands, from the repository root:

    python -B 05-knowledge/results/external_pin_sections_20261024/theory/check.py
    python -O -B 05-knowledge/results/external_pin_sections_20261024/theory/check.py
    python -B 05-knowledge/results/external_pin_sections_20261024/theory/audit.py
    python -O -B 05-knowledge/results/external_pin_sections_20261024/theory/audit.py
