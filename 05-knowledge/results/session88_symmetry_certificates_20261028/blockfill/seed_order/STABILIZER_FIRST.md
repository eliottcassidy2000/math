# Stabilizer-first seeds: a small saving and a concrete failure of monotonicity

Session88 bounded follow-up. **PROVED** ordering observation; **FINITE-EXACT**
complete-output comparisons through order eight. No runtime or novelty claim.
This experiment does not change the preceding published algorithms.

## Ordering observation

For each pointed state j in the complete bank, the full pointed stabilizer
H_j is already available from the canonical deletion parent's automorphisms
that stabilize its attachment mask. In a fixed underlying child U, let m_j be
the size of j's vertex orbit. Orbit-stabilizer gives

    |H_j| = |Aut(U)| / m_j.

Process all unfilled pointed states in decreasing |H_j|, breaking ties by the
existing named key. The first visited state in each ordinary block maximizes
|H_j| within that block, hence minimizes m_j. This requires no advance knowledge
of the child's full automorphism group or of its block.

This says **smallest orbit**, not **minimum adaptive classification cost**.
Point stabilizers at different vertex orbits need not be nested, and the
adaptive algorithm must discover each pointed type, not merely recover the
full group. Its next queried vertex and the labeled presentation also matter.

## Exact experiment

`check_seed_order.py` changes only the order in which first block rows are
constructed. It calls the unchanged adaptive row engine. Because the first
seed may change, it deliberately retains the original least named state as
each output representative, then orders those representatives exactly as before.
Consequently the comparison checks equality of all six original output fields:
complete Q, state bank, ordered raw child bank, all vertex-indexed transporter
sets, all full automorphism groups, and parent alignment data.

| Child order | Original adaptive queries | Stabilizer-first queries | Changed block seeds |
|--:|--:|--:|--:|
|1|0|0|0|
|2|1|1|0|
|3|3|3|0|
|4|9|8|2|
|5|38|37|3|
|6|244|242|7|
|7|2,593|2,586|25|
|8|47,391|47,377|87|

At eight, every changed seed has a strictly smaller root orbit. Across all
6,880 blocks, 15 save one query, **one costs one extra query**, and 6,864 are
unchanged: net saving 14 queries, about 0.030%. The total canonicalizer requests
including the unchanged raw-mask bank/preparation are 106,657. Stabilizer sizes
and ordering have additional costs; no matched timing benchmark was performed.

## Explicit hostile: a smaller root orbit costs more

One order-eight block has the named row

    ((131074,103),5), ((155648,2),3).

Its full automorphism group has order 15. The original seed `(131074,103)` has
orbit size five and stabilizer order three. It constructs the row with **one**
new root classification. Stabilizer-first chooses `(155648,2)`, orbit size three
and stabilizer order five, and needs **two** classifications.

The reason is directly visible in the fixed adaptive engine's first query:

| Seed | Literal child code | Type at first queried vertex 0 | Full group after first query | Known pointed types after first query |
|---|--:|---|--:|--:|
|(131074,103)|19922946|(155648,2)|15|2|
|(155648,2)|221478976|(155648,2)|15|1|

Both first queries recover the full group. From the first seed, the query also
discovers the other root type, so propagation completes the row. From the
smaller-orbit seed, the first query stays within the already-known type; the
five-vertex orbit's type still needs one classification. Thus **full-group
recovery is not yet complete named-root classification**. This supplies a
specific failure mechanism, not just an abstract warning about incomparable
stabilizers. It does not rule out a different adaptive query policy.

For a contrasting positive block, the seed changes from `(67858,56)` with
orbit seven and stabilizer three to `(85298,127)` with orbit one and stabilizer
21, reducing two queries to one.

## Reproduce

From this directory in the repository:

```sh
python3 -B check_seed_order.py
python3 -B -O check_seed_order.py
```

For staging, append `--baseline-root /path/to/repository`. The checker needs no
pre-existing JSON receipt, checks complete-output equality and every seed's
orbit-stabilizer relation, reconstructs the original per-block query totals,
and reproduces the hostile's first-query explanation. Normal and optimized
runs produced byte-identical JSON. The expected query table is checked only
after construction. The inherited baseline remains hash-verified.

This is the final bounded ordering test; no further heuristic expansion or
unmeasured acceleration is claimed.
