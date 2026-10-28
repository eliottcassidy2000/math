# Audit and reproduction

**Same-agent mathematical audit and alternate exact computations.**
No independent-agent review, formal proof-assistant result, or global
conjecture closure is claimed. The source universe and failure boundaries
are stated in [THEORY.md](THEORY.md).

From the repository root, with Python3 and the standard library:

```powershell
python3 -B 05-knowledge/results/kernel125_closed_carriers_20261028/kernel125.py
python3 -B 05-knowledge/results/kernel125_closed_carriers_20261028/strong_ports.py
python3 -B 05-knowledge/results/kernel125_closed_carriers_20261028/root_minimum.py --max-order 9 --degree-first
python3 -B 05-knowledge/results/kernel125_closed_carriers_20261028/root_minimum.py --max-order 9 --degree-first --neighbor-degrees
python3 -B 05-knowledge/results/kernel125_closed_carriers_20261028/sequence_channels.py --large 200
python3 -B 05-knowledge/results/kernel125_closed_carriers_20261028/legacy_audit.py
python3 -B 05-knowledge/results/kernel125_closed_carriers_20261028/audit.py
```

The older minimum-root variant is reproduced by omitting --degree-first.
Its saved receipt predates the addition of degree-only instrumentation:
zero degree counters and the explicit false degree_first field were absent.
Its bank is byte-identical to the predecessor; the degree-first bank has
different root representatives, checked by an exact class bijection.
The later neighbor-degree counters are also absent from the two earlier
receipts; omitted counters mean zero in those earlier algorithms. The
auditor normalizes these instrumentation-only differences. Both improved
selectors are replayed through8 under optimized Python.
Timing fields are descriptive and are excluded from semantic comparisons.
The new code does not save a duplicate compressed bank; reproduction checks
against the existing bank only after constructing its answer.

| Receipt | Explicit universe | Positive / hostile or alternate control |
|---|---|---|
| kernel125.json | all125 phases, all32768 binary rows, all64 linear subspaces, all125 actual hub attachments,3125 five-cell phase tuples | Direct binary responses agree iff relative coordinates agree; native guard separates the diagonal fibre; actual converse includes polarity; native complement leaves T5. |
| strong_ports.json | all5040 permutations of P7,21 actual matching attachments,21^3 four-cell relative configurations | Literal tournament automorphisms vs semilinear action; nonnormal but separating promised orbit; whole-field collision; native complement exits family. |
| root_degree.json / root_neighbor.json | all parent-mask orbits through9,191536 final classes for each selector | Independent Burnside totals; every prior automorphism-order histogram; full child-group calculation through7; exact full-bank bijections;81/27 control. |
| root_minimum.json | same universe, original pointed-key order | Entire order-nine bank byte-identical to predecessor; same hostile. |
| sequence_channels.json | m=0..30,50,100,200; three scalar channels | Independent rational partition sums through30; independent ordinary future-response counts at every executed order; native self-converse census through8. |
| legacy_repairs.json | original and repaired producers; corrected functions0..30 | Cyclic pair-action average through30; literal ordered-pair tracing on all cycle types through8. |
| audit.json | saved outputs, source hashes, alternate short executions | Full receipt consistency, explicit witness checks and documented links. |

The previous F64/F256 field calculations are inherited and source-pinned,
not rerun as new tests here. The actual Paley test shows why a normality
hypothesis may be weakened on a promise set, and why extending the
operation language can invalidate that weakening.

Review gates covered: all iff directions; named versus unmarked outputs;
phase count versus unlabelled class count; chart group versus its order;
global converse versus native complement; polarity under actual converse;
linear invariant quotients versus valid observer quotients; nonabelian
composition order; exact weights versus real approximations; integer
divisibility; and finite verification versus all-order proofs.

The new class generator inherits the exact parent-isomorphism oracle.
Its colour refinement only selects buckets; a literal bijection preserving
every arc is required. The reference class bank is audit-only input.
The scalar engine inherits the existing contracted graph builder; sharing
that builder is not an independent proof of its older implementation.
Its direct rational partition comparison is the alternate path.

Source hashes are recorded in [source_pins.json](source_pins.json).
The legacy audit additionally records the full pre-repair Git commit.
The historical control is intentionally reproducible even though its
producer was false. No historical output file is overwritten.

The shared legacy file received the same parity repair concurrently in
Session87. Rebase retained its stronger front warning and our explicit
minimal witness and repair lineage. The merged source was re-audited.
