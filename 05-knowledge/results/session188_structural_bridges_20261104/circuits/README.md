# Session 188 circuit branch: exact constant-region surplus

**PROVED local incidence identity; FINITE-EXACT regression checks.**

[PROOF.md](PROOF.md) derives

    |constant gates| + |their boundary gates| = 2r + f + E,

where r is the induced constant-region source count, f its one-parent node
count, and E its distinct-child fanout excess. Equality is exactly a disjoint
full binary merge forest ending in distinct unary boundaries. The result
extends Session 118's constant-antichain lemma without assuming independent
subcircuits or counting a shared gate twice.

The deletion consequence uses a free final literal, as the pinned elimination
candidate explicitly does. Strict B2 has a sharp possible one-gate terminal
NOT correction; the full bound applies when the restricted output depends
on at least two ordinary primary inputs. No asymptotic coefficient improvement
or external candidate verification is claimed.

## Reproduce

    python3 verify.py > receipt.json
    python3 -O verify.py > receipt_optimized.json
    cmp receipt.json receipt_optimized.json

Both runs passed, with identical output. The complete topology universes
contain 64,938 closed constant regions in live distinct-pin DAGs through six
gates and 29,202 in live repeated-pin DAGs through five. The Boolean census
contains every one- or two-gate presentation on three inputs with all 16
truth tables: 6,240 live presentations, each tested under all 18 coordinate
restrictions fixing at least one and leaving at least one input. It checks
11,712 nonempty constant-region cases with nonconstant output, including
5,856 cases needing the strict terminal correction. These are presentation
counts, not counts of distinct Boolean functions.

[source_pins.json](source_pins.json) records the exact recovered sources.
[AUDIT.md](AUDIT.md) lists the mathematical boundary checks.

The remaining frontier is unchanged: one must force additional f or E (or
some different surplus) frequently in minimum circuits with a small affine
constant-flat dimension. The exact accounting says where a surplus goes;
it does not prove that such a surplus exists.

## Sequential continuation

[SEQUENTIAL_BUDGET.md](SEQUENTIAL_BUDGET.md) proves that constant-region
fanout surplus cannot be reused: its cumulative value is bounded by the
initial gate-only fanout potential. A sharp k-stage alias ladder shows why
unary-boundary fanout is different: the same one-unit branch point can move
backward through k boundary deletions before being consumed once.

    python3 verify_sequential.py > sequential_receipt.json
    python3 -O verify_sequential.py > sequential_receipt_optimized.json
    cmp sequential_receipt.json sequential_receipt_optimized.json

Both runs pass. The checker exhausts all literal reroutes through five-gate
live distinct-pin topologies, verifies ladders through k=32, and compares
literal Boolean functions exhaustively at k=1,2,3.
