# Session132: prime labels versus square-sum graph order

**PROVED:** exact degree/edge/arrival-root formulas; connectedness threshold
14; 2-connectedness threshold 31; girth thresholds 15 and 30.
**FINITE-EXACT:** all edge lists and invariants of G_1..G_97; Hamiltonian path
and cycle counts through order 32; explicit cycles at every prime 37..97.

The supplied lists have 25 entries but only 24 distinct primes: 71 occurs
twice and 67 is missing. No silent correction is made.

Two interpretations are included, rather than guessed:

- `prime_order_first25.md`: p_i labels G_{p_i};
- `rank_order_first25.md`: p_i labels G_i, for i=1,...,25.

The user's A exactly equals the disconnected prime-order graphs, explained
by G_n becoming connected at 14. Other supplied groups do not classify the
verified structural regimes. In the rank interpretation the nontrivial path
orders 15,16,17,23,25 carry prime labels 47,53,59,83,97.

See `PROOFS.md` for proofs, boundary witnesses, the map/loss ledger, and
limits. `graphs.json` contains exact edge lists; `witnesses.json` contains
checkable paths, cycles and the spanning open-ear certificate for G_31.
Run `render_atlas.py` to generate SVG contact sheets of the actual 25 graphs
under both interpretations. Rendered SVG/PNG files are intentionally excluded
from this text-only Git release; the frozen PNG previews are delivered separately.

## Portable reproduction

Only Python's standard library is required for generation and verification.
From this directory (or use absolute script paths from any directory):

```sh
python3 generate.py
python3 render_atlas.py
PYTHONDONTWRITEBYTECODE=1 python3 check.py --output /tmp/s132-prime-normal.json
PYTHONDONTWRITEBYTECODE=1 python3 -O check.py --output /tmp/s132-prime-optimized.json
cmp /tmp/s132-prime-normal.json /tmp/s132-prime-optimized.json
python3 verify_packet.py
```

The independent verifier regenerates the SVG contact sheets in a temporary
directory, so it runs directly from this text-only release. It does not import
the graph producer: it constructs edges
by testing every vertex pair and counts Hamiltonian paths/cycles with an
unpruned memoized recurrence. The producer uses square-root reflection
matchings and degree/connectivity-pruned search. Counts are for undirected
paths modulo reversal and cycles modulo rotation/reversal. Explicit gates
remain active under `python -O`; `--test-failing-gate` is a negative control.

The source-pinned Session131 odd-root channel theorem is inherited and
credited. No new all-order Hamiltonian existence theorem is claimed, and
no publication was performed by this lane.
