# Shared-C3 observers and sandwich resources

The all-size statement and inheritance boundaries are in [THEOREM.md](THEOREM.md). This packet extends a known square/diagonal centre law by an actual two-cell global-observer construction, capacity-sensitive off-diagonal restriction, full kernel transports, two different work counters, and explicit common-target maps.

Main distinction: for each fixed pair of distinct native prefix keys, the two local slots retain their own normalized resources through all three C3 phases. One global observer realizes both local partners. The local occurrence work is6; global observer work is3. Diagonal pairs still admit a global observer but double-use a local resource. Every selected global presentation reaches one common strong self-converse target class, so the growing resource count is not a distinct-neighbor count.

The square/blue profile and nonstrong-source looplessness are explicitly credited to THM-2454 and THM-4729. The values1,3,9,1729 and the failed power-of-three extrapolation are also inherited and independently verified. No numerical harmonic is transported between different tournament orders.

## Reproduction

Run from this directory:

    python3 verify_manifest.py
    python3 check.py > /tmp/sandwich-normal.json
    python3 -O check.py > /tmp/sandwich-optimized.json
    cmp /tmp/sandwich-normal.json normal.json
    cmp /tmp/sandwich-optimized.json optimized.json
    cmp /tmp/sandwich-normal.json /tmp/sandwich-optimized.json

Only the Python standard library is used. All proof inputs are vendored text archives with original paths and immutable commit provenance in `SOURCE_PINS.json`. No runtime imports from a repository, network access or Git are required. Checks use explicit exceptions, not removable assertions. `MANIFEST.json` records every durable file except itself; Python runtime `__pycache__` files are excluded. The verifier rejects missing, altered and unlisted durable files.

## Exact checked universes

- Every labelled R on1..4 owners,75 cases. Each check directly enumerates all sandwich paths and automorphisms as an independent product-law control.
- Complete R paths, full adjacency-pruned automorphism enumeration and every normalized R word for C3, cyclic C5, Paley7 and Paley11. H is also recomputed by subset DP.
- Every pair of R path-orbit representatives, and every central phase, for the C3/C5/Paley7 sandwich profiles. Paley11's large sandwich profile is derived from the proved pair bijection, not falsely reported as a full pair census.
- Native r=1..4, global orders19,35,51,67: all binary key pairs and three phases, actual matrix restrictions, both local resources, global resource quotient, both costs, literal source/partner phase transports and common-target maps.
- A labelled order10 local-capacity hostile for arbitrary R, with a genuine pair of distinct R path orbits that become reflection-paired after appending C3.
- Empty-wing hostile showing why the sandwich closing-cycle proof requires nonempty R.

The normal and optimized receipts are byte-identical. `NORMAL_O_EQUIVALENCE.json` and `PORTABILITY_REPLAY.json` record execution and relocated manifest-only reproduction.
