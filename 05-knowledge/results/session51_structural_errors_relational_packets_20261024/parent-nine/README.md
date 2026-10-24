# Canonical parent growth and order-ten prediction

Starting from the authenticated order-eight action-equipped atlas, this packet
constructs191536 ordinary order-nine parents without nine-vertex canonicalization,
then evaluates97213472 pointed order-ten extensions by color refinement alone.
It obtains U10=97336068 and predicts a10=(U10-5508)/10=9733056.
The5508 correction is independently proved in the pinned structural packet.

- `build_nine.cpp`: invariant singleton selection, residual parent-only marked
  keys, full automorphism recovery by parent transporters.
- `trace_ten.cpp`: refinement-only target trace; FULL_REFINEMENT disables its
  safe singleton early exit for a complete independent-path replay.
- `parents-eight.txt`: pinned complete input with actual actions.
- `parents-nine.txt`: generated complete parent atlas with actual actions.
- `THEOREM.md`: selection, action completeness, target trace, guards and limits.
- `reproduce.py`: compiles and reruns both builder modes and both trace modes.
- `verify_packet.py`: checks hashes, output equality, arithmetic and boundaries.

All guards use exceptions and survive optimized builds. The production sources
contain neither the order-nine output count nor the order-ten target count.
Only eight-vertex parents are canonicalized; no target-order-ten isomorphism,
deletion or automorphism computation occurs.

Input provenance is explicit: the original certified JSON atlas and its
completeness/action certificate are archived in `sources/`. The verifier checks
the JSON SHA against that certificate and compares every text-input row and
action against the JSON, not merely its class count or labelled mass.

Reproduction (a separate empty directory is recommended):

    PYTHONDONTWRITEBYTECODE=1 python3 -B reproduce.py --output-directory /tmp/order-ten-replay
    PYTHONDONTWRITEBYTECODE=1 python3 -B verify_packet.py
    PYTHONDONTWRITEBYTECODE=1 python3 -B -O verify_packet.py

The full two-mode replay takes several minutes. Compiled binaries and Python
caches are not deliverable payloads and are excluded from the manifest.
