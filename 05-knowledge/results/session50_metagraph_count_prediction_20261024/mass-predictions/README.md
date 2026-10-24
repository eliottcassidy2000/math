# Session50: finite odd-mass moment prediction

Read `THEOREM.md` for proofs and scope. The main result converts a short vector of fixed-path fibre/collision moments into an exact tournament class count, using factorial divisibility and rational reciprocal bounds. It is not a closed n-only enumeration formula.

Two blind order-eight certificates (plus a short hand-checkable upper polynomial in `hand_certificate.json`):

1. Eight full moments, no support bound: rational interval plus residue modulo32768 isolates6880.
2. Twelve modular moments, plus N=2^21: residue modulo2^22 alone isolates6880.

Run with Python3.10+ and the standard library:

    python3 predictor.py --moments eight_moments.json --output fresh-blind.json
    python3 predictor.py --moments twelve_residues.json --modular-universe-size 2097152 --output fresh-modular.json
    python3 check.py --output fresh-normal.json
    python3 -O check.py --output fresh-optimized.json
    cmp fresh-blind.json blind_prediction.json
    cmp fresh-modular.json modular_prediction.json
    cmp fresh-normal.json normal.json
    cmp fresh-optimized.json optimized.json
    python3 verify_packet.py

`predictor.py` never reads an atlas. Its eight-moment mode receives only a JSON list of integers. Its modular mode receives only a list of residues and the explicit universe cardinality. All calculations use integers and exact rational arithmetic. The predictor does not certify that arbitrary user-supplied numbers are genuine moments; provenance is checked separately by `check.py`.

`check.py` verifies byte-pinned historical inputs, regenerates the frozen moments, runs the blind predictor, and only then validates against held-out atlas totals. It checks orders3..8, exact polynomial identities, modular divisibility, half-step Newton denominators, the n-only h=1 recurrence, and positive/hostile controls. It uses no bare assert. The upstream order-eight per-row H/Aut audit is inherited and credited, not silently rerun here. A separate independent Session50 audit recomputes those original rows.

`verify_packet.py` verifies all manifest bytes, copies only manifested files into a genuinely new temporary directory, and reruns normal, optimized, eight-moment, and modular predictions there. Historical prose/data snapshots have `.txt` extensions and are not imported as executable code. The source pins name exact paths and commitf26bc37f4e77f72a7ce4e4f9eabf2b271acc3390.
