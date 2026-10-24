# Independent regular-eleven score-mass audit

PASS: a separately written fixed-labelled score recurrence, exhaustive tiny
score-vector controls, and full ordinary-group recomputation certify the
atlas-free closure: 48,251,508,480 labelled regular-eleven tournaments;
1,205 rigid and 1,223 total ordinary regular classes; P10=13,333;
E_regular11=133,210. This is a score-count/group-mass sidecar, not bare-metagraph
sufficiency and not an unrestricted order-eleven count.

See AUDIT.md for the proof and precise inherited-completeness boundary.

Run:

    python3 independent.py
    python3 -O independent.py
    python3 hostiles.py
    python3 -O hostiles.py
    python3 crosscheck.py
    python3 verify_audit.py
    python3 replay.py --output physical-copy-replay.json

Only files listed by MANIFEST.json are needed. The relocation receipt is an
external witness and intentionally not included in its own digest graph.
