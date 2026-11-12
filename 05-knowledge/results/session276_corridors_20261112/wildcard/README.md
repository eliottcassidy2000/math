# Exact target languages and infinite-source readiness

Read [PROOF.md](PROOF.md) for target-capped automata and the complete all-order
three-reservoir fibre grammar. The [all-k hierarchy](../hierarchy/PROOF.md) extends the reservoir theorem
to every fixed labelled apex tournament and separates the target-controlled
transient degree from the apex-controlled polynomial growth degree.

Read [READINESS.md](READINESS.md) for the exact
binary-menu obstruction and a physical old-arc branch repair.

All counts concern named two-apex words, not tournament isomorphism classes.
Chan–Pak remains open. The all-height classifications depend on the inherited
positive transfer; the finite-core formula also uses the cited Busch floor.

Run from the repository root:

```
python3 -B 05-knowledge/results/session276_corridors_20261112/wildcard/capped.py 21 /tmp/COUNT_DFA21.json
python3 -B 05-knowledge/results/session276_corridors_20261112/wildcard/fibres.py 63 > /tmp/FIBRES63.json
python3 -B 05-knowledge/results/session276_corridors_20261112/wildcard/compare.py
python3 -B 05-knowledge/results/session276_corridors_20261112/wildcard/readiness.py
python3 -B 05-knowledge/results/session276_corridors_20261112/wildcard/selector.py
python3 -B 05-knowledge/results/session276_corridors_20261112/wildcard/insertion.py
```

Repeat with `python3 -O -B` for optimized-mode controls. The executable checks
raise exceptions rather than relying on Python assertions. `CLOSURE21.json`
contains a nondeterministic timing field; compare its mathematical fields,
not timing. Other JSON receipts are deterministic.
