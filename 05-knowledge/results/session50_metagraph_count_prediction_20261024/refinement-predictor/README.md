# Blind refinement predictor at order eight

**FINITE-EXACT:** A000568(8) is obtained as `(55068 - 28) / 8 = 6880`
from complete ordinary parent representatives through order seven, their actual
permutation actions, attachment-subset orbits, and directed color refinement.
The prediction stage performs no child canonicalization, no child automorphism
search, and no deletion/re-rooting isomorphism tests.

- [THEOREM.md](THEOREM.md) states the exact input contract and count mechanism.
- [VALIDATION.md](VALIDATION.md) records the verified boundary and controls.
- `predictor.py` is the only blind prediction entrypoint.
- `freeze_inputs.py` and `validate_inputs.py` are explicitly separate
  preparation/validation stages. They can use lower-order isomorphism tests.
- `parents-through-seven.json` is the pinned complete input atlas.
- `receipt.json` records the result and operation counts, not an input.
- `source-pins.json` and `archive/` retain literal historical source bytes;
  historical Markdown is stored with `.md.txt` suffixes.

Reproduce the blind stage:

```sh
python3 predictor.py --receipt fresh.json
python3 -O predictor.py --receipt fresh-optimized.json
```

Reproduce the entire checked package from an actual isolated manifest copy:

```sh
python3 replay.py --manifest-only
```

The arbitrary-order statement remains a conditional correction identity, not
an all-order scalar recurrence. This entrypoint rejects target order nine.
The frozen holonomy theorem has explicit new strong order-nine errors, so the
order-eight correction must not be extrapolated. No timing or speedup claim
is made.
