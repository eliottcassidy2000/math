# Independent positive-controller replay

These files are unchanged copies of the Session137 main auditor's separate
implementation and review. They import no producer module.

```
python3 -B audit_width_universality.py > check.json
cmp check.json width_universality.normal.json
python3 -B -O audit_width_universality.py > check.optimized.json
cmp check.optimized.json width_universality.optimized.json
```

Both relocated replays were executed and matched the original receipts.
`audit_width_states.py` is imported here only for its independent profile
engine. Its separate Phase1 command-line entry point expects an additional
Phase1 witness snapshot; that entry point is not part of this reproduction.
The universality checker and its imports are fully self-contained as shipped.
