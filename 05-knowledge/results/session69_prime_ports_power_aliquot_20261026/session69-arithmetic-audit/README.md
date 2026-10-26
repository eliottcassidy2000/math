# Session69 arithmetic audit

Read `REPORT.md` for the verdict, proof review, primary-source scope, and limits.

Self-contained independent computation:

```sh
python3 check_independent.py > independent.normal.json
python3 -O check_independent.py > independent.optimized.json
cmp independent.normal.json independent.optimized.json
```

`--negative-control` must exit nonzero in both modes with
`negative-control canonical first-six`. The supplied negative receipts record
exit code 1 in both modes; analogous producer aliquot receipts record its
target-six completeness failure.

The copied `power.*.json` and `aliquot.*.json` receipts are fresh executions
of the producer checkers, not merely copied original output. Validate the
producer manifests, all receipts, and all 999 complete inverse fibers using:

```sh
python3 verify_inputs.py --power PATH_TO_SESSION69_POWER \
  --aliquot PATH_TO_SESSION69_ALIQUOT > verification.json
```

This writes only this audit's `INPUTS.json`; producer files are read-only.
The independent checker imports no producer code. The final cross-comparison
deliberately calls the inherited inverse algorithm after the independent
complete-fiber digest has already been computed.

`SHA256SUMS` covers this packet, excluding itself. No push was attempted by
this reviewer; the parent owns integration.
