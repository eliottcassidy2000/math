# Aliquot inverse-fiber research packet

Read [REPORT.md](REPORT.md) for proved bounds, current literature status, exact computations, and explicit open boundaries.

**Concrete scoped result:** Chen–Zhao's established lower-density theorem, evaluated with certified integer arithmetic at an enlarged parameter, gives

`liminf U(X)/X > 0.060452936134095`.

This is a stronger instantiation than their displayed0.0602757 parameter, not a claim of literature-record novelty, natural-density existence, or a second odd untouchable.

**Reproduce:**

```sh
python3 scripts/verify.py
python3 -O scripts/verify.py
python3 scripts/density_lower_bound.py
python3 scripts/verify_density_monotonicity.py
```

The complete fiber and prefix checks use proved source bounds. All dependencies are Python standard library. Independent audit is separately supplied to the parent. [Sources](sources/PRIMARY_SOURCES.json); [fiber/census output](results/verification.json); [density enclosure](results/density_lower_bound.json).
