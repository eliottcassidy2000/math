# Session105 — exact integer projective A000568 recursion

Read [INTEGER_PROJECTIVE.md](INTEGER_PROJECTIVE.md). Entry point:

```python
from integer_response import count_packed
value, statistics = count_packed(100)
```

This combines Session98's fixed-degree quotient with Session87's integer
arithmetic and packed degree-bucket implementation. It is a single-target
counter, not a tournament/metagraph generator. `a_0=1`; all odd permutation
cycle types and their centralizer weights are included.

* `integer_response.py`: transparent tuple reference, cached transfer variant,
  3-tail ablation, packed full-factorial and packed odd-local realizations.
* `baseline98.py`, `baseline87.py`: byte-identical inherited source snapshots.
* `verify.py`: deterministic independent checks, active under Python `-O`.
* `benchmark.py`: all methods, actual stronger prefix baseline, repeated timings.

The projective integer realization improves the Session98 Fraction engine.
It is not uniformly faster than Session87's all-coefficients packed engine;
retain that implementation when measured speed or a whole prefix is wanted.
