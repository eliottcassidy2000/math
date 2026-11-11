# Physical additive and quadratic controllers

- [PROOF.md](PROOF.md): one-way trades, all-length exact digit selector,
  greedy inverse, no-unit-deletion-bank theorem, conditional interval criterion.
- [EXPANSION.md](EXPANSION.md): mixed-bank parity resolution and causal
  forward/converse product gates, including a literal paid multiplication.
- `replay.py` / `RESULTS.json`: initial controls.
- `expand.py` / `EXPANSION_RESULTS.json`: expansion controls.
- `carrier.py`: relevant helper functions inherited from Session264 fusion code.
- `search*.py`, `binary_search.py`, `fiber.py`: bounded exploratory discovery
  scripts. They are not needed to reproduce any proof certificate, and their
  negative outputs are not global nonexistence claims.

From the repository root:

```sh
python3 -B 05-knowledge/results/session275_exact_selection_20261111/switches/replay.py
python3 -B 05-knowledge/results/session275_exact_selection_20261111/switches/expand.py
```

No finite coverage-range extension is claimed. The constructions select exact
integers in explicitly described infinite families and control their
interactions. Chan–Pak and eventual all-odd coverage remain OPEN.
