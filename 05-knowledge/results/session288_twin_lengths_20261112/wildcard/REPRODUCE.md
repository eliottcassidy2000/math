# Reproduction

No third-party packages for the proof certificate and literal witnesses:

```sh
python3 05-knowledge/results/session288_twin_lengths_20261112/wildcard/verify.py
python3 -O 05-knowledge/results/session288_twin_lengths_20261112/wildcard/verify.py
cmp 05-knowledge/results/session288_twin_lengths_20261112/wildcard/REPLAY.json 05-knowledge/results/session288_twin_lengths_20261112/wildcard/REPLAY.optimized.json
```

This computes the complete finite two-port skeleton, compares literal
Hamilton/path-cover counts for sixteen chain blowups, and independently
counts the15- and17-vertex arbitrary-module witnesses. It checks the ray
criterion for48 bounded sample directions; the all-direction statement is
proved symbolically, not inferred from those checks.

Optional discovery replay requires `numpy` and `numba` (tested numpy2.5.4,
numba0.68.1), and the inherited census blob pinned in SEARCH_PROVENANCE.json:

```sh
python3 05-knowledge/results/session288_twin_lengths_20261112/wildcard/screen.py
python3 05-knowledge/results/session288_twin_lengths_20261112/wildcard/modules.py
```

The loader uses the repository copy when present and otherwise reads the
pinned Git blob without modifying the repository. It verifies the source
hash. `screen.json` is a discovery census, not a needed assumption of the
new all-order gate or universal module identity. `module_pairs.json` gives
one representative for each detected ordered profile pair, not every
isomorphism pair. Both scans inherit completeness of the earlier ordinary
class inventory; no new full census claim is made.
