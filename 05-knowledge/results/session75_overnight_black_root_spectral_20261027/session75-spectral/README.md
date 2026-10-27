# Session75 spectral exception budgets

Read [RESULT.md](RESULT.md) first. Self-contained all-shape proof, sharp-leading
Paley arc debt, separate rank/block and native-clock budgets, and exact
join/C3 edit distances **3 at Paley7 / 11 at Paley11**. Independent parent
audit is separate; no ROOT positivity or native reachability closure is claimed.

* `check.py`: standard-library exact verifier, independent Bellman partition
  implementation, field27 and tower spectra, unbalanced trees, native XOR,
  clone loss and unequal pair-weight controls.
* `certificate.json`: complete subset lower bounds and attaining labelled
  parses; `generate_certificate.py` independently regenerates it.
* `receipt.normal.json`, `receipt.optimized.json`: identical successful replays.
* `test_corruption.py`: four rejection controls, also identical under `-O`.
* `SOURCE_PINS.json`, `sources/`: exact inherited definitions and mechanisms.
* `SHA256SUMS`: file integrity manifest (not a mathematical audit).

From this directory:

```sh
python3 generate_certificate.py > /tmp/session75-spectral-regenerated.json
cmp certificate.json /tmp/session75-spectral-regenerated.json
python3 check.py > /tmp/session75-spectral-normal.json
python3 -O check.py > /tmp/session75-spectral-optimized.json
cmp /tmp/session75-spectral-normal.json receipt.normal.json
cmp /tmp/session75-spectral-optimized.json receipt.optimized.json
python3 test_corruption.py
python3 -O test_corruption.py
sha256sum -c SHA256SUMS
```
