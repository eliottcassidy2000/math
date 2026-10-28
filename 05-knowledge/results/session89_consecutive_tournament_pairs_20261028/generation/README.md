# Consecutive tournament-pair generation

See [ALIGNED_GENERATION.md](ALIGNED_GENERATION.md) for the proof, precise output
contract, nested/independent/overlapping-stream distinctions, measurements and
scope boundaries.

- `check_pairs.py`: three complete-bank algorithms; literal independent audit
  through core order 4; matched full-output enumeration through core order 5.
- `check_lift.py`: independent full self-converse lift-bank audit through order 7.
- `verify.py`: pinned replay of deterministic results, hostile controls and lift.
- Running the benchmark produces `results_n5.json`; raw measurement receipts are not required inputs.
- The checkers generate their own local JSON outputs; none are required committed inputs.

From this folder inside the repository:

```sh
python3 -B verify.py
python3 -B -O verify.py
```

From an external staging folder, add `--baseline-root /path/to/source-repository`.
Only the unchanged byte-pinned Session87 and Session88 files are imported. No
baseline source files or saved result receipts are modified or needed as input.
