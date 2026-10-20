# Blue maxima continuation

This add-only packet proves two complete new finite maxima and a separate all-order lower bound:

- `M_6=1185`, with 12 maximum-blue classes at order 13.
- `M_7=8295`, with 24 maximum-blue classes at order 15.
- `liminf M_k/k! >= sqrt(e)`, with explicit one-sided error bounds.

Start with the plaintext `REPRODUCE.txt` runbook. Read `BLUE_MAXIMA_K6_K7.md` and `SQRT_E_LOWER_BOUND.md` for the mathematical statements and proofs. The general exact maximum and a matching asymptotic upper bound remain open. Records beyond `k=7` are lower bounds only.

## Reproduce the exact finite results

Use an isolated copy of this directory. Python 3 and a C++17 compiler supporting unsigned 128-bit integers (GCC or Clang) suffice. No third-party Python packages are needed. Keep validation enabled. The producer verifiers also work under `python3 -O`; the separately authored `independent_audit.py` must run without `-O`.

```sh
g++ -O3 -std=c++17 enumerate_reduced.cpp -o /tmp/blue_reduced
g++ -O3 -std=c++17 replay_subcubes.cpp -o /tmp/blue_subcubes
g++ -O3 -std=c++17 enumerate_walsh.cpp -o /tmp/blue_walsh
g++ -O3 -std=c++17 replay_zeta.cpp -o /tmp/blue_zeta
/tmp/blue_reduced 6 > exact_k6.json
/tmp/blue_subcubes 6 A_reps_k6.txt reduced_k6.bin > replay_k6.json
/tmp/blue_walsh 6 reduced_k6.bin > walsh_k6.json
/tmp/blue_walsh 7 > exact_k7.json
mkfifo blue_replay.fifo
/tmp/blue_zeta 7 A_reps_k7.txt blue_replay.fifo > replay_k7.json &
replay_pid=$!
/tmp/blue_walsh 7 - blue_replay.fifo > replay_walsh_k7.json
wait "$replay_pid"
rm blue_replay.fifo
python3 -O verify_maxima.py
python3 -O verify_local_design.py
python3 independent_audit.py
```

The two `A_reps` files are the complete representative lists from the retained exact JSON records; the Python verifier independently checks their underlying orbit coverage via the JSON. The array producer and consumer communicate native signed 32-bit integers and must run on the same machine. The retained portable certificates are JSON; the generated `reduced_k6.bin` is not needed in the published packet.

The order 15 comparison streams 3.8 GB through the FIFO but does not store that amount on disk. Peak working memory is under 1 GB. The producer and independent replay took roughly 12 and 19 seconds respectively in the original run; timings depend on hardware. The complete Python orbit/witness verifier takes longer than the C++ transforms but remains bounded by 456*7! permutations plus the small witness enumerations.

The programs stop with a nonzero exit code on a mismatch. `verify_maxima.py` additionally checks full labelled means and second moments against an independent exact rational implementation, histogram identities, independent tournament-orbit coverage, direct signed-word witnesses, and matching-completion orbits. A `--skip-orbit-coverage` option exists for a quick nonexhaustive rerun; it should not be used for a full audit.

## Optional nonexhaustive search reproduction

```sh
g++ -O3 -std=c++17 search_local.cpp -o /tmp/blue_local
g++ -O3 -std=c++17 search_extensions.cpp -o /tmp/blue_extend
/tmp/blue_local 8 100000 > local_k8_deep.json
/tmp/blue_local 9 10000 > local_k9.json
/tmp/blue_local 10 20000 > local_k10.json
/tmp/blue_extend 6 314 3618 > extension_k7.json
/tmp/blue_extend 7 268 39375 > extension_k8.json
```

These deterministic searches use seed 20261020. Their forward/backward subset-DP gradient is checked against every one-bit finite difference at initial test states. `verify_maxima.py` independently recounts the reported witnesses in Python. These are not exhaustive global searches and must never be described as determining `M_8`, `M_9`, or `M_10`.

## Validation records

- `exact_k6.json`, `exact_k7.json`: complete reduced-census summaries and full weighted labelled histograms.
- `replay_k6.json`: complete independent subcube/DP agreement at `k=6`.
- `walsh_k6.json`: complete Walsh/DP agreement at `k=6`.
- `replay_k7.json`: all 956,301,312 ordinary-zeta/Walsh values agree.
- `replay_walsh_k7.json`: the corresponding second complete Walsh run.
- `verification.json`: independent orbit, moment, word, and completion checks.
- `witness_n13.json`, `witness_n15.json`: explicit maximum-blue tournaments and automorphisms.
- `local_design_verification.json`: exact probabilities on all 384 local signed permutations, an explicit affine 16-point K4 design, rational lower bounds, and quantitative termwise checks.

Independent review is APPROVED in `INDEPENDENT_REVIEW.md`. The separately authored `independent_audit.py` and its results, plus the reviewer's fresh `independent_replay_k6.json` and `independent_replay_k7.json`, are included. No remote writes are performed by the producer or verifier.
