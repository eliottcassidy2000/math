# Exact A000568 queries: integrated recurrence and verified cache

**Delivered:** exact Python and C++/GMP single-coefficient engines, a fully
checked table through n=250, and reproducible same-machine benchmarks.
**Not delivered:** a new polynomial-time algorithm, an instantaneous arbitrary-n
formula, or a solution of Chan–Pak.

The GF/projective response quotient is inherited from Session98, including its
maximality theorem. This implementation combines it with degree-local
odd-factorial scaling, provably exact negative shifts, half-defect packing,
a late gauge switch and unit-tail contraction. The Python version additionally
contracts {1,3} by a precomputed Laurent polynomial with one weight multiplication
per state. See `PROOF.md` for the exact invariants and complexity bound.

## Practical results

Interleaved medians of three runs, same machine, in seconds:

| Target | new C++ projective | new C++ raw | new Python projective+terminal | new Python raw single | inherited Python prefix | historical C/GMP partition |
|---:|---:|---:|---:|---:|---:|---:|
|80|0.00714|0.00754|0.01159|0.01153|0.01977|0.02093|
|120|0.03064|0.03356|0.07992|0.08698|0.10895|0.76498|
|160|0.21109|0.20287|0.45966|0.49385|0.59156|not run|
|200|1.09560|1.10662|2.12423|2.22354|2.50235|not run|

At 250, a single additional run measured C++ projective 8.713s, C++ raw 9.045s,
Python projective+terminal 12.584s, Python raw-single 12.944s, and inherited
Python prefix 13.955s. These are single-run observations, not medians.

These are practical query comparisons, not all identical outputs: the inherited
Python baseline returns the entire prefix while the other engines return one
term. The C++ raw/projective comparison has the same requested output and
language. It shows that projective state compression is not a robust wall-time
win by itself (it loses slightly at 160). The roughly 25x speedup over historical
compiled partition enumeration at 120 combines the inherited response DP and
implementation, not a newly discovered mathematical formula. The best recovered
Python baseline is much stronger than historical partition code.

See `receipt.json` for all samples, source hashes, compiler commands embodied
in the runner, transition counts, and the additional single-run target250
comparison. The shared host is not an isolated benchmark machine. C++ timings
in the `seconds` arrays include launching the process and capturing its full
integer output; the nested C++ telemetry reports kernel time separately.
Python timings exclude decimal conversion. Comparisons therefore do not hide
C++ startup or output time.

The isolated Python ablation (`ablation.json`) measured, at 200, raw response
2.320s, projective direct 2.290s, and projective plus {1,3} terminal 2.232s
(medians of three). The terminal alone saved about 2.5% in this small sample;
this is not a large or asymptotic acceleration.

## Near-instant bounded queries, honestly accounted for

`prefix.json` stores **all 251 exact values a(0),...,a(250)** as decimal strings.
Every coefficient was compared between the independent degree-local
odd-factorial engine and Session239's tuple/factorial response engine. Building
the table took 17.38s in the first engine plus 45.87s for independent verification.
The file is 730,950 bytes, including indentation, and is checksum-checked on load.

Loading, checking the hash and parsing all values took about 0.017–0.023s.
A fresh Python CLI process querying a(250), including interpreter startup,
file verification/parsing and capturing the full decimal output, took roughly
0.05–0.08s. Filesystem caches were **not** flushed. Warm list indexing is tens
of nanoseconds in the saved loop measurement, but this excludes formatting
or transmitting a large integer; do not call that end-to-end query latency.
The table costs O(sum_(n<=N) bitlength(a(n)))=O(N^3) bits asymptotically; a warm
lookup returns an existing integer reference in O(1), while emitting it costs
at least its output length. This is preprocessing, not a new enumeration
complexity result. Querying beyond 250 falls back to exact fresh recurrence and
can be slow.

## Reproduce

From the repository root:

```sh
D=05-knowledge/results/session270_chanpak_enumeration_20261111/enumeration
python3 "$D/cache.py" 100                 # verified table query
python3 "$D/cache.py" 250                 # largest saved exact value
python3 "$D/cache.py" 120 --fresh         # Python recurrence, no table
python3 "$D/cache.py" --time-cache        # fresh-process query timing
python3 "$D/cache.py" --build 250         # rebuild and verify every coefficient
python3 "$D/benchmark.py"                 # independent controls and comparisons
python3 "$D/ablation.py"                  # isolate the Python terminal contraction
g++ -O3 -std=c++17 "$D/counter.cpp" -lgmpxx -lgmp -o /tmp/a000568-s270
/tmp/a000568-s270 200                     # projective integer recurrence
/tmp/a000568-s270 200 -1                  # same-language raw-response ablation
```

The benchmark retrieves unchanged historical baselines with `git show` at
pinned commit `f2e03e29358d7a6e566f27f75794ec601d0449f5`; it needs no network.
Compiled binaries are deliberately not committed. Python requires only its
standard library; C++ requires a compiler and GMP development libraries.
Direct independent odd-partition controls cover every n<=30, every tested
late switch and both compiled modes. Larger results agree between all engines
run at those targets and with the independently checked whole-prefix cache.

## Final source provenance

`receipt.json` preserves the actual timed Python source hash `75eeba6c...`.
After those timings, one explanatory comment was corrected from
`b_j=j*q_1-q_j` to `b_j=(j*q_1-q_j)/2`; executable statements did not change.
The final file hash is `63f8cab0...`, also used by the isolated ablation.
`verify_final.py` reverses exactly that comment, reconstructs the full timed
source hash, compares parsed Python ASTs, and reruns final Python/C++ controls.
`FINAL_CHECKS.json` records the outcome and hashes **every final Python script**
and the final C++ source. `SHA256SUMS` additionally covers all delivered files.
Historical baseline hashes remain pinned to the declared repository commit,
not silently relabelled as final implementation hashes. Redundant text-output
copies were removed; JSON receipts retain their results.

Run `python3 "$D/verify_final.py"` to repeat the final packaging controls.
