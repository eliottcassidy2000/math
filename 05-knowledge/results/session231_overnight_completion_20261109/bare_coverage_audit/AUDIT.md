# Independent audit: fixed-core bare singleton contexts

Verdict: PASS for the exact construction-specific identities and complete target-bounded census. This is not an unrestricted Chan–Pak spectrum certificate, nor an asymptotic density theorem.

## Independent computations

`audit.py` imports no producer code. It builds each prefix state by a singleton recurrence, composes the fixed core once, then explores a binary suffix tree using only three quantities (A,T,E). Appending c maps these to (2A+45c,T+c(A+E),45c). Empty suffixes are tested separately. Every suffix ending in zero is tested at its last edge; a branch is pruned only when its current A exceeds N−45, so no qualifying future final-zero suffix can occur. The prefix bound is 122248. All arithmetic uses exact Python integers.

At N=100000000, independent results are:
- 916809 normalized k=0 representations.
- 75076 distinct positive odd targets through 10000000.
- 907984 distinct positive odd targets through 100000000.
- All ten reported largest internal gaps match. Largest: 854 missing odds between attained 987933 and 989643; second: 754 between 77196595 and 77198105.
- The independently produced 6250000-byte union bitset is byte-for-byte identical to the producer's UNION.bin. Encoding: bit floor(n/2), least-significant bit first in each byte. SHA-256: 871cca7430a860219c3658c8672a18728cbecd01989f28e0ad3af4431bb99c77.

`controls.py` independently counts Hamilton paths by recursive outgoing-path extension, rather than the producer's subset-endpoint iteration. All 129 graphs with total context length at most four pass, including nonnormalized prefixes and suffixes. Leading-zero/trailing-one normalization passes. The closed scalar formula is separately compared to a literal pair-sum definition of f for every binary context and split through total length eight.

`verify_receipts.py` independently reconstructs the seven-state core (45,1440,878,577,332,25,14), verifies all 180 retained witnesses (including all 160 k=1 witnesses) by actual absorption and transfer, checks order, verifies membership in the independent union, verifies the bare formula through k=3, and checks all four endpoint examples. RESULT.json and RESULT.normal.json are byte-identical.

## Proof audit

Deleting a leading unselected prefix vertex deletes a global source; deleting a trailing selected suffix vertex deletes a global sink. Every Hamilton path uses those vertices first or last respectively, giving a bijection. Iteration gives precisely the specified normal form. The nonempty cut and complement in the strongly connected core, together with the selected first prefix and unselected final suffix when present, make each normalized graph strong.

The scalar context formula has the correct empty-suffix convention c_0=1 and final selected-suffix correction. Since T counts the subfamily of A-covers having the additional selected terminal condition, T≤A. For a normalized suffix ending zero, F=2A_before−T_before+45≥A_before+45. Earlier suffix letters at least double A, proving the stated b and P bounds. For an empty suffix the singleton pair count satisfies f(P)≤P, yielding the stated upper bound on P. The denominator is strictly positive: q−G−45=818+3513600J.

The bare target lower bound follows from deterministic Hamilton-path insertion when vertices are added. The producer's length loop includes every length permitted by 2^(b−1)A≤N−45; exact filtering removes its harmless superset. At this cap only k=0 and k=1 pass the bare lower bound. Python integers avoid arithmetic overflow. Explicit check() raises, rather than optimization-disabled assert, protect the producer controls.

Absorption X^j into the prefix and Y^(k−j) into the suffix preserves the entire ordered, cut-labelled tournament, not just its count. Thus the all-k union equals the k=0 arbitrary-context union at every scale. The interval decomposition is additional organization, not additional occupied targets. Endpoint corrections are −20(P mod 2) at the first step and +14(1−c_b) at the final step; for k=1 they occur together. Therefore full unit spacing holds exactly for even P and c_b=1. All hostile numerical examples match.

## Preregistration and scope

The supplied preregistration specifies the fixed-core context universe, exact target-dependent bounds, positive and hostile controls, distinct-target reporting, and absence of a global spectrum claim. It transparently retains and corrects the erroneous proposed suffix-zero positive control. Its filesystem modification time precedes RUN.log's modification time, consistent with the stated order. Those mutable files/timestamps are not an independent immutable record of when the hypothesis or correction was first written. The current source was modified after the initial run, and later receipt/bitset additions were disclosed. Accordingly the audit verifies the present mathematics and reproducible outputs, not historical preregistration priority.

The percentages are finite counts among all positive odd numbers up to each cap. They neither prove sparse asymptotic density nor exclude another core, arbitrary surrounding blocks, or unrestricted tournaments.

## Optimization-safe gate repair

Integration identified that the initial independent scripts used Python `assert` for correctness checks. These have all been replaced with explicit exception gates. Each executable includes `--false-gate`, which deliberately invokes a failing gate before any work. All three fail with exit code 1 and the expected RuntimeError under both normal Python and `-O`. All three full normal/optimized runs produce identical output. The RESULT.json and CONTROLS.json receipts remain byte-identical to their preserved pre-repair snapshots, and the union still matches the producer byte-for-byte. GATE_REPAIR.json records these checks, and per-script logs retain the outcomes. Original source snapshots are retained as nonexecuted `.pre-gates.txt` files.

The receipt verifier first looks for sibling `../bare_coverage`, then falls back to `../session231-bare-coverage`. An isolated publication-layout test passed. The imported controls module emits no stdout; the receipt verifier now prints one JSON object. Its direct execution still prints its own control receipt.
