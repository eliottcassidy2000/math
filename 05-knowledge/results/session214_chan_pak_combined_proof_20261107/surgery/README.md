# One-exit surgery: obstructions and one-arc escapes

**Independent audit: PASS.** The [class-exclusion audit](../surgery-audit/SURGERY_AUDIT.md) and [one-arc construction audit](../surgery-audit/ONE_FLIP_AUDIT.md) approve the proofs and finite certificates below. These signed reports supersede the provisional "pending independent audit" labels preserved in the frozen proof files.

The [initial proof](PROOF.md) shows that the least odd omission beyond 7 and 21 in the class "transitive chain plus at most two apexes, one apex beating at most one chain vertex" is **7,477**, at every order. This class permits relocating the binary comb's unique exit by reversing old arcs, so it goes beyond Session213's induced-retention barrier.

The [addendum](ADDENDUM.md) classifies the case of two nontrivial strong components and proves the next all-order class omission, **27,371**. Through 27,371, the only missing positive odd values are 7, 21, 7,477 and 27,371. Arbitrary dominance products of this class still omit prime 7,477, but recover 27,371 as 101 times 271.

The [positive surgery controls](ONE_FLIP_CONTROLS.md) attain both new omitted values, 7,477 and 27,371, by **one old-arc reversal of an actual Session213 comb**, retaining its first exit and adding a second:

- B^-(8,14,0): H=3,187; reverse x_0->a to obtain H=7,477.
- B^+(9,179,1): H=24,693; reverse x_5->a to obtain H=27,371.

One is the exact minimum number of old-arc reversals when choosing any original comb and keeping its vertex set. This does not cover additions, deletions or arbitrary rebuilding, and these two controls are not an all-target recurrence.

The unrestricted Chan–Pak spectrum remains **OPEN**. The exclusions are class-relative, not ambient forbidden values; the larger two-exit class is not proved complete.

## Verification

Author normal and optimized runs agree. The independent audit uses a separate path-cover convolution census for the class exclusions and literal Hamiltonian-word enumeration for the one-arc controls. The cutoff proofs, not finite extrapolation, give the all-order exclusions.

`counting.py` contains the run transfer, quadratic cut interpolation, Gray traversal and direct subset counter. The transfer and quadraticity mechanisms are inherited and explicitly credited. Source pins are in `SOURCES.json`, `ADDENDUM_SOURCES.json` and `ONE_FLIP_SOURCES.json`.

## Reproduction

Run from this directory with Python 3; the author checkers use only the standard library:

    python3 verify.py > /tmp/surgery-7477.json
    cmp RESULT.json /tmp/surgery-7477.json
    python3 -O verify.py > /tmp/surgery-7477-opt.json
    cmp RESULT.json /tmp/surgery-7477-opt.json

    python3 addendum_verify.py > /tmp/surgery-27371.json
    cmp ADDENDUM_RESULT.json /tmp/surgery-27371.json
    python3 -O addendum_verify.py > /tmp/surgery-27371-opt.json
    cmp ADDENDUM_RESULT.json /tmp/surgery-27371-opt.json

    python3 one_flip_verify.py > /tmp/surgery-one-flip.json
    cmp ONE_FLIP_RESULT.json /tmp/surgery-one-flip.json
    python3 -O one_flip_verify.py > /tmp/surgery-one-flip-opt.json
    cmp ONE_FLIP_RESULT.json /tmp/surgery-one-flip-opt.json

The signed audit reports give the independent C++ census and Python literal-enumeration commands. `SHA256SUMS` and `PUBLICATION_FILES.json` pin the fourteen release payload files.
