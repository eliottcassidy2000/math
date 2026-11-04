# Independent audit: prime-five quotients with triple children

**PASS.** The proof and exact finite certificate establish the stated fifteen-old-name class. No all-order, arbitrary ambient, or full reconstruction claim is inferred.

## Mathematical review

The maximal-child-run decomposition is a bijection for every block-uniform physical query. For each child, its successive nonempty runs are filled by an ordered path cover, in their occurrence order. Consequently the complete profile, rather than the internal labeling or only its Hamilton-path count, is sufficient. A three-name tournament has precisely one of the two complete profiles `(0,1,6,6)` and `(0,3,6,6)`. All eight actual labeled triples are independently classified here by permutation/cut enumeration. Thus all `8^5` physically labeled child choices reduce to all 32 profile patterns, with no cubic or quartic contributions removed.

We independently test all 1,024 five-name quotients for primeness and strong connectivity. Exactly 264 qualify. Explicit simultaneous vertex maps place each in exactly one of the three stated classes, with sizes 120, 120, and 24. For each map, the same map carries the opposite quotient to the opposite representative, all 32 profile patterns round-trip, and the actual physical query is transported to the corresponding child names. The compared tournaments are never independently relabeled. The receipt retains every map.

The three-representative reduction is therefore valid. All 96 profile comparisons have a nonzero integer difference at the advertised query; their minimum absolute differences agree with the producer. The assertion is exact integer separation at a single transported physical query, not a modulo-four statement or a retained full two-port coefficient table.

## Independent implementations and controls

`verify.py` does not import producer code. `count.cpp` supplies two independent recurrences:

- A weighted suffix recursion on run states, multiplying the first, second, and third visits of child i by `h_i`, `6/h_i`, and 1. These weights telescope to its complete cover count. A terminal state contributes one and can continue, accounting correctly for every possible number of runs.
- A literal backward subset/first-vertex Hamilton-path recurrence on the actual seventeen-vertex augmented tournament. It has no quotient or cover-profile assumptions.

All 192 response counts are checked by both methods. We additionally change the block labels and choose independently among the six transitive or two cyclic actual child orientations, then repeat the literal physical counts at the transported query. This gives 384 physical counts across all 96 patterns, not merely the producer's eighteen spot controls. All eight actual labeled child orientations occur. The audit reconstructs all 192 multilinear profile-polynomial coefficients independently by interpolation on the `h_i` values 1 and 3.

The C++ unsigned 64-bit recurrences are exact: their counts are bounded by `17! = 355687428096000`, below `2^64`. State transitions strictly increase the number of visited symbols, so the suffix recursion is acyclic. Explicit Python checks and C++ exceptions remain active with Python optimization.

Normal and optimized audit runs agree byte for byte. Fresh producer normal and optimized replays also match the frozen producer receipts byte for byte.

## Reproduction and pins

Run from this directory:

    python3 verify.py > normal.json
    python3 -O verify.py > optimized.json
    cmp normal.json optimized.json

An already installed `g++` is required; the executable is compiled in a temporary directory. No installation is performed. The source receipt is used only as the expected certificate and to obtain its three representative row tuples. All classifications, count recurrences, relabelings, and interpolations are independently implemented.

Frozen producer proof: `971a34905a1703fa2b9e9d6d64ea3718d366c05f1757389b87739ca03decf200`.
Producer verifier: `cf6af36e2405ce356a4e64d5e810b8793a1fcba16aa52b57d5b9084b3011e094`.
Both producer receipts: `c1df3507138216ddf95147350c9bb64a5a0fe9040770534e08f812d692ba0a75`.
Full source pins and all 96 audit records are retained in both audit receipts.

Publication allowlist: `AUDIT.md`, `verify.py`, `count.cpp`, `normal.json`, `optimized.json`, and `SHA256SUMS`. The identical producer replay copies need not be published a second time.
