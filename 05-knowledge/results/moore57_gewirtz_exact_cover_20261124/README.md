# Moore(57,2) follow-up: sharpen the 56-match test

**Status (2026-11-24):** no unrestricted sharply transitive Gewirtz row or nonexistence certificate was found. The Moore graph remains **OPEN**. This follows the [transport and flag session](../moore57_switching_flags_20261124/README.md); now the working object is an exact permutation cover, not a degree estimate.

- **FINITE-EXACT:** A sharp 56-set in the full degree-56 Gewirtz group has trivial left and right multiplier groups. No order-3, 5, or 7 left–right symmetry survives. [Proof and orbit-cover replays](SYMMETRY.md).
- **PROVED conditional, using finite class inventories:** Every sharp set can be moved and conjugated to include the identity and a fixed element of one of **two** derangement classes. The exhaustive pair of search problems has 21,172 and 20,292 candidates (down from 40,920). [Exact reduction and optional solver](SEARCH.md).
- **PROVED with exact rank certificate:** The normalized target J−I is in the *signed integer* span of the derangement matrices; even without permitting nonderangements, modular weight obstructions alone cannot decide the sharp-set question. [Proof](gewirtz-obstruction.md).
- **PROVED reduction:** After a single sharp row is found, juxtaposing an independently chosen second row is not enough. The next test is one exact cover on *two* Latin tables sharing the same unknown contact at each branch. The generator and [Hoffman–Singleton calibration](moore-paired-branches.md) are included.

The symmetry tests are genuine restrictions on the candidate group, not an assumption that every graph must have these symmetries. The reusable lesson from tournament switching is to take the quotient on a pair of endpoints, not just on symbols or path counts. Here left and right endpoint relabelings move the entire set of 56 matchings, and an equivariant cover must tile **cell orbits**, not merely hit all positions in aggregate.

## Reproduce

Plain Python 3.10+ for the proof replays, GNU C++ for the involution clique pass, optional NumPy/SciPy for bounded unrestricted MILP. Run from a writable directory for the C2 data files.

    python3 odd-stabilizers.py
    python3 gewirtz-sharp-c7.py
    python3 gewirtz-sharp-c2-data.py
    g++ -O3 gewirtz-sharp-c2.cpp -o verify-c2
    ./verify-c2
    python3 gewirtz-obstruction-check.py
    python3 seeded-cover.py
    python3 moore-paired-branches.py

The next live check is one of the two reduced unrestricted 0–1 covers. Do not infer infeasibility from a timeout, a partial packing, a regular-subgroup obstruction, or the signed integer result.
