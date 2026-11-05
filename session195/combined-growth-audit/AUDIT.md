# Independent audit of the combined four-coordinate growth grammar

Session 195, 2026-11-05. **PASS after two presentation corrections, both resolved.** The full `{7,21}` Hamiltonian-path spectrum remains open.

## Frozen source and correction status

`SOURCE_SHA256SUMS` identifies the reviewed producer snapshot. The corrected proof SHA256 is

    740979ceedc7621232d5b668662f871bfb24a397f3f9c708a764445ab3fbe8f6

Two issues found during review were corrected before this snapshot:

1. The displayed rank matrix has continuations as rows and inputs as columns. Its original caption called it input-by-continuation. The corrected caption now agrees with the entries; the determinant remains -16.
2. The justification of `R>=H` now uses the valid injection obtained by splitting each Hamiltonian word immediately before the root. Not every arbitrary two-path cover can be joined into a Hamiltonian word.

Neither correction changes a recurrence, theorem conclusion, or computational result. No unresolved mathematical correction remains.

## Mathematical audit

- **Physical A and K births:** the new vertex has precisely the old root as an old out-neighbor. The path split gives `h'=h+R-e`. Changing the mark to the new vertex gives A; retaining the old mark forces the final endpoint segment x,r and gives K. The K two-cover partition correctly removes exactly d covers with singleton rooted second path.
- **Neutral sink D:** the old-root endpoint count becomes zero, and the new sink can be appended to either path of a rooted cover, including its previously empty first path. Thus R doubles without changing h or d.
- **General gadget:** an actual reversed seam is the only possible gadget-to-input transition. The seam/no-seam partition gives all three rooted statistics without identifying fragments or losing empty-bin cases. Deleting the new root restores an ordered join, giving `d'=H(G-b)h`.
- **Physical cone and positive coordinates:** `0<=e<=d<=h<=R` and `R+e>=h+d` are supported by explicit injections. All A/K/D transforms and all four general-gadget coordinates agree algebraically with the original state. In particular, deleting the root from an L-cover proves `B>=L`; this is the required nonnegative fourth-coordinate coefficient.
- **Four-dimensional linear minimality:** the actual continuation matrix has determinant -16. This excludes lower-dimensional homogeneous linear continuation signatures for this grammar. It does not exclude nonlinear encodings or assert a fixed-dimensional state for arbitrary colored gluing contexts.
- **Arbitrary-depth bounded completeness:** R is nondecreasing, every non-D operation outputs h at least the input R, and D leaves h unchanged. Therefore a prefix with R above the target can contribute only its already recorded current h, followed by neutral D operations. Terminal outputs with h within the target but R above it are explicitly retained as values. Deduplicating equal four-coordinate states preserves every continuation. The finite saturation is complete without any graph-order or depth cutoff.
- **K-spine and gcd statements:** the closed form follows by induction from K and then one A. Every generator is integer-linear, so coordinate gcd divides every descendant count. These statements constrain this marked grammar, not all witnesses for the same scalar H.

The proof accurately separates structural interface closure, finite exact support, and the unresolved all-order spectrum question.

## Independent implementation controls

`verify_independent.py` was written separately and imports no producer code. It computes ending and starting subset tables, constructs actual graph operations, and regenerates the full gadget signature set.

It verifies:

- all 1,099 labelled tournaments of orders one through five and all 5,405 choices of root;
- exact equality with the supplied set of 65 **eight-entry** signatures, including the root-deletion coordinate;
- 855 physical A/K/D transitions on every rooted labelled input through order four;
- 1,885 physical gadget transitions on every rooted labelled input through order three and every inventory gadget;
- the old-triple collision and its distinct H55/H59 continuations;
- the determinant -16 continuation matrix and the H9 coordinate-gcd hostile;
- the explicit singleton-to-H163 and singleton-to-H247 derivations, every intermediate state, final labelled rows, and final root.

The two final witness states were additionally checked by literal permutations and literal rooted path cuts, independently of subset dynamic programming:

    H163: (163,70,168,75), order 8, root 7;
    H247: (247,112,250,115), order 9, root 8.

The independent tuple-based Python saturation at target 1,000 gives exactly 21,235 retained states and misses precisely 7,21,67,173. Normal and optimized Python runs have byte-identical stdout and receipts.

## Full census replay and code inspection

The frozen C++ saturation code was inspected and compiled independently. Its packing is safe because it enforces a target at most 65,535 and packs a state only after all four coordinates are bounded by R and R by the target. Arithmetic before packing is signed 64-bit. The parser's actual inventory is hash-bound and matches the independently regenerated physical signature set.

Replays establish:

- Combined bank, target 1,000: 21,235 retained states, agreeing with independent Python.
- Combined bank, target 20,000: 2,559,611 retained states, 9,996 attained odd counts, exactly 7,21,67,173 missing.
- A/K/D alone, target 32,768: 2,046,988 retained states and exactly the 29 omissions listed in the producer proof.

The full high-bound runs use the reviewed producer C++ implementation, while the lower-bound tuple saturation and physical/inventory checks use independent code. The all-depth completeness follows from the proof above, not an unexplained finite graph census.

## Reproduction and release

In the published layout, run from `session195/combined-growth-audit`:

    python3 verify_independent.py --source ../combined-growth --bound 1000 > normal.out
    python3 -O verify_independent.py --source ../combined-growth --bound 1000 > optimized.out
    cmp normal.out optimized.out

For the C++ replay, compile `../combined-growth/census.cpp` and pass its `inventory.txt` at the documented targets. The local audit used an exact frozen copy of those sources.

Publish `AUDIT.md`, `verify_independent.py`, `receipt.json`, `SOURCE_SHA256SUMS`, `cpp-1000.json`, `cpp-20000.json`, `cpp-near-32768.json`, and `SHA256SUMS`. Frozen source copies, compiled binaries, and stdout copies are local audit support. Supplemental order-six gadgets are outside this baseline result and must retain their own audit scope.
