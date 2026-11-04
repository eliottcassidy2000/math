# Independent review: transitive-module distance reconstruction

**PASS.** Reviewed the frozen Session189 collision packet whose `PROOF.md`
has SHA-256 `dcaad9182ddea6a9520a1fd15deaa9f32bb60b1f18e9e54dae52b2dc3619d824`.
`frozen_source_manifest.json` pins all fourteen producer files. The review
approves the corrected all-order statements and their stated boundaries.
The earlier extra factorial is not present in the reviewed proof or verifier.

## Mathematical review

1. **The prefactor is r+1, not (r+1)r!.** Cutting r excursion transitions
   and the distinguished two-cover boundary leaves r+2 increasing runs.
   The distinguished boundary has r+1 possible positions; the other run
   labels already name the remaining ordered slots. Inclusion-exclusion
   on those labels is exhaustive, and another factorial double-counts
   their placements. At top degree r=n-2, the corrected coefficient is
   `2(n-1)!`, matching the unrestricted permutation/cut count directly.

2. **The coefficientwise comparison is valid.** Consecutive subtraction
   gives the forward difference of
   `(s+1)^(n-d-2)((s+2)^(d-1)-s^(d-1))`. For d>=2 its ordinary coefficients
   are nonnegative and its constant is positive. The surjection identity
   for finite differences pays every coefficient; d=1 gives equality.

3. **The ambient transfer keeps joint supports and ordered excursions.**
   Each of the two endpoint-constrained paths meets the module. There is
   exactly one external prefix, one external suffix, and a collection
   of excursions between module names. Squarefree multiplication enforces
   disjointness. The factor `tau^r` assigns ordered excursions to the
   selected transition slots; it is compatible with the corrected run
   count and requires no additional division or multiplication by r!.
   `theta_0=H(R+z)>0` follows by cutting at z and is valid in arbitrary
   ambients, not only ordered or strongly connected ones.

4. **Reconstruction and the threshold are sound.** Summing deletion
   coordinates gives the stated S_d, with `S_1=S_2<S_3<...`. The unique
   maximizing pair is the pair of physical endpoints. Either endpoint's
   scores determine all other positions except distances one and two;
   the opposite endpoint distinguishes those at m>=5. For m<=4 every
   available distance after deletion is at most two, so this particular
   battery has exactly the stated failure boundary.

5. **The collision reduction retains its hypotheses.** I read the pinned
   Session188 hereditary-monotonicity, highest-node, normal-form and
   all-cover-restoration sources. Full response equality supplies scalar
   equality, hence complete hereditary ordered-cover equality on the
   original modular tree. Other highest nodes are restored at actual
   queries whose ports imitate the selected core. Saturated bidegree
   (2,2) prevents fixed exterior bits from adding aliases. Reversing each
   top-degree cover path is valid because both ports are internal and
   the port-port arc is unused. Other children can then be restored by
   their full unmarked profiles. The resulting comparison has the same
   exterior and differs only in the transitive child order. Applying the
   proved reconstruction theorem is therefore legitimate. This yields
   binary child choices, not exclusion of the remaining reversal cube.

An equivalent check of the top-degree reversal step is available directly
from the run proof: for an arbitrary module interior, its symmetric
endpoint-weighted two-cover polynomial is invariant under whole converse,
by reversing both internal paths and exchanging their slots. The same
nonnegative exterior excursion coefficients apply. This agrees with the
producer's more general all-cover transfer argument.

## Independent computations

The producer's normal and optimized verifier receipts were reproduced
in a separate directory and are byte-identical to each other and to the
frozen receipt. No producer file was edited during review.

The new `literal_bank.cpp` does **not** use Hamilton subset DP or Boolean
finite differences. It visits every permutation of the old names and the
two ports, tests literal old arcs and the fixed x->y edge, and expands
the actual port-contact factors. It independently reproduces all 44 full
named coefficient vectors, including every zero:

- Total coefficients: 38,888.
- Cubic coordinates: 13,392.
- Quartic coordinates: 20,784.
- Lower-degree coordinates: 4,712.

Every saved physical query value is independently reconstructed from those
literal coefficients. This strengthens the producer's two complete literal
controls to every graph in the retained bank. It is still that declared
finite bank, not a census of all tournaments of the displayed orders.

The producer's exterior universe of at most three names cannot exercise
two nonempty excursions: each excursion needs at least two vertices.
This is an evidence limitation, not a flaw in the all-order proof. The
review adds direct capacity controls with genuinely positive higher terms:

- m=5 with four exterior names: theta=(11,17,4).
- m=6 with six exterior names: theta=(51,229,216,36).

For each, the identity is checked on the original transitive order, whole
converse, and an internal adjacent swap. Theta is independently counted
by literal exterior permutations and ordered bins, not by the producer's
squarefree polynomial multiplication. All physical rows and shared-name
coordinates remain in `audit_receipt.json`.

Finally, direct shared-quartic counts verify the quotient-converse/profile
restoration on a nontransitive C3 quotient with child sizes 5,3,2, an actual
internal permutation of the first child, and conversed other children.
The checked whole orders are 10,12,13, including a non-ordered exterior.
These controls support the structural reduction; they are not asserted
to be examples of full-response equality.

The independent `audit.py` passes normally and under `python -O` with
identical stdout. It uses explicit exceptions, not removable assertions.

## Reproduction

Let SOURCE be the directory containing the fourteen frozen producer files.
From this audit directory run:

    g++ -std=c++17 -O3 literal_bank.cpp -o literal_bank
    ./literal_bank < bank_input.txt > literal_coefficients.txt
    python audit.py "$SOURCE" > audit.normal.stdout
    python -O audit.py "$SOURCE" > audit.optimized.stdout
    cmp audit.normal.stdout audit.optimized.stdout

The producer replay receipts are retained as `normal_receipt.json` and
`optimized_receipt.json`. To regenerate them without modifying the frozen
producer directory, copy its fourteen files to a temporary directory and
run its verifier there, directing stdout to the corresponding audit files.
The producer verifier writes deterministic data alongside itself.

Publishable audit files: `REVIEW.md`, `audit.py`, `literal_bank.cpp`,
`bank_input.txt`, `literal_coefficients.txt`, `frozen_source_manifest.json`,
`audit_receipt.json`, `normal_receipt.json`, `optimized_receipt.json`,
`audit.normal.stdout`, and `audit.optimized.stdout`.

No unrestricted reconstruction, all-pairs orbit injectivity, actual response
collision, or all-prime theorem has been inferred. The packet is approved
at its stated scope.
