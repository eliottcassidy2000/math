# Independent audit: induced pruning and the all-order H541 obstruction

Status: **PASS, no correction required.** Source files remain unchanged.

## Frozen source binding

All source files are in `session197-hybrid-constructions/induced-pruning/`:

- `PROOF.md`: `a5ec2d6842289f75e833881e4bb1c5e0ab715c38446a3486b5119b06c64bdee0`.
- `core.py`: `3204eeb1e9d2b05cac9bcc21c98f25697d10f925d6ff743393e945ce4e73615c`.
- `pruning.py`: `1869476a5c45595da0b434f5c3eb6e4e832f82ab857bc3ae9f842c12d242f409`.
- `census.cpp`: `d6a6ff83b0b96c451ee22e656bfa2b84745734998732f215d881f0eb5c944751`.
- `verify.py`: `72e20c2e050523d5144f92112358ceda68efc68c2c791c16549560eda66ad1c6`.
- `VALIDATION.json` and `VALIDATION_OPTIMIZED.json`: `216da0765f2ced584da9e3273dfbe50ab80d5c3e19182247a2db18776af411df`.

## Proof audit

The induced-pruning argument is valid for the exact full language stated. A deleted original port can be represented by a kept physical vertex that is no longer an original port. The representative was continuously exposed in the simulation, so this does not re-expose a forgotten vertex. It must receive default incidence at later kept births. The explicit translation of masks and forgetting indices through the current simulated port list is essential and is implemented correctly. Original and simulated list positions need not agree after skipped births. The invariant proves every retained final original port remains exposed. The empty bijection handles k=0.

The SCC reduction is also valid for arbitrary nonempty seeds: deleting only new singleton final strong components preserves H and leaves the whole seed. Nontrivial component counts are odd and at least three. Their number is at most Omega(N), and multiplying Busch's inequalities gives a total nontrivial order at most Omega(N)+floor(log_5(N^3)). For a nontrivial strong seed, its whole vertex set lies in one such component, allowing subtraction of the seed order. This addresses arbitrary singleton padding before applying an order bound.

The decision theorem covers all finite words, including adaptive choices, because any successful word has a shortened witness. The fixed-k alphabet and logarithmic depth give a polynomial in the numeric N, not its bit length. The stated SCC-based capped counting procedure also makes the seed-order dependence polynomial. The even-target and N=1 cases are correctly separated.

For the H179 seed with exactly ports [0,1,4], primality of 541 forces one nontrivial final SCC. The integer inequality 5^11 <= 541^3 < 5^12 restricts that component to order at most 12, so at most four kept births suffice. Induced Hamilton-count monotonicity makes every H>541 prefix rejection safe. The theorem does not make 541 an unrestricted missing tournament count or apply to other seed interfaces.

## Source replay

Normal Python and `python -O` runs both reproduced the frozen validation receipt byte for byte. This rechecks 193,536 exhaustive short pruning controls, 3,000 long controls, 200 SCC/order-bound controls, the independent Python depth-two census, and the full source C++ census. The source's early return when an induced subset has more than 541 Hamilton paths is sound by the same monotonicity.

## Independent full census

`check_pruning_census.cpp` was written independently for this audit. It reconstructs the physical births directly, visits the operations in a different order, and computes scalar H by the odd-cycle-collection formula rather than either Hamilton endpoint recurrence.

For each graph, it enumerates each directed odd cycle once by fixing its least vertex. If more than 270 odd cycles are present, H is already at least 543. Otherwise it computes weighted disjoint cycle collections using

    F(S) = F(S-v) + 2 * sum_{C contains v, C subset S} multiplicity(C) F(S-C),
    F(empty)=1.

All terms are nonnegative, so saturation above 541 is safe. Exact adjacency arrays are cached, but word prefixes are not deduplicated. The checker tests all 99,088 word-edge graphs, representing 12,545 distinct adjacency arrays, and reproduces every histogram in the source receipt. Accepted prefix counts are exactly 1, 36, 440, 5,716 and 79,768. No target 541 occurs. Optimized and unoptimized library-checking builds produce identical receipts.

Independent checker `check_pruning_census.cpp` SHA256: `1d001e5fa49089e7cc57f228551ff816527587a75e3d09436c0c2831a2d9371e`.
Independent receipt `pruning-independent.json` SHA256: `fba396469949fd375de6de67dde8683c4d82c492c9429490bbe8d8142cb3de67`.

Reproduction:

    g++ -O2 -std=c++17 check_pruning_census.cpp -o /tmp/check_pruning_census
    /tmp/check_pruning_census

The proof establishes the all-order reduction; the complete independent finite census supplies the required bounded obstruction. Both parts pass.
