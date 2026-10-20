# Generator-free circulant naturality

20 October 2026. **Four scoped theorems; independent audit PASS.**

The [proof](GENERATOR_FREE_NATURALITY.md) provides a generator-free augmentation projection for prime orders p≡3 mod4, an exact abstract-cyclic naturality obstruction for p≡1 mod4 with p>5, a Frobenius obstruction in the even-degree finite-field multiplier-type model, and a first-order nondefinability result for this particular projection.

These are different notions: a uniform arithmetic algorithm is not the same as a map natural after forgetting field addition, and naturality is not first-order definability. The prime-power model is not automatically the full sign-even graph problem; the q=9 rook graph gives an explicit hostile.

The [independent audit](review/INDEPENDENT_AUDIT.md) confirms all four statements, checks the circuit lower-bound source, and supplies independently implemented exact controls. The producer's 183,074 checks and the independent 176,751 checks agree in normal and optimized Python. Frozen final proof hash is recorded in both checksum manifests.

Run `python3 verify_generator_free_naturality.py --output verification.json`, then repeat with `python3 -O`. Independent controls are in review/. No historical priority claim or universal nondefinability claim for every prime-field bijection is made.
