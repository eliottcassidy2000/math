# Publication status and audit navigation

Read AUDIT.md for the independent verdict and final supplement checks. THEOREM.md is a clearly marked corrected publication edition of the audited AP-only snapshot; its original false suggestion of repeated path realizations across frames is removed. The original producer README below is retained as historical context, not as a separate final verdict. All-frame results are in bare-extension/THEOREM.md. SOURCE_SHA256SUMS.txt preserves the historical manifest; SHA256SUMS describes this publication edition.

---

# One affine square plane, exactly 6(q−9) disjoint AP frames

Read THEOREM.md for definitions, proof, representative independence, all factors, and the non-AP exclusion. This is an exact characteristic-three counting result conditional on a fixed square-valued affine plane, not an integer existence theorem or a bare-biclique census.

Reproduce using the standard library:

    python verify.py
    python -O verify.py

Both runs check q=81,729,2187 witnesses and generate results.json. The F81 pass independently enumerates every unordered pair of ambient affine lines and checks all 31,104 fully ordered frame-to-magic-layout images, each of the 432 layouts having exactly 72 valid gauges. Runtime requirements use explicit exceptions and survive optimization.

Source pins in sources/ are copied local Session 57 mathematical packets, not mutable imports. SOURCE_PINS.json records source origin and SHA256; SHA256SUMS records the whole deliverable. No old source, repository, or remote was modified.
