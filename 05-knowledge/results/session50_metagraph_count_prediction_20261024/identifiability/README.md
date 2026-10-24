# Session 50: metagraph count identifiability

This standalone packet separates exact count reconstruction, arithmetic extraction, and a genuine next-order prediction requirement.

Main new results:

1. A necessary-and-sufficient realization theorem for the abstract odd-fibre partition layer on the **actual** complement/reflection cube.
2. A bounded n=6 hostile with native totals, native local signature palette, native loop totals and connected layers, but S=12 versus20 at equal V,N,B,M2. It is explicitly non-native, not a conflicting tournament census.
3. Same-uncolored-support countermodels matching every fixed finite degree-moment truncation in the abstract relaxation.
4. Merged-degree congruences extracting S modulo8 from M2, and modulo2^(j+2) from M_(2^j).

Read [THEOREM.md](THEOREM.md) for proofs and exact scope. Historical parity and black-twin results are credited rather than rediscovered. `SOURCE_PINS.json` records exact upstream commit, blob, byte length and SHA-256; historical Markdown is archived unchanged with `.txt` suffixes.

Reproduce using only the standard Python library:

    python native_reference.py --check
    python reproduce.py --check
    python verify.py
    python -O verify.py
    python verify_manifest.py

The normal and optimized verification reports are identical. `PORTABILITY.json` records replay after copying only manifest-listed files into an unrelated temporary directory. No Git, repository checkout, network, package installation or neighboring session is needed for replay.

`native_reference.py` independently enumerates fixed-path words through n=6 and canonicalizes over all permutations within equal-score blocks. The resulting native merged-degree histograms agree with the archived independent historical census. The model checker inspects all1024 words of each saved n=6 hostile without importing the producer. It also emits concrete score-multiset witnesses showing why each hostile partition is not native.

The exploratory integer-program search is not a verification dependency. Its bounded witness is stored as explicit integer multiplicities in `reproduce.py` and checked exactly. The finite higher-moment identities use rational Lagrange interpolation, not floating-point fitting.
