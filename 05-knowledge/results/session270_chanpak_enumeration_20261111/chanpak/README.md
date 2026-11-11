# Compile Hamilton-path targets through one million

**FINITE-EXACT:** Every positive odd target≤1,000,000 except7 and21 has a literal tournament witness on at most16 vertices. The family is the inherited marked two-apex/transitive-spine grammar G(2,0), not a new family. The compact atlas contains499,998 witnesses. All were replayed;626 received independent subset-DP checks. Sizes are minimum within this fixed-core word language, not global minima. Strongness is not required.

    python3 compile_target.py 999999

This returns a15-vertex adjacency witness and its color word. Load/decompress the atlas once for repeated direct-address requests. The finite compiler is fast; the infinite Chan–Pak problem remains OPEN.

**PROVED:** A36-state nonnegative transfer supports arbitrary recursive color words. Exact invariant closure yields a minimal17-dimensional rational representation and a nonzero Hankel minor certificate, upgrading the inherited empirical rank17. See `rank.json`.

**PROVED:** For any fixed k-core C, let A_m sum Hamilton counts over all2^(km) cross orientations to an m-vertex transitive spine. With y=exp(2^k z), its EGF is

    Σ_m A_m z^m/m! = (y+1)^2/4 · Σ_r q_r(C)((y−1)/4)^(r−1).

It implies a small rational ordinary GF, constant-coefficient recurrences, and triangular recovery of the cover vector. Two-core mass is(12^m+3·8^m+3·4^m)/8 for m≥1. This is mass over presentations, not support or an unmarked tournament census. Averaging trades an exponential number of orientations for smaller individual sensor graphs.

`PROOF.md` provides exact statements, source lineage, limitations and reproduction commands. `atlas.json` gives binary layout, checksum and minimum-order histogram; `atlas.b85` is deterministic zlib/Base85 of2MB packed data. `census14.progress` records the length-by-length coverage frontier. `SOURCE_PINS.json` pins inherited mechanisms. No changes outside this directory are part of this lane.
