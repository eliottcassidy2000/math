# Session 133: independently replayable low-H tournament atlas

**Checkpoint 5811bb4 correction:** the H15 polynomial hostile was a trailing-
zero comparison error. The corrected H9 control normalizes coefficients. See
[CORRECTION.md](CORRECTION.md); all atlas counts are unchanged.

`ATLAS_PROOF.md` states scope, completeness proof, primary citations and
inheritance. `build_atlas.py` regenerates `atlas.json` with every ordinary
class through order six, actual adjacency, H, automorphism order, SCC sizes,
literal cycle-packing profile and self-converse flag. `atlas.optimized.json`
is byte-identical. The low-count strong atom inventory is complete for H<=21
at every order by Busch's bound.

This atlas is a verified input for new FRACTRAN sequence compilers, not a new
claim that 7 and 21 are excluded or that all other odd counts are attained.
The latter completeness conjecture remains open.
