# Independent self-converse supplement audit

**PASS.** The separate Session197 self-converse supplement is correct as stated. This report covers the explicit eleven- and fourteen-vertex graphs, their exact Hamilton counts and automorphism groups, the reflection-fixed word counts, and the sharp minimum **order** of eleven in the self-converse nonrigid `7p` sector. It does not assert a minimum Hamilton count or an unrestricted spectrum theorem.

## Frozen inputs

All inputs are in `session197-orbit-wildcard/`:

- `SELF_CONVERSE.md`: `79c2fa74478e485df184c7b4b9831b08f14c67008d42bcea5c31afec0ff553a9`.
- `verify_self_converse.py`: `6af2f980b72d638e5cb167035fdf9842e5d1643f8db93d88a53508da1f099698`.
- `self-converse-normal.json`: `c35d8ba1b4da4157d5521ca545be0010cd68b1d73e85803762c5440347beedd7`.
- `self-converse-witness.json`: `6c52d4b98277eb0e83927ba33ac0bab2b5e508f17e38894d07704c080ab3df17`.

The original fourteen-vertex packet was separately audited in `REPORT.md`. Its source files are unchanged.

## Independent computations

`check_self_converse.py` reconstructs the eleven-vertex graph directly from the five-vertex quotient and the carousel, rather than importing the author's construction or checker. It imports only the independent graph-audit utilities. The resulting outgoing rows are

    [398,412,440,496,481,451,391,1536,1664,127,639].

Walk inclusion-exclusion, rather than either of the author's subset endpoint recurrences, gives `H=56483=7*8069`. Trial division through 89 and a recursively checked Lucas certificate both prove 8069 prime. Literal quotient-run enumeration gives `(11,17,4)` and independent ordered-cover enumeration gives `(175,1638,6678)`, reproducing the substitution count. Exact automorphism backtracking gives precisely the seven carousel rotations, with 82 search nodes. Strong connectivity and every arc of the displayed anti-involution are checked directly.

For both graphs the reflection-fixed Hamilton words are counted in two independent ways: a dynamic program on selected involution pairs, and inclusion-exclusion over pair-colored walks. These give **111** fixed words at order eleven and **6111** at order fourteen. The center condition for the odd graph and the joining edge for the even graph are imposed explicitly. The displayed witness words are checked literally. The fourteen-vertex dihedral relation and Burnside count give exactly **1084789** word orbits.

The original fourteen-vertex count `15144269`, primality of its quotient `2163467`, exact automorphism group of order seven, and module-indecomposability are already independently certified by the primary graph audit. This supplement does not substitute a fitted recurrence or a residue test for any graph construction.

## Minimum-order proof review

The cited S196 `orbit-spectrum/PROOF.md` at commit `cf3b44b9a442c3cff0f28cab924e73c4c00bb00c`, Sections 5–6, was read directly. Its inherited minimum order ten and reduction to full automorphism group of order seven have the needed scope, including the exceptional small prime cases.

At order ten a group of order seven has one orbit of size seven and three fixed vertices. An anti-automorphism normalizes the automorphism group, so it preserves both of these odd-sized sets. A nontrivial self-converse tournament has an anti-involution: the automorphism group has odd order, hence an anti-automorphism of order `2m`, with `m` odd, has an odd power of order two that still reverses arcs. Such an involution has a fixed vertex in each preserved odd-sized set. An anti-automorphism cannot fix two vertices because it would reverse their connecting arc while fixing its endpoints. This rules out order ten, and the explicit order-eleven graph attains the bound.

The argument is used only for a graph of order ten; the harmless degenerate one-vertex convention for “anti-involution” does not affect this proof. No source correction is required.

## Reproduction

Run `python check_self_converse.py` and `python -O check_self_converse.py` in this audit directory. The normal and optimized receipts agree byte for byte.

- Independent checker SHA-256: `22c2a0103c97a62180edf1bd406b8e21a5620902c920fb4d5afd25436499f2f3`.
- Receipt SHA-256: `1c44ecc311ad2b597bdc7f5f13c01d4a9783bb935f02c1fe47a011dd530a6866`.

No source file was changed during the audit.
