# Odd-transitive kernel span and local noncancellation

Generic local module-converse detectability is already incoming E4. This packet adds a restricted explicit mixed-quartic battery, the full span/jet description, and a constructive matching witness; it does not claim the first local noncancellation theorem.

The proof packet establishes:

1. The selected odd-transitive mixed-quartic kernels span all polynomials in `t(1+t)` through the expected degree.
2. For n physical exterior names, their response functional has only `floor((n+1)/2)+1` relevant Taylor moments, with an explicit inverse-jet certificate.
3. A maximum matching in the directed cross-edge graph supplies an all-order signed nonzero moment for every nonempty physical exterior.
4. Consequently, if `m >= 2 floor((n+1)/2)+3`, reversing an odd transitive module A of size m is detected by a named mixed quartic. Arbitrary groups contained in one incidence side cannot hide all of the detection.
5. Mixing the two incidence sides can cancel all these moments: the contracted directed triangle gives a minimal exact counterexample to unrestricted grouping.

This is a **local-A-converse theorem**. It is not an unrestricted collision theorem, not a quotient-only root theorem without restoration, and not an arbitrary-fixed-ambient source-family result. The domination condition counts physical exterior names, not quotient vertices.

The new proof uses the exact two-channel local-converse identity from `../session189-odd-transitive/PROOF.md`. The finite controls do not replace the proof. Independent review is pending.

## Reproduction

Run `python verify.py` from any directory. Standard library only. It verifies the kernel ranks/reflection for odd orders 3 through 31, exhausts all 33,866 physical exterior/root-incidence cases through exterior order 5, checks the exact matching formulas and moment inversion, and compares 42 literal physical quartic differences in seven cases with the derived formula.

`receipt.json` contains the exact controls and source hashes. `verify.stdout` records the normal run. An optimized run is recorded separately when available.
