# Session223 3-adic research lane

Main result: **ALL_K.md** proves the full fixed-length two-J valuation law at every precision 3^k. **PROOF.md** gives the shorter mod-3 derivation and the length-composition resolvent argument. Neither claims the height-truncated M_2 factor.

Independent executable controls:

- `python3 verify.py`: exact mod-3 projection plus literal enumeration through length9 and integer transfer through length1000.
- `python3 exact_higher.py`: full finite modular component graphs at q=3,9,27, exact rational projection, comparison to the analytic formula.
- `python3 exact_81.py`: same check at q=81, with1,056,321 normalized states.
- Repeating each command with `python3 -O` produces identical receipts. All checks use explicit active exceptions rather than disabled assertions.

Receipts: VERIFY.txt / VERIFY_OPT.txt; HIGHER.txt / HIGHER_OPT.txt; HIGHER81.txt / HIGHER81_OPT.txt.

Independent proof reviews and independent graph-based controls are in `../audit/`. Exploratory scripts are not proof dependencies.
