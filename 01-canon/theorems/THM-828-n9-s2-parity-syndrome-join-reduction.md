---
id: THM-828
title: The n=9 raw-S2 collision join is supported on a 14-dimensional XOR syndrome code
status: RESERVED + PROVED LOCAL SYNDROME REDUCTION; exact R8 join in progress
source: codex-2026-07-15-S13/referee
depends_on: [THM-809, THM-818, THM-825]
related: [THM-801, THM-814, HYP-6880]
---

# THM-828 — the `n=9` raw-S2 join has only `2^14` possible differences

This number reserves the exact XOR-syndrome reduction now driving the full
THM-818 computation.  For an apex-zero upper pair `(u,v)`, put `D=u xor v`.
Equality of every raw four-state `S2` histogram forces even `D` parity on the
low and high halves of each of the six nonfixed layers, even parity on the
fixed layer, and `D_apex=0`.  These thirteen parity equations plus the apex
equation have independent rank fourteen in the 28 upper bits, leaving exactly
`2^14=16,384` possible differences, including zero.

Because every nonfixed layer has size at most three and the fixed layer has
only three free positions, these parity conditions are also existence-exact:
every syndrome word is realized by at least one pair of layer words with the
same raw histogram.  Actual metagraph collisions still require the three face
kernel conditions, upper colour, and equality for the particular glued base;
the full join and its survivor census are not claimed complete in this stub.

The final theorem will contain the elementary `F_2^2` cycle proof, an explicit
rank-fourteen basis, the closed `16,441,671,680` count of all apex-zero ordered
raw-S2-equal pairs before face constraints, the difference-indexed `R_8`
algorithm, exact witnesses, and the resulting collision genealogy.
