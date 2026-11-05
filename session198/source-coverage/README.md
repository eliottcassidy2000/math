# Certified source-pair coverage and a sharper C3 inverse

The proof gives five exact scalar C3 forms instead of eight, physically sharp bilinear product constants, and the sharp common source cutoff U<=floor(2H/15), improving U<=floor(H/3).

Run `python3 verify.py` and `python3 -O verify.py` from this directory. These standard-library runs reproduce `VALIDATION.json`, `SUPPORT.json`, and `WITNESSES.json`.

At target bound 5,001:

- Strong C1 has 2,471 represented values; C2 has 2,454; C3 has 2,466.
- C1 union C2 misses exactly 7,21,67 among the odd values from3 through5,001.
- Adding C3 leaves exactly 7,21. This is a finite-exact support statement, not universal coverage.
- The C3 source bank shrinks from24,975 profiles to6,339; the candidate evaluations shrink from758,909 to94,591, with exactly the same scalar image.
- All126,946 source recipes in the C1 bank are checked by a separate recurrence implementation. Eighteen selected actual witness graphs are counted independently for both H and C.

`PROOF.md` contains the all-order arguments and the inherited dependency contract. `SUPPORT.json` contains explicit physical source recipes for each supported value in each class. `WITNESSES.json` contains selected actual adjacency rows, roots, and recipes. `VALIDATION_NORMAL.json` and `VALIDATION_OPTIMIZED.json` must match.

The universal representability goal remains open. Positive coordinate matrices and finite support do not imply that arbitrary nonnegative profiles are physically realizable or that arbitrary odd targets admit descent.
