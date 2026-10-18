# Pre-registration (S44, third round): out-of-sample test of THM-4710's closed forms

Written and committed BEFORE any computation at the values below. The closed forms were fitted on N = 9..48
(`witt_rigidity_interpolation_9_48.json`).

**Predictions.**
- **P9.** For every N = 49, ..., 64, the forced degree-N generators (unique admissible lowest-weight vectors,
  generators at Li 3, 5, 7, 9) satisfy the following, in Cayley coordinates:
  - `G_5 = A_5 - ((5N-1)/(4N+1)) A_3 A_2`;
  - `G_7 = A_7 - (28N-33)/(6(4N+1)) A_5A_2 - (48N^3+224N^2+185N+123)/(2(N+4)(4N+1)(4N+3)) A_4A_3
          + (560N^3+1856N^2-1419N-477)/(24(N+4)(4N+1)(4N+3)) A_3A_2^2`;
  - the `A_7 A_2` coefficient of `G_9` is `-3(3N-8)/(2(4N+1))`;
  - the other six `G_9` coefficients are the rational functions recorded in `witt_rigidity_interpolation_9_48.json`;
  - `dim S_L = 1` at L = 3, 5, 7, 9.
- **P10.** Vir_(4,13) (N = 18, c = -217/26, not used anywhere before) has, as its degree-18 generator at Li degree 5,
  `G_5(18) = A_5 - (89/73) A_3 A_2`, and has no degree-18 generator at Li 7. Its generators sit at Li 0, 3, 5.

**Runs.**
- `witt_rigidity_interpolation_20261018.py 49 64`, compared with the formulas by `check_witt_closed_forms_20261018.py`.
- `semiinvariant_generators`-style extraction for (4,13).
