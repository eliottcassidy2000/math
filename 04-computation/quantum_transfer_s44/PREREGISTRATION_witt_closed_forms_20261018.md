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

## Outcome (recorded after the runs): P9 and P10 PASSED

- **Timing.** The pre-registration was committed in 22048565d9 at 15:50:16 and pushed. The first run produced N = 49
  at 15:56.
- **Deviation in method, forced by memory.**
  - The pre-registered PBW-basis run passed at N = 49 (dims 1, 0 mismatches). It then died with MemoryError at
    N = 50. `vir_ss.block(2N+L, N)` enumerates every partition of `2N+L` into parts `>= 2`, about `10^8` at N = 50,
    and the stuck process was killed at 16 GB.
  - The test was completed by `check_witt_closed_forms_stable_20261018.py`. It is the same recursion in the faithful
    stable-range coordinates `mu`, partitions of `L`, where `d` and `D_n` are exact formulas in `N`; it is valid for
    `N >= L`.
  - That implementation was first validated: on N = 9..49 it reproduces the PBW-basis results exactly ("validation
    OK").
  - The predictions themselves were not changed.
- **P9 PASSED.** For every N = 49..64, `dim S_L = 1` at L = 3, 5, 7, 9, and every Cayley coefficient of `G_5`, `G_7`
  and `G_9` equals the closed form fitted on N = 9..48. There are 0 mismatches.
- **P10 PASSED.** Vir_(4,13) (N = 18, c = -217/26) has degree-18 lowest-weight generators exactly at Li 0, 3, 5, none at
  Li 6 or 7, and its Li-5 generator is `A_5 - (89/73) A_3 A_2 = G_5(18)`.
- Outputs: `check_witt_closed_forms_stable.out` and `.json`. The PBW-basis partial output is
  `check_witt_closed_forms.out`.
