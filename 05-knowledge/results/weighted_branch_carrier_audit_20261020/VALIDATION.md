# Validation receipt

- Exhaustive independent switching-state enumeration at n=2,4,6,8: PASS.
- Exact rational conjugacy-class centralizer sum: PASS class by class.
- n=8 weight counts 19,14,16 at 1/2,1/4,1/8: PASS.
- Correct marked fixed count 7/31, unmarked coincidence count 9/31: PASS.
- Independent matrix implementation, 12,114 commuting-isomorphism tests: PASS.
- D8 hostile: complete S8 stabilizer check and explicit exchanged markings: PASS.
- C++ optimized/replay and address/undefined sanitizer comparison: identical output.
- Python ordinary/-O marked checks and comparison controls: identical output.
- Partition branch totals n=1..16: reproduced, including exceptional n=2.
- Optional terminal-involution flip on all 2-elements: FINITE-EXACT 9/49 fixed at n=8; no denominator-clearing theorem.
- Root mathematical and scope review: PASS after the explicit denominator-clearing caveat.

The full report specifies universes, encodings, controls, losses, source provenance and reproduction commands. An involution alone does not clear denominators 4 or 8.
