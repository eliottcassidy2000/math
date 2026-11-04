# S188: close the fixed-quotient module return route

The new result is an all-order contextual no-return theorem, not a global spectrum obstruction.

Replace address24's mixed-port homogeneous pair by **any tournament of order at least2**, retain the fixed external quotient and its external ports, and choose any nonempty internal port set. This cannot attain either next constant -103 or -101, even after arbitrary marked ordered contexts and arbitrary fixed integer clock shifts.

The proof uses exact decorated path-cover coefficients and two elementary splitting injections. Every module of order at least4 already has c<=-106. The remaining two- and three-vertex modules reduce to the explicitly checked old small candidates. This is not a bounded-search extrapolation.

A positive side result gives an explicit growing transitive-module strong family preserving old names and arcs, with all seven assembly coordinates. Its joint defect is exactly35*2^(m-1), the sharp minimum over nonempty transitive module-port sets. The proof explains why the nonzero3^m two-cover channel does not survive in the actual attachment count.

- `PROOF.md`: full proof, exact scope, sources and remaining obligation. Section7 is the strongest general-module result.
- `carrier.py`, `inherited.py`: literal source-derived carrier and attributed subset-DP/attachment helpers.
- `skeleton.py`, `family_coeffs.py`: symbolic transitive joint/full-signature certificates.
- `verify.py`, `receipt.json`, `normal.stdout`, `optimized.stdout`: transitive, direct-attachment and independent-permutation checks.
- `module_coeffs.py`: nonnegative decorated skeleton coefficients for arbitrary modules.
- `module_extension.py`, `module_extension.receipt.json`, `module_extension.normal.stdout`, `module_extension.optimized.stdout`: all75 labelled modules through order4,1,023 nonempty marked modules,285 exact singleton formulas and570 injection checks.

Reproduce with `python verify.py`, `python -O verify.py`, `python module_extension.py` and `python -O module_extension.py`; standard library only. Both normal/optimized pairs byte-match. No native Hamilton-complement program, reusable two-branch closure, global coverage or historical novelty claim is made.
