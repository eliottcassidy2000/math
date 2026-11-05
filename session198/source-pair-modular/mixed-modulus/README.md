# Reviewed-source C2/C3 all-modulus transfer

The proof is in `PROOF.md`. It constructs actual source pairs for every
parity-compatible residue modulo every positive integer, separately with two
or three Hamilton cycles. It does **not** establish exact odd-target coverage.

Run from a copied packet directory:

```
python3 verify.py > VALIDATION-normal.json
python3 -O verify.py > VALIDATION-optimized.json
cmp VALIDATION-normal.json VALIDATION-optimized.json
```

Both modes write `VALIDATION.json`. The verifier is standalone Python standard
library code and constructs literal graphs instead of importing source-oracle
code. Dependencies are frozen proof snapshots; their SHA256 values appear in
the validation receipt. The all-order proof, not the finite checks, establishes
the congruence theorem.
