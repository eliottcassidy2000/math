# HYP-7064 — lift the 1,895 bad seven-shift directions

**Status:** RESOLVED by `THM-910` (codex-2026-07-16-S18).

`THM-908` assigns every nonzero projective residue direction `r in P^4(F_7)` the
signed affine-line ceiling

```text
C(r)=max_z sum_j L_6(z+jr).
```

Exactly `906/2801` directions have `C(r)<=32` and hence satisfy
`-F_6<=32/343<0.097` without using any speed magnitude.  The remaining `1,895`
directions have ceiling `34,36,48,50`, or `62`.  Every direction with five distinct
nonzero residues is already closed more sharply at `25/343`.

The successful lift is simpler than the proposed discrepancy/B3/B4 split.  A zero
residue supplies an exactly uniform quotient coordinate; this plus three unordered
pair-ray rows closes all zero-containing directions.  With no zero residue, the affine
determinants `r_j z_i-r_i z_j` have expectations governed by the labeled pair-ray law;
twenty rational orbit certificates close the remaining `675` directions.

This quotient preserves the signed residue-six kernel ceiling but destroys quotient
speeds, additive relation lattices, and wall chronology.  Those are mandatory sidecars,
not optional refinements, on the `1,895` residual directions.

The exact verifier also compares the signed ceiling to the positive box-certificate
ceiling.  Their fifteen-stratum representative tournaments are transitive but differ by
nine edge flips, quantifying the cancellation lost by a positive majorant.
