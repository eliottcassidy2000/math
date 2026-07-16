# HYP-7064 — lift the 1,895 bad seven-shift directions

**Status:** OPEN after an exact finite residue reduction (codex-2026-07-16-S18).

`THM-908` assigns every nonzero projective residue direction `r in P^4(F_7)` the
signed affine-line ceiling

```text
C(r)=max_z sum_j L_6(z+jr).
```

Exactly `906/2801` directions have `C(r)<=32` and hence satisfy
`-F_6<=32/343<0.097` without using any speed magnitude.  The remaining `1,895`
directions have ceiling `34,36,48,50`, or `62`.  Every direction with five distinct
nonzero residues is already closed more sharply at `25/343`.

The proposed next lemma is a lift-sensitive improvement on the bad directions: the
within-cell trajectory `z(u)=sec(eu/7)` cannot remain on a maximizing affine coset for
enough measure unless the quotient speeds carry a short additive relation.  In the
nonresonant case, discrepancy should lower the average below `33`; in the resonant
case, `THM-906`'s B4 corner sum and `THM-907`'s B3 tail give a finite relation census.

This quotient preserves the signed residue-six kernel ceiling but destroys quotient
speeds, additive relation lattices, and wall chronology.  Those are mandatory sidecars,
not optional refinements, on the `1,895` residual directions.
