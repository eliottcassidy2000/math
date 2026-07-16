# HYP-7065 — projective pair-ray completion of the limiting residue-six sign

**Status:** RESOLVED by `THM-909` (codex-2026-07-16-S18).

`THM-908`'s raw seven-shift ceiling closes `906/2801` projective residue directions.
Restoring exact quotient pair data closes every remaining direction:

- a zero-coordinate uniform-marginal bound closes `1470/1505` zero directions;
- three unordered pair-ray certificates close the remaining `35`;
- the `675` bad nonzero directions reduce to twenty permutation/projective orbits;
- twenty affine-invariant labeled pair-ray certificates close those orbits, with worst
  bound `3240009/48020000=0.06747...`.

Together these give the universal limiting estimate

```text
-F_6 <= 32/343 < 0.097.
```

The solver-free verifier stores all literal rational weights and checks every pointwise
row and exact ray endpoint.  It also verifies the labeled ray interpolation on all
`3,540` ordered speed pairs through `60`; the universal interpolation law itself is the
Fourier theorem already proved before symmetrization in `THM-891`.

The quotient uses projective residue directions as vertices and determinant classes
`r_j z_i-r_i z_j` as pair sidecars.  It preserves the seven-shift signed kernel and the
labeled pair-ray law, but destroys wall chronology and higher relation lattices.  This
does not contradict the original pair-marginal obstruction: the decisive operation is
signed averaging before the pair projection.
