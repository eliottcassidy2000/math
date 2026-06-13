# HYP-1780: Residue rank stratifies the next obstruction

**Status:** OPEN meta-hypothesis; one new Lean residue formalization.
**Session:** kind-pasteur-2026-05-30-S3

## Statement

Many project obstructions are not controlled by the raw size of a support
family, but by the smallest surviving residue after a natural projection.
The next obstruction in a family appears when that residue has just enough
rank, parity, or independence to survive cancellation.

In slogan form:

```text
large support family + projection = residue;
first obstruction = smallest nontrivial residue with the wrong shape.
```

## Evidence

- OCF reduces Hamiltonian-path parity to the mod-2 residue `H(T) % 2 = 1`.
  The new Lean theorem `H_mod_two_eq_one_from_ocf` records this residue
  explicitly.
- THM-354 says good-cut height is not a path-coordinate statistic; its residue
  is the strong-component defect `n - #SCC(T)`.
- THM-351/352/353 say Boolean-cube transport rows have forced conservation;
  the interesting information is the off-diagonal escape residue.
- H=63 single-core classes are exact projection kills: deleting the core
  vertex leaves no odd-cycle residue.
- THM-025 is a near-kill: almost all odd cycles die after deleting one vertex,
  but the small residue still has a top-heavy independence profile.
- Paley/Interval comparisons show that equal support shadows can have different
  multiplicity and disjointness residues.

## Predictions

1. First non-real-rooted `I(Omega,x)` examples should cluster at small but
   nonempty deletion residues, not merely high total odd-cycle count.
2. H-spectrum gaps beyond parity should be explainable as impossible residue
   profiles in the OCF conflict graph, not just as missing numeric values.
3. Bucket-transport features should be organized by residual rank of the
   off-diagonal transport matrix after the row checksum is removed.
4. Path-homology ghost-cycle failures should resemble near-kill odd-cycle
   residues: most structure dies under the projection, but one low-rank
   survivor remains outside the expected boundary image.

## Next Tests

- Use the new S355 deletion-residue rank and top-heavy alpha-profile features
  in the `omega_*` Tournament TDA block to scan n=9 real-root failures.
- For good-cut transport, replace raw `Delta g` with `Delta #SCC` and compare
  against spine/ribs/sea escape mass.
- For H-spectrum work, collect candidate alpha-vectors by residue rank and
  see whether forbidden values correspond to absent low-rank profiles.
- For path homology, compare ghost-cycle examples against the THM-025
  near-kill signature.

## Related

- HYP-1779: support-residue calculus.
- HYP-1785: deletion-residue rank filter.
- THM-025, THM-344, THM-351, THM-352, THM-354.
- `07-reflections/residue-feedback-loop.md`.
- `04-computation/residue_rank_probe_s355.py`.
