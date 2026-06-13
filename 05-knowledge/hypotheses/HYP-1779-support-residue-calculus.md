# HYP-1779: Support-residue calculus organizes projection phenomena

**Status:** OPEN meta-hypothesis; THM-354 gives one new confirmed instance.
**Session:** codex-2026-05-30-support-residue

## Statement

Many tournament invariants in this repository are best understood as residues
of support families under a projection or quotient:

```text
support family -> projection/forgetful map -> surviving residue geometry.
```

The projection itself is usually lossy. The useful invariant is the defect:
which certificates are killed, which survive, and what independence, parity, or
homology structure remains among the survivors.

## Confirmed Instance

THM-354 proves that the apparently base-path-dependent good-cut count is really
the strong-component defect:

```text
g(T, P) = n - #SCC(T)
```

for every Hamiltonian path `P`. Thus the good-cut projection descends to
tournament isomorphism and complement-merged classes because SCC count does.

## Evidence Across Threads

- H=63 single-core classes: deleting the core vertex kills all odd cycles, so
  `Omega=K31` and `H=1+2*31`.
- THM-025 real-root failure: deleting the near-core vertex leaves a tiny
  two-cycle residue with enough disjointness to support the unique independent
  triple.
- Paley versus Interval: support shadows and multiplicities separate; Interval
  wins through disjointness residues amplified by THM-143's co-occurrence
  gradient.
- Bucket balance: internal half-lines are paired by a fixed-point-free
  involution; cross-lines are the quotient escape residue.
- Transfer diagonal: `M[a,a]` is the even-minus-odd position residue after
  projecting Hamiltonian paths to the position of `a`.
- Path homology ghost cycles: the old projection asks whether a through-vertex
  cycle survives as a genuine old class or dies into a boundary.

## Testable Predictions

1. Non-real-rooted `I(Omega,x)` examples at n>=9 should be enriched by large
   maximum deletion-loss fraction plus a small nonempty residue with top-heavy
   independence profile.  S355 turns this into explicit features:
   `omega_near_kill_vertices`, `omega_near_kill_rank2_vertices`, and
   `omega_max_loss_residue_rank`.
2. Any future `Delta g != 0` transport theorem should reduce to a change in
   SCC count, not to a special property of the base-path tiling coordinates.
3. In circulant Paley/Interval comparisons, THM-143 slope profiles should
   predict the first Walsh degree where Interval's disjointness advantage
   overtakes Paley's multiplicity advantage.
4. HYP-408-style ghost-cycle failures should resemble THM-025 near-kills:
   almost all relevant support dies under projection, but a small structured
   residue remains outside the expected boundary image.

## Related

- THM-354: good-cut count equals `n - #SCC(T)`.
- THM-344: n=8 H=63 complete-Omega classes.
- THM-025: real-rootedness failure at n=9.
- THM-143: interval co-occurrence linearity.
- THM-350: finite unordered bucket balance layer.
- HYP-1780: residue rank stratification.
- HYP-1785: deletion-residue rank filter.
- `07-reflections/support-residue-calculus.md`.
- `07-reflections/residue-feedback-loop.md`.
