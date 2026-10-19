# Torus-filling geometry extensions, 2026-10-19

## Results

Three connections produced precise results beyond the existing source audit.

1. **Seifert gauge torsion.** For arbitrary orders `a,b`, the specified
   two-move gauge quotient is `Z + Z/gcd(a,b)`. An exact primewise formula
   counts admissible labelled classes at each signed clutch value `p`.
   Three distinct admissible classes at `(a,b,p)=(5,5,10)` have the same
   cyclic presentation cokernel. Coprime orders are exactly the boundary
   where a `|p|=1` class can occur.
2. **Exterior monodromy versus cusp index.** For a rank-two square-zero
   shear with Smith factors `s|t`, exterior degrees one and three retain
   `(s,t)`, but degree two retains only `(s,s)`. The mapping-torus middle
   homology therefore also loses `t`. Under the earlier translation-fan
   hypotheses the toric cusp has `st` components and Euler `2st`.
   The family `diag(1,d)` keeps degree-two Smith data constant while these
   geometric counts grow. Rank alone also rules out injective top-degree
   specialization when `d>1`.
3. **Conductor deformation holonomy.** The full normal-crossing cotangent
   presentation gives the line bundle `T^1=Ext^1(Ω,O)` on the conductor.
   Degree `-1-1+2=0` on its three normalized curves leaves two multiplicative
   gluing holonomies. Both must equal one for a smooth-total-space
   smoothing. For the manuscript's displayed toric maps they do equal one:
   every map fixes the global character `t`, hence the dual conormal frame.
   This is an equivariant certificate conditional on the claimed local
   unramified quotient; it does not verify that quotient's free/proper
   analytic existence.

The PDF is byte-identical to the repository's August 24 artifact, so this
is an extension session rather than a newly revised-source intake.
SHA-256: `283bba102dd1d5dc346af81b28145bdaaea6654398d5032e76e97bafb9a858f2`.
The complex-S6 headline remains a manuscript claim under audit.

## Inheritance pass

Closest proved mechanisms: the coprime Seifert quotient in the earlier
Hopf/S6 intake; THM-3955 and THM-3957's conductor cotangent presentations;
THM-3991's lattice-index/factorial Euler theorem.

Canonical hostiles: reversed branch orientations change integral homology
without changing Euler counts; raw gauge coordinates need not be genuine
moduli; normalized differentials erase conductor-supported information.

Corrected near miss: MISTAKE-480, which forbids declaring the raw `p` kernel
an analytic moduli sidecar before computing the actual gauge quotient.

Least-used sidecars recovered here: the second integral shear Smith factor,
noncoprime gauge torsion, and the actual line-bundle transition products on
the conductor rather than only their normalized restrictions.

## Concept board and portfolio

- Anchor: global torus-filling presentation, represented by integral twist
  lattices. Test: which marked gauge classes share the cyclic clutch group?
- Niche: integral nilpotent monodromy, represented in all exterior degrees.
  Test: which degree retains the finite index of the shear image?
- Wildcard: nonnormal conductor deformation theory, represented by Ext^1
  and graph holonomy. Test: degree-zero versus trivial line bundle.
- Boundary: periodic A2 fan quotient, represented by lattice index and
  component strata. Test: Euler and top-degree specialization dimensions.
- Literature bridge: Mumford degeneration and small toric resolutions,
  represented by square-to-triangle subdivisions. A global product-Tate
  realization of the hostile family is proposed, not claimed proved.

The first discovery, noncoprime torsion, sharpened the quotient ledger:
finite presentation torsion is real, but still not analytic moduli. The
second discovery showed that even saturated degree-two monodromy can miss
integral geometric size, forcing a degree-one sidecar. The third discovery
turned a normalization loss into a positive local gluing certificate:
actual transition characters, rather than only degrees, settle its two
holonomies from the displayed toric maps.

## Artifacts

- `seifert_gauge_torsion_count.md`: unrestricted lattice/count theorem and
  proof, examples, exact loss ledger, and scope.
- `monodromy_smith_cusp_index.md`: exterior Smith theorem, Wang-sequence
  corollary, conditional cusp count and specialization obstruction.
- `conductor_smoothing_holonomy.md`: local Ext calculation, triple-point
  formula, seminormal graph Picard calculation, and actual toric certificate.
- `verify_geometry_extensions.py` and `.out`: exact finite independent
  controls. There are 277,729 Seifert residue cases, 961 gauge Smith
  matrices, 500 exterior matrices, four mixed-basis controls, five endpoint
  complexes, 50 translated A2 charts, and a symbolic shear-character check.

All files are outside the shared checkout for the designated publisher.
No Git operations or checkout writes were performed by this worker.
