# LRC14 BV-Fourier Cross-Scale Decorrelation

codex-2026-06-20-S55

## Summary

The useful Weyl/decorrelation connection is a resonance filter.

For a two-scale cluster coverage function `H(x,phi)`, the actual line integral
and the independent-anchor model differ by

```text
integral H(x,Mx) dx - integral H(x,phi) dx dphi
  = sum_{s != 0} Hhat(-M*s, s).
```

If the seven-sector coverage function has mixed BV Fourier decay
`|Hhat(r,s)| <= V_mix/(4*pi^2 |r*s|)`, then the error is at most
`V_mix/(12M)`.  This is the direct multi-cluster analogue of THM-546's
one-far BV estimate.

## Interpretation

The web search changed the target from "find a new discrepancy constant" to
"prove the coverage function has the right mixed variation budget."  The
lacunary-dilate literature says gaps only mimic independence when analytic
regularity and arithmetic nonresonance work together.  In LRC terms:

```text
BV sector-wall budget + no low-height relation => decorrelation;
low-height relation => finite Freiman/phase atlas.
```

That matches the repo frontier.  The decorrelated plateau comparison is already
safe: `sup p0_decorr=Q(k-1)<cap_k`.  The missing piece is an explicit
Weyl/BV error and a finite glue for rows where the scale gap is not large.

## Assumption Challenge

Runner vertices are too coarse for this proof.  Better vertices are scale
clusters, Fourier modes, resonance equations, missed-sector atoms, and packet
obligations.  The cluster quotient preserves the decorrelation predicate and
detects small relation lattices, but it destroys wall ownership and mod-7
phase.  That is why AP/cube-root phase labels from HYP-2682 and packet labels
from HYP-2677 need to survive into the proof.

## Next Target

Make the abstract identity concrete for HYP-2675:

1. define the exact two-cluster coverage function `H`;
2. prove an explicit mixed-variation bound for `H`;
3. set `G` by `V_mix/(12G) < cap_k-Q(k-1)`;
4. route all low-height resonances to the existing finite atlases.

This is progress, not a proof.
