# HYP-1805: Kernel pinning principle

**Status:** EXPLORATORY meta-hypothesis.

## Statement

Many repo obstructions should be attacked by first identifying a kernel whose
second moment is already a pinned local obstruction:

```text
compressed shadow -> kernel/autocorrelation -> pinned local statistic
```

Then the remaining proof work is an exclusion/incidence problem, not a
shadow-to-locality problem.

## Evidence

HYP-1804 gives the model example for circulant tournaments:

```text
Var_{v != 0} J_3(0,v)
  = E(S)/(p-1) - (p^2 - 2p + 5)/16.
```

Thus Fejer/IPR concentration, additive energy, and pinned triangle localization
are the same one-dimensional axis.  The hard bridge begins at `J_5`, where
repeated-vertex exclusion diagrams appear.

Other nearby examples have the same shape:

- Lonely Runner: forbidden intervals form a finite endpoint-cover kernel; the
  obstruction is endpoint protection incidence.
- Caccetta-Haggkvist: out-neighborhood growth is a return kernel; the
  obstruction is the first return residue.
- Endpoint transfer: support shadows are insufficient; private pivots and
  collision hypergraphs are the incidence lift.

## Prediction

For several open repo questions, the next simplification will come from
rewriting a global statistic as the variance or second moment of a pinned
kernel, then classifying the illegal coincidence diagrams.

Immediate candidates:

- `J_5` pinned cycle profiles for circulant tournaments.
- Path-homology boundary-rank variance.
- Real-root Newton defects for `I(Omega,x)`.
- Endpoint-protection graphs for Lonely Runner tight instances.

## See

- `07-reflections/kernel-residue-trick-atlas-2026-05-30.md`
- `07-reflections/fejer-kernel-wild-session-2026-05-30.md`
- `04-computation/fejer_kernel_wild_session.py`

