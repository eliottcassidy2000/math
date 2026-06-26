# LRC14 Fiber-Zipper Convergence

codex-2026-06-26-S188

## What Changed

HYP-3023 showed that exact magnitude kills route mixing in the automatic
fibers, but exact magnitude is close to packet identity.  This pass tested a
smaller convergence path on the first target word `MFCMMCCFFFCCC`.

The useful split is:

```text
residue-terminal fiber
  -> Erdos-Turan residue discrepancy bins
  -> Henselian unit rule
```

Erdos-Turan bins are not the final certificate; they are the contraction
gauge.  They reduce target mixed fibers from `27` to `6` and max mixed size
from `30` to `4`.

The Henselian unit rule is stronger than expected on this target: unit-root
and denominator-unit data at `p=2,3,7` gives `0` mixed route fibers before
exact magnitude is attached.

## Assumption Challenge

The vertices are not runners or raw automatic words.  They are zipper teeth:
ET discrepancy, Henselian unit roots, q/unit-excess lanes, exact magnitude,
barcode shadows, and packet labels.

The quotient preserves the LRC predicate by retaining enough data to keep
theorem routes pure.  It destroys raw row identity and much of exact magnitude.

## Next Pull

Turn the target-fiber readout into a lemma: simple unit roots lift, singular
unit roots carry local debt, and the remaining denominator-unit classes route
to q-witness, AP/GW boundary, petal, or covering.  Then run the same zipper on
the full HYP-2963 bank.
