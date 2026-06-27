# LRC14 Apex-Aperture Comb Trichotomy

Created: 2026-06-24T16:02:03Z
Author: codex-2026-06-24-S154
Hypothesis: HYP-2967

The latest proof move: stop treating the many-14-multiple branch as live.
THM-571 already closes `|M14|>=7`.  The actual covering-strictness residue is
`1<=|M14|<=6`.

I added an exact local certifier:

```text
04-computation/lrc14_apex_aperture_comb_certifier_codex_s154.py
05-knowledge/results/lrc14_apex_aperture_comb_certifier_codex_s154.out
```

At a unit apex `a/14`, the non-14-multiple core is closed-safe.  If one side has
a strict core aperture and the few `14q` danger combs do not cover it, the
midpoint of an uncovered gap is a rational strict witness.

Full S151-bank audit:

```text
live_low_m14_rows = 18909
aperture-comb-certified = 12548
all-apertures-first-order-blocked = 4661
all-apertures-comb-saturated = 1700
```

Rebased over the incoming HYP-2964 Moon-core skeleton, HYP-2965 boundary-gap
packet bridge, and HYP-2966 NORK pinch-template atlas, this does not prove
LRC14, but it gives a sharper local residue:

```text
remaining bad row
  => full unit-support/AP-GW-skeleton core
     or comb-saturated tiny aperture.
```

The first branch should go through HYP-2960 skeleton/source labels or the
HYP-2908/THM-572 state lift.  The second branch should become a finite-comb
inequality: tiny aperture saturation forces a scale relation, hence either
scale separation or a bounded finite-core atlas.

Tournament note: the six-unit aperture tournament has only one achieved
isomorphism class, the transitive class.  The local scalar aperture quotient is
a proof gate, not a final classifier.
