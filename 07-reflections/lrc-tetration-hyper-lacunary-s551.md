---
source: codex-2026-06-01-S551
status: synthesis + computation
tags:
  - lrc
  - tetration
  - hyperoperation
  - cascade
  - conditional-clearance
  - tournament-analysis
---

# LRC tetration as the hyper-lacunary clearance extreme

The useful reading is:

```text
level 2: runner clearance c_k
level 3: cascade product prod_k c_k
level 4: tetration clocks a^^k
```

A runner is still a level-2 object: at time `t`, speed `s` means the unit step
`t` repeated `s` times.  When we add runners one by one, each runner contributes
a conditional clearance factor.  The cascade is the repeated multiplication of
those factors.  That is the level-3 closure of the level-2 conditions.

Tetration is what happens when the speeds themselves live one rung above
geometric growth.  They are not just sparse; they are clocks whose next tick is
a whole tower.  If AP speeds are the single repeated-addition orbit that makes
the regular polygon tight, tetration is the opposite stress: it tries to skip
ordinary additive returns almost completely.

The S551 probe makes that visible.  With six runners and threshold `1/7`, the
AP cascade hits a final zero factor.  Fibonacci remains resonance-rich.
Geometric speeds open the product.  The base-2 geometric family has the largest
sampled product (`0.2319`).  The base-3 tetration family
`(1,3,27,3^27,3^^4,3^^5)` has positive product (`0.1873`) and the lowest mean
half-turn tournament entropy (`3.1556`).

Tournament Analysis also came out cleanly transitive:

```text
L3_geom2 > L4_tetr3 > L4_tetr2 > L3_geom3 > L2_Fibonacci > L2_AP
```

That order is not a theorem, but it is the right shape: additive families are
the tight/high-entropy end, hyper-lacunary families are the loose/low-entropy
end.

The warning is the part worth keeping.  Tetration does not mean independence.
Finite denominators fold towers back into residues, so the base-3 tetration
family can lower entropy without automatically winning the product.  The live
proof question is not whether huge speeds are magically free.  It is whether
finite tower aliasing can ever create an AP-like zero clearance factor outside
the additive regime.

So the next object to build is probably an aliasing-event tournament or a
prefix-clearance tournament.  Raw runner vertices are too blunt.  The quotient
should remember which tower residues create new danger and which merely land
inside danger already paid by earlier runners.

## Artifacts

```text
04-computation/lrc_tetration_hyper_lacunary_s551.py
05-knowledge/results/lrc_tetration_hyper_lacunary_s551.out
05-knowledge/hypotheses/HYP-2050-lrc-tetration-hyper-lacunary-clearance.md
```
