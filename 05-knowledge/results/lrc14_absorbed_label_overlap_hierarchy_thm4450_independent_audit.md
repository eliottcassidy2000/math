# THM-4450 independent audit: overlap tiers and component addresses

**Status: ACCEPT / VERIFIED-EXACT; LRC(14) OPEN.**

Two clean-room referees independently checked complementary parts of
THM-4450 without importing the primary program.

The arithmetic referee rebuilt the Bernoulli expressions and a separate
rational-wall integrator. It reproduced the four class-uniform tenth-order
constants `1/63,1/70,1/70,1/77`, every subcritical and equality ratio, and
all four primitive sharp controls. It independently checks `184` equal-radius
and `158` mixed-radius literal pairs, three clock-two identities, three
clock-four fibre identities, every entry-threshold sum, and the mixed
`1/28` statistic. The mixed finite universe is exactly `pq<=49`: the
oscillation bound equals `1/28` at product `49` and is strict from `50`.

The decoder referee reconstructs the quotient walls directly. It reproduces
all five primitive open quotient-`y` component lists and their scale pullbacks
at `t=1,5,7`. It then builds every closed safe component, including isolated
points, for the four structured hostile bodies. Literal full-set containment
and the component-address criterion agree in all twenty ray checks. The
trapped-component and safe-endpoint vectors match the primary output, while
direct witnesses confirm that every hostile row is safe. This verifies that
the endpoint sieve is necessary only and that a single chosen component is
not a faithful state.

The coordinate distinction is explicit: primitive cells and beta widths are
quotient-`y` quantities. Their physical-`x` pullbacks are twice as numerous
and half as wide.

## Replay and hashes

```text
python -B 04-computation/lrc14_absorbed_label_overlap_hierarchy_thm4450_cleanroom.py
python -O -B 04-computation/lrc14_absorbed_label_overlap_hierarchy_thm4450_cleanroom.py
python -B 04-computation/lrc14_absorbed_label_overlap_hierarchy_thm4450_independent.py
python -O -B 04-computation/lrc14_absorbed_label_overlap_hierarchy_thm4450_independent.py
```

Normal, optimized, and frozen outputs are line-identical and end in `PASS`.

```text
04-computation/lrc14_absorbed_label_overlap_hierarchy_thm4450_cleanroom.py
  832d3012823d374046964affd3028275bd125e259e4f8183ed769255856079e3
05-knowledge/results/lrc14_absorbed_label_overlap_hierarchy_thm4450_cleanroom.out
  40bd7e00e62e0e2a87c591116e217d95a8e2087373202494f1bc6a4fc1743f98
04-computation/lrc14_absorbed_label_overlap_hierarchy_thm4450_independent.py
  174c5f2d0db6552595949283e9817c92aea9dd04a2c5517c3611e54f7839b940
05-knowledge/results/lrc14_absorbed_label_overlap_hierarchy_thm4450_independent.out
  ff175e08f3f0c2b4c230bb202c0fab4d53f6e5300bb5fb9afcfabd92c1df0b36
```

The independent computations audit finite atlases and controls. The
set-theoretic inclusion-exclusion, compact-versus-open equality gate, and
two-direction component decoder remain analytic parts of the proof.
