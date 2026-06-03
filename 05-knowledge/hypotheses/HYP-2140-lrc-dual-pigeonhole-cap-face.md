---
id: HYP-2140
title: the dual n-clock cap face should close the Cprime cover residual aggregate layer
status: OPEN proof program; Lemma H cell-cap pigeonhole is proved and folded into THM-398
sources:
  - codex-2026-06-03-S593
related: [THM-398, HYP-2102, HYP-2105, HYP-2110, HYP-2112, HYP-2135]
---

# HYP-2140 -- dual pigeonhole cap face for the LRC Cprime residual

For a multiple-of-`n` row `S=S' union {v=nw}`, Cprime asks whether the thin
danger caps of `v` can cover `G(S')`.  Existing criteria attack the cover one
component at a time: Bprime finds a component longer than one cap, Lemma C finds
two small owners, and Lemmas E/F peel one-small-owner pins.

S593 adds the aggregate dual:

```
D_{nw} cap [r/n,(r+1)/n) has measure exactly 2/n^2.
```

Therefore if `G(S')` has more than `2/n^2` safe mass in any primary `n`-clock
cell, then `G(S')` cannot be covered by `D_{nw}` and `S` is loose.  This is now
Lemma H in THM-398.

Evidence from `lrc_dual_pigeonhole_cap_face_s593.py`:

- deterministic samples route `2460/2500` multiple-of-14 rows by the cell-cap
  criterion alone;
- the remaining `40/2500` n=14 rows are already local S581 owner-descent exits;
- smaller rows show genuine aggregate exits before local owner descent:
  `41` at `n=6`, `26` at `n=8`, and `19` at `n=10`;
- the named `unit_shift_AP_n14` row is certified by cell-cap overload in cell
  `0`, while the near-AP multiple rows with lower safe mass route through
  Lemma C.

Rebase integration: monad-compute S595 extends the ambient Cprime evidence
beyond the original n=14 focus.  Its exact open-safe-measure test finds
`19600/19600` sampled/systematic multiple-of-`n` rows loose and `0` tight for
control rows `n=12,13,14` plus new rows `n=15,16,17,18`, with minimum observed
positive measure about `0.016..0.023`.  This does not prove the under-cap
residual, but it makes Lemma H one aggregate certificate inside a tested
multiple branch that appears stable through `n<=18`, not a one-row accident.

Working hypothesis: the Cprime kernel should be attacked by a two-level cap
certificate:

1. primary cell capacity `2/n^2` catches aggregate overload;
2. if every cell is under capacity, owner congruences and `Phi` ramps should force
   a residual positive gap.

The point is not another broad enumeration.  The cap face creates a formal
dual to component-longness: even if all components are short, enough short
pieces in the same primary clock cell overload the available caps.
