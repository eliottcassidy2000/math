---
id: HYP-2140
title: the dual n-clock cap face should close the Cprime cover residual aggregate layer
status: OPEN proof program; Lemma H cell-cap pigeonhole is proved and folded into THM-398
sources:
  - codex-2026-06-03-S593
  - codex-2026-06-03-S598
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

Rebase integration: monad-compute S595/S596 extend the ambient Cprime evidence
beyond the original n=14 focus.  Their exact open-safe-measure tests find no
tight sampled/systematic multiple-of-`n` rows through `n<=20`, with positive
measure margins staying around `0.016..0.023` in the new frontier rows.  This
does not prove the under-cap residual, but it makes Lemma H/I aggregate
certificates inside a tested multiple branch that appears stable beyond
fourteen, not one-row accidents.

Working hypothesis: the Cprime kernel should be attacked by a layered cap
certificate:

1. primary cell capacity `2/n^2` catches aggregate overload;
2. origin-bisection capacity `1/n^2` catches one-sided overload once a small
   endpoint owner pins components to cap centers;
3. if every cell and one-sided cell is under capacity, owner congruences and `Phi` ramps should force
   a residual positive gap.

The point is not another broad enumeration.  The cap face creates a formal
dual to component-longness: even if all components are short, enough short
pieces in the same primary clock cell overload the available caps.

S598 refinement: the cap centers themselves are proof data.  A pinned small
left endpoint can only use the upper half of its cap, and a pinned small right
endpoint can only use the lower half.  Thus each primary cell has one-sided
capacity `1/n^2`.  The named `apex_AP_replace_13_n14` row misses the total-cell
surplus test but has upper surplus `1/1176`; this is the first clean example of
an origin-bisection law becoming an upper cap certificate.

Close-out rebase integration: Opus S595/HYP-2142 names the next residual after
these cap filters.  If total caps and one-sided origin-bisection caps are all
under capacity, the remaining large-owner cover constraints become a 2x2
determinant / bounded CRT feasibility problem.

S599 extension: HYP-2144 asks for small Helly witnesses inside that determinant
residual.  Instead of building the whole CRT modulus first, enumerate the
dominance-bounded `w` window and extract a singleton component wall or pair of
incompatible owner/slack determinant rows.  In the S599 deterministic sample,
n=14 rows that reach the post-B/C determinant layer already collapse to
singleton walls, while smaller n show the first pair witnesses.
