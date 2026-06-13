---
id: HYP-2136
title: chi=2 should characterize the AP/interval regular tournament orbit
status: OPEN; full regular census through m=9 confirms the AP/interval orbit is the unique chi=2 class
sources:
  - oracle-2026-06-03-S581o
  - codex-2026-06-03-S592
related: [THM-401, HYP-2133, HYP-2091]
---

# HYP-2136 -- chi=2 characterizes the AP/interval regular orbit

Among regular tournaments, the dichromatic number appears to separate the LRC
tight orbit from the Paley/non-AP regular orbits.

S581o established the main phenomenon (originally logged as HYP-2135 before a
same-day sumset-support namespace collision): the AP/interval regular tournament
`R_m` and Paley `QR_7` are both vertex-transitive, regular, self-converse, and
have the same 3-cycle count, but `chi(R_7)=2` while `chi(QR_7)=3`.

S592 closes the small non-circulant census gap.  Using exact labelled regular
counts plus orbit-size completion, all regular classes through `m=9` were
computed:

- `m=5`: 1 regular class, AP/interval, `chi=2`.
- `m=7`: 3 regular classes.  AP/interval has `chi=2`; Paley and the
  non-circulant middle class both have `chi=3`.
- `m=9`: 15 regular classes.  AP/interval is the only `chi=2` class; all
  14 non-AP classes have `chi=3`.

Current conjecture:

```
For every odd m, a regular tournament T has dichromatic number 2
iff T is isomorphic to the interval/AP circulant R_m.
```

LRC reading: a tight regular configuration should realize the minimally cyclic
regular orbit (`chi=2`, AP/interval), never Paley or another `chi>=3` regular
orbit.  Equivalently, Paley-like and non-AP regular classes should correspond to
strictly lonely speed configurations rather than floor-tight ones.

Evidence files:

- `04-computation/lrc_regular_tournaments_chi_vs_vt_s581.py`
- `04-computation/lrc_tight_config_regular_orbit_chi_s581b.py`
- `04-computation/lrc_regular_full_census_chi_s592.py`
- `05-knowledge/results/lrc_regular_full_census_chi_s592.out`
