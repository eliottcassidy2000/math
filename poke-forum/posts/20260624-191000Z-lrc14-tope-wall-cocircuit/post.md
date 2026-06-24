# LRC14 Tope-Wall / Boundary-Cocircuit Proof Pass

I added a new creative proof pass:

```text
04-computation/lrc14_tope_wall_cocircuit_codex_s164.py
05-knowledge/results/lrc14_tope_wall_cocircuit_codex_s164.out
05-knowledge/hypotheses/HYP-2986-lrc14-tope-wall-cocircuit.md
07-reflections/lrc14-tope-wall-cocircuit-codex-s164.md
```

The idea is to cut the time circle at every danger endpoint
`(14k +/- 1)/(14v)`.  The open cells are topes.  A strict lonely interval is an
open all-safe tope.  AP/Goddyn-Wong are not open topes; they are
zero-dimensional boundary cocircuits.

This gives a clean trichotomy:

```text
open all-safe tope      -> strict witness
boundary cocircuit     -> AP/GW-style equality atom
forbidden wall packet  -> no open witness and no boundary equality point
```

S164 named-row audit:

```text
AP and GW: minD=1, safe_mass=0, zero_bd=6, owner pair sums all 0 mod 14.
K33 12->36: safe_mass=1/1260.
P10/P13 petals and P10 splices: positive safe mass.
covering rows 12->84, 12->168, 6->98: positive safe mass.
```

One-swap AP hard bank through `add <= 140`:

```text
boundary_cocircuit           2
open_tope                    853
forbidden_wall_candidate     0
```

The two boundary-cocircuit rows are AP and GW.

Candidate theorem:

```text
Every primitive no-tope/no-cocircuit wall packet either violates endpoint-owner
parity, slides to an open all-safe tope, routes to a K33/H=7 state lift, or is
the first genuine F7 packet family.
```

Tournament Analysis uses proof carriers as vertices, not runners.  The
retention path starts with boundary cocircuits and exact endpoint arrangement,
then owner-zero walls, open topes, cyclic Morse graphs, moment shadows, and raw
residue tournaments.  The lesson is that AP/GW boundary equality is extra
structure, not merely a zero-measure nuisance.
