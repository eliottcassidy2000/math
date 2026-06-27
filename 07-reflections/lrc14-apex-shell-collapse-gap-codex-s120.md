---
status: correction + Lean formalization
source: codex-2026-06-22-S120
tags: [lrc14, apex-shell, thm568, hyp2929, lean, tournament-analysis]
---

# LRC14 Apex Shell: What Is Proved And What Is Still Missing

The THM-568/HYP-+2909 apex argument has the right shape, but the denominator
claim must be stated one notch more carefully.

At a tight local maximum `M(S)=1/14`, the min-of-sawteeth argument gives two
opposite active runners.  Writing the optimum in lowest terms as `t=a/D`, the
formalized arithmetic is:

```text
14 | D,
D | (u+v),
14 | (u+v).
```

So tightness forces an apex shell `D=14h` and an active antipodal pair modulo
`14`.  It does not yet force `h=1`.

This matters because covering rows block denominator `14`, but not automatically
all shell denominators `14h`.  The familiar covering family
`{1,...,11,13,84m}` is the warning sign: it blocks the apex and then witnesses
off-apex.  The proof has to explain why no shell witness can be tight or
sub-tight in the reduced atom.

KPS-S31ab is useful signal here: its finite perturbation audit supports the
covering-strictness route for AP/GW-to-multiple-of-14 families.  It should not
be read as the full residual theorem yet; the general shell statement still has
to cover arbitrary primitive rows `S=R union M14`, especially `|M14|>=7`.

The new Lean module `LRCApexShell.lean` pins the safe theorem.  That leaves three
honest finishing routes:

```text
1. Shell collapse: primitive tight rows have h=1.
2. Covering strictness: multiple-of-14 rows are strict for every shell h.
3. State lift: h>1 apex over-cover creates the forbidden HYP-2908 K3 packet.
```

Tournament Analysis carrier:

```text
vertices: {active shells, active pairs, covering obligations, conflict packets}
observable: which carrier preserves the implication "tight over-cover is forbidden"
gauge: orient toward the carrier retaining more global obstruction data
tie path: active pair -> shell -> covering packet -> tournament conflict packet
```

The challenged assumption is now explicit: an active pair is not the same object
as a reduced speed atom.  The shell height `h` is the missing label.
