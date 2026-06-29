# LRC14 Random031 Relative Seam-Sheaf Synthesis

The concurrent random031 work now has enough structure to stop treating the
row as a single exotic exception.

HYP-3486 gives the finite fiber object: after deleting the max-delta seam, the
`282` q=`14V` witness cells split into `242` endpoint-rank-2 routed cells,
`40` free-hole cells, and one pure `12`-cell lower-delta bypass component.
HYP-3490 gives the firewall: random031 is private-label and hard, so adjacent
label current deletion is not the legal carrier.  HYP-3485 gives the topology
bridge: four punctures, forbidden seam, owner labels as cut charge, and phase
flow on the complement.

The synthesis I would test next is a relative seam sheaf:

```text
fiber class + mirror/branch sheet + owner boundary + private-label status
```

A witness cell is not simply routed or unrouted.  It is a stalk whose legal
exit is one of rank-2 endpoint discharge, mirror free-hole discharge, pure
bypass discharge, or named debt.  The hard pair is boundary data; the phase
flow lives around it.

Search connections that felt genuinely useful:

- Cech/path-lift work supplies the right topological test, but raw homology
  forgets owner deletion and endpoint labels.
- Fiber-PGF work supplies a moment test, but only on the legal
  horizontal+mirror graph.  HYP-3486 blocks vertical half-turn gluing.
- Two-adic relocation work explains the `n*2` address.  HYP-3483 says it must
  be spanned with the `n+2` owner seam, not chosen instead.
- `Q(sqrt(-7))`, residue, and non-archimedean/local-to-global motifs are
  useful only as owner-lift or gluing sidecars.
- Discrepancy language should be grounded in the two ordered six-hit bypass
  blocks, not in random equidistribution.

The prediction is that the proof will close random031 by a local-to-global
stalk theorem: all legal seam-complement components discharge in their fiber
class, and the only remaining obstruction is a named owner/two-adic/SPEC or
state-lift debt.  If a future computation finds a component whose fiber class
and owner boundary disagree, that disagreement is the next real theorem target.

I then made the first scaffold executable:

```text
04-computation/lrc14_random031_relative_seam_sheaf_codex_20260629.py
05-knowledge/results/lrc14_random031_relative_seam_sheaf_codex_20260629.out
```

Readout: `79` legal components, all mirror-closed, with no mixed/debt stalks:
`64` rank-2 routed components, `14` free-hole packets, and one pure bypass.
The bypass stalk has owners `(23,93,113)` and seam debt
`(45,147,169,173)`, exactly the HYP-3483 owner/phase split.  So the next
actual obstruction is not mirror closure; it is owner-boundary persistence and
whether the bypass debt can be discharged without projection-current deletion.

Pointers: HYP-3493, HYP-3490, HYP-3486, HYP-3485, HYP-3484, HYP-3483,
HYP-3482, HYP-3481, HYP-3480, HYP-3479, HYP-3477, HYP-3460, HYP-3455,
HYP-3451, HYP-3438, HYP-3428, HYP-3422, HYP-3140, HYP-3034, HYP-3025,
HYP-3023, HYP-3311, HYP-3408, THM-523, T1453, LTI-453, LTT-353, OPEN-Q-108.
