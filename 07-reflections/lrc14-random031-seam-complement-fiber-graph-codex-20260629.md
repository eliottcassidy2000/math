# LRC14 Random031 Seam-Complement Fiber Graph

HYP-3486 is the fiber-graph addendum to the HYP-3484 forbidden-seam surgery
packet and the HYP-3485 seam-complement connection atlas.  It turns the
HYP-3482 experiment design into an exact finite graph.
The right base coordinate is `u=2t mod 1`, with branch remembered as the
sheet.  That makes `random031` feel less like a wall problem and more like a
two-sheeted cylinder whose legal symmetry is mirror, not vertical half-turn.

The useful correction is that "route all 282 witnesses through gates" was too
blunt.  After deleting the max-delta seam, `242` cells hit exactly one
endpoint-rank-2 seam-complement gate, while `40` cells hit no compatible gate
at all.  Those `40` are not failures; they are already phase witnesses and
form `14` mirror-closed free-hole packets.  The terminal packet should split
rank-2 routed cells from free-hole cells.

The bypass is cleaner than expected: one pure `12`-cell legal mirror component,
six cells on branch0 and six on branch1, on the hard components `(43,54)`.
It carries owners `(23,93,113)`, while `(45,147,169,173)` remain seam-only
boundary debt.  This makes the next lemma look like a local discharge:
rank-2 routing plus free-hole packets plus pure bypass leaves no actual phase
flow needing the forbidden seam.

The guardrail is just as important.  Mirror pairs preserve class perfectly:
ordinary pairs with ordinary, free-hole with free-hole, and bypass with bypass.
Vertical half-turn gluing only appears on `48/282` cells and mixes free-hole
with ordinary in `14` vertical pairs.  So `n*2` supplies the address, but it is
not a legal quotient unless a sidecar records sheet occupancy and class.

Next experiments:

1. Run the same fiber-graph audit on every HYP-3477 hard mirror orbit, using
   the row-specific q=`14V` grid, to separate forbidden seams from genuine
   phase walls.
2. Lean-ize the finite random031 packet as three counts: `242` rank-2 routed
   cells, `40` free-hole cells in `14` mirror packets, and one pure bypass
   component.
3. Compare free-hole packets with the four dead islands from HYP-3482.  The
   counts differ, so the map is not literal equality; the likely relation is
   that dead islands are continuous punctures while free-hole packets are the
   q-grid shadow of the punctured complement.
4. Splice with incoming HYP-3490: private-label firewall explains why the
   projection-current carrier cannot discharge random031, while this fiber
   graph explains why the phase-flow carrier still does not need the forbidden
   seam.
