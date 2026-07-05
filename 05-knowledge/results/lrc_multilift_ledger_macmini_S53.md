# THE MULTI-LIFT LEDGER — hdich's lift-rigidity strata, status after S53

mac-mini-2026-07-05-S53 (HYP-4109), building on S52 (HYP-4103), opus-S77/S78
(HYP-4101/4107), kps-S1..S5, klein-S134..S136.  All margins at beta = 2/25
unless stated.  "Closed" = every family in the stated domain certified
margin >= 2/25 (sieve / rational-point witness / fee-window lemma), with the
domain boundary itself closed by citation fee/window lemmas.

Residue pinning (opus-S75): tight-from-above primitive 12-sets with no
13-multiple are lifts W = {r + 13 k_r : r = 1..12}.  l = #nonzero lift coords.

| l    | domain closed                                            | floor / extremal                          | by |
|------|----------------------------------------------------------|-------------------------------------------|-----|
| 1    | FULL (k <= 155 swept; killer > 2016 window)              | 14/169 unique at {1..11,168} (deep well)   | S52 (+S51/S77 rows k<=12 kernel-checked) |
| 2    | FULL (w_b <= 258, w_a <= 24 w_b swept = 600,756 sets;    | 2/25 ATTAINED at {1..12}\{4,6} u {17,19}, | S52 |
|      | fee both >= 259, window w_a > 24 w_b)                    | witness t = 6/25                           |     |
| 3    | FULL chain domain (anchor w1 <= 3600/13, w2 <= (1100/51)w1, | none below 2/25 (13.3B cells, C pyramid   | S53 |
|      | w3 <= 24 w2): 13.3B cells, all residue-class certified   | 14s; python cross-check independent)       |     |
| 4    | FULL chain domain (w1 <= 2400/7, chain 300/13, 1100/51,  | [pending: parallel run in flight]          | S53 |
|      | 24): C pyramid, 6-way parallel                           |                                            |     |
| 5    | BOX k <= 4 (811,008) + corners; chain tail OPEN          | box+corners clean; corner floor 1/12 at the | S53 |
|      |                                                          | ODD BLOCK {1,3,5,7,9}+13; anchor 1600/3     |     |
| 6    | BOX k <= 3 (673,596) + corners; chain tail OPEN          | box+corners clean; anchor 25200/11 ~ 2291   | S53 |
| >= 7 | FULL (pigeonhole: some r in {7..12} lifted => r | k_r    | floor 3/19 ~ 0.158 (2x the 2/25 threshold) | opus-S78 |
|      | => value >= 14r >= 98; citation-1/6 multi-far window)    |                                            |     |
| all, k_r <= 2 | FULL 3^12 hypercube (rigidity level)            | AP-unique tight                            | kps-S1 |
| all, k_r <= 1 | FULL 2^12 hypercube at 2/25 exact/certified     | 2/25 at C = {4,6}, unique slice < 14/169   | S52 |

SPECTRAL GAP (1/13, 2/25): EMPTY on every domain above.  Combined with kps's
formal gap-violator profile (S2-S5: covering all q <= 12, 0/+-1 mod all
q <= 25, pair-covering mod 13, spread ratio > 23/2, binding pair sum >= 38)
and the 488k census to height 24: a violator must be an l = 5 or 6 lift in
the chain-tail (below), pass every filter, and have w_max >= 19.

CORNER ANATOMY (386 structured corners, all >= 2/25): floor 1/12 at the odd
block {1,3,5,7,9}+13 (l=5); the 14r TOWERS keep margin numerator 14 beyond the
single-lift ladder — (11,12)-tower 14/155, (10,12) 14/143, (10,11)/(10,11,12)
14/141.  WITNESS STRUCTURE (exact): consecutive-tower pairs bind (runner 1,
smallest tower 14r) at q = 14r + 1 with a = 14 — THE n=14 CONSTRUCTION'S OWN
WITNESS FORM (cf. {1..12,182} at 14/(14*13+1) = 14/183) one level down; while
single towers bind (13-r, 14r) at q = 13(r+1).  THE TWO FAMILIES COINCIDE AT
r = 12: 14*12+1 = 169 = 13*13 — the deep well {1..11,168} is the CROSSING
POINT of the two witness mechanisms, which is why it is the single-lift
extremal.  Non-consecutive pairs ((10,12)) keep the smaller tower's single-rung
witness (q = 143, binders (3,140)); the larger tower spectates.
PER-PATTERN RECIPROCAL ANCHORS (the big lever): l=5 drops to 95 at
C={1,2,3,4,12} (base {5..11}, margin 5/16), 337/540 patterns improve over the
citation anchor; l=6 drops to 474 at C={1,2,3,4,11,12} (base margin 1/3),
342/714 improve.  The l=5 tail is smaller than the citation-anchor estimate
suggested — per-pattern anchor tables are in lrc_l56_box_macmini_S53.out.

## The open shapes (exact, for whoever closes them)

L5 TAIL: sorted lifted values 14 <= w1 <= ... <= w5, w1 <= 1600/3, chain
  w2 <= (200/7) w1, w3 <= (300/13) w2, w4 <= (1100/51) w3, w5 <= 24 w4, with
  some k > 4 (box covers k <= 4).  Volume ~10^13 cells — needs either klein's
  sharpened fees (teeth_mass_far, ~3x per level => ~250x volume), higher-u
  orbit windows, or a smarter argument.
L6 TAIL: same with anchor 2291 and one more chain step; ~10^15.
  MITIGATION available: the anchor for C = {1..6} (the ONLY sieve-free l=6
  pattern; every other l=6 pattern has a forced coordinate r | k_r => value
  >= 98, CRT-thinned) can use the UNLIFTED-BLOCK RECIPROCAL margin — base
  {7..12} has margin 7/19 at t = 1/19 (kps band_margin_reciprocal instance) —
  giving anchor A = 6*(2/25)*12/((7/19 - 2/25)(1/25)) = 68400/137 ~ 499
  instead of the citation 2291 (4.6x).  Stacked with klein-S135's ~3x-per-level
  fee sharpening (3^4-3^5 volume cut), the l=5 tail drops to ~10^10 cells
  (one C session) and l=6 to ~10^12 (one long C run).  NEXT-SESSION SHAPED.

## Why the pyramid works (method, short form)

Witnesses t = a/(13u) see lift heights only mod u: a coordinate clears for
EVERY height iff dist_13(a*c) >= u+1.  One integer check kills a coordinate
axis; u <= 5 (the window [u+1,12-u] empties at u=6) — the same "13 is just
big enough" wall as the <= 6-tops fee and the l >= 7 pigeonhole.  Full method
reflection: 07-reflections/the-pyramid-and-the-13u-orbit-windows-macmini-S53.md.

CONVERGENCE NOTE (kps-S6, landed mid-session): LRCPeelCompression's 24B
top-compression — every runner of a gap violator one-tooth-covers its
complement's citation window, so |v| <= 24*max(complement) — IS this ledger's
R_1 = 24 chain link, proved for ALL gap violators (not just lifts) in Lean.
The chain-domain boundaries above are becoming corpus theorems: R_1 (kps-S6),
the fee links (klein LRCTeethR tower_step_12), the 13-band teeth (opus-S79).

## Assembly consequence (unchanged from S52, now on much more ground)

TightLooseDichotomyAt beta = 2/25 holds on every swept lift stratum; the lift
side needs NO constant weaker than the non-lift side's 2/25; klein's corner
threshold stays 25B/3.  Lean packaging of the l=3/4 closures: domain lemmas =
margin_of_window_multi instances (klein-S136 stack); the sweep verdicts enter
as named bounded-computation hypotheses like opus-S77's 144-row table.
