# AMM 12592 — Lane D1: the bulk rule (design sweep on the T6 junk flow)

Session: boxeph multifront, 2026-08-03. Exact int arithmetic throughout; no
floats in any decision. Script: `04-computation/amm12592_bulkrule_design_
sweep_boxeph.py` (extends the certified fast engine
`amm12592_transient_fast_junkflow_boxeph.py`). Ledger:
`amm12592_bulkrule_design_sweep_boxeph.out`, traces
`amm12592_bulkrule_trace_R*_D0*_<variant>_boxeph.json`, nu-profiles
`amm12592_bulkrule_nuprof_R512_D0{4,5}_plain_boxeph.json` and `..._D04_train_...`.

## 0. Setup: the generalized clamp is licensed by T2

T2's conjugacy proof uses only linearity and the [x^0]-matching condition:
it holds for ANY admissible in-box choice u_t (even, in the asymmetric c-box
[-2C(d-1,t), +2C(d-1,t-1)]), not just the nearest-point clamp. The engine
freedom per cell is exactly the interval I_t = [w_t - hi_t, w_t - lo_t] for
the junk jn_t = w_t - u_t: at an overflow cell the SIGN of jn_t is forced and
only the magnitude is free (can only take MORE junk, up to the full box width
2C(d,t) beyond the minimum); at an in-box cell junk of either sign can be
injected up to the distances to the box edges. Death iff jn_d != 0 (cell d
box is [0,2]; no choice can evade it), capture iff jn = {} after feed stops.

Key structural fact used by every variant: with next-kernel spread
S = 1 + delta' (K = (1,1) or (1,2,1)), the next-row load at output cell c is
w'_c = jn_c + 2 jn_{c-1} + jn_{c-2} (S=2 case) and is FINALIZED by its
lowest contributor jn_{c-S}; feed lands only at cells {0,1}, never at a
finalized output. So a single DESCENDING pass (front -> cell 0) can choose
each jn_t with full knowledge of the two cells above, i.e. exact one-row
lookahead at zero cost. An ASCENDING pass leaves the front outputs
uncontrolled (control lands below, where it is useless).

## 1. Variants implemented (all exact, all in one engine)

- plain: nearest-point clamp = rule A. REGRESSION-CERTIFIED bit-identical
  to the certified fast engine (outcome + full per-row trace + minus2 +
  front0) at (R,D0) = (256,0),(256,1),(512,4),(512,5).
- desc (V2+V3): descending pass; at each cell, if some jn_t in I_t makes the
  finalized output land INSIDE the next row's box (absorbable = zero junk
  there), take the feasible point nearest the plain overflow; else take the
  I_t endpoint minimizing the next overflow. VENT: cells above the front
  get an alternating halving pre-cancellation ladder (prefs
  (-1)^k excess/2^k) that fires only when the T9a edge transport of the
  front cell would overflow the next box.
- **desc1 (THE FINALIST)**: desc with the vent ladder truncated at depth 1
  (single vent cell tb+1, pref -excess/2). Reproduces every desc closure
  AND closes 256 at D0=0 (the depth-2 ladder pref +excess/4 at tb+2 breaks
  exactly that one closure — vent-shape sensitivity, see sec. 4).
- descw32/128/512: desc restricted to a window of that depth below the
  front (plain below) — separates bulk action from low-cell action.
- asc: ascending control (shows orientation matters).
- train / train16 / train3 / tr15 / tr54 / tr43 / tr74 (V1): desc + TRAIN
  GUARD in the top layer (t in [front-K, front], K = 16 or 48): when the
  chosen junk has the alternation-correct sign (opposite to the cell above)
  but magnitude below rho0 * |jn_{t+1}|, boost it (within I_t) to restore
  envelope ratio rho0; never flips a forced sign. rho0 in {3/2, 2, 3, 5/4,
  4/3, 7/4}. fill/fill15: additionally fill absorption holes inside the
  train with alternation-phased junk.

Verification battery on every run: per-cell box/parity asserts; streamed
independent identity spots at x = 2 and x = 3 (sum_i x^i Delta_i(x) =
(1-x)^{R-1}, computed from the u-ledger, not the flow state) + ballot debt
= (R-2)/2; for closures at R <= 1024 the full blocks are rebuilt
(ballot + u) and passed to the hostile referee `verify_witness`
(admissibility + full polynomial identity via t-substitution + ballot laws).

## 2. RESULTS: exact D0* thresholds (every value pinned by an adjacent
die/close pair; closures spot-checked, finalists referee-verified)

| R    | plain (known) | desc1 (finalist) | tr15   | best of family | gain |
|------|---------------|-----------|-----------|----------------|------|
| 128  | 0             | 0         | 0 (cap 81)| 0              | 0    |
| 256  | 1             | **0** (cap 154, REFEREE ok, witness sha 5950bd42) | 0 (cap 153) | **0** | 1 |
| 512  | 5             | 5 (4 dies @123) | **4** (cap 319, REFEREE ok; 3 dies @119) | **4** | 1 |
| 1024 | 15            | **14** (cap 628; REFEREE ok + witness sha 81391349 on the desc run, desc1 closes identically; 13 dies @245) | >= 16 (15 DIES @542!) | **14** | 1 |
| 2048 | 38            | **36** (cap 1271; 35 dies @503) | 37 (36 dies @1051) | **36** | 2 |
| 4096 | 89            | **87** (cap 2537; 86 dies @1019; pinned with desc1) | not run | **87** | 2 |
| 8192 | 192           | **<= 188** (desc1, cap 5064, spots ok; desc closes 189 too) | — | **<= 188** | >= 4 |

(At 512/1024/2048/4096 desc and desc1 have identical outcomes, die rows,
capture rows, minus2 and T6b on every pair tested; the two variants differ
only at 256 D0=0.  A desc1 probe at R=16384 D0=400 — the exact point where
plain DIES at row 4055, refuting the 25/1024 saturation — was launched at
close-out; its trace flushes to
`amm12592_bulkrule_trace_R16384_D0400_desc1_boxeph.json`.  If it closes,
gain >= 1..16 at 16384 and the desc1 curve stays strictly below plain
through seven doublings.)

All closures have minus2 = (R-2)/2 exactly (T8 debt paid), capture at
~0.62R, spots x2/x3/debt all True. Witness file (with SHA-256):
`04-computation/amm12592_witness_R512_bulk_tr15_D0_4_boxeph.json`
(sha256 311975f6..., 1.19 MB), referee ALL-PASS. desc 256/512 closures and
tr15 512 closure and desc 1024 D0=14 closure referee ALL-PASS
(ok/adm/identity all True, debt 511 = (R-2)/2).

Reading of the curve: best-of-family normalized 1024*D0*/R =
0 (256), 8 (512), 14 (1024), 18 (2048), 21.75 (4096), <= 23.5 (8192) —
increments 8, 6, 4, 3.75, <= 1.75 — the same declining-increment LINEAR
signature as plain (4, 10, 15, 19, 22.25, 24), shifted down. The local
bulk-rule class does NOT produce flat or O(log R) slack on this data; it
uniformly dominates plain (desc/desc1 never worse anywhere measured) and
shaves the threshold by 1..4+ units with the absolute gain GROWING in R
(0,1,1,1,2,2,>=4 — lower 8192 probes not yet run, so >=4 is a floor) while
the relative gain stays ~2%: consistent with a boundary-layer effect on an
unchanged slope, though the growing absolute gain leaves a slope reduction
not yet excluded. eps_inf estimate for the family: <= 0.023 vs plain
~0.0245 (both still linear on this data).

## 3. MECHANISM (the real yield of this lane)

Exact anatomy of the die/close dichotomy (nu-profiles; nu := (1+sigma)j =
pairwise-sum residue = non-alternating content):

1. **Alternation signature.** In every CLOSE trace the front-region junk
   envelope stays strictly sign-alternating with per-cell magnitude ratio
   ~2-4 growing inward, and the front recedes. In every DIE trace, at the
   dichotomy rows a SAME-SIGN BLOCK nucleates a few cells behind the front
   (R=512 D0=4 plain: rows 10-15; birth sequence at row 27 in the train
   trace: two adjacent same-sign cells -> 5-block -> block pair
   (+block, -block) -> marching spike). The marching front value is the
   block's leading value, FROZEN (T9a edge transport), e.g. 292 bits
   (plain) / 323 bits (train) at R=512.

2. **Immortality threshold.** Once a same-sign block's value exceeds the
   local cap scale, no in-box rule can dissolve it: absorption removes at
   most ~2C(d-1,t) per cell-row, caps along the march diagonal DECAY
   geometrically (T9b: tau > 1/phi), so the total remaining absorption
   budget along the diagonal is a bounded multiple of the first cap, while
   the block value is hundreds of bits above it. All shaping must act
   BEFORE the block crosses cap scale — i.e. in the top layer where junk =
   load - cap is a difference artifact with |junk| << cap (full sign
   control). This is why the dichotomy is decided in rows ~20-70 and why
   only the top ~2-3 cells are absorbable per row (excess/cap ratio grows
   ~x1.5 per cell downward from the front: the boundary layer is O(1)-thin
   at fixed D0, and D0 shifts it multiplicatively).

3. **Ratio-2 marginality law (new, exact).** For an alternating train with
   inward envelope ratio rho: kernel (1,2,1) maps the cellwise value to
   (1-rho)^2 x (own cell), kernel (1,1) to (1-rho) x — so rho = 2 is
   EXACTLY magnitude-preserving ((z+1)^2 resolvent), rho < 2 decays,
   rho > 2 grows with coherent phase FLIP on (1,1)-rows (dominance from the
   interior side). This explains: (i) rho0 = 3/2 guard beats rho0 = 2
   (which beats 3) at R=512 — sub-doubling envelopes are the annihilating
   regime; (ii) the guard's failure mode at scale: maintained ratio-2
   trains never decay, junk L1 GROWS post-feed (R=1024 D0=15 tr15: peak
   1366 bits vs 982 plain, death at row 542 where plain closes — the guard
   can BREAK a plain closure), and defects nucleate where the boosted layer
   meets the natural envelope.

4. **Where the gain comes from.** desc's gain (1..3 units) is exactly the
   absorb-window boundary layer: choosing u inside the box (not at the
   nearest point) absorbs the front taper ~1-2 cells deeper per row and
   pre-cancels the T9a edge transport (vent). Window-depth tests (32/128/
   512/full) change nothing at 512 — confirming only the front boundary
   layer matters for the die and low-cell shaping is irrelevant to it
   (consistent with B1 steering-invariance).

5. **Negative (hazard-disciplined).** Within the class "one-row-lookahead
   local clamp rules on the T6 flow", no variant tried achieves sublinear
   slack; deaths of sub-threshold variants are NOT evidence about epoch
   feasibility. The class's structural obstruction is the immortality
   threshold: any rule that lets a same-sign block cross cap scale near the
   front loses; preventing nucleation requires controlling the interplay
   of feed-repumped loads with the (1,1)/(1,2,1) Beatty kernel word over
   MANY rows — i.e. scheduling across rows (choosing WHERE the residue
   accumulates over the whole feed phase), not cellwise-local choices.

## 3b. THE SUPER-BLOCK LAW (new, exact; the sharpest result of this lane)

Mode `block` of the sweep script: per row, take the longest same-sign junk
run within 40 cells of the front; call it SUPER when its MINIMUM |value|
exceeds the ENTIRE remaining cap tail sum_{t' > lead} 2C(d-1, t') (the
total absorption budget strictly beyond the run's deepest cell). Exact
outcomes on all near-threshold pairs (three rules, four scales):

| run | outcome | first super row | super rows after | first run>=5 |
|-----|---------|------|------|------|
| 512 D0=4 plain | DIE @121 | 13 | ALL (109) | 11 |
| 512 D0=5 plain | CLOSED | never | 0 | 82 (sub-cap) |
| 1024 D0=14 plain | DIE @250 | 17 | ALL (232) | 18 |
| 1024 D0=15 plain | CLOSED | never | 0 | 160 (sub-cap) |
| 2048 D0=37 plain | DIE @508 | 27 | ALL (482) | 25 |
| 2048 D0=38 plain | CLOSED | never | 0 | 320 (sub-cap) |
| 1024 D0=13 desc1 | DIE @245 | 16 | ALL (228) | 18 |
| 1024 D0=14 desc1 | CLOSED | never | 0 | 187 (sub-cap) |
| 512 D0=3 tr15 | DIE @119 | 13 | ALL (106) | 15 |
| 512 D0=4 tr15 | CLOSED | never | 0 | 57 (sub-cap) |
| 4096 D0=88 plain | DIE @1014 | 31 | ALL (984) | 28 |
| 4096 D0=89 plain | CLOSED | never | 0 | 630 (sub-cap) |

Perfect 12/12 separation: **death <=> a super block forms** (and it forms
by row 13/17/27/31 at R = 512/1024/2048/4096 — a DECLINING fraction of
R, LONG before the visible T9' march); once super,
super on EVERY later row (no recovery ever observed, 2141 super rows
total); closures contain long same-sign runs too, but always sub-cap. The
first-super row is the true point of no return; the die row is just its
delayed consequence. This (i) explains why all local shaping fails (by
first-super time the block already dwarfs every box it will ever meet),
(ii) gives the E-lin scaling limit its order parameter — sup over rows of
(block min value)/(cap tail) as a function of eps = D0/R — and (iii)
upgrades the immortality statement to a concrete provable target:
"super-block persistence under every admissible clamp" (the kernel adds
same-sign contributions inside the block faster than the sub-tail boxes
can drain; T9b gives the geometric tail). Ledgers:
`amm12592_bulkrule_blockledger_R*_D0*_*_boxeph.json`.

## 4. Rule non-monotonicity + spec sensitivity (new hazard entries)

Rule-family thresholds are NOT monotone in D0 or in variant strength:
train at 512 dies at D0=4 (row 169) but at D0=3 dies LATER (row 170);
tr15 closes 512 at D0=4 but dies at 1024 D0=15 where plain closes. Any
future claimed threshold for a designed rule must be pinned by explicit
die/close adjacent pairs (as done here), never inferred from monotonicity.

Spec sensitivity (caught by re-run): extending the vent ladder from depth 1
to depth 2 (adding pref +excess/4 at cell tb+2) FLIPPED the 256 D0=0
closure to a death at row 63 with everything else unchanged. The finalist
spec desc1 is pinned in `VARIANTS` (vlad=1) and all table closures above
reproduce under the committed script; the near-threshold dynamics is
chaotic in the rule parameters even though every individual run is exact
and deterministic.

## 5. Handoff to D3

- The exact spec of the best rules is in `VARIANTS` + `choose_row` of the
  sweep script (desc = absorb-window + vent ladder; tr15 = + guard(16, 3/2)).
- Exact witnesses saved as files: R=256 D0=0 (desc1, sha 5950bd42),
  R=512 D0=4 (tr15, sha 311975f6), R=1024 D0=14 (desc, sha 81391349);
  R=2048 D0=36 / R=4096 D0=87 / R=8192 D0<=188 closures have traces +
  x2/x3/debt spots (witness JSONs regenerable deterministically from the
  pinned specs).
- The march/stall dichotomy is now a statement about same-sign-block
  nucleation vs cap scale in the front boundary layer. The scaling-limit
  question (D2/E-lin) should be phrased in these variables: block value /
  cap ratio at nucleation as a function of eps = D0/R.
- Attainment: with best-of-family, T(n) <= n + 1 + floor(gamma* n) + 14
  for critical n <= 2047 (1024-epoch slack 14 replacing 15) — the desc 1024
  D0=14 closure passed the full hostile referee (ok=True) and both identity
  spots; via the THM-3329 assembly this improves the standing attainment
  constant from 15 to 14.

## 6. Status labels

- PROVED: generalized-clamp conjugacy (T2 extension — proof identical,
  choice-independent); finalized-output structure of the descending pass;
  ratio-rho kernel action law ((1-rho)^2 / (1-rho) cellwise multipliers).
- VERIFIED-exact: regression vs certified engine (4 configs, bit-identical);
  D0* table above with adjacent die/close pairs; referee ALL-PASS for
  desc1 R=256 D0=0, tr15 R=512 D0=4, desc R=1024 D0=14 (witness files with
  SHA); x2/x3/debt spots on every closure listed; super-block ledger 12/12.
- MEASURED (no proof): die/close alternation signature; block-birth
  sequence; gain pattern 0,1,1,1,2,2,>=4; eps_inf(family) <= 0.023.
- CONJECTURED: immortality threshold as an exact lemma (block value >
  geometric cap-tail along the T9b diagonal implies death for every
  admissible clamp rule) — this looks provable with the T9 machinery and
  would upgrade the class-negative from measured to proved.
