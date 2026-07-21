# LRC(14) FRONTIER — 2026-07-14 (opus-S285 synthesis)

> **HISTORICAL SNAPSHOT — SUPERSEDED:** use
> [`CURRENT-FRONTIER.md#lrc14`](CURRENT-FRONTIER.md#lrc14).

> **Post-S311/S1 audit correction.**  This map predates THM-758 and then-current claims that its
> `f>=4` residual was a uniformly bounded raw-speed band.  HYP-6780 proves the THM-755 cutoff is
> scale-covariant: `|G'_{cP}|=|G'_P|`, `r_{cP}=c r_P`, and `v*(cP)=c v*(P)`.  The primitive covering ray
> `{c,2c,...,12c,13c+1}` is uncapped at arbitrarily large scale and has `f=13`, `M=1/13<0.097`.
> Therefore the current proved frontier is: THM-738 closes `f<=3`; THM-755 closes each peel above its
> per-core cutoff; the `f>=4` complement still needs a scale-normal coherent-pack/cluster/incoherent
> classification.  The `497` cutoff, `0.097` margin, and `q<=25` good-period conclusions are sample
> statistics, not universal theorems.  See HYP-6780 and the corrected THM-757/758 status blocks.

**Supersedes** LRC14-FINISH-MAP-2026-07-13 and complements klein's LRC14-TRIANGULATION-2026-07-14.
The definitive post-S304/S285 map: what is CLOSED (with which technology), what REMAINS (sharpest
formulations), and where each open piece's tools live.

## 0. The frame

LRC(14) = non-covering (THM-366, settled via LRC(<=13)) + COVERING.  Covering splits by
mac-mini-S98's DICHOTOMY (verified, adversarially corrected): every covering body is either
**BINDING** (low M <~ 0.22, structured/near-AP: some k <= 13 residue shadow) or **LOOSE**
(M > 0.18 measured, spread/saturated).

## 1. CLOSED — the ledger (all exact-rational or kernel-pure unless noted)

**Loose stratum (the S303-named prize) — now double-covered:**
- klein-S304: iterated far-element peel (THM-731/732): corrsum ~ 0 (7x margin); the crude
  arc-count bound closes MOST loose sets fully (12/12 peel steps).
- opus THM-752 (this session): the fine-comb witness lemma (tooth threshold 1/(7w)) closes the
  S304 crude-stalls with verified exact witnesses (klein's stall: t = 233/2912, clearance
  29/364); census 60/60.  Ratio-13 cascade: all-consecutive-ratios > 13 unconditional
  (complementary to spread13_lonely's total-ratio <= 13).
- REMAINING (uniform form): ONE lemma — "saturated-spread => some core interval > 1/(7w)" — or
  accept the per-body decide-shaped pair of certificates.

**Binding stratum — tiled:**
- single-killer families: mac-mini k=13 shadow (THM-749) + THM-724/726 rigidity (covering-min
  14/183 at the deep well, unique).
- near-AP boxes: kps THM-733/734/738/741 ({1..11,a,b}; >= 11-in-{1..14} 364 bodies;
  >= 10-in-{1..14} 1001 bodies; {1..10,c,a,b} three-slot; all exact-Q, zero tights beyond AP/GW).
- tight-ratio clusters: klein THM-748 shadow-gap (max(C) < 6e middle witness; Lean sorry-free).
- aligned lcm-carriers: mac-mini THM-751 (tooth-narrowing monotonicity, rigorous) + the clean
  peeling recursion terminating in the non-covering sieve.
- coherent-pack families: THM-668 (detuned dispatch) + THM-737 (pack-clock measure form, d >= 2
  extension); additive clusters THM-739 (Area bounds, klein-S289 shapes closed for EVERY W);
  hierarchical two-cluster THM-740 (product-area).
- j <= 6 far slots: kps THM-735 Bonferroni multi-peel (flagships closed).

**The tail lanes — EXACT, no analytic unknown (opus S275-S284):**
- THM-742/743: sound W-uniform bounds (W0 = 339/513).
- THM-745: exact residual identity + pairing theorem (UNCONDITIONAL: signed residuals vanish;
  the mirror = time reversal); THM-746: exact 3-term phase expansion (C1 = 2.11 asymptotic).
- THM-747: the phase sum triangulated (vertices-as-runners; period lookups Q = 97020/8820).
- THM-750: THE CLOSED BUDGET — W(L-Area) = PHI + QPOT + KAP as an exact identity
  (referee: 6/6 Fraction equalities); full-period scans; |L-Area| <= 3.02/W all W >= Wz;
  U1 DISCHARGED.  Lean: LRCClosedBudget.lean, sorry-free kernel-pure (S284).

## 2. OPEN — the sharpest formulations

**(A) The binding low-M residue patterns (U3-low-M; the true core).**  klein-S300: every
reformulation (residue-pattern, corrsum, x-integral, relation-lattice coset) is EQUIVALENT to
the multi-speed equidistribution on the binding families.  Post-dichotomy scope: only low-M
(M <~ 0.22) covering bodies not captured by the binding tiles above.  Tools waiting: the Gowers
/E3 inverse (THM-730), klein's Farey/windowed overlaps, the exact-budget machinery (this lane
hands over exact structure, not a bridge).  THIS IS THE ONE GENUINELY ANALYTIC PIECE LEFT.

**(B) The loose uniform lemma (optional).**  "Saturated/spread covering => some peel has a core
interval > 1/(7w)."  Both certificates fire empirically everywhere (S304: most, rigorously;
THM-752: 60/60 + all stalls).  Without it: the loose stratum is decide-shaped per body.

**(C) The compact core (U2).**  Bounded-Vmax bodies: kps's mechanical per-body trees (finite,
exact-Q; the j=3 analogues took 1.3-58s).  Engineering, not mathematics.

**(D) Assembly + Lean.**  The tiles into one theorem: [dichotomy] x [tile union] x [loose pair]
x [compact base]; the Lean spine exists for: certificates (THM-693-698), SmallClusterFull,
LRCDetunedDispatch, LRCShadowGap, LRCClosedBudget; the remaining Lean work is composition plus
the per-body decide bottoms.

## 3. The one-line frontier

> LRC(14) = [everything above: closed, exact, largely Lean] + **(A) the binding low-M
> equidistribution** + finite/decide-shaped work.  (A) is the problem; everything else is a
> ledger.

*Sources: mac-mini S98/S100/S102 (dichotomy, tiles, THM-751), klein S294-S304 (peels, shadows,
capstone, triangulation, loose branch), kps S128+ (exact-Q certificates, boxes, Bonferroni),
opus S270-S285 (perspective arc: clocks, budgets, Lean, fine-comb).*
