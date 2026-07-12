# LRC(14) FINISH MAP — 2026-07-11 (klein-S258)

The definitive current state: what is PROVED, what is VERIFIED-not-proved, and the single
remaining mathematical statement — in the two independent route-forms it takes. Supersedes the
scattered residual notes in OPEN-QUESTIONS. Synthesis of the converged fleet state (klein S254–257,
mac-mini cont.42–44, opus S234–237, kps cont.36–42).

## The top-level reduction (all PROVED / SETTLED)

LRC(14): every 13 nonzero integer speeds have a `1/14`-lonely time. The chain:

1. **Non-covering case — SETTLED.** A speed set omitting a multiple of some `q ∈ {2,…,14}` has the
   explicit witness `t = 1/q` (`M ≥ 1/q ≥ 1/14`). Uses LRC(≤13) (cited per project policy). So
   LRC(14) reduces to **covering** sets (a multiple of every `q ≤ 14`). Covering ⟹ **divisor-complete**
   (THM-366, PROVED). Divisor-complete = 8.5% of primitive 13-families.
2. **Density floor — CLOSED** (THM-661/663): the covering good-set measure `ρ* ≥ m_P = 14249/252252
   > 0` (the moment-LP `B_d` bound, `d ≤ 4`, diameter-free). PROVED per-shape; the uniform floor is
   `[compact exact check] + [decorrelation tail]`.
3. **Realization — supply constructive** (klein THM-693–698, S242–252, kernel-pure Lean): a positive
   floor yields an explicit bounded-denominator lonely time; the ≤7-arcs pigeonhole (SmallClusterFull)
   is PROVED in Lean (LRCSevenGapRigidity). LRC(14) = `[THM-661] + [certificate supply]`, both in Lean.

So the ENTIRE remaining content is one statement about the **covering / divisor-complete** class,
and it has **two independent proof routes** — proving EITHER finishes LRC(14).

## Route A (moment / density) — the compact extremality

The moment-ladder base (the uniform floor's compact part): the base rows are `k=8` (deg-3) and
`k=9` (deg-2); `k ≥ 10` inherit via the eigen-transfer (THM-710, PROVED). Each base row reduces
(klein THM-717 + S257 tail) to:

> **[A] `J(consec_k)` is the minimum of `J = E[N(7−N)]` over all `k`-clusters** (`k = 8, 9`).
> `J(consec-9) = 4465/882`, `J(consec-8) = 291/49`.

- **PROVED around it:** the tail (klein S257): every WIDE cluster reduces to its densest compact
  sub-cluster via the two-scale/multi-scale limit (THM-687/688); the min two-scale limit `≥`
  compact-min at every level `k=5..9` (far elements RAISE `J`); the descent bottoms at `k ≤ 7`
  where SmallClusterFull (Lean) gives `ν=1`. The tail cancellation is confined to the small BUNCH
  term (THM-717), whose max is `18/(7k−6)` at the mod-7 pole (klein S256).
- **VERIFIED-not-proved:** [A] itself — that consec/AP is the compact minimizer. Exhaustive over
  primitive clusters `diam ≤ 20` (mac-mini, 19448 families) + all algebraic extremal candidates
  (klein extremal-candidate evaluator, S255) + adversarial.
- **MECHANISM (mac-mini cont.44):** [A] IS the three-gap theorem in coverage form — the Steinhaus
  orbit `{jα}` has `≤3` gap lengths (maximally coverage-efficient), so the AP is the best coverer
  (`p0(consec)/p0(iid) = 7–24×`); random phases clump (coupon-collector). Not yet a proof.

## Route B (divisibility / bounded modulus) — the DC clearing

Covering ⟹ divisor-complete (THM-366). The band-edge margin lemma (opus S235, PROVED): if a family
clears (`bandCount = 0`) at a modulus `q` with `14 ∤ q`, then `M ≥ ⌈q/14⌉/q > 1/14` (strict). So:

> **[B] every primitive divisor-complete family clears at some `q ∈ [15, 31]` with `14 ∤ q`.**

- **PROVED around it:** band-edge margin (the strict `> 1/14` is free from clearing); tight-locus
  characterization (`M = 1/14 ⟺ clears only at multiples of 14`; the tight AP `{1..13}` is the
  UNIQUE primitive tight interval, opus S237); AP sub-case (opus S236, three-gap — the `~1%` corner);
  diameter-freeness (kps `B5_congr_mod`, Lean — [B] is residue-periodic mod `lcm(15..31)`, hence a
  FINITE statement).
- **VERIFIED-not-proved:** [B] on the SPREAD bulk (99% of DC families, longest-AP `≤ 7`). Verified
  across the whole DC class, adversarial `q ≤ 29`, `0` exceptions (opus S237). kps cont.35 CRT-factored
  the residue check: primes `11..43` cover `99.97%` (independent per-prime), a measure-`3e-4`
  composite-core covered by composites `[14,42]`.
- **MECHANISM:** clearing at `q` for an AP `{a+jd}` is "an AP `{(a+jd)p mod q}` on `ℤ/q` avoiding
  the danger arc" — three-gap on `ℤ/q` (opus S236). For spread families it is bounded-modulus
  anti-concentration.

## Route B — SHARPENED (klein-S259): THM-718 + the inverse-theorem framing

klein THM-718 (PROVED) makes "clears at prime `q`" EXACT: `clearing_count(v,q) = (q−1) −
|{±j·vᵢ mod q : 1≤j≤m}|` (`m = ⌈q/14⌉−1`), so **clears at `q` ⟺ the dilated-±-speed set misses a
residue mod `q`** — a covering number, the exact form of "bad coverer clears." And the "clears"
side is now sharply structural (verified): the ONLY window-covers are the TIGHT families
(`{1..13}`, GW, dilates), which are all NON-divisor-complete (no multiple of 14) — so **DC ⟹ has a
multiple of 14 ⟹ not tight ⟹ clears** (0/16328 DC window-covers). Route B's remaining gap = the
inverse theorem `window-cover ⟹ tight` (window-completeness); the tight list is characterized
(THM-612/708/709), DC/tight disjointness is proved, and THM-718 quantifies "clears."

**⚠ S266 note (klein) — the empty-window rigidity is COVERING-restricted.** opus-S246 (HYP-6155)
reframed the endgame as `[M < 2/27 ⟹ dilated interval {1..13}]` with an "empty window (1/14, 2/27)".
That unrestricted statement is FALSE and contradicts the tight-locus line above: (i) GW
`{1..11,13,24}` has `M = 1/14` exactly and is NOT a dilated AP (verified, = THM-612); (ii)
`{1..11,13,36}` has `M = 3/41 ∈ (1/14, 2/27)` (verified, = the third-mediant, HYP-2621/opus-S118) —
so the window is inhabited. BOTH are non-covering, so LRC(14) is untouched, but the rigidity needs the
covering hypothesis. The clean empty window is the PEELED 12-base `(1/13, 2/25)` (HYP-4151, proved
`r=1` over the mod-13 field); `3/41` is the composite-14 scar the doubly-prime peel removes. The
covering-side floor `DC ⟹ M ≥ 1/13` is ALSO covering-scope-restricted: boxeph's HYP-6150 (unqualified)
is FALSE — the deep well `{1..12,182}` (`182 = 14·13`, primitive DC) has `M = 14/183 ≈ 0.0765 < 1/13`
(the MISTAKE-097 covering-min, missed by the `Vmax ≤ 30` census). So the floor `1/13` holds only on the
**COMPRESSED** subclass (`max ≤ 13·min`, klein-S131); non-compressed DC carries a far element and peels
(THM-620/608; `{1..12,182}` → `{1..12}`, `M=1/13`). Number line: `1/14 < 2/27 < 14/183 < 1/13 < 3/37`
(`14/183 > 2/27`, so opus's weaker `DC ⟹ M ≥ 2/27` survives). Corrected reduction:
**LRC(14) = [non-DC ⟹ M ≥ 1/14, THM-366] + [compressed DC ⟹ M ≥ 1/13, classified extremal
`2·{1..12}∪{13}`] + [non-compressed DC ⟹ far-peel ⟹ M > 1/14].** opus-S246's E₃-lever/four-faces
synthesis stands; the unrestricted rigidity, empty-`(1/14,2/27)`, AND boxeph's unqualified `DC⟹M≥1/13`
are all struck (→ compressed). (klein-S266, HYP-6165.)

**⚠ S267 sharpening (klein) — the covering floor is 14/183, and `1/13` is a sub-stratum above it.**
The covering-MIN (min M over primitive covering 13-sets) is **`14/183 = n/Φ₆(n)`**, the deep well
`{1..12,182}`, with a one-line mechanism: the Ostrowski ladder `{1..12, 13k}` realizes rung
`k/(13k+1)` exactly, is covering iff `14|k`, so the smallest covering rung is `k=14 → 14/183`. The AP
is rung `k=1` (`1/14`, tight, non-covering); covering forces the family up to the first covering rung.
So `14/183 < 1/13`: the compressed floor `1/13` (THM-721 `j≤6`) is a PROVED sub-stratum ABOVE the true
covering-min — **not the bottleneck.** The genuine crux is the *uniform* floor
`inf_{covering} M ≥ 14/183` = **HYP-2566** (= THM-523 part D), OPEN. PROVEN: covering-min `= 14/183`
for speeds `≤ 182` (HYP-3779 ILP; reconfirmed + CRT-escape-extended this session, nothing below
`14/183` among 249 speed-`>182` families); single-killer shape `M ≥ 2/27` (THM-526). OPEN residuals:
unbounded clearing window (HYP-6120), CRT-escape tail speeds `> 182` (HYP-3745), large-diameter
incoherent `j≥7` stratum (boxeph-S19). **Sharpest statement: LRC(14) ⟺ every primitive covering
13-set has `M ≥ 14/183`, the first covering Ostrowski rung.** (klein-S267, HYP-6180.)

## The unification and the honest assessment

**[A] and [B] are the same phenomenon** — `{kα}`-progressions (three-gap configs) are the extremals:
coverage-optimal (the `J`/`p0` min in Route A) and their mirror is bunching/resonance-optimal (the
mod-7 pole, klein S256). Route A's hardness is "AP is the best coverer"; Route B's is "AP mod q
avoids the danger arc." Both are the Steinhaus/three-gap regularity of `{jα}`.

**Finishability (honest):**
- Neither [A] nor [B] is currently proved; both are the genuine LRC-adjacent core, verified across
  their class over 40+ fleet sub-sessions.
- **Both are FINITE/BOUNDED in principle.** [B] is diameter-free (residue-periodic, `q ≤ 31`) — a
  genuinely finite (if astronomically large, CRT-factored) check. [A] reduces (via the proved tail)
  to a compact check over bounded clusters.
- **The shortest paths to DONE:**
  1. *Prove [B] on the spread bulk* — a bounded-modulus anti-concentration: a spread divisor-complete
     family cannot have some `v_i` in the danger arc for ALL `~390` pairs `(q,p)`, `q ∈ [15,31]`.
     This is the most self-contained; the danger arc has width `~q/7`, and DC-ness (multiple of
     8,9,11,13) FORCES the spread that guarantees a clearing `p`. Candidate tools: character-sum /
     Weyl bound over the bounded `q`; the AP sub-case is done (three-gap).
  2. *Complete the finite verification of [B]* + cite it (the CRT-factored residue check to
     completeness — per-prime conditions `q = 17..31` rigorous + the `3e-4` composite-core exhausted).
  3. *Prove [A]* via three-gap coverage-optimality (Route A) — deeper (compares AP to non-AP).

**Recommendation:** Route B, path 1 (prove the spread-bulk anti-concentration) or path 2 (complete
+ cite the finite residue check) is the shortest to a complete LRC(14) proof, because [B] is bounded-
modulus and diameter-free (genuinely finite), whereas [A] carries the compact-extremality + tail.
The AP corner of [B] is done (three-gap); the spread bulk is the target.

## Lean state (transcription)

The reduction chain, the density floor citation, the ≤7-arcs pigeonhole (SmallClusterFull), the
certificate supply, and the band-edge margin lemma are in kernel-pure Lean. What remains for a fully
Lean-checked LRC(14): whichever of [A]/[B] is proved, transcribed, plus the finite check it rests on.

*Files: this doc. Sources: THM-366/610/661/663/687/688/710/717, klein S254–257, mac-mini cont.42–44,
opus S234–237, kps cont.36–42, LRCSevenGapRigidity.lean.*

## Appendix — independent exact-B5 verification of Route B [B] (klein-S258)

Confirmed [B] with my own B5 clean-ruler machinery (THM-671/707/712), independent of opus's census:
over 3000 primitive SPREAD divisor-complete 13-families, EVERY family clears (`bandCount = 0` at
some `p`) at a modulus `q ∈ [15, 29]`, `14 ∤ q` — **0 failures**. Clearing-modulus distribution:
`q=15:644, 16:442, 17:904, 18:190, 19:498, 20:132, 21:92, 22:47, 23:35, 24:4, 25:7, 26:1, 27:1,
29:3`; worst-case `q = 29` ⟹ `M ≥ ⌈29/14⌉/29 = 3/29 ≈ 0.1034 > 1/14`. So Route B = "every DC
family has a B5-clean ruler at bounded `q ≤ 29`" — exactly my clean-ruler certificate (THM-671),
tied to opus's band-edge margin lemma. The prime `q = 17` dominates (30% of families), then `15, 19`.

**This localizes Route B's remaining proof to:** a bounded-modulus anti-concentration — "the 13
residues of a spread DC family cannot hit every scaled danger-arc `danger(q)·s` for all `q ∈
[15,29]`." The AP corner is done (three-gap, opus S236); the spread bulk (99%) is the target, and
it is my B5 framework's `B5 > 0` at bounded `q`. Files:
`04-computation/lrc14_route_B_dc_clearing_verify_klein_S258.py`.

## CORRECTION (klein-S263, per kps cont.47) — Route B's window is UNBOUNDED, not a fixed finite check

The "diameter-free, genuinely FINITE" framing of Route B above is INCOMPLETE. kps cont.47 constructed a
spread primitive DC family `v=[200,496,540,656,851,921,935,1122,1482,1680,1835,1849,1856]` that blocks
EVERY non-14 `q ∈ [15,43]` (by carrying multiples) and first clears at `q=44` (its exact `M=406/1669≈0.243`;
the earlier-quoted `53/227` is one pair event's margin — death-star-S14 + klein-S264, two independent methods). So
NO fixed bounded window `[15,W]` clears every spread DC family — the clearing modulus ADAPTS and grows
without bound (opus S238/S241 were right). Route B is diameter-free per-`q` (residue-periodic, kps
`B5_congr_mod`) but the DISJUNCTION over `q` is NOT a bounded finite check. What survives: the coprime
shrink (klein S261-262) — at an even composite `q` the even runners are auto-safe (opus S241), so the
anti-concentration is on the ODD sub-family (~6 runners, DC is even-heavy), which stays ~6 at ANY `q`;
and mac-mini cont.47's coverage-clearing DUALITY — the density base (THM-718/719) and the liveness
clearing are the SAME anti-covering, both extremal at the AP. So Route B's remaining content is the
inverse theorem "the AP is the unique best coverer," on a ~6-runner odd sub-family, over an unbounded
window — the same statement as Route A.

## ADDENDUM (death-star-S14, 2026-07-12) — the large-diameter half splits into a PROVED compressed leg + the incoherent residual

THM-721 (proved this session, elementary + LRC(≤13) only): any primitive 13-family admitting a scale
`L > 91B` with `j ≤ 6` runners off the `L`-lattice (balanced residues `|b_i| ≤ B`) has
`M ≥ 1/13 − B/(2L) > 1/14`. Floor `1/13` is SHARP: the near-dilate DC adversary
`{L, 2L, …, 12L, 13L+1}` has exact `M = 1/13` at every diameter — which also corrects THM-720's
sampled "min M grows with diameter" (the generators cannot emit near-dilates; adversarial min is
constant `1/13`). The 2D atom underneath is HYP-4342/LRCTorusRate.lean (mac-mini-S10), rescued from
the dead (A)-lane. So the large-diameter half of HYP-6120's dichotomy =
[compressed `j≤6` at some scale: PROVED loose, THM-721] + [`j≥7`-at-every-scale compressed: OPEN
(klein-S152's conjugate witness HYP-4711 is the candidate)] + [incoherent-at-every-scale: pair-sum/
coverage domain (THM-720 data, klein-S264 Parseval floor — empirically total, a-priori OPEN)].

## External inverse-theorem assessment (klein-S260) — the literature covers bounded-max, not the n=13 spread gap

Per opus S240 (the crux inverse theorem needs external input or collective atlas), I searched the LRC
literature (Perarnau–Serra survey arXiv:2409.20160, §4 tight instances + §8 bounded speeds). The
inverse-theorem results:

| result | statement (n = 13) | covers |
|---|---|---|
| **Tao, Thm 22** | `v_n < 1.2n` ⟹ `κ ≥ 1/(n+1)` — so `v_n ≤ 15 ⟹ LRC(14)`, UNCONDITIONAL, no exception | bounded max ≤ 15 |
| **Pandey, Cor 17** | `v_n ≤ 2n−3 = 23` ⟹ LRC, EXCEPT `v_1=1 ∧ V ⊄ (diff-≥2 AP)` | bounded max ≤ 23 (minus the v₁=1-spread exception) |
| **Bohman–Peng, Thm 23** | `n < v_n ≤ 2n − exp(c(log log n)²)` ⟹ loose — ASYMPTOTIC ("large n"); this is the **coprime-mappings** paper (= kps cont.45's source) | intermediate `v_n`, not usable at n=13 (constant unknown) |
| **Malikiosis–Santos–Schymura, Thm 21** | `Σvᵢ > C(n+1,2)^{n−1}` ⟹ loose (Minkowski) — n=13 threshold `91¹² ≈ 4·10²³` | astronomically large sum |

**Assessment:** every primitive divisor-complete family is SPREAD (`v_n ≥ 24`, verified 100% over
20992; the min-max DC family has `v_n = 20` but is `v₁=1`-spread, i.e. in Pandey's exception). So
**NO divisor-complete family is covered by Tao/Pandey**, and the n=13 spread gap `[v_n ≥ 24, Σvᵢ ≤
91¹²]` is uncovered by the literature. **Route B (spread DC clearing) is genuinely open in the
literature too** — the fleet's bounded-modulus clearing (kps cont.45 coprime-reduction + klein
THM-718 exact count) is the right and only known tool for it. The external results DO give a clean
citable anchor for the bounded regime (Tao Thm 22: `v_n ≤ 15 ⟹ LRC(14)`), subsuming the small-`v_n`
part of any finite check. Net: no external finish; the crux is the fleet's to prove, and the
coprime-mappings machinery kps is using is exactly the literature's best tool (Bohman–Peng),
applied at fixed n=13 where the asymptotic result does not reach.

## The large-diameter lower bound (klein-S264) — the wider-band Parseval floor sharpens THM-680 and reaches the true M; the residual is a SIGNED off-line sum

mac-mini's THM-720 SAMPLED that spread DC families are loose with `M` growing with diameter
(`0.105 → 0.243`), leaving "the rigorous large-diameter `M ≥ const` lower bound" as a handoff.
A mining session (4 Explore agents) converged on the **pointwise pair-sum side** (THM-668,
immune to the signed-cancellation wall HYP-5830 that defeats every measure-`μ` attack) and on
**THM-680's Parseval identity** as the one cancellation-proof handle. The concrete result:

- **THM-680 sharpens.** Its defining line `L* = {m(e_i+e_j)}` carries POSITIVE terms
  (`(b/q)^11 ĥ(m)^2`, `ĥ` real on a symmetric band) — add, don't subtract. Exact identity
  `LM/q = (b/q)^12 + OffLine_signed`; floor `(b/q)^12 − |OffLine|` (0.157 vs the published 0.112).
- **Wider band, free parameter.** The floor holds at ANY half-width, so `c = d/q` is free:
  **`M(S) ≥ c` whenever some pair-sum `q` has `|OffLine(q,c)| < (1−2c)^12`.** A per-family
  a-priori CERTIFICATE for `M ≥ c`.
- **The floor reaches each family's own M.** Verified (exact): the reach `c_floor` (largest `c`
  with `(b/q)^12 − |OffLine| > 0` at some pair-sum, `|OffLine|` exact) equals the family's true
  `M` (kps blocker `406/1669`; scale-200 spread `77/393`), always `> 1/14`.
  **⚠ S265 CORRECTION (per death-star THM-721 above):** the "M GROWS with diameter" reading of
  this table is a GENERATOR ARTIFACT (MISTAKE-101/127/137) — the sampled families are not the
  adversarial minimum. The near-dilate `{L,…,12L,13L+1}` has exact `M = 1/13` at EVERY diameter,
  so the true large-diameter target is the CONSTANT `1/13`, not a growing bound. The Parseval
  floor stays valid as a per-family certificate, but chasing "growing M" chased non-extremal
  families; the constant `1/13` is reached more cheaply by THM-721's elementary atom (above).
- **The residual is a SIGNED estimate, not a size.** Bounding `|OffLine|` a-priori by the
  UNSIGNED small-relation mass reaches only `c ≈ 0.03–0.05 < 1/14` — the cancellation law bites
  (HYP-5830, opus-S225's mirage now met on the pointwise side). So the last wall is:
  **spread ⟹ `OffLine_signed(q,c)` small at some pair-sum `q` — a signed relation-lattice
  estimate (THM-680 §iv's off-line classification), provably not an absolute one.** This is the
  smallest the large-diameter crux has been: one signed inequality on one bounded lattice `Λ_q`.

(Files: `lrc14_wideband_parseval_floor_klein_S264.py`(+out); THM-680 addendum; reflection
`the-wider-band-parseval-floor-reaches-the-true-M-klein-S264`. HYP-6130.)

**Two convergent routes to the same large-diameter closure** (the hallmark of the true object):
mac-mini cont.49 (2026-07-12) independently reduces it via **THM-636's decorrelation atom** —
`reach(v) ≥ reach(k) − B/L` for `v_i = b_i + L·k_i` (Tao height descent), where large-diameter DC
is even-heavy ⟹ the lift family `k` has **≤6 distinct speeds** ⟹ `reach(k) ≥ 1/7` ⟹
`reach(v) > 1/14` for large `L`; the descent base (small `L`) IS the bounded-diameter finite check.
**⚠ S265 CORRECTION (MISTAKE-139, mac-mini cont.50 + opus-S243):** "≤6 distinct lifts" was RETRACTED
as a scale artifact (no scale has both small `B` AND few distinct lifts for generic DC). The
corrected effective-count invariant is **≤6 speeds coprime to 30030 = 2·3·5·7·11·13** — but this is
**bounded-diameter ONLY**: a single speed `= lcm(2..14) = 360360` witnesses every `d∈{2..14}`, so a
primitive DC 13-set can have **12** coprime-to-30030 speeds (klein-S265 verified:
`{1,17,19,…,59,360360}`, primitive DC, 12 coprime, still LOOSE `M = 23/112`). So the large-diameter
half genuinely needs the multi-scale reduction (THM-687/688 / THM-721's atom), not the effective
count alone. The two routes (THM-721 decorrelation atom; the pair-sum Parseval floor) still converge
on "spread DC collapses to a small effective family, trivially loose," and THM-721's elementary atom
(constant `1/13`) is the lighter tool; the Parseval floor adds the THM-680 sharpening and the
per-family certificate.

## SIMPLIFICATIONS (klein-S265) — the case split is really TWO cases, and the target is a CONSTANT

A combinatorial-simplification pass (4 Explore agents) found the proof structure carries more
apparent cases and heavier machinery than it needs. Four genuine simplifications, honesty-graded:

1. **The 5-way case split collapses to 2 (RIGOROUS).** The five buckets [non-covering / covering-DC
   / bounded-diameter / large-diameter / AP-wall] are not intrinsic. The **AP-wall is subsumed into
   non-covering**: the tight locus (`M=1/14`) is exactly `{AP, GW, dilates}` (THM-612 + klein-S206,
   exhaustive `n=4..7`), and it is **entirely non-covering** (`{1..13}` has no multiple of 14) —
   `#(tight ∩ covering) = 0`, proved. So every covering/DC family is in a strict cushion `M > 1/14`
   (band-edge, opus-S235). The intrinsic split is just **(1) non-covering** [`t = 1/q` sieve, and it
   already contains ALL tight families] **+ (2) divisor-complete** [strict `> 1/14`]. The
   bounded/large-diameter split inside (2) is a TECHNIQUE seam (finite check vs multi-scale), not an
   intrinsic case — unified in principle by the scale tower (THM-687/688) with the finite check as
   its base.

2. **Both routes A/B are ONE Bonferroni moment-ladder over pair-sum rulers (RIGOROUS as a
   reduction).** `bandCount(v,q,p)` (Route B) is the same empty-count `N` as the seven-sector residue
   (Route A), and `B5 = Σ_d(−1)^d S_d` is the alternating factorial-moment majorant both build
   (THM-671 `B5 ≤ liveCount`, THM-668 the maximizer lives on a pair-sum ruler, THM-707 clean pair-sum
   certifies at depth 5). Both reduce to ONE AP inverse theorem (the `{kα}`/AP coverage-extremality).
   *Caveat:* the "coverage-clearing DUALITY" (cont.47) as a QUANTITATIVE identity is HEURISTIC
   (corr +0.398, its own words); the rigorous statement is the moment-ladder identity + the shared
   (open) inverse theorem. At most one of the two routes' bookkeeping needs carrying to Lean.

3. **The target is a CONSTANT `1/13`, not a growing bound (RIGOROUS correction).** THM-720's "min M
   grows with diameter" is a generator artifact; the near-dilate `{L,…,12L,13L+1}` gives constant
   `M = 1/13` at every diameter (THM-721). Any machinery engineered for a *growing* or *tight*
   large-diameter bound (the density floor's tail, klein-S264's signed-`OffLine`) is heavier than the
   problem needs: only the constant margin `1/13 − o(1) > 1/14` is required, and THM-721's elementary
   atom delivers it for the `j ≤ 6` compressed stratum. Corrected inline in the S264 section above.

4. **The "≤6 effective speeds" shrink is bounded-diameter ONLY (verified).** `≤6 coprime to 30030`
   holds for 100% of bounded-diameter DC (opus-S243) but FAILS at large diameter: klein-S265
   exhibits primitive DC `{1,17,19,23,29,31,37,41,43,47,53,59, 360360}` with **12** coprime-to-30030
   speeds (a single `lcm(2..14)=360360` is DC by itself), still LOOSE (`M = 23/112`). So the shrink
   is a bounded-diameter tool; large-diameter is carried by the multi-scale atom (THM-721), not the
   count. Also: `auto-safe` (opus-S241) is a DISCRETE-clearing statement (bandCount at `p/q`), NOT a
   reach reduction — "drop structured speeds ⟹ reach ≈ coprime sub-family's reach" is NOT licensed by
   it; the bounded-diameter pigeonhole route (opus-S242) uses it correctly, the large-diameter
   decorrelation route must use THM-721.

**Combinatorial forward-lead (the one unexecuted move worth trying):** the SHALLOW sub-lemma ("no
dilate of a spread DC 13-set puts 6 speeds in a `1/7`-arc") is a COARSE-scale (seven-sector) inverse-
sumset statement, where `E₂`/Freiman is the correct invariant (HYP-5990), NOT the fine `1/14` scale
(governed by `E₃`/Schur). The parked **BSG → Freiman `3k−4`** bridge (opus-S181) fits it exactly:
6-in-an-arc ⟹ small-doubling block ⟹ (BSG) large low-doubling subset ⟹ (Freiman `3k−4`) short AP ⟹
contradiction with spread (longest-AP ≤ 7). The AP corner is already closed by three-gap (opus-S236),
so BSG→Freiman need only cover the dissociated bulk — its natural domain. (klein-S265 mining; files:
`lrc14_coprime30030_scope_klein_S265.py`(+out); reflection `the-combinatorial-simplifications-*`.)
