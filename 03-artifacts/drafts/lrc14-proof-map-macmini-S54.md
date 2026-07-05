# THE LRC(14) PROOF MAP — every leg, its Lean name, its status

**mac-mini-2026-07-05-S54 (HYP-4117).** The consolidated obligation tree as of
2026-07-05 ~15:30, after opus-S81's MISTAKE-105 correction and today's census.
Status codes: **GREEN** = kernel-pure Lean in the corpus ([propext,
Classical.choice, Quot.sound]); **COMP** = exact computation, zero violations,
not yet a Lean object; **OPEN** = mathematics not yet discharged either way.
Owner policy: LRC(<=13) enters as the named citation hypothesis, not a sorry.

---

## The top surface (GREEN, two equivalent forms)

    lrc14_of_dichotomy_and_corner  (cite) (hdich) (hcorner) : LRC14Statement
        [klein-S133, LRCHcompSurface.lean]
    lrc14_of_spread_dichotomy_and_corner : ... (dichotomy needed on SPREAD
        bases only, ratio > 23/2)   [kps-S2, LRCHarmonicGate.lean]

with `cite : LRCUpTo13` (citation), and the two remaining named predicates:

    hdich   : TightLooseDichotomyAt (2/25)   -- every primitive compressed
              covering family's argmax-peel base is c*{1..12} up to sign
              OR has margin 2/25 at a point
    hcorner : CornerLonelyAt (2/25)          -- loose base + killer <= 25B/3

Glue all GREEN: hcomp_of_primitive (scale/sieve split), argmax peel,
tight_free_rider' (tight side consumes subset-of-c*AP-up-to-sign, no
permutation glue), gcd_killer_of_primitive.

---

## Branch 1: hdich = the n=12 spectral-gap dichotomy

The tight side is DONE (free-rider, GREEN).  What remains is exactly:
**no primitive 12-set has M in (1/13, 2/25) unless it is {1..12}** (then
c-scaled by the peel's gcd).  Every constraint below on a hypothetical
violator is GREEN unless marked:

### 1a. The formal violator profile (kps S2–S7 + klein + opus, all GREEN)
  - ratio gate: spread (w_max > 11.5 w_min)          [tightLooseDichotomyAt_of_spread]
  - covering every q <= 12 in every direction         [not_loose_dvd]
  - pinned 0/±1 mod every q <= 25, every direction    [not_loose_near_unit]
  - pair-covering mod 13 (no 13-mult case)            [pair_pinning_13]
  - binding pair >= 38 (=> w_max >= 19)               [gap_forces_big_pair, ω-parity]
  - 24-compression: w_max <= 24 w_2nd                 [peel_height_bound / gap_compressed_24]
  - order-statistic ladder on the top seven:
      w_(l) <= C_l w_(l+1), C_1 = 24, C_2..C_6 ~ 64.7/69.2/85.7/133/573
                                                      [gap_ladder_rung, LRCGapLadder]
  - certificate completeness: loose <=> bounded integer search (s <= B/(2Δ)+1)
                                                      [cert_of_margin, LRCCertCompleteness]
  - exact-M grid attainment (THM-592 in Lean)         [maximizer_on_grid, LRCGridAttainment
                                                       + klein-S137 Stage A in flight]

### 1b. The lift neighborhood (full residue system mod 13)
Residue pinning (GREEN, residue_pinning_13/residue_injOn_13): tight-from-above
no-13-mult families are lifts {r + 13 k_r}.  Strata:
  - l=1: GREEN — all 144 rows kernel-checked [LRCLiftRigidityRows]; floor
    14/169 at {1..11,168}; floor-level cutoff 2016 swept (COMP, S52).
  - l=2: COMP — full structural domain 600,756 sets (S52); floor 2/25
    ATTAINED at {1..12}\{4,6} ∪ {17,19}; the ≥2/25 witness is GREEN
    [block46_margin, LRCLiftFloorRows] + six ladder rows GREEN.
  - l=3: COMP — full chain domain, 13.3B cells, two independent
    implementations with exact level-count agreement (S53).
  - l=4: COMP-in-flight — chain domain, 6-way parallel, zero unresolved so
    far (S53 workers; post-pass committed).
  - l=5/6: COMP on boxes (k<=4 / k<=3) + 386 corners; chain TAILS **OPEN**
    (exact shapes + sharpened anchors in the S53 ledger; one C session each).
  - l>=7: legs 1–2 GREEN [LRCLiftPigeonhole]; leg 3 rebuilt: opus-S81 GAP
    DESCENT (spread tops dodged at any count) + the BOUNDED RESIDUAL SWEPT
    CLEAN (S55: 46 patterns, 7,071,570 sets, all witness-certified >= 2/25);
    remaining: the high ratio-cluster leg (opus's).

### 1c. The non-lift shapes (doubled ± pairs, 13-multiples)
Outside every lift sweep; constrained by 1a only.  TODAY (S54): the census
[25,40] — 129M leaves, 20.2M cheap-filtered, 58,178 full-profile survivors,
**every one witness-cleared >= 2/25 (zero hard cases)**, including all
32,764 doubled-pair shapes and 25,321 with-13-mult shapes.  With kps-S3's
[1,24] census: **the spectral gap is empty to height 40** (B=48 extension
running).  Status: COMP per height window.

### 1d. The remaining crux (OPEN — now NAMED and FINITE-SHAPED, S55)
**No residue filter can anchor height — proven** (S55 pole-necessity corollary:
the profile is periodic on CRT-frozen rays; the floating 7-cluster passes at
every frozen scale).  The mechanism must be witness-side, and the data says it
is: **THE Q50 CONJECTURE** — every profile-passing family has a margin-2/25
witness at modulus q <= 50.  Evidence: 511,947/511,947 census survivors
(first-q <= 45); 100% at composite-CRT moduli; 12,095/12,095 CRT-lifted
families at heights to 1.16e14 (first-q <= 35).  THE REDUCTION:
    gap emptiness = Q50 (finite template verification, single scale)
                  + cluster descent (opus-S81, between scales).
Remaining work: (i) the per-q template check at q = 26..33 (enumerate mod-q
templates consistent with pinned divisors, verify clearing-dilation
existence — bounded, engineering-shaped); (ii) opus-S81's high ratio-cluster
leg (the bounded residual is swept CLEAN: 7.07M sets, S55); (iii) compose in
klein's assembly language.

## Branch 2: hcorner = CornerLonelyAt 2/25

  - Geometry GREEN: IntervalEscape (one-tooth containment), THM-619 bands,
    two-pin lemma (prose+COMP), klein's window stack [LRCTeethR, tower_step_12].
  - COMP: 1,568 structured loose bases swept through bands+pins+window —
    1,566 EMPTY, survivors cleared exactly (S48/S49, THM-620).
  - **OPEN**: (i) the loose-base space enumeration argument (composition note:
    far base elements peel first, re-entering the pipeline at smaller scale)
    as a theorem; (ii) the band-certificate decidable tables in Lean
    (schema named by klein-S134; certificates exist in the .out files).

## What "finished" requires, minimally

  1. The 1d crux: one absolute-height mechanism (or per-scale census to a
     PROVEN cutoff).  This is THE mathematical gap.
  2. l=5/6 tails + l>=7 descent residual (bounded, sweep-shaped, in motion).
  3. l=2/3/4 strata as named decidable hypotheses or kernel tables in Lean
     (pattern exists: LRCLiftRigidityRows did it for 144 rows; the 10^10-cell
     strata need the pyramid-witness-stream packaging instead).
  4. hcorner's two OPEN items (composition + tables).
  5. Assembly: plug 1–4 into lrc14_of_spread_dichotomy_and_corner.

Nothing on this list is analytic; items 2–5 are engineering-shaped.  Item 1
is the only place where new mathematics might still be required — and the
census keeps saying the gap is empty everywhere anyone looks.

## The three recurring walls (why the constants are what they are)

  - 2ρl < 1 (fee-mean ceiling, MISTAKE-105) = l <= 6 at ρ = 1/13.
  - the 13u orbit window [u+1, 12-u] empties at u = 6 (pyramid, S53).
  - >= 7 lifted coords must hit a unique-multiple position (pigeonhole, S78).
  All three: "13 is just big enough."  The same largeness that makes the
  rigidity TRUE (q=5 contrast, opus-S74) makes it CHECKABLE.
