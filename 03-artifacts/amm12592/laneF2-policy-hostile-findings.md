# Lane F2 (G4) findings: policy-independence of the band freeze — REVISED

Session: death-star-2026-07-30-coinC2, lane G4 (repo lane F wave 2).
Script: `04-computation/amm12592_laneF2_policy_hostile_deathstar.py`
(modes: witness / beam / prefixgreedy / symdfs / asymdfs / selftest; selftest
PASS: scaled integer S == 2 D b^Lmax == descended laneD ledger value at all 9
biases, exact). Frame: THM-2966 normal form; depth law
d_m = floor(g1*m/g2) + D0, default gamma = 1/2, D0 = 0;
PROVED necessary envelope (E_m): |D_m(p)| <= (p^{m+1}+q^{m+1})/2 at every
row-end m, checked exactly at rational biases. Search bias set "9-bias" =
{1/2, 1/3, 2/3, 2/5, 3/5, 1/4, 3/4, 1285/2181, 8847357/11821757} (the two
certificate biases included). All asserted numbers are exact integer/Fraction
computations; floats only for display/scores.

STATUS SUMMARY (details and mechanisms below):

1. The lane-D "frozen-residual law" is a POLICY artifact in its original
   form, and lane C2's "M=11 wall" was a SEARCH-ORDER artifact: anticipatory
   policies blow far past both. (EVIDENCE, with exact re-verified witnesses.)
2. The surviving corridor's mechanism is (near-)ANTISYMMETRY
   (complementary-mirror schemes, delta_side1 = -delta_side0), which kills
   D(1/2) identically — the p=1/2 pinch that killed greedy — plus exact
   D=0 resets at the parity-free dyadic checkpoints M = 2^r - 1
   (THM-2976 T1 clock, confirmed in vivo). (EVIDENCE + PROVED parity frame.)
3. The freeze does not die; it MOVES: with only the 9 biases controlled, the
   stranded band value re-emerges as a frozen residual at p near 1/2
   (grid-death at m ~ 21 at p = 15/32); augmenting the controlled set with
   near-1/2 biases postpones it. The live policy-independent question is the
   PRECISION-vs-SURVIVAL law (measured below). (EVIDENCE.)

## 1. Witness anatomy at (gamma, D0) = (1/2, 0), M = 10 (witness mode)

Deterministic DFS (abs-first ordering + sound relaxation prune, 1562 nodes)
reproduces the lane-C2 witness; re-verified exactly (box + parity + E_m at 9
biases, all 10 rows). Anatomy:

- rho_m := max over 9 biases of |D_m(p)|/env_m(p) climbs to the wall:
  0, 0.500, 0.668, 0.894, 0.914, 0.957, 0.988, 0.9990, 0.9984, 0.9977:
  this witness ENVELOPE-SURFS: it survives at margin ~1e-3, with deep cells
  saturated (mean |delta|/cap = 0.83-1.0 per row).
- Binding-bias march: 1/2 (rows 1-2) -> 3/5 (3-5) -> 1285/2181 = cert bias A
  (6-7) -> 1/3 (8-10). The certificate bias is *transiently binding* exactly
  in the handoff rows.
- rho(1/2) locks onto 1/2 EXACTLY for rows 2-5 (D_m(1/2) = -env/2 =
  -2^{-m-2}): the survivor keeps a half-budget geometric residual at the
  symmetric point instead of zeroing it.
- Ledger split at p = 1/2 (laneD Ledger replay): the witness strands
  POSITIVE band value (+1.5e-3 at m=10, band coeffs {7:-120, 8:+84, 9:+228})
  cancelled in VALUE by negative cone mass (-1.2e-3); greedy strands
  same-sign band value {7:-106, 8:-138, 9:-106} it can never cancel. Forced
  band parity (beta_M restriction) PASS for both — parity is
  choice-independent (PROVED), the SIGNS are the whole game.
- Witness diverges from laneD zeroing-greedy at row 3 = m* - 1:
  anticipation must begin BEFORE band birth (m* = 4).
- laneD zeroing-greedy violates (E_m) first at m = 8 (its frozen residual
  3.33e-3 crosses the halving envelope 2^{-m-1}); g3's per-cell
  N(1/2)-zeroing greedy violates from m = 4. Policy ladder so far:
  zeroing greedies die at m*..2m*; the DFS witness reaches 10 at margin 1e-3.

## 2. Prefix-greedy continuation: the delay-not-reroute law (prefixgreedy)

Replay witness rows 1..P, then laneD zeroing greedy; first (E_m) violation
row V(P), exact, 9 biases, continuation to M = 30:

    P    : 0   2   3   4   5   6   7   8   9   10
    V(P) : 8   8   6   6   7   8   8   9   10  11

- For P >= 7: V(P) = P + 1 — death IMMEDIATELY after the anticipated prefix
  ends, with residual frozen thereafter (D(1/2) identical at M=20 and M=30).
- For P in {3,4,5}: V(P) < V(0) — an anticipatory prefix is WORSE for a
  zeroing continuation than greedy's own prefix: the witness's rows are
  tuned for future anticipation, not for value-zeroing.
- Lane G3's exhaustive extend result sharpens this beyond any policy: NO
  continuation of the M=10 witness survives M=11 (INFEASIBLE for clamped
  prefixes back to row 7): surviving TO M does not put you on a corridor
  THROUGH M. Anticipation with horizon H buys survival ~ H and nothing more.

## 3. Beam policy family: the M=11 wall falls (beam mode)

Beam search over per-cell transport choices, PROVED envelope as pruning
oracle, score = projected worst-bias rho (max-scoring; 'sum' = LP-guided
variant), width W = anticipation capacity; W=1 is lookahead-0 rho-greedy;
dedup + rand-keep for diversity. Every completed row-M beam state is an
exact witness for the (E_1..E_M) system at the searched biases
(independently re-verified: box + parity + envelope).

- W=1 (lookahead-0, max-score): survives to row 14, dies in row 15; tail
  rho grows x1.5/row while binding at p=2/3 (frozen residual at 2/3, env
  ratio 2/3 per row: 0.21 -> 0.31 -> 0.48 -> 0.74 -> dead). Even
  lookahead-0 with the RIGHT SCORE beats every zeroing policy (14 vs 8).
- W=256: survives ALL 20 rows (M=20 >> m*=4), best rho at 9 biases 0.02,
  nowhere near pinch. RE-VERIFIED EXACTLY at M=13 and M=20. So:
  * the lane-C2/G3 node-cap at M=11 was an artifact of abs-first DFS
    ordering (it explores the pinched near-symmetric corridor first);
  * M = 11..23 are FEASIBLE for the 9-bias envelope system: the M=20
    (9-bias) and M=23 (13-bias superset, W=256 augmented run) survivors
    are independently re-verified, and a witness's PREFIX is a witness
    for every smaller M — the first exact witnesses past the former
    frontier, found by policy search, not by exhaustive search;
  * the M=10 "pinch" (rho 0.998, pinned saturated cells) was an artifact
    of that same ordering: the abs-first DFS finds the WORST witness.
    (The in-class restricted DFS modes, symdfs/asymdfs, inherit the same
    ordering pathology and node-cap at M=13 — beam search, not DFS, is
    the right witness extractor in this landscape.)
- Mechanism: the W=256/M=13 witness is mirror-ANTISYMMETRIC on 43/55
  (m,k)-cells (choice list dumped and re-verified; the M=20 survivor is
  antisymmetric on 77/120 as band-vs-cone value play takes over):
  delta_side1 = -delta_side0, i.e. v_{m,k} = binom(d_m,k) - w_{m,k}
  (complementary-mirror labeling). Then D(p) = -D(1-p), so D(1/2) == 0
  IDENTICALLY — the strongest constraint (smallest envelope) is killed
  structurally, not fought. This is the same antisymmetry as the classical
  ratio-2 middle-pair trick (THM-2160), rediscovered by the beam at
  gamma = 1/2 < 1. Formally: complementary-mirror schemes need exactly
  f(p) - f(1-p) = (p - 1/2) with f(p) = sum p^m q W_m(p) — a ONE-SIDED
  antisymmetric flow problem (new reduction worth its own lane).
- Dyadic checkpoint exploitation (question (b)): beam bestrho hits 0
  EXACTLY at rows 3 and 7 (M = 2^r - 1, beta_M == 0: parity-free levels,
  THM-2976 T1 verified in-script for both depth laws) and near-0 (0.002)
  at row 15: the corridor RESETS at the binary-clock checkpoints. The
  (D3) clock is not construction-dead after all — it is exactly where the
  surviving corridor re-zeroes.

## 4. The freeze MOVES to p near 1/2: dense-grid audit (dense mode)

The 9-bias W=256 survivor at M=20, audited exactly on the dense grid
p = j/64 (all rows): ZERO violations through m=14, then a frozen residual
at p = 15/32, 17/32 (nearest grid points to 1/2) with rho DOUBLING per row:
m=15: 0.03 -> 0.06 -> 0.12 -> 0.21 -> 0.41 -> 0.80 (m=20) — extrapolated
grid-death at m ~ 21. Mechanism: antisymmetry kills D(1/2) but the band
(directions theta in (1/3, 2/3) at gamma=1/2) strands value in the ODD GERM
of D at 1/2 — D(1/2+eps) ~ eps * D'(1/2) — and the envelope at fixed eps
halves per row. The band IS the germ of D at the symmetric point: lane D's
freeze, restated at the accumulation point of the uncontrolled biases.

Augmenting the search set with near-1/2 biases {15/32, 17/32, 31/64, 7/16}
(13 biases, W=192, M=20): survives all 20 rows, dense j/128 audit ZERO
violations, worst grid rho 0.355, and D'(1/2) now DECAYS (6e-5 at m=10 ->
3e-7 at m=20; normalized x2^{m+1}/(m+1): 0.01-0.03, not frozen). The policy
controls whatever it is told to control; the stranded content moves to the
next uncontrolled scale. But best search rho then grows GEOMETRICALLY
(x1.5-2.0 per row) from an onset row, and this growth rate is
WIDTH-INDEPENDENT (W=64/192/256 tails nearly identical per-row values at
gamma=1/2: rho(18) = 0.078/0.081; rho(20) = 0.282/0.294) — structural, not
a search artifact. Controlling near-1/2 biases bought only a constant
factor (~2.7x =~ 1.4 rows) over the 9-bias run's dense-grid trajectory at
p = 15/32: augmentation does not change the growth RATE.

## 5. The death law is the A-CLOCK (anti-checkpoints A_M = 2^r - 1)

New clock, dual to THM-2976 T1's M-clock, machine-checked in-script
(`clock` mode) for all five depth laws below: at levels with
A_M = M + d_M + 1 = 2^r - 1, EVERY band position is forced ODD — because
binom(2^r - 1, o) is odd for all o (Lucas, all-ones), and the correction
(1 + x^{M+1})(1+x)^{d_M} in beta_M is supported on o <= d_M and o >= M+1,
disjoint from the band [d_M + 2, M - 1]. (PROVED; this sharpens lane D
sec. 2's parenthetical into the driving clock.) So the A-clock levels
plant |band| forced-odd UNTOUCHABLE coefficients all across the band.
Anti-checkpoints m_ac (band width w in parens):

    gamma=1/2 D0=0: 4(0), 20(8), 84(40)
    gamma=1/3 D0=0: 2, 5(2), 11(6), 23(14), 47(30)
    gamma=1/4 D0=0: 2, 5(2), 24(16), 50(36)
    gamma=3/5 D0=0: 4(0), 9(2), 19(6), 39(14)
    gamma=1/2 D0=2: 3(0), 8(0), 19(6), 40(16)

Measured death rows (beam, augmented 13-bias set, PROVED envelope oracle):

    gamma=1/2 D0=0: W=64 dies in row 22; W=256 survives 23 (rho 0.26,
                    x1.9/row tail -> projected death ~25)  [m_ac = 20]
    gamma=1/3 D0=0: W=64 dies in 23; W=128 dies in 25; W=512 dies in 26
                    with rho(25) = 0.986 — pinned to the wall [m_ac = 23]
    gamma=1/4 D0=0: W=64 dies in row 28 (M=30 rerun; rho tail 0.116 at
                    24 -> 0.203 -> 0.384 -> 0.697 -> dead)   [m_ac = 24]
    gamma=3/5 D0=0: W=64 dies in 21                          [m_ac = 19]
    gamma=1/2 D0=2: laneD zeroing greedy first violates at m = 16
                    (= 2 x m_ac,1 = 2x8, mirroring D0=0 greedy's 8 = 2x4);
                    beam W=48 dies in row 20 (rho(19) = 0.68, x1.7/row
                    tail from 15; binding 2/5 and 1285/2181) [m_ac = 19]
                    — D0=2 moved m_ac from 20 to 19 and the death row
                    followed it down (D0=0 W=64 died in 22).

- Death row = m_ac + (0..4) across all five depth laws, while band birth
  m* ranges 2.5..9 and does NOT correlate: (1/2, D0=2) opens its band
  LATEST (m* = 9) yet dies EARLIEST (row 20, m_ac = 19); gamma=1/3 opens
  at m* = 3 and dies at 23-26. The death law is the A-clock, not the
  band-birth law.
- The gamma-INDEPENDENT M-clock reading ("all die ~21 after the last
  bridgeable checkpoint M=15") is REFUTED: gamma=1/4 survives to 27 while
  gamma=3/5 dies at 21 at identical width and bias set.
- Every death lands 0-4 rows after the first MATURE anti-checkpoint (the
  first m_ac with m_ac >~ 19, where the envelope 2^{-m} is already tiny);
  earlier anti-checkpoints (bw <= 6 at small m) are survivable because
  cone budgets still dominate. No run at any gamma < 1 or any width came
  near the NEXT anti-checkpoint (84 / 47 / 50 / 39 all unreached).
- Width buys rows only logarithmically: at gamma=1/3, x8 width bought +3
  rows (23 -> 26); at gamma=1/2, x4 bought ~+2. The geometric rho-growth
  rate after onset is width-invariant.
- Historical alignment: laneC2's gamma=1/3 DFS node explosion at M=11
  (74M nodes) sits exactly on that law's m_ac = 11; greedy deaths at
  m = 8 (D0=0) and m = 16 (D0=2) are 2x the first anti-checkpoint.

## 6. Verdict

Lane D's frozen-residual law is NOT policy-independent AS STATED: "no
assignment survives past band birth m*" is FALSE (exact re-verified
witnesses to M = 20-26 >> m* at four rates; the M=10/11 'pinch' and wall
were DFS-ordering artifacts — abs-first ordering explores the doomed
near-symmetric corridor first, and the 'pinned saturated' M=10 witness it
finds is a dead-end corridor, not the corridor).

The defensible policy-independent restatement (sharpest we can currently
support, each clause with its mechanism):

  (i)   ZEROING policies (greedy family) die at 2 x m_ac,1 with lane D's
        PROVED mechanics (band strands same-sign value).  [PROVED + EVIDENCE]
  (ii)  ANY finite-bias-controlling policy survives far past m* by
        (near-)antisymmetry (kills D(1/2) structurally), exact D=0 resets
        at the parity-free M-clock checkpoints M = 2^r - 1, and band-vs-
        cone VALUE cancellation.  [EVIDENCE: exact witnesses]
  (iii) What actually kills every tested policy is the A-CLOCK: at
        A_M = 2^r - 1 the whole band is forced odd (PROVED); after the
        first mature anti-checkpoint the best achievable worst-bias rho
        grows geometrically at a width-independent rate, and death follows
        within 0-4 rows (+ log W).  [EVIDENCE, 5 depth laws]
  (iv)  Bridging to the NEXT anti-checkpoint (spacing x2 in m) would need
        the post-reset residual to be ~2^{-m_ac} — precision DOUBLING per
        epoch, i.e. a doubly-exponential precision wall in the number of
        epochs. No policy achieved even one bridge.  [EVIDENCE + mechanism]

Consequences for HYP-9061:

- C* = 2 is NOT in doubt from finite survival alone: survival to 20-26 is
  policy-achievable but epoch-bridging never was. The construction-gate
  reading of (27) revives only WEAKLY: any gamma < 1 construction must
  solve the antisymmetric one-sided flow f(p) - f(1-p) = p - 1/2 (new
  reduction, sec. 3) across infinitely many all-odd band levels.
- For the LOWER-bound direction this lane hands the dual a concrete
  choice-independent target: the forced-odd band vector at A_M = 2^r - 1
  (PROVED, choice-independent) is the natural per-epoch obstruction; a
  proof that its band value cannot be VALUE-cancelled by in-cone mass at
  band-direction biases (theta in (gamma/(1+gamma), 1/(1+gamma)), i.e.
  p in the band-dual interval around 1/2) at rate faster than 2^{-m}
  would convert the measured death law into C* = 2 — or, at the
  certificate rate, into the two-ray entropy comparison of (27). The
  binding biases observed at death (15/32, 17/32, 31/64, 7/16 — and
  1285/2181 at gamma=3/5's death row) are exactly band-dual biases.
- The A-clock epochs are dyadic in A (levels 2^r - 1), i.e. dyadic in
  T(m) = A_m: the epoch structure lives on the DEADLINE axis. Note
  d_m = 2^rho - 1 with A_m = 2^r - 1 is (D3)'s deep-corner timing; the
  full-band-odd statement here is its all-positions strengthening.

## 7. Reproduction

    python3 amm12592_laneF2_policy_hostile_deathstar.py selftest
    python3 amm12592_laneF2_policy_hostile_deathstar.py clock --g1 1 --g2 3
    python3 amm12592_laneF2_policy_hostile_deathstar.py witness --M 10
    python3 amm12592_laneF2_policy_hostile_deathstar.py prefixgreedy --M 10 \
        --prefixes 0 2 3 4 5 6 7 8 9 10
    python3 amm12592_laneF2_policy_hostile_deathstar.py beam --M 20 \
        --width 256 --dense 64
    python3 amm12592_laneF2_policy_hostile_deathstar.py beam --M 23 \
        --width 256 --extra "15/32,17/32,31/64,7/16"
    # + the gamma/D0 sweeps of sec. 5 (--g1/--g2/--D0/--width as tabled)

All witnesses re-verified by independent exact replay (box + parity +
envelope); beam completions ARE witnesses. Beam deaths are policy deaths,
NOT infeasibility proofs: exhaustive infeasibility past M = 10 remains
open at every rate (the one exhausted statement: no continuation of the
specific M=10 abs-first witness survives M=11 — lane G3).
