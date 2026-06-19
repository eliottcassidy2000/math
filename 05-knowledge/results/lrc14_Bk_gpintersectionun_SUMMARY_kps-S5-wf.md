# LRC(14) gp-intersection-uniform angle — adversarial verification (kps-S5-wf)

All numbers EXACT (fractions.Fraction). mu_exact cross-validated vs 20000-pt grid (40 shapes, both
thresholds, max disagreement 0; gap=theta breakpoints complete).

## Claim (0) m_P  [CONFIRMED EXACT]
min meas(G_P) per |P|=0..10: 1, 6/7, 66/91, 55/91, 1979/4004, 2243/5880, 3029/10780,
45107/229320, 2479/17640, 10601/114660, 14249/252252.  ALL MATCH.
m_P = 14249/252252 at |P|=10, P={1,2,3,5,7,8,9,11,12,13}.  MATCH.

## Claim (1) via-max refutation  [PRIMARY WITNESS CONFIRMED EXACT]
k=7 P={1,2,3,6,12,13} E={0,2,3,4,5,6,8}: meas(G_P)=515/1092, mu_2/7=13/35, rho*_2/7=0. ALL MATCH.
(My guessed E for the k=8,10 secondary witnesses gave rho*!=0; the prompt only specified P there,
not E. The k=9 guess hit 0. So the zero-locus is E-specific; the headline k=7 witness is exact.)
"All reconstructed S lonely" realizability NOT independently checked here (prompt flags consistency-only).

## Claim (2) k<=7 pigeonhole mu_1/7=1  [CONFIRMED (a.e.)]
Exhaustive bounded spread: k=3..7, min mu_1/7 = 1 over 36/120/330/792/1716 shapes. MATCH.
(k=7 maxgap<=1/7 only on measure-zero shifted-grid set — standard.)

## Claim (3) k>=8 consecutive mu_1/7 & union floor  [CONSTANTS EXACT; minimizer VERIFIED not PROVED]
Consecutive mu_1/7: k=8 691/735, k=9 247/294, k=10 38/49, k=11 1381/2205, k=12 13823/24255,
k=13 477/1078.  ALL MATCH.
thr_k = 1 - min_{|P|=13-k} meas(G_P): 3637/5880, 2025/4004, 36/91, 25/91, 1/7, 0.  ALL MATCH.
consec mu_1/7 >= thr_k at every k.  Union floor min[meas(G_P)+mu_1/7(consec)-1] = 1891/5880
at k=8 P={1,5,7,8,9}.  MATCH.
HUNT (consecutive minimizes mu_1/7, min>=thr):
  k=8 exhaustive spread<=12 (792): min=691/735=consec.
  k=9 exhaustive spread<=13 (1287): min=247/294=consec.
  k=10 (1698 struct+rand, spread<=2k): min=38/49=consec.
  k=11 (2072): min=1381/2205=consec.
  k=12 (2636): min=13823/24255=consec.
  k=13 (3197): min=477/1078=consec.
  Targeted near-grid/mod-7 + exhaustive subsets of {0..k+2}, k=8..13: min=consec, 0 violations.
  TOTAL mu_1/7<thr violations found: 0.
Status: VERIFIED computationally; "consecutive minimizes mu_1/7 for ALL integer E" remains a
REDUCTION/conjecture, not a theorem (no proof for unbounded shapes). Slack >=0.32 everywhere.

## Claim (4) quasi-independence R_1/7  [CONFIRMED EXACT]
min R_1/7 over consecutive E = 67053/84241 at k=11 P={1,12}. MATCH (exact).
Per-k: k=8 30570/33859, k=9 726/893, k=10 203/247, k=11 67053/84241, k=12 36253/41469, k=13 1.

## Surviving gap (precise)
1. The k>=8 branch is a REDUCTION conditional on the WEAK 1/7-spread bound mu_1/7(E)>=thr_k for ALL
   integer co-offset sets E (k=8..13). VERIFIED (consecutive minimizes, >=0.32 slack, exhaustive
   k<=9 + bounded/structured/near-grid k>=10) but NOT PROVEN for unbounded/all integer shapes.
   Strictly easier than via-max B(k); favoured by "large spread raises mu" (confirmed: mu rises
   toward 1 as spread grows).
2. Whole angle rests on the UPSTREAM global-witness implication rho*_1/7(P,E)>0 => M(S)>=1/14
   (THM-527/HYP-2581e), assumed here; finite-Vmax O(1/Vmax) discrepancy uniformity NOT re-derived.
3. Realizability "all reconstructed S lonely" is consistency (sampled Vmax), not a proof.
