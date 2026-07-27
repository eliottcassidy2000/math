# LRC(14) gap (1/14, 3/41) — adversarial engineering ledger (scratch)

Engine: exact breakpoint-moduli M (moduli {v_i+v_j} u {|v_i-v_j|} u {2v_i}),
integer arithmetic, Fraction comparisons. Referee: 7/7 canon gates
(1/14, 1/13, 2/27, 1/14, 3/41, 8/105, 14/183) + 12/12 random families vs
full q<=2*vmax sweep.

Constraint stack used (all canon-derived):
covering 2..13; v_max >= 65 (THM-1290+S320, no later extension found);
spread > 12; pair-sum admissibility/D=M*s (THM-1268/1269);
pinning W(q) <= ceil(3q/41)-1 at every q (>= 3/41 escape otherwise);
window [43,48] blocked by construction (S323 law in reverse); <= 91^12.

## Assembled candidates (reached exact-M or a fold certificate)

| family | exact M | witness | verdict class |
|---|---|---|---|
| {4,8,...,44,52,96} = 4*{1..11,13,24} | 1/14 | 1/56 | tight dilate (non-primitive) |
| {4,8,...,48,52,276} = 4*({1..5,7..13}+{69}) | 3/40 | 14/320 | dilate-lift kill q=4*80 |
| {4,8,...,48,52,300} | 7/86 | 15/344 | dilate-lift kill |
| {3,22,...,50,360} | 11/69 | 34/69 | pack-balance destroyed |
| {1..11,13,1944} (round-2 sole primitive 1-tail) | 162/1949 | 812/1949 | fold saturation at q=h+5 -> 1/12 |
| {1..12,65} | 5/66 | 5/66 | corner exact |
| {1..12,66} | 1/13 | 1/13 | corner exact |
| {1..11,13,65},{...,66},{...,67} | 1/12 | 1/12 | corner exact |
| single-far {1..11,13}+w, w in [56,3000] | min 1/13 (w=60) | 27/65 | fold saturation |
| single-far {1..12}+w | min 5/66 (w=65) | 5/66 | fold saturation |
| single-far {1..5,7..13}+w | min 3/40 (w=69) | 14/80 | fold saturation |
| 13-term APs (two-anchor star rays), all (a,d) d<=219 | min ~97/206 =.4709 | - | stack catastrophic |

NOTHING BELOW 3/40 = 0.075 EVER ASSEMBLED. Gap ceiling is 3/41 = 0.07317.

## Kill-stage statistics by round

R1a survivor one-swaps (~270k tried): 4 pass [14,48] stack; 0 in gap.
R1b two-scale smalls+pack P>=65 (122k): 0 pass stack; deaths at q=17..23 (d=1 band).
R2 12-core + CRT tail (732 cores, h<=2600): every q individually fixable for
   all 732, but joint CRT admits tails for only 3 cores (23 tails);
   1 primitive assembled -> fold-killed; rest non-primitive (g=2,4).
R3 11-core + two tails (4137 cores, h1<=1100, h2<=5200): FINAL 2,602,378
   single-residue passes, ZERO joint CRT pairs, zero assembled.
R4A load-2 doorway, 12-core+tail (65 cores, h<=8500): 73 admissible tails,
   ALL non-primitive (g=4: 38, g=2: 35) -> collapse to dead single-far cores.
R4B load-2, 11-core + two tails: S1: 37 tried, 0 assembled; S2: 1 tried,
   non-primitive g=4. Dilation self-closure confirmed.
R5 pack-at-top (twin-peak, k>=7-far, 10 menus x 14 patterns x N<=6500):
   0 of ~900k pass the [14,60] tau-CRT. The one-translation-per-modulus
   freedom cannot satisfy ~40 simultaneous unit-cover conditions.
R6 window-covered-by-inclusion (11 duty blocks at scale 78..192):
   0 exact covers with 4-7 fillers (pool [1,600]); greedy leaves 24-26 of
   188-284 unit-pairs uncovered; per-filler useful contribution ~10 moduli
   vs 26-32 active moduli needing ~2-3 contributions each.
R7 beam over irregular hybrids (swaps from duty pool [15,300]): depth-1
   complete (141,530 states): the ONLY zero-bad-[14,60] states reachable are
   the dilates 4*{1..11,13,24} (g=4) and 4*{1..13} (g=4, vmax 52) — i.e.
   the tight orbit itself; nothing primitive, height>=65, stack-clean exists
   within one swap of the entire escape atlas. Depth 2+ stopped (beam best
   already = tight dilates; every deeper state inherits a strictly positive
   bad-count or the same collapse).

## The four unavoidable witness classes (every death is one of these)

W1. d=1-band units (q in [17,27], esp. 19/23/25/26/27): uncovered q needs
    {±v} to cover ALL phi(q)/2 unit pairs; 13 slots minus duties cannot.
    Kills all low-structure and pack-translation designs (R1b, R5, R6).
W2. escape-window [43,48] depth-3 arcs (S323 absolute window): kills every
    height-<=55 pinning survivor (atlas); blocking it by inclusion costs
    6+ slots (R6) or CRT residues whose joint lattice has gap >= 2x10^3
    (load 2) to 10^16 (load >= 7) vs fold-safe windows of length <= ~10^3.
W3. fold moduli q = h ± v (and 2h, c*q' for dilates): the 12-core's Farey
    comb. W_V(q) = W_C(q) exactly (distance-class collapse), and
    W_C(q) >= q*M(C) - H/2 with M(C) >= 1/13 (12 classes, LRC(13) settled)
    -> hard ratio bound v_max <= 533*H/4 + H ~ 134*v_(2); empirically the
    comb dies far earlier (last safe q = 233 for {1..12} vs bound 1866) and
    max safe-run <= 10 < fold span 2H >= 24: comb-run pigeonhole.
W4. dilation collapse: designs forced into g*Z lattices reduce to primitive
    cores at height/g whose fold combs were already exhausted (single-far
    curves >= 3/40 for all w >= 56).

## The two closed mini-theorems (exact, referee-checked)

T-A. No completion of the cores {1..12} or {1..11,13} lies in the gap at ANY
   height: h <= 64 by THM-1290; h in {65,66(,67)} exact values 5/66, 1/13,
   1/12; h >= 67(68): the 24 fold moduli span an interval of length 2H+1 >=
   25 inside [55, q*], but SAFE(C) has maximal run 10, and SAFE(C) = empty
   for q > 233 (checked exhaustively to the saturation cap q* = 1866/640).

T-B. Ratio bound: any 13-speed family with M < 3/41 has
   v_max <= H/(2(M(C)-3/41)) where C = family minus v_max, H = max(C),
   M(C) >= 1/13. In particular v_max < 134*v_(2): no lonely top scale.
