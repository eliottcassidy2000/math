# ANGLE H — Adversarial re-verification of the LRC(14) reduction chain (mac-mini-2026-06-18-S?)

Independent engine (`lrc14_adversarial_chain_macmini_0618sH.py`), exact rational.
All canon anchors reproduce; every link tested. **Verdict below.**

## Anchors reproduced EXACTLY (independent engine, no shared code with canon scripts)
- mu_{2/7} consec: 1, 19/21, 9/14, 4/7, **83/210** (k=3..7). [confirms THM-530's correction of THM-527-F: consec_7 = 83/210, NOT 13/35]
- mu_{1/7} consec: 1, 691/735, 247/294, 38/49, 1381/2205, 13823/24255, 477/1078 (k=7..13). ALL match THM-531-A.
- meas(G_P) min row psz=1..10: 6/7, 66/91, 55/91, 1979/4004, 2243/5880, 3029/10780, 45107/229320, 2479/17640, 10601/114660, **14249/252252**. m_P confirmed.
- meas_S7(consec_8) = 481/1470 = 0.32721; M7(8)=20160/823543=0.02448. Match THM-532-A.
- Scale-invariance of meas_S7: meas_S7(dE)=meas_S7(E) confirmed exact for d=2,3,5.

## LINK-BY-LINK VERDICT

| Link | Claim | Verdict |
|---|---|---|
| (1) slow-fast change of var (THM-527-A) | rho_K -> rho*, O(1/Vmax) error, Vmax<=V0 finite check | **HOLDS empirically but GLUE UNWRITTEN.** Good-period count >=22 in all sweeps (never 0). BUT the "Vmax<=V0 finite check" for k>=3 is **asserted, never enumerated in any script** — the only enumerated finite cores are the k=2 slice (<=62) and general covering samples. The O(1/Vmax) error is not made a rigorous LOWER bound anywhere. |
| (2) global-witness gap>1/7 => M>=1/14 (THM-530/531, "upstream/assumed") | exact, boundary | **HOLDS.** Half-width 1/14 teeth, gap>2*(1/14)=1/7 leaves a safe phi: arithmetic correct. End-to-end: every reconstructed S=P∪{Vmax-e} has exact M(S)>=1/14 (6 families x 24 Vmax each). rho*>0 in every case. No counterexample to sufficiency. |
| (3) N(E) subset S7(E) | maxgap<=1/7 => every sector hit, up to meas zero | **HOLDS, no boundary leak.** 0 violations over random exact x. Boundary x=1/7 (maxgap=1/7 exactly): all 7 sectors STILL hit (half-open [j/7,(j+1)/7) convention saves it; a point at j/7 lands in sector j). |
| (4) k<=7 pigeonhole mu_{1/7}=1 | k<=6 always, k=7 a.e. | **HOLDS.** mu_1/7=1 for all shapes k<=7 (consec, perforated, dissociated). k=7 tight locus {maxgap<=1/7} is measure-zero (0/200000 sampled) => rho*=meas(G_P), benign. |
| (covering used?) is meas(S7)<=cap_k / mu_1/7>=thr_k claimed for ALL integer E? | find breach | **CLAIMED FOR ALL INTEGER E, and NO BREACH FOUND.** k=8 (dangerous): exhaustive box B<=20 (>77k subsets via earlier run) + 8000 random large-spread: max meas_S7 = 0.32721 = consec, cap_8=0.38146, **0 breaches**. k=9..13 boxes: 0 breaches. So covering/primitivity is genuinely NOT needed for the S7/mu_1/7 crux — the bound is universal over integer E. (Covering IS used upstream, only to guarantee G_P nonempty / S is a covering 13-set.) |

## THE WEAKEST LINK (the actual gap, in priority order)

1. **THM-532-D crude product bound `corr <= C*W` DOES NOT CLOSE (confirmed exact).**
   corr(consec_8) = 0.3027 (true), but C*=0.0395, W(consec_8)=9.73 => C*W = 0.384 > margin_8 = cap_8 - M7(8) = 0.357. The worst ratio C* and max weight W are at DIFFERENT shapes, so the product overshoots. **The high-relation-height certificate is not a closed proof**; it needs the finite low-height residual, which is **described but NOT enumerated**. This is the headline honest gap and it is REAL.

2. **The finite-Vmax glue (THM-527-A gap iii) is the load-bearing unwritten piece.**
   "Vmax <= V0 is a finite check" — for the k>=3 S3 residual this finite set is **never enumerated or checked in any script**. The whole reduction is asymptotic (Vmax->inf); the passage to a witness at FINITE Vmax rests on (a) rho*>0 + (b) a uniform O(1/Vmax) bound + (c) a finite check, of which only (a) is established (and even (a) is the OPEN crux mu_1/7>=thr_k / HYP-2602). Empirically the good-period count is large (>=22), but that is verification, not proof.

3. **The crux itself, HYP-2602 (mu_1/7(E) >= thr_k for all integer E, equiv. consecutive minimizes mu_1/7 / meas_S7<=cap_k), remains CONJECTURAL** — verified exhaustively in boxes and 8000+ random shapes with ZERO failures, but no proof. THM-532 reduces it to (1)+a finite residual; neither is closed.

## NET
No NEW gap or off-by-one was found in the *stated* logic: links (2),(3),(4) and the union bound are sound; the anchors are all exactly correct; the S7/cap reduction is not breached and is genuinely universal over integer E (covering not needed for it). **The chain's weakness is exactly where the canon says it is — but slightly worse than billed:** the two "finite checks" (the low-height S7 residual in THM-532-D, and the small-Vmax check in THM-527-A) are **named but not actually enumerated in code anywhere**, and the one written certificate (corr<=C*W) provably overshoots at k=8. So the residual is NOT just "bookkeeping": it is an unwritten finite enumeration PLUS an unproven universal extremal inequality (HYP-2602). LRC(14) remains open; the reduction is honest about being open, and its load-bearing glue is unwritten rather than wrong.
