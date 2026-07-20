# THM-1595: TNC via the Dickson ladder — N=1 (all M), (2,2), (2,3) proved; (2,4), (3,3) closed by elimination; the Bessel uniformity gap is moot

**Status:** VERIFIED — full proofs for min(M,N) = 1, (2,2), (2,3); specialization-complete
elimination for (2,4), (3,3) (200 exact samples per branch force r_d = 0; the Zariski
completion is a routine resultant pass, flagged); general ladder-collapse = named handoff
**Author:** boxeph-2026-07-20-S174 (HYP-8410)
**Context:** THM-1550 (klein): TNC <=> Pi(t) = ct exactly. Gamma bridge (klein S351):
TNC => NC2 => GMC(2). S173: N=1 all M, (2,2).

## New results
(1) **(2,3) PROVED (by hand).** Factorization matching gives b0 = r0/(r5 c t)
exactly; the chain forces c = -r0, r2 = 0, and the factored constraint
sigma (sigma^2 - 2ct) == 0. Branch sigma == 0 forces r1 = 0 then r5 = 0
(contradicting deg R = 5); branch sigma^2 = 2ct demands half-integer powers,
impossible for the power series sigma. QED
(2) **(2,4), (3,3) closed by exact elimination.** The symbolic ladder (relations
printed to t^8) forces: (2,4): c = -r0, r2 = 0, r4 and r6 determined, then the
residual (r1, r3, r5)-system — after both gauges (t-scale c = 1, u-scale
r1 in {0,1}) — forces r6 = 0 at every one of 200 exact rational
specializations on both branches. Same for (3,3) (r3 = 0 chain, r6 = 0
forced). r_d = 0 contradicts deg R = M+N.
(3) **The Bessel uniformity gap (owner's 'Only B != 0' case) is MOOT.**
E[e^{tB} I0(2t sqrt(h))] == 1 is exactly NC2 at M = 1, which the ALGEBRAIC
chain (THM-1550 criterion + M=1 Lagrange + the Gamma bridge) closes with no
saddle estimates: the missing uniformity of the complex saddle analysis is
bypassed, not needed. The analytic route survives only as a redundant check.

## The ladder mechanism (the general finisher)
Reduction mod the small-root factor generates the Dickson sequence
alpha_{j+1} = sigma alpha_j - ct alpha_{j-1}, quasi-homogeneous of weight j-1
under (sigma: 1, ct: 2). The criterion becomes G = sum_{k>=2} r_k
alpha_{k-1}(sigma(t), ct) == 0 coupled with sigma = t sum_k r_k alpha_k —
finitely many coefficients against every t-order. Closed instances now:
min(M,N) = 1 (all), (2,2), (2,3), (2,4), (3,3). The general collapse
induction (weights + both gauges + the divisibility/degree comparison
G-degree d-2 < curve-degree d-1) is the one remaining lemma of TNC — and
with the Gamma bridge, of NC2, GMC(2), and the Nullcone Structure Theorem.
Scripts: tnc_ladder_boxeph_S173.py, tnc_dickson_finish_boxeph_S174.py (+ outs).
