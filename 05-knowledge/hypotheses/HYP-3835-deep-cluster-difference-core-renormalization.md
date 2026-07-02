# HYP-3835: the DEEP-CLUSTER residual renormalizes to a DIFFERENCE CORE one scale down; the AP is the fixed point of the difference flow; the tiling times are continuous-Fraenkel points

**Status:** OPEN (program) -- first predictions CONFIRMED (two-scale factorization exact-limit verified;
worst pattern = consecutive = AP difference core; D-zeros = exact tilings; 3.8x empirical margin)
**Instance:** opus-2026-07-01-S32
**Script:** `04-computation/lrc14_simultaneous_peel_lemma_opus_20260701.py`, part (P) (+ `.out`)

## The residual (from HYP-3834)

After the simultaneous peel, the r=2 residual is: 11-cores containing a gap-free cluster
F = {N + c_1, ..., N + c_j}, j >= 7, N -> infinity, ratios bounded, plus a compact part C_low with
<= 4 elements. The union bound provably dies at j = 7 (danger measures sum to 1: a 7-cluster CAN tile).

## The renormalization claim

As N -> infinity the fast variable u = Nt decouples from the slow offsets tau_i = c_i * t, and

    meas(L_{C_low u F})  -->  Int_{L_{C_low}} D_F(t) dt ,   D_F(t) = meas_u( n_i (A - c_i t) ),

i.e. the deep cluster acts as ONE effective fast runner carrying its DIFFERENCE PATTERN {c_i - c_1} as a
frozen offset configuration driven by the slow time. The cluster is replaced by its difference core one
scale down. Iterating: the scale tower TERMINATES (each step trades j deep elements for a (j-1)-element
difference core + 1 effective runner).

**The AP is the fixed point:** differences of a consecutive AP are again an AP. So the extremal deep
cluster should be AP-like at every scale -- the same reason the census extremizers (pentagon = AP minus
2, tight locus {AP, GW} of THM-523) are what they are. The census extremizers and the deep-cluster
extremizers are ONE family seen at heights h=1 and h->infinity; dilation-invariance meas(L_{gC}) =
meas(L_C) is the exact symmetry connecting them.

## Verified today (exact/grid, script part P)

1. **Convergence** {N..N+10}: meas -> limit integral 0.158059; diff +0.0035 (N=30) -> +0.00002 (N=120).
2. **Two-scale factorization** {1,2,3,4} u {N..N+6}: predicted Int_{L_low} D_7 = 0.1041; actual
   0.1054 (N=240). Factorization quantitatively right; all values > 1/36 (3.8x margin).
3. **Worst pattern = consecutive** (difference core = AP {1..6}), 0.1059, vs even-spread 0.123,
   dilated-AP-shadow 0.198, lacunary 0.192, random 0.215. Prediction CONFIRMED at j=7, N=210.
4. **The tiling zeros:** D_7(t) = 0 EXACTLY at t = k/7 -- at those slow times the seven danger arcs
   tile the circle (total measure exactly 7 * 1/7 = 1). An exact cover by interval-APs with distinct
   speeds = a CONTINUOUS FRAENKEL / disjoint-covering-system configuration (Fraenkel's conjecture
   territory -- NEW connection, not previously in the repo). The zeros are isolated and non-degenerate
   (tilings break under perturbation since arcs move at distinct speeds), so Int_{L_low} D stays
   positive; the whole question is the QUANTITATIVE floor of that integral.

## What would close the deep-cluster residual (the honest list)

(a) An effective convergence rate meas -> Int D with error O(1/N) uniform over bounded-ratio patterns
    (= HYP-3787's O(1/w) equidistribution, now in exactly the right coordinates);
(b) a positive floor for Int_{L_{C_low}} D_F over: compact C_low (<= 4 elements, census-able) x
    difference patterns c (a j-1 <= 10-element integer pattern -- INDUCTION on the tower);
(c) the observation that (b) is again an LRC-type average -- the residual of LRC(14) at depth 1 is a
    SHIFTED LRC with effectively 8 runners (7-cluster + observer): n drops 14 -> 8; the tower loses
    6 runners per level and terminates.

The 1.16x tightness of the census (pentagon = 1.161 x 1/36) vs the 3.8x deep-cluster margin says the
binding case is COMPACT (height 1), consistent with the fixed-point picture: the flow contracts toward
the AP at height 1, where the census already holds.

## Boundary falsification probe (same session): SURVIVED

The dangerous region is the boundary between census (1.16x tight) and deep clusters (3.8x loose).
Probed (`lrc14_deep_cluster_boundary_minimax_opus_20260701.py`, exact rationals, ~7500 configs):
(A) consecutive 7-cluster x ALL 715 compact 4-cores x heights 14..120: min 0.0662 (N=19);
(B) pentagon-difference drop-pattern clusters: min 0.0818;
(C) mixed compact-k + (11-k)-cluster, k=0..3: min 0.0635 at {1,3,4}+{20..27} (k=3, N=20);
(D) barely-gap-free two-scale clusters (ratios 2-5): min 0.1718.
GLOBAL boundary min = 0.063457 = 1.967 x pentagon = 2.284 x 1/36. NOTHING dips below the height-1
census extremizer; the worst direction is clusters DESCENDING toward compact height (N=20 binding,
large-N looser) -- exactly the contraction-toward-the-fixed-point picture. Non-monotone resonance dip
around N=19-23 (cluster commensurate with compact max 13) noted; C_low = {1,3,4,11} recurs as the
worst compact 4-core (a near-divisor-chain -- worth a look).

-> HYP-3834, MISTAKE-090, THM-522, THM-523, HYP-3787, OPEN-Q-108; Tao 1701.02048 (bounded-speed
   reduction -- this renormalization is its quantitative local form); Fraenkel disjoint covering
   systems (new lead, see INVESTIGATION-BACKLOG).
