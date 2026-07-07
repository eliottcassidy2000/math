---
source: opus-2026-07-07-S144
status: audits + exact closed forms + two refutations (per-q minimality in primitive
  normalization; the mu = M ladder collapse) + one surviving new law (deep-resonance-tail
  AP-maximality, T_5/T_6) + the V_r invariant profile. Concurrent with mac-mini-S54
  (HYP-5147) and klein-S173 (HYP-4971), reconciled below.
tags:
  - lonely-runner
  - LRC14
  - per-q-windows
  - density-floor
  - three-gap
  - motzkin
  - circular-coloring
  - validity-audit
---

# The invariant resonance profile, and the ladder separation

**opus-2026-07-07-S144 (HYP-5137).** Owner mandate: understand the LRC(14) history and
validity deeply, find what we've been missing, extend forgotten factoids, aim at how/why.
Three deliverables: (1) a validity audit of the per-q window program (HYP-5117) that
converged, mid-session, with two concurrent fleet corrections — plus the pieces they did
not have; (2) the **V_r invariant profile**, exact at the AP, with a surviving new
extremal law (deep-resonance-tail maximality); (3) the resolution of the S141/S143
**μ = M hunt**: the homomorphism-ladder collapse is **refuted** at |S| = 4 with
hand-checkable certificates, and the ladder is exactly rigid up to |S| = 3.

## 1. The per-q audit — what three agents found in one day, reconciled

kps-S72 (HYP-5117) proposed: μ_{1/7}(E) = Σ_{q=2}^{6} W_q(E), the AP minimizes each W_q,
each W_q is σ-odd residue data ⟹ the first σ-odd crossing of the σ-even floor. Three
audits landed within hours (none aware of the others until mid-session pulls):

- **mac-mini-S54 (HYP-5147):** derived the true local-window formula
  `W_q(E) = Σ_p Σ_± max_s (G_s/q − 1/7)_+ / D_s^±` and showed kps's measured inequalities
  were **attribution-saturated** (the nearest-rational cell fills with generic good mass at
  high μ); for the *true* windows the AP does **not** minimize (random families' huge drift
  rates give W_2 = 0.0067 ≪ 0.0649). Robust reformulation: the window-sum **lower bound**.
- **klein-S173 (HYP-4971):** the crude-cap version is *provable in a paragraph*:
  μ_{1/7}(E) ≥ 146/(35·span) for every 13-set — an elementary diameter floor discharging
  span ≤ 73, independently matching the roof-route ledger.
- **This session:** two pieces the others did not have, plus the exact AP-side constants:
  - **The q = 1 ledger.** The attribution is a Voronoi partition *including q = 1* (the
    origin resonance), which every statement of HYP-5117 omitted. Exactly:
    **W_1(AP_13) = 1/7** — 32.3% of μ(AP_13) — and the AP (minimal diameter) is the
    W_1-*maximizer* shape among its perturbations, not minimizer. Any per-q program must
    carry an origin-budget clause: at k = 13 the AP's origin windows are worth exactly
    one whole 1/7-threshold of measure. (Cute and true: the AP_13 origin window half-width
    (6/7)/(k−1) equals **1/14, the LRC threshold itself**, exactly at k = 13.)
  - **Primitive near-affine adversaries.** kps's escape hatch (ii) was "fix the primitive
    representative." That does not save the claim even for the *attributed* windows:
    spread-AP-with-one-bump families — primitive, gcd of differences 1, **not** affine
    images — inherit the dilated profile and undercut individual attributed windows while
    paying elsewhere: `{1,3,…,23,26}` has W_3 = 0.0782 < 0.0857 = W_3(AP_13);
    `{1,4,…,34,38}` has W_2 = 0.0494 < 0.0649. (μ-excess +0.11 in both.) The dilation
    action pushes window mass 2-adically along the q-axis (q=2 → q=1, q=4 → q=2, …), and
    that transport **survives perturbation into primitive families**. Per-q minimality is
    false in every normalization; only the sum is minimal.
  - **Exact AP windows (the roof's closed form).** On the Farey-k grid the roof gives, for
    each node p/q (q ≤ 6) with Farey-k neighbor denominators q_L, q_R:
    window mass = (7−q)/(7q) · [1/(q_L−q) + 1/(q_R−q)]. Summing:
    **W_q(AP_13) = 1/7, 5/77, 3/35, 8/147, 23/294, 4/245 (q = 1..6), summing to
    477/1078 exactly** — the canon constant, now decomposed. The identity behind it:
    the *cluster-gap closing rate* at node p/q **is** the Farey flank minus q
    (e.g. at 1/3 the three gaps close at rates {11,11,8} on the right and {10,10,10} on
    the left; binding rates 8 = 11−3 and 10 = 13−3 are exactly the flanking denominators
    of 1/3 in F_13 minus 3). mac-mini's drift-rate formula and THM-637's Farey
    combinatorics are the same object; this is the exact-rational bridge.
  - **Why lift families do not collapse windows** (the audit's first surprise): for
    `{1..12, 13+60m}` the far element's huge rate does *not* shrink the q ≤ 6 windows —
    a fast phase **transits** a big gap (splitting it briefly, both pieces often still
    > 1/7) rather than closing it. The window ends when the **slow bulk** spreads, not
    when the fast element moves. This is the far-peel/simultaneous-peel mechanism
    (HYP-3900) appearing *inside each resonance window* — and it is why the per-q window
    masses of lift families stay ≥ the AP's (verified m = 1, 5, 20; mod-420 variant too).

**Net state of HYP-5117 after the day:** the decomposition idea survives as klein's
provable lower-bound assembly (crude caps) + mac-mini's exact formula for the structured
side; the *extremal per-q lemma is dead in all normalizations* (three independent
refutation routes); the σ-odd "crossing" was an artifact of a non-invariant coordinate
(the q-label) — the honest invariant content stays σ-even, which reaffirms kps-S67's
grading against kps-S72's optimism.

## 2. The V_r invariant profile, and the law that survives

The dilation problem is not bookkeeping; it is the statement that *the q-label is not a
function of the configuration*. The label that **is**: r(x) = the number of circular gaps
of {frac(e·x)} exceeding 1/7. Define

> **V_r(E) = meas{x : exactly r gaps > 1/7}**, r = 1..6 (seven gaps > 1/7 cannot sum ≤ 1),
> so μ_{1/7}(E) = Σ_r V_r(E), and **every V_r is fully affine-invariant** (gap multisets
> are translation/dilation/reflection-invariant). Verified: AP_13, 2AP−1, 3AP−2 have
> identical profiles to 4 decimals.

Near a q ≤ 6 node the phase vector is *exactly* linear in δ, so the r-profile of each AP
window is sliced by the sorted closing rates, giving exact rationals:

> **V_r(AP_13) = 17/99, 41/539, 94/1155, 149/2695, 65/1386, 6/539** (r = 1..6), summing
> to 477/1078 exactly (numerics agree to 1e-6). Tails: **T_4 = 5497/48510,
> T_5 = 563/9702, T_6 = 6/539**, where T_j = Σ_{r≥j} V_r.

The empirical shape (part-2/part-3 runs, structured zoo + 200 random + multi-start
hill-climbs, k = 13 and k = 8, MISTAKE-119 discipline):

- **T_4-maximality FAILS**: mod-7-structured families beat the AP at k = 13
  (0.1245 > 0.1133); {1..7, far} beats at k = 8.
- **T_5- and T_6-maximality SURVIVE everything thrown at them**:
  > **The deep-resonance-tail law (new, open):** among k-element integer families, the AP
  > maximizes T_5 and T_6 — the measure of times with ≥ 5 simultaneous gaps > 1/7,
  > i.e. phases confined to < 2/7 of the circle in ≥ 5 clusters.
  (At k = 8, {1..7, far} ties T_5 within numeric resolution — the far limit collapses to
  AP_7 + a transiting phase, a flat direction, honest caveat.)
- The perturbation data reads as an **exchange rate**: near-affine bumps destroy resonant
  mass (V_6 → 0, ΔT_4 ≈ −0.02) but create bulk mass (V_1..V_3) at ~5–6× the rate, netting
  Δμ ≈ +0.11. **R2 in this frame = "bulk creation beats resonance destruction at rate
  > 1"** — the AP minimizes μ *because* it is the unique family whose good measure is
  maximally deep-resonant (my S134 "pure window, zero bulk," now quantified per level).

Why this is a better structured-side target than per-q W_q: (i) affine-invariant by
construction — the entire dilation pathology of §1 vanishes; (ii) the T_6 event is
*finite and arithmetic*: 13 phases in ≤ 6 arcs of total width < 1/7, pairwise separated
by > 1/7 — a 6-block partition/simultaneous-clustering event whose probability is
governed by the gcd structure of within-block differences (common divisors cluster
simultaneously on positive measure; incommensurate differences only on product-small
sets). That is the σ-odd mechanism kps wanted, now attached to an invariant statement.
(iii) The three-gap quantization (mac-mini-S15's frame, factoid C3) is the *why* of the
AP's dominance: three-gap forces the AP's window gaps to be *equal*, and equal gaps
maximize #{gaps > θ} at fixed sum — the AP tops the r-count pointwise inside its windows.
Schur-flavored, but on the gap vector at fixed x, not on the family step vector (which
kps-S63 rightly killed).

## 3. The ladder separation: μ = M is FALSE from |S| = 4, exactly

S141 built the homomorphism ladder (witness ⟹ circular coloring ⟹ independent set:
LRC(14) ⟹ GRAPH-14 (χ_c ≤ 14) ⟹ MOTZKIN-14 (μ ≥ 1/14)) and asked whether the rungs
collapse (χ_c = 1/M?, μ = M?); S143 proved |S| = 2 collapse and left the hunt running.
This session replaced the period-scan engine (N ≤ 320 cap, O(300 DPs/set)) with the right
object: **μ(S) is the maximum cycle mean of the S-avoiding window graph** (finite graph ⟹
the Motzkin density is attained by a cycle — Cantor–Gordon periodicity for free, no
period cap), decided per set by one Bellman–Ford positive-cycle test at mean = M, with
the witness independent set certifying μ ≥ M per set. Complete sweep, exact both sides:

> - **|S| = 3, max ≤ 14: μ = M on all 325 primitive sets.** With χ_f = 1/μ ≤ χ_c ≤ 1/M
>   this *pins* **χ_c = 1/M for every 3-element set in range** — the ladder is rigid
>   through |S| = 3.
> - **|S| = 4, max ≤ 14: 26 separations μ > M.** First: **S = {1,3,4,5}, M = 2/9 (t = 4/9),
>   but {0,2} mod 8 avoids all differences at density 1/4 > 2/9.** All certificates
>   machine-verified; this one hand-checkable in a minute.
> - **|S| = 5, max ≤ 12: 67 separations.** Total 93.
> - **χ_c = 1/M (GRAPH-LRC's bold form) is refuted at the same set:** {0,2},{1,3},{4,6},
>   {5,7} mod 8 properly 4-colors Cay(ℤ, ±{1,3,4,5}), so χ_c ≤ 4 < 4.5 = 1/M.

So the S141 hunt's "μ = M on all 8 test sets" was a small-sample accident (all its test
sets happened to be collapse cases), and the Haralambis-expectation flag was right: the
fractional/Motzkin relaxation has genuine slack from four generators on. **Consequences
for the program:** the LP/fractional route cannot prove LRC(14) by itself — a 13-set
could conceivably have M < 1/14 ≤ μ, invisible to MOTZKIN-14 — *unless* one proves a
**tight-locus collapse lemma** (μ = M on near-tight covering sets). Notably the S141 data
already showed the four repo extremal families (GW, prim-sat, parity record, deep well)
sitting exactly at the fractional bound: the relaxations appear tight precisely at the
binding configurations. That refined conjecture — collapse at the tight locus, slack in
the bulk — is now the honest form of the ladder question, and it rhymes with everything
else in this problem: exactness lives only on the extremal shell.

**Where the slack comes from (factoid wiring).** Every |S| = 4 exception contains additive
triples (sum-closure: {1,3,4,5}: 1+3=4, 1+4=5; {2,3,5,8} and {3,5,8,11/13}: Fibonacci
runs; {1,11,12,13}: 1+11=12, 1+12=13; …). Additive triples a+b=c in S are exactly
**triangles in the distance graph** G_S (x, x+a, x+a+b), i.e. the minimal vectors of the
relation lattice (factoid D1: kissing number = additive energy shells). Short relations
create local chromatic structure that a rotation (one global witness) cannot exploit but
a combinatorial periodic set can — the Motzkin LP wins exactly where the relation lattice
is rich. The old "Fibonacci is the covering-min's foil" (C1) returns with a sharper job:
Fibonacci-closed sets are the *maximal ladder-slack* sets, the opposite pole from the
AP's slab-realizable optimum. (Observation-level; the exception census is in the .out.)

## 3b. The centerpiece: the GW family separates — the tight locus splits

Testing the tight-locus collapse directly (part-4 run) gave the session's sharpest result:

> **μ(GW) = 1/13 EXACTLY, while M(GW) = 1/14.** The Goddyn–Wong family {1..11, 13, 24} —
> THE tight non-AP instance of LRC(14) (THM-612, CASE-tight-locus-has-GW) — admits the
> periodic avoiding set **{0, 12} mod 26** (all differences ≡ ±12 mod 26, never in GW;
> hand-checkable) of density 1/13, and no periodic set does better (max-cycle-mean test at
> 1/13: no positive cycle). So χ_f(G_GW) = 13 < 14 = 1/M: **the fractional bound is a full
> rung loose at the second tight family.**
>
> Same split at k = 4: tight {1,2,3,4} collapses; tight **{1,3,4,7} — the Lucas sequence,
> sum-closed (1+3=4, 3+4=7) — has M = 1/5, μ = 1/4 exactly** ({0,2} mod 8). Both
> separations are exactly **one unit-fraction rung**: μ = 1/(1/M − 1).
>
> Meanwhile prim-sat (M = 1/13), parity record (M = 1/12), every AP prefix k ≤ 8, and all
> 13 single-move neighbors of GW **collapse** (μ = M exactly). GW is an isolated
> separation point in its own neighborhood: **tightness-for-M and separation-for-μ
> co-occur exactly at the non-AP tight instance.**

Two readings. (i) For the program: the tight locus {AP, GW} — which every extremality
argument must treat as twin equals — is *split* by the fractional relaxation: the AP is
LP-faithful (its witness slab is Motzkin-optimal), GW is LP-slack (a combinatorial
periodic set beats every rotation slab by one rung). Any route to LRC(14) through
μ/χ_f-type relaxations is therefore structurally blind at GW; conversely, a proof
technique that *distinguishes* slabs from general periodic sets is exactly what the
GW-side needs. The sum-closure mechanism is visible in both separating tight instances:
GW's defining substitution 12 → 24 creates the closure **11 + 13 = 24** (and 24 = 2·12,
the 2-adic doubling move again; the optimal set's period 26 = 2·13 is its shadow); the
k=4 separator is the Lucas run. (ii) For the S141 record: "the tight families sit AT the
fractional bound" was a lower-bound certification (witness slab, μ ≥ M) misread as
equality — the S144 engine shows equality holds at AP/prim-sat/parity-record and *fails*
at GW. Cite-check before any external claim: GW entered the literature from circular-flow
theory (Goddyn–Wong 2006, BGGST 1998), and χ_f of these distance graphs may be known
there or in Liu–Zhu; the named open question left here is **χ_c(G_GW) ∈ [13, 14] —
which?** (If χ_c(G_GW) = 14, the circular rung stays faithful at the tight locus where
the fractional rung fails — it would single out χ_c as the right graph invariant for
LRC; if χ_c < 14, every graph relaxation is blind at GW.)

And the tournament half is standing right there: sum-closure = directed 3-cycle structure
= the same odd-cycle objects the OCF/H(T) side is built from (a session on the "sum
tournament" of a speed set — orient a+b → c — may make this precise; the μ > M sets
should be the odd-cycle-rich ones, the μ = M sets the transitive-like ones).

## 4. Corrections filed / absorbed

- HYP-5117's per-q extremal lemma: refuted as stated (this session's primitive
  near-affine adversaries; convergent with mac-mini-S54's true-window refutation and
  attribution-saturation diagnosis). INDEX updated under HYP-5137 with cross-refs.
- HYP-4972 (S141): the μ = M collapse hope is now resolved *negative* (93 separations);
  the |S| ≤ 3 rigidity + tight-locus collapse is what survives. INDEX updated.
- The S141/S143 `witness()` helper carried the pre-S143 divisor-closure bug (raw pair
  sums, not divisors — returned None on all-odd sets); fixed in-module, matching the
  documented M_exact fix. The recurring pattern's third appearance; engine authors: the
  candidate q-set must always be divisor-closed.
- MISTAKE-123 (kps-S73, absorbed): all (A′)-side bars in this note are stated without
  reliance on the circulated positivity-T_k table.

## Ledger

- EXACT/NEW: W_q(AP_13) = 1/7, 5/77, 3/35, 8/147, 23/294, 4/245 (Farey-flank closed
  form; identity: closing rate = flank − q); V_r(AP_13) = 17/99, 41/539, 94/1155,
  149/2695, 65/1386, 6/539; T_6(AP_13) = 6/539; both sum exactly to 477/1078.
  **μ(GW) = 1/13 exactly ({0,12} mod 26); μ({1,3,4,7}) = 1/4 exactly ({0,2} mod 8).**
- REFUTED: per-q W_q AP-minimality in primitive normalization (near-affine adversaries);
  T_4 AP-maximality; μ = M collapse (93 exact separations); χ_c = 1/M at |S| = 4;
  **the tight-locus collapse itself, at GW** — the surviving equality locus is the
  AP-side (AP prefixes, prim-sat, parity record, GW's whole single-move neighborhood).
- SURVIVES (new, open): T_5/T_6 deep-resonance-tail AP-maximality; |S| = 3 ladder
  rigidity (μ = M ⟹ χ_c = 1/M, all 712 sets through max = 18); the AP-locus collapse +
  one-rung law μ = 1/(1/M − 1) at the sum-closed tight instances; χ_c(G_GW) ∈ [13,14].
- Files: `lrc_per_q_audit_opus_S144.py`, `lrc_per_q_nearaffine_and_Vr_opus_S144.py`,
  `lrc_resonant_tail_law_opus_S144.py`, `lrc_mu_eq_M_maxcycle_opus_S144.py`,
  `lrc_tight_locus_collapse_opus_S144.py` (+.outs);
  `witness()` fix in `lrc_graph_interpretation_ladder_opus_S141.py`.
- Builds on: THM-637 roof (the closed forms are its Farey combinatorics), kps-S72
  (the program audited), mac-mini-S54 + klein-S173 (concurrent, reconciled), kps-S67
  (σ-grading, reaffirmed), boxeph-S1 (anchor reduction), S134/S135/S141/S143 (own thread),
  factoids C1/C3/D1/D2 + HYP-3900 far-peel (wired in above).
- Does NOT prove LRC(14), R2, or (A′). It kills two false hopes precisely, replaces them
  with one invariant target (deep-resonance tail) and one refined conjecture (tight-locus
  collapse), both with exact AP-side constants.
