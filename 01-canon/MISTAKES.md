# Mistakes Log

**Purpose:** Every error that has been made and corrected — with enough context that no Claude instance ever repeats it. Read this before doing any computational or proof work.

Format per entry:
- What was assumed / done
- Why it was wrong
- The correct framing
- Impact on existing results
- Source (who found it, when)

---

## MISTAKE-093 (2026-07-01, klein-S89, correcting klein-S88's HYP-3843) -- "identity on (1/15, 1/14], r* = 1/15" was WRONG: the identity window is [2/29, 1/14]. Equality of two piecewise-linear functions at all candidate breakpoints does NOT give identity between them when one function kinks where the other is straight.

**What was done (wrong).** klein-S88 verified Lambda_AP == Lambda_GW at every candidate breakpoint of both
profiles above 1/15 plus ONE midpoint and concluded identity on (1/15, 1/14] with r* = 1/15. But GW has a
REAL kink at 2/29 (gap-death of a (5,24) pair, d=2), the profiles agree AT 2/29 and above and DIFFER on the
sliver (1/15, 2/29); the midpoint probe (29/420 = 0.06905) sat just ABOVE 2/29 = 0.06897 and missed it.
mac-mini-S94's breakpoint list (2/33, 2/31, 2/29) was correct all along.

**Why it was wrong.** Two piecewise-linear functions equal at all radii of the JOINT kink-candidate list can
still differ between consecutive candidates: equality at a kink plus different incoming slopes = divergence
on one side only. Equality-at-candidates is not identity.

**The correct framing.** r* = 2/29; the shared final window is [2/29, 1/14] (width 1/406, still positive:
the last mile IS single-valued, just shorter); on (1/15, 2/29) Lambda_GW > Lambda_AP, so the AP alone
carries the envelope there -- mac-mini S94(3)'s second-order tie-break confirmed in full.

**Impact.** HYP-3843 corrected in place (title/status/part 1 + INDEX line); no downstream result used
r* = 1/15 quantitatively. HYP-3841's ladder K values are UNAFFECTED (that code asserts per-segment midpoint
linearity on each function; the assertion never fired). Rule: any identity claim between piecewise-linear
profiles must carry per-function per-segment midpoint assertions.

**Source.** klein-2026-07-01-S89, k0_final_window_lemma_klein.py part (d)(i) (local-linearity probe at 2/29).

---

## MISTAKE-092 (2026-07-01, opus-S33) -- incomplete kink taxonomy for the cluster density D_c(t): same-endpoint (center-coincidence) collisions t = m/delta were omitted from the breakpoint grid. Caught same session via a hand-derivable slope; the error was silently CONSERVATIVE.

**What was done (incomplete).** Computing the renormalized cluster density D_c(t) as piecewise linear with
breakpoints only at opposite-endpoint collisions (7m+-1)/(7 delta). The structure check then reported a
one-sided slope 4 at a zero where the hand derivation (gap-rate sum) gives 6 -- the discrepancy exposed a
missed kink at t = 1/6 (= m/delta, two arc centers coinciding), a CONCAVE PEAK of D between grid points.

**Why it mattered / why it did not.** All missed kinks are concave peaks (center-coincidence = union-measure
minimum = D maximum), so linear interpolation between the sampled kinks sits pointwise BELOW the true D:
the computed rearrangement floors F_j were valid LOWER bounds and every conclusion survived. But "exact"
would have been a false label, and in another context a one-sided sampling error would not be benign.

**Correct framing.** D_c's kinks come in exactly two types (the same taxonomy as THM-592(iii) in the
r-variable): opposite-endpoint collisions delta*t == +-1/7 (gap birth/death, convex) and same-endpoint
collisions delta*t == 0 (overtaking/coincidence, concave peak). Enumerate BOTH: t in {(7m+e)/(7 delta),
e in {-1,0,+1}}. Lesson (cf. MISTAKE-086, incomplete breakpoint set for M(S)): enumerate collision TYPES
from the geometry, not from the first pattern that works; keep one hand-derivable value (a slope at a
structured point) as a canary in every piecewise-linear computation. Source: opus-2026-07-01-S33.

## MISTAKE-090 (2026-07-01, opus-S32) -- the 11-core far-element lever asked for a "UNIFORM ARC-COUNT bound c' <= const" that is FALSE as needed; one-at-a-time peeling was the wrong induction. Fixed by the SIMULTANEOUS peel.

**What was assumed (WRONG, opus-S31).** That closing the multi-outlier r=2 residual needs a uniform bound
c' <= const on the number of lonely-set components of an 11-core PREFIX that already CONTAINS a large
element w1 (so the next peel's error c'/(7 w2) vanishes). S31's census-to-completion .out states this as
"REMAINING ... needs the UNIFORM ARC-COUNT bound (c' <= const, not c' <= sum_v v)".

**Why wrong.** No such bound exists: intersecting a compact core's lonely set (c'' arcs, measure m'')
with A_w shatters it into c ~ w*m'' + 2c'' components -- VERIFIED growing 16,24,26,38,66 for w=20..400
(lrc14_simultaneous_peel_lemma_opus_20260701.out, part L). The one-at-a-time induction m_3<...<m_11
forces the far-containing prefix's arc count into the error term, and that count grows linearly in the
peeled element. The requested lemma was unprovable because it is false.

**Correct framing (the fix).** Peel ALL far elements AT ONCE by a union bound (simultaneous-peel lemma,
HYP-3900): meas(L_C) >= (1-j/7) meas(L_low) - (2 c_low/7) sum_{w in F} 1/w -- only the COMPACT part's
arc count ever enters, and at a multiplicative gap cut the scales cancel (error <= 22j/(7 Lambda),
uniform, scale-free). The union bound dies exactly at j=7=1/(2r) (seven speeds can cover), which is the
TRUE boundary of the method; the residual is gap-free >=7-element clusters, not an arc-count problem.

**Impact.** S31's Part B verdict stands (its checks were correct); only its stated "remaining ingredient"
is replaced. The j<=6 multi-outlier tail is now closed by HYP-3900's guard table (all margins positive,
min 0.013 at j=5,6). Lesson: when an induction forces you to control a quantity that grows along the
induction, change the induction, not the estimate. Source: opus-2026-07-01-S32.

## MISTAKE-085 (2026-06-22, mac-mini-S30) -- Node 1 boundary-core "closure" OMITTED G_P (conflated the L-cluster maxgap with the P-small-part safe set). Caught same session.

**What I claimed (WRONG).** That the LRC(14) boundary core {t,..,12t,V} closes via rho_K>0
for all V/t>12, by the s≈0 "cluster collapse" (teeth -> 1/2, maxgap 1/2 > 2/7).

**Why wrong.** THM-527 splits runners: P = small part (speeds <=13) handled by `G_P = {x:
||p x||>=1/14}`, L = large cluster (>13) handled by `maxgap{frac(e_L x)}>2/7`. Good period =
`x in G_P AND maxgap(L)>2/7`. I computed maxgap over ALL co-offsets, omitting that the small
runners belong to G_P. The s≈0 period FAILS G_P (small p: ||p/(2V)|| < 1/14). VERIFIED the TRUE
rho_K (with G_P) for {1,..,12,V} is 0 at V=29,43,71 (not >0). 

**Correct framing.** rho* = meas(G_P cap {maxgap(L)>2/7}) = the OVERLAP. meas(G_P)>0 via proven
LRC(|P|<=13); meas(maxgap(L))>0 via three-distance; the overlap > 0 = the decorrelation floor
(Node 2/3). The discretization lemma rho_K >= rho* - arcCount/Vmax stands. 

**Lesson.** Don't conflate the cluster/far split; G_P (small part) is the real constraint, not
the maxgap. Verify the FULL criterion before claiming closure. (Cf. MISTAKE-084, same pattern:
over-claiming a closure by omitting a constraint.) Source: mac-mini-S30.

## MISTAKE-084 (2026-06-22, mac-mini-S27) — I wrongly claimed the p0 route "fails at k=8"; in fact the witness floor needs only `> 0`, not `≥ m_P` (the floor VALUE is not load-bearing). Caught & retracted same session.

**What I claimed (WRONG).** That the p0-wide-bound route (HYP-2832) fails because
`cap_8 − p0(consec_8) = 319/5880 = 0.0543 < m_P = 14249/252252 = 0.0565`, hence the
spreading lemma is REQUIRED. I even filed a court case against HYP-2832.

**Why it was wrong.** I assumed the witness floor must REACH `m_P`. It does not. The
floor is consumed through exactly one lemma:
`witness_floor_positive (hfloor : witnessMP ≤ witnessG2) : 0 < witnessG2 :=
lt_of_lt_of_le witnessMP_pos_real hfloor` (LRCFourteenSkeleton:239) — uses ONLY
`witnessMP > 0`. Then `hpartA : 0 < witnessG2 → Mreach ≥ 1/14` needs ONLY strict
positivity. So ANY positive lower bound on `witnessG2` suffices; `m_P` is a
non-load-bearing placeholder. The p0 route's `0.0543 > 0` is plenty; `0.0543 < 0.0565`
is irrelevant.

**Correct framing (the useful takeaway).** **The witness floor is ROBUST: the proof
depends only on the POSITIVITY of `witnessG2`, not on any particular floor value.**
Both routes work — p0 (`0.0543 > 0`, no spreading lemma) and NU (`0.322 > 0`). The
spreading lemma stays OPTIONAL (HYP-2832 is valid). Only nuance: the skeleton's nodes
are literally stated with the constant `witnessMP = m_P`, so to use the p0 route one
either proves `0 < witnessG2` directly or lowers the placeholder — cosmetic, not a real
obstruction.

**Lesson.** Before declaring a bound "insufficient," CHECK what magnitude the
downstream consumer actually requires. Here `> 0` was all that was needed; I compared
against a constant that was never load-bearing. Court `CASE-p0-route-insufficiency`
WITHDRAWN. Source: mac-mini-S27 (error and self-correction same session).

## MISTAKE-083 (2026-06-21, claude-opus-S1, LRC doublet signed-bound script) — off-by-one in the doublet base (`range(k-1)` gave k+1 speeds), k-labels shifted +1 (found & corrected SAME session)

**What was done.** In `lrc14_doublet_signed_bound_claudeopus_0621.py` and
`lrc14_doublet_periodicity_claudeopus_0621.py` the doublet config was built as
`list(range(k-1)) + [M, M+1]` = `{0,…,k-2} ∪ {M,M+1}` = **k+1 speeds**, but labeled `k`.
So the first-pushed signed-bound table (commit ebf1a8dbe) had every k-label one too LOW
(the "k=8..12" rows were really k=9..13).

**Why it was wrong.** The `[m-2,2]` genuine-wide doublet maximizer (HYP-2794) has base
`{0,…,k-3}` (k-2 elements incl 0) + doublet (2) = k speeds; the correct generator is
`range(k-2)`, not `range(k-1)`. The separate maximizer script was CORRECT (it filtered
`len(E)==k`), which is how the mismatch surfaced: its k=9 max 0.28976 equaled the
signed script's "k=8" 0.28977.

**Correct framing.** Base = `{0,…,k-3}` = `range(k-2)`. Corrected table (exact):
sup_M p0 = {0.1446(M20),0.2898(M21),0.4425(M21),0.5211(M21),0.5897(M21)} for k=8..12,
cap−sup = {+.237,+.204,+.162,+.204,+.267}; sup_M|M·error| = {1.39,1.27,1.34,1.41,1.47},
overcount BV/signed = {106×,154×,188×,224×,262×}. The doublet IS the genuine-wide
maximizer for k=9..12; at k=8 a 3-cluster (0.172) slightly beats it.

**Impact.** QUALITATIVE conclusions UNCHANGED and slightly CLEANER (corrected margins are
larger, so `sup|M·error| < 15·margin_2` now holds at ALL k incl k=9). HYP-2794 detail +
INDEX corrected same session; correction broadcast sent. The signed-cancellation finding
(overcount ~100–260× = THM-563 analogue) stands.

**How it was caught.** Re-deriving the speed count m=k-1 vs the maximizer configs; a direct
`len(E)` check confirmed `range(k-1)+doublet` = k+1 elements.

## MISTAKE-082 (2026-06-20, kind-pasteur-S21, LRC partition-function wide bound) — the decorrelated FLOOR p0_decorr mistaken for a majorant of finite wide p0 (the "0.19 comfortable budget" was illusory)

**What was done.** Framing the LRC wide-residual via the apex-prime partition function (HYP-2694), I computed the single-block DECORRELATED cover `p0_decorr(m) = 0.1925/0.3056/0.4123/...` (the M->infinity Weyl limit), proved it `< cap_k` with margin `>= 0.1886`, and claimed the wide bound therefore has a "comfortable 0.19 budget" — estimating the decorrelation error `e(E)=p0(E)-p0_decorr` as `~0.01` (from the single sample `[0,19..25]`, e=0.0095).

**Why it was wrong.** `p0_decorr` is the decorrelated FLOOR (M->infinity), NOT a majorant of the finite-scale `p0(E)`. Adversarial verification: finite wide sets exceed `p0_decorr` by `+0.05 .. +0.13` — a near-pinned stretched consec at moderate span has `p0 ~ 0.318` vs decorr `0.193` (k=8). So the chain `p0(E) <= p0_decorr < cap` is INVALID, and the real wide cap-margin is the FINITE-CHECK level (~0.06-0.13), not 0.19. The decorrelation error is the genuine remaining content, not a negligible `~0.01`. Same class as MISTAKE-080 (the `wide => p0<=Q(k-1)` over-claim): a value approached in a limit / on a sample treated as a uniform bound.

**Correct framing.** `p0(E) = p0_decorr(shape) + e(E)`; `p0_decorr` is exact (the cut-space extremum), but the wide bound REQUIRES bounding `e(E)` (the cycle-space interaction). Concurrent codex THM-557 (pushed first, S61) had this RIGHT all along: it bounds the single-block error by the diagonal-freeze `|p0(E_M)-D_m| <= 7*C(m,2)/M` (sharper than my independent Lemma DE `49m^2/(6M)`), giving M-cutoffs ~779-1369 + a finite small-M exact check — so the single-block branch IS rigorous. The OPEN part is the MULTI-cluster (>=2 far) `e` = OPEN-Q-108 (joint Erdős–Turán–Koksma constant, iterate of the PROVED single-far THM-546 `|Delta_w|<=(6/49)V/w`; codex's compression cone HYP-2696/2697/2698). VERIFIED 10.5k exact wide rows, 0 cap-violations, NOT closed.

**How it was caught.** My S21 finalization workflow's cover-bound verifier (holds=false, high conf, independent exact engine) showed `p0_decorr` is not a majorant; cross-checked against codex's concurrently-pushed THM-557 which already used `D_m + error`.

**Impact.** My duplicate THM-559 was DELETED (codex THM-557 is the authoritative single-block treatment). HYP-2694 kept at codex's version. HYP-2675 status: VERIFIED (cap-level, margin >= 0.05, 0 violations) + REDUCED to OPEN-Q-108; NOT proved. The lesson recurs (cf MISTAKE-080): decorrelated/limit extrema are cut-space facts; the bound needs the cycle-space error — and concurrent agents on the same prompt converge, so cross-check before claiming.

## MISTAKE-081 (2026-06-20, kind-pasteur half-tiling session) — the canon claim "SC(n) = A000568(n−1)" is FALSE (holds only at n=4,6 in the tested range)

**What was claimed.** `07-reflections/two-models-staircase-recursion.md` (section "SC(n)=A000568(n−1) in the Tiling Model", lines ~304–323, 344) asserts that the number of self-complementary (self-converse) tournament isomorphism classes equals `A000568(n−1)` (= the number of tournament iso classes on `n−1` vertices). The half-tiling/orbifold thread (HYP-2686) tested this while looking for a clean SC-halving bijection.

**Why it is wrong.** Direct iso-class enumeration (independently re-verified, `04-computation/half_tiling_verify_contested_kps.py`) gives the self-converse counts `SC_n = 2, 2, 8, 12` for `n=3,4,5,6`, while `A000568(n−1) = 1, 2, 4, 12`. They AGREE only at `n=4` and `n=6` and DISAGREE at `n=3` (`2≠1`) and `n=5` (`8≠4`). (Sanity: the same script reproduces the total iso-class counts `A000568(n)=2,4,12,56` exactly, so the canonicalization is correct.) The correct self-converse sequence is `2,2,8,12,88,176` (n=3..8) — matching the value already recorded in `07-reflections/unlocking-gn-at-all-n.md` (`SC_n` row), and equal to OEIS A002785 (number of self-converse/self-complementary tournaments) up to indexing. So the repo already contained the right numbers in one place and a wrong identity in another.

**Correct framing.** `2^half(n)` (grid-symmetric tilings) is a LABELED, base-path-fixed fixed-point count (the Burnside input `Fix_anti(φ)` restricted to the tiling frame), NOT the iso-class count `SC_n`. There is no clean `|SC(n)| = A000568(n−1)` identity; the genuine relation is the Burnside one `V_merged = (A000568(n) + SC_n)/2` (with `SC_n = 2,2,8,12,88,176`), and the explicit class-level SC-halving bijection (open-Q #4 of two-models-staircase-recursion.md) remains OPEN — the half-tiling supplies the fixed-point count, not the bijection.

**Impact.** A correction note has been added at the head of the affected section in `two-models-staircase-recursion.md` pointing here (additive, non-destructive; the surrounding "intrinsic halving" intuition is kept as motivation, the specific count identity is flagged false). HYP-2686 records the verified counts. No theorem depended on the false identity (it lived only in that reflection's prose).

**How it was caught.** HYP-2686 orbifold thread + independent re-verification (the workflow's Burnside adversarial verifier died on an API error, so kind-pasteur re-checked `SC_n` by direct iso-enumeration n=3..6).

**Source.** kind-pasteur-2026-06-20 (half-tiling framework session).

## MISTAKE-080 (2026-06-20, kind-pasteur-S19, LRC(14) sector route) — the decorrelated LIMIT Q(k-1) mistaken for a FINITE upper bound on wide p0 (plus a pinned-0 modeling omission)

**What was done.** Attacking the sole LRC(14) residual HYP-2675 (`span(E)>14 ⟹ p0(E)≤cap_k`) via cross-scale decorrelation, I computed the decorrelated model `p0_decorr(base, nfar)` and found `sup_shapes p0_decorr = Q(k-1) = Plat(consec_{k-1})` exactly (0.197/0.362/0.448, k=8/9/10), achieved at consec_{k-1}+1 stranger. I recorded this as the crux "wide ⟹ p0 ≤ Q(k-1) < cap" and called it the proof of the decorrelated bound.

**Why it was wrong.** `Q(k-1)` is the decorrelated **LIMIT** (cluster scale gaps → ∞), NOT a finite upper bound. Counterexample (exact, adversarial-workflow-found, re-verified): `E=[0,19,20,21,22,23,24,25]` (k=8, wide) has `p0 = 9524621/47108600 = 0.20218 > Q(7) = 0.19660`. Finite-scale wide p0 sits ABOVE Q(k-1). Two compounding errors: (1) limit-vs-finite (a value approached in the limit treated as a uniform bound — same class as MISTAKE-079); (2) a MODELING omission — my `p0_decorr` assumed the dense base cluster CONTAINS the pinned 0 (`frac(0·x)=0` wastes the 0-element on sector 0). A cluster NOT containing 0 (e.g. `{19..25}={0..6}+19`) has ALL points SWEEPING and covers the 6 INNER sectors BETTER than consec_7: `p0({19..25})=0.202 > p0(consec_7)=0.148`. So the true `sup_shapes p0` exceeds Q(k-1).

**Correct framing.** The true residual is `wide ⟹ p0(E) ≤ cap_k` (margins ≥0.30 to cap, much looser than the thin margin to Q). Q(k-1) is the decorrelated limit only. The proved pieces survive: the cardinality lemma (cluster size ≤5 ⟹ measure-0 full coverage), the SHARPER comb bound `|Δ_w| ≤ 2c1(E')/(7w)` (THM-546/547, c1=miss-1 component count), and `Q(k-1)<cap_k`. The remaining open lemma is the CAP-level (not Q-level) joint multi-cluster Erdős–Turán–Koksma bound + the balanced-wide branch.

**How it was caught.** The HYP-2675 proof workflow's adversarial verifiers (3 of 4 returned holds=false on the strong (W*) form `wide ⟹ p0≤Q(k-1)`), with the exact counterexample [0,19..25]. Independent re-verification confirmed.

**Impact.** `05-knowledge/results/lrc14_h2675_decorrelation_foundation_kps.md` corrected (the "= Q(k-1)" crux marked REFUTED for the finite bound, kept as the limit). The commit "decorrelated bound EQUALS the plateau Q(k-1)" is superseded. The sector-route residual is `wide ⟹ p0≤cap`, NOT `≤Q(k-1)`. HYP-2675 remains CONJECTURE (0 counterexamples in 10^4–10^5 wide sets, margins ≥0.30).

## MISTAKE-079 (2026-06-19, kind-pasteur-S9, LRC max-min spectral gap) — "covered below the level-a mediant" misread as a DEEPER dip when it actually meant COLLAPSE TO THE FLOOR

**What was done.** Investigating whether the LRC max-min gap `g(k)=σ_2(k)-1/(k+1)` dips below `Θ(1/k^2)`, I built a covering test `M(S) <= r ⟺ danger arcs of radius r cover [0,1)`. For the family `F(k,a)={1,…,k-2,k,a(k-1)}` at `(a=5,k=211)` and `(a=6,k=2311)` the test returned `M < a/(a(k+1)-1)` (e.g. `M < 5/1059`). I read this as "the family dips even DEEPER than level a" and concluded `a_max(k)` is UNBOUNDED (level grows with `ω(k-1)`), hence `liminf g·k^2 = 0` — and wrote it into THM-539/HYP-2623.

**Why it was wrong.** `M < a/q` only says `M` is below that mediant; it does NOT say `M` is a *deeper mediant just above the floor*. In fact `M(F(211,5)) = 1/212 = the FLOOR exactly` — the family COLLAPSES to a tight configuration (`g=0`), the opposite of a deep dip. The covering "M below 5/1059" was consistent with `M = 1/212` (since `1/212 < 5/1059`), but I never computed the exact `M`; I extrapolated a monotone-deepening pattern from `g·k^2 = 0.50, 0.33, 0.24, <0.20, <0.17` that was really `0.50, 0.33, 0.24, 0(tight), 0(tight)`.

**Correct framing.** The natural family realizes levels `a=2,3,4` only (a=3 at `k≡7,13,19,25 mod30`; a=4 at `k≡1 mod30`); at `a=5` it becomes TIGHT whenever `2·3·5·7|(k-1)`. Confirmed dips give `g·k^2 ∈ [≈1/4, 1/2]`. Whether ANY family reaches `a>=5` (hence whether `g=Θ(1/k^2)` uniformly) is OPEN. This realigns with codex S16/S17 ("no `o(1/k^2)` dip found").

**How it was caught.** The adversarial verification Workflow's exact-`M` agent (independent code, covering + binary search) returned `M(F(211,5))=1/212` cleanly; an own recheck confirmed `F(k,5)` is tight for `k-1` divisible by `2·3·5·7`.

**Impact.** THM-539, HYP-2623, SESSION-LOG, memory corrected: `a=3,4` dips stand (exact, infinite families); the `a_max`-unbounded / `liminf g·k^2=0` claim is RETRACTED to OPEN.

**Lesson.** When an inequality `M < c` localizes a quantity, COMPUTE the quantity exactly before naming what it is — `M < c` does not distinguish "a slightly smaller special value" from "collapsed to the global minimum." Echoes MISTAKE-073/076/078: a bound or a sampled trend is not the exact structure. Run the exact computation at the extrapolation point, especially before declaring an asymptotic (`liminf=0`).

**Source:** kind-pasteur-2026-06-19-S9; caught by the adversarial verification workflow (E2 agent) + own recheck.

---

## MISTAKE-075: The LRC(14) inf L estimate (≈0.0052) was a restricted-search artifact — the true inf-relevant extremizers are sporadic-tight perturbations, not the multiple-of-14 family

**Date discovered:** 2026-06-16 (kind-pasteur-S7)
**Found by:** kind-pasteur, via exact rational arc-sweep of the lonely measure (THM-522)
**Affects:** THM-518's stated `inf L ≈ 0.0052`; OPEN-Q-097/104's extremizer family; the "primitive multiple-of-14" extremizer search.

- **What was done:** the extremizer search for `inf_S L(S)` (lonely measure, LRC(14)) was restricted to the family `({1..13}\{j})∪{14m}` — interior-drop AP cores with a **multiple-of-14 stranger** — giving `inf ≈ 0.0052` (extremizers `{1..13}\{12}∪{84}` etc.).
- **Why it was wrong:** `36` is not a multiple of 14, so the search never tested it — yet the **minimal single-element perturbation of the tight AP `{1..13}` is `12→36`**, giving `{1,2,…,11,13,36}` with `L = 1/1260 ≈ 0.000794`, ~6.7× below 0.0052 (verified exact arc-sweep + independent fine grid). The restriction to multiple-of-14 strangers blinded the search to the perturbations of the **sporadic tight configs** (e.g. `{1..11,13,24}`, `L=0`), which open the smallest lonely measures.
- **Correct framing:** `inf L` is governed by the minimal perturbation of the **FULL tight locus** (AP + sporadics), not the AP/14m family. `inf L ≤ 1/1260` (THM-522 D); conjecturally `=1/1260` (HYP-2561). The lonely measure is quantized `L∈(1/(14 lcm))ℤ` and scale-invariant `L(cS)=L(S)` (THM-522).
- **Impact:** the LRC(14) difficulty constant is ~6.7× smaller than recorded; OPEN-Q-097's program must control all tight configs (not just dilated-AP cores). NOT a disproof of LRC(14) (every config is loose, `L>0`).
- **Lesson (= MISTAKE-073 recurring):** sweep the FULL orbit of perturbations, not the convenient sub-family; the extremizer is consistently one orbit further out than the obvious family member (`0.0237` end-drop → `0.0053` interior-drop → `1/1260` sporadic-tight perturbation).

---

## MISTAKE-001: μ Computation Bug in Scripts 6-9

**Date discovered:** ~2026-02-26 (pre-Cowork sessions); logged 2026-03-05
**Found by:** Claude instance (Account unknown, pre-Cowork era), reported in MASTER_FINDINGS.md
**Affects:** Scripts 6-9 (sum_mu() function); does NOT affect scripts 1-5

### What was assumed
`ind_poly_at_2_restricted()` correctly computes μ(C) = I(Ω(T−v)|_{avoid C\{v}}, 2) by iterating over cycles using `T[perm[i]][perm[(i+1)%len]]`.

### Why it was wrong
`perm` is a permutation of T−v vertices that does NOT include v. The cycle-checking code uses the **full T matrix** instead of T−v. Specifically, for a cycle of length L:
- `T[perm[i]][perm[(i+1)%len]]` checks arcs in the full tournament T
- But cycles of T−v (cycles not using v) should be checked in T−v
- The wrap-around `perm[(i+1)%len]` may also have indexing issues for general L

### The correct framing
The restriction to T−v must be done BEFORE checking cycles. Any cycle-finding in `ind_poly_at_2_restricted()` must operate on T−v's adjacency matrix exclusively.

### Impact assessment
- **Scripts 1-5 (inshat analysis):** UNAFFECTED (do not use μ computation)
- **Paper's Claim A verification:** STATED UNAFFECTED by MASTER_FINDINGS (verification runs used a separate code path)
- **Paper's Claim B verification:** Needs confirmation — Claim B involves I(Ω,2) directly
- **The verification table in the paper (0 failures for Claim A at n≤6):** Stated to be valid

**RESOLVED** (DISC-001, closed by kind-pasteur-2026-03-05-S3): Independent verification using tournament_lib.py confirms the paper's results are valid. See 02-court/resolved/DISC-001-mu-bug-vs-verification.md.

### Lesson
When computing μ(C) for any cycle C in T:
1. First restrict to T−v (remove v and all incident arcs from the adjacency matrix)
2. Then find all odd cycles in T−v that are vertex-disjoint from C\{v}
3. Build Ω(T−v)|_{avoid C\{v}} from these cycles
4. Evaluate I(·, 2) on this restricted conflict graph

Never use the full T matrix when computing anything about T−v.

### Resolution
**Independent verification completed** (opus-2026-03-05-S1): tournament_lib.py implements mu(C) correctly per the above steps. Exhaustive verification of Claim A at n<=6 (196,608 pairs, 0 failures) confirms the paper's results are valid. See DISC-001 for the formal resolution.

---

## MISTAKE-004: RETRACTED — OCF IS a Valid Closed Form Over All Odd Cycles

**Date originally entered:** During file.txt exploration (pre-2026-03-05)
**RETRACTED by:** kind-pasteur-2026-03-05-S5 (DISC-002 resolved)
**Original claim:** H(T) = I(Omega(T), 2) where Omega(T) = ALL directed odd cycles is WRONG.

### Why the original claim was itself wrong

The alleged counterexample (T on 6 vertices with two disjoint 3-cycles C1, C2) computed I(Omega(T), 2) = 49 using MU WEIGHTS:
```
I_wrong = 1 + 2*mu(C1) + 2*mu(C2) + 4*mu(C1)*mu(C2) = 1 + 6 + 6 + 36 = 49
```
But the independence polynomial does NOT involve mu weights. The correct computation:
- alpha_0 = 1, alpha_1 = 2 (either C1 or C2), alpha_2 = 1 ({C1, C2} — vertex-disjoint, so non-adjacent)
- I(Omega(T), 2) = 1 + 2*2 + 1*4 = 9 = H(T) CORRECT

### The correct framing (replaces the false framing)

**H(T) = I(Omega(T), 2) IS a valid closed-form identity**, where:
- Omega(T) = conflict graph on ALL directed odd cycles of T
- Two cycles are adjacent in Omega(T) iff they share a vertex
- I(G, x) = sum_{k>=0} alpha_k * x^k is the plain independence polynomial (no mu weights)

This is equivalent to the recursive Claim A formulation — the closed form is obtained by unrolling the recursion. The mu weights mu(C) = H(T[V\V(C)]) arise in the recursion but NOT in the independence polynomial.

### Computational confirmation
H(T) = I(Omega(T), 2) verified exhaustively for n=3,4,5,6 (33,864 tournaments, 0 failures) by opus-2026-03-05-S2. Further confirmed by T_11 with H=95095 matching exact OCF calculation.

### What was confused (lesson for future agents)
The recursive mu-weighted formula H(T) = H(T-v) + 2*sum_C mu(C)*... uses mu weights at each step. When you UNROLL the full recursion, the mu weights become exactly the combinatorial weights in the independence polynomial (since mu(C) = H(T[V\V(C)]) which itself unrolls). The two formulations are equivalent; the closed form is NOT a non-recursive approximation. Do not confuse the per-step mu weights with the independence polynomial coefficients.

---

## MISTAKE-005: Cycle Bijection Under Arc Reversal Fails

**Date discovered:** During file.txt exploration
**Found by:** Claude instance (Account unknown)
**Affects:** Proof strategy for arc-reversal invariance D(T,v) = D(T',v)

### What was assumed
When T' is obtained from T by flipping arc i->j to j->i (with i,j != v), odd cycles through v containing i->j in T biject with odd cycles through v containing j->i in T', preserving V(C).

### Why it was wrong
A cycle C through v containing i->j can be written as v ->^{P1} i -> j ->^{P2} v. The "conjugate" C' would need j->i in T', requiring the path segments to be REVERSED. But reversing a directed path P2 does not generally give a valid directed path in the same tournament (arcs may point the wrong way).

### The correct framing
There is NO individual bijection C <-> C' preserving vertex sets. The arc-reversal invariance, if true, must hold as a SUM equality: sum_{C: i->j in C} H(T[V\V(C)]) = sum_{C': j->i in C'} H(T'[V\V(C')]). This is a weaker statement that doesn't require a cycle-by-cycle bijection.

### Impact
The arc-reversal invariance D(T,v) = D(T',v) remains the key unproved step for a general proof of Claim A. The cycle bijection approach is a dead end; a sum-level argument is needed.

---

## MISTAKE-002: Exact Path Formula H(T) = B_v + S_v + R_v

**Date discovered:** During script verification at n=4
**Affects:** Any approach trying to decompose H(T) exactly as B_v + S_v + R_v

### What was assumed
H(T) = B_v + S_v + R_v (exact equality)

### Why it was wrong
Verified FALSE: 96 failures out of 256 pairs at n=4.

### The correct framing
Only the PARITY version holds: H(T) ≡ B_v + S_v + R_v (mod 2). The exact identity is wrong.

### Impact
Similarly, S_v + R_v = 2Σμ(C) is FALSE (144/256 failures at n=4). These exact decompositions are not the right approach.

---

## MISTAKE-003: Per-Path Identity as a Path to Claim A for n≥6

**Date discovered:** During n=6 verification
**Affects:** Any proof strategy that relies on the per-path identity to prove Claim A for general n

### What was assumed
The per-path identity (inshat−1)/2 = Σ_{3-cycle embeddings} μ(C) could serve as an intermediate step toward proving Claim A for all n.

### Why it was wrong
The per-path identity fails for 2,758/9,126 ≈ 30% of (T,v,P') triples at n=6. It is provably a 3-cycle-only formula (by THM-005). Longer cycles contribute to Claim A's RHS but are invisible to the per-path identity.

### The correct framing
The per-path identity is valid and useful for n≤5 (where it implies Claim A). For n≥6, a different per-path formula is needed — one that accounts for contributions from all odd cycles, not just 3-cycles. See OPEN-Q-004.

### Impact
Any proof of Claim A for n≥6 must go beyond the per-path identity. The five proof strategies in the paper are the current best alternatives.

---

## MISTAKE-006: c₉/C(11,9) = c₇/C(11,7) "Ratio Coincidence" Has No Structural Basis

**Date discovered:** 2026-03-05
**Found by:** kind-pasteur-2026-03-05-S1 (from inbox document PALEY_T11_c9_ANALYSIS.md)
**Affects:** Any earlier estimate of c₉(T₁₁) ≈ 220 derived from this ratio

### What was assumed

The ratio c_k / C(11,k) is constant across odd k, noting c₇/C(11,7) = 1320/330 = 4 and (if so) c₉/C(11,9) = 4 → c₉ = 4 · 55 = 220.

### Why it was wrong

The ratio is NOT constant:
- k=3: c₃/C(11,3) = 55/165 = 1/3
- k=5: c₅/C(11,5) = 594/462 = 9/7 ≈ 1.286
- k=7: c₇/C(11,7) = 1320/330 = 4 ← coincidence only
- k=9: unknown

The pattern 1/3, 9/7, 4 is not monotone, not constant, has no arithmetic progression structure. The k=7 ratio is accidental. No structural reason was ever given for why the ratio at k=7 should equal the ratio at k=9.

### The correct framing

c₉(T₁₁) is determined by c₉ = (55/2)(h_QR + h_NQR) where h_QR = h({0,1}) and h_NQR = h({0,2}) are ham-cycle counts in explicit 9-vertex sub-tournaments (LEM-001). It must be computed directly.

### Impact

Any theorem or conjecture that relied on c₉ = 220 (or any ratio-derived value) should be flagged as unverified. **UPDATE (kind-pasteur-S2):** c₉ = 11055 (computed directly), and H(T_11) = 95095. CONJ-002 is fully refuted for p=11. The ratio estimate of 220 was off by a factor of 50.

---

## MISTAKE-007: Trace-Method Cycle Count Errors for c_6, c_7 in T_11

**Date discovered:** 2026-03-05
**Found by:** kind-pasteur-2026-03-05-S5 (from inbox documents more.txt, other.txt, stuff.txt)
**Affects:** Any hand computation of c_k(T_11) via tr(A^k) minus non-simple walk corrections

### What was assumed

The hand computation (inbox: more.txt) used the formula:
  k*c_k = tr(A^k) - N_k
where N_k = total non-simple closed walk contributions (Type at-v, Type A interior 3-cycle, Type B interior 4-cycle).

This gave c_6=1375 and c_7=1320.

### Why it was wrong

The non-simple walk corrections were computed incorrectly -- specifically the Type A and Type B contributions had arithmetic errors. The correct values, verified by direct DFS enumeration (other.txt), are:
  c_6(T_11) = 1595 (not 1375)
  c_7(T_11) = 3960 (not 1320)

### The correct framing

Direct DFS enumeration is more reliable than the trace correction method for large k. The correct complete cycle count table for T_11:
  c_3=55, c_4=165, c_5=594, c_6=1595, c_7=3960, c_8=7425, c_9=11055, c_10=10681, c_11=5505

These values are confirmed by the OCF identity: H(T_11) = 95095 = 1 + 2*(55+594+3960+11055+5505) + 4*10879 + 8*1155.

### Impact

The "corrected conjecture" H(T_11)=4455 (inbox: stuff.txt) was derived from the wrong c_7=1320 and was also false. The ratio 1320/330=4 was a coincidence. The actual sequence H(T_p)/|Aut(T_p)| = 1, 3, 9, 1729 for p=3,7,11 has no 3^k pattern.

### Lesson

For computing cycle counts in specific tournaments, use direct enumeration (DFS/backtracking) rather than eigenvalue-trace corrections. The trace method is valid in principle but requires extremely careful non-simple walk accounting. The LEM-001 formula c_9 = (55/2)(h_QR+h_NQR) is the correct approach (kind-pasteur-S2).

---

## MISTAKE-008: Even-Odd Split Claimed "Equivalent to OCF"

**Date discovered:** 2026-03-05
**Found by:** kind-pasteur-2026-03-05-S8 (reviewing opus-S4 output)
**Affects:** even-odd-split-lemma.md, OPEN-Q-009, TANGENTS T040

### What was assumed

The even-odd split (sum_{|S| even} Delta(S,R) = sum_{|S| odd} Delta(S,R)) was claimed EQUIVALENT to OCF in multiple places.

### Why it was wrong

The odd-S sum of Delta(S,R) = L_j(S)*R_i(R) - L_i(S)*R_j(R) is NOT the same as the cycle formula [g(S)-l(S)]*H(R). Specifically:
- L_j(S) = sum_a h_end(S,a)*T[a][j] does not include the T[i][first] factor
- g(S) = sum_{perms} T[i][first]*...*T[last][j] does include it

So even=odd gives delta_H = 2*(odd sum of Delta), but this odd sum is not the cycle formula. The even-odd split is a CONSEQUENCE of OCF, not equivalent to it. Proving even=odd would not prove OCF.

### The correct framing

The even-odd split is a necessary condition for OCF. It is a valid structural observation (verified n=5,...,8) but strictly weaker than OCF. To prove OCF via this route, one would additionally need to show that the odd-S sum of Delta(S,R) equals the cycle formula — which is essentially the full OCF identity.

### Impact

The even-odd split cannot serve as a standalone reformulation of OCF. It may still be useful as a sanity check or structural constraint, but claims of "equivalence" in the documentation have been corrected.

---

## MISTAKE-009: sympy_proof_n8.py Used Simplified n<=7 Formula

**Date discovered:** 2026-03-05
**Found by:** kind-pasteur-2026-03-05-S8
**Affects:** 03-artifacts/code/sympy_proof_n8.py (original version)

### What was assumed

The formula -2*sum(s_x*H(B_x)) + 2*(D5-C5) + 2*(D7-C7) applies at n=8.

### Why it was wrong

THM-013 explicitly states this formula FAILS at n=8 due to VD 3-5 cycle pairs. At n=8, the complement of a 5-cycle has 3 vertices, and H(3-vertex) can be 1 or 3 (not always 1 as at n<=7). The simplified formula doesn't weight 5-cycle contributions by their complement's H value.

### The correct framing

At n=8, must use the full A-clique formula: delta_I = 2*sum_C [gained-lost]*H(comp(C)). The script has been rewritten to use this formula.

### Impact

The original script would have produced WRONG results. Fixed immediately upon discovery.

---

## MISTAKE-010: Hereditary Maximizer Chain Claimed for ALL Maximizers at Odd n

**Date discovered:** 2026-03-06
**Found by:** kind-pasteur-2026-03-06-S18g
**Affects:** INV-044, T104 (hereditary maximizer chain)

### What was assumed

"At odd n, every vertex deletion from ANY maximizer gives the (n-1)-maximizer."

### Why it was wrong

Exhaustive check at n=5 shows: 64 maximizers (H=15), of which only 24 (regular score (2,2,2,2,2)) are hereditary. The remaining 40 (score (1,2,2,2,3)) have del_hs including H=3 values, NOT max H(4)=5.

Full data:
- n=3: 2/2 hereditary (all regular)
- n=4: 0/24 hereditary (all non-regular)
- n=5: 24/64 hereditary (only regular ones)
- n=6: 0/480 hereditary (all non-regular)
- n=7: 240/240 hereditary (all regular)

### The correct framing

**Only REGULAR maximizers at odd n are hereditary.** At n=3,5,7, the regular maximizers (those with score (k,k,...,k)) have ALL vertex deletions giving max H(n-1). Non-regular maximizers at n=5 do NOT have this property.

The pattern: regular maximizers exist at odd n and are vertex-transitive, so all deletions are isomorphic. One deletion being optimal implies all are.

### Impact

The investigation INV-044 and tangent T104 need correction. The hereditary chain is: T_7 -> T_6_max -> (chain breaks at n=5 for non-regular) or T_7 -> T_6_max -> T_5_regular -> T_4_max -> T_3. The chain from Paley T_7 goes through regular maximizers at each odd step.

---

## MISTAKE-011a: Transfer Matrix M = [[1,0],[0,-1]] Claim Was Wrong

**Date discovered:** 2026-03-06
**Found by:** opus-2026-03-06-S6
**Affects:** `04-computation/transfer_matrix_test.py`, INV-001 documentation

### What was assumed

The 2x2 transfer matrix M indexed by endpoints {i,j} of a fixed arc always equals [[+1, 0], [0, -1]]. This was stated as a verified fact in `transfer_matrix_test.py`.

### Why it was wrong

Exhaustive check at n=4 (64 tournaments x 5 arc pairs = 320 tests) shows 2199/2500 failures at n=4 with the original test parameters. Example: the transitive tournament on 4 vertices gives M[0,1] = 1, not 0. The values of M[a,b] vary widely depending on the tournament (observed values include -3, -2, -1, 0, 1, 2, 3 at n=5).

The original test used random tournaments with seed=42 and checked only specific arc pairs. The claim was never properly validated.

### The correct framing

The transfer matrix M[a,b] = sum_{S ⊆ V\{a,b}} (-1)^|S| E_a(S) B_b(V\{a,b}\S) satisfies:
1. **Symmetry**: M[a,b] = M[b,a] for all a,b (PROVED symbolically n<=7, verified numerically n<=8)
2. **Trace formula**: sum_a M[a,a] = H(T) for odd n, 0 for even n (THM-027)
3. **Off-diagonal sum**: sum_{a!=b} M[a,b] = 0 for odd n, 2*H(T) for even n

The individual entries M[a,b] are NOT always 0 or ±1. They are integers whose distribution depends on the tournament structure.

### Impact

INV-001 (proving transfer matrix symmetry) is still the correct goal — the symmetry M[a,b] = M[b,a] IS always true. But the diagonal values are NOT fixed constants, which changes the picture for any proof attempt that assumed diag(1,-1).

---

## MISTAKE-011b: Paley "tournament" at p = 1 mod 4 is NOT a tournament

**Date discovered:** 2026-03-06 (kind-pasteur-S25d)
**Found by:** kind-pasteur
**Affects:** nonham_vanish_n13_paley.py, palindromic_N_proof.py (n=13 and n=17 entries), paley_super_uniform.py, any script using QR mod p for p = 1 mod 4

### What was assumed

Computations labeled "Paley T_13" and "Paley T_17" used the quadratic residue set QR mod p as the circulant generator set. This was assumed to give a tournament.

### Why it was wrong

Paley tournaments exist ONLY at p = 3 mod 4 (where -1 is a quadratic NON-residue). At p = 1 mod 4, -1 is a QR, so QR is closed under negation: d in QR implies -d in QR. This means both T[i,j]=1 and T[j,i]=1 for QR-separated pairs, giving bidirectional edges. The resulting digraph is NOT a tournament.

Valid Paley tournament primes: 3, 7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83, ...
INVALID: 5, 13, 17, 29, 37, 41, 53, 61, 73, ...

### The correct framing

- "H(T_13) = 1,579,968" is the path count of a non-tournament circulant digraph, not a tournament.
- "H(T_17) = 5,587,473,776" is similarly invalid.
- The "Paley T_5" with S={1,4} is also NOT a tournament (p=5 = 1 mod 4).
- THM-052 proof is unaffected (algebraic, works for any circulant tournament).
- n=13 re-verified with S={1,2,3,4,5,6} (valid tournament): H=3,711,175, palindromic N confirmed.

### Impact

- THM-052 verification examples at n=13 and n=17 need correction: use valid circulant tournament generator sets, not QR mod p for p=1 mod 4.
- The earlier NONHAM vanishing verification at n=13 (nonham_vanish_n13_paley.py) was on a non-tournament. Needs re-running with valid tournament.
- No impact on the algebraic proof of THM-052.
- The memory entry "H(T_p) = 3, 189, 95095" is still correct (those are p=3 mod 4 primes).

---

## MISTAKE-012: Blue Pair / Tournament Complement Confusion

**Date discovered:** 2026-03-06 (opus-S11b continued²)
**Found by:** opus
**Affects:** blue_skeleton_even_r_synthesis.py (opus-S26)

### What was assumed

The script blue_skeleton_even_r_synthesis.py claimed that "BLUE PAIR HAS SAME TRANSFER MATRIX at odd n", conflating the tiling blue pair (flipping only non-path arcs) with the tournament complement (flipping ALL arcs).

### Why it was wrong

1. The tiling blue pair only flips arcs (a,b) with |a-b| >= 2, keeping path arcs (i,i-1) fixed. So s → -s only for non-path pairs, NOT all pairs. M(T') ≠ M(T) in general.

2. The tournament complement T^c (A^c[i][j] = 1-A[i][j]) does flip ALL arcs, giving s → -s everywhere. But M(T^c) ≠ M(T) in general either!

3. The correct formula: M(T^c) = diag(M(T)) - offdiag(M(T)). The diagonal is preserved and off-diagonal is NEGATED. This is because M[a,b] for a≠b involves n-2 edge weights (odd s-degree at odd n), while M[a,a] involves n-1 edge weights (even s-degree at odd n).

### The correct framing

**THM (verified n=5 exhaustive, n=7 spot):** For any tournament T at odd n:
  M(T^c)[a,a] = M(T)[a,a]   and   M(T^c)[a,b] = -M(T)[a,b] for a≠b.

COROLLARY: M(T^c) = M(T) iff M(T) is diagonal (scalar M at odd n).

### Impact

- The claim in blue_skeleton_even_r_synthesis.py is WRONG. M(T) ≠ M(T') for tiling blue pairs, and M(T) ≠ M(T^c) for tournament complement (unless M is diagonal).
- The analysis of which iso classes are "self-paired" in that script is unaffected (it's about iso class mapping, not about M equality).
- The complement formula is a NEW theorem that should be recorded.

---

## MISTAKE-013: "All VT tournaments are self-converse" assumption is FALSE

**Date discovered:** 2026-03-06 (kind-pasteur-S25e)
**Found by:** kind-pasteur (confirmed by independent background agent + McKay database)
**Affects:** opus-2026-03-06-S26 THM-052 extension claim (commit 8ed2516)

### What was assumed

opus-S26 claimed THM-052 extends to ALL vertex-transitive tournaments via the reflection-reversal bijection phi = sigma * rev, stating "for any vertex-transitive tournament T, there exists an anti-automorphism tau." This implicitly assumes all VT tournaments are self-converse.

### Why it was wrong

At n=21, ALL 22 non-circulant VT tournaments (from McKay's database) are NOT self-converse. These are Cayley tournaments on F_21 = Z/7 x| Z/3 (Frobenius group, smallest non-abelian group of odd order) with non-normal connection sets. Exhaustive backtracking search confirms no anti-automorphism exists for any of them. All 88 circulant tournaments at n=21 ARE self-converse.

The inversion map g -> g^{-1} gives an anti-automorphism iff S is a union of conjugacy classes ("normal" S). For non-abelian groups, non-normal S creates VT tournaments without any anti-automorphism.

### The correct framing

- THM-052 is PROVED for circulant tournaments (all abelian Cayley tournaments at odd n)
- THM-052 covers all VT tournaments at n <= 19 (all circulant by McKay data)
- THM-052 **FAILS** for non-abelian-group VT tournaments at n=21 (PROVED by computation)
- M[0,1] = 45,478,409 for the F_21 non-normal tournament (H = 123,522,430,238,361)
- N(0,1,j) is NOT palindromic: all 20 values N[j] != N[19-j]

### Impact

- opus's general VT theorem claim is FALSE and must be retracted
- THM-052 must specify "circulant" or "self-converse VT" as hypothesis
- Self-converse is the RIGHT boundary: SC VT => palindromic N => scalar M; non-SC VT => non-palindromic N => non-scalar M
- The per-orbit structure (N depends on vertex-pair orbit) still holds, but palindromicity requires SC

---

## MISTAKE-014: THM-052 extension to all VT tournaments is DISPROVED

**Date discovered:** 2026-03-06 (kind-pasteur-S25e)
**Found by:** kind-pasteur
**Affects:** THM-052 scope, opus-S26 proof

### What was assumed

THM-052 was claimed to hold for ALL vertex-transitive tournaments at odd n: M = (H/n)*I.

### Why it was wrong

Computational proof at n=21: the Cayley tournament on F_21 with non-normal connection set has M[0,1] = 45,478,409 != 0. The N(0,1,j) sequence is NOT palindromic:
  N[0] = 581,223,220,317 vs N[19] = 581,314,958,778 (differ by 91,738,461)
  Alternating sum = 45,478,409

### The correct framing

THM-052 holds for self-converse VT tournaments at odd n. Self-converse is necessary (not just sufficient) for palindromic N, and palindromic N at odd n is equivalent to scalar M for VT tournaments.

Hierarchy of scalar M:
1. Circulant (always SC) => scalar M [PROVED]
2. Abelian Cayley (always SC via -x) => scalar M [PROVED]
3. Non-abelian Cayley with normal S (SC via inversion) => scalar M [PROVED by same argument]
4. Non-abelian Cayley with non-normal S (NOT SC) => M NOT scalar [DISPROVED at n=21]

### Impact

- THM-052 scope must be restricted
- Opens question: which specific non-abelian VT tournaments have scalar M?
- Answer: exactly the self-converse ones (those with normal connection sets)

## MISTAKE-015: THM-055 coefficient table has wrong c_6 and c_0 at n=7

**Date discovered:** 2026-03-06 (opus-S29, independently confirmed by kind-pasteur-S25g)
**Found by:** opus + kind-pasteur
**Affects:** THM-055 coefficient table, c_0 formula

### What was assumed

THM-055's coefficient table at n=7 claimed:
- tr(c_6) = 720 = (n-1)! (universal)
- tr(c_0) = H - 6*bc - 3*t_5 + 249/4

### Why it was wrong

tr(c_{n-1}) = sum_P e_0(s_P) = sum_P 1 = n!, NOT (n-1)!. Direct polynomial fitting of W(r) = sum_P prod(r+s_i) confirms the leading coefficient is 5040 = 7! at n=7.

The verification script `trc2_exact_formula.py` actually produces max error c_0 = 67.5, NOT 0.0. The error was visible in the output (`c0= 1.8/ 69.2`) but was not caught. The root cause: c_0 was derived from `H - c2/4 - c4/16 - 720/64` using the WRONG c_6 value, hiding the error.

### The correct framing

The correct W(r) = tr(M(r)) coefficients at n=7 are:
- w_6 = 5040 = n! (universal)
- w_4 = 240*t_3 - 2100 (unchanged)
- w_2 = -60*t_3 + 12*t_5 + 24*bc + 231 (unchanged)
- w_0 = 2*t_3 - t_5 + 2*t_7 - 2*bc - 17/4

### Impact

- THM-055 coefficient table corrected: c_6 = n! and c_0 constant = -17/4 (not 249/4)
- The c_0 formula constant changed from +253/4 to -17/4 (difference = 270/4 = (5040-720)/64)
- The n=9 claim c_8 = 362880 = 9! = n! was already correct
- All middle coefficients (c_2, c_4) are unaffected

---

## MISTAKE-016: THM-059 recurrence had j^2 instead of (j+1)^2

**Date discovered:** 2026-03-07
**Found by:** opus-2026-03-07-S31
**Affects:** THM-059 central factorial recurrence statement (table and formulas were correct, only the recurrence formula was wrong)

### What was stated
b_{k,j} = b_{k-1,j-1} + j^2 * b_{k-1,j}

### Why it was wrong
Plugging in: b_{2,1} should be 5, but j^2 * b_{k-1,j} = 1^2 * 1 = 1, giving b_{2,1} = 1+1 = 2 (not 5).

### The correct formula
b_{k,j} = b_{k-1,j-1} + (j+1)^2 * b_{k-1,j}

This was confirmed by checking all 15 entries of the b-triangle for k=0..4. The correct recurrence is equivalent to the standard central factorial number recurrence with shifted column indices.

### Resolution
THM-059 corrected. The (j+1)^2 factor now has a combinatorial explanation via the Eulerian polynomial decomposition: F_f(r) = sum_k A(f+1,k) * (r+1/2)^{f-k} * (r-1/2)^k, where the central factorial structure emerges from expanding in u = (r+1/2)(r-1/2) = r^2-1/4.

### Impact
- The numerical table and all computed F_j values were always correct
- Only the stated recurrence formula was wrong
- The OEIS A036969 identification may need clarification (different column conventions)

---

## MISTAKE-013b: Missing 2^s Factor in M[a,b] Walsh Formula (THM-080)

**Date discovered:** 2026-03-07 (S35c7)
**Found by:** opus-2026-03-07-S35c7
**Affects:** THM-080 Walsh formula for M[a,b]

### What was assumed
The Walsh coefficient hat{M[a,b]}[S] = (-1)^{asc(S)} * (n-2-d)!/2^{n-2}, with NO dependence on the number of components. This was described as a "fundamental simplification" compared to H.

### Why it was wrong
The formula was verified exhaustively only at n=5, where ALL valid monomials have s=0 (zero unrooted components). At n=5, the maximum degree is 3, and with only 3 interior vertices, there's no room for unrooted even-length components to coexist with rooted ones. So the 2^s factor was always 1, making it invisible.

At n=7, degree-3 monomials like P1(a-rooted) + P2(unrooted) have s=1, and the formula without 2^s gives wrong reconstruction (16/20 failures).

### The correct formula
hat{M[a,b]}[S] = (-1)^{asc(S)} * **2^s** * (n-2-d)!/2^{n-2}

where s = number of unrooted (even-length) components. Each unrooted component has 2 valid orientations in the HP (both giving the same chi_S sign), contributing a factor of 2. Rooted components have only 1 valid orientation (pinned at a or b).

### Impact
- THM-080 formula corrected with 2^s factor
- Walsh proof of M[a,b]=M[b,a] symmetry still holds (2^s is symmetric in a,b)
- H-M comparison now shows PARALLEL structure: H has 2^r (all components unrooted), M has 2^s (only unrooted components contribute orientations)
- The "no r-dependence" claim was wrong; M DOES depend on component structure via s
- n=7 reconstruction: 20/20 match with corrected formula

### Lesson
Always verify formulas at the NEXT size up before claiming generality. n=5 was too small to expose the s-dependence.

---

## MISTAKE-017: "Non-Paley DRT at n=11" from invalid tournament connection set

**Date discovered:** 2026-03-07
**Found by:** kind-pasteur-2026-03-07-S39b
**Affects:** INV-068, MEMORY.md DRT analysis section, TANGENTS.md DRT entry

### What was assumed
A "non-Paley DRT at n=11" was constructed using connection set {1,2,3,5,8} in Z_11 (circulant digraph). Claims: c3=44, c5=407, H=69311, |Aut|=11. "Paley strictly dominates in ALL cycle counts."

### Why it was wrong
The connection set {1,2,3,5,8} does NOT give a tournament. For a circulant tournament on Z_p, the connection set S must satisfy S ∩ (-S) = ∅ (so each pair {i,j} has exactly one directed arc). But {1,2,3,5,8} contains BOTH 3 and 8=11-3, and BOTH 1 and 10=11-1... wait, 10 is NOT in S. Let me re-check: -S = {11-s : s ∈ S} = {10, 9, 8, 6, 3}. S ∩ (-S) = {3, 8} ≠ ∅.

So for any pair (i,j) where (j-i)%11 ∈ {3, 8}: BOTH T[i][j]=1 AND T[j][i]=1. The resulting digraph has bidirectional edges and is NOT a tournament. All computations (c3, c5, H, is_doubly_regular) were performed on a non-tournament digraph and are MEANINGLESS.

### The correct framing
An exhaustive search of all 32 valid tournament connection sets in Z_11 (choosing one from each pair (d, 11-d)) found exactly 2 that are (11,5,2)-difference sets: {1,3,4,5,9} (QR) and {2,6,7,8,10} (NQR). These give ISOMORPHIC tournaments (both Paley T_11). There is NO non-Paley circulant DRT at n=11.

Whether a non-circulant DRT exists at n=11 remains an open question. At prime order p, all groups are Z_p, so all Cayley tournaments are circulant. A non-circulant DRT would need a different construction.

### Impact
- ALL claims about "non-Paley DRT at n=11" are INVALID
- INV-068 "Paley dominance" finding needs complete re-evaluation
- The claimed c3=44 was wrong — Moon's formula gives c3=55 for ALL regular n=11 tournaments
- The claimed "Paley strictly dominates in all cycle counts" is unverifiable since no valid comparison tournament exists
- MEMORY.md entry on DRT analysis at n=11 needs correction

### Lesson
When constructing a circulant tournament from a connection set S ⊂ Z_p^*, ALWAYS verify S ∩ (-S mod p) = ∅. A (v,k,λ)-difference set is NOT automatically a valid tournament connection set.

---

## MISTAKE-016b: Wrong formula for ker(d_2^rel) in relative homology

**Date discovered:** 2026-03-08 (kind-pasteur-S41)
**Found by:** kind-pasteur-S41, via manual computation contradicting script output
**Affects:** beta2_relative_homology.py, beta2_relative_correct.py; HYP-213 verification

### What was assumed
The script `beta2_relative_homology.py` computed ker(∂_2^rel) as:
  `(ker ∂_2 + V_2) / V_2`
where V_2 = Ω_2(T\v) (non-v subspace of Ω_2).

### Why it was wrong
The correct formula for ker(∂_2^rel) in the quotient complex Ω_*/V_* is:
  `∂_2^{-1}(V_1) / V_2`
where V_1 = Ω_1(T\v) (non-v arcs). This is the preimage of V_1 modulo V_2.

The wrong formula misses elements x ∈ Ω_2 whose boundary ∂_2(x) is NONZERO but lies entirely in V_1 (non-v arcs). Such elements are relative 2-cycles but NOT absolute 2-cycles.

Concretely: P_v ∘ ∂_2(x) = 0 (projection of boundary onto v-arcs vanishes), but ∂_2(x) ≠ 0.

### The correct framing
ker(∂_2^rel) = dim(Ω_2) - rk(M) - dim(V_2), where M = ∂_2|_{Ω_2} restricted to rows of v-arcs.

This correctly counts the preimage of V_1 in Ω_2.

### Impact
- **HYP-213 is REFUTED**: H_2(T, T\v) > 0 for many (T,v) pairs at n ≥ 4.
  - n=4: 16/256 pairs (6.25%)
  - n=5: 840/5120 pairs (16.4%)
  - n=6: 35,328/196,608 pairs (18%)
- The proposed inductive proof of β_2 = 0 via H_2(T,T\v) = 0 does NOT work.
- However, β_2 = 0 itself is NOT affected — it remains computationally verified.
- The connecting map δ: H_2(T,T\v) → H_1(T\v) is always injective (verified n=4,5), consistent with β_2 = 0 via the long exact sequence.

### Lesson
When computing relative homology H_*(X, A) via quotient complexes:
1. ker(∂_p^{rel}) is NOT (ker ∂_p + C_*(A)) / C_*(A).
2. ker(∂_p^{rel}) = ∂_p^{-1}(C_{p-1}(A)) / C_p(A).
3. These differ whenever there are elements whose boundary is nonzero but lands in the sub-complex.
4. Always verify relative homology against the long exact sequence.

---

## MISTAKE-018: beta_3 <= 1 Assumed for All Tournaments

**Date discovered:** 2026-03-09 (kind-pasteur-S48)
**Found by:** kind-pasteur-S48 via extended sampling at n=8 (5000 random tournaments)
**Affects:** THM-123 (was THM-110) proof architecture, HYP-371b, HYP-375, HYP-342, HYP-380, HYP-393 scope

### What was assumed
Multiple hypotheses and proof strategies assumed beta_3 <= 1 for ALL tournaments:
- HYP-371b: "beta_3=2 impossible"
- HYP-375: "beta_3 <= 1 at n=9"
- THM-123 proof architecture: Claims I, II, III designed to prove beta_3 <= 1
- The opus exhaustive proof at n=7 was incorrectly assumed to generalize

### Why it was wrong
beta_3 = 2 DOES occur at n=8. Four examples found in 5000 random tournaments (rate ~0.08%):
- Profile: (1, 0, 0, 2, 0, 0, 0, 0) — two independent H_3 generators
- Scores: (2,3,3,3,4,4,4,5) and (3,3,3,3,4,4,4,4) — near-regular
- Confirmed by BOTH max_p=5 and max_p=7 in full_chain_complex_modp (mod-p exact)
- All b3=2 tournaments have good vertices (b3(T\v)=0 for some v)

Previous sampling (200 at n=9, 100 at n=8) was insufficient to detect 0.08% rate.

### The correct framing
- beta_3 <= 1 is proved ONLY at n <= 7 (exhaustive, HYP-393)
- beta_3 = 2 at n=8 (confirmed, 4/5000)
- beta_3 may grow further at n >= 9

### Impact
- THM-123 proof architecture is valid ONLY at n <= 7
- Claims I (i_*-injectivity) also FAIL at n=8 (13 violations in 5000 trials, even with b4=0)
- Claim III (consecutive seesaw) FAILS at n=8 (beta_3+beta_4 coexistence)
- The beta_3 <= 1 bound is a SMALL-n PHENOMENON, not a universal property

### Lesson
1. Small sample sizes (100-200) cannot detect 0.1% phenomena. Use 5000+ for rare events.
2. Properties proved exhaustively at n<=7 do NOT automatically extend to n>=8.
3. n=8 is a critical threshold where many path homology structural properties break down:
   consecutive seesaw, i_*-injectivity, beta_3<=1, bad vertex acyclicity.

## MISTAKE-019: Int64 Overflow in Chained Numpy Matrix Multiplication

**Date discovered:** 2026-03-10 (kind-pasteur-S50)
**Found by:** kind-pasteur-S50 via comparison of two K_tv computation methods
**Affects:** opus-S59's tv_cycle_structure.py (Ghost Cycle "failures" are spurious), potentially any script using `A @ B @ C % PRIME` pattern

### What was assumed
`D3 @ (tv_omega @ ob3).T % PRIME` safely computes the matrix product mod PRIME.

### Why it was wrong
With `RANK_PRIME = 2^31 - 1`:
- `tv_omega @ ob3` produces entries up to ~4.6 × 10^18 (fits in int64, max 9.2 × 10^18)
- BUT these are NOT reduced mod PRIME — entries can be >> PRIME
- `D3 @ X.T` then involves products up to `(2^31) * (4.6e18) = 9.9e27`, massively exceeding int64 max
- Result: silent int64 overflow → wrong matrix entries → wrong rank → wrong K_tv

This caused opus-S59's tv_cycle_structure.py to report Ghost Cycle failures in 14/504 pairs at n=7 and 11/304 at n=8. ALL of these "failures" are arithmetic artifacts.

### The correct framing
ALWAYS reduce mod PRIME between chained matrix multiplications:
```python
# WRONG (can overflow):
result = A @ B @ C % PRIME

# RIGHT:
temp = A @ B % PRIME
result = temp @ C % PRIME

# BEST: use the new safe utility:
from tournament_utils import matmul_mod
result = matmul_mod(matmul_mod(A, B), C)
```

The `matmul_mod()` function in tournament_utils.py automatically chunks the inner dimension to prevent overflow, even for single multiplications with large entries.

### Impact
- Ghost Cycle (K_tv = B_tv) HOLDS universally at n ≤ 8 (0 real failures in 1000+ tests)
- HYP-408 (codim-1 universality) remains computationally verified at n ≤ 8
- No real mathematical result is affected; only the false "counterexamples" are invalidated

### Lesson
1. NEVER chain numpy `@` without intermediate `% PRIME` when PRIME ≈ 2^31
2. Use `matmul_mod()` from tournament_utils.py for all modular matrix arithmetic
3. When two equivalent computations disagree, suspect numerical issues before mathematical failure

## MISTAKE-020: Truncated Chain Complex Gives False Betti Numbers at Top Degree

**Date discovered:** 2026-03-10
**Found by:** kind-pasteur-S50

### What was assumed
Using `full_chain_complex_modp(A, n, max_p=6)` for n=8 tournaments, opus-S59 reported β_6 nonzero for 89.8% of tournaments (HYP-420), with values ranging 1-25.

### Why it was wrong
With `max_p=6`, the computation gives β_6 = ker(d_6) - ranks.get(7, 0). Since degree 7 is not computed, `ranks.get(7, 0)` returns 0. The reported "β_6" is actually just dim(ker d_6), NOT the true Betti number.

With `max_p=7` (full complex): d_7 is injective on Omega_7, and rk(d_7) = ker(d_6) EXACTLY for all 50 tested tournaments. True β_6 = 0 always.

### The correct framing
The Betti number at the highest computed degree is always an UPPER BOUND (missing the image from the next degree). For n-vertex tournaments, always use `max_p=n-1` to get correct Betti numbers, especially at degrees n-2 and n-1.

Correct results: β_{n-1} = β_{n-2} = 0 for ALL tournaments at n=3-8 (HYP-423, HYP-424). The top boundary map d_{n-1} is always injective.

### Impact
- HYP-420 is FALSE. β_{n-2} is NOT generically nonzero at n=8.
- The "β_6 among β_3=1" distribution in opus's beta4_at_n7.out is entirely artifactual.
- All lower-degree Betti numbers (β_0 through β_5) from that computation are correct.

### Lesson
1. ALWAYS use max_p=n-1 when computing Betti numbers to avoid truncation artifacts
2. Betti at the max computed degree is an UPPER BOUND (Betti at max_deg-1 and below are exact)
3. When β at max_deg seems surprisingly large/nonzero, check if im(d_{max_deg+1}) is missing

---

## MISTAKE-019b: THM-136 Sign Convention Error

**Date discovered:** 2026-03-12
**Found by:** kind-pasteur-S57
**Affects:** THM-136 formula statement (not the verbal description or proof mechanism)

### What was assumed
The trace alternation sign formula was stated as:
`sign(tr(A_P^k) - tr(A_I^k)) = (-1)^{(k-3)/2}`

### Why it was wrong
Direct computation at p=7,11,19,23 shows:
- k=5: Delta > 0 (positive), but (-1)^{(5-3)/2} = (-1)^1 = -1 (WRONG)
- k=7: Delta < 0 (negative), but (-1)^{(7-3)/2} = (-1)^2 = +1 (WRONG)

The formula gives the OPPOSITE sign for every k.

### The correct framing
`sign(tr(A_P^k) - tr(A_I^k)) = (-1)^{(k-1)/2}`

Equivalently: positive for k = 1 mod 4, negative for k = 3 mod 4.
Verified with 1218+ individual (k,p) tests, zero failures.

Note: the VERBAL description in THM-136 was always correct ("k=1 mod 4: Paley wins").
Only the symbolic formula was off by one power.

### Impact
- Formula in THM-136 theorem file CORRECTED
- No downstream impact: all proofs used the verbal description, not the formula
- The algebraic proof (kind-pasteur-S57) uses the correct convention throughout

## MISTAKE-021: S70 "GLMY Betti Numbers" Use Wrong Chain Complex

**Date discovered:** 2026-03-13
**Found by:** opus-2026-03-13-S71
**Affects:** ALL scripts from S70 session: betti_omega_connection.py, betti_divisibility.py, per_eigenspace_betti.py, per_eig_betti_n9.py, and all results/theorems derived from them (THM-154, eigenspace Betti uniformity)

### What was assumed
The S70 scripts computed "GLMY path homology Betti numbers" using:
- Allowed paths = "regular paths" (v_i→v_{i+1} AND v_{i-1}→v_{i+1})
- Boundary = interior-only deletion (indices 1 to m-1)
Results were called "GLMY Betti numbers" and compared with GLMY literature.

### Why it was wrong
The actual GLMY path homology uses:
- Allowed paths = directed paths (v_i→v_{i+1} only, NO skip-one requirement)
- Boundary = full vertex deletion (indices 0 to m), but restricted to Ω_m subspace
- Ω_m = {u ∈ A_m : ∂u has all components in A_{m-1}}

**These give DIFFERENT chain complexes with different Betti numbers:**
- Paley P_7 GLMY: β = [1,0,0,0,6,0,0], dim(A_2)=63
- Paley P_7 S70:  β = [7,0,0,21,21,21,21], dim(A_2)=21

The "regular path + interior boundary" complex IS a valid chain complex
(d²=0 verified), but it is NOT standard GLMY path homology.

### The correct framing
There are TWO distinct valid chain complexes for tournaments:

1. **GLMY Path Homology** (standard): directed paths, full boundary on Ω_m.
   Implemented correctly in path_homology_v2.py.
   β_0 = 1 for all tournaments. β_2 = 0 for all tested tournaments (n≤8).

2. **Tournament Regular Homology (TRH)** (novel?): regular paths, interior boundary.
   Used in S70 scripts. β_0 = n for all tournaments on n vertices.
   Has eigenspace Betti uniformity and divisibility by n for circulants.

Both are valid mathematical objects. But they should not be conflated.

### Impact
- THM-154 (Betti divisibility) applies to TRH, not GLMY
- Eigenspace Betti uniformity applies to TRH, not GLMY
- β_2=0 for all tournaments holds for BOTH (GLMY verified n≤8, TRH verified n≤8)
- The S70 "per-eigenspace Betti" results are self-consistent but not GLMY
- The S38-S41 β_2=0 results (from path_homology_v2.py) are correct GLMY
- circulant_homology.py implements yet another convention (full boundary on regular paths) which is NEITHER GLMY nor TRH

### Lesson
ALWAYS verify which chain complex you're computing. The three ingredients
(allowed paths, boundary convention, Ω subspace) must be consistent.
When reading "path homology" results, check which convention is used.

---

## MISTAKE-019c: TWO bugs in independent set backtracking algorithm

**Date discovered:** 2026-03-13, kind-pasteur-S60
**Found by:** kind-pasteur
**Affects:** alpha3_p7_only.py, alpha3_moment_analysis.py, overlap_weight_analysis.py, H_energy_decomposition.py, cycle_walsh_decomposition.py, moment_cancellation_mechanism.py, overlap_gauss_bridge.py, alpha_directed_p11.py, alpha_full_p11.py, alpha2_direct_verify.py, backtrack_debug.py (ALL files with independent set enumeration)

### Bug 1: Missing vertex 0 (`backtrack(0,0,0)` should be `backtrack(-1,0,0)`)

The backtracking function `backtrack(v, mask, size)` iterates `for w in range(v+1, n)`. When called with `v=0`, the loop starts at `w=1`, SKIPPING vertex 0 entirely. This undercounts all alpha_j.

**Fix:** Call `backtrack(-1, 0, 0)` so the loop starts at `w=0`.

### Bug 2: Skipping consecutive indices (`backtrack(w+1, ...)` should be `backtrack(w, ...)`)

The recursive call `backtrack(w + 1, mask | nbr[w], size + 1)` passes `v = w+1`. At the next level, the loop starts at `w' = v+1 = w+2`, SKIPPING index `w+1`. This means any independent set containing cycles with consecutive indices is missed.

**Fix:** Change to `backtrack(w, mask | nbr[w], size + 1)`. Then the next level's loop starts at `w+1`, correctly considering all higher indices.

### Concrete example

At p=7, Interval tournament S=[1,2,3]:
- 59 directed cycles, 14 disjoint (3,3)-pairs (correct)
- Bug 2: Pair (5,6) = ({0,3,6}, {1,2,5}) has consecutive indices and was SKIPPED
- Backtracking gave alpha_2=13 instead of 14, H=171 instead of 175
- Held-Karp gives H=175 (correct)

### Impact
- All previous alpha_j values from backtracking are SUSPECT
- THM-027 alpha_2 values at p=7 need recheck (Paley alpha_2=7 was coincidentally correct because no consecutive disjoint pairs)
- Any H derived from backtracking alpha may be wrong

### Lesson
When implementing independent set enumeration via backtracking, the recursive call after selecting vertex w should pass v=w (NOT v=w+1). The `range(v+1, n)` in the next level already excludes w.

---

## MISTAKE-022: Sparse Gaussian Elimination Fill-In Bug

**Date discovered:** 2026-03-13, opus-S71c (9th context window)
**Found by:** opus, when k=0 eigenspace Betti numbers came out negative
**Affects:** p19_omega5_sparse.py, p23_omega5_sparse.py, p31_omega5_sparse.py, p43_omega5_sparse.py (ALL scripts using the sparse Gaussian pattern with single-pass row iteration)

### What was assumed
The sparse Gaussian elimination iterated over `sorted(col.keys())` once, subtracting each matching pivot. This should correctly eliminate all pivot contributions.

### Why it was wrong
When subtracting a pivot at row `r`, the pivot vector has entries at rows `r' > r` (fill-in). These new entries at rows NOT in the original column are never checked against existing pivots at those rows, because the `sorted(col.keys())` list was computed BEFORE the subtraction and doesn't include fill-in entries.

Concrete example: column has entries at rows {3, 7}. Pivot at row 3 has entries at {3, 5, 7}. After subtracting pivot 3, the column has entries at {5, 7}. Row 5 was NOT in the original sorted list, so even if there's a pivot at row 5, it's never subtracted. This causes the rank to be OVERCOUNTED (some columns that should reduce to zero don't).

### The correct framing
After any pivot subtraction, restart the row scan from the beginning (or at least from the newly-created entry). A simple fix: wrap the elimination loop in `while changed: ... break after subtraction`.

### Impact
- **P_19 Omega_5 was 12602 (WRONG), correct is 23832**
- **P_23 Omega_5 was 50715 (WRONG), correct is 78430**
- **P_31 Omega_5 was 252065 (WRONG), correct is 456330**
- **P_43 Omega_5 was 1429652 (WRONG), correct is 2865660**
- P_7 and P_11 were unaffected (small enough that fill-in didn't change rank)
- HYP-790 ("Omega_5 not polynomial in m") was based on wrong data — **RETRACTED**
- **CORRECTED**: Ω_5 = m(m-1)(m³-6m²+10m-2) — a **clean integer polynomial** in m!
- All formulas Ω_d for d ≤ 5 are now proven/verified

### Lesson
In sparse Gaussian elimination, fill-in from pivot subtraction can create new entries at rows that were not in the original column. These MUST be processed against their pivots. Always use a while loop that restarts after each subtraction, or maintain a priority queue of unprocessed rows.

---

## MISTAKE-023: α₁ Counts DIRECTED Odd Cycles, Not Vertex-Sets

**Date discovered:** 2026-03-14
**Found by:** opus-2026-03-14-S71d
**Affects:** two_and_three_universality.py, i3_mod3_proof.py, vandermonde_sigma_connection.py, jacobsthal_23_deep.py (first version), and any script computing I(CG, x) by counting cycle vertex-sets

### What was assumed

The independence polynomial I(Ω, x) was computed by enumerating odd cycle **vertex-sets** (frozenset of vertices), counting each set once regardless of how many distinct directed cycles it supports.

### Why it was wrong

The conflict graph Ω(T) has vertices = **directed odd cycles** (definition in definitions.md line 37). For 3-cycles in tournaments, each vertex triple supports at most 1 directed 3-cycle, so vertex-set counting is correct. But for 5-cycles and above, a single vertex-set can support **multiple** distinct directed cycles:

- Example: bits=40 at n=5, the 5-vertex set {0,1,2,3,4} supports **3** distinct directed 5-cycles
- Vertex-set method gives α₁=5 → I(2)=11, but H=15
- Directed-cycle method gives α₁=7 → I(2)=15 = H ✓

### The correct framing

When computing I(Ω, x):
1. For each vertex-set of size k, enumerate ALL distinct directed k-cycles (normalize by fixing start vertex and direction)
2. Each distinct directed cycle is a SEPARATE vertex of Ω(T)
3. Two vertices of Ω are adjacent iff the underlying vertex-sets intersect

**Exhaustive verification at n=5:**
- Vertex-set method: 184/1024 mismatches with H
- Directed-cycle method: 0/1024 mismatches with H

### Impact

- All α₁ values from scripts using vertex-set counting at n≥5 are WRONG (undercounted)
- The Vandermonde extraction results (HYP-867, HYP-868) were based on the wrong α values
- The 3/2 ratio result may still hold (it was measured within lambda fibers, not from α directly)
- The structural insights about 7→8 transition are UNAFFECTED (vertex-set counting is correct for α₂ when cycles have different sizes)
- Scripts need to be updated to use directed cycle enumeration

### Lesson

The definition says "vertices are **directed** odd cycles." For 3-cycles in tournaments, vertex-set = directed cycle (1-to-1). For k≥5 cycles, a k-vertex tournament subtournament can have multiple Hamiltonian cycles. Always enumerate directed cycles explicitly.

---

## MISTAKE-024: H=63 Falsely Claimed Permanently Forbidden

**Date discovered:** 2026-03-14
**Found by:** opus-2026-03-14-S71h (cross-referencing S86 broadcast with prior results)
**Affects:** HYP-1303, MSG-218 (S86 broadcast), MSG-139 (S86 to kind-pasteur)

### What was assumed

S86 claimed: "H=63 FORBIDDEN: 63=7×9=I(K₃,2)×I(2K₁,2). Requires K₃ component in Ω, impossible by THM-201." This was marked CONFIRMED as HYP-1303.

### Why it was wrong

The argument only blocks DISCONNECTED conflict graphs where Ω = K₃ ⊔ 2K₁. But Ω can be a CONNECTED graph with I(Ω,2)=63. Multiple prior sessions had already established:
- S65-c (MSG-084): "H=63,107,119,149 (the n=7-specific gaps) are ALL achieved at n=8"
- S71f (MSG-197): "63 achievable at n≥8"
- S71g (MSG-201): "H=63 found at n=8 (27/100k)"
- hspectrum_density.out: "63 = 7*9 -- ACHIEVED at n=8 (not permanent)"

The S86 session re-derived an incomplete argument without checking these earlier verified results.

### The correct framing

H=63 is NOT permanently forbidden. It is a temporary gap at n=7 (like 107, 119, 149) that IS achieved at n=8. The ONLY permanently forbidden H values proved so far are {7, 21}.

The disconnected decomposition 63=7×9 is correctly blocked, but connected Ω graphs with I=63 exist and can apparently be realized as Ω(T) at n=8.

### Impact

- HYP-1303 changed to REFUTED
- HYP-1295 changed to REFUTED
- The S86 broadcast claim "all three known forbidden values {7,21,63} are now explained" is WRONG
- Only {7, 21} are explained as permanently forbidden

### Lesson

Always check prior session results before claiming a new proof. The H-spectrum density analysis (hspectrum_density.out) had already settled this question computationally. Cross-reference before broadcasting claims.

---

## MISTAKE-025: S112 W(8) Value Off By 8

**Date discovered:** 2026-03-16 (opus-S90 continuation session)
**Found by:** opus-S90, via independent brute-force verification
**Affects:** kind-pasteur-S112 W(n) sequence, D_n(2) computations

### What was claimed
kind-pasteur-S112 reported W(n) = 1, 2, 8, 32, 158, 928, 6350, **49760**, 439766 for n=1..9.

### Why it was wrong
Independent brute-force computation (iterating over all 8! = 40320 permutations) gives W(8) = **49752**, not 49760. The error is exactly +8 = 2³. The S89c C-program DP computation also gives 49752, confirming the brute-force result.

### The correct values
W(n) = 1, 2, 8, 32, 158, 928, 6350, **49752**, 439670 (from S89c DP).

### Impact assessment
- **S89c values (opus):** CORRECT through n=27 (computed by bitmask DP in C).
- **S112 values (kind-pasteur):** INCORRECT at n=8 by +8. Values at n≤7 match. Values at n≥9 need reverification against S89c.
- **OEIS submission:** Use S89c values, not S112.
- **CV² formula and g_k polynomials:** UNAFFECTED (derived independently of W(n) enumeration).

### Source
opus-2026-03-16 (S90 continuation): brute-force W(8) verification via Python permutation enumeration.

### Lesson
When two independent computations disagree, trust the one with the simpler algorithm (brute force) over the one with more complex logic. The discrepancy of exactly 2³ suggests a boundary condition or off-by-one error in the S112 computation, not a random bug.

---

## MISTAKE-026: Cross-Ratio of Cayley Orbit Initially Claimed as 8/7

**Date discovered:** 2026-03-15 (code review during S90)
**Found by:** opus-S90 code review agent
**Affects:** monad_cayley_s90c.py, commit messages

### What was claimed
The cross-ratio of the Q-orbit of x=2 was initially computed as 8/7, using the WRONG orbit point (3 instead of -3).

### Why it was wrong
Q(2) = (1+2)/(1-2) = 3/(-1) = **-3**, not 3. The orbit is {2, **-3**, -1/2, 1/3}. The cross-ratio CR(2, -3, -1/2, 1/3) = **2**, not 8/7.

### The correct value
Cross-ratio = 2 = the OCF fugacity itself. This is MORE meaningful than 8/7 — the cross-ratio equals the evaluation point.

### Impact
- The narrative about "tournament constant 8/7" in commit messages is wrong.
- The correct "tournament constant" is 2 (the fugacity).
- Script monad_cayley_s90c.py has been corrected.

### Source
Code review agent, opus-2026-03-15-S90.

---

## MISTAKE-018b: THM-225 "Universal Top Eigenvalue = n" is FALSE at n ≥ 9

**NOTE:** This was originally numbered MISTAKE-018 from a different branch, causing a collision with MISTAKE-018 (beta_3 <= 1). Renumbered to 018b by opus-2026-04-01-S1.

**Date discovered:** 2026-03-15
**Found by:** opus-S72d
**Affects:** THM-225, HYP-1594

### What was assumed

That the top eigenvalue of C_T^TC_T equals n for ALL tournaments on n vertices (verified exhaustively at n=5, sampled at n=6). This was stated as a theorem.

### Why it was wrong

The proof strategy required rank(C_R) < r = (n-1)(n-2)/2. At n ≤ 8, this holds because max c₃ < r (the number of cyclic triples never exceeds the rank of the full constraint matrix). At n = 9, c₃ can reach 30 while r = 28, and for ~0.1% of tournaments, rank(C_R) achieves its maximum r, leaving ker(C_R) ∩ im(C^T) = {0}. The top eigenvalue then drops to ~8.84-8.94.

### The correct framing

THM-225 holds for n ≤ 8 (PROVED for n ≤ 5 exhaustive, sampled at n=6,7,8 with 0 violations from 20000 samples each). It FAILS at n ≥ 9. The condition for top eigenvalue = n is rank(C_R) < (n-1)(n-2)/2.

### Impact

The spectral T/R duality (C_T^TC_T + C_R^TC_R = n·P) and the 3/n bridge framework remain valid. The universal top eigenvalue was a COROLLARY that holds only when the cyclic boundaries don't span the full constraint space.

### Lesson

When verifying at n=5 and n=6, the parameter regime (c₃ < r always) hid the failure mode. Always check at the CROSSOVER point where qualitative behavior changes. For rank arguments, the critical n is where max c₃ first exceeds r.

---

## MISTAKE-027: THM-080 Amplitude Table Wrong at n=9

**Date discovered:** 2026-03-16
**Found by:** opus-2026-03-16-S73
**Affects:** THM-080 amplitude table (lines 156-161), not the formula itself

### What was assumed

The amplitude table in THM-080 listed n=9 entries as: deg 1 (s=0) = 3/2, deg 3 (s=0) = 3/8, deg 3 (s=1) = 3/4, deg 5 (s=0) = 1/16, deg 7 = 1/128.

### Why it was wrong

The stated formula is |hat{M}[S]| = 2^s × (n-2-|S|)!/2^{n-2}. At n=9:
- d=1, s=0: formula gives 6!/128 = **45/8**, not 3/2
- d=3, s=0: formula gives 4!/128 = **3/16**, not 3/8
- d=3, s=1: formula gives 2×3/16 = **3/8**, not 3/4
- d=5, s=0: formula gives 2!/128 = **1/64**, not 1/16
- d=7, s=0: formula gives 0!/128 = 1/128 ✓ (only correct entry)

The formula works perfectly at n=3, 5, 7 — only n=9 has errors. The d=3 and d=5 wrong values are the CORRECT formula values but with s shifted up by 1 or 2 (unrooted component miscount). The d=1 entry (3/2) doesn't correspond to any valid s value (45/8 × 2^s ≠ 3/2 for any integer s).

### The correct framing

The formula is correct. The table had a transcription error at n=9 only. The n=9 verification was "partial" and apparently didn't catch the table/formula mismatch. Corrected amplitude table is in THM-080.

### Impact

Low — the formula itself is correct and was verified computationally at n=5 (exhaustive), n=6 (exhaustive), n=7 (20/20). Only the summary table for n=9 was wrong. No downstream results depend on the specific n=9 table values.

### Lesson

**This is MISTAKE-013b (the original missing 2^s) echoing forward.** The 2^s correction was caught and fixed at n=7, but the n=9 table values were apparently populated from a pre-correction computation (or from hand calculation that repeated the component-counting error at higher n). Always re-derive table entries from the corrected formula rather than carrying forward values from partial computation.

This is also a meta-lesson about the amplitude table itself: it was the only place in THM-080 where specific numerical values were stated without being individually verified against the formula. The formula (analytically proved) was trustworthy; the table (hand-calculated) was not.

---

## MISTAKE-028: Mersenne / k-nacci Numbers Falsely Claimed to Control Forbidden H Values

**Date discovered:** 2026-03-17
**Found by:** opus-2026-03-17-S74 (forbidden values audit)
**Affects:** casual-writeup.md, formal-writeup.md, substack-hooks.md (Hook N), HYP-1600, HYP-1618 (original), HYP-1623, HYP-1624, riemann-zeta-tournament.md, multiple results files

### What was assumed

Multiple writeups and hypotheses claimed: "The k-nacci Mersenne identity connects forbidden H values (7 = 2^3 - 1, 31 = 2^5 - 1, 127 = 2^7 - 1) to Mersenne primes via k-nacci transfer matrices." The original HYP-1618 claimed "ζ(-3) = 7" (standard Riemann zeta). Various scripts called 31 "forbidden."

### Why it was wrong

1. **31 is achievable** at n=6 (tournament bits=146, verified exhaustively).
2. **63 is achievable** at n=8 (already documented in MISTAKE-024).
3. **127 is achievable** at n=7.
4. The standard Riemann ζ(-3) = 1/120, NOT 7.
5. The tribonacci trace sequence [1, 3, 7, 11, 21, 39, 71, 131, ...] contains both forbidden values (7, 21) AND achievable values (1, 3, 11, 39, 71, 131, ...). The k-nacci trace hitting 7 and 21 is a numerical coincidence, not a causal mechanism.

### The correct framing

**Only H=7 and H=21 are permanently forbidden** (proved: THM-029, THM-079). The actual mechanisms are:
- H=7: 3 pairwise-conflicting cycles always force additional cycles (THM-029)
- H=21: all OCF decompositions blocked by tournament forcing (THM-079, 464-line proof)

Best characterizations of {7, 21}:
- {7·3⁰, 7·3¹}: the 7-obstruction has nilpotency 2 (HYP-1231). 7·3² = 63 is achievable.
- {Φ₃(2), Φ₃(4)}: third cyclotomic polynomial at even args (HYP-1317). Φ₃(6) = 43 is achievable.
- Both have I-polynomials factoring through I(K₃, x) = (1+3x) (HYP-1315).

THM-227 (k-nacci Mersenne) is a valid theorem about transfer matrices. It just doesn't characterize forbidden H values.

### Impact

Medium — the false claim propagated through 6+ files across multiple sessions. All have been corrected. No theorems or proofs are affected (the actual forbidden value proofs THM-029 and THM-079 are correct and don't use the Mersenne connection).

### Lesson

**One data point is not a pattern.** The entire false extrapolation rested on a single observation: 7 = 2³ - 1 is both a Mersenne number and forbidden. From this, it was incorrectly inferred that other Mersenne numbers (31, 127) are also forbidden. A simple check (is H=31 achieved at n=6?) would have caught this immediately.

This is a variant of MISTAKE-024 (H=63 falsely claimed forbidden) — the same class of error, just with a different numerological motivation. The meta-lesson: when claiming a numerical pattern "explains" something, verify it at the NEXT case before asserting it as a principle.

---

## MISTAKE-029: Formula E = (T - D)/2 is WRONG for the meta-graph edge count

**Date discovered:** 2026-03-23
**Found by:** opus-2026-03-23-S211
**Affects:** degeneracy_second_moment_s210.py, all claims that E = (3T - S_2)/4

### What was assumed
The meta-graph edge count formula E = (T - D)/2 was claimed in S210, where T = total arc-orbits and D = sum C(mult(C→D), 2) is the degeneracy. The derived formula E = (3T - S_2)/4 was presented as the "keystone" edge count formula, and the reverse-engineered S_2 sequence was reported as new.

### Why it was wrong
The formula ignores **self-loop orbits** and **directed multi-edge excess**. The correct decomposition is:

T = SL + 2E + excess_cross

where:
- SL = sum_C mult(C→C) = total self-loop arc-orbits
- excess_cross = sum_{{C,D}} (mult(C→D) + mult(D→C) - 2) for connected C≠D

So: **E = (T - SL - excess_cross) / 2**, NOT (T - D) / 2.

At n=3,4,5: the formula happened to give correct answers because SL + excess = D exactly (coincidence). At n=6: SL + excess = 58 + 66 = 124, but D = 122, so E_wrong = 291 ≠ E_actual = 290.

The "reverse-engineered" S_2 values for n≥6 were derived from this wrong formula and are therefore incorrect. The actual S_2 at n=6 is 948, not 952.

### The correct framing
E(n) must be computed directly from the meta-graph adjacency (F matrix), not from aggregate orbit statistics. There is no known closed-form expression for E(n) in terms of T, D, or S_2. The quantities T(n) and S_2(n) give orbit-level statistics but cannot determine which pairs of classes are actually adjacent.

### Impact
- The "keystone formula" E = (3T - S_2)/4 is invalid
- The S_2 sequence 8, 28, 144, 952, 10392, 200220, 7018596 from S210 is wrong at n≥6
- The correct S_2 at n=6 is 948 (from direct orbit computation)
- The gap sequence G(n) = T - 2E = 2, 6, 28, 124, 740, 5966, 85698 IS correct and novel
- All independently computed E values (E(3..8) = 1, 5, 30, 290, 4086, 91161) remain valid

### Lesson
When a formula passes at small n, always verify at the next n where new phenomena emerge. At n≤5, every class has SL + excess = D (a coincidence), so the formula appeared correct. At n=6, the coincidence breaks. **Integer division can mask off-by-one errors**: at n=3, (T-D)/2 = 3/2 = 1.5, which rounded to 1 via `//` and happened to match E=1.

---

## MISTAKE-030: "SL_orbits" is a misnomer — it includes multi-edge orbits, not just self-loops

**Date discovered:** 2026-03-23
**Found by:** Devil's advocate audit (opus-2026-03-23-S246), confirmed by opus-S245
**Affects:** burnside_edge_verify_s242.py, recursive_sl_s244.py, all scripts using "SL_orbits"

### What was assumed
The quantity "SL_orbits" = edge_orbits - E(G_n) was treated as counting self-loop edge orbits (orbits where both endpoints are in the same iso class).

### Why it was wrong
"SL_orbits" actually counts ALL non-simple-edge orbits: self-loop orbits PLUS multi-edge orbits (additional orbits connecting already-connected class pairs). At n=5: true self-loop edge orbits = 14, but "SL_orbits" = 20. The difference of 6 is multi-edge orbits.

At n=3,4: multi = 0, so the values coincide — masking the error (same pattern as MISTAKE-029).

### The correct framing
- **gap_orbits** = edge_orbits - E = self_loop_orbits + multi_orbits (RENAME from "SL_orbits")
- **self_loop_orbits** = #{S_n-orbits on {T, T^e} with T ≅ T^e} = 2, 5, 14, ... (computed via Burnside)
- **multi_orbits** = #{edge orbits connecting already-counted class pairs} = 0, 0, 6, ...

### Impact
- The recurrence search for "SL_orbits" in S242/S244 was wasted effort on a DERIVED quantity (= T/2 + (n-2)! - E). Any pattern found is just a pattern in E in disguise.
- The formula edge_orbits = T/2 + (n-2)! is CORRECT and independently valuable.
- Future work should target E(G_n) directly, not the gap.

### Lesson
**Name quantities precisely.** "SL_orbits" was never defined as "self-loop edge orbits" — it was defined as "edge_orbits - E" and then ASSUMED to count self-loops. The assumption failed at n=5. Always verify definitions against direct computation before building analysis on them.

---

## MISTAKE-031: Tiling complement ≠ tournament complement

**Date discovered:** 2026-03-24
**Found by:** Devil's advocate audit (kind-pasteur-S20ex)
**Affects:** wiggly_metagraph_deep_s20ev.py, aw_precision_s20ew.py, wiggly_patterns_s20eq.py, unified_weights_s20et.py

### What was assumed
Flipping all TILE bits (`mask ^ ((1<<m)-1)`) gives the complement tournament.

### Why it was wrong
The tiling model has m = C(n-1,2) tiles (non-base-path arcs). Flipping tile bits only reverses these arcs, leaving the n-1 base path arcs unchanged. The true tournament complement reverses ALL C(n,2) arcs. These produce different tournaments.

### Impact
V_merged was wrong at n>=5: got 9 (should be 10) at n=5, 41 (should be 34) at n=6. All spectral analysis, W/A comparisons, and eigenvector correlations in affected scripts were computed on the WRONG quotient graph. Corrected in wiggly_corrected_s20ex.py.

### Lesson
When working in the tiling model (fixed base path), always compute the complement via the ADJACENCY MATRIX (reverse all arcs), not via bit flipping on tiles.

---

## MISTAKE-032: Grid-symmetric fraction formula was wrong

**Date discovered:** 2026-03-24
**Found by:** Devil's advocate audit (kind-pasteur-S20ex)
**Affects:** CLAUDE.md (line about "Grid-sym fraction = exactly 2^{-(n-2)}")

### What was assumed
The fraction of grid-symmetric tilings is 2^{-(n-2)} for all n.

### Why it was wrong
The correct formula is 2^{(floor((n-1)/2) - C(n-1,2))/2}, giving exponents 0, -1, -2, -4, -6, -9 for n=3..8. The formula 2^{-(n-2)} gives -1, -2, -3, -4, -5, -6 which only matches at n=5,6.

### Impact
Claims about blue fraction being exactly 2^{-(n-2)} per class are wrong. The correct formula accounts for the number of fixed tiles on the anti-diagonal of the staircase.

### Lesson
Always verify formulas at multiple n values, not just the first few where coincidences can mask errors.

---

## MISTAKE-033: Confused complement-tiling with complement-tournament in blue/black analysis

**Date discovered:** 2026-03-24
**Found by:** User correction (opus-S295)
**Affects:** three_graphs_s295.py, wiggly_vs_lines_s275.py reasoning

### What was assumed
Blue/black lines were modeled as connecting T to T^op (the tournament complement, flipping ALL C(n,2) arcs including base-path arcs). This led to the claim that blue/black "lives outside Q_m" and has zero cross-class edges in the merged meta-graph.

### Why it was wrong
In the tournament-tiling-explorer, a blue/black line connects a TILING to its COMPLEMENT TILING = flip all m tiles (XOR with 2^m - 1). This stays INSIDE Q_m. The complement tiling gives a tournament where all non-base-path arcs are reversed but base-path arcs are PRESERVED. This is NOT T^op (which reverses ALL arcs).

The complement tiling IS at Hamming distance m in Q_m. It gives a different labeled tournament that may be in a DIFFERENT iso class. Blue/black lines DO create cross-class edges in both the unmerged and merged meta-graphs.

### The correct framing
- **Complement TILING** = flip all m tiles = bits XOR (2^m - 1). Stays in Q_m. THIS is what blue/black lines are.
- **Complement TOURNAMENT** (T^op) = flip all C(n,2) arcs. Leaves Q_m (changes base-path arcs). NOT the same as complement tiling.
- Blue/black lines ARE in Q_m, they ARE waggly lines (at distance m), and they DO connect different iso classes.
- The S295 analysis incorrectly modeled blue/black by computing T^op instead of the complement tiling.

### Impact
The three_graphs_s295.py blue/black weight matrix is WRONG. The claim "blue/black is purely diagonal" is FALSE. Must recompute using the correct definition: complement tiling = XOR all tile bits.

### Lesson
ALWAYS check definitions against the actual explorer behavior. The tiling complement (flip tiles) and tournament complement (flip all arcs) are different operations. In the tiling model with fixed base path, flipping all tiles does NOT give T^op.

---

## MISTAKE-034: Band-limitedness at m/2 does NOT hold at n=5

**Date discovered:** 2026-03-25
**Found by:** kind-pasteur-2026-03-25-S1
**Affects:** h-is-band-limited.md (opus-S306), OPEN-Q-040 item 2

### What was assumed
"H is EXACTLY zero in upper Walsh spectrum (k >= m/2). PROVED at n=5,6." (From OPEN-Q-040 and h-is-band-limited.md reflection.)

### Why it was wrong
At n=5 (m=6): the Walsh degree of H is 4 = 2*floor((5-1)/2). Since m/2 = 3, the Walsh degree EXCEEDS m/2. There are 7 nonzero Walsh coefficients at weight 4, and alpha_4 = sum of Walsh coefficients at weight 4 = 0.375 != 0.

Additionally, complement symmetry H(t) = H(~t) FAILS in the tiling model because flipping all tile bits is NOT T^op (base-path arcs are fixed). This means odd-weight Walsh coefficients are nonzero (17 at n=5, 907 at n=7).

### The correct framing
- Walsh degree = 2*floor((n-1)/2) for ALL n >= 4 (THM-260, proved via THM-076)
- Band-limitedness at m/2 holds for **n >= 6** (since 2*floor((n-1)/2) < C(n-1,2)/2 iff n >= 6)
- At n=4,5: Walsh degree exceeds m/2 — NOT band-limited at midpoint
- In the tiling model, both odd and even Walsh weights can be nonzero

### Impact
The "upper half vanishes" claim in h-is-band-limited.md needs correction at n=5. The main qualitative finding (H is low-frequency, concentrated in lower Walsh spectrum) is correct and gets STRONGER as n grows. The asymptotic ratio degree/m -> 0 still holds.

### Lesson
When making claims about "all n," verify at the boundary cases (smallest n). The n=5 case is special because m = C(4,2) = 6 is comparable in size to n-1 = 4. For n >= 6, the quadratic growth of m dominates the linear growth of the Walsh degree.

---

## MISTAKE-035: "G_n is a DAG under H-gradient" — False Claim Propagated Across Repo

**Date discovered:** 2026-04-01
**Found by:** opus-2026-04-01-S1 (systematic audit)
**Affects:** CLAUDE.md, OPEN-QUESTIONS.md, 4 reflection files, paper draft, ~20 agent messages, gn_merged_cascade_s221.py (hardcoded output), local_gradient_s186.py (hardcoded CONFIRMED)

### What was claimed
"The meta-graph G_n is a DAG under H-gradient (0 downhill edges, verified n=3..7)" — CLAUDE.md line 326 (pre-fix). OPEN-QUESTIONS.md claimed "HOLDS at n=3..8."

### Why it was wrong
THREE distinct errors compounded:

1. **Trivially true claim conflated with nontrivial property.** For ANY undirected graph with a real-valued function on vertices, orient edges by function value → the result is always a DAG (modulo level edges). This was explicitly noted in `meta_graph_deep_s181.py` lines 366-368 but the insight was never propagated. The REAL nontrivial question is about **level edges** (same H, different class).

2. **Level edges exist from n=5 onward.** G_n level edges: 0, 0, 1, 15, 136 for n=3..7. G_n/Z_2 level edges: 0, 0, 1, 5, 71 for n=3..7. The graph is NOT a strict DAG from n=5 onward.

3. **H-decreasing edges exist at n=7.** `merged_n7_deep_s20co.out` shows: G_7 has uphill=2988, downhill=962, level=136. G_7/Z_2 has uphill=1633, downhill=419, level=71. The "downhill" count here reflects edges where the class with more neighbors (higher index) has LOWER H — these are real H-reversals in the metagraph. `gap_inventory_s196.py` correctly listed this as REFUTED.

4. **Hardcoded output bugs.** `gn_merged_cascade_s221.py` line 487 prints "DAG: 0 H-decreasing edges (all n)" unconditionally, even though its own data (line 68 of output) shows "DAG: Y, Y, N, N, N, N" for n=3..8. `local_gradient_s186.py` prints "CONFIRMED: all negative-DeltaH flips stay in-class" unconditionally even when the script found counterexamples.

### The correct framing
- G_n has a **strong H-gradient**: most edges increase H. The ratio uphill/(uphill+downhill) is 100% at n≤6 (for the nontrivial edges), and ~76% at n=7.
- G_n is NOT a strict DAG from n≥5 (level edges) and has H-decreasing edges from n≥7.
- The level edge fraction stays small (~3-5%) and may decrease asymptotically.
- The H-gradient is a useful organizing principle but not an absolute law.

### Impact
- CLAUDE.md, OPEN-QUESTIONS.md, 4 reflection files, paper draft all corrected (opus-2026-04-01-S1).
- Every new agent session was reading this false claim and propagating it.
- The `unlocking-gn-at-all-n.md` file listed H-DAG as a "Proved Structural Law" — it was not proved and is not true.

### Lesson
**Three compounding failures:** (1) A trivially-true observation was mistaken for a nontrivial theorem. (2) The discoverer of the triviality (meta_graph_deep_s181.py) did not propagate the correction. (3) Later scripts hardcoded "CONFIRMED" messages that print regardless of results. When a claimed property is trivially true, that's a red flag that you're measuring the wrong thing.

---

## MISTAKE-036: Diameter conjecture diam(G_n) = n-2 is WRONG

**Date discovered:** 2026-03-23
**Found by:** kind-pasteur (gap_inventory_s196)
**Affects:** the-isomorphism-class-graph.md, merged-metagraph-invariants.md, multiple broadcast messages

### What was claimed
"Diameter of G_n is n-2" — conjectured based on n=3 (diam=1), n=4 (diam=2), n=5 (diam=3).

### Why it was wrong
At n=6: diam=4 = n-2 (still holds). At n=7: diam=**7** ≠ 5 = n-2. At n=8: diam=**8** ≠ 6. The actual growth is closer to quadratic (~n²/4), not linear. The diameter-is-feedback-arc-set.md reflection explains: diam ≈ max FAS count difference, which grows quadratically.

### The correct values
diam(G_n) = 1, 2, 3, 4, 7, 8 for n=3..8.

### Impact
- `merged-metagraph-invariants.md` self-contradicts: says "CONFIRMED" at line 84 and "REFUTED" at line 172.
- `the-isomorphism-class-graph.md` still lists "Prove diameter = n-2" as an open problem.
- Multiple broadcast messages from S170, S177, S305 assert or propose proving diam=n-2.

### Lesson
Patterns that hold for 4 consecutive values (n=3..6) can still fail at n=7. Always test at the next case before conjecturing.

---

## MISTAKE-037: H-convexity conjecture is FALSE

**Date discovered:** 2026-03-23
**Found by:** kind-pasteur-S20ch
**Affects:** gap_inventory_s196.py line 176

### What was claimed
That the H-landscape on G_n is "convex" — a specific technical condition about H values along paths in the metagraph.

### Why it was wrong
Refuted at n=6 by kind-pasteur-S20ch. Specific counterexample documented in gap_inventory_s196.py.

### Impact
Low — this was a tentative conjecture, not widely propagated.

### Lesson
Convexity-like properties in combinatorial spaces are fragile and should be tested thoroughly before conjecturing.

---

## MISTAKE-049: SC(n) = A000568(n-1) — Fabricated Identity

**Date discovered:** 2026-05-07 (oracle session)
**Found by:** oracle-2026-05-07
**Affects:** `07-reflections/product-graph-sc-spine-fractal-dimensions.md`

### What was assumed
The reflection claimed SC(n) = A000568(n-1), "verified n=2..10," with a table showing SC(3)=1, SC(5)=4, SC(7)=56, SC(8)=456, SC(9)=6880 — all matching A000568(n-1).

### Why it was wrong
The correct SC values from THM-283's Burnside formula are SC(3)=2, SC(4)=2, SC(5)=8, SC(6)=12, SC(7)=88, SC(8)=176, SC(9)=2752, SC(10)=8784. These do NOT match A000568(n-1) except at n=4,6 (coincidences). The previous session's code had a bug that produced wrong SC values, and two coincidental matches (n=4 and n=6) created a false pattern.

### The correct framing
The true identity is **SC(2m) = A(m, 4)** where A(n,q) = Σ_{odd λ of n} q^{c(λ)}/z(λ) is the q-deformed tournament count. A(n,2) = A000568(n) and A(m,4) = SC(2m). This is proved algebraically via the doubling bijection λ → 2λ, which gives c(2λ)=2c(λ)+K and z(2λ)=2^K·z(λ), so 2^{c(2λ)}/z(2λ) = 4^{c(λ)}/z(λ).

### Impact
- Medium: the false identity was in a reflection file only, not in canon theorems.
- The CORRECT identity (SC(2m)=A(m,4)) is new and provably correct.
- The correct SC values are already in THM-283 and anti_aut_integration_s20ci.out.

### Lesson
Two coincidental matches in a sequence identity are not verification. Always run the sequence through at least n=8 where the values diverge significantly. The Davis/SC partition Burnside formula should be the canonical source for SC values, not ad-hoc code.

---

## MISTAKE-050: H=63 Reintroduced as a Universal Lean Theorem

**Date discovered:** 2026-05-29
**Found by:** opus-2026-05-29-S8
**Affects:** `04-computation/lean/TournamentH7/H63.lean`, `HSpectrum.lean`, `SUBMISSION.md`, `OPEN-Q-055`, HYP-1754

### What was assumed

The Lean formalisation introduced a theorem/axiom `H_ne_sixtythree` claiming
H(T) ≠ 63 for every tournament T, citing exhaustive n≤7 evidence.
`HSpectrum.lean` bundled this into a universal forbidden trio {7,21,63}.

### Why it was wrong

This repeats MISTAKE-024. H=63 is already known to be achievable at n=8.
The S8 audit re-verified a concrete n=8 counterexample from
`h63_verify.out`:
- H(T)=63 by Held-Karp DP
- H(T)=63 by direct enumeration of all 8! vertex permutations
- Ω(T) has 31 directed odd cycles and is the complete graph K31
- Therefore I(Ω(T),2)=1+2·31=63, matching OCF

### The correct framing

H=63 is a temporary n≤7 gap, not a permanent forbidden value.
The Lean theorem is now demoted to:
`H_ne_sixtythree_le_seven (hn : n ≤ 7)`.
The universal forbidden bundle is {7,21}; the finite n≤7 bundle is {7,21,63}.

### Impact

HYP-1754 is REFUTED. OPEN-Q-055 has been corrected. Any document saying
"universally forbidden {7,21,63}" should be treated as stale unless it explicitly
means n≤7.

### Lesson

Finite exhaustive evidence must carry its finite quantifier into Lean. A theorem
with no `n≤7` hypothesis turns computational evidence into a false universal
axiom. Also: H=63 unlocks in the simplest possible OCF shape, Ω=K31, so the
old disconnected-factor obstruction was measuring the wrong graph shape.

---

## MISTAKE-051: Universal TRRT Revived Despite THM-025 Counterexample

**Date discovered:** 2026-05-29
**Found by:** opus-2026-05-29-S8 during repo scour
**Affects:** OPEN-Q-047, OPEN-Q-051/052/053 priority labels, INV-189/INV-186, HYP-1729

### What was assumed

Newer notes revived the Tournament Real-Rootedness Theorem (TRRT): for every
tournament T, I(Ω(T),x) has all real negative roots. The revived entries cited
small samples at n=9,10 with zero failures and treated the Hermite-Biehler
program as a route to a universal theorem.

### Why it was wrong

Canon THM-025 already disproves universal real-rootedness at n=9. The explicit
counterexample has score sequence [1,1,3,4,4,4,6,6,7] and
I(Ω,x)=1+94x+10x²+x³. Newton's inequality fails at k=2:
10² < (3/2)·94·1, so the polynomial has non-real roots.

### The correct framing

Real-rootedness is proved for n≤8 via claw-freeness and is common in samples,
but it is not universal. The right open problem is to characterize the
real-rooted subclass and locate the THM-025 failure inside any
Hermite-Biehler/interlacing framework.

### Impact

OPEN-Q-047 is retitled as a characterization problem. The HB lemmas are
downgraded from "critical to prove universal TRRT" to "important for the
real-rooted subclass program." HYP-1729 is marked REFUTED as a universal
theorem.

### Lesson

Sampling cannot override a canon counterexample. Before reviving a conjecture,
search `01-canon/theorems/` and `MISTAKES.md` for explicit disproofs.

---

## MISTAKE-052: THM-390 claimed twice in one day (codex-S547 vs codex-S548)

**Date discovered:** 2026-06-01
**Found by:** monad-reviewer-2026-06-01 (QC startup audit)
**Affects:** `01-canon/theorems/THM-390-*`, HYP-2036, HYP-2038, definitions.md,
TANGENTS.md, results/INDEX.md, hypotheses/INDEX.md, reflections, SESSION-LOG

### What happened

Two **distinct, both-PROVED** LRC theorems were independently filed under the same
id THM-390 on the same day:
- codex-2026-06-01-**S547** — `lrc-padic-zero-branch-cover-core` (committed fa44a9d):
  the denominator-sieve semantics (`z_q=0 ⇒ t=1/q` witness) and the minimum AP
  open cover `U_n={u: 2u≥n}` of size `floor(n/2)`.
- codex-2026-06-01-**S548** — `lrc-zero-branch-star-core-peeling` (committed 2264cf3):
  a single q-grid zero-branch star has empty strict endpoint-protection core, with
  explicit peel layers `|C|·m_s`.

S548 did not notice S547 had already taken THM-390 (concurrent sessions, both under
the `codex` line). The collision made every `THM-390` reference ambiguous — HYP-2036
in particular cited both theorems under the one number.

### Why it matters

Ambiguous canon ids break `depends_on` graphs and citations: a reader cannot tell
which theorem a reference means. This is the same class of issue as the
THM-013/THM-082 filename collisions (resolved as THM-012b / THM-084) and the
MISTAKE-018/018b renumber.

### Resolution

First claimant keeps the number. **S547 cover-core stays THM-390; S548 star-peeling
renumbered to THM-391.** File renamed, all star-pointing references updated
(definitions.md, TANGENTS, results/INDEX, hypotheses/INDEX, HYP-2036 [now depends on
both], HYP-2038, two reflections, the verifier script, SESSION-LOG entry). Both
proofs were independently re-derived and are correct (see verification notes in each
theorem file). Historical inbox/broadcast messages left as-is.

### Lesson

Before filing a new `THM-N`, run `ls 01-canon/theorems/ | grep THM-N` to confirm the
id is free — especially in concurrent multi-agent sessions where two agents may pick
the same "next" number on the same day. The repo still carries older unresolved id
collisions (THM-260×3, THM-338×2, THM-336/337 dups); those are latent debt that
should likewise be renumbered when next touched.

---

## MISTAKE-053: Systemic HYP-number collisions — five `HYP-N` reused in one 30-hour LRC burst

**Date discovered:** 2026-06-02
**Found by:** monad-reviewer-2026-06-02 (QC startup audit)
**Affects:** HYP-2050, HYP-2052, HYP-2058, HYP-2061, HYP-2063 (and their INDEX
entries, files, reflections). This is MISTAKE-052 (the THM-390 collision)
repeating at scale for the `HYP-*` namespace.

### What happened

Between 2026-06-01 and 2026-06-02, three concurrent agent lines (`opus`,
`oracle`, `codex`) ran the LRC@14/n=17 frontier in parallel and each picked the
same "next" HYP number within **3–12 minutes** of one another. Five collisions:

| HYP | First claimant (UTC) | Second claimant (UTC) | Both have a file? |
|-----|----------------------|------------------------|-------------------|
| 2050 | codex-S551 tetration 20:53 | oracle-S549o Lean 20:56 | only codex |
| 2052 | opus-S551 sieve-no-completeness 21:11 | oracle-S552 loneliness-spectral-gap 21:21 | **BOTH** |
| 2058 | oracle-S553o almost-lonely 15:03 | opus-S556 proof-lite-and-tension 15:21 | only opus |
| 2061 | oracle-S555o pinch-time-pigeonhole 17:41 | codex-S558 small-pinch-shield 17:54 | only codex |
| 2063 | opus-S559 2q-tight-tuple-apex 18:03 | codex-S559 n17-prime-gate 18:15 | **BOTH** |

### Why it matters

Same as MISTAKE-052: an ambiguous id breaks `depends_on`/citation graphs — a
reader cannot tell which hypothesis "HYP-2061" means. THM-396 already
`depends_on: HYP-2059, HYP-2060`, and HYP-2059's INDEX entry chains into HYP-2061,
so the ambiguity reaches a canon theorem's dependency closure.

### Resolution (this session)

- **HYP-2063 (both-file collision, newest):** fully renumbered. First claimant
  opus keeps `HYP-2063` (2q-apex); codex's n17-prime-gate → **HYP-2069**. File
  renamed, INDEX/SESSION-LOG/TANGENTS updated, 0 stray refs remain.
  **Caution — the frontier is a live race:** my first reassignment to `HYP-2064`
  immediately collided *four ways* — a rebase mid-session pulled in oracle-S557o
  (gap-bound), codex-S560 (gate-skip-transfer, has file), and monad-researcher-S560
  (A000568-Burnside), all independently filed under `HYP-2064` within hours. I moved
  codex-n17 clear of the contested 2050–2068 band to **HYP-2069**. The residual
  three-way `HYP-2064` (oracle-S557o / codex-S560 / monad-researcher-S560 — not my
  artifacts) is left to its owning sessions + the cleanup session, banner-flagged in
  the INDEX (suggest HYP-2070/2071 by first-commit timestamp). monad-researcher-S560
  already self-flagged it as a known 3-way collision in its SESSION-LOG entry.
- **HYP-2052 (both-file collision, older, 16 refs):** documented but **NOT yet
  renumbered** — the reference web is dense and a botched mass-rename would create
  more inconsistency than it removes. Canonical assignment: opus-S551
  `lrc-sieve-no-finite-completeness` is first claimant and keeps `HYP-2052`;
  oracle-S552 `lrc-loneliness-spectral-gap` is the duplicate and should be
  renumbered (suggested **HYP-2065**) in a focused future cleanup. Until then,
  always disambiguate by the file slug, not the bare number.
- **HYP-2050 / 2058 / 2061 (single-file collisions):** the idea that owns the file
  keeps the number (minimizes churn); the file-less duplicate (always an `oracle`
  index/reflection entry) is latent debt — suggested reassignments HYP-2066
  (oracle almost-lonely, ex-2058), HYP-2067 (oracle pinch-pigeonhole, ex-2061),
  HYP-2068 (oracle Lean-formalization, ex-2050). Disambiguate by slug meanwhile.

### Lesson

The MISTAKE-052 lesson ("`ls | grep` before filing") was logged for `THM-*` but
not adopted for `HYP-*`, and the failure rate is far higher because HYP numbers
advance many times per day across ≥3 concurrent lines. **Reserve the id first
(Step 5c checkpoint) before doing the work**, and `grep "HYP-N" 05-knowledge/hypotheses/INDEX.md`
+ `ls 05-knowledge/hypotheses/ | grep HYP-N` immediately before `finish_session`.
A sub-300-second reservation push at session start would have prevented all five.
Latent renumber debt remaining: HYP-2052 (both-file), and the three single-file
oracle duplicates above.

**Additional pre-existing two-file HYP collisions found in the same audit** (older
than this 24h window — full latent debt list for the future cleanup session):
- HYP-1969: `lrc-h-phase-plateau` vs `lrc-proof-route-currencies`
- HYP-1992: `lrc-n18-observer-source-gate-battlefield` vs `lrc-rapidity-formal-group-bridge`
- HYP-1995: `lrc-exact-gap-race-wall-ledger` vs `lrc-twin-roots-of-unity-bridge`
- HYP-2009: `lrc-polygon-outside-inside-arcs` vs `resonance-debt-conjecture`
- HYP-2040: `lrc-conditional-clearance-wedge-transitivity` vs `lrc-n4-measure-gap-unique-tight`

These confirm the collision rate has been chronic across the whole LRC era, not a
one-off. The cleanup session should resolve all of them by first-commit-timestamp
and rebuild a contiguous HYP index.

## MISTAKE (oracle-2026-06-03-S576o): pinch-M with a gcd(m,C)=1 filter gives SPURIOUS LRC counterexamples
When computing the loneliness radius M(S)=max_t min_i ||v_i t|| as a max over PINCH times
t=m/(v_a+v_b) (HYP-2075: the optimum is a pair-sum pinch), you MUST range over ALL
m=1..C-1, NOT only the coprime m (gcd(m,C)=1). The optimal pinch need not be in lowest
terms: e.g. S=(1,4,5) has M=1/3 attained at t=2/6 (pair-sum C=1+5=6, m=2, gcd=2), which a
coprimality filter drops, yielding a false M=2/9 < 1/4 -- a spurious "counterexample" to the
PROVEN LRC(4). Symptom: bounded-speed censuses report min M < 1/n at small n where LRC is a
theorem. Fix: drop the gcd filter (evaluate every m). Caught in lrc_even_ladder_selfconverse_proof_s576.py.

## MISTAKE (monad-compute-2026-06-03-S4): minH_strong(m)=m²−5m+9 is a 4-point coincidence; true value at m=7 is 25 not 23
HYP-2180 (opus-S599s) fit the strong-tournament Hamiltonian-path minimum minH_strong(m)=3,5,9,15 (m=3..6, exhaustive) to the quadratic m²−5m+9 and used a *near-transitive scan* to assert minH_strong(7)=23. EXHAUSTIVE enumeration of all 2^21 tournaments on 7 vertices (reversal-halving, `strong_H_spectrum_m7_exhaustive_monad_s4.py`) gives **minH_strong(7)=25, not 23** — and 23 is not a strong-tournament H-value at m=7 at all. The quadratic matched m=3..6 only by coincidence (same trap as MISTAKE-028/036: a pattern holding for 4 values then failing). ~~The CORRECT law is the known **Busch (2006) recurrence p(n)=p(n−1)+p(n−2)+1** for the minimum number of Hamiltonian paths in a strong tournament, giving 3,5,9,15,25,41,67,….~~ **[SUPERSEDED by MISTAKE-055, monad-compute-2026-06-06-S5/S6:** this recurrence is a MIS-CITATION of Busch. Exhaustive iso-class enumeration gives minH_strong(8)=**45** (not 41) and minH_strong(9)=**75** (not 67). The recurrence p(n)=p(n−1)+p(n−2)+1 itself fits only m≤7 then breaks at m=8 — the very same trap it was logged to correct. Busch's TRUE published sequence is **3,5,9,15,25,45,75,125,225,…**, which the exhaustive computation reproduces exactly.]** Everything else in HYP-2180 survived the exhaustive check: 7,21,63 are NOT strong values at m=7; 35=7·5 and 49=7² do fill in; only {7,21} are permanent H-gaps (63 achievable at n=8). Lesson (again): fit a candidate closed form only after it is verified at the FIRST genuinely new case, and trust a near-transitive scan for nothing more than a lower bound.

## MISTAKE-054: Incremental 3-cycle counter swapped i-beats-j / j-beats-i (under-pruning)

**Date discovered:** 2026-06-04 (monad-compute-2026-06-04-S4)
**Found by:** monad-compute, via ground-truth disagreement with the direct-count engine
**Affects:** the FIRST version of `h21_finite_check_v2_monad_s4.py` (the DFS-pruned
extension `extend()`); FIXED before any result was reported.

### What happened
The fast engine v2 builds each new vertex's orientation by DFS, accumulating the
new 3-cycles `{j, i, new}` incrementally and pruning when partial `c_3 > CAP`.
The triple's out-degrees were coded as
  `dj = ij + (1-nj)`,  `di = ji + (1-ni)`
i.e. vertex `j`'s out-degree used `ij` (i beats j) instead of `ji` (j beats i),
and symmetrically for `i`. Because the cycle test requires BOTH `di==1` and
`dj==1`, this is **not** a harmless relabel — with `nj`/`ni` attached to the
wrong term it tests a different condition, so some true 3-cycles were not counted.

### Symptom
v2 reported MORE iso-classes with `c_3<=10` than the direct-count engine v1
(m=7: 453 vs 339; m=9: 17,667 vs 2,575). Both engines still reproduced A000568
with the cap removed (the bug only affects the *capped* count), which is why the
no-cap self-validation did not catch it.

### The fix / correct framing
For triple `{j, i, new}` with `j<i`:
  `dj = ji + (1-nj)`  (j beats i? + j beats new?),
  `di = ij + (1-ni)`  (i beats j? + i beats new?),
  `dn = ni + nj`.
3-cycle iff `dj==di==dn==1`. After the fix, v2 matches v1 EXACTLY for m<=10 and
runs ~10x faster.

### Lesson
A no-cap / total-count self-check does NOT validate threshold/pruning logic.
Always cross-check a fast pruned engine against a slow direct-count engine on the
ACTUAL filtered quantity (here `#{iso classes with c_3<=10}`), not just the total.

## MISTAKE-055: Busch (2006) strong-min recurrence mis-cited as p(n)=p(n−1)+p(n−2)+1 (gives 41,67); true minH_strong is 3,5,9,15,25,45,75

**Date discovered:** 2026-06-06 (monad-compute-2026-06-06-S5/S6)
**Found by:** monad-compute, via exhaustive iso-class enumeration at m=8 and m=9
**Affects:** the MISTAKE-(2026-06-03-S4) entry above; HYP-2180; HYP-2271's "Busch-type" reduction; opus-S699j/k's strong-min(8)≤45 search bound; any downstream use of the 41/67 values.

### What was assumed
The prior monad-compute session (2026-06-03-S4), while correcting an *earlier* bad fit (the quadratic m²−5m+9), asserted that the minimum number of Hamiltonian paths in a strong tournament obeys the recurrence p(n)=p(n−1)+p(n−2)+1, giving 3,5,9,15,25,**41**,**67**,… and attributed this to Busch (2006).

### Why it was wrong
That recurrence matches the EXHAUSTIVE values 3,5,9,15,25 (m=3..7) but BREAKS at m=8 — the identical "holds for several values then fails" trap (cf. MISTAKE-028/036/054) the entry was written to warn against. EXHAUSTIVE enumeration of all non-isomorphic strong tournaments (generated by canonical augmentation, validated against A000568 = …,456,6880,191536 for n=7,8,9) gives:

  minH_strong(m) = 3, 5, 9, 15, 25, **45**, **75**   for m = 3..9   (NOT …25,41,67)

opus-S699j/k's non-exhaustive reversal-search bound strong-min(8) ≤ 45 was therefore TIGHT (=45), not loose; and strong-min(9)=75.

### The correct framing
Busch, "A Note on the Number of Hamiltonian Paths in Strong Tournaments" (Electron. J. Combin. 13 (2006), #N3) proves the minimum equals Moon's (1972) upper bound, with sequence **3, 5, 9, 15, 25, 45, 75, 125, 225, 375, 625, …** (n≥3). The exhaustive computation reproduces this EXACTLY through m=9. Empirically the data satisfies p(n)=3·p(n−2) for every step except n=7 (25 vs 27); the asymptotic growth is ~(√3)^n. (Do NOT re-fit a clean recurrence here without checking against Busch's closed form.)

### Impact — POSITIVE for the program
- HYP-2271 (opus-S699j/k) reduced the delta-field polarization / "7,21 never H" theorem to the lower bound **strong-min(m) ≥ 22 for all m≥7**. Busch (2006) proves the minimum is 25,45,75,… (strictly increasing, ≥25 for m≥7) FOR ALL n ⟹ the reduction is CLOSED by a published theorem, not just "Busch-type, to be proven".
- {7,21} are confirmed absent from the strong H-spectrum exhaustively for m≤9 (7,21,35 below strong-min; 49,63 ARE strong values at m=8). Combined with strong-component multiplicativity H=∏H(Cᵢ), this verifies the phantom-volume theorem (only {7,21} are durable forbidden H, genus-2 multiplicative semigroup) for all tournaments whose strong components have ≤9 vertices, and reduces the all-n statement to the published Busch bound.

### Lesson
When citing a literature recurrence, verify its VALUES against the first genuinely new exhaustive case before propagating it. The "41/67" recurrence was adopted as the fix for a bad fit and itself silently inherited the same coincidence failure mode. Exhaustive iso-class enumeration (via gentourng/nauty-style canonical augmentation — here a pure-Python canonical-augmentation generator validated by A000568) makes m=8,9 cheap (6880 / 191536 classes) where labeled enumeration (2^28 / 2^36) is not.

---

## MISTAKE-056: Signed-LRC worry-set "split" claimed first at n=14 — it is first at n=8

**Date discovered:** 2026-06-06
**Found by:** monad-explorer-2026-06-06-S708b
**Affects:** opus-S699 reflection `signed-lrc-theory-sign-is-a-cut-and-the-worryset-splits-s699.md`, HYP-2262 (the "MAIN RESULT" narrative), and the broadcast MSG-001 ("n=14 is the FIRST n whose C=2n−1 admits a doubled-speed shell-partner"). Does NOT affect the theorems T1–T4 (all correct).

### What was claimed
"Through n=7 every tight (M=1/n) config is shell-partner-free; it FAILS at n=14 (V*=(1..11,13,24), shell-partner 3+24=27). n=14 is the FIRST n whose C=2n−1 admits a doubled-speed shell-partner (24=2·12)."

### Why it is wrong
S699 verified n=4,5,6,7 (shell-partner-free) and then jumped straight to the *known* n=14 frontier, never checking n=8,10,12. But **n=8 already carries shell-partner tight configs.** Exhaustively (exact M, and independently the S592 floor test), the n=8 worry-set has 3 floor-tight primitive configs and **two carry a shell-partner**:
- `(1,2,3,4,5,7,12)` = AP_8 with 6→12, where 12=2·6≡−3 (mod 15), shell-partner (3,12), 3+12=15=2·8−1. M=1/8.
- `(1,4,5,6,7,11,13)`, shell-partner (4,11). M=1/8.

The first is the SAME "double the (n−2) speed" mechanism as n=14's V* (double 12→24). So n=8 (C=15=3·5) is the first n whose C admits a doubled-speed shell-partner tight config — not n=14.

### The correct framing
"tight ⟹ no shell-partner" holds for n≤7 and FAILS first at **n=8**. The V*-type (shell-partner-carrying tight) configs form the infinite **Family II** = AP_n with (n−2)↦2(n−2), floor-tight ⟺ **n≡2 (mod 6)** = every even n with 3∣(2n−1) = {8,14,20,26,…} (verified exact n≤29). n=14 is special only as the first such n whose C is a pure prime power (3³). The shell-partner is always (3, 2(n−2)). See HYP-2281 / reflection `the-worryset-split-is-at-n8-shell-transversality-as-the-gauge-invariant-s708.md`.

### Impact
- The "split exists / is finer than M" conclusion STANDS — only the "first n" is corrected (8, not 14).
- POSITIVE: gives a minimal, SOLVED (LRC(8) is true) laboratory for the prime-2×prime-3 doubling mechanism that recurs unsolved at n=14; the (3,24) carry attack should be prototyped on n=8's (3,12).
- Also reframes the carrier as a purely UNSIGNED, gauge-invariant property: "carries a shell-partner" ⟺ "S mod 2n−1 is not a shell-transversal" (HYP-2281 L1–L2).

### Lesson
When a property is verified up to n=k and then claimed to "first fail" at some larger known-frontier n=N, CHECK every n in (k,N). The interesting frontier (here n=14, C=3³) is rarely the *first* instance of a phenomenon; the first instance (n=8, C=3·5) is usually smaller, more tractable, and already solved.

## MISTAKE-057: THM-427 + HYP-2294 + T765 triple-claimed by two concurrent monad-explorer-S3 instances

**Date discovered:** 2026-06-06 (monad-explorer-2026-06-06-S3, the gcd-torsion lane, at close-out)
**Found by:** monad-explorer (self, on post-push `ls` of theorem dir)
**Affects:** `01-canon/theorems/THM-427-*`, HYP-2294 (INDEX), T765 (TANGENTS), and the two-tower reflection/script. Same class as MISTAKE-052 (THM-390 dup) / MISTAKE-053 (HYP-* dups).

### What happened
Two DISTINCT, both-good LRC results — both responding to the same dispatched seed's "find the unifying statement" — were filed by two concurrent `monad-explorer-2026-06-06-S3` instances under the SAME ids THM-427 / HYP-2294 / T765:
- **gcd-torsion lane** (commit 63ed166, 2026-06-07 01:38:09 UTC): composite-LRC cell-leak `= N_i·n − g·W_i(g)`, a function of `gcd(r,n)=n/ord(r)`.
- **two-tower lane** (commit dba3832, 2026-06-07 01:46:44 UTC): the clock ℤ/n × shell ℤ/(2n−1) coprime-CRT witness group.

The two-tower commit landed ~8.5 min later, when the gcd-torsion THM-427 was already on origin — it did not rebase-detect the taken id (the live-race failure mode of MISTAKE-053).

### Resolution
First claimant keeps the number (gcd-torsion, earlier commit + already on origin). The two-tower lane renumbered: **THM-427→THM-428, HYP-2294→HYP-2295, T765→T766**. Theorem file `git mv`'d; self-references flipped in the two-tower theorem file, its reflection, its script+`.out`, and the shared INDEX table-row / TANGENTS entry. 0 stray refs remain (the two results are complementary: gcd-torsion = mod-n leak face, two-tower = the mod-n ⟂ mod-2n−1 CRT product — they reinforce, not conflict).

### Lesson
The MISTAKE-053 fix ("reserve the id at Step 5c BEFORE the work; `ls 01-canon/theorems | grep THM-N` immediately before finish") still was not adopted. Sub-300s reservation pushes at session start would have prevented this. When two agents share a machine-name line (`monad-explorer`) and a date, the `[machine]-[date]-S[N]` id does NOT disambiguate concurrent instances — both became "S3". Consider a per-instance random suffix when a line is run in parallel.
---

## MISTAKE-058: a THIRD concurrent monad-explorer-S3 lane (signed-pairwise) also hit THM-427/HYP — the collision was 3-way, not 2-way

**Date discovered:** 2026-06-06 (monad-explorer-2026-06-06-S3, the **signed-pairwise** lane)
**Found by:** monad-explorer (this session), at session-end rebase — after MISTAKE-057 (the two-tower
lane) had already documented the *gcd-torsion vs two-tower* pair.
**Affects:** completes MISTAKE-057. The same window saw a THIRD distinct LRC result claim `THM-427`.

### What happened
MISTAKE-057 recorded TWO concurrent `monad-explorer-2026-06-06-S3` lanes colliding on `THM-427`. There
was a **THIRD**: this signed-pairwise lane (`THM-427-signed-pairwise-floor-is-a-maxcut-LRC`,
`Gstar ≥ 1/(2 r_min)`), committed 20:51:32 -0500 — after gcd-torsion (20:37) and two-tower (20:46).
Three distinct, all-good LRC theorems under one id `THM-427`, all from the same instance name.

### Resolution (first-come keeps the number; consistent with MISTAKE-057)
- `THM-427` → **gcd-torsion** (first). `HYP-2294`, `T765` → gcd-torsion.
- `THM-428`, `HYP-2295`, `T766` → **two-tower** (second; self-renumbered, MISTAKE-057).
- `THM-429`, **`HYP-2296`**, T764-update → **signed-pairwise** (third, this lane): file `git mv`'d,
  id + `signed_lrc_rmin_bound_monad_s3.py` docstring updated; HYP renumbered 2295→**2296** (2295 is
  the two-tower's), references flipped in THM-429, the reflection, INDEX, TANGENTS, SESSION-LOG. My
  already-pushed *commit messages* still say THM-427/HYP-2295 (immutable history); the canon files are
  **THM-429 / HYP-2296**.

### Lesson
Even after a collision is "resolved," re-check before finishing: a 2-way resolution can be incomplete
if a third concurrent instance is in flight. And renumber by first-commit author-date end-to-end
(here gcd-torsion < two-tower < signed-pairwise ⟹ 427/428/429, 2294/2295/2296). The deeper fix
remains MISTAKE-053's: reserve ids at Step 5c before doing the work; when a `[machine]-[date]` line is
run ≥3-way in parallel, the `S[N]` suffix does not disambiguate — use a per-instance random tag.
The Step-5c "reserve the id first, `ls | grep` before filing" rule (MISTAKE-053) must run **even
against your own instance id** — concurrency can duplicate the *session name*, not just the number.
A one-line reservation push at session start (claiming THM-N/HYP-N as honest stubs) would have
prevented all three. When three files share `THM-N`, resolve by first-commit author-date, not by
who notices last.

**ADDENDUM (same session, on rebase):** it was a THREE-way race, not two. A *third* concurrent
`monad-explorer-2026-06-06-S3` filed THM-427 = "signed pairwise floor is a max-cut LRC"
(commit 20:58 -0500, 20 min after the gcd-torsion claim; it also forward-referenced HYP-2294 for
its asymptotic question). First-claimant rule again: gcd-torsion keeps THM-427; the max-cut lane →
**THM-429**, its HYP-2294 forward-ref → **HYP-2296** (free). Three independent S3 instances,
three THM-427 claims, all within 20 minutes — the strongest evidence yet for per-instance id
suffixes when one agent line is fanned out in parallel.

---

## MISTAKE-059: "Exactly 3-to-1" inferred from a count ratio without checking the map (caught + corrected same session)

**Date discovered:** 2026-06-07 (monad-explorer-S6, self-caught)
**Affects:** THM-436 ADDENDUM (2″) as first checkpointed; HYP-2305; reflection `the-icosahedral-fifteen-s6.md` (all corrected before any agent built on them)

### What was assumed
The commutator map {60 oriented overlapping cyclic-triangle pairs on a 5-set} → {20 three-cycles of A_5}
was stated as "**onto and exactly 3-to-1**" — inferred purely from `60 / 20 = 3`, and dressed up as the
icosahedral **face-vertex flag** incidence (`20` faces × `3` vertices = `60` flags, flag→face uniformly
3-to-1).

### Why it was wrong
The fibers are **not uniform**. Direct enumeration (`04-computation/icosahedral_flag_fibers_monad_s6.py`)
gives fiber sizes `{2 (×3), 3 (×14), 4 (×3)}` (sum 60 over 20 three-cycles). The `3`-to-`1` holds only
on average. So `60 = |A₅|` is the group order, NOT a flag count, and the commutator covering is NOT the
icosahedral flag map.

### The correct framing
What is actually true and robust: **every one of the 60 oriented overlapping pairs has a commutator of
cycle-type a 3-cycle** (conjugation/inversion-invariant ⇒ order-convention-independent), and the 60
commutators are **onto all 20** three-cycles. That type-uniformity — not any multiplicity-uniformity — is
the real content of "A_n perfect realized by overlapping triangle pairs."

### Lesson
A matching TOTAL (`60 = 20·3`) does not certify a uniform MAP. When a count coincides "too cleanly" with
a known structure (here, icosahedral flags), verify the **fibers / the map**, not just the cardinality.
This is the project's own "too clean ⇒ test it" rule applied to itself; the test refuted the clean story
and left the honest one (type-uniformity + a 15-fold canonical bijection) standing.

---

## MISTAKE-060: THM-438 "bigon-trees ARE the Catalan count" — top order is a SIGNED cactus cancellation, not a +1-per-tree count

**Date discovered:** 2026-06-07 (monad-explorer-2026-06-07, deep-research / analytic lane, 3rd session)
**Found by:** monad-explorer, while attempting the "small remaining write-up" (the +C_k sign) flagged OPEN in THM-438's Honest-status section
**Affects:** THM-438 Part B proof MECHANISM + error term; the reflection `the-paley-cluster-integrals-are-catalan-numbers-tree-walks-and-the-moment-method.md` ("Patterns with any non-bigon cycle ... are O(p^{k+1/2})"; "the top order is an all-bigon graph ... a tree of bigons ... counted by C_k"); the reflection's stated O(1/sqrt p) convergence. **Does NOT affect** the STATEMENTS A_{2k}=C_k p^{k+1} or R(p)->e — both stand (verified).

### What was assumed
The leading order p^{k+1} of the cluster integral `A_{2k} = sum_{distinct x_0..x_{2k}} prod chi(x_{i+1}-x_i)`
is reached ONLY by all-bigon coincidence patterns; a **tree** of k bigons maximizes V=k+1; each such bigon-tree
(= Euler tour of a plane tree) contributes **+1**, so the leading coefficient is literally the Catalan count
C_k. "Patterns with any non-bigon cycle ... are O(p^{k+1/2})." Error term O(p^{k+1/2}).

### Why it was wrong (verified exactly, `04-computation/paley_cluster_cactus_census_monad.py`)
Three things are false:
1. **Bigon-trees do NOT each contribute +1.** Via the partition-lattice Moebius inversion
   `A_{2k} = sum_sigma mu(0,sigma) M_sigma`, a bigon-tree pattern `sigma` carries Moebius weight
   `mu(0,sigma) = prod_blocks (-1)^{|B|-1}(|B|-1)!`, which is NOT 1 when a vertex is visited >=3 times.
   The bigon-tree leading coefficient (sum over non-crossing edge-pairings of `prod_v (b_v-1)!`) is
   **1, 3, 13, 69, 421, 2867** (k=1..6) = **OEIS A088368** (g.f. `A=sum n! x^n A^n`, `a(n)~e*n!`) —
   NEITHER C_k NOR (2k-1)!!. At k=2 bigon-trees give **3**, at k=3 they give **13** (census confirms).
2. **Even cycles DO reach the top order p^{k+1}.** The single 2k-cycle pattern (`x_0=x_{2k}`) equals
   `tr(M^{2k}) = (-p)^k(p-1) ~ (-1)^k p^{k+1}` — the SAME order as bigon-trees, not O(p^{k+1/2}).
   It enters with `mu=-1` and SUBTRACTS. More generally every **even cactus** (connected graph whose
   biconnected blocks are all even cycles, total half-edges k) contributes at p^{k+1}.
3. **The Catalan number is a signed cancellation.** Census:
   `k=2: bigons(+3) + 4-cycle(-1) = 2 = C_2`;  `k=3: bigons(+13) + {bigon+4cyc} + {6cyc} = 5 = C_3`.

### The correct framing
**Closed form (PROVED via Gauss-sum inversion `chi(w)=g^{-1} sum_t chi(t) omega^{tw}`, verified exactly):**
```
M_sigma = (-1)^k * p^{V-k} * F(sigma),    F(sigma) = sum over F_p-flows t on G_sigma of prod_e chi(t_e),
```
V = #blocks, flow space = cycle space (dim m = 2k-V+1). A pattern reaches p^{k+1} iff F reaches full order
p^m; those are exactly the **even cacti**. The leading coefficient of A_{2k} is the **signed sum over even
cacti** `sum mu(0,sigma) * lead(M_sigma) = C_k` — an inclusion-exclusion that converts the all-pairings
overcount (A088368, ~e*n!) into the **non-crossing** count C_k. This is the genuine free-probability /
moment-method content: the two-point Gauss spectrum's even-cycle terms are PRECISELY the corrections that
take Gaussian-style pairings to the semicircle's non-crossing pairings.

**Error term:** `A_{2k} = C_k p^{k+1} + O(p^k)` (NOT O(p^{k+1/2})). Verified: `(A_4-2p^3)/p^2` is STABLE
(~ -7.1..-7.8 -> ~-8), while `/p^{2.5}` drifts to 0. Hence R(p)-e has relative correction **O(1/p)**,
resolving the reflection's stated O(1/sqrt p) vs the close-out's "favors 1/p" IN FAVOR OF 1/p.

**Part C simplifies:** R(p)->e needs **NO Weil bound**. The only V=2k no-leaf pattern is the single
2k-cycle = `tr(M^{2k}) = (-p)^k(p-1)` (elementary); V<2k is trivially `O(p^{2k-1})=o(p^{2k})`.

### Impact
- THM-438 Part B mechanism CORRECTED (addendum added). Statements A_{2k}=C_k p^{k+1}, R(p)->e UNCHANGED.
- Part C upgraded: fully elementary (no Weil).
- Error term corrected p^{k+1/2} -> p^k; convergence rate of R->e pinned to 1/p (feeds HYP-2307 #2).

### Lesson
The project's own rule again: a clean final count (C_k) reached by a clean-sounding mechanism
("bigon-trees = Catalan") does not certify the mechanism. The Moebius weights and the
equal-order even-cycle patterns were invisible at the level of "count the bigon-trees." Always
decompose the inclusion-exclusion and check which patterns share the leading order — here the
cancellation `A088368 -> C_k` is the actual phenomenon, and it is the free-probability fingerprint
the moment-method slogan was pointing at.

---

## MISTAKE-061: THM-438 — the top-order patterns are NOT "even cacti"; they are the larger class of EVEN-SERIES patterns (even theta graphs included)

**Date discovered:** 2026-06-07
**Found by:** monad-explorer-2026-06-07 (deep-research / analytic lane, 4th session)
**Affects:** THM-438 ADDENDUM and MISTAKE-060 (the *characterization* of which coincidence
patterns reach the leading order `p^{k+1}`). Does NOT affect the Catalan law `A_{2k}=C_k p^{k+1}`
itself (re-confirmed here, rigorously) nor `R(p)->e`.

### What was assumed (MISTAKE-060 / THM-438 ADDENDUM)
"`M_sigma` reaches the top order `p^{k+1}` **iff** `F(sigma)` reaches full order `p^m` — exactly
the **even cacti** (connected, all biconnected blocks even cycles)." The census then grouped the
leading coefficient as bigon-trees (+A088368) corrected by even-cycle **cacti** down to `C_k`.

### Why it was wrong
`F(sigma) = sum_{flows} prod_e chi(t_e)` reaches full order `p^m` iff the flow-form product
`P(s) = prod_e ell_e(s)` is a **perfect square** (then `chi(P)=chi(Q^2)=+1` off the zero locus,
so `F ~ +p^m`). `P` is a perfect square iff **every series-class of edges has even size** (each
distinct flow-line occurs an even number of times). The even cacti satisfy this — but so do
**even theta graphs** (two vertices joined by three even paths; biconnected block is NOT a single
cycle) and, generally, all "even series-parallel" 2-connected patterns. These are NON-cacti yet
reach `p^{k+1}` and MUST be counted. Verified (`04-computation/paley_cluster_theta_check_monad.py`):
at `k=3` the `V=5, m=2` top-order patterns are **6 even-cacti{2,4} + 1 even-theta(2,2,2)** — the
even theta (mu=+1) was invisible to the "even cacti" census (it sat in the `(6,)` biconnected
bucket, silently cancelling the single 6-cycle, so the *total* still came out right).

### The correct framing (VERIFIED k<=4; the `g` step PROVED)
Let `c0 = lim A_{2k}/p^{k+1}`. Then
```
c0 = (-1)^k * sum_{rho : connected, EVERY series-class even}  mu(0,rho) * g(rho),
```
and `g(rho) := lim F(rho)/p^m = +1` for EVERY such pattern. **`g==+1` is PROVABLE:** within each
series-class the closed Euler walk passes straight through the degree-2 internal vertices, so all
edges of the class get the SAME orientation sign `s in {+1,-1}`; the class is even, so
`prod_{e in class} sign_e = s^{even} = +1`; hence `P = (prod sign_e) Q^2 = +Q^2` and
`g=chi(P)=+1`. Therefore the entire character/Gauss-sum content collapses and
```
($$)   sum_{rho : even-series pattern}  mu(0,rho)  =  (-1)^k C_k        (number-theory-FREE).
```
RIGOROUSLY CONFIRMED `c0 = 2, 5, 14 = C_2, C_3, C_4` by clean Richardson (`1/p`) extrapolation of
the exact flow-Moebius value (`04-computation/paley_cluster_topterm_monad.py`) — this also REPLACES
the prior slowly-converging census (which read `1.56, 2.77, 3.11` at `p<=23` and only *looked* like
it might reach `5`). The breakdown:
```
k=3:  bigon-trees(m=3) +13,  (m=2: cacti+theta) -9,  (m=1: 6-cycle) +1   = 5 = C_3
k=4:  bigon-trees(m=4) +69,  (m=3) -72,  (m=2) +18,  (m=1: 8-cycle) -1   = 14 = C_4
```
(bigon-tree sub-sums `+13, +69` = OEIS A088368, the all-pairings overcount, as before).

### Impact
- THM-438 ADDENDUM-2 added: Catalan law `A_{2k}=C_k p^{k+1}` RE-CONFIRMED (rigorous, k<=4), error
  `O(p^k)` unchanged, `R(p)->e` unchanged.
- The MECHANISM is corrected a SECOND time: top-order class = **even-series patterns** (perfect-square
  flow product), strictly larger than even cacti. `g==+1` is proved, reducing handoff #1 to the
  clean number-theory-free Moebius identity `($$)`.
- Free-probability reading SHARPENED: the random skew-Rademacher matrix gives `C_k` *directly* from
  non-crossing pairings (each `+1`, no factorials); the deterministic Paley Moebius expansion
  over-counts to A088368 in the bigon sector and the even cacti + even thetas + ... cancel it back to
  `C_k`. The equality is Wigner quasirandomness; `($$)` is its exact combinatorial fingerprint.

### Lesson
MISTAKE-060 corrected the *value* mechanism (bigon-trees -> A088368 -> C_k) but inherited a wrong
*support*: "even cacti." A pattern can saturate the flow character-sum without being a cactus — any
even series-parallel skeleton does. When a leading-order census gives the right TOTAL, that does not
certify the per-class STRUCTURE: a missing pattern (the even theta) can hide inside a coarse bucket,
cancelling against another, leaving the total correct and the story wrong. Characterize the support
by the actual saturation condition (perfect-square flow product = even series-classes), not by the
most familiar sub-family.

---

## MISTAKE-062: even-series pattern count is NOT OEIS A215257 — a 5-term coincidence that breaks at k=6

**Date discovered:** 2026-06-07 (monad-explorer, 8th session)
**Found by:** monad-explorer-2026-06-07 (deep-research, 8th session)
**Affects:** THM-438 ADDENDUM-3 point (2); HYP-2308; the reflection `the-drt-engine-is-S-squared-equals-J-minus-nI-the-catalan-is-genus-zero.md`; INDEX/SESSION-LOG entries asserting "even-series count = A215257"

### What was assumed
THM-438 ADDENDUM-3 (5th session) identified the number of EVEN-SERIES patterns of the
path `[0..2k]` (= the unsigned support of `(**)`) as **OEIS A215257**: the values for
`k=1..5` are `1, 3, 13, 67, 383 = A215257(k+1)` (indecomposable deque-sortable
permutations). The recursion script hardcoded the *predicted* next value `A215257(7)=2345`
for `k=6` but NEVER actually computed `k=6` (its `KMAX` default was 5).

### Why it was wrong
A direct exhaustive count at `k=6` (fast integer enumerator
`04-computation/paley_starstar_triangle_fast_monad.py`, cross-validated against the original
SVD test `04-computation/paley_starstar_crosscheck_monad.py` with **0 disagreements** over
all `Bell(13)=27.6M` partitions exhaustively at `k<=5` and a 300k sample at `k=6`) gives
```
   even-series count, k=1..6  =  1, 3, 13, 67, 383, 2351.
```
The OEIS b-file gives `A215257(7) = 2345 != 2351`. An OEIS search for
`1,3,13,67,383,2351` returns **no results** — the unsigned even-series count is not (yet)
a catalogued sequence. The A215257 match was a **5-term small-number coincidence**.

### The correct framing
- The unsigned even-series pattern count is `1, 3, 13, 67, 383, 2351, ...` (computed,
  rigorous through k=6), NOT A215257, and presently matches no OEIS sequence.
- This does NOT touch any headline result. The Moebius-SIGNED sum
  `(**) S_k = sum_{even-series} mu(0,sigma) = (-1)^k C_k` is independently re-verified
  exhaustively at `k=6` (`S_6 = 132 = C_6`), as is the cycle-rank triangle row
  `t(6,m) = 1, 45, 560, 2626, 4845, 2867` and the loop equation
  `S_k = -sum_{i+j=k-1} S_i S_j`.
- If anything the refutation SHARPENS the thread's thesis: the *unsigned* count is so
  unstructured it is not even a known sequence, while the *signed* sum is the cleanest
  possible (Catalan). "The Catalan is a cancellation, not a count" is now literal.

### Impact
- THM-438 ADDENDUM-3 (2) corrected (see ADDENDUM-6). HYP-2308 / INDEX A215257 cells updated.
- The "indecomposable deque-sortable permutations" bijection program (ADDENDUM-3/4 handoff)
  is moot — there is no A215257 bijection to find because the counts differ.

### Lesson
A 5-term OEIS hit is weak evidence — A215257 and the even-series count share five terms by
chance. NEVER hardcode an "expected next" OEIS value as if computed; compute it. Generic
divergence of two integer sequences after a short common prefix is the default, not the
exception (cf. MISTAKE-006 ratio coincidence, MISTAKE-010 small-n pattern break).

---

## MISTAKE-063: THM-438 ADD-9 wrongly "refuted" `A088368(m) ~ e*m!` — the original asymptotic is CORRECT (Kotesovec/OEIS); ADD-9 sampled only the pre-peak rising side

**Date discovered:** 2026-06-07 (monad-explorer, 12th session)
**Found by:** monad-explorer-2026-06-07 (deep-research, 12th session), via direct OEIS lookup of A088368
**Affects:** THM-438 ADDENDUM-9 point (6) ("CORRECTION FLAG"); any reflection/INDEX/SESSION-LOG
line asserting "A088368(m) ~ m!(m+2)/2, NOT e*m!"

### What was assumed (the erroneous "correction")
ADD-9 point (6) flagged the long-standing claim "`A088368(m) ~ e*m!`" (the diagonal
`t(k,k)` of the cycle-rank triangle = the all-pairings overcount) as **NOT supported by the
data**, observing `A088368(m)/m! = 1, 1.5, 2.17, 2.875, 3.51, 3.98, 4.45` (m=1..7) is
"monotonically increasing PAST e ≈ 2.718", and proposed instead the empirical
`A088368(m) ≈ m!(m+2)/2` (ratio `(m+2)/2`).

### Why it was wrong
The asymptotic `a(n) ~ e * n!` is an **established OEIS result** for A088368 (Vaclav
Kotesovec, Apr 10 2019; verbatim on the A088368 page). Computing the ratio with the
**true** b-file values (ADD-9 also had a transcription slip: `A088368(7) = 21477`, not
"22417") shows `a(n)/n!` does NOT increase monotonically — it OVERSHOOTS e, **peaks at n=8
(≈ 4.359)**, then **strictly DECREASES** back toward e:
```
   n:        2     3     4     5     6     7     8(peak) 9    10    12    16    20
   a(n)/n!: 1.50  2.17  2.88  3.51  3.98  4.26  4.36   4.32  4.19  3.85  3.36  3.14  ... -> e
```
ADD-9 sampled only `m<=7` — entirely on the **rising side, before the peak** — and mistook
the slow, overshooting approach for divergence. ADD-9's `(m+2)/2` fits the rising side only
and diverges (it predicts 11.0 at n=20, where the true ratio is 3.14 and falling).
Verified: `04-computation/paley_starstar_diagonal_noncrossing_monad.py`.

### The correct framing
- `A088368(m) ~ e * m!` is CORRECT. The diagonal of the cycle-rank triangle grows like
  `e * m!` (this is the "wild end" of the bridge polynomial; equivalently `h_m(m) =
  A088368(m)/m! -> e`).
- A088368 = **"number of partitions of [n] into sets of NONCROSSING LISTS"** (Callan,
  arXiv:0711.4841), G.f. `A(x/F(x)) = F(x)` with `F(x) = sum n! x^n`. It is a named,
  closed-form object — the diagonal is NOT "uncatalogued" (only the OFF-diagonal columns are).
- The original `~e*m!` slogan (ADD, ADD-8) should be RESTORED; ADD-9 point (6) is retracted.

### Impact
- THM-438 ADD-9 point (6) retracted; ADDENDUM-10 records the restoration + the noncrossing-
  lists identity. No headline result was ever affected (`A_{2k}=C_k p^{k+1}`, `R(p)->e`,
  `(**)`, column rationality all stand).
- The "tame<->wild bridge" of ADD-8/9 now has explicit asymptotic endpoints: `h_m(m) -> e`
  (wild/A088368) and `h_m(-1) -> 0` (tame/Mersenne, super-exponentially).

### Lesson
The MIRROR of MISTAKE-062. There, a sequence MATCH was over-trusted from 5 small terms;
here, an asymptotic was REFUTED from 6 small terms. **A factorial-scale ratio that is still
changing at n<=8 has converged to nothing** — slow asymptotics routinely overshoot and turn
around (here the turn is at n=8). Before declaring an `~ c*n!` claim false from a finite
ratio table, (i) check OEIS for an established asymptotic, and (ii) extend the ratio far
enough to see whether it is still rising — a monotone prefix is not a limit.

---

## MISTAKE-064: misread Erdős Problem 64 as "even cycle" (2k) when it is "power-of-2 cycle" (2^k); proved the trivial settled statement and wrongly framed it as the open problem

**Author:** opus-2026-06-08-S708 (caught by the user same day)
**Severity:** framing/attribution error (the math proved is correct but answers the WRONG question)

### What happened
The user asked to "Work On: Does every finite graph with minimum degree at least 3 contain a
cycle of length 2𝑘 for some 𝑘 ≥2?" — **Erdős Problem 64 = the Erdős–Gyárfás conjecture (1995):
every graph with min degree ≥3 contains a cycle whose length is a POWER OF TWO** (`2^k`: 4, 8,
16, …). This is **OPEN and falsifiable.** I read "2𝑘" as "2·k" (an even number ≥4), proved the
classical longest-path result *"min degree ≥3 ⟹ an even cycle of length ≥4"* (TRUE and settled),
and **wrongly presented it as the user's problem** in THM-443/HYP-2313/S708.

### Why it is wrong
"Min degree ≥3 ⟹ even cycle" is settled-true, so it CANNOT be an open falsifiable problem — that
alone should have flagged the misreading. The real condition is **multiplicative/2-adic** (length
exactly a power of 2), enormously stronger than "even." Petersen (girth 5) has no 4-cycle but an
8-cycle; a counterexample would be a min-degree-3 graph avoiding 4, 8, 16, 32, … *all at once*.

### Fix
- THM-443 **rescoped**: it correctly proves the EVEN-cycle lemma (classical) — relabeled to state
  explicitly that this is NOT Erdős 64. The power-of-2 problem is treated honestly as OPEN in
  THM-444/S709 (verified computationally on small/structured graphs; no counterexample; proven for
  cubic planar by Heckman–Krakovski; Markström computational searches).
- HYP-2313's parity-covering lens stands as a lens but its even-cycle leg is downgraded to the
  weaker statement; the power-of-2 ("dyadic") version is the open one.
- The S707 Pfaffian/even-dicycle bridge applies to EVEN cycles (Pólya), NOT to power-of-2 cycles.

### Lesson
**When a problem is stated to be OPEN, first check that your proof does not settle it — if it
does, you have misread the problem.** A one-character misread (`2^k` vs `2k`) flips a trivial
exercise into a famous open conjecture. Parse the *difficulty claim* as a constraint on the
interpretation.

---

## MISTAKE-065: Tile-Bit Negation Under Path Reversal — T^op Along a Reversed Path is a Grid Transpose WITHOUT Negation

**Date discovered:** 2026-06-09 (caught in-session by branch-C verification before propagating)
**Found by:** kind-pasteur-2026-06-09-S1 (hand derivation wrong; computational branch corrected it)
**Affects:** THM-447(5) original wording, HYP-2335 original wording, T767 original wording (all corrected in place)

### What was assumed
In the canonical frame of the skew-Sylvester double D(T) (path, twin, reversed path), the copy-2
tile block was claimed to be the "grid-transposed NEGATED" copy of T's tiling x — reasoning:
copy 2 is T^op (all arcs reversed), so its tile bits must be complemented.

### Why it was wrong
Copy 2 is traversed along the REVERSED base path. Reversing the path also reverses the
upper/lower convention of every tile, which complements each bit a second time. The two
negations cancel exactly: t(X,Y) = x(σ_n(X,Y)) — grid transpose with NO negation. Verified in
100% of 1098 framed + 1096 labeled cases (n=3..6).

### The correct framing
The single negated Sylvester copy in the tile schema [[H,H],[H,−H]] lives in the CROSS block:
σ-partner cross tiles (i,j) ↔ (j,i) carry complementary bits A[i][j] vs A[j][i]. See THM-452.

### Lesson (same family as MISTAKE-033)
Arc reversal (T^op), path reversal, and tile-bit complementation are THREE involutions that
compose in non-obvious ways. T^op + reversed path = grid transpose, NOT complement-transpose.
Whenever a claim involves "negated copy" at the tiling level, track ALL involutions explicitly
— two of them silently composing to the identity is the recurring trap. (See also the broader
NUMBER-conflation hygiene reference [[polysemous-constants-bridges-traps-and-homonyms]] (klein-S7):
the "2" homonym above is one of several constants — 2, n, 7, 14, 28, 6, 3 — that wear arithmetic vs
dimensional hats; run the PERSISTENCE TEST before treating any numeric coincidence as structure.)

### Impact
None propagated: corrected in THM-447(5-CORRECTED), THM-452(1), HYP-2335 status, T767 note,
all within the same session.
## MISTAKE-066: Erdős 592 finite-bridge direction stated BACKWARDS in first tree-grid script

**Date:** 2026-06-09
**Found by:** mac-mini-2026-06-09-S1 (same session, caught while writing the pysat version)
**Affects:** `04-computation/erdos592_treegrid_dichotomy_macmini_s1.py` original docstring
(corrected in place); draft doc §3.2 v1 (corrected)

### What was assumed
"An infinite witness for ω^n ↛ (ω^n,3) restricts to finite witnesses on every finite
grid, so the negative relation implies Q(n,t) SAT for all t; UNSAT at some t is evidence
for the positive direction."

### Why it was wrong
An infinite witness only guarantees that no FULL-TYPE (infinite) independent set exists.
A finite binary subgrid is not of full type, so nothing forces any finite subgrid to
contain a blue edge: the restriction of an infinite witness to a finite grid can be
empty. The implication as stated is unsupported in BOTH directions at the finite level.

### The correct framing (THM-453 part D)
The true bridge runs the other way and is a compactness statement:
- Q(n,t) SAT for ALL t ⟹ (König) a triangle-free graph on the full grid with no
  independent binary subgrid ⟹ ω^n ↛ (ω^n,3) with a STRONG witness.
- Hence positive relations FORCE finite cutoffs: ω^n → (ω^n,3) ⟹ R(n,2) < ∞.
  Specker at n=2 forces R(2,2) < ∞ even though SAT witnesses persist through t=10:
  the cutoff is real but large.
- A negative ordinal relation does NOT formally imply Q(n,t) SAT for all t
  ("strong witness" is a priori stronger than "witness").

### Lesson
For infinite-to-finite shadows, check WHICH quantifier compactness actually transports.
"Kill all full-type sets" has no finite trace on a single finite configuration; only
the universal finite family (all t at once) carries ordinal content, and it does so in
the direction finite-SAT-everywhere ⟹ infinite negative. Cf. MISTAKE-064 (parse the
statement before proving the wrong one).

---

## MISTAKE-067: incomplete subgrid-verifier falsely certified Q(n,t) SAT — R(2,2) is actually 5, not >14

**Date:** 2026-06-09 (same session, caught by structure-reading the "witnesses")
**Found by:** mac-mini-2026-06-09-S1
**Affects:** `erdos592_treegrid_pysat_macmini_s1.py` (find_independent_binary_subgrid),
`erdos592_invariant_quotient_macmini_s1.py` (imports it), results files
`erdos592_treegrid_pysat_*.out`, `erdos592_treegrid_push_*.out`,
`erdos592_invariant_quotient_*.out` (their SAT lines beyond the corrected thresholds);
fixed + fully re-run in `erdos592_verify_fix_macmini_s1.py`

### What happened
The CEGAR loop's final certificate ("no independent binary subgrid exists in the model")
used a backtracking search that committed to the FIRST consistent subtree under each
chosen child and never explored alternative subtrees of the same child. The search was
therefore incomplete: it could return "none found" when an independent subgrid existed,
falsely certifying SAT. Q(2,t) was reported SAT through t=14 ("R(2,2)>14"); with a
complete generator-based search, **Q(2,5) is UNSAT: R(2,2) = 5 exactly**
(Q(2,4) SAT with 35 edges).

### How it was caught
Reading the structure of the t=11 "invariant witness": its printed B_1 visibly missed
rectangles like {2,3}×{4,5}, which a genuine witness must hit — the certificate and the
object contradicted each other on inspection.

### What was NOT affected
- All UNSAT results (solver-side, sound): the n=1 calibration, and the corrected runs.
- Labs 1–2 avoidability conclusions: those used POSITIVE results of exists_S_free_grid,
  which constructs explicit grids whose pairwise pattern-checks happen at append time —
  sound. (exists_S_free_grid has the same incompleteness on its FALSE answers; its
  False-derived "caps" are lower bounds only — none were used for headline claims.)

### Lessons
1. In a CEGAR loop the FINAL "no counterexample" check is a proof obligation — it must
   be a complete decision procedure, not the same heuristic used to generate clauses.
2. Read the witness, not just the verdict: printing the (R, B_g) structure exposed the
   hole immediately (cf. MISTAKE-015's "the error was visible in the output").
3. Greedy-commit backtracking (commit to first consistent subtree, never revisit) is a
   recurring trap in tree searches — both grid searchers this session had it.

---

## MISTAKE-068: Cycle-Anchored Subset DP Reused for Longest-Path Problems

**Date discovered:** 2026-06-09 (caught in-session by branch-II self-check P3; never propagated)
**Found by:** kind-pasteur-2026-06-09-S2 branch II (blowup spectrum lab)
**Affects:** any script computing longest paths with a min-vertex-anchored subset DP

### What was assumed
A subset DP that anchors each cycle at its minimum-label vertex (correct for cycle enumeration:
every cycle can be rooted at its minimum) was reused for LONGEST PATHS by anchoring paths at
their minimum vertex.

### Why it was wrong
A path's minimum-label vertex can be INTERIOR — such paths are never generated when extending
only from the anchor. Symptom: longest path of P3 reported as 2 vertices; downstream, the
circumference law c(G[K2]) = 2p showed 275 fake "beaters" (true count after fix: 56).

### The correct framing
Cycles: min-anchoring is safe (rotation moves the minimum to the start). Paths: must allow
two-sided extension from the anchor, or run DP from every start vertex. Cycle spectra computed
with the anchored DP were NEVER affected.

### Lesson
"Anchor at the minimum" is a CYCLE symmetry (rotation), not a PATH symmetry (paths only have
reversal). Check which group acts before porting a canonical-form trick between objects.

---

## MISTAKE-069: "First Power-of-2 Cycle in Enumeration Order" Reported as "Smallest"

**Date discovered:** 2026-06-09
**Found by:** kind-pasteur-2026-06-09-S2 branch III (double-checked: own DFS checker ≡ networkx)
**Affects:** S710's Erdős-64 verification table (THM-446 context): "McGee → C16" and the cage list

### What was assumed
S710's cage battery reported, per graph, "the" power-of-2 cycle found — implicitly the smallest.

### Why it was wrong
The search reported the FIRST 2^k cycle encountered in enumeration order, not the smallest k.
McGee (girth 7) was reported as "→ C16" but in fact contains 34 EIGHT-cycles. (Petersen has 15
eight-cycles, not the 10 sometimes quoted.) The Erdős–Gyárfás CONCLUSIONS are unaffected (a
2^k cycle exists either way), but any downstream use of "McGee is C8-free" would be wrong —
e.g. it would corrupt the girth ladder of THM-457 (the true min order of a girth-7-or-more
cubic C8-free graph is > 46 by SA floors, not 24).

### The correct framing
When verifying "∃ cycle of length in S", always record the FULL dyadic profile (which members
of S occur), or at minimum the smallest. Exact anchor census now in canon (THM-457(1)).

### Lesson
"Found one" ≠ "found the smallest". Enumeration order silently masquerades as minimality.

## MISTAKE-070: Circular Proof Inventoried as PROVED (|Aut| | H "by tiling-count integrality")

**Date discovered:** 2026-06-10
**Found by:** kind-pasteur-2026-06-10-S1 Thread B (adversarially verified: independent canon sweep + line-by-line re-read of every cited source)
**Affects:** `gap_inventory_s196.out` item 15 ("|Aut| divides H for all tournaments: PROVED (opus S182)"); downstream confidence in THM-048 Step 3 and CLAUDE.md's tiling-fibration line

### What was assumed
That the universal divisibility |Aut(T)| | H(T) had been PROVED (opus-S182), so later artifacts cited it freely: THM-048 Step 3 says "by orbit-counting" with no proof; CLAUDE.md line ~351 states "Tilings * |Aut| = H for every iso class (orbit-stabilizer on tiling fibration)".

### Why it was wrong
S182's argument (`aut_H_deep_s182.out` Part 5) derives the divisibility from the integrality of the tiling count H/|Aut| — which is the SAME statement (it even hedges: "orbit size divides |Aut|, and different orbits may have different sizes"). The S20bt tiling formula behind the CLAUDE.md line was verified only at n=4,5 and implicitly ASSUMED freeness ("|Aut| are related by automorphisms"); `forbidden_tiling_counts.out` line ~241 even flags freeness as an assumption. The fact WAS true — but no proof existed anywhere in canon, while the inventory said PROVED.

### The fix
**LEM-003** (kind-pasteur-2026-06-10-S1): for ANY finite digraph, an automorphism fixing a directed Hamiltonian path's arc set is the identity (a directed Ham path's arc set determines its vertex sequence — the unique in-degree-0 source anchors an induction), so Aut acts freely on directed Ham paths, all orbits have size exactly |Aut|, and |Aut| | H. One paragraph, zero tournament structure. Exhaustively verified n≤6 (all 2^10+2^15 labeled tournaments, explicit orbits for all 3432 masks with |Aut|>1) + independently re-verified with different machinery. Honest boundary: FAILS for Hamiltonian CYCLES (no distinguished start): C3 has 1 Ham cycle and |Aut|=3 ∤ 1; circulant RQ5 has BOTH its Ham cycles rotation-fixed.

### Lesson
"X is an integer because it counts something" requires proving the count is well-defined — here, that all orbit sizes equal |Aut| (freeness), which was the entire content. An asserted proof citation can be as hollow as an extrapolated pattern (cf. MISTAKE-028/036/055). When an inventory says PROVED, spot-check the cited argument before building on it.

## MISTAKE-071: "Verified exhaustive n=4,6" That Checked Only One Maximizer Class (HYP-2312)

**Date discovered:** 2026-06-11
**Found by:** mac-mini-2026-06-10-S2 (det census) + independent adversarial re-verification (fresh exhaustive enumeration at n=6)
**Affects:** HYP-2312 ("the H-maximizing tournament has |Pf| = 1"), THM-442 section (3), and any plan to restrict A038375 extremal searches to Pf = ±1 tournaments

### What was assumed
That at n = 4, 6 the H-maximizing tournament has minimal Pfaffian |Pf| = 1, recorded as "VERIFIED exhaustive n=4,6", and conjectured for all even n (HYP-2312).

### Why it was wrong
At n = 6 the maximum H = 45 is attained by TWO iso classes (240 labeled each, both score (2,2,2,3,3,3)): one has |Pf| = 1 but the other has |Pf| = 7 (det S = 49; same in the det(I+2A) convention, so not a convention artifact). The earlier verification evidently examined only one maximizer class. At n = 8 the six H = 661 classes have |Pf| ∈ {1, 9, 17}.

### The fix
HYP-2312 amended to the EXISTENTIAL form: "at least one H-maximizing tournament has |Pf| = 1" — true at n = 4, 6 (exhaustive) and supported at n = 8. The existential form still justifies searching Pf = ±1 tournaments for the max VALUE of H, but such a search misses maximizer CLASSES. Universal form: REFUTED at n = 6.

### Lesson
"The maximizer" is a set, not an element. Any claim of the form "the extremal object has property P" must be checked on ALL extremal classes (argmax can be plural), and the verification record should state how many maximizers were examined.

---

## MISTAKE-072: "dim_nonspec(H) = n−5 / the overlaps are spectral shadows" extrapolated from n≤8 to all n

**Date discovered:** 2026-06-15
**Found by:** monad-explorer-2026-06-15 (n=9 carrier-dimension chain)
**Affects:** THM-505 (original dimension section), HYP-2513, reflection `the-zeta-function-and-the-ocf-read-complementary-halves` (closing paragraph)

### What was assumed
A carrier-dimension probe at n≤8 found that the non-spectral content of `H` is exactly the simple-cycle vector `(c6,…,c_n)`, of dimension `n−5`, with every overlap defect (`p33, TQ, Q44, TF`) a spectral function of it. This was stated for general `n`: "the simple cycles are the genuine hidden coordinates; the overlaps are their spectral shadows."

### Why it was wrong
At **n=9** the simple-cycle counts do NOT determine `H`. Nested-carrier chain over 130 000 cospectral samples: `sig→+(c6,c7,c8)→+c9→+Q44→+T333` splits `14804→482→24→1→0`. There are explicit cospectral witnesses with identical `(c6,c7,c8,c9)` but different `H` (each satisfying `ΔH = 4ΔQ44 + 8ΔT333`), and even `(c6,c7,c8,c9,Q44)` leaves `H` split. So `dim_nonspec(H) = 6 > n−5 = 4` (carriers `{c6,c7,c8,c9,Q44,T333}`), and BOTH `Q44` and `T333` are INDEPENDENT non-spectral carriers — not spectral shadows. The break occurs exactly when the triple level `α_3` (and the `(3,5)`-pair structure) gains room; n=8 was the last size where the higher-correlation configs were pinned by the simple counts for lack of room.

### The fix
THM-505 / HYP-2513 / the reflection all amended: dim = `n−5` ONLY for `n ≤ 8`; at `n ≥ 9` the non-spectral content is a TOWER of cycle-correlation rungs (simple counts; pair-overlaps `Q44`; triple-packings `T333`; …), each rung independent of those below once `n` gives it room. `H` is universal-linear in the full carrier set but not a bounded-degree polynomial in the simple cycles alone past n=7.

### Lesson
The repo's own refrain — "patterns that hold at n=3,4,5 often break at n=6 or n=7" — applies one notch higher here: a dimension/structure pattern verified exhaustively/by-heavy-sampling at `n≤8` STILL broke at `n=9`, precisely at the `n` where a NEW combinatorial level (here the triple disjoint-cycle packing `α_3`, needing `3·3=9` vertices) first has room. When a carrier/dimension count is conjectured "for all n," always test it at the first `n` where the next independent-set / correlation level can appear, not just one size up.

### Addendum 2 (monad-explorer-2026-06-15-S4): even the PACKING basis over-counts *for H* — H's own non-spectral dim is `⌊n/3⌋`, not `A000009(n)−3`
The S3 partition law `dim = #{λ:odd≥3,Σλ≤n}−3 = A000009(n)−3` (now identified as the named
sequence: the cumulative restricted-partition count equals partitions of n into distinct/odd
parts, by the one-line GF identity `Σ_{s≤n}[x^s]Π_{odd≥3}1/(1−x^k)=[x^n]Π_{odd≥1}=q(n)`)
correctly counts the rank of the **individual packing carriers** `(N_λ)`. But it does **not**
measure the non-spectral dimension of `H`. `H = I(Ω,2) = 1 + Σ_{j≤⌊n/3⌋} 2^j α_j` depends only
on the **level-sums** `α_j = Σ_{|λ|=j}N_λ` — it never sees the split of a level into its
length-types. Since `α_j=0` for `j>⌊n/3⌋`, **`dim_func(H)(n) ≤ ⌊n/3⌋` (LINEAR), PROVED**, and
`< A000009(n)−3` for n≥8. VERIFIED n=8: carrier basis `{c7,D33,D35}` rank 3 (the S3 number) but
level-sum basis `{c7, D33+D35}` rank 2 with `H` in span — `D33,D35` are independent functions
yet `H` reads only their sum. **Third lesson:** the chain of corrections (trace basis 6 →
packing basis 5 → level grading ⌊n/3⌋) shows that "the non-spectral dimension of `H`" must be
measured against *what `H` is a function of*, not against any convenient carrier list — and `H`,
being the scalar `I(Ω,2)`, factors through the coarse level-marginals `α_j`. The
`A000009(n)−3` law is true and beautiful, but it is the dimension of the finer **packing
vector**, not of `H`. See THM-505 TWO-DIMENSIONS section + reflection
`H-reads-only-the-level-grading`. (`04-computation/ocf_two_dimensions_monad.py`)

### Addendum (monad-explorer-2026-06-15-S3): the "6" was itself a basis over-count → the real law is a partition function
The corrected `dim_nonspec(H)=6` at n=9 above is itself **basis-dependent and over-counts by one**. It used the TRACE basis `{c6,c7,c8,c9,Q44,T333}`, but `c8` and `Q44` enter `H` **only through their sum** `c8+Q44 = D35` (the disjoint-(3,5)-pair count): the closed form's `+4c8+4Q44` is `4·D35` after `D35 = c3c5 − W8 + c8 + Q44`. In the basis-independent OCF *packing* basis `N_λ` (rank of the within-class carrier-delta matrix), the intrinsic dimension is **5, not 6** — carriers `{c7, c9, D33, D35, T333}`. The general law is `dim_nonspec(H)(n) = #{partitions of s≤n into odd parts ≥3} − 3` = 1,2,3,5,7,9,… (verified rank 3,5,7 at n=8,9,10). **Second lesson:** a "dimension" measured by a *nested chain over a fixed (trace) carrier list* counts basis vectors, not degrees of freedom — over-complete bases inflate it. Measure intrinsic dimension by the RANK of the carrier-delta space in the natural (here OCF packing) basis. See THM-505 growth-law section + reflection `the-non-spectral-dimension-of-H-is-a-partition-function`.

## MISTAKE-073 (2026-06-15, THM-503 infimum) — searched a SLICE, reported it as the infimum
THM-503 reported `inf L ≈ 0.0237` over primitive multiple-of-14 LRC configs, "attained at the near-tight cores {1,…,12}∪{14m}". This searched only the **end-drop** family (drop the LAST runner of the tight AP {1..13}). A broad search over **interior-drop** configs ({1..13}\{j}∪{14m}, j interior) goes ~4× lower: `{1..11,13,84}` has L≈0.00535, `{1..13}\{6}∪56` has L≈0.00561 (both verified loose, M=7/89 & 2/23, L stable to Q₀=24000). The true inf L ≈ 0.0053, not 0.0237; the minimizing drop is the MIDDLE (j=6), not the end. The qualitative inf>0 survived, but the value/extremizer/margin were wrong. Lesson: when reporting an extremum over a symmetric-ish search space, **sweep the full orbit of perturbations, not the first slice** — here "drop a runner of the tight AP" has 13 inequivalent positions and the end (most natural to try) is the WORST extremizer, the middle the best. Echoes the census-horizon lesson (MISTAKE-018) and the parameterize-the-exception lesson: the extremal object is rarely the most obvious member of the family. Fix: THM-503 CORRECTION addendum + HYP-2520.

## MISTAKE-074 (2026-06-15, THM-518 / OPEN-Q-104) — a Fourier-TRUNCATED Riesz ratio reported a false looseness certificate
Probing the Bedert Riesz certificate `∫M·R/∫R < 1` for the LRC(14) extremizer `{1..13}\{6}∪56`, I computed `∫M·R` via the Fourier side `Σ_v Σ_{|k|≤14} s(k) R̂(vk)` (truncating the sinc sum at `Kmax=14`) and got ratio **0.9506 < 1** — flagged as "CERTIFICATE!". The EXACT direct-grid integral `∫M·R = (1/Q)Σ_a M(a/Q)R(a/Q)` (Q=30030, no truncation) gives ratio **1.064 > 1** — NO certificate. The tail `Σ_{|k|>14} s(k)R̂(vk)` is NOT negligible: `s(k)~1/k` decays slowly and `R̂(vk)` is nonzero on the dense relation sum-set, so truncating at `Kmax=14` dropped a positive contribution that flips the verdict. Lesson: for a `∫(periodic)(Riesz)` pairing where the test function `M=Σ1_danger` has slowly-decaying (`1/k`) Fourier coefficients, **compute the integral DIRECTLY on a fine grid — never truncate the Fourier/sinc sum** unless the tail is rigorously bounded. The direct grid is exact and cheap here; the Fourier shortcut silently fabricated a certificate. (Self-caught the same session by running the direct-grid verification BEFORE writing any canon claim — the discipline "verify a too-good certificate by an independent exact method" is what saved it.) Fix: HYP-2547 caveat + THM-518 §C(1); all OPEN-Q-104 ratios now direct-grid.

## MISTAKE-076 (2026-06-17/18, THM-526 / HYP-2580a) — "criterion C verified universal" was a sampling artifact; C is NOT necessary
The LRC(14) proof-program reduced "M(S)≥1/14 for covering 13-sets" to the **criterion** `C(S): ∃v∈S, W(S\{v})>1/(7v)` (arc-width lemma: C ⟹ M≥1/14, genuinely PROVED). HYP-2580a/Angle-F then claimed `C(S)` holds **universally** ("verified on ~12,000 covering sets / 4000/4000, 0 failures"), making "prove C for all covering S" the target. **C is NOT universal.** Exact counterexample (found by two independent verifiers, reconfirmed directly): **S\* = {1,2,3,5,7,8,9,10,11,12,13,38,42}** — primitive, covering, case S3 (k=2, Vmax=42) — has **all 13 ratios W(S\*\{v})·7v < 1** (max 429/532≈0.806), yet **M(S\*)=2/23 ≥ 1/14** (a global witness, not any single-removal arc). So "prove C(S) for every covering 13-set" is a **false target** — it can never close LRC(14). Lesson: a SUFFICIENT criterion validated only by sampling can be silently non-necessary; the C-failure locus here is sparse (≈1 in thousands) and slipped through every Monte-Carlo sweep. Before elevating "criterion holds on all X" to a proof target, **search adversarially for the criterion's failure set** (it is cheaper to refute universality than to prove it), and keep the direct quantity (here M itself, via the global witness / finite check) as the fallback. The PROVED pieces (S1, S2, the k=2 slice, cluster-collapse Lemma A) never used C-universality and stand. Fix: THM-526 correction section + HYP-2580a marked REFUTED + HYP-2581; a correct closing lemma must remove ≥2 runners or use a global witness. Echoes MISTAKE-073 (searched a slice, called it the infimum) and the census-horizon lesson (MISTAKE-018).

## MISTAKE-077 (2026-06-18, μ_θ exact engine) — naive order-change breakpoints UNDER-count the gap measure
Computing μ_θ(E)=meas{x: maxgap{frac(e_i x)}>θ} exactly requires breakpoints where the cyclic ORDER of the orbit points changes (collisions (e_i−e_j)x∈ℤ) AND, WITHIN each order-cell, the analytic crossings where an affine gap equals θ. An engine that used ONLY order-change breakpoints (sampling each cell's midpoint) gave a WRONG, smaller value for μ_{1/7}(consec_9) — 4829/5880 instead of the correct 247/294 — because a gap can cross θ strictly inside an order-cell. The fix (engine A / the workflow's mu_theta): on each order-cell every cyclic gap is affine in x, so {maxgap>θ}=∪_t{affine_gap_t>θ}, each a sub-interval solved exactly (s·x+c>θ); take the union length. Cross-validated against a brute grid with denominator dividing 7·lcm(differences) (which cannot be wrong). Lesson: for a piecewise-linear extremal (max/min of affine pieces) the level-set {f>θ} has breakpoints BOTH at the piece-swaps AND at the θ-crossings inside pieces — sampling cell-midpoints silently drops the latter. Always solve the affine=θ crossings, or cross-check against a provably-fine rational grid. Affects any μ_θ/EWLB computation. → THM-530, HYP-2602.

## MISTAKE-078 (2026-06-19, LRC(14) wide-spread bound) — "verified on canonical families + loose ceiling" mistaken for "uniform theorem"; the envelope tail DIVERGES
Closing the LRC(14) residual to the wide-spread cap bound (meas(S7(E))≤cap_k for span>B), I argued (HYP-2611b, kps-S9) that the route was "structurally complete, remaining is engineering": every tested wide regime sits ≤~0.21 ≪ cap (margin 0.17), the tight case is the exact finite check, so the wide bound is "loose." This conflated **two different things**: (i) the inequality is TRUE and loose on the tested families (~40k exact shapes, 0 violations) — correct; (ii) a UNIFORM-over-all-wide-E PROOF is easy — WRONG. The S10 assembly + 3 verifiers showed the support-6 correction tail eps(B), bounded by the free per-coordinate envelope Σ_{m,7∤m} c1/|m| (c1=0.697303), **DIVERGES HARMONICALLY** (partial sum to 10⁵ = 7.42). So no finite tail bound follows from the envelope; a successive-minima/Minkowski lattice-point count over the support-6 relation lattice (|K(n)|≤c1⁶/(λ₁···λ₆)) is REQUIRED and was never executed. Lesson: a LOOSE margin + exhaustive sampling on structured families is strong EVIDENCE but is NOT a uniform inequality — the uniform statement can still need a genuine (here lattice-geometry) argument that the per-term envelope, summed freely, cannot supply (harmonic divergence). When claiming "the remaining work is engineering," check that the summation/quantifier over the INFINITE part actually converges with an EXECUTED bound, not just that each instance is safe. Echoes MISTAKE-076 (sampling ≠ universality) and MISTAKE-073 (a slice ≠ the infimum). Fix: HYP-2611b corrected; the residual is HYP-2608(a) the wide-spread bound = the Minkowski tail count, OPEN. Everything else (k≤7, glue G1, per-E bound, caps, bounded finite check to span 16) is gap-free.

## MISTAKE-078 AMENDMENT (2026-06-19, CASE-thm538 resolved) — the "support-6 / ≥6-body" mechanism was ALSO wrong, not just the uniformity
MISTAKE-078 flagged that "verified-on-families ≠ uniform theorem" (correct, stands). But its framing also leaned on THM-538's "support-6 floor ⟹ the wide correction is a ≥6-body object" — and that mechanism is itself FALSE (CASE-thm538-support6-floor-zero-padding, conceded): the support-6 floor holds only for the active-coordinate sum Q, NOT the zero-padded kernel K that appears in the measure (short support 2–5 relations DO contribute; the AP's correction is support-3-dominated, not support-6). The CORRECT structure (HYP-2646): K(n)=D7(n mod 7)/∏n_j, correction = Σ_c D7(c)S_c(E) conditionally convergent, ruled by support-6 relation DENSITY R6 (HYP-2645), with the convergent representation being the finite x-cell integral / far-element plateau (HYP-2644/2610), NOT the box-truncated lattice sum. Double lesson: (1) a "verified exhaustively, max 5e-17" check can silently compute a DIFFERENT object (Q, active coords) than the one in the statement (K, zero-padded) — always confirm the verified quantity IS the claimed quantity on a case where they should differ; (2) a clean algebraic vanishing ((1−1)^{6−|U|}=0) can be an artifact of dropping a |T|-dependent factor — re-derive keeping every factor. The wide-spread bound's real content is R6 density + far-element decorrelation, unaffected.


## MISTAKE: L7 'CLOSED' overstated (mac-mini-2026-06-21-S13, caught by the S13 rigor-audit workflow)
**The error:** S13's broadcast / HYP-2738 / reflection claimed the LRC(14)-S3 sector route was "CLOSED" and the audit "ALL PASS". The adversarial audit (Thread A, completed after close-out) found this OVERSTATED. Actual: L7 is REDUCED, not closed.
**Invalid step (L7 Step 4):** the closure (lrc_q108_L7_closure_kps.md) cites the finite-f1 convergence |p0(B u far)-p0_inf| = O(1/f1) as "= THM-546 (PROVED)". FALSE: THM-546 peels ONE far element with the REMAINDER BOUNDED. In the L7 limit BOTH far elements grow (f2 = gamma*f1), so a single THM-546 peel of f2 leaves E'=B u {f1} UNBOUNDED (V = Theta(f1)) => bound (6/7)/gamma ~ 0.43 = O(1), NOT O(1/f1). The proved D_{p,q} <= 14/p bounds the discrepancy of the LIMIT LAW mu_{p,q}, a DIFFERENT object from the finite-f1 convergence RATE. The closure conflated them.
**True status:** the O(1/f1) rate IS true (verified broadly: |err|*f1 bounded ~<= 0.75) but NOT proven; it needs a JOINT 2D Erdos-Turan/Koksma bound peeling the BOUNDED BASE from the FAST FAR-PAIR (not a single peel of f2). Verified-not-proved also: the r>=3 -> pairwise reduction, base-size domination, and "consec maximizes meas(S7)" (HYP-2602, open as a theorem).
**Lean caveat (gap #7):** delsarte_bound_k8/k9/k11 prove the per-shape q0 <= L_y; the content is the readout identity (moment-LP functional = dual covector), and the bound is then trivial (q3,q6>=0). It does NOT formalize L_y <= cap (the extremality / cap content).
**Lesson:** "reduced to a finite atlas + a numerically-verified rate" is NOT "closed". Do not cite THM-546 (single-far, bounded remainder) for the two-far joint limit. The genuine remaining input is the joint 2D ET-Koksma rate. -> HYP-2730, HYP-2738, THM-546, lrc_q108_L7_closure_kps.md.

## MISTAKE-86: computing M(S)=max_t min_s||st|| with an INCOMPLETE breakpoint set underestimates M and reports FALSE tights

**What happened (kind-pasteur-2026-06-28-S256):** A census probe over single-speed replacements of the AP used candidate breakpoints t = k/(2 s_i) and k/(s_i - s_j) only, OMITTING t = k/(s_i + s_j). This made M_of UNDERESTIMATE M(S) (it missed the true maximizer, which often sits at an s_i+s_j crossing), so ~23 single-replacements were reported as "tight" (M=1/14) when they are actually LOOSE (M>1/14). With the complete breakpoint set the census collapsed correctly to the single extra tight set GW (12->24).

**The fix:** f_S(t)=min_i ||s_i t|| is piecewise linear; its local maxima (the candidates for max_t) occur where some ||s_i t|| has a kink (t=k/(2 s_i)) OR where two terms cross. Two terms ||s_i t||=||s_j t|| cross at BOTH t=k/(s_i - s_j) AND t=k/(s_i + s_j). The full candidate set is {k/(2 s_i)} U {k/(s_i - s_j)} U {k/(s_i + s_j)}. Sanity anchors: M(AP)=1/14, M(GW)=1/14, M({1..11,13,26})=1/12. Always include the s_i+s_j family; cross-check with a fine rational grid.

**Lesson:** A "tight" verdict from an exact-rational M computation is only as good as the breakpoint set. Incomplete breakpoints give a one-sided error (underestimate M) that fabricates tights -- the most dangerous direction for a census.

## MISTAKE-087 (2026-06-30, mac-mini-S47) -- the construction n/Phi_6 was assumed to be the COVERING-MIN for n>=7 from a heuristic + a restricted scan; it is NOT (beaten exact at n=7,8,9)

**What happened:** HYP-3701 (my own, S42) inferred a "transition at n=7": drop-2/mediant 2/(2n-1) is the covering-min for n<=6, and the construction {1,..,n-2,n(n-1)} = n/Phi_6(n) takes over as the covering-min for n>=7. The inference rested on (a) the PG(2,6)-failure heuristic (the first projective plane fails at q=6=n-1) and (b) opus's "107-set scan confirmed 14/183" -- a scan of NEAR-CONSTRUCTION variants. The whole subsequent arc (HYP-3703 tiling-optimality, 3704 three-routes, 3717 three-gap, 3722 observer-escape; the Kershner/Eisenstein/A2 framing) built on "convergent = covering-min for n>=7."

**The refutation (exact, dense-grid cross-checked):** smaller-M primitive covering (n-1)-sets exist at every n=7,8,9:
- n=7: {1,2,5,6,7,8}, M=2/13=0.15385 < 7/43=0.16279 (t=4/13)
- n=8: {1,4,5,6,7,11,16}, M=2/15=0.13333 < 8/57=0.14035 (t=8/15)
- n=9: {1,3,4,5,7,11,18,32}, M=4/33=0.12121 < 9/73=0.12329 (t=29/33)
So there is NO transition at n=7; the sub-convergent (mediant at n=7,8) keeps beating the construction. opus's 14/183 is a restricted-family min, NOT the global covering-min.

**Why it slipped through:** the beaters are SPREAD-STRUCTURED (a speed of 32 ~ 3.5n at n=9), so near-construction perturbation scans and low-speed exhaustion both MISS them. The PG-failure heuristic was a post-hoc fit to n<=6, not a mechanism.

**Not fatal to LRC:** all candidates (mediant 2/(2n-1), 4/33, convergent) are > 1/n; the covering-min being smaller just makes the floor margin tighter (~1/(n(2n-1)) vs ~1/n^2). The LRC floor M>=1/n is untouched.

**Lesson:** (1) a "covering-min" claim must be tested against SPREAD-structured covering sets (large speeds), not only near-construction perturbations or low-speed exhaustion. (2) An elegant structure on a particular covering set (the construction's hexagonal AP / three-distance / Eisenstein) does NOT make that set extremal -- do not conflate "a beautiful covering" with "the minimal covering." (3) A transition inferred from a number-theoretic coincidence (PG(2,6)) needs exhaustive confirmation at the first n past the alleged transition. -> HYP-3725, HYP-3701, CASE-convergent-not-covering-min.

## MISTAKE-088 (2026-06-30, mac-mini-S52) -- claimed a clean "n>=12 => covering-min = 1/(n-1)" from an ILP whose speed bound V was BELOW the construction scale n(n-1)

**What happened:** The covering-min ILP (HYP-3731/3732) with speed bound V=72 returned M_prim=1/(n-1) for n=12,13,14,15 (clean!), and I wrote it up as a transition at n=12, the LRC14 hard core = 1/13, and a pinning of HYP-2566's looseness as 1/(n(n-1)). 

**The error:** V=72 is BELOW the construction's largest speed n(n-1) (=156 at n=13, 182 at n=14). The ILP literally could not see the construction {1..n-2, n(n-1)}, whose M = n/Phi_6(n). And n/Phi_6 < 1/(n-1) for ALL n (n^2-n < n^2-n+1), with the construction a valid primitive covering set. So M_prim(n) <= n/Phi_6 < 1/(n-1) ALWAYS -- 1/(n-1) is NEVER the covering-min. The ILP's 1/(n-1) was just the best LOW-SPEED primitive set; the clean "n>=12" regime was a search-bound artifact. (klein-S36 had flagged exactly this under-resourcing in the same-numbered HYP-3731.)

**The fix / what survives:** EXACT (V=4n suffices, klein-confirmed): n=7..11 covering-min 2/13,2/15,4/33,4/37,3/31, depth a(n)=2,2,4,4,3, the Stern-Brocot ray [0;n-1,k] frame (floor k=1, covering-min k=a(n), construction k=n). RETRACTED: the n>=12 clean regime, the n=12 transition, LRC14=1/13, the HYP-2566 pinning. For n>=12 the covering-min is <= n/Phi_6, exact value OPEN (needs V ~ n(n-1)).

**Lesson:** before declaring a computed extremum (min/max), ALWAYS evaluate the KNOWN constructions as upper/lower bounds and check the optimizer's search range covers their scale. Here a single check -- "is the known construction n/Phi_6 below my ILP's answer?" -- would have caught it instantly (it is, by 1/(n^2-n+1)/... a hair). A clean formula emerging right at the edge of the search range is a red flag for under-resourcing, not a discovery. -> HYP-3732, HYP-3731, HYP-2566, THM-523.

## MISTAKE-089 (2026-06-30, klein-S52, caught mid-session before canon) -- searching "covering-min escapes" over ALL sets instead of COVERING sets returns the tight/mediant minimizers as PHANTOM escapes

**What happened:** working the lowness lemma's large-speed residual, I searched sets S = {1..n-2}\{k} U {2 speeds} for the min M, imposing only "kills resonances k, n-1." The search returned M = 1/n (GW) at n=14 and M = 2/(2n-1) at n=10,12,16 -- all BELOW n/Phi_6 -- which looked like the lowness lemma FAILING (a low-M set missing a core speed).

**The resolution:** those minimizers are NOT covering sets. A covering set (THM-523) must contain a multiple of EVERY q in {2,...,n}; GW = {1..11,13,24} and the drop-2 mediant sets miss a multiple of n, so they are THM-523's TRIVIAL class (killed by the q=n witness with M=1/n exactly), NOT things the covering-min ranges over. Re-running with the full covering constraint (mult of every q in {2..n}) gives the genuine single-drop covering-escape minima 5/43, 2/21, 7/89 (n=10,12,14), all > n/Phi_6 as expected (and 4/29 < 8/57 at n=8, the real small-n failure, MISTAKE-087).

**Lesson:** the covering-min and the LRC floor 1/n live on DIFFERENT set-classes. Any covering-min / lowness search MUST impose the covering predicate (mult of every q in {2..n}) up front; otherwise the tight `1/n` and mediant `2/(2n-1)` minimizers (which sit BELOW n/Phi_6) masquerade as escapes and appear to refute the lowness lemma. The large-multiple-forced lemma (HYP-3763) itself is unaffected -- it holds for any low-M set, and GW does carry the forced multiple 24=12*2. -> HYP-3763, THM-523, HYP-3701.

## MISTAKE-091 (2026-07-01, mac-mini-S93, caught in-session before canon; renumbered from 090 -- opus-S32 claimed 090 concurrently) -- applying the tight-set slope formula to NON-PRIMITIVE sets under-counts witnesses by the dilation factor; formula-only censuses without direct cross-verification

**What happened:** the THM-593 tight-set slope formula `c_S = (2/q) sum_{u unit} 1/v_max(u)` silently assumes the argmax of `f_S(t)=min_v||vt||` is exactly the `phi(q)` unit fractions `a/q`. Applied in a cross-modulus census to `{2,4,6,8}` (= 2*AP, non-primitive) at q=5 it returned `5/12` -- but the lonely measure is dilation-invariant (`m_{cS} = m_S` via `u = ct`), so the slope must be `5/6`. A dilated set has `c*phi(q)` witnesses (at `a/(cq)`); the formula misses `(c-1)*phi(q)` of them. The first landscape output also reported this phantom `5/12` as a "beater."

**The resolution:** (1) always normalize to gcd 1 before applying witness-count-sensitive formulas (dilation invariance of measure/slope is the free cross-check: `slope(cS) == slope(S)` catches the bug instantly); (2) the formula is in general a LOWER bound `c_S >= (2/q) sum 1/v_max(u)`, with equality iff the tight set is "clean" (argmax = unit fractions) -- all primitive tight sets found at q=5..16 are clean, verified by DIRECT exact measurement in the last linearity cell; (3) never publish a formula-only census: the corrected landscape recomputes every slope directly and prints `formula=direct` per row.

**Lesson:** a closed form derived under a witness-structure hypothesis must carry that hypothesis explicitly, and any census built on it needs an independent direct computation per entry (here: exact measure in the last Farey cell). Dilation/scale invariances are free consistency tests -- run them on every new invariant. Secondary workflow lesson: `timeout N python | tee out` loses ALL buffered output on timeout with exit 0 (tee masks the kill); use `python -u` + per-section flush for long sweeps. -> THM-593, HYP-3840.
