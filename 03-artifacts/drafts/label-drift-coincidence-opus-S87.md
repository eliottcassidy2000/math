# The label-drift formula and the co-incidence closure (HYP-4156)

**opus-2026-07-05-S87.** The two concrete remaining sliver pieces (HYP-4146) worked to
proofs and exact tables. Numerics: `label_drift_opus_S87.py.out`,
`coincidence_check_opus_S87.out`.

## 1. THE LABEL-DRIFT FORMULA (proved mechanism, exact lengths)

At a resonance t = p/q + ε the cluster offsets sit at (d_j·p/q + d_j·ε): a q-grid
configuration whose points drift at INTEGER rates d_j. Order the offsets circularly;
gap k drifts at rate (d_next − d_k). The union covers (2/13-net persists) exactly while
every gap ≤ 2/13; the desert in ε is the REGIME-CHAINED feasibility interval: piecewise
linear between offset-collision breakpoints (coincident labels split and re-order —
essential at q ≥ 8 where labels repeat mod q). Verified exact: predictions/measurements
ratio 1.015–1.08, converging as O(1/W) (the tooth-width correction).

**First-regime closed form** (single-regime resonances, e.g. q = 7): the desert is
slack·(1/D₊ + 1/D₋), slack = 2/13 − 1/q, D± the minimal positive/negative binding
drifts of the label permutation. Consecutive clusters at p = ±1/7: drifts +1 (six gaps)
and −6 (wrap): length = (1/91)(1 + 1/6) = **1/78**.

**Worst-case table (exhaustive over difference patterns, the maximizer is ALWAYS the
consecutive pattern):** q=7: 1/78, 1/156, 1/130 for p = ±1, ±2, ±3; q=8..11 tabulated in
the .out (max 9/130 ≈ 0.069 at q=11, p=1). These are now EXACT inputs for the assembly.

## 2. THE FAST-PHASE FATTENING LEMMA (proved)

If the net condition fails throughout an interval of length 1/v_1 (an uncovered θ-arc
persists), the reference phase θ(t) = −v_1·t meets the arc within that interval: θ moves
at −v_1, offsets at +d_j ≥ 0, relative speed ≥ v_1, so the relative position sweeps a
full circle in ≤ 1/v_1. Hence the true desert is contained in the (net-ok ε-set)
fattened by intervals of transient failure each ≤ 1/W: **desert ⊆ formula-interval +
O(1/W)** — exactly the measured convergence. (The subtlety: the arc can close before θ
arrives — but then the net is restored and the formula-interval continues; only
PERSISTENT failure breaks deserts, and persistent failure is caught in one wrap.)

## 3. LOCALIZATION (the remaining lemma, status honest)

Claim: components > K/W occur only at census resonances (numerics: 100%, K ≈ 5).
Proved here: the PUMPING form — a covered interval of length ℓ with ℓW/2 exceeding the
handover-state count forces a simultaneous approximation |v_j·T − n_j| ≤ ε (state
pigeonhole on the deterministic handover walk; kps's pair-walk alternation is the c = 2
instance where the state argument is trivial). The pumping constant is astronomical
(state space ~ c!·(1/ε)^c); the SHARP constant needs the c-walk cascade (kps's offered
Stern–Brocot recursion) — the one remaining analytic ingredient, and it is the SAME
ingredient the fleet needs elsewhere: klein/kps's "equioscillation / AP-residue rigidity
uniform in r" and mac-mini's alignment-band confinement (THM-619's per-component residue
bands) are this lemma in other coordinates. One rigidity, three names.

## 4. THE CO-INCIDENCE CLOSURE (exact computation, all strata)

Per stratum |B| = #unlifted (citation 1/(|B|+1), window 2δ = 2(1/(|B|+1) − 1/13)/12,
danger zones = census resonances for c ≤ 12 − |B| with the worst-case radii of §1 +
window half-width + 1/134 fast-phase safety):

| stratum | cite | window | zones | PASS |
|---|---|---|---|---|
| B=5 (l=7) | 1/6 | 0.0150 | q=7 only | **792/792 — CLOSED** |
| B=4 | 1/5 | 0.0205 | q≤8 | 484/495 |
| B=3 | 1/4 | 0.0288 | q≤9 | 176/220 |
| B=2 | 1/3 | 0.0427 | q≤10 | 25/66 |
| B=1,0 | 1/2,1 | 0.07–0.15 | q≤12 | saturated |

**The original sliver (l = 7, five unlifted) is EMPTY** — every 5-base keeps a
margin-(1/6−) point clear of every possible q=7 desert — conditional only on §3's
localization. The B ≤ 4 failures over-grant the adversary: the table places ALL
resonances' deserts simultaneously, but ONE cluster has ONE difference pattern, and
simultaneous resonance at many q forces near-CONSECUTIVE patterns — which are the
S85-exact families (M = a/(2a+r−1) machinery, loose outright) up to perturbation. The
failing 109 shapes at B ≤ 4 are therefore the finite worklist for the per-family
simultaneity refinement, not open territory.

Pending refinements flagged: q ∈ {14, 21} short resonances for c = 7 (expected pass —
the Farey moat at those denominators is comparable); the B ≤ 4 simultaneity table.

## 5. The synthesis (what the repo history said)

THM-619 (mac-mini): each good-set component imposes a residue band on the killer —
the DUAL of the net condition (there the killer must align to the base's components;
here the cluster must align to the window). klein-S140: M < 2/25 forces AP residues —
the same rigidity in the loose branch. kps walks: the c = 2, 3 sharp cases. The S85
concurrent session: consecutive combs saturate M = v_min/(v_min+v_max) — the extremal
objects of ALL these statements are the same consecutive/AP configurations, met from
four directions. The missing analytic center is ONE lemma: **coverage/extremality
persisting over an interval forces AP structure at an explicit rate** — the c-walk
cascade (or equivalent) proves it; everything downstream is already exact tables.


## 6. S89 addendum: the degeneracy stratification (the correct scope of the closure)

The resonance census refines by the difference-gcd type g = gcd(d_2, …, d_c):

* **g = 1 (generic):** q = 14 needs 2 | all differences, q = 21 needs 3 | — so a
  generic cluster resonates ONLY at q = 7. The S87 792/792 computation is therefore
  the RIGOROUS closure of the generic type at |B| = 5 (its zone set was exactly right).
* **g ≥ 2 (degenerate):** deserts appear at p/(7g′)-families (verified: 2AP deserts at
  odd p/14, length 1/156 + O(1/W); 3AP at p/21, 1/234 + O(1/W); consecutive does NOT
  desert at 1/14 — measured 0.0015 = an ordinary short component). These clusters have
  D ≥ 6g (universal deserts ≤ 1/(13·6g)) and all values in one residue class mod g —
  strong family constraints. Zone-avoidance tables saturate here NOT because of deserts
  (negligible) but because 42+ windows of irreducible width 2δ cover the circle: the
  right tool is the DILATION/AFFINE TOWER (a g-degenerate cluster is a scaled copy of a
  smaller-scale cluster; if g divides every cluster value the fold t → gt drops the
  cluster to scale W/g = a smaller stratum of the same assembly), plus the same-parity
  family constraints. The degenerate types join the near-consecutive STRUCTURED
  WORKLIST: finite, recursive, S85-exact-adjacent.

Net: generic clusters CLOSED by tables (792/792); degenerate clusters = a finite
structured recursion. The B ≤ 4 saturation (S87) persists for the same window-geometry
reason and lands in the same structured worklist.
