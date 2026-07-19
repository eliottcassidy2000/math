# The LRC(14) campaign read as history: strategy families, the triage law, the wall count, and the tournament arc

**Instance:** opus-2026-07-19-S399 (owner-directed synthesis session)
**Mandate:** synthesize the full history of the repo's LRC(14) work; find common patterns and
mistakes, similarities and differences between strategies, and opportunities to leverage
big-picture understanding for critical insights and reframing of assumptions.

**Relation to the two prior owner-commissioned reviews — read those too, this does not repeat them:**
- `death-star-2026-07-18-S58` (SESSION-LOG:~1252): the five failure genera (sampling artifact,
  MAX-not-MEAN, dilation bites, raw-count-grows, sufficient-not-necessary).
- `07-reflections/the-lrc14-near-completion-history-and-the-certificate-rung-ladder-boxeph-S130.md`:
  the fourteen near-completion episodes, the mistake-percentage taxonomy, and the
  **certificate-rung-ladder mechanism** (every "closes 95%" tool is a bounded-denominator
  certificate; extremality *means* blocking or saturating every such certificate; hence the
  residual is one invariant object under seven names).

What is new here: (§1) the **full-arc era map** from first contact 2026-05-30 — the prior reviews
start at 07-08/07-11; (§2) a **strategy-family comparison table** across the whole campaign;
(§3) the **five-axis triage law**, a predictive instrument extracted from ~175 MISTAKES entries and
~25 strategy-level refutations, validated retrodictively against every major death; (§4) an
adjudication of **how many walls there actually are** (answer: two and a half, with the proof of
which named residuals are the same object); (§5) the **tournament arc** — the project-mandate
thread, from "transitive tournaments are unfaithful" to THM-1240's *forced nontransitivity*;
(§6) a cross-N data point verified this session; (§7) unabsorbed incoming ideas; (§8) reframings.

---

## 1. The full era map (2026-05-30 → 2026-07-19)

Fourteen phases. Sources: two full SESSION-LOG sweeps (72,690 lines, both write-regions — note the
log has a top prepend region AND a bottom append region used by codex/claude/monad-explorer;
neither alone is chronological), MISTAKES.md complete, hypotheses INDEX complete, all 17 LRC14
navigation docs.

| # | Dates | Phase | Lead(s) | Outcome |
|---|-------|-------|---------|---------|
| 0 | 05-30 | First contact: residue/endpoint-protection calculus, "finite anti-Bohr boundary theorem" | codex | Absorbed later; THM-357 trichotomy |
| 1 | 06-03..06 | Depth-distribution/resonance-lattice; **Vitali wall** (no finite-moment invariant decides tightness) | claude | Pinch oracle THM-369; tight census {AP, {1..11,13,24}} |
| 2 | 06-06..13 | Cross-problem lenses (LRC@19, Hadwiger-Nelson, band ladders) | claude, claudebox, kps | Framing only |
| 3 | 06-13..16 | **Lonely-measure L / singular series era** | kps, mac-mini | DIED by reframe: L has no floor (signed cancellation), L=0 necessary-not-sufficient |
| 4 | 06-16..17 | **The gap reframe** THM-523: object is M(S)=max min, counterexamples must be covering sets | mac-mini-S3 | The permanent backbone; "the disproof search PRODUCED the proof reduction" |
| 5 | 06-18..21 | μ_{1/7} witness floor ("the 1/7 pivot"), sector/Delsarte/wide-bound analytics | kps, mac-mini | Crux crystallized as coverage extremality |
| 6 | 06-22..28 | First Lean skeleton; carrier sprawl; **apex-7 frames** (Lee-Yang, Chebyshev, Beurling-Selberg, Borsuk-Ulam, cyclotomy) | kps, mac-mini, codex | Frames "DESCRIBE, do not prove"; binding obstruction is **2-adic not 7-adic**; polynomial method dies at composite 14=2·7 |
| 7 | 06-29..30 | **Covering-min crack via the observer lens**: three-gap mechanism, deep well {1..12,182}, 14/183=n/Φ₆ | klein arrives, mac-mini, opus | The first structural rigidity result |
| 8 | 07-01..02 | Consolidation; **Lean swarm** (8-module DAG); owner policy: LRC(≤13) settled as citation | fleet | Formal spine born |
| 9 | 07-03..04 | Far-peel, tight-locus formalization, confinement | kps, opus, mac-mini, klein | "LRC(14) ≤ LRC(13) + compressed case" |
| 10 | 07-05..06 | (G) spectrum-gap campaign; **Route 2 collapse** (MISTAKE-116/117: J–K citation was about accumulation points; finite covering system provably impossible) | opus-S130 audit | (C) survives as true math; the *bridge* was broken, not the statement |
| 11 | 07-07..08 | **Route-1 density-floor endgame**: all six per-k legs closed; THM-663 covering assembly; "one analytic item left" | mac-mini, klein, boxeph, death-star, monad arrive | First of the fourteen near-completion episodes at fleet scale |
| 12 | 07-09..13 | Grand assembly in Lean (foundational-axioms-only, modulo cite+hB5); covering-min rigidity proved; Ostrowski/Farey coordinates; FINISH-MAPs ("prove EITHER route") | fleet | Walk-backs begin: window unbounded, Q_s premise later refuted |
| 13 | 07-13..17 | Peel machinery; frame program; residue-six closure; certificate v2; B5 funnel proves CANONICAL_LONELY; 7-wall Lean chain | mac-mini, boxeph, death-star, opus | "Closed modulo a band" episodes; Q_s=Θ(r²) refutation kills the 07-13 premise |
| 14 | 07-17..19 | **The convergence**: hard stratum q0>14; tight locus proved = exactly TWO families; r=3,4,5,6 clustered strata closed; sojourn 2/21 (THM-1203); n=12 Hamming program to c=35; **(D,s)/Farey coordinates** (slack, D=M·s); **mod-p spread ladder** (13/19/23); comb/blocker-cycle package (THM-1197..1250); kernel splits (S58f-i); three history audits | everyone | The residual named precisely; see §4 |

Two structural observations the phase table supports:

**(a) The two great reframes were both object-corrections, not technique-upgrades.** Phase 4
(L → M: mean-type → max-type) and phase 7 (family statistics → observer-lens/three-gap at the
marked origin). Each unlocked more than every technique refinement combined. The current
candidate for a third object-correction is §5's (the cycle-bearing owner/blocker graph as
first-class object).

**(b) Every "one item left" was declared at a phase boundary** — 06-28 ("exactly ONE new theorem"),
07-08 ("one analytic item"), 07-09 ("sole obligation hB5"), 07-13 ("one inequality, prove either
form"), 07-16 ("two named uniformities"). boxeph-S130 chronicles fourteen of these. The
arithmetic reason ("the sliver": the free pigeonhole floor 2/29 is 96.6% of target, so every
certificate family gets close before stranding) is boxeph's §4; the *sociological* reason visible
in the full arc is that **phase boundaries are where syntheses are written, and syntheses
compress the residual into the current frame's vocabulary** — the compression IS the overclaim.
Rule 3 of boxeph-S130 (type the residual) is the correct fix and should be adopted verbatim.

---

## 2. Strategy families: what each was, proved, and where each stalled

The campaign ran ~13 distinguishable strategy families. Table columns: what the method bounds
(mean/measure/max/exact), which symmetries it respects (T=translation, D=dilation), certificate
class (single-modulus SM / cross-modulus XM / analytic AN / exact-enumeration EX), headline
result, stall mechanism.

| Family | Bounds | Sym | Class | Proved | Stalled on |
|---|---|---|---|---|---|
| 1. Lonely measure L / singular series | mean | T | AN | — | L has no floor; signed cancellation; wrong object (phase-3 death) |
| 2. Depth/moment invariants | mean | T,D | AN | Σk·p_k=2n/(n+1) | **Vitali wall**: all finite moments blind to tightness |
| 3. q-witness sieve / bounded-denom certificates | max | D | SM | non-covering case ENTIRELY (+ tight locus lives here) | rung ladder: blocked by divisibility engineering; covering modulus unbounded (MISTAKE-116) |
| 4. Density floor μ_{1/7} (per-k legs, tent/window/energy/PZ) | measure | D (after primitivity peel) | AN | all six legs k=8..13; m_P floor; THM-663 | realization/finite-Vmax glue; means-vs-tails corrections en route |
| 5. Fourier / soft-Weyl / Q_s / Delsarte | mean+abs | T,D | AN | Q_s=O(r) at k=13 core; exact means/variances | absolute bounds lose to **signed cancellation**; Q_s=Θ(r²) at the resonant peel; Delsarte in-principle impossible (lonely set has measure ZERO at tight families, THM-1185) |
| 6. Freiman/additive (doubling, energy, BSG/PFR) | invariant | **T**,D | AN | E₃ inverse THM-730 | **translation-invariance blindness** (THM-1225); "M<1/13 supplies NO additive energy — the crux is the Diophantine→energy bridge" (boxeph-S104) |
| 7. Covering-min / observer lens | max | D | XM | deep well unique 14/183 (THM-724/726/883); lowness lemma | global form was the subject of the one substantive court case; feeds but does not close the inverse theorem |
| 8. Strata/horn decompositions (r=2..6, carrier charts) | max | D | EX+SM | r=3 (THM-1094), r=4 (THM-1097), clustered r=5 (THM-1214), r=6 strata (THM-1212), five-comb (THM-1198) | averaging step provably halts at r<3.5 (THM-1140); each closure re-spawns a renamed residual (HYP-7744→46→48→50 chain) |
| 9. (D,s)/Farey/slack coordinates | exact | D | EX | M=D/s pinch; slack=nD−s∈ℤ; slack-0⟺extremal; D≥4 in (1/14,3/41); D=M·s (THM-1269, ex-1245); tight locus = 2 families (THM-1120/1142) | "bound D" ⟺ bounded primitive speeds — which is the inverse theorem again (§4) |
| 10. Mod-p spread ladder (13, 19, 23) | max | **sensitive to T** — correct side of the triage | SM (necessary conds) | antipodal-spread lemmas, kernel-pure; ledger rung `gap_regime_mod19_spread`; mod-23 near-bijection pin (HYP-7880) | necessary-not-sufficient; "the obstruction is irreducibly cross-modulus" (boxeph-S126) |
| 11. Comb/slow-gap owner-incidence program | measure+exact | D | AN+EX | five-comb dual density; 2/21 sojourn (THM-1203); blocker-cycle + address compression (THM-1233..1250) | the **oriented germ lift / handoff-debt** lemma — a transport statement, categorically NOT a bounded certificate (§4) |
| 12. Lean formalization ladder | — | — | — | 569 LRC modules; grand assembly foundational-axioms-only; **the only asset class never walked back** | kernel-pure ≠ non-vacuous (MISTAKE-136/146/154/155/186) |
| 13. Tournament/observer bridges | order data | — | — | phase-clock walk theorems; Paley bridge THM-640; source-marked classes (THM-381) | transitive projections erase sign/phase/mass — until THM-1240; see §5 |

**Similarities across strategies the table makes visible:**
- Families 1, 2, 5 all died the same death at different altitudes: an **averaging step**
  (mean over t, over moments, over relation supports) that integrates across the certificate
  rungs and therefore cannot see the single saturated rung carrying the extremal family.
  kind-pasteur's THM-1140 later made this *exact* for family 8: the union/averaging step
  requires 7−r > r, halting at r<3.5 — the one place the campaign proved *why* an averaging
  method halts rather than just observing that it did.
- Families 3, 8, 10 are all rung machines (boxeph-S130 §4); their closures are real and
  permanent (the entire non-covering case, all clustered strata) but their residual is invariant.
- Families 7, 9, 11 are the three that produced **rigidity** statements (unique deep well;
  two-family tight locus; forced blocker cycle) — all three are exact/XM methods that broke
  translation invariance and worked on primitive orbits. This is the triage law working (§3).

**Differences that mattered:**
- Family 4 vs 5: both analytic, but the density floor bounds a *measure of a good set* (a tail
  event) while soft-Weyl bounds a *mean of a signed sum*. The first closed its legs; the second
  was structurally refuted. Tail vs mean, again.
- Family 9 vs 3: both Diophantine, but (D,s) is an *exact identity frame* (every attained M has
  the form, no bound claimed), so it cannot be "blocked" — it converts questions instead of
  answering them. Its danger is the opposite: conversions can be mistaken for progress
  (episode #13 in boxeph's chronicle).
- Family 11 is the only current line whose frontier statement is *not* a bounded-denominator
  certificate and *not* a global-uniformity claim — it is a local-to-global transport/descent
  lemma. This is why §4 counts it as a genuinely separate wall.

---

## 3. The five-axis triage law (the predictive instrument)

Extracted from all ~175 MISTAKES entries, the ~25 strategy-level refutations, and the three
rigidity successes. **A method can in principle decide near-floor structure at n=14 only if it
satisfies all five:**

1. **Breaks translation invariance.** The floor is not translation-invariant (translate ladder
   {1+k,..,13+k}: doubling fixed, M runs 1/14→5/22 — THM-1225). Kills: doubling, additive
   energy, Freiman dimension, direction-set/Kakeya invariants, any conv-hull/difference-set
   statistic. (Positive examples: q0, residues mod p, D.)
2. **Respects dilation invariance — i.e., works on primitive orbits.** M is
   dilation-invariant; covering, slack, and speed-size are not. Any claim quantified over raw
   coordinates dies by dilation transport (MISTAKE-073/075/091/126/128/140/156/170 chain —
   "third time dilation has bitten"). Positive form: state everything on primitive families;
   normalize gcd FIRST, then rederive non-invariant predicates (covering!) afterward.
3. **Bounds a max/tail, not a mean.** Existence of a witness is a sup fact. E[maxgap] is not
   AP-minimized; fine-scale truncated means are vacuous; L had no floor; the Vitali wall says
   ALL finite moments are blind. (MISTAKE-129 "existence is a MAX"; death-star HYP-4777.)
4. **Tolerates signed cancellation — tracks phase, or is exact.** Absolute-value envelopes
   diverge or lose: certificate v4 (near-orthogonality invisible to absolute bounds), soft-Weyl
   13-fold product, |S|=Θ(R), PZ-on-the-moat. The successes are exact identities (D=M·s,
   folded identity THM-965, additive-triangle deletion THM-1203) or explicitly signed objects
   (mirror-pair parity THM-953).
5. **Is cross-modulus adaptive or exactly enumerative — never a fixed finite certificate
   family.** Any fixed modulus set/window/denominator bound is defeated by CRT/lcm engineering
   (MISTAKE-116, THM-1105/1110/1115, escape families ≡ AP mod lcm(2..Q₀) clearing only at
   nextprime(Q₀)). What survives: complete finite strata enumerated exactly (r≤6 charts,
   Hamming banks), or statements whose modulus adapts with the family.

**Retrodictive validation** (the law's content is that these were *predictable*): every entry in
the campaign's refuted-routes lists violates a named axis — counting lemma (4), Q_s=o(r²) (3,4),
uniform Bonferroni ledger (4: sign flips within a residue class), Delsarte LP (3: measure-zero
lonely sets at tight families), doubling/Freiman detector (1), n=13 rigidity transfer (not an
axis violation — an n-specific fact, THM-1220, the one genuine surprise), finite covering
systems (5), E[maxgap] program (3), sampled-margin uniform claims (5, in the instrument rather
than the math). The three rigidity successes (deep well, two-family tight locus, blocker cycle)
satisfy all five. I could not find a counterexample to the law in either direction in the corpus.
**Proposed practice:** before any new route, write one line per axis. It is a two-minute check
that would have pre-killed most of the fourteen episodes.

(Relation to boxeph-S130 §4: the rung ladder explains why axis-5 violations *specifically* fail;
the triage law is the union of that mechanism with the four analytic axes.)

---

## 4. How many walls? Two and a half.

The final-week documents name at least eight "residuals". The equivalence web, with the canon
that proves each identification:

**Wall A — the inverse-theorem complex (one object, seven names):**
- n=12 AP-uniqueness (HYP-7310, "Tao's optimistic conjecture") = the open half of boxeph
  THM-1017 (LRC(14) ⟺ [M<1/13 covering ⟹ AP core]) — by THM-1017 itself.
- death-star's "5% rational-time-evasive Freiman core" (HYP-7750/7895) — *named to be*
  HYP-4382 = the 12-set-uniqueness core, in S58i's own words.
- opus's "bound D" (THM-1268/1269, formerly opus-numbered 1240/1245 — renumbered this session
  after first-push collisions with codex's THM-1240 and kind-pasteur's THM-1245): D=M·s makes
  bounding D ⟺ bounding the active-pair sum ⟺
  bounded primitive speeds near the floor — which is exactly "near-floor primitive families are
  structured" = the inverse statement. The rung equations (D≥4 ⟹ s=55, 69, 83, …) are Wall A
  in (D,s) coordinates, not a separate finite problem: residue feasibility at any fixed rung
  never sees the global maximizer (opus-S397 "the residues are never the obstruction"), so the
  decision is analytic — the axis-5 lesson one level up.
- the certificate-blocking class (boxeph-S130 §4), the compact/non-dilated/alignment/escape-tail
  cores — identified there.
- Depth-minimal numeric shadows: 3/38 (n=12 gap) and 4/55 (n=14 interval).

**Wall B — the phase-transport complex (codex's comb program):** the oriented private-germ
lift around the finite blocker word / first-failed-continuation handoff-debt alternative
(PROOF-MAP top banner), with its seven-wall sibling (the connected G_gt branch with phase
retained, HYP-7870). This is *not* a bounded-denominator certificate and *not* a global
uniformity claim; it is a local-to-global transport with an explicit finite address alphabet
(THM-1248's digit words, contraction 1−A > 4691/5503716, THM-1250's lcm debt). The guardrail
theorems (THM-1239 resonant blocker, THM-1242 mixed q=15 sunflower) show the adversary still
adapts — but along *phase*, not modulus. Whether Wall B secretly reduces to Wall A is **open
and important**: death-star-S58 cont.9 found the r=6 tail resists three independent method
classes "confirming from three sides that the residual is genuine low-complexity arithmetic
rigidity, the same wall the fleet's inverse-theorem route has circled" — suggestive of identity;
but no theorem identifies them, and the comb program's objects (owner words, holonomy, lcm
debts) have no current translation into family-structure statements. **Adjudication: treat as
distinct until a reduction is written.** Consequence: the fleet has two independent shots, and
should protect Wall-B work from being re-aimed at Wall A by synthesis pressure.

**The half — n=12 sporadic emptiness** (the H5/H6 banks, next face c=36, plus the "metric
reverse-content law"): explicitly "must not be counted as closed by the direct LRC(14)
progress" (S78 banner). It is inverse-flavored but runs on its own finite-face ladder with its
own closure mechanism; on current evidence it is a *sub-wall of A* that has a working
assembly line, hence "half".

**What this means operationally:** (i) any new "residual" named in a future session should be
typed as A, B, or A-half — or justified as genuinely new (extends boxeph rule 3); (ii) Wall A
progress claims must state which of its seven coordinates they advance and whether the advance
survives translation to the others; (iii) Wall B has the only frontier lemma that is
*falsifiable in finite terms* right now ("lift the germ or pay new lcm debt") and deserves
sustained specialist attention rather than fleet-wide swarming — swarms produce rung machines.

---

## 5. The tournament arc (the project-mandate thread)

The repo's mandate is parity in tournaments; LRC entered through tournament analysis
(observer-source tournament THM-381: *the observer is lonely iff the marked vertex is a
source*; phase clocks THM-373; pressure lifts THM-377). Then, across phases 12–14, at least six
independent verdicts of the form "the faithful object is NOT a runner tournament" (codex-S73,
S75 update, THM-1219, THM-1226, UNIFIED-SYNTHESIS §0, TARGETS doc). It would be easy to read
the arc as "the tournament lens failed."

That reading is wrong, and the precise statement matters:

- Every unfaithfulness verdict concerns the **transitive order tournament** (runner order,
  determinant-sign order). Transitive = an order = zero cycle structure = no information beyond
  a permutation. The verdicts say: the proof-bearing data (sign, phase, mass, incidence) lives
  in what the transitive projection quotients away.
- The campaign's own endgame then *forced nontransitivity*: THM-379 (any nonempty
  owner-compatible endpoint core has a directed owner cycle), and now THM-1240 — every
  six-cover induces a blocker functional graph with a directed cycle of length 2..6, "the
  first genuinely nontransitive tournament forced by the hard slow-gap geometry" — with
  strictly positive product holonomy in the canonical gauge (S82 banner). The Wall-B frontier
  lemma is literally a statement about transporting structure around a directed cycle.
- So the arc is: **the tournament lens returns exactly when cycles appear, and the obstruction
  is carried by cycle holonomy** — which is the repo's core theme (Rédei/OCF: the interesting
  invariant lives on odd cycles, not on orders) resurfacing in LRC clothing. The observer-lens
  memory ("structures seen via the marked origin") cracked covering-min in phase 7; the
  cycle-holonomy lens is its successor object.
- Concrete, non-decorative next step: the blocker cycles have length 2..6; the holonomy
  positivity is a signed/parity statement on cycles. The repo owns exact technology for signed
  invariants on tournament cycles (Rédei involution transplant LEM-020, mirror-pair parity
  THM-953). Whether an involution-pairing argument can force the germ lift around odd blocker
  cycles is a well-posed question that only this repo would ask. Filed in the backlog rather
  than oversold here.

Also on-mandate: THM-640's Paley bridge ("the density-floor M-minimizer is the H-maximizer";
"composite-14 is why LRC(14) is hard in tournament language") and HYP-3805 (Paley heptagon as
LRC extremal object) remain the two cleanest tournament↔LRC structural identifications; both
are diagnostic, neither is load-bearing yet.

---

## 6. A cross-N structural fact (verified this session)

Script + frozen output: `04-computation/lrc14_dms_identity_crosscheck_opus_S399.py`,
`05-knowledge/results/lrc14_dms_identity_crosscheck_opus_S399.out`.

death-star-S59b discovered `{1..29,31,120}` at M=4/127 (N=32 runners) this morning. Verified
here: the maximizer is t*=55/127, the active straddling pair is (7,120), s=127, **D=M·s=4
exactly, slack n·D−s = 1** — the same slack-1 structure as {1..11,13,36} (pair (5,36), s=41,
D=3, slack 1) and the slack-0 structure of the extremals. Three consequences for the synthesis:

1. The (D,s)/slack frame is **structural across N**, not n=14 numerology.
2. **Rung realizability is N-dependent**: the D=4 slack-1 rung is realized at N=32 but has no
   known realization at N=14 (the 4/55 target). So "(1/14, 3/41) is empty" — if true — is an
   N-specific fact needing an N-specific mechanism, consistent with THM-1220 (n=14 the unique
   non-rigid n in 10..18 via 12↔24). Cross-N census (death-star's S59 line) is the right
   instrument and is genuinely underused.
3. The near-floor families across N so far are all "AP + one or two engineered defects with a
   small-slack active pair" — the inverse-theorem shape, again, now with cross-N evidence.

---

## 7. Unabsorbed incoming ideas (from the inbox deep read, MSG-1599–1850, 07-14 → 07-19)

Coverage note: MSG-1665–1850 read in full (via two digests + direct read of 1790–1850);
MSG-1599–1664 covered only via filenames and later cross-references — that ~92-file window
(klein S298–S312, mac-mini S95–S108, codex THM-773–830 arc) is the one honest gap.

**The headline — the single highest-leverage unworked direction in the whole window:**

1. **Mine the Sungkawichai–Trakulthongchai LRC(11–13) proofs for their equality case**
   (boxeph-S114, MSG-1823, proposed as candidate (2) for discharging INV): *"if the proof pins
   M(C)=1/13 ⟺ AP, HYP-4382 follows and the residual collapses — the most plausible route,
   living inside a proof not in this repo."* Grep confirms ZERO mentions afterward. Since
   LRC(≤13) is settled by owner directive and cited freely, the equality-case *structure* of
   those proofs is fair game — this converts Wall A into a literature-reading task. If their
   sieve pins the tight locus, the n=12 AP-uniqueness may already be implicit in a paper the
   repo already trusts. One session: fetch, read the equality analysis, report exactly what
   their argument pins.

**The rest, ranked:**

2. **The mod-23 near-bijection pin** (boxeph-S130, HYP-7880): 12-sets in the uniqueness gap
   are mod-23 near-bijections (slack 1). Canon contains 2/23 only as an attained value;
   nobody has used p=23 as a constraint. One Lean file by the LRCMod19Spread template.
   Axis-5 caveat: it carves rungs, it cannot decide them. Pending handoff: boxeph→opus asks
   for the mod-17/19/23 CRT stack on 4/55 realizers *before any further sampling*;
   death-star (MSG-1850) sharpens the aim to non-single-far D=4 realizers.
3. **The metagraph transport** (boxeph-S110, MSG-1818, follow-up (b)): 183=|PG(2,13)|; the
   deep-well AP is the *transitive pole* and the Singer difference set the *regular pole* of
   the same object. Make the LRC-config → tournament-iso-class map precise so G_n's proved
   transitive-class isolation applies to M. The only place in the window where the tournament
   project and LRC(14) actually met — feeds §5 directly. Untouched.
4. **The blocking-budget count** (boxeph-S130 §7): 13 speeds service ~15 simultaneous
   certificate constraints; prove a constraints-per-speed counting bound forcing near-AP
   structure. The Diophantine→additive bridge in finite form. Untouched.
5. **The signed cross-orthogonality estimate** (boxeph-S49): the named successor to the
   certificate line — "absolute bounds end where signed cancellation begins." §5 resonance:
   the repo's parity technology is the natural toolkit. Unattempted.
6. **Function-field port of deficit→margin** (death-star-S58e/g, MSG-1842/1845): in the
   function-field model the archimedean carry vanishes — the natural place for a clean
   Schur-deficit→margin inequality. Distinct from boxeph's abandoned S91 FF attempt (that
   targeted difference-closure). Not run.
7. **Height-bound citation frame** (mac-mini-S120, MSG-1799): Erdős/Jacobsthal tight families
   with v_max=2n−Θ(log n), plus Pomerance's "n<v_max<2n−c·log²n ⟹ not tight" — sharpens
   HYP-7450 from the literature side. Never referenced again.
8. **Cross-N D-graded gate census** (death-star-S59b, MSG-1850, fresh): extract the full D=4
   gate, test the predicted next opening (N=61?), prove the first-gap side as a *band*
   phenomenon. Converges with §6.
9. **The hybrid tree×moment bound** (kind-pasteur c66/c67): pairwise-only moment LP is
   *worse* than the spanning tree — the bounds are non-nested; a hybrid may beat both. Its
   r=6 motivation evaporated, but the inequality question stands alone. Unworked.
10. **The formalization manifest as idle assembly line** (MANIFEST-2026-07-17): ten
    decide-shaped one-session items + mac-mini's three built rungs (FragmentationCount
    repaired after MISTAKE-138, TieSplitWalk, KillerBudget). Unclaimed rungs whenever no
    agent picks one up.
11. **Engineering-mandate dormancy** (standing violation of CLAUDE.md's equal-priority
    mandate): boxeph recommended the engineering pivot three times (MSG-1810/1811/1823);
    grep shows no engineering-deliverable letters since March. Concrete revival targets in
    §8.7.
12. **PROOF-MAP staleness** (boxeph-S130 warning, still true; handoff unassigned): the map
    has not absorbed the mod-19 rung, the clustered floor, the rational-time floor, or the
    D-bound reduction. The MISTAKE-183 staleness channel stays open until an editor folds
    them in.

**Pending handoffs answered this session (opus's mail):** kind-pasteur c84's request to
re-scope the opus branch-2 write-ups is now DONE — THM-1215 carries the detection-floor
banner (the search negatives are vacuous below 1/14; the 1200/1200 grid certificates are the
citable evidence). The THM-1245 and THM-1240 first-push collisions are resolved: opus's files
renumber to THM-1269 and THM-1268. boxeph's CRT-stack request: acknowledged, deliberately not
run this session (§9); it belongs to whoever works the 4/55 rung next, with the §4 Wall-A
caveat attached.

---

## 8. Reframing opportunities, ranked (honest status on each)

1. **Adopt the five-axis triage as a standing gate** (§3). Cost: one paragraph per proposed
   route. Would have pre-killed most recorded deaths. Status: meta-practice, adoptable today.
2. **Type every residual against the wall census** (§4, extending boxeph rule 3): A / B /
   A-half / justified-new. Makes "seven names, one object" structurally impossible to repeat.
3. **Protect Wall B as an independent program.** It is the only frontier lemma that is not a
   rung machine; its failure mode (phase-adaptive adversary) is different in kind; do not
   re-aim comb sessions at Wall A vocabulary. Status: allocation decision for the owner.
4. **Exactification bounty**: every measured constant in canon gets converted to an exact
   rational or an identity (the S128c85 model: 0.164 → 31807/194040; the hypotheses sweep
   found "an identity that explains a frozen constant essentially never dies"). Cheap,
   compounding, and each exactification is Lean-ready.
5. **Cross-N before n=14-specific bets** (§6): the rung-realizability map over N is finite
   work per N and converts "is (1/14,3/41) empty?" from a single hard question into an
   interpolation question with data. Status: death-star has started; needs a dedicated census.
6. **The cycle-holonomy lens** (§5): attempt the germ-lift lemma with involution/parity
   technology on odd blocker cycles. Status: speculative, on-mandate, well-posed; one
   exploratory session to determine if the pairing even typechecks.
7. **Engineering deliverables from the campaign's instruments** (dual mandate): (a) the
   **instrument validation gate** (any search must rediscover known extremals from inside its
   sampling range before its negatives count — kind-pasteur cont.84's practice, generalizes to
   any computational-math QA); (b) the exact-rational sweep + straddling-pair analyzer (this
   session's script is a 60-line seed); (c) the certificate-rung profiler (boxeph-S130's
   script); (d) the Lean mod-p spread template (LRCMod19Spread → parametrize by p). These are
   concrete `04-computation` library targets aligned with the mod_rank/PyPI roadmap.

## 9. What this session did and did not do (scope honesty)

Did, beyond the document: the §6 verification (script + frozen out); the THM-1215
detection-floor re-scope banner (accepting kind-pasteur's audit — the session's one canon
correction, of opus's own prior work); the THM-1240/1245 → THM-1268/1269 renumbers (first-push
precedent, opus loses both); backlog entries for §7's leads.

Did NOT: no new mathematics beyond §6; no new HYP claimed (all factual claims sourced to
existing canon). The m=1 CRT feasibility check on the 4/55 rung was considered and deliberately
NOT run — the §4 Wall-A analysis predicts it is feasible-and-undecisive as a decision
procedure, and running it anyway would have been the axis-5 mistake with fresh paint
(episode #15). As a *search-space carver* for future 4/55 work it retains value — that is how
boxeph's pending request should be consumed, by whoever works the rung next.

## 10. Cross-links

death-star-S58 synthesis (SESSION-LOG:~1252) · boxeph-S130 reflection + HYP-7880 ·
THM-523/369/381/379/373/377/640 · THM-663/671 · THM-724/726/883 · THM-1002/1008/1017/1028/1043 ·
THM-1094/1097/1120/1140/1142/1149/1158/1171/1185/1198/1203/1210/1212/1214/1220/1225 ·
THM-1230/1235 · THM-1268/1269 (ex-opus-1240/1245) · codex THM-1240/1248/1250 · kind-pasteur
THM-1245 (witness law) · death-star THM-1284/1272 (cross-N) · HYP-7310/4382/7750/7870/7880 ·
MISTAKE-116/117/129/136/
138/140/146/154/155/156/160/161/162/166/170/173/174/182/183/186 · LRC14-PROOF-MAP.md ·
LRC14-UNIFIED-FRONTIER-SYNTHESIS-2026-07-18 · LRC14-FORMALIZATION-MANIFEST-2026-07-17 ·
`04-computation/lrc14_dms_identity_crosscheck_opus_S399.py` (+ frozen .out).
