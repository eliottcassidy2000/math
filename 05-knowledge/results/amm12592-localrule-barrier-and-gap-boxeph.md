# AMM 12592 — Lane E3: the local-rule barrier and the existence gap (the honest golden statement)

Session: boxeph multifront, 2026-08-04 (after D1 bulk sweep, D2 E-lin, D3
golden bound, and the Estimate-E hostile referee; all citations T1–T9',
Theorems A/B/C/S, G1–G3, THM-3329 to those notes). All computations exact
int/Fraction; no floats in any decision; no numpy; sympy unused.

Script: `04-computation/amm12592_localrule_barrier_gap_boxeph.py`
(stages A record+transportation point, B lemma certificates, C class sweep,
D schedule beam). Outputs:
`05-knowledge/results/amm12592_localrule_barrier_gap_boxeph.{out,json}`.

Status labels: **PROVED** (complete argument), **VERIFIED-exact** (exact
computation at stated scales), **OBSERVATION** (rigorous point with stated
scope caveat), **CONJECTURED / OPEN** (precise statement, open).

---

## 0. Headline

1. **RECORD CORRECTION (the gap has moved).** The task premise
   "exact-floor feasibility at R >= 256 with D0 = 0 is OPEN; D0*(256) = 1
   for every method tried" is stale. The D1 sweep closed **R = 256 at
   D0 = 0** (rule desc1) on 2026-08-03; the witness
   (`amm12592_witness_R256_bulk_desc1_D0_0_boxeph.json`, sha256
   5950bd42…) passed the D1 referee AND the Estimate-E hostile referee's
   fresh streaming verifier, and was re-verified again this session
   (stage A: sha + verify_witness ALL-PASS). "D0*(256) = 1" was a
   statement about plain rule A, not about feasibility. **The exact-floor
   existence frontier is now R = 512, D0 = 0 (OPEN, with D0 = 1, 2, 3
   also open; best known at 512 is D0 = 4, tr15, referee-verified).**
2. **NEW CERTIFICATE: transportation integer point at (256, 0).** The
   desc1 witness translates exactly to the nonnegative transportation
   form `f_{i,t} = (C(d_i,t) - delta_{i,t})/2`: all 58837 cells integer in
   `[0, C(d_i,t)]` (8.45% boundary-saturated), and the (**) identity
   `sum_i x^i F_i = T_R` holds at x = 2, 3, -1 exactly (stage A; with
   block admissibility + the epoch identity this is an integer point of
   the THM-2966/3002 polytope). So the transportation LP relaxation at
   the exact floor profile is feasible at R = 256 **with an integer
   point** — the polytope chain 8..128 (LP note) extends to 256.
3. **Formalization finding (OBSERVATION, load-bearing).** "Clamp-local
   rule" cannot be made a robust complexity class by window-locality
   alone: a width-2 rule that reads its local cap values can reconstruct
   (d, t) — and generically (R, i) — and thereby replay an arbitrary
   precomputed global schedule (Section 2.2). A provable barrier must
   therefore quantify over either (a) a finitely-parameterized concrete
   rule space (Section 2.3), or (b) reachable STATES (the G3 mechanism
   route, Section 2.5). Both precise conjectures are stated; neither is
   proved; the evidence ledger is assembled.
4. **New PROVED lemmas (small but structural).**
   (L1) crossover/front complementarity `(t_2 - 1) + F(0) = d_0`;
   (L2) front marginality: the dominant front-cell load equals HALF the
   c-box width exactly, `C(R-2-F_0, d-F_0) = C(d, F_0)`;
   (L3) boundary-layer rate law
   `(A_{t-1}/A_t) / (C(d,t-1)/C(d,t)) = (R-1-t)/t` → `g/(1-g) =
   1.48748…` per cell of depth (closed form for D1's measured "~x1.5");
   (L4) the influence cone: junk/choices at cells in `[c-2k, c]` are the
   ONLY influences on cell c after k rows — feed cannot touch the front
   region for `~0.2R` rows, so the nucleation battle is fought entirely
   inside an `o(R)`-cell boundary layer (proves the causal content of
   B1's steering invariance and D1's window-invariance measurements).
   Plus the new limit constant: front-cell overflow/box-width →
   `kappa = (3-4g)/(2(1-g)) - (1-g)/g in (0.0839816, 0.0839819)` at
   D0 = 0 (VERIFIED-exact table 128..8192; limit PROVED given L2/L3).
5. **Bounded exact searches at the frontier (hazard-labeled).** All 21
   D1-class rules die at (512, D0 = 0..3), 84/84 runs — die rows
   107..174 (asc degenerately at 1). A beam
   over cross-row action schedules (12-action alphabet, protected pure
   lanes, G3-order-parameter scoring) closes the (256, 0) control at
   capture 153 but ALL-DIES at (512, D0 = 0..3) within its width/budget.
   **Search negatives prove nothing about feasibility** (the record
   itself shows why: beams and LP failed at (256,0) where desc1 later
   closed). (512, 0) remains genuinely OPEN in both directions.
6. **Synthesis (Section 4): three nested programs to C\* <= 1 + gamma\***
   with exact status of every route, plus the floor-side honesty note
   (post-THM-3024-demotion the general-class lower bound is open, so the
   golden VALUE is currently pinned only for balanced-block schemes).

---

## 1. The record, corrected (exact-floor feasibility ledger)

Exact-floor = profile `d_i = floor(g(R+i))`, D0 = 0, g = gamma*.

| R | exact-floor status | constructions | verification |
|------|------|------|------|
| 8..64 | CLOSED | direct beam / slim witnesses; rule constructions | referee'd (multiple sessions) |
| 128 | CLOSED | three+ independent: direct beam + banded repair; LP-parity-rounding v5; parameter-free rule (family-hunt); plain rule A itself (capture 81) | all referee'd |
| 256 | **CLOSED** (2026-08-03, D1) | desc1 (descending absorb-window + depth-1 vent) | witness sha 5950bd42…; D1 referee ALL-PASS; Estimate-E fresh referee ALL-PASS; **this session: re-verified + transportation integer point (stage A)** |
| 512 | **OPEN at D0 = 0,1,2,3**; CLOSED at D0 = 4 (tr15) | 21 rules die at D0<=3 (stage C); schedule beam all-dies (stage D) | tr15 D0=4 witness sha 311975f6…, twice referee'd |
| >=1024 | open below the D1 curve (14/36/87/<=188 …) | — | — |

Notes. (i) The 256 closure is NOT a plain-rule fact: plain rule A dies at
(256, 0) at row 61; the vent freedom (junk injection at in-box cells,
licensed by the T2 generalized-clamp conjugacy) is what closes it. (ii)
The task-brief's "LP/beam did not close 256 at D0 = 0" is true of those
methods; the false step was reading method thresholds as a feasibility
threshold. (iii) Vent-shape sensitivity (D1 sec. 4: depth-2 vent FLIPS
the (256,0) closure to death) means the frontier is razor-thin in rule
space; nothing monotone can be inferred from rule failures.

## 2. Part (1): the clamp-local barrier

### 2.1 The flow and the full class (recalled, PROVED elsewhere)

In T2/T6 junk coordinates, a continuation is a choice, at every row i and
cell t, of an even in-box value `u_t in [-2C(d-1,t), +2C(d-1,t-1)]`
(junk `jn_t = w_t - u_t`; death iff `jn_d != 0`; capture iff `jn = {}`
post-feed). By the T2 bijection, **the set of surviving-to-capture choice
sequences IS the set of epoch witnesses at the profile.** Hence:

**No-full-class-barrier principle.** Any "barrier" proved for ALL
admissible continuations at slack `< cR` is an epoch-INFEASIBILITY
theorem. Since a (256, 0) witness exists, no such barrier holds at
R = 256; and if the golden program (GG + S') is true, none holds at slack
o(R) for any large R. A local-rule barrier is therefore necessarily a
statement about a RESTRICTED class, and the restriction carries all the
content.

### 2.2 Why window-locality alone is vacuous (OBSERVATION)

Say a rule has window W and lookahead L if its choice at (i, t) is a
fixed function phi of the loads/junk in cells `[t-W, t+W]`, the cap
values in that window, and the Beatty letters `delta_i..delta_{i+L-1}`,
the same phi for every R (R-uniformity).

*Observation (near-collapse of the class).* From two adjacent cap values
`2C(d-1,t), 2C(d-1,t+1)` in its window, phi recovers (d, t) exactly
(binomials are strictly monotone in n at fixed k; the ratio pins k
against n). From d and the slack convention D0(R), the row index
n = R + i solves `floor(g n) = d - D0(R)` with at most 2 solutions
(Beatty preimage), and dyadic R is the unique power of two in (n/2, n];
residual ambiguity (rare collisions across epochs) would additionally
have to match the entire load window, which contains epoch-specific
~Theta(R)-bit integers. So — up to input collisions never observed in
any recorded run — a width-2 rule can identify (R, i, t), reconstruct
the (deterministic) global state offline, and replay ANY precomputed
global schedule. Window-locality without information- or
description-restrictions is thus not a smaller class in any provable
sense. (This is why the barrier is formalized below over (a) a concrete
parameterized rule space and (b) reachable states — not over "all local
functions".)

### 2.3 The concrete rule space and the parameterized barrier conjecture

**Definition (the D1 rule space RS).** A rule in RS is an instance of the
`choose_row` schema of the D1 sweep engine: parameters
`(orientation in {desc, asc}; vent depth V >= 0; vent truncation vlad in
{1, 2, ..., inf}; window Wwin in {0, 1, ..., inf}; guard (K, rho0 in Q>1)
or none; fill in {0,1})`, acting per row by descending (or ascending)
pass with exact one-row kernel lookahead, in-box absorption targeting,
vent pre-cancellation ladder, and optional alternation guard. Every D1
variant (plain, desc, desc1, descnv, descw*, train*, tr*, fill*) is an
instance; RS is the class for which exact D0* data exists.

**Conjecture E3-BC(param) (the barrier, parameterized form).** There is
a universal `c0 > 0` (data suggest `c0 ~ 0.02`) such that every fixed
rule p in RS has `liminf_R D0*_p(R)/R >= c0`; i.e. no rule in RS closes
dyadic epochs at sublinear slack. CONJECTURED, not proved.

**Evidence (all exact, die/close-pinned, referee'd at the records).**
- plain: D0* = 0/1/5/15/38/89/192/[401,416] at R = 128..16384,
  normalized → ~0.0245 with declining-increment linear signature.
- best of family (desc1/tr15): 0/0/4/14/36/87/<=188 through 8192 — the
  SAME linear signature shifted ~2%; the desc1 (16384, 400) probe was a
  pure 1/row march projected to die at row 4071 (referee sec. E).
- 21 rules x D0 = 0..3 at R = 512: 84/84 die (stage C; die rows 107..174
  apart from asc's degenerate row-1 deaths; train-guard variants last).
- Spec sensitivity and D0-non-monotonicity (D1 sec. 4) — no monotone
  coupling argument can prove the barrier: desc1 closes (256,0) where
  plain dies; tr15 closes (512,4) where plain dies; adding vent depth
  can flip a closure to a death. Any proof must survive this chaos.

### 2.4 The proved skeleton a barrier proof can stand on (PROVED, new)

Let `d = d_0 = floor(gR) + D0` in the T4b window, `F_0 = R-2-d` (initial
front), `t_2 = 2d-R+3` (G2 ratio-2 crossover). Stage-B certificates:
21/21 grid points (dyadic 128..8192, D0 in {0, 1, ceil(R/32)}).

- **(E3-L1, complementarity).** `(t_2 - 1) + F_0 = d`. *Proof:*
  `(2d-R+2) + (R-2-d) = d`. The crossover cell and the front are exact
  binomial complements in d.
- **(E3-L2, front marginality).** Dominant front-cell load
  `A_{F_0} = C(R-2-F_0, d-F_0) = C(d, 2d-R+2) = C(d, F_0)` = exactly
  half the c-box width `2C(d, F_0)`. *Proof:* `R-2-F_0 = d` and L1.
  The front cell is born at box scale — the taper is marginal, which is
  why O(1) slack shifts decide the die/close dichotomy.
- **(E3-L3, boundary-layer rate).** Exactly
  `(A_{t-1}/A_t)/(C(d,t-1)/C(d,t)) = (R-1-t)/t`; at the front this is
  `(d+1)/F_0 -> g/(1-g) in (1.4874828, 1.4874887)`. The forced
  load-to-box ratio grows by ~x1.49 per cell of depth: control by ANY
  in-box choice is confined to an O(1)-deep taper at fixed D0
  (closed form for D1's measured "x1.5 per cell; boundary layer
  O(1)-thin").
- **(E3-kappa, taper constant; limit PROVED from L2/L3 + T4,
  finite table VERIFIED-exact).** Front-cell overflow / box width at
  D0 = 0 equals `w_{F_0}/width - (1-g)/g` with
  `w_{F_0}/width -> (3-4g)/(2(1-g))`, giving
  `kappa -> (3-4g)/(2(1-g)) - (1-g)/g in (0.0839816, 0.0839819)`;
  measured 0.084038 at R = 8192 (table in .out; at eps = 1/32 the ratio
  drops to ~0.0622). The initial front OVERFLOWS by a constant ~8.4% of
  its box: nucleation pressure at the front is Theta(box), not o(box).
- **(E3-L4, influence cone).** For any admissible continuation, the junk
  at cell c after k further rows is a function only of the row-i junk
  and the choices at cells in `[c - 2k, c]` (kernel support {0,1,2},
  choices act cellwise). *Proof:* induction on k. **Corollaries.**
  (i) Feed (cells {0,1}) cannot influence cell c before row `(c-1)/2`:
  the front region is feed-free for `~F_0/2 ~ 0.2R` rows — nucleation
  (measured rows 13..57) is decided by the transported initial data and
  the choices in the `[F_0 - 2k, F_0]` band alone.
  (ii) B1's steering-invariance and D1's window-depth invariance
  (descw32 = descw512 outcomes) are causal necessities, not accidents.
  (iii) For EVERY schedule, local or global, the first-super event at
  row k is decided inside a 2k-cell boundary layer: with first-super at
  o(R) rows, the battle is fought on an o(R)-cell strip. Any barrier
  proof is a statement about this strip-game.

### 2.5 The mechanism conjecture (the provable-looking target)

**Conjecture E3-BC(mech) (forced nucleation).** There exist `c > 0` and
`k(R) = O(log R)` such that for every dyadic R and `D0 <= cR`: every
admissible continuation whose choices, at every row and every cell of
the front band `[F(i) - 40, F(i)]`, deviate from the plain clamp by at
most a bounded multiple of the local box width (all RS rules qualify),
reaches a G3 super-pair state — hence dies, by the PROVED G3 theorem —
by row k(R). CONJECTURED.

Support and why it looks provable: the initial band data is exactly
alternating and convex (G2, PROVED/certified), so by G1 pure transport
creates no defects; defects are forced ONLY by the clamp taper (E3-L2:
the front is at box scale; E3-kappa: it overflows by a constant
fraction) and grow at ratio `(rho-1)^{1+delta}` in the expanding band
`t > t_2` (G1-d); bounded-multiple deviations can displace but not
cancel the defect mass budget (G1-c: transported L1 = 2x one-signed
D^2-mass); once a same-sign front pair crosses `Theta'`, G3 kills every
continuation. The 12/12 super-block separation, the G2-d defect-birth
dichotomy (front defects at rows 10/15 in dying runs), and the
first-super rows 13/17/27/31/39/57 (declining fraction of R) are its
finite instances. What is missing is the accounting lemma "bounded-
multiple band deviations cannot keep the band D^2-mass below box scale
for more than O(log R) rows at D0 <= cR" — the vent-depth flip at
(256, 0) shows the constant-fraction regime is delicate exactly at the
frontier, so the lemma's constants must (and, by E3-kappa/L3, can)
come out of the closed forms, not from sweeps.

**Honest labels.** A proved E3-BC(mech) would be a real theorem
("greedy-like schedules need linear slack") and would explain the whole
D0* landscape; it would NOT bound epoch feasibility (the T2 bijection
warning of 2.1) nor non-band schedules. E3-BC(param) would follow for
all RS rules with bounded vent/guard multiples.

## 3. Part (2): the existence gap at the frontier (512, D0 = 0)

### 3.1 What is now certified

- **(256, 0): integer point of the transportation polytope** (stage A;
  headline 2). Consequence for the program: the exact-floor profile
  itself stays integer-feasible through R = 256; the golden assembly
  needs `D0*(R) = O(1)`-ish behavior of FEASIBILITY (not of any rule),
  and the only frontier evidence against it is absence of a 512 witness.
- **(512, 0..3): all 21 RS rules die** (stage C table, exact); **the
  schedule beam** (12-action alphabet {plain, desc, desc1, desc2, desc3,
  descnv, tr15, tr43, tr54, tr74, train16, fill15}, protected pure
  lanes so the beam dominates every pure rule, two-key ranking by the
  G3 order parameter (same-sign-pair bits minus front-cap bits) and by
  (front, L1)) **all-dies at D0 = 0, 1, 2, 3** — last state dead at
  rows 109/111/113/202, ~226-288-state frontiers per row (best pure
  deaths among all 21 rules: 109/111/114/174; the beam alphabet lacks
  'fill', whose 114 at D0 = 2 outlives the beam's 113 by one row). At
  D0 = 3 the cross-row mixing outlives every pure rule by 28 rows
  (174 -> 202) and still dies — the G3 no-return signature: schedule
  gains are pre-nucleation only. Its (256, 0)
  control CLOSES (capture 153, one row before desc1's 154), so the
  machinery can find floor closures when they are reachable from its
  alphabet. VERIFIED-exact as searches; **no feasibility conclusion**
  (first beam attempt also pruned the known (256,0) winner until pure
  lanes were protected — a live demonstration of search-negative
  worthlessness, kept in the ledger per hazard discipline).
- **LP relaxation at (512, 0): status OPEN.** No real point known; the
  r128 LP note's mechanism ("per-cell real LP optimum IS the clamp")
  means greedy real-LP = rule A which dies, so LP feasibility at 512
  needs either a genuine LP solve (~1.6·10^5 variables, out of budget
  for exact rational simplex here) or a structured construction. The
  known necessary-condition families do not exclude it: (ARCH) at
  m = 512 passes at gamma* (the floor-audit ladder's largest refuted
  rate at m = 512 is 508/319 < gamma*), and the ballot cone /
  x-evaluation capacities are satisfiable at the floor profile (prior
  notes). No exclusion known.

### 3.2 Candidate structural exclusions (none proved; listed as targets)

- "No (512, 0) witness has the ballot normal form / enters the
  E-attractor by row r" for explicit r: decidable by the S_L endgame
  algebra on the reachable-set boundary; the r128-endgame-algebra hunts
  show the method (9 exhausted hunts at 128 before the LP closure).
- "Every (512, 0) witness leaves the NI class by row 15": plausible via
  the G3 noseed chain run backward from the plain death at 107 — would
  be the first PROVED "vent is necessary" statement.
- Full exclusion at (512, 0) (infeasibility) would refute exact-floor
  universality and force per-epoch slack >= 1 somewhere in the golden
  assembly — allowed by all current theory (slack O(1) is as good as 0
  for every C* consequence via LIFT).

### 3.3 What a (512, 0) witness would give

Exact-floor feasibility through 512, the pattern 0, 0, 0(?) continuing
at the frontier, and the sharpest form of the golden program (Section 4
route R0) alive at another doubling. Via THM-3329 it would also improve
the finite attainment for critical n <= 1023 from
`T(n) <= n + 1 + floor(gamma* n) + 4` (max slack over epochs 8..512 is
tr15's 4) to `T(n) <= n + 1 + floor(gamma* n) + O(1)` with per-epoch
slack 0 (epochs 8..256 already close at D0 = 0). Not claimed; the
search stands open.

## 4. Part (3): what would close C\* <= 1 + gamma\* (route ledger)

Three nested upper-bound programs (strongest first), then the floor side.

- **R0 (exact-floor universality).** Every dyadic epoch closes at
  D0 = O(1). STATUS: true through R = 256 (D0 = 0!); frontier (512, 0)
  open; no obstruction known. Gives `C* <= 1 + gamma*` with constant
  per-epoch slack (LIFT makes any O(1) enough). The transportation form
  (integer points now at 8..256) is the natural vehicle; what is
  missing is a rounding/construction theorem — and E3-BC says it will
  not come from RS rules.
- **R1 (GG + S', the golden gap).** o(R)-slack schedules keeping
  `sup_i m_i/Theta'_i < 1` through the feed phase, then S'. STATUS:
  open; per E3-BC the schedule must be non-local (cross-row); the
  concrete candidate is the LP-parity-rounding pipeline with the ballot
  invariant (closed R = 128 at the floor) upgraded from search to
  theorem: parity is free (proved), the ballot cone is the one global
  scalar cut (proved), so the target is a bounded-rounding-debt lemma
  against the (**) polytope's forced faces. Acceptance test: the exact
  G3 order parameter (margin profile in D3 sec. 6).
- **R2 (E-lin + S(eps)).** PROVED chain modulo Hypothesis S(eps):
  `C* <= 1 + gamma* + eps`; S verified exactly for dyadic
  128 <= R <= 16384 at eps = 1/32, 1/16; open at 32768 (probe died at
  1.9%, inconclusive — referee sec. B). Nearest unconditional
  milestone: prove S(1/32) → `C* < 1.6292382`. The proof shape is
  pinned (D2 sec. 5 envelope program; G1 calculus gives the feed-phase
  decay engine; C-A/C-F give the post-feed drain; only the ordered
  extinction bookkeeping and the alternation-integrity envelope are
  missing).
- **R3 (doubling / cancellation constructions).**
  `q^{R-1} = q^{R/2-1}·q^{R/2}` composition: blocked by profile-floor
  mismatch (doubling-structure note's obstruction record); no current
  path; would give R0 directly if solved.
- **R4 (amortized epoch decompositions).** Average per-epoch slack
  -> 0 with rare expensive epochs: no mechanism known that beats the
  per-epoch story; subsumed by R0/R1 in all current formulations.

**Floor-side honesty.** After the hostile floor audit, THM-3024's
general-class promotion is demoted: `C*_general >= log_5(5 phi^2)` is
OPEN (the transportation relaxation alone yields no general floor); the
golden floor stands for balanced-block schemes modulo one unwritten
Stirling-transfer lemma. So closing any of R0-R2 would pin
`C* = log_5(5 phi^2)` only against the block class; the general lower
bound needs the deadline-bounded-routing ingredient (or the HYP-9061
dual-side program). Upper-bound work and this note are unaffected.

## 5. Hazards honored

- Rule and search negatives (stages C, D; the 84/84 deaths; the beam
  all-dies) are NEVER used as feasibility evidence; the note's own
  record-correction (256, 0) is the standing counterexample.
- The class-collapse observation (2.2) is labeled OBSERVATION with its
  collision caveat; no theorem is claimed there.
- E3-BC(param) and E3-BC(mech) are labeled CONJECTURED; the proved
  items (L1-L4, kappa limit, stage A/B certificates) are separated.
- Quantifiers: all barrier statements are per-rule/per-class liminf
  statements over dyadic R; no non-dyadic claims; T4b window respected
  (all stage-B grid points have d_0 in the window).
- The beam dominates every pure rule in its alphabet BY CONSTRUCTION
  (protected lanes) — but only relative to that alphabet and width;
  no exhaustiveness claim.
- Witness verifications cited are the referee chain of record (D1
  quick-stage + Estimate-E fresh streaming verifier + this session's
  stage A); no new unverified closure is claimed anywhere.

## 6. Status ledger

- **PROVED (new):** E3-L1 complementarity; E3-L2 front marginality
  (front load = half box width, dominant term); E3-L3 depth-rate law
  with limit g/(1-g); E3-L4 influence cone + corollaries (causal
  steering-invariance, feed-free front phase ~0.2R rows, o(R)-strip
  reduction of the nucleation game); kappa limit formula from L2/L3+T4.
- **VERIFIED-exact (new):** stage A record re-verification + the
  (256, 0) transportation integer point (58837 cells, identity at
  x = 2,3,-1); stage B 21/21 lemma grids + taper table (kappa ->
  0.08398, eps = 1/32 value ~0.0622); stage C 84/84 rule deaths at
  (512, 0..3); stage D beam: (256, 0) control closure at capture 153,
  all-die at (512, 0..3).
- **OBSERVATION:** window-locality class collapse (2.2).
- **CONJECTURED / OPEN:** E3-BC(param) (barrier over RS); E3-BC(mech)
  (forced nucleation for bounded-deviation band schedules); exact-floor
  feasibility at (512, 0) — the frontier — and LP-relaxation
  feasibility at (512, 0); the structural-exclusion targets of 3.2.
