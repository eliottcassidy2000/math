# The Leader Ledger: the metric lives ON the walls — a conservation law for the nearest-runner process, and the Ostrowski climb as a stopping time

*boxeph-2026-07-12-S21. Owner prompt: make mathematical progress on LRC(14); be creative in
exploring past threads, especially unrelated ones. This session's seed was an "unrelated"
thread: mac-mini-S57's winding-tournament reflection, which ends on the honest verdict that
the tournament frame captures the ORDER of the runners and loses the METRIC ("the metric
lives at the optimum"). The leader ledger is the missing half: the metric lives ON THE
WALLS — every wall-crossing of the phase-clock movie carries an exact depth, and the depths
obey a conservation law. THM-722 (proved), LEM-025 (proved), HYP-6280 (data + refutation).*

## The object

Runners = distinct speeds, not all odd. At each time t there is a **leader** — the runner
nearest the observer's position 0 — with signed phase φ(t) ∈ (−1/2, 1/2). The function
|φ| = f is the distance-to-loneliness whose maximum is M(S). Three elementary facts
(THM-722, proof in canon):

1. φ rises at the leader's speed and only ever FALLS by jumping +x → −x at **sum-handoffs**
   ((v_i + v_j)t ∈ Z) of depth x = f(t). Difference events pass leadership continuously.
2. **The ledger balances:** ∫₀¹ v_λ dt = 2 Σ_handoffs x_h. Total climb = total fall.
3. The circle partitions into H⁺ = #handoffs **chains**, each with speed-mass
   x_in + x_out ≤ 2M, each crossing zero exactly once (a leader landing); if S has an even
   element, **H⁺ is even** (ι: t ↦ −t fixes exactly the chains through 0 and 1/2, pairing
   the rest).

Everything is exactly computable in rationals; verified on 50+ families with zero failures
(`lrc14_leader_ledger_boxeph_S21.py`).

## What the ledger showed (exact, this session)

- **The AP's witnesses are the units mod 14.** {1..13} has H⁺ = 58 handoffs, all depths
  1/q for q ∈ [14, 25]; exactly six reach 1/14 — at t = a/14, gcd(a,14)=1: the φ(14) = 6
  unit witnesses (matches codex-S175's "first_unit_witness=1/14(6)" from the other side).
- **The deep well's ledger IS the Ostrowski climb.** {1..12,182}: the first 14 handoffs are
  t = k/183, pair (1,182), depth k/183, k = 1..14 — a linear staircase on the (1,182)
  sum-ruler (spacing 1/183 = 1/Φ₆), cut at k* = 14 by runner 12's lander: 183 − 12k < k ⟺
  k > 183/13, i.e. **k* = ⌊183/13⌋ = 14 is a stopping time**. mac-mini cont.56's "omit the
  distance-1 lander" and the +1 in Φ₆ = 13·14+1 (the ruler offset) become one floor
  computation.
- **One one-line lemma is tight at ALL FOUR extremals.** LEM-025: with v = min, f = max,
  B = second max, q = v + f: M ≥ v·⌊q/(B+v)⌋/q. Tight (bound = exact M) at the AP (1/14),
  every Ostrowski rung {1..12,13k} (k/(13k+1)), the deep well (14/183), and the compressed
  extremal 2·{1..12}∪{13} (1/13). The whole bottom of the covering M-spectrum sits ON the
  (min,max)-ruler climb. (GW is the exception that probes the rule: its witness 3/14 rides
  the (1,13) ruler, and the lemma is not tight there.)
- **Two-line, citation-free closure of the covering {1..12}∪{f} stratum:** covering forces
  182 | f, so q = f+1 ≡ 1 (mod 13) and LEM-025 gives M ≥ (f/13)/(f+1) ≥ 14/183 for all
  f ≥ 182. (Statement previously known via mac-mini's exact family formula + LRC(13); the
  ledger proof needs neither.)
- **Equioscillation is graded, and the covering rung does NOT reach it.** Witness-chain
  ratio (x_in+x_out)/2M along the ladder k = 1,2,3,7,14: 0.78, 0.75, 0.83, 0.93, 0.964 —
  climbing toward but never reaching 1. Contrast: the AP's zero-chain and the compressed
  extremal 2·{1..12}∪{13} equioscillate EXACTLY (ratio 1). **My conjecture "covering ⟹
  max-chain-mass ≥ 28/183" is REFUTED by the deep well itself** (2639/17751 < 2716/17751):
  the covering-min extremal is NOT chain-equioscillating. Honest surprise; logged in
  HYP-6280.
- **Ledger efficiency η = ∫v_λ/(2H⁺M) decays along the ladder** (0.785 → 0.534): the deep
  well pays 220 handoffs where the AP pays 58 — covering families waste depth on shallow
  handoffs. The average-depth bound M ≥ ∫v_λ/(2H⁺) is a factor ≈ 1.9 below M at the deep
  well, so the conservation law ALONE cannot reach the crux (recorded limit, not a route).

## The threads it joins (the creative-mining ledger)

- **Winding tournament (mac-mini-S57, THM-373):** difference-handoffs are wall-crossings of
  the runner phase clock; the leader FSM is a metric OBSERVABLE of the tournament movie
  that the iso-class census forgets. Order (tournament) + depth ledger (metric) = the full
  object.
- **ι-parity (klein-S270, the Rédei aesthetic):** "odd-modulus loneliness comes in ι-pairs"
  gets a geometric home — witnesses are deepest handoffs, handoffs pair under ι except the
  two fixed chains, H⁺ is even. Same proof shape as the project's involution-pairing
  parity arguments (fixed points carry the parity).
- **Chebyshev equioscillation (mac-mini-S40 thread):** the ledger quantifies it (η, chain
  ratios) and shows the covering-min extremal is *not* equioscillating in the chain sense —
  the two-point binding pair (x_in = x_out = M) is the AP/compressed phenomenon, not the
  deep well's.
- **Rotation orbits (mac-mini cont.56) / Ostrowski (klein-S269):** ruler depths ARE
  rotation-orbit closest approaches; the ladder rungs are stopping times of the climb; the
  ledger is the dynamical object whose maxima the witness-search formula enumerates.

## Where this could bite the crux (open, for the fleet)

The crux in ledger coordinates: **covering ⟹ the deepest handoff is ≥ 14/183.** The ledger
adds structure to "deepest": depth accumulates only along sum-rulers, linearly at rate
v_pair/q per step, and is cut only by landers. So the crux = "a covering family cannot cut
every ruler's climb below 14/183." The AP cuts all climbs at 1/14 — with landers at every
k/14 — and is non-covering precisely because that lattice is full. The quantitative handle:
cutting a (v,f)-ruler climb at step k needs a runner u with ||u·k/q|| < vk/q, i.e. a lander
within vk of the ruler point — a dense-lander requirement that fights the divisor-
completeness pattern. This is the lander-exclusion count (klein-S270's honest target (a))
with an explicit process attached. NOT claimed proved; the deep well shows the minimal
sustained climb (14 steps at the slowest rate 1/183).

## Engineering note (mandate)

The ledger is an exact interleaving auditor for rotating/phased systems: given integer
rates, it outputs the full schedule of "who is nearest the shared origin", handoff times on
the (rate_i + rate_j) rulers, depths (safety margins), and the conservation check — O(Σv)
events, exact rationals, pure Python. Candidate uses: TDMA/token-passing collision margins,
cam/gear mesh phasing, polling schedules. Seed: `lrc14_leader_ledger_boxeph_S21.py`
(leader_ledger() is dependency-free).

*Files: THM-722 (canon, proved), LEM-025 (inside THM-722), HYP-6280,
04-computation/lrc14_leader_ledger_boxeph_S21.py,
05-knowledge/results/lrc14_leader_ledger_boxeph_S21.out. Credits: mac-mini-S57/cont.55/56,
klein-S267/S269/S270, opus-S249–252, kps cont.51–55, THM-668-mac-mini, THM-373, codex-S175,
Kravitz (pair-sum reduction precedent, backlog klein-S253).*
