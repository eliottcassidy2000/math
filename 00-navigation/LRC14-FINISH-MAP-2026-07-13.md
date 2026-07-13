# LRC(14) FINISH MAP — 2026-07-13 (klein-S284)

**The definitive current state.** Supersedes `LRC14-FINISH-MAP-2026-07-11.md` (which predates both the
completed covering-min rigidity, THM-724/726, and the fully-reduced density route, klein S273–S283).
Synthesis of the converged fleet state: covering side (mac-mini THM-724/726/730, kps cont.55–68, opus
S253–266), density side (klein S272–283), Lean spine (klein S242–252, kps).

**One-line status.** LRC(14) is reduced to a single equidistribution/cancellation inequality, which
appears in two independent, now-unified forms — one per route. A vast rigorous skeleton surrounds it on
both sides; the remaining inequality is genuine harmonic analysis, and every *elementary* tool has been
tried and is exhausted on both routes (klein S283, opus S266).

---

## 1. The top-level reduction — PROVED / SETTLED

LRC(14): every 13 nonzero integer speeds have a `1/14`-lonely time.

1. **Non-covering ⟹ done.** A speed set omitting a multiple of some `q∈{2,…,14}` has the explicit
   witness `t=1/q` (`M≥1/q≥1/14`). Uses LRC(≤13) (cited per project policy, SETTLED). **THM-366/523.**
   So LRC(14) reduces to **covering** sets (a multiple of every `q≤14`) ⟹ **divisor-complete (DC)**.
2. **Realization — supply constructive, in Lean.** A positive good-set floor yields an explicit
   bounded-denominator lonely time; the `≤7`-arcs pigeonhole (SmallClusterFull) is PROVED in Lean
   (LRCSevenGapRigidity), and the certificate supply is kernel-pure Lean (klein THM-693–698).

So the entire remaining content is one statement about the **covering / divisor-complete** class, with
**two independent proof routes — proving EITHER finishes LRC(14).**

---

## 2. Route B (covering / divisibility) — STRUCTURE PROVED, one analytic inequality left

**The covering-min rigidity is now a THEOREM.** The minimum of `M` over primitive covering 13-sets is
`14/183 = n/Φ₆(n)` (`Φ₆(14)=183`), attained **uniquely** at the deep well `{1..12,182}` (`182=14·13`):

- **THM-724 (PROVED):** every primitive covering *single-killer* 13-set has `M≥14/183`; deep well is the
  unique minimizer (balance lemma + shallow-witness counting + prime-13 tightness). Kernel-pure Lean
  (kps cont.60/61/66, LRCSingleKillerLadder).
- **THM-726 (PROVED):** every primitive covering 13-set with `≥2` far outliers has `M≥1/13>14/183`
  (far-element monotonicity + 64317-config finite check + THM-724). ⟹ deep well is the **unique global
  covering-min**.
- **Structural skeleton complete (kps cont.68):** single-killer (proved+Lean); multi-killer =
  separate-13/14 (`core≤11`, enumerated+reduced, cont.58/59); the 2-runner extremal unifies frames but
  introduces no new cases.

So `LRC(14) ⟺ inf_{covering} M ≥ 14/183` (the first covering Ostrowski rung; AP is rung `k=1=1/14`,
non-covering). The **combinatorial extremal step is also PROVED:**

- **THM-730 (PROVED, elementary):** the AP uniquely maximizes additive triples — `T(A)≤C(k,2)`, equality
  iff `A` is a dilated AP (the E₃/Schur inverse, opus-S182's target).

**The single remaining piece [B]:** the `|core|=1` smooth-body discrepancy — equivalently the
*resummation* linking the Schur/E₃ deficit to loneliness (mac-mini-S78: `L=(6/7)^{13}(1+\text{corrsum})`,
`L>0 ⟺ \text{corrsum}>−1`; "Schur deficit ⟹ `L>0`" **is** the covering case, not yet proved). opus
(S262/266): `ε_v` is an **alternating multi-linear (Gowers-type / E₃)** cancellation of a distinguished
core arc against the product good-set `∏_w(1−1_{D_w})`; the completion identity handles the *bilinear*
(pairwise) part cleanly, but the signal is `≥2`-way. **Elementary tools exhausted (opus-S266).**

---

## 3. Route A (density / moment) — FULLY REDUCED, one analytic inequality left

The uniform good-set floor's compact base rows are `k=8` (`Φ≤cap₉=1979/4004`) and `k=9` (`J≥432/91`);
`k≥10` inherit via the eigen-transfer **THM-710 (PROVED)**. The whole route is now a rigorous chain
(klein S272–283) bottoming on one cancellation:

- **k=8 row CLOSED** modulo the tail constant (klein-S272): `Φ` majorant PROVED (deg-3 LP), compact
  `d≤25` EXHAUSTIVE (THM-719), tail = the explicit `Φ_∞` transfer.
- **The tail constant reduces (rigorously) to a 2nd moment.** The far-element two-scale error
  `Error=Φ(E'∪w)−Φ_∞(E')` satisfies `|Error| = |S|/w`, and (THM-727) `S=Σ_s Σ_ℓ (−1/2πiℓ)U_s(ℓw)ĝ_s(ℓ)`
  with `U_s` the `R_s`-endpoint exponential sum; (THM-728) `U_s^{e'}(N)=e'·χ̂(N mod e')` a **1-D DFT of a
  swing-derivative**; (THM-729) the 2nd moment `Q_s=Σ_ℓ|U_s(ℓw)|²/ℓ²` **equals** `(2πw)²·[Riemann-
  discrepancy of the autocorrelation of `1_{R_s}` at the `w`-grid]` — a **1-D discrepancy**.
- **Cauchy–Schwarz:** `|S|=O(√Q_s)`, so **`Q_s=O(r)` ⟹ Error→0 ⟹ the row closes** (`r=#R_s`-arcs`=O(diam)`;
  peel `w=d≥diam` ⟹ `Error=O(1/√d)→0`).
- **RIGOROUS around it:** crude `Q_s≤4π²r²/3` (BV Fourier); the closed-form diagonal `diag=Σ_i
  4π²{w w_i}(1−{w w_i})=O(r)`; and the *downgrade* — **any** power-saving `Q_s=O(r^{2−ε})` suffices (S281).
- **VERIFIED-not-proved:** the sharp `Q_s=O(r)` (klein-S280, `Q_s/r∈[1.0,1.7]` to `diam 199`).

**The single remaining piece [A]:** `Q_s=o(r²)` — the oscillatory cancellation of the width-weighted
Weyl sum `Σ_i w_i e(−ℓw c_i)` (arc midpoints `c_i` under `×w`), equivalently `offdiag` = the
convex(`B₂`)-minus-lattice balance being `o(r²)`. **Elementary tools exhausted (klein S281–283):** large
sieve gives `O(r³)`; Montgomery–Vaughan / width-weighted give `O(r²)`; `offdiag≤0` is **refuted** (one
counterexample); `B₂`-convexity anti-correlation is overwhelmed by lattice-straddles.

---

## 4. THE UNIFICATION — both routes bottom on the SAME estimate

The two remaining pieces [A] and [B] are the **same kind of object** (klein-S279):

| | Route B (covering) | Route A (density) |
|---|---|---|
| distinguished element | a **core arc** `D_v` | a **swing offset** `e'` |
| the product set | good set `∏_w(1−1_{D_w})` | cover set `∏_{e''}(1−1_{·∈∪A})` |
| the object | `ε_v=\mathrm{Cov}(D_v,G')` | the endpoint sum `U_s^{e'}` / `Q_s` |
| bilinear part | completion identity `Cov(D_v,D_w)≤1/(3vw)` — clean | derivative gain `sin(πn/7e')` kills `n=0` — clean |
| the crux | `≥2`-way (Gowers / E₃) cancellation | oscillatory (Weyl / convex-minus-lattice) cancellation |

Both are **a distinguished element correlated against a product of arc-indicators**, bilinear-clean,
with the signal in the higher-order equidistribution. The two differ in *what they need*:

- **Route B is TIGHT** (`14/183`, no slack): needs the **sharp** multi-linear (Gowers/E₃) cancellation.
  Naive Erdős–Turán is `~700×` too weak; averaging is provably insufficient (mac-mini-S76). Harder.
- **Route A has SLACK** (box extension): needs only **any** power-saving `Q_s=o(r²)` — a *soft*,
  1-linear equidistribution bound. Lower-order and more tractable, but the same category.

So LRC(14) is one equidistribution cancellation, and Route A is the softer face of it.

**⚠ S285 sharpening (klein) — it is literally ONE lattice, three cosets.** The unification is exact, not
merely "same kind": both faces are sums of `Ĝ(a)` over the speed **relation lattice** `L={a : a·E'=0}`.
Density's `f̂(ℓw)=Σ_{a·E'=ℓw}Ĝ(a)` is the **`ℓw`-coset**; covering's `ε_v=Σ_{h≠0}b_hĝ(−hv)` is the
**zero-coset** `L`; the old support-6 kernel (THM-538, MISTAKE-078) is the **same** relation-lattice sum
(the conditionally-convergent one — no Minkowski count closes it; convergent form = the finite `x`-cell
integral HYP-2645 = the `Q_s`/`ε_v` we keep landing on). Density's slack is a **coset advantage**: the
`ℓw`-coset is bounded away from `0`, so its terms are uniformly high-order — which the covering zero-coset
lacks. (klein-S285, HYP-6455.)

---

## 5. What is PROVED vs OPEN — the ledger

**PROVED (rigorous, much of it Lean):** the top reduction (THM-366/523); the covering-min rigidity and
deep-well uniqueness (THM-724/726); the Schur/E₃ inverse (THM-730); the density eigen-transfer (THM-710);
the density Fourier/DFT/autocorrelation reduction (THM-727/728/729) with crude `Q_s≤4π²r²/3`,
closed-form `diag=O(r)`, and any-ε-suffices; the k=8 majorant + `d≤25` exhaustive box (THM-714/719);
the realization/certificate supply and SmallClusterFull (Lean).

**VERIFIED-not-proved (extensive computation, 0 exceptions):** the sharp density `Q_s=O(r)` (S280); the
covering good-set positivity across the DC class (opus); the deep well as the numerical global min.

**OPEN — the single inequality, two forms:**
- **[B]** `inf_{covering} M ≥ 14/183` via the smooth-body / resummation cancellation (`Schur deficit ⟹
  L>0`) — sharp multi-linear (Gowers/E₃).
- **[A]** `Q_s=o(r²)` (any power-saving) — soft oscillatory (Weyl) cancellation of the arc midpoints.

**Proving EITHER finishes LRC(14).** Elementary methods are exhausted on both; the remaining work is
genuine harmonic analysis / equidistribution — a real analytic task (or an external result), not a
further elementary reduction. Route A is the softer target.

---

## 6. Concrete next moves

1. **The soft Weyl bound [A]** — `Q_s=o(r²)` for the arc midpoints under `×w`. Softest target; any
   power-saving closes it. Needs the convex-minus-lattice balance or an equidistribution/Weyl estimate
   beyond the large-sieve family.
2. **The sharp multi-linear bound [B]** — the E₃/Schur resummation `Schur deficit ⟹ L>0`
   (mac-mini/opus). Harder (tight), but the combinatorial half (THM-730) is done.
3. **Import external equidistribution** — both faces are standard-flavoured; a targeted lemma
   (van der Corput / Weyl for structured dilated point sets, or a Gowers-inverse input) may discharge
   one side.
4. **Lean consolidation** — the proved skeleton (both routes) is largely kernel-pure; the finish is
   `[proved skeleton] + [one cited analytic inequality]`.

*Supersedes LRC14-FINISH-MAP-2026-07-11.md. Sources: THM-366/523/710/714/719/724/726/727/728/729/730,
HYP-2566/6350/6415/6440/6445, klein S272–284, mac-mini S70–78, kps cont.55–68, opus S253–266.*
