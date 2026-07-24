---
source: kind-pasteur-2026-07-23-S132 (Opus 4.8)
status: MECHANISM IDENTIFIED (strong numerical+structural evidence) + a creative program to push it to 1/14.
  The snippet is a SOUND Riesz-weighted variational lower bound on the lonely-runner GAP (not a log-energy
  measure certificate — my S130 worry was about the wrong functional). This reframes LRC(14) as a
  "certifiable concentration" problem and yields seven concrete attack angles. Computationally grounded.
tags: [lrc, lonely-runner, snippet, variational-bound, gap, free-energy, certifiable-concentration, program, ideas]
related: [kps-S129, kps-S130, kps-S131, THM-518, THM-515, klein-S404/405, macmini-S169, THM-252]
---

# The snippet is a variational gap bound — and LRC(14) as certifiable concentration

**kps-2026-07-23-S132.** Owner: leverage the snippet toward FINISHING LRC(14). I found what the snippet
mechanically IS, which turns the abstract "log-energy" talk into a concrete, sound, and extensible method.

## 1. Identification (strong evidence): a SOUND variational gap bound
For any probability measure μ, `gap(S) = max_τ min_v‖vτ‖ ≥ ∫ min_v‖vτ‖ dμ`. Take μ = R/∫R with the
**Riesz product** `R(τ) = ∏_{v∈S}(1 + a·cos 2πv(τ−τ*))`, `a<1` (so R>0), centred at the lonely time τ*.
Computed for `S={1,…,13}`, τ*=1/14:

| a (amplitude) | ∫g·R/∫R | note |
|---|---|---|
| 0.55 | 0.04530 | |
| **0.59–0.60** | **0.04570–0.04578** | **= snippet X = 0.045725 to ~0.1%** |
| 1.00 | 0.04721 | |

`a=0.6 ⟹ ρ=(1−√(1−a²))/a = 1/3 ⟹ 2·arctanh(⅓)=log2` — mac-mini's magic amplitude. **The snippet's certified
value is (to 0.1%) this variational bound.** Consequences that settle prior open questions:
- **SOUND.** `∫g dμ ≤ max g = gap` for any probability μ. This is an *average of g*, not the entropy
  `∫M·log R` I flagged unsound in S130 — that worry was about the wrong functional. No sign problem.
- **LOSSY ⟹ WIDER-GAP, not the conjecture.** For `{1,…,13}` it gives 0.0457 vs the true gap 0.0714 (66%).
  A band-limited μ is not a point mass, so `∫g dμ < gap` strictly. The snippet therefore proves
  `gap > 1/25`, a *uniform improvement on the union bound `1/26`* — it structurally **cannot reach 1/14**.
  (Resolves the S131 fork: reading (a). If the LLM claimed the full conjecture, this step is insufficient.)
- **Why logs.** `∫g·R/∫R` in closed form runs through the tent's odd-harmonic Fourier series against R's
  relation structure; the snippet packages it as `2·arctanh` (=rapidity) values of two rationals. (Exact
  closed form still to pin — the tent gives 1/k² rationals, so the two logs likely come from a
  free-energy/`log∫e^{βg}R` refinement of the plain average, which is ≥ the average and still ≤ gap.)

## 2. The master reframe: LRC(14) = CERTIFIABLE CONCENTRATION
`sup_μ ∫g dμ = gap` **exactly** (μ → point mass at τ*). So LRC(14) is *entirely* about the CERTIFIABLE
class: **over nonnegative trig-polynomial measures whose `∫g·R/∫R` is evaluable by the artanh sandwich, is
the sup ≥ 1/14?** Two obstructions, both computationally exhibited this session:
- **(O1) Lossiness.** Certifiable (band-limited Riesz) μ can't fully concentrate; degree-Σv=91 gives ~0.047
  for `{1,…,13}`, short of 0.0714.
- **(O2) Centering.** The bound is fragile to τ*: centring μ at 1/14 for a drop-core `{1..13}\{6}∪14` gives
  **0.034 < union bound** (its true peak is elsewhere); centring at its own argmax recovers 0.042. **For a
  general config τ* is unknown a priori** — this is the crux the pure-AP case hides.

## 3. Seven creative angles to close 0.047 → 0.0714 (the "factor-2" difficulty, Tao)
**A. Certifiable-concentration calculus.** Characterise the loss `gap − sup_{deg≤D}∫g dμ` vs degree D and
find the D that reaches 1/14. Non-uniform amplitudes + 2nd harmonics + products-of-Fejér push the class
toward point masses. *Testable now.*
**B. τ*-location via the three-distance theorem.** The optimal lonely time sits at a rational `a/q` fixed by
the config's continued-fraction / Stern–Brocot arithmetic (three-distance structure of `{vτ}`). Reduce
"find τ*" to a CF search, then centre μ there — makes the bound config-adaptive and kills (O2).
**C. Primal–dual.** `sup_μ ∫g dμ = gap` is an LP; its dual is a *covering/packing potential* certificate.
Strong duality ⟹ dual optimum = gap. Construct the DUAL witness (a decomposition of g forcing gap ≥ 1/14)
instead of μ — often more directly artanh-certifiable.
**D. Curvature matching = why Σv².** Loss `= ∫(gap−g)dμ`; near τ*, g is a KINK with slopes = binding speeds,
aggregate curvature `½Σv²` (klein-S405). Optimal μ-bandwidth `~1/√(Σv²)` matches the peak and minimises loss
— **this is what the snippet's `2457=3Σv²` weight encodes.** Make bandwidth config-adaptive via Σv².
**E. Ground-state energy / reflection positivity.** `gap = −min_τ H`, `H(τ)=−min_v‖vτ‖` a cosine-coupled
lattice Hamiltonian. Deploy reflection positivity + chessboard estimates (standard for such H, unused in the
repo's LRC) for a ground-state lower bound.
**F. Finite reduction to gap-near-extremizers.** The gap is minimised at `{1,…,13}` and dilates (=1/14).
Near-extremizers form a structured finite family (gap analog of THM-518's resonant strangers). Then
`gap ≥ 1/14 − ε ∀S` reduces to a snippet-style certificate per class — a FINITE check. The snippet is one.
**G. Odd-harmonic bridge → Lean-ready.** Tent (odd harmonics, 1/k²) and artanh (odd harmonics, 1/k) share
support; a generating-function identity certifies `∫g·R/∫R` float-free. Every numerical loneliness bound
becomes machine-checkable (opus-S2 engine). The snippet is the prototype.

## 4. Why this beats the stalled routes (the payoff)
THM-518: the Riesz-product **measure** certificate (`∫M·R/∫R<1`) STALLS on AP-cores (low additive dimension).
The variational **average** `∫g·R/∫R` is a *different functional*: it is powered by the L²-energy `Σv²`
(large on APs, angle D), so it is **strong exactly where the measure-certificate is weak**. That
energy-vs-additive-dimension decoupling (klein-S405) is the structural reason this line can go past the
stall. It is sound, beats the union bound, and is second-moment-driven — a genuinely new tool for `inf L>0`.

## 5. Honest assessment + next experiment
The snippet-method is sound but caps ~0.047 for `{1,…,13}` — a wider-gap result, NOT the conjecture. Closing
to 1/14 needs angles A+B+D (better certifiable measures, config-adaptive centre & bandwidth) and/or C/E, then
F to finitise. Files: `/tmp/{variational,match,floor,freeenergy}.py`.

## 6. Experimental results (this session) — the sharp picture
Ran the experiment. Two clean facts for `{1,…,13}`, τ*=1/14:
- **The Riesz product PLATEAUS at 0.0472** (per-frequency amplitude optimisation drives every `a_v→0.99`
  and stops). It is a *weak concentrator* — its low-frequency factors (`v=1,2`) are broad and cap the
  concentration. **So the snippet is a weak, low-degree instance giving the wider gap ~0.047, and by itself
  CANNOT approach 1/14.**
- **A Fejér concentrator CLIMBS to the gap:** `∫g dμ` = 0.051 (D=91), 0.065 (D=500), 0.0689 (D=1500),
  0.0705 (D=5000 = 98.7% of 1/14). So `sup_μ ∫g dμ = gap` is realised — the variational *route* is not capped;
  **the Riesz product is simply the wrong measure.** Use Fejér/Jackson kernels, not Riesz products.

**The sharp limit (crucial).** `∫g dμ < gap` *strictly* for every non-point-mass μ. `gap({1,…,13}) = 1/14`
*exactly*. Therefore the variational method proves `gap > 1/14 − ε` (any ε, large enough D) but **never the
sharp `gap ≥ 1/14`** for the extremiser or any near-extremiser whose gap `= 1/14 + o(1)`. This is exactly the
"exact-value" wall (THM-518). **Consequence — the division of labour to finish LRC(14):**
- **Bulk configs** (gap bounded away from 1/14): a high-degree Fejér variational bound, config-adaptively
  centred (angle B) and certified float-free (angle G), gives `gap > 1/14`. *Feasible via this route.*
- **Near-extremisers** (gap `→ 1/14`): the lossy variational bound fails here; they need **exact/rigidity**
  input — prove the only configs with `gap ≤ 1/14 + δ` are dilates of `{1,…,n}` (a stability theorem, angle F),
  which have `gap = 1/14` known exactly. This is the true crux, and where THM-518's resonance analysis lives.

**Net:** the snippet is a sound but weak wider-gap tool. Its real gift is the *method* — a sound,
second-moment-powered, artanh-certifiable variational gap bound that (with a proper Fejér concentrator +
config-adaptive centring) reaches `1/14 − ε` for the bulk, reducing LRC(14) to a **rigidity/stability
statement about near-extremisers**. That reduction is the concrete path to the finish. Files added:
`/tmp/{optimize,opt2,fejer}.py`.
