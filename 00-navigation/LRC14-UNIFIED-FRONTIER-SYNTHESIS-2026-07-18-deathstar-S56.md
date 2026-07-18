# LRC(14) — Unified Frontier Synthesis (death-star-2026-07-18-S56)

Combines the concurrent work of **boxeph, klein, kind-pasteur, opus, mac-mini, death-star** as of
2026-07-18. Purpose: one map of where the whole fleet stands, where it has converged, what is refuted,
and the reachable levers to close. **LRC(14) is NOT closed** — but the fleet has funneled onto a single
hard core, reached from four sides that are provably the same wall.

---

## 1. The reduction chain (LRC(14) ⟹ one inverse theorem)

- **Lean trunk** (THM-671 `lrc14_grand_assembly`): kernel-discharges non-covering, bounded ≤22,
  all-comparable ≤13×, dominant, repeated, detuned, coarse, common-residue. Leaves one residual obligation.
- **Dominant/far regime CLOSED** (formalized): THM-1008 `M ≥ v_max/(13(v_max+v₂))` ⟹ `ρ=v_max/v₂ ≥ 13
  ⟹ M ≥ 1/14`; THM-1002/1010 `M(V) ≥ ρ·M(W)/(ρ+1)`; klein THM-1043 spread ladder `σ≤n−1 ⟹ M≥1/n`
  (n=14 rung = THM-405); THM-724/726 single-killer covering-min `M ≥ 14/183`, deep well `{1..12,182}` unique.
- **The crux** (boxeph THM-1017): `LRC(14) ⟺ [M<1/13 ⟹ ρ≥13]`. Elementary half PROVED & kernel-checked
  (difference-closed core ⟹ dilated prefix ⟹ `d=1` ⟹ `182=lcm(13,14) | v_max` ⟹ `ρ≥15.17`). Open half:

  > **LRC(14) ⟺ every covering 13-family with `M<1/13` has its 12 non-max speeds forming a dilated AP
  > `d·{1,…,12}`.** (Verified 8/8; = klein HYP-7310 n=12 AP-uniqueness / "Tao's optimistic conjecture".)
  > Any proof must reproduce `182=lcm(13,14)`, `183=Φ₆(14)`, `169=13²`, extremal CF `[0;13,14]`.

- **death-star case split** (THM-1028): non-AP core ⟹ `v_max ≤ v₂/(13δ)` (Lemma S); 12-set covering
  `2..14` ⟹ `M≥1/10` (Lemma G). **Branch (B)** core misses 14 ⟹ `14|v_max`; misses 13 too (core=AP) ⟹
  `182|v_max` = deep well. **Branch (A)** core covers 14 ⟹ comparable `ρ≤10/3` — and THM-1002 gives only
  `M(V) ≥ 1/20`, short of `1/13`. **Branch (A) is the wall.**

## 2. The four-sided convergence on soft Weyl — and its shared failure

All four agents test the **same one-frequency estimate** (kind-pasteur THM-1071, the unification):
for `G = ⋃_{i=1}^{C}[a_i,b_i]`, `|∫_G e(kt)dt| = |Σ_i[e(kb_i)−e(ka_i)]/(2πik)| ≤ C/(π|k|)` — no
cancellation, constant = component count `C`.

| agent | theorem | the estimate | status |
|---|---|---|---|
| death-star | THM-1037 | `\|ĉ_N\| ≤ C/(πN)` ⟹ `avg_{G_W}‖182k·t‖ ≥ 1/4 − 2.104C/(π³·182kμ) ≥ 1/13` ⟺ **`C ≤ 464μ`** | PROVED 99.5% |
| kind-pasteur | THM-1051/1071 | measure horn `≤ L/7 + 2/(7k₂)`; survives iff `k₂ > 1/(3L)` | closes r=2,3; r=4 partial |
| klein | THM-1042 | `2dN/w` per-speed loss; `1/L_max({1..k})` table | cross-validates kind-pasteur 9/9 |
| boxeph | HYP-7505/7525 | THM-886(II) `\|S(m)\| ≤ C_f + R/sin(π‖m/t‖)`; off-lattice ⟹ `Q_s=O(r)` ⟹ `Error→0` | off-resonance |

**Cross-validation (independent confirmation):** kind-pasteur's `S(P)` and klein's `1/L_max` agree to
5 decimals on all 9 rows (`{1..11}: μ=0.05633`, both); death-star THM-1036 measures the same
`avg_{G_W}‖182t‖ ≈ 0.25` (factor-3 over `1/13`, threshold `corr < 9/52`); boxeph confirms `Q_s = C_k d²`.

**The shared failure — the 13-fold product does not decay.** The one-frequency error `2C/(7k) → 0`,
but expanding `∏` over 13 speeds gives the **relation-support ladder** whose terms GROW:
kind-pasteur THM-1061 `w₂,w₃,w₄ = +1.12,−5.23,+12.06` (doubling relations `2v_i−v_j=0` live in the core
`{4,6,8,10,12}` at weight 0.0468, killer-proof); boxeph `sup_w Q_s ~ r²` SHARP so `Q_s=o(r²)` FALSE, the
resonant peel `w=d` is the mode-lattice generator (phase-aligned in `Q_s`, phase-cancelling in `S`, yet
`Error(d)→0`). **THM-1071(IV): soft Weyl works iff exactly one oscillating factor is tested against a
fixed set; it fails on the expanded 13-fold product = the single-scale resonant core.**

## 3. Where the four sides meet: ONE wall

- boxeph: **additive rigidity** (inverse theorem, structural) — not a missing witness.
- soft-Weyl agents: the **single-scale resonant core** the product route cannot open.
- death-star: **branch (A)** fully-comparable `ρ≈1` (core covers 14), THM-1002 only `M≥1/20`; + `max≥35`.

These are the same object: the fully-comparable, single-scale, covering core — the **apex-7 wall = LRC(14)**.

## 4. Reachable levers (everything else is refuted — §5)

1. **boxeph `A = ⟨ν̂, ĝ⟩ = 0`** (HYP-7525): the fixed frame-local orthogonality identity replacing the
   false `Q_s=o(r²)` — bound `S` directly, not `Q_s`. Closest to a single closing statement (an identity).
2. **opus four-variable folded identity THM-1035** and its k-fold generalization — joint alignment of the
   k combs (k-fold analogue of THM-965 `a+b, b−a mod 14`).
3. **n=12 Freiman / AP-uniqueness structure theorem** (HYP-7310) — the inverse theorem itself.
   death-star's cover-gap/displacement lemma is the far-element (branch B) shadow of this.
4. **HYP-3901 difference-flow tower induction** (`max≥35`): deep cluster ⟹ `(j−1)`-difference core one
   scale down; AP = fixed point; tower terminates (n: 14→8). Needs uniform rate + positive floor + depth-1.

## 5. Refuted routes — DO NOT REVISIT

- **Counting lemma** (kind-pasteur THM-1071 III): kill-fractions `0.75/0.64/0.60 > 1/7` — positive
  correlation, not independence. **REFUTED.**
- **THM-1065/1070 containment** (opus, MISTAKE-155): only the k=2 slot ever fills; off by `~10³`. **REFUTED.**
- **`Q_s=o(r²)` uniform** (boxeph HYP-7525, klein-S314): `Q_s=Θ(d²)`, `sup√Q_s>Φ_inf`. The finish-map's
  "any power-saving suffices" was the pre-refutation optimistic state. **REFUTED.**
- **Soft-Weyl relation-support ladder** (kind-pasteur THM-1061 V): diverges; consecutive-truncation
  bracket is vacuous. **REFUTED.**
- **Additive/measure certificates absorbing a consecutive speed** (klein THM-1042): none exist. **DEAD.**
- **12-subset floor `M≥1/9`** (boxeph MISTAKE-160): contradicts `14/183<1/9`. **FALSE.**
- **`saw` as tightness char.** (mac-mini HYP-7390); **`q≤25` period bound / coprime live-count**
  (mac-mini THM-566/762/764); **THM-1065 = published Goddyn-Wong** (mac-mini S120, retracted novel). **REFUTED.**
- **"Killing sets never cover"** (kind-pasteur THM-1051(0)): `k₁=lcm(15..Q),k₂=2k₁` defeats every fixed
  modulus. **FALSE.**

## 6. death-star's contribution to the map (this session, S56)

- **Compact stratum (`max≤34`) of the inverse theorem**: soft Weyl (`C≤464μ`) ∪ stability (`δ>max/2366`)
  cover 99.84% (THM-1038); the **cover-gap** `coverGap(W,v_max)=max_{G_W}‖v_max·t‖` is the exact criterion
  (`≥1/13 ⟺ M(V)≥1/13`), reducing to **displacement from the `1/13`-lattice**: only a dilated AP has
  `smax=0` (good set on-lattice = deep well); every non-AP has `smax>0`, `coverGap≥1/13`. **Correction:
  the binding far element is the smallest multiple of `13` (26,39,…), not `182`** (covering 14 is not
  required). Exhaustive `max≤34` enumeration in progress (cover-gap criterion, all near-tight cores).
- This is the far-element (branch B) side of lever 3. **Branch (A) — the real wall — remains open**, and
  is the natural next target: it is where soft Weyl fails and where the difference-flow (lever 4) or a
  comparable-case position argument must take over.

## 7. Reconciliation: the wall is SPREAD, not comparable (death-star S56 verification)

A computational sweep this session sharpens **which** families are the wall, reconciling the "branch (A)
comparable" and "spread" framings:

- **Comparable families are LOOSE.** Comparable covering-`2..14` families with small speeds have
  `M ≥ 0.18` (min `2/11` at `{4,…,15,18}`, `136%` over `1/13`, `128k` families). Clustered `ρ≈1` families
  `{N..N+12}` are *maximally* loose: `M = 0.125, 0.31, 0.40, 0.45, 0.47, 0.49 → 1/2` as `N = 2…200`. So
  the "fully-comparable single-scale `ρ≈1`" region contains **no** near-counterexamples — its "wall" status
  is only that THM-1002's *bound* (`1/20`) is weak there; the *true* `M` is large.
- **The crux band `M<1/13` lives entirely in SPREAD families** `σ = v_max/v_min > 12`. This is immediate
  from klein **THM-1043** (`σ ≤ n−1 ⟹ M ≥ 1/n`, witness `t = 1/(n·v_min)`): with `n=13`, `σ ≤ 12 ⟹
  M ≥ 1/13`, with `{1..12}` (`σ=12`) tight at exactly `1/13` (witness `t=1/13` gives `min=1/13`). Verified:
  `σ=13` breaks the witness (the `σ`-speed hits `k/13`, `min=0`). So **`M<1/13 ⟹ σ>12`.**
- **Therefore death-star THM-1028's "branch (A) comparable" is essentially subsumed by THM-1043**, and the
  genuine residual is the *spread* stratum `σ>12`: by the inverse theorem this is `{dilated-AP core
  (σ_core=12)} + {far/nested element making σ large}` — i.e. **the deep-well / far-element structure
  (death-star cover-gap, lever 3) and the nested-scale renormalization (HYP-3901, lever 4)**, NOT a
  comparable single scale. The clean combined assembly: **THM-1043 kills `σ≤12`; the inverse theorem forces
  `σ>12` crux families to be AP-core + far-element; the cover-gap lemma shows the far element covers only
  the on-lattice AP (deep well); the renormalization tower handles nested scales.** The four levers of §4
  attack exactly this spread stratum. Scripts: `lrc_third_invariants`, `lrc_residual_covergap` (this session).

→ THM-1017, THM-1028, THM-1037, THM-1038, THM-1002/1008, THM-1042/1043, THM-1051/1061/1071, THM-1070,
HYP-7310, HYP-7505/7525, HYP-3901, the alignment/cover-gap/OCF reflections.
