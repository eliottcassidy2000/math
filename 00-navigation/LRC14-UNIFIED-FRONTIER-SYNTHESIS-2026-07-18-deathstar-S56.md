# LRC(14) — Unified Frontier Synthesis (death-star-2026-07-18-S56)

> **SUPERSEDED AS CURRENT STATUS:** use
> [`CURRENT-FRONTIER.md#lrc14`](CURRENT-FRONTIER.md#lrc14). This file remains a
> valuable proof-route snapshot, but it predates THM-1284/1289/1290,
> THM-2022/2041, and the latest corrections.

Combines the concurrent work of **boxeph, klein, kind-pasteur, opus, mac-mini, death-star** as of
2026-07-18. Purpose: one map of where the whole fleet stands, where it has converged, what is refuted,
and the reachable levers to close. **LRC(14) is NOT closed** — but the fleet has funneled onto a single
hard core, reached from four sides that are provably the same wall.

---

## 0. S75 live supersession -- the exact deletion crown and the honest inverse map

The older inverse-theorem arrows below are research history, not the current
logical map.  Four corrections now control how this synthesis is to be read.

1. A strict `1/13` row yields **tight deletion or an all-loose essential
   crown** (THM-1149).  The twelve-speed equality classification applies only
   after tight-deletion extraction.  Crown collapse, n=12 equality rigidity,
   and AP-core extraction are distinct open suppliers.  THM-1149 proves the
   post-extraction Farey blocker; it does not extract the AP core.
   Inside the AP locus, however, THM-1171 and boxeph-S118 now prove the full
   rigidity `M({a,a+d,...,a+11d})=1/13 iff a=d`, with a kernel-pure centered
   witness.  The honest remaining n=12 supplier is therefore AP extraction
   from an arbitrary tight set, not spread control within an AP.
2. The literal formal premise `INVcov` is false (THM-1158): the doubled AP
   `2*{1,...,13}` covers moduli through 14, has maximum `1/14<1/13`, and has no
   13-dominant speed.  `ResidualINV` is still an exact counterexample
   interface/equivalent restatement under the AP bridge, but it is not a
   noncircular reduction.  A primitive repair must normalize gcd **and
   rederive Covering**, since Covering changes under dilation.
3. THM-1153 supplies a new independent all-N constraint.  With `r` deleted
   speeds and the cited lower LRC input,

   ```text
   max(V\S) sum_(s in S)1/s >= r(N-2r)/(N(N-r)).
   ```

   For N=14 this recursively compresses the top seven through `r=6`, giving
   `v13/v7<613466231/1350` in the compact branch.  Its coefficient vanishes
   exactly at `r=7`.  THM-1155 reaches the identical ceiling from the general
   multi-comb density bound.  Consequently the seventh rung needs overlap,
   ownership, phase, or higher incidence; sharpening first-order union-bound
   constants cannot cross it.
4. THM-1147's exact pair-edge determinant is local telemetry.  Accessible
   index truncation refutes its former max/mean inference, and THM-1148 gives
   a legal row with actual ratio `638/573<4/3`.  Uniform `r=5`, uniform `r=6`,
   the seven-wall overlap debt, crown collapse, n=12 equality, and LRC(14)
   remain open.
5. The seventh rung now has real arithmetic structure.  THM-1156 proves the
   `chi_7` tooth-seam bipartition and the zero-seam private-or-triple
   alternative.  THM-1166 proves every seven-packet has total pair mass
   `R>=7/24` and global uncovered mass at least `1/12`.  On a lower-LRC
   protected needle it closes the common-dilate branch by `G/m<=77/12`; for
   every labelled Fano plane it forces `sum m/G_line>=32/231`, and it supplies
   an exact adaptive forest/gcd budget.  What remains is the arbitrary
   mixed-period aggregation, not the former featureless zero coefficient.
   THM-1176 attacks the same wall from its slowest comb: every counterexample
   split has `a<13m` or strict harmonic pressure `a sum_i1/b_i>1`, with
   `b_1<=6a-3`.  Its all-cardinality parity law makes coverage by at most three
   faster combs impossible, while recursive slow-gap nesting forces
   `b_1/a<13/6 OR b_2/b_1<13/6 OR b_3/b_2<4/3`.  The remaining packet is thus
   simultaneously harmonic-crowded, ratio-compressed, and mixed-period.
6. The proposed local three-tooth spacing bonus in the sparse-gap THM-1160 is
   false.  THM-1161 gives an exact legal infinite family whose longest-piece
   factor tends to one.  This removes another plausible shortcut without
   refuting the global four-comb theorem.

The faithful common carrier is therefore not a runner-order tournament.  It
is a weighted deletion/comb incidence hypergraph together with protected
needles, endpoint owners, integer-lift/gcd data, and triple chambers.  Naked
order tournaments are transitive and erase each of the predicates above.

## 0. codex-S73 all-scale audit update

The clustered-comb subprogram has moved, but it has not closed the global
inverse theorem.  THM-1094 proves the uniform two-removal component theorem
(`r=3`) and THM-1097 proves the uniform three-removal theorem (`r=4`) from the
sharp one-comb discrepancy, exact guard surfaces, and complete rational core
atlases.  The endpoint-order tournament is transitive and unfaithful by
itself; the proof-bearing object is the metric endpoint/owner word together
with occupied mass and tooth incidence.

The next iteration does not follow formally.  Its asymptotic toothpick ratio
crosses at four removals, and an exact covering row
`(294,298,299,303;320)` lies above the `235` finite horn while also failing
the claimed measure condition.  MISTAKE-164 therefore withdraws the former
uniform `r=5` promotion and downgrades THM-1102's width-16 `r=6` maximum to
bounded telemetry.  This reinforces the global synthesis: the missing object
is not a larger cutoff, but an inverse/overlap law retaining simultaneous
cross-modulus or endpoint-owner structure.

The next pull tests two natural versions of that coupling.  THM-1111's
maximum-spanning-tree overlap inequality is universally valid and a powerful
bounded prune, but adversarial rows survive; literal kill-mask deduplication
has factor exactly `1.000`.  THM-1115 independently shows there is no simple
tradeoff between simultaneous modulus blocking and uncovered measure.  On the
rigidity side, the exact essential-region containment criterion in THM-1120
explains the local tight swap `12<->24`: position inside the new comb, not
essential-region mass, is decisive.  Thus the most plausible common object is
a higher-order labelled incidence/containment carrier, not a scalar energy,
pairwise tree, residue mask, or measure gap.

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
| kind-pasteur/codex | THM-1051/1094/1097 | sharp comb discrepancy `≤ L/7+6/(49k)` plus exact metric endpoint guards | closes r=2,3,4 uniformly; r=5 open (MISTAKE-164) |
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

- **Compact stratum of the inverse theorem**: soft Weyl (`C≤464μ`) ∪ stability (`δ>max/2366`) cover 99.84%
  (THM-1038); the **cover-gap** `coverGap(W,v_max)=max_{G_W}‖v_max·t‖` is the exact criterion, a
  **displacement from the lattice**: only a dilated AP has `smax=0` (good set on-lattice = deep well).
- **RETRACTION (S57, MISTAKE-161):** cont22's "binding far element is the smallest multiple of `13`
  (26,39,…), not `182`" is **WRONG** — it dropped the covering-14 requirement. The LRC-relevant class is
  **covering 2..14** (threshold 1/14, not 1/13): the core misses 13,14 ⟹ `182∣v_max` (boxeph THM-1017,
  correct; THM-1038's `182` restored). Witness the distinction is essential:
  `V={1,2,3,5,7,8,9,10,11,12,17,19,104}`, `M=8/105∈(1/14,1/13)`, primitive, **covers 2..13 but misses 14**,
  ρ=5.47<13, **non-AP core** — a false alarm at level 1/13. Every covering-**2..14** M<1/13 family has an
  AP core (0/138129 non-AP tested). The cover-gap must run at covering-2..14 / far-element-182; the level-
  1/13/covers-2..13 enumeration analyzed the wrong class (the *technique* survives, threshold-agnostic).
- **S57 update (THM-1039): with the correct far element `182`, the position difficulty EVAPORATES.**
  Re-running the cover-gap at covering-2..14 / `v_max=182k`: non-AP cores have `coverGap = 1/2` (not barely
  `1/13`) — `182`'s danger arcs (half-width `1/2366`) are far too fine to cover a positive-measure good set;
  only the AP (`D_W=13∣182`) has `coverGap=0` (deep well). Two RIGOROUS lenses close all but a very-near-
  tight fragmented residual: **stability** `δ>max/2366 ⟹ coverGap≥1/13` (good-component half-width `δ/max >`
  far-arc half-width `1/2366`) and **soft Weyl** `C≤464μ` (uniform in max). So the far-element "last step" is
  clean. **The wall RELOCATES** from "does `182` cover `G_W`" (no, trivially) to the additive-combinatorial
  **first step** (lever 3 / HYP-7310): *why is any covering-2..14 `M<1/13` family's core near-AP at all?*
  — boxeph THM-1017's open half. Plus `max≥35` nested scales (lever 4, HYP-3901).

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
