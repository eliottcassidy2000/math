# A dictionary of maps on the loop: the runners are multiplication, the OPUC recursion is Blaschke, the Dedekind cocycle connects them — and the census is the quotient

*kind-pasteur-2026-07-01-S8. The owner's cue: "pushing Verblunsky to the unit circle is a nice recursive metaphor for LRC (runners on a loop) — find creative functions between points on a loop, a whole dictionary; they may operate group-like." They do. The circle `T = ℝ/ℤ` carries a tower of natural group actions, every LRC object is an orbit/invariant/cocycle of one of them, and the 11-speed census is literally the quotient by the dilation×complement group. Verified computationally; then synthesized into the census next step.*

## The dictionary (and which entries are groups)

| # | map | formula on `T=ℝ/ℤ` | composition | structure | LRC role |
|---|---|---|---|---|---|
| 1 | rotation `R_a` | `t ↦ t+a` | `R_a∘R_b=R_{a+b}` | **group** `(T,+)` | the runner **flow** |
| 2 | dilation `M_v` | `t ↦ v·t` | `M_v∘M_w=M_{vw}` | monoid `(ℤ,×)`; **group `(Z/N)*`** on `N`-torsion | the **runners** themselves |
| 3 | affine `A_{v,a}` | `t ↦ v·t+a` | semidirect | **group** `Aff = T⋊(Z/N)*` | a runner with a phase |
| 4 | inversion `ι` | `t ↦ −t` | `ι²=id` | **group** `ℤ/2` | antipode / complement (Borsuk–Ulam) |
| 5 | doubling `D` | `t ↦ 2t` | `D^k=M_{2^k}` | monoid `(ℤ_{≥0},+)` | the "2" of `14=2·7`, dyadic tower |
| 6 | Blaschke `b_a` | `z ↦ (z−a)/(1−ā z)` | Möbius∘Möbius | **group** `SU(1,1)/±` | the **OPUC/Verblunsky** recursion |
| 7 | Gauss `G` | `x ↦ {1/x}` | natural extension | groupoid (modular flow) | `t*=[0;13,14]`, three-gap |
| D | character `χ_k` | `t ↦ e(kt)` | `χ_k∘χ_l=χ_{k+l}` | **dual group** `Ẑ=ℤ` | Fourier; Ramanujan sum = `Σ` unit characters |
| c | sawtooth `((·))`, `s(h,k)` | 1-cocycle of `SL₂(ℤ)` | reciprocity | **cocycle** | the connector; value `→ −1/12 = ζ(−1)` |

Three of these carry the whole story, and all three are verified:

**Runners are map #2, and the atoms are its group.** `M_v∘M_w = M_{vw}` is literally "run at speed `v`, then at speed `w` = run at speed `vw`." On the `N`-torsion points the invertible runners are exactly `(Z/N)*` — a genuine abelian group (verified closed/invertible), and the LRC atoms of an AP core (HYP-3793) are one of its orbits. The group *type* varies: `(Z/10)*≅ℤ/4` (cyclic, gen 3), `(Z/12)*≅(ℤ/2)²` (Klein four!), `(Z/14)*≅ℤ/6` (cyclic). Scale-invariance (THM-522) is nothing but the action of `M_c`.

**The OPUC recursion is map #6 — Blaschke composition — and "push to the unit circle" is its boundary dynamics.** Verified: `b_{a₁}∘b_{a₂}` is again a disk automorphism (the maps form the group `SU(1,1)/±`), and iterating a fixed `b_a` sends the orbit of `0` to the boundary — `0 → 0.5 → 0.80 → 0.93 → 0.98 → 1` for `|a|=0.5`, faster as `|a|→1`. This is the owner's metaphor made exact: the Szegő recursion *is* Blaschke composition, and the extremizer's `|α_k| → 1` (HYP-3793) is the orbit hitting the boundary = the loose measure collapsing onto its `(Z/N)*` atoms. **Tight = atomic = the Verblunsky orbit reaching the circle.**

**The Dedekind sum (map c) is the connector between the additive circle (#1) and the modular/Möbius world (#6/#7)** — which is *why* the census collar-sum (additive, over units) evaluates to a Dedekind sum (`−1/12 = ζ(−1)`, HYP-3773). Verified reciprocity `s(h,k)+s(k,h) = −1/4 + (h/k+k/h+1/hk)/12`. And a clean gem fell out:

> **`s(2,7) = 1/14` and `s(3,7) = −1/14`** — the LRC(14) threshold is the apex-prime-7 Dedekind quantum. For `k=7` prime, all `s(h,7) ∈ (1/14)ℤ`; the threshold `1/14 = 1/(2·7)` is the minimal nonzero value. Same source `14 = 2·7`, two faces: the danger band and the Dedekind quantum coincide.

## The census is the quotient by the loop-group

Because dilation `M_c` and complement `ι` are symmetries (verified: `M_2,M_3,M_5` fix `M(S)=1/10` and the primitive class), **the space of loose 11-cores is a quotient**, and the `(Z/N)*` label is the invariant that names each orbit. So the census reduces to enumerating orbit representatives — one per modular class. Doing that over a broadened pool (**6846** loose 11-cores: all 2-drops of `{1..16}`, 1-drop+far to 60, random to speed 22):

- **Global minimum `= 313/9702 = 0.03226` at `{1..13}\{6,10}` — the full orbit `(Z/10)*`. Every core clears `1/36`.** (Last turn's extremizer survives a 90× larger search.)
- **Organized by class `N`, every modular class clears `1/36`:**

| `N` | `(Z/N)*` | `φ` | min `meas` (full-orbit) | ≥ 1/36 |
|---|---|---|---|---|
| **10** | ℤ/4 | 4 | **0.03226** | ✓ |
| 6 | ℤ/2 | 2 | 0.04413 | ✓ |
| 8 | (ℤ/2)² | 4 | 0.05335 | ✓ |
| 12 | (ℤ/2)² | 4 | 0.05633 | ✓ |
| 9 | ℤ/6 | 6 | 0.06655 | ✓ |
| 7,5,4,3,18 | … | … | 0.082–0.163 | ✓ |

- **The near-tight locus is the dichotomy `{full unit orbit, partial two-clash}`** (bottom-30: 10 full-orbit, 19 two-clash, 1 other). The binding competition is `(Z/10)*` (`0.03226`) vs the two-clash `(Z/19)` core `{1,2,3,5,7,8,9,10,11,12,13}` (`0.03238`) — a **0.4% gap**. This is the 11-speed echo of THM-523's tight-locus dichotomy `{AP dilate, Goddyn–Wong sporadic}`: the full orbit is the "AP," the two-clash is the "GW."

## What this buys, honestly

**The group structure is real and verified**: runners compose multiplicatively into `(Z/N)*`, OPUC composes into `SU(1,1)`, and the Dedekind cocycle connects them, with the threshold `1/14` appearing as the apex-prime Dedekind quantum. The census is genuinely a quotient by the dilation×complement group, and the `(Z/N)*` label organizes it into finitely many modular classes, each of which clears `1/36` in a broad search — with the infimum pinned at `(Z/10)*` and a two-clash `(Z/19)` competitor a hair behind.

**It is not yet a proof.** The census is a broad *search* (6846 cores, up to speed 22), not an exhaustive bound over all primitive 11-cores; and the dichotomy is *dominant, not perfectly clean* (one "other" core in the bottom-30). But the group picture sharpens the target into something finite: prove the `{full-orbit, two-clash}` dichotomy is exhaustive on the near-tight locus (THM-522's scale-invariance already bounds representatives to scale-1, THM-560 gives the enumeration template), then evaluate the two families in closed form — the full-orbit one as a `(Z/N)*` collar sum minimized at `N=10`, the two-clash one as a Diophantine collision minimized at `N=19` — and confirm both clear `1/36`. That is now a two-family closed-form check, not an open search.

## Convergence (the frontier arrived at the same dictionary independently)

Fetching after this work, three agents had just pushed the same picture from different sides — a strong triangulation that the loop-group structure is the right frame:

- **opus-S13** built the *same* dictionary from the **position-space** (runner-cloud) side and explicitly named my S7 work its "time-space dual." Their `lrc_circle_dictionary` groups it as **AGL(1) {rotation, dilation=RUNNER, reflection, affine}**, **PSL₂(ℤ) {Gauss, mediant, Möbius}**, **Szegő {Verblunsky vertex-insertion}**, **analytic {sawtooth (ι-odd), distance (ι-even)}** — matching my #1–#7 with the same "runners = dilation part of AGL, units = invertible dilations" reading. Their signature result is the **harmonic law**: the AP `{1..n−1}` at `t=1/n` has `|α_j| = 1/(n−1−j)` *exactly* (verified n=5..20), with `α_0 = 1/(n−1)` the runner-centroid (their observer lens) = the covering-min ceiling, and `α_{n−2}=1` the terminal atomic value. This is the **position-space** companion to my **measure-space** OPUC thermometer: opus's `|α_j|` are the AP runner cloud's balance tower; mine are the lonely *measure's* moments — same Szegő machine, the two sides of the moment problem.
- **mac-mini-S78 (HYP-3794): the open core is FINITE** — MSS 2024 (arXiv:2411.06903) caps any LRC(14) counterexample's speeds at `C(14,2)^12 = 91^12`. So "unbounded is the difficulty" is obsolete; the census over modular classes is *doubly* finite (by MSS *and* by the `(Z/N)*` organization).
- **klein-S68 (HYP-3800)** and **mac-mini-S77 (HYP-3792)** independently reduced the far-coupling to finite **modular** invariants — the phase-residue `p(w)=nw mod Φ₆ ∈ ℤ/Φ₆` and the safe-band `(ℤ/q)` frame. Same multiplicative-group backbone.

So the whole frontier has converged on: *the LRC(14) open core is finite and modular*, living on `(Z/N)*`, `ℤ/Φ₆`, `ℤ/q`. My net-new pieces in that shared frame are **(i) the Dedekind-cocycle connector with the threshold identity `1/14 = s(2,7)`**, and **(ii) the census synthesis by modular class** (the `{full-orbit, two-clash}` dichotomy and the per-`N` table clearing `1/36`).

— Related: HYP-3793 (atoms = `(Z/N)*`, moments = Ramanujan sums, ζ as unit-restricted denominator — this reflection extends it), opus-S13 (circle-map dictionary + harmonic law `|α_j|=1/(n−1−j)`, the position-space dual), mac-mini HYP-3794 (open core finite via MSS), klein HYP-3800 / mac-mini HYP-3792 (phase-residue / safe-band modular frames), HYP-3773 (margin = 2-fold Dedekind sum → `−1/12 = ζ(−1)`), THM-522 (scale-invariance = the `M_c` symmetry), THM-523 (`{AP, Goddyn–Wong}` dichotomy — the two-family echo), THM-560 (census template), HYP-3762 (Gauss/CF map, three-gap), OPEN-Q-108. Scripts: `04-computation/lrc14_loop_map_dictionary_kps.py`, `lrc14_census_dichotomy_synthesis_kps.py` (+ `.out`). Extends HYP-3793; converges with opus-S13 — **no new HYP reserved** (halting the collision cascade the agents flagged).
