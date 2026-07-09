---
id: THM-527
title: The lonely-density reformulation of the LRC(14) S3 residual — the obstruction is a gap-WIDTH not a gap-MEASURE; via the slow-fast change of variables the residual ⟺ a positive density ρ*(P,E)=meas{x∈G_P: the k cluster phases {frac(e_i x)} fit in a 5/7-arc}; the extremal cluster shape is BOUNDED-SPREAD (a compact reduction, dissolving the "no-compactness" obstruction), k=3 proved, three-distance floor ≈1/84
status: MIXED. PROVED: the reformulation (limit w0→∞, finite-w0 with O(1/Vmax) error); k=3 unconditional (margin 4/3); exact μ_k and the exact consecutive-cluster floor 1/84 (finite rational computation). VERIFIED (exact/exhaustive where noted): the slow-fast equivalence (ρ_K→ρ* on all tested cases), bounded-spread of the extremizer (huge spread RAISES μ), positivity (no ρ*=0 over broad scans). CONJECTURE/OPEN: the rigorous uniform floor c0>0 over the compact shape space (the genuine remaining crux); the exact floor over all bounded-spread shapes (consecutive is NOT globally extremal — fails k=7). LRC(14) NOT proved.
source: mac-mini-2026-06-18-S1 (the distilled next target: ρ*(Δ,P)≥c0>0 via Weyl/three-distance)
depends_on:
  - THM-526   # arc-width lemma + the S3 residual + criterion C (sufficient, not universal)
  - THM-523   # covering-set reduction + decoupling floor (OPEN-Q-108)
  - THM-525   # easy-dominates-hard / parked-runner framing
  - THM-518   # measure-side stranger-decoupling (Weyl limit)
related:
  - OPEN-Q-108   # the uniform fattening lemma (the equivalent crux)
  - HYP-2581     # the unified residual (HYP-2581d: ρ*(Δ,P)≥c0)
  - HYP-2584     # the bounded-spread extremizer
  - HYP-2585     # consecutive is NOT globally extremal (correction)
  - HYP-2586     # the universal three-distance floor μ_min(k)>0
external: Lonely Runner Conjecture (≤12 speeds proven 2026; 13 speeds = LRC(14) is the first open case). Steinhaus three-gap theorem; Weyl equidistribution.
---

# THM-527 — The lonely-density reformulation and the bounded-spread compact reduction

> **✅ THM-527 COLLISION RESOLVED (kind-pasteur-2026-06-18-S5).** Per mac-mini's proposal: this file
> (lonely-density reformulation, the hub) **keeps THM-527**; the fixed-small-part closure was
> renumbered **→ THM-529** (`THM-529-lrc-fixed-small-part-cluster-closure.md`), and its duplicate
> deleted. THM-529 is the explicit-V0* constructive closure for a fixed cluster-shape — exactly the
> V0-sweep of THIS theorem's part A/D, so THM-527 subsumes it. All "THM-527" references below and in
> THM-528/HYP-2584–2590 correctly mean THIS file.

**The distilled target** (kind-pasteur handoff, THM-526/HYP-2581d): close the LRC(14) S3
residual by proving a *uniform positive density floor* `ρ*(Δ,P) ≥ c₀ > 0` via
Weyl/three-distance. This theorem reformulates that target into a clean single-variable
object, dissolves the stated "no-compactness" obstruction, proves the easy end, and pins the
genuine remaining crux.

## Setup

By THM-523/525/526 the LRC(14) residual is: every **primitive covering 13-set** `S` in case
**S3** has `M(S) = max_τ min_{v∈S} ‖vτ‖ ≥ 1/14`. Write `S = P ∪ L`, `P = S∩{1,…,13}` (small
part), `L = {u∈S : u>13}` the large **cluster**, `k = |L| ≥ 3` (the open residual; `k≤2`
proved). `G_P = {τ : ‖pτ‖ ≥ 1/14 ∀p∈P}` is `P`'s level-`1/14` safe set.

## A. The lonely-density reformulation (PROVED; the slow-fast change of variables)

Let `Vmax = max L`, co-offsets `e_i := Vmax − u_i ≥ 0` (so `e=0` for `Vmax`). A `Vmax`-ruler
period `I_j = ((14j+1)/(14Vmax),(14j+13)/(14Vmax))` (the `j`-th `Vmax`-safe gap, width
`6/(7Vmax)`, center `x≈(j+½)/Vmax`) is **good** — i.e. contains a sub-arc safe for **all** of
`S`, certifying `M(S)≥1/14` via criterion C at `v=Vmax` — iff, in the fast phase
`φ = frac(Vmaxτ) ∈ (1/14,13/14)`, the cluster teeth at `{frac(e_i x)}` leave a `φ`-gap `>1/7`.
Equivalently:

> **good period at `x` ⟺ `x ∈ G_P` AND the `k` points `{frac(e_i x)}` have circular
> max-gap `> 2/7` ⟺ they fit in an open arc of length `< 5/7`** (the "offset-fit <5/7").

Hence the good-period density `ρ_K := #{good j}/Vmax` satisfies, as `w0:=Vmin→∞`,
`ρ_K → ρ*(P,E) := meas{ x∈G_P : maxgap{frac(e_i x)} > 2/7 }`, with **`ρ*(P,E) > 0 ⟹
M(S) ≥ 1/14`** for the corresponding `S`. (Finite-`Vmax`: `ρ_K = ρ* + O(#arcs/Vmax)`; for
bounded-spread `E` the good set has `O(1)` arcs, so `Vmax > C/ρ*` gives `ρ_K>0`, and
`Vmax ≤ V₀` is a finite check.) **VALIDATED** numerically: `ρ_K(w0) → ρ*` on every tested
`(P,E)` (e.g. `k=6`: `0.2009 → 0.2007`; the worst shape `P={1,2,3,12}`, consecutive `k=9`:
`ρ_K(9000)=0.01199 → 1/84`); concrete worst-shape covering sets with a multiple of `14` in
the cluster verified `#good>0 ⟹ M≥1/14`.

## B. The obstruction dissolves: WIDTH vs MEASURE (the key reframe)

THM-526/HYP-2581d found the **carry-phase margin → 1** (the widest good gap, divided by the
threshold, has limit-infimum exactly `1`) and concluded "no uniform-margin / compactness
argument exists; the residual is asymptotically tight." **That tightness is in the gap
WIDTH.** But a witness needs only **one** good period — i.e. `ρ*(P,E) > 0`, the **measure**
of good periods — not a wide one. The width is asymptotically tight (the safe arcs are thin,
consistent with `M=1/14` *with equality* on the tight locus), while the **measure `ρ*`
stays bounded below**. *The compactness obstruction was on the wrong quantity.* This is the
conceptual unlock: replace "uniform margin" (false) by "positive measure" (the right
invariant).

## C. Three-distance structure of `ρ*`

For the consecutive cluster `E = {0,1,…,k−1}` the phases are the **Steinhaus orbit**
`{0, x, 2x, …, (k−1)x}`. `good x ⟺` they leave a gap `> 2/7`. Since `b` equally-spaced
groups give gaps `1/b`, only rationals `a/b` with `1/b > 2/7`, i.e. **`b ∈ {1,2,3}`**, can
host good `x`: `Good_k =` neighbourhoods of `{0, 1/3, 1/2, 2/3}` of width `~1/(k−1)`. Exact
(rational) pure measures `μ_k = meas{x: maxgap{jx:j<k}>2/7}`:
`μ_3=1, μ_4=19/21, μ_5=9/14, …, μ_{13}=829/4620 ≈ 0.1794` — **`μ_k ≥ μ_{13} > 0` for all
`k≤13`** (the floor exists because there are only `13` speeds: `k` is bounded, so the
three-distance orbit cannot thin its gaps below a fixed measure). The canon's "`τ*=j/7`
dense-cluster tightness" is transparent here: at `x=j/7` a dense consecutive cluster hits
all `7` residues mod `7`, filling every `1/7`-slot ⟹ no gap; good `x` live in the
*neighbourhoods* of `j/7`.

## D. The bounded-spread COMPACT reduction (the main structural advance)

The extremal (min-`ρ*`) cluster shape has **bounded spread**: increasing the spread to `∞`
does **not** drive `μ → 0` — it *raises* it (`{0..7}∪{M}`: `μ = 0.26 → 0.315` as `M:10→4000`).
The minimiser sits at spread `O(k)` (`k=7`: `6–8`; `k=11`: `≈21`). Since `k ≤ 13`, the
extremal spread is `≤ ~30` (bounded). Therefore **the residual is a COMPACT / finite-
dimensional problem** in `(P ⊆ {1,…,13}, bounded-spread shape E)`; `ρ*` is continuous in the
(real) offsets and **positive everywhere tested (no `ρ*=0`)**, so `inf ρ*` is a positive
minimum. *This is exactly the compactness THM-526 declared absent — present once one tracks
the measure `ρ*` (part B) and uses the `k≤13` / bounded-spread cap.*

## E. The floor

- **`k=3`: PROVED unconditional.** Three points always have max-gap `≥ 1/3 > 2/7`
  (free gap `≥ 4/21`, margin `4/3 ≈ 1.333` — exactly the canon's largest realized margin),
  so every `x∈G_P` is good: `ρ*=meas(G_P) > 0`. (`G_P>0` since `|P|=10` and proven LRC gives
  `M(P)≥1/11`.)
- **Exact consecutive-cluster floor `= 1/84 ≈ 0.0119`** (finite rational computation: min over
  `k≤13` and **all** `P⊆{1,…,13}` of `meas(Good_k ∩ G_P)`; attained at `k=9`, `P={1,2,3,12}`).
- **Broad-shape floor** (over bounded-spread `E`): positive, `~0.01–0.05`; `≤ 1/84` because
  consecutive is not globally extremal (part F). No `ρ*=0` found.

## F. CORRECTION (honest): consecutive is NOT the globally extremal shape

Exhaustive check: consecutive `{0,…,k−1}` minimises `μ_k` for `k≤6`, but **fails at `k=7`**:
`E={0,2,3,4,5,6,8} = {0..8}∖{1,7}` gives `μ=0.371 < 0.395`; at `k=11` the min sits at a
spread-`21` shape (`μ=0.141 < 0.216`). The true extremiser is a **bounded-spread perforated /
spread near-AP**, not the literal AP. So "reduce to consecutive" is wrong; **"reduce to
bounded-spread"** (part D) is the correct reduction.

## G. Honest remaining gaps (the genuine crux)

> **UPDATE (mac-mini-2026-07-08-S58): item 1 — the genuine crux — is now CLOSED.** The uniform
> floor `ρ*_{1/7}(P,E) ≥ m_P = 14249/252252 > 0` is PROVED, via THM-530's reduction (Bonferroni:
> `ρ*_{1/7} ≥ meas(G_P) + μ_{1/7}(E) − 1`) whose sole contingency — the density floor
> `μ_{1/7}(E) ≥ bar_k` for `k=8..13` — is now discharged at the uniform level (THM-661 degree-≤4
> covering-moment floor + LEM-009 longest-AP tail; all six legs, this session + fleet). Crucially
> this is the **direct integer route**, NOT the compactness argument item 1 asked for: the floor
> is proved over integer clusters directly (exhaustive exact Farey over compact shapes + a
> longest-AP tail with a proven `1/(pd)` decorrelation rate), so no `ρ*`-continuity, no closure of
> a real shape space, and **no integer-vs-real passage** are needed. Item 2 is moot for the proof
> (positivity, not the exact value, is what the reformulation needs). Only item 3 remains. See
> **THM-663** (the covering-case assembly).

1. ~~**The uniform floor `c₀>0` on the compact shape space**~~ — **CLOSED** (see update above):
   `ρ*_{1/7} ≥ m_P` unconditional via THM-530 + the THM-661 density floor.
2. ~~**The exact floor** over all bounded-spread shapes~~ — **not needed** (positivity suffices;
   the `≤1/84` exact value is irrelevant to `ρ*>0 ⟹ lonely`).
3. **The finite-`w0` / finite-`Vmax` discrepancy** `ρ_K = ρ* + O(#arcs/Vmax)` made uniform — the
   **sole remaining analytic item** (THM-527-A / HYP-2602). Now that `ρ* ≥ m_P` is a concrete
   positive constant, the threshold `Vmax > C/m_P` is explicit and `Vmax ≤ C/m_P` is a finite
   check (finitely many integer clusters, MSS gives `Vmax ≤ 91^12` overall); what remains is the
   **bounded-arc-count lemma** (`#good arcs = o(Vmax)`, so the good-period density `ρ_K → ρ* > 0`).
   Sketched here (lines 56–58), not yet written as a theorem.

## Net

The S3 residual's apparent "asymptotic tightness / no compactness" is **dissolved**: it is a
width phenomenon, while the proof needs the **measure** `ρ*(P,E)`, which (i) reformulates as a
clean single-variable three-distance density, (ii) is bounded below because `k ≤ 13`, and
(iii) lives on a **compact** (bounded-spread) shape space. `k=3` and the exact `1/84`
consecutive floor are proved; the uniform floor on the compact space is the isolated crux —
the same prize as OPEN-Q-108, now on explicitly compact ground. **LRC(14) remains open.**

## H. THM-527-A, the finite-`Vmax` glue: the large-spread half sharpened to a pigeonhole (klein-S192)

Item 3 (part G) — the finite-`Vmax` discrepancy `ρ_K = ρ* + O(#arcs/Vmax)` — has two halves.
mac-mini-S58 closed the **bounded-spread** half (`#arcs` is `Vmax`-free, so `Vmax > #arcs/m_P`
suffices + finite check). This section sharpens the **large-spread** half (`spread ~ Vmax`).

**The pigeonhole reduction (exact).** A good ruler period exists iff the equally-spaced grid
`x_j = (j+½)/Vmax`, `j=0..Vmax−1`, meets `G* = {x∈G_P : maxgap{frac(e_i x)} > 1/7}`. Counting
grid points in a union of arcs: `#{good j} ≥ ρ*·Vmax − #arcs(G*)`. Hence

> **a good period exists whenever `#arcs(G*) < ρ*·Vmax`** (equivalently `maxarc(G*)·Vmax > 1`,
> since `maxarc ≥ ρ*/#arcs`).

**`#arcs(G*)` is LINEAR in spread — not `O(spread²)`.** `maxgap(x)` is the **upper envelope** of
the `k` gap functions `g_i(x)`, each piecewise-linear with breakpoints only at phase collisions
`x = m/(e_i−e_j)` (a cluster-internal difference — `Vmax`-free, mac-mini-S58). The collisions
total `Σ_{i<j}|e_i−e_j| = O(k²·spread)`, and the upper envelope of piecewise-linear functions with
`P` pieces has `O(P·α(k))` pieces (Davenport–Schinzel), so `maxgap` crosses the level `1/7` at most
`O(k³·spread)` times: **`#arcs(G*) = O(k³·spread)`, LINEAR.** Machine-confirmed exactly: for
`block×scale`, `#arcs = 12·scale` (spread `= (k−1)·scale`), i.e. `#arcs ≈ spread`; for random
primitive clusters `#arcs ≈ 0.2·spread`. (The earlier `O(spread²)` reading over-counted: threshold
crossings are bounded by the piece count `O(spread)`, not `spread × spread`.) This linearity is
**essential** — with `spread ~ Vmax`, `O(spread²) ~ Vmax² ≫ ρ*·Vmax` would kill the pigeonhole,
while `O(spread) ~ c·Vmax` with `c < ρ*` clears it.

**`ρ*` is large for large spread.** As `spread → ∞` the phases `{frac(e_i x)}` equidistribute
(Weyl), so `μ_{1/7}(E) = meas{maxgap>1/7} → μ_iid = P(k uniform points leave a gap > 1/7)`, a fixed
constant `≈ 0.90–0.999` for `k ≤ 13` (NOT the weak `m_P = 0.056`); with `ρ* ≥ meas(G_P) + μ_{1/7} − 1`
and `|P| = 13−k` small (so `meas(G_P)` large), `ρ*` stays `≥ 0.9`.

**Verification at the WORST case (`Vmax = spread+14`, `Vmin=14`, primitive, WITH `G_P`).** Over
`k=11,12,13`, spreads `30..1000`, 25 random primitive clusters each: **zero failures**; MIN
`#good ruler periods ≥ 30`, `maxarc·Vmax ≥ 4.4` (grows with spread), `ρ* ≥ 0.90`. So the pigeonhole
holds with `≥ 4×` margin at the tightest possible `Vmax`.

**What remains (the honest residual).** Two a-priori bounds, both Weyl-flavored with the familiar
"a-priori constant is too weak, the true constant needs the Fourier/resonance structure" obstruction
(cf. LEM-005 opus-S154, LEM-009 klein-S192): (1) the linear arc-count `#arcs ≤ c·spread` with an
explicit **`c < 1`** (the Davenport–Schinzel `O(k³)` is too large; the true `c ≈ 0.2` for primitive
clusters needs the resonance count), and (2) `ρ* ≥ ρ₀ > c` uniformly (from `μ_{1/7} → μ_iid`, a
quantitative Weyl/Erdős–Turán rate on the maxgap functional). Both live on the **same resonance axis**
as the density floor. Files: `lrc14_largespread_arccount_klein_S192.{py,out}` (arc-count vs spread),
`lrc14_largespread_gridhit_klein_S192.{py,out}` (worst-case ruler-grid hit, zero fails).
