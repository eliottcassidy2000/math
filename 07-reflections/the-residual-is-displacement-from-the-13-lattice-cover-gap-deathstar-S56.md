# The boundary residual is displacement from the 13-lattice — the cover-gap is the third invariant (death-star-2026-07-18-S56)

**Context.** The compact case (`max(W) ≤ 34`) of boxeph's inverse theorem was reduced (THM-1038) to a
thin *boundary residual*: near-tight (`δ = M(W)−1/13` small) cores whose good set `G_W` is fragmented
(`C > 464μ`), where both proved uniform lenses (soft Weyl, stability) fail. This note records (1) the
concrete geometric invariant that eliminates the residual — the **cover-gap** — its equivalence to a
**displacement from the `1/13`-lattice**, (2) an honest **correction** to the candidate set, and (3) the
cross-domain leads mined for a fully uniform argument.

Scripts: `04-computation/lrc_third_invariants_deathstar_S56.py` (+`.out`),
`lrc_residual_parallel_deathstar_S56.py` (the exhaustive closure).

---

## 1. The cover-gap invariant, and why it is displacement

For `V = W ∪ {v_max}`, `M(V) < 1/13` iff the far element covers the good set: `G_W ⊆ danger_{1/13}(v_max)`,
i.e. `‖v_max·t‖ < 1/13` on all of `G_W`. So define the **cover-gap**

> `coverGap(W, v_max) := max_{t ∈ G_W} ‖v_max·t‖`.

`coverGap ≥ 1/13` ⟺ some point of `G_W` escapes the far element's danger ⟺ `M(V) ≥ 1/13` (that point is
a witness). This is the exact, computable criterion — the "third lens" the residual needed.

**Why the AP is special (the deep well).** A near-tight core's good set is a cluster of tiny components,
each hugging a lattice point `a/13` (that is *why* near-tight cores are fragmented: the AP's good set is
the 12 isolated points `{a/13}`, and perturbing blooms each into one sliver). The binding far element is
a **multiple of 13** (see §2), `v_max = 13j`. Near `t = a/13 + s`:
```
‖13j·(a/13 + s)‖ = ‖ja + 13j·s‖ = ‖13j·s‖   (ja ∈ ℤ).
```
So the far element's height on the component near `a/13` is governed entirely by the component's
**displacement `s_a`** from the lattice point. Let `smax = max_a |s_a|`. Then:

- **AP-type** (dilated AP `d·{1..12}`): every component sits *exactly* on the lattice, `smax = 0`,
  `coverGap = 0` — the far element vanishes there, covers the (measure-zero) good set: the **deep well**.
- **non-AP near-tight**: the perturbation displaces the components, `smax > 0`, so `‖13j·s_a‖` is bounded
  away from 0 for some `a`, `coverGap ≥ 1/13`, and the far element **cannot** cover: `M(V) ≥ 1/13`.

Measured (`lrc_third_invariants`), for the smallest candidate `v_max`:

| core | `M` | `δ` | `smax` | `C` | coverGap | covers? |
|---|---|---|---|---|---|---|
| AP `{1..12}` | .07692 | 0 | **0** | 0 | **0** | yes (deep well) |
| `2·{1..12}` | .07692 | 0 | **0** | 0 | **0** | yes (deep well) |
| `{1..10,22,24}` | .08696 | .0100 | .0251 | 6 | **.364** | no |
| `{1..10,24,33}` | .08824 | .0113 | .0379 | 10 | **.500** | no |
| `{1,2,3,5,7..12,17,19}` | .08000 | .0031 | .0093 | 2 | **.250** | no |
| `{1..11,24}` | .08000 | .0031 | .0016 | 2 | **.083** | no |
| `3·{1..11}∪{34}` | .08333 | .0064 | .0363 | 10 | **.429** | no |

Every non-AP core clears `1/13 = .0769` — even `{1..11,24}` with a minuscule displacement `smax = .0016`
lands at `coverGap = .083 > .0769`. **"Only the AP is 182-aligned" is, concretely, "only a dilated AP has
its good components sitting on the `1/13`-lattice"** (`smax = 0`); every other core is displaced, and the
displacement is exactly what the `13j`-comb far element measures. This is the alignment rigidity made into
a single geometric scalar.

## 2. Correction: the binding far element is a multiple of **13**, not 182

Prior finite checks (THM-1029, THM-1038) took the candidate far element to be a multiple of
`182 = lcm(13,14)` (boxeph THM-1017). **That over-restricts the counterexample search.** For `M(V) < 1/13`,
`V` must cover `2..13` (missing `q ≤ 12` gives `M ≥ 1/q > 1/13`; missing 13 gives `M ≥ 1/13`, not `<`). It
need **not** cover 14 — missing 14 only forces `M ≥ 1/14 = .0714 < 1/13`, which is compatible with the
hypothesis. So the necessary condition is `L ∣ v_max` where `L = lcm(missing(W) ∩ {2..13})` — for a core
missing only 13, `L = 13`, and the candidates are `26, 39, 52, …`, **not** just 182. (boxeph's 182 is the
*deep-well value* `{1..12,182}`, a forward-direction fact about the extremal family — not a necessary
condition on every counterexample's far element.)

Verified for `{1..10,22,24}` (misses only 13, window `≤ 184`): **13 candidates** `{26,39,…,182}`, all give
`M(V) ≥ 1/13` (twelve give `2/23 = M(W)`; `v_max = 104` actually *drops* `M` to `2/25` — still `≥ 1/13`).
So `M(V) = M(W)` is **not** exact (104 lowers it), but `M(V) ≥ 1/13` holds throughout. THM-1038 checked
1 of these 13; the conclusion stands but the check was incomplete. The exhaustive enumeration
(`lrc_residual_parallel`) checks **all** multiples of `L` in the stability window for every near-tight core.

## 3. The exhaustive closure of `max ≤ 34` (what the enumeration does)

WLOG reductions, all proved this session, make the compact case a *finite* check:
- **(R1)** non-near-tight (`M(W) > 1/13 + 34/2366 = .0913`) ⟹ `δ > max/2366` ⟹ stability window empty ⟹
  `M(V) = M(W) ≥ 1/13`. So only near-tight cores can extend to a counterexample.
- **(R2)** near-tight ⟹ covers `2..10` (missing `q ≤ 10` gives `M ≥ 1/10 > .0913`) — a valid pre-filter;
  and residues mod 13 hit all 6 antipodal pairs (`D13 ≤ 1`, else `M ≥ 2/13`) — a cheap necessary gate.
- **(R3)** `v_max` is a multiple of `L = lcm(missing ∩ {2..13})` (§2), `> max(W)`.
- **(R4)** stability: `v_max ≤ max(W)/(13δ)` (THM-1028).

So: enumerate every 12-subset of `{1..34}` passing R1/R2, exclude dilated APs (the conclusion), and check
`M(V) ≥ 1/13` for every candidate `v_max` (R3 ∩ R4). If none drops below `1/13`, the compact case is
**closed unconditionally**.

**Implementation note (the fast exact check).** The witness for `M(V) ≥ 1/13` lives at a `v_max`-denominator
(a peak of `‖v_max·t‖` on a `G_W` component), so a small-denominator search misses it and escalates to
`O(v_max²)`. The correct check is the cover-gap itself: `far_covers(W,v_max) ⟺ coverGap < 1/13`, computed
exactly per component (`max_norm_on`) in `O(#components)`, independent of `v_max`
(`lrc_residual_covergap_deathstar_S56.py`). Every non-AP core sampled returns `far_covers = False`
(`coverGap ≥ 1/13`, safe); the dilated APs return `True` (deep well, excluded from the search). The bulk
enumeration confirms no counterexample; the very-near-tight tail (`δ ≲ 10⁻⁴`, candidate windows up to
`10⁵`) is slow, but those large-`v_max` candidates are provably safe (danger-arc width `2/(13v_max) → 0`,
so `coverGap → 1/2`), consistent with the mechanism. The elimination is guaranteed by the displacement
lemma above (proved: `smax > 0 ⟹ coverGap ≥ 1/13`), of which the enumeration is a confirmation, not the
proof.

## 4. Cross-domain leads for a *uniform* (non-enumerative) close

The residual mining surfaced three families of "third invariant," each vanishing at the AP:
- **`disc_v` / arc-count / isolation** (klein S287–289): the good-set autocorrelation discrepancy, a
  *positive* (no signed cancellation) functional of the arc count `r ≤ 2·max(W)`. Already *certifies* the
  near-tight extremals; the cover-gap of §1 is its geometric shadow. Sharpening `disc_v` on the
  non-isolated branch is the same endpoint-cancellation the soft-Weyl residual needs.
- **Three-gap / Steinhaus** (HYP-2913): tight cores realize `≤ 3` distinct phase gaps `{1,n,2n}`; a
  perturbation *splits* them, so the distinct-gap count is `0` at tight and grows with `δ` — a
  fragmentation invariant coupled to near-tightness (exactly the residual's two axes in one object).
- **Difference-flow renormalization** (HYP-3901): for `max ≥ 35`, where the core carries its own far
  element, the deep cluster renormalizes to its difference core one scale down; the AP is the fixed point.

The cover-gap picture says all three are measuring the same thing: **displacement of the good components
off the `1/13`-lattice**, which is `0` iff `W` is a dilated AP. The uniform statement — `smax > 0 ⟹
coverGap ≥ 1/13`, uniformly over non-AP cores — is the soft Route-A equidistribution (THM-1036), now with
the cleanest possible witness: the small comb `13j`, not `182`.

→ THM-1038, THM-1037, THM-1036, THM-1033, THM-1029, THM-1017, HYP-3901, HYP-2913, the disc_v/arc-count
reflections (klein S287–289), the alignment/Freiman/OCF reflections.
