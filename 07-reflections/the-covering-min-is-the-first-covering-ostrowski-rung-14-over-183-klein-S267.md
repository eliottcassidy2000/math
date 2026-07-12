# The covering-min is the first covering rung of the Ostrowski ladder — 14/183, not 1/13

*klein-2026-07-12-S267. Owner directive: keep sharpening the LRC(14) crux via hypothesis
investigation. Three fleet sessions in three days had posted three different "DC floor" numbers —
kps cont.51 said 1/12, boxeph-S20 said 1/13, my own S266 deep well said 14/183. I set out to
reconcile them and pin the actual covering-minimum. They are nested box artifacts, the true answer
is 14/183, and there is a clean one-line reason why.*

---

## The three floors are one ladder, boxed at three scales

The claimed floors, each verified exact, and each the minimum of its own search box:

| claim | value | search box | what it actually is |
|---|---|---|---|
| kps cont.51 (HYP-6175) | 1/12 = 0.0833 | Vmax ≤ 22 exhaustive / ≤32 hill-climb | the **two-block** near-interval stratum floor (`{1,2,3,4,10..18}`) |
| boxeph-S20 (HYP-6150) | 1/13 = 0.0769 | Vmax ≤ 30 exhaustive (10.97M) | the **compressed** subclass floor (`2·{1..12}∪{13}`) |
| klein-S266 deep well (HYP-6165) | 14/183 = 0.0765 | needs Vmax = **182** | the **global covering-min** |

Ordered by structure: `deep-well 14/183 < compressed 1/13 < near-AP 3/37 < two-block 1/12`. Each box
ends just short of the next-deeper extremal, so each session reported its box's floor as "the" floor
— the exact MISTAKE-141 genus, three times over, one layer deeper each time. My hunt (deep-well
variants, base+outlier, 15k structured covering families) found **nothing below 14/183**; the
follow-up CRT-escape search (249 covering families with a speed > 182) found nothing below it either
(min there `28/365`). The covering-min is 14/183, and this independently reconfirms the fleet's own
ILP certificate (HYP-3779: lazy-cut infeasible for all speeds ≤ 182) — the "spread beater below the
deep well" (HYP-3764) is real only for n=7..11 and was retracted by its author at n≥12.

## Why 14/183 — the covering constraint is a rung selector

The clean reason sits in one computation. The **Ostrowski ladder** `{1,…,12, 13k}` realizes the rung
`M = k/(13k+1)` **exactly** (verified k=1..16, every match on the nose):

```
k :   1     2     3    ...   13     14      15
M : 1/14  2/27  3/40  ...  13/170  14/183  15/196
cov:  no    no    no   ...   no    *YES*    no
```

The family is **covering iff `14 | k`** (its far element `13k` must supply a multiple of 14). The
smallest such `k` is **14**, giving `M = 14/183`. That is the whole mechanism:

> The AP `{1..13}` sits at rung `k=1` (`M=1/14`, tight, **non-covering**). The covering constraint —
> "you must carry a multiple of 14" — forbids the low rungs and forces the family up the ladder to
> the **first covering rung, `k=14`**, whose far element is `13·14 = 182` and whose value is
> `14/(13·14+1) = 14/183 = n/Φ₆(n)` (since `Φ₆(14) = 14²−14+1 = 183 = 13·14+1`).

So `14/183 = n/Φ₆(n)` is not a coincidence of the deep well's arithmetic — it is *forced* as the
first rung of the tight ladder that clears the covering demand. The AP and the deep well are the two
named ends of the same Stern–Brocot ray (mac-mini's Ostrowski-ladder reflection); the covering
constraint is exactly the projection off the tight end onto the first covering rung.

## What this sharpens about the crux

**The covering floor is 14/183, and `1/13` is a red herring for the bottleneck.** My own S266 framing
— "compressed DC ⟹ M ≥ 1/13, the hard half" — is true but is a *sub-stratum* statement: `14/183 <
1/13`, and the family that achieves the covering-min (the non-compressed deep well) sits **below**
the compressed floor. The genuine covering bottleneck is `inf over primitive covering 13-sets of M`,
conjectured `= 14/183 > 1/14`, which is **HYP-2566** (= THM-523 part D), and it is **open** — "uniform
looseness IS the covering case of LRC(14), not lighter" (HYP-4067). The compressed `1/13` result
(THM-721, `j≤6`) is a proved island *above* the real floor; the deep well and the large-diameter
incoherent families around it are what live at 14/183.

The honest proof map, after this session:

- **PROVEN:** non-covering ⟹ `M ≥ 1/14` (elementary, THM-366/523). Single-killer shape `{1..12,m}`:
  covering ⟺ `182|m`, `M ≥ 2/27`, min `14/183` (THM-526, gap-free). Covering-min `= 14/183` for all
  speeds ≤ 182 (HYP-3779 ILP). Deep-well loneliness Lean-certified (HYP-4065). Band-edge:
  clear at `14∤q` ⟹ `M ≥ ⌈q/14⌉/q > 1/14` (opus-S235, *conditional on clearing*). Compressed `j≤6`
  ⟹ `M ≥ 1/13 − B/(2L) > 1/14` (THM-721).
- **OPEN (= HYP-2566, the covering case):** the *uniform* `inf_{covering} M > 1/14`. Three linked
  residuals: the **unbounded clearing window** (band-edge margin `→ 1/14` as diameter `→ ∞`,
  HYP-6120); the **CRT-escape tail** (speeds `> 182`; my S267 search extends the empirical certificate
  but proves nothing); the **large-diameter incoherent stratum** (`j≥7`, boxeph-S19's un-shrunk core).

So the single sharpest statement of the LRC(14) crux is now:

> **LRC(14) ⟺ every primitive covering 13-set has `M ≥ 14/183`** (equivalently `> 1/14`), the
> covering-min being the first covering rung `14/183 = n/Φ₆(n)` of the Ostrowski ladder, uniquely
> attained by the deep well `{1,…,12,182}`. Proven for speeds ≤ 182 and on the compressed and
> single-killer strata; the uniform bound over all speeds is HYP-2566, open.

## The shape of it

Four sessions, and the same lesson keeps arriving one layer deeper: a "floor," freshly reported from
an exhaustive box, is the floor *of the box*, and the true extremal is a structured family just past
its edge — the AP at Vmax 13, the compressed near-dilate at 24, the deep well at 182. The fix is never
a bigger box; it is finding the *rung* the constraint selects. Here the covering constraint turned out
to be a rung selector on a ladder that was already in canon, and the covering-min fell out as
`n/Φ₆(n)` with a one-line reason. The number was known; what was missing was that it is *forced*, and
that the `1/13` everyone (me included) had been calling the crux sits strictly above it.

*Files: `04-computation/lrc14_covering_min_hunt_klein_S267.py` (+out), and the ladder / CRT-escape
verifications inline. HYP-6180. Reconciles kps HYP-6175 (1/12), boxeph HYP-6150 (1/13), klein-S266
HYP-6165 (14/183) as nested box artifacts; reconfirms HYP-3778/3779 (covering-min = 14/183, ILP);
sharpens klein-S266 (the covering bottleneck is 14/183 = HYP-2566, not the compressed 1/13). Connects
[[the-covering-min-is-an-ostrowski-ladder-and-the-ap-and-deep-well-are-its-ends]],
[[the-empty-window-needs-the-covering-restriction-GW-and-3-41-klein-S266]], THM-526, THM-721,
HYP-2566/4067 (uniform looseness = the covering case).*
