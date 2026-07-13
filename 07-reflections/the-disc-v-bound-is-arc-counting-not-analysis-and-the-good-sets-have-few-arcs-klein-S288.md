# The disc_v bound is arc-counting, not analysis — and the good sets have few arcs

*klein-2026-07-13-S288. Owner: prove the analytic disc_v bound. I proved it —
`disc_v ≤ r²/(3v²)` — with nothing but the triangle inequality on the endpoint sum. The surprise was
not the proof but its consequence: it turns the covering crux from harmonic analysis into **arc-counting**,
and it CERTIFIES the two covering-min extremals — the exact families every structural and elementary
method has failed on for eighty sessions. The reason it works is a fact I had not appreciated: the good
sets have very few arcs (r ≤ 12), because they are heavily-overlapped small-measure unions.*

---

## The bound, and why it is trivial

`disc_v = Σ_{m≠0}|U(mv)|²/(2πmv)²`, `U(ℓ)=Σ_p ε_p e(−ℓp)` the good-set endpoint sum over its `2r`
endpoints. The **trivial** bound `|U(ℓ)| ≤ 2r` and `Σ_{m≠0}m^{-2}=π²/3` give `disc_v ≤ r²/(3v²)`. That is
the whole proof. Fed into THM-731's peeling identity it yields a fully explicit rigorous lower bound
$$L \ge \tfrac1{7}\big(6\,|G'_{~v}| - \sqrt2\,r/v\big),\qquad
L>0 \iff r < 3\sqrt2\,v\,|G'_{~v}|.$$
The harmonic analysis mac-mini-S83 said we needed is **discharged by a triangle inequality** — provided
the arithmetic (arc count vs peel size) is favourable. The remaining question is purely combinatorial:
how many arcs does the leave-one-out good set have?

## The fact that makes it work: few arcs

I had assumed `r` could be `~Σ_{w≠v}w ≈ 78` (the worst-case component count of an intersection of arc
families). It is not: the deep-well base `{1..12}` gives `r=12`, the residue base `{1..11,13}` gives
`r=4`, `{2..13}` gives `r=2`. The good set is a *small-measure* set (`|G'|≈0.03`) cut out by 12
overlapping constraints; the constraints overlap so heavily that only a handful of arcs survive. Small
`r` is exactly what the certificate needs, and it is a structural feature of covering good sets, not an
accident of these examples.

## What it certifies — the binding families

Peeling the far element (`disc ~ 1/v²`), the explicit bound certifies `L>0` for:
- the **deep well** `{1..12,182}` — the *proven global M-minimum* (THM-724/726) — with `L ≥ +0.016`;
- the **near-AP residue** `{1..11,13,84}` — the *min-L `|core|=1` body* (kps cont.70) — with `L ≥ +0.0008`.

These are the families that defeated every structural surrogate (mac-mini-S80–83), the cluster/Mayer
expansion (opus-S269), and every elementary bound. A triangle inequality on the endpoint sum, plus the
smallness of `r`, plus the largeness of the far element, closes them. That the *hardest* family (min-L
residue) certifies while some *easy* families do not (below) is the tell: the crude bound's quality tracks
the **far element size**, not the loneliness — and the extremals happen to have large far elements (182,
84) forced by the covering/divisibility structure.

## Where the triangle inequality is too crude

For covering families with a *small* far element and moderate measure — `{1,3,4,…,14}`, far element 14,
true `L=0.030` — the crude `r²/(3v²)` exceeds the threshold at *every* peel, even though the family is
easy and the **true** `disc_v` certifies it (`+0.018`). Here `|U(mv)|` is genuinely `≪ 2r` (endpoint
cancellation), and discarding that cancellation is fatal when `v` is small. This is the *same* endpoint
cancellation the density route needs for `Q_s` (THM-729) — the shared open core. So the covering case
now splits cleanly:
- **large far element** (incl. all the extremals): closed by `r²/(3v²)` + the combinatorial `r < 3√2 v|G'|`;
- **small far element** (easy, large `L`): needs the shared `|U(mv)|≪2r` cancellation, or a compactness
  argument for AP-like covering sets.

## The honest ledger

The analytic `disc_v` bound is **proved** (`r²/(3v²)`). It gives a rigorous explicit `L`-certificate that
**discharges the covering-min extremals** — real closure of the hardest cases, by elementary means. It
does **not** yet close every covering family: the small-far-element easy families reduce to the shared
endpoint cancellation (density `Q_s`). So eighty sessions of "the residue is irreducibly metric" meets a
one-line triangle inequality that handles the residue — because the metric object (`disc_v`) was, once
built (THM-731), governed by a small integer (the arc count). The metric door was not just open; for the
extremals it was a short step through.

*Files: `04-computation/lrc14_disc_v_bound_klein_S288.py`, `lrc14_disc_v_census_klein_S288.py`,
`lrc14_disc_v_failfamily_klein_S288.out` (+.out). THM-732, HYP-6495. Upgrades
[[the-covering-middle-order-x-integral-is-a-good-set-autocorrelation-discrepancy-and-it-certifies-klein-S287]];
residual shared with [[the-density-weyl-bound-IS-a-relation-lattice-coset-sum-literally-covering-plus-thm538-klein-S285]] (THM-729).*
