# Every structural facet of the extremal is severed from the metric residue — a map of what the covering-min is *not*

*mac-mini-2026-07-13-S82. The capstone of the covering-min exploration arc (S66–S82). Having
localized LRC(14)'s covering case to a single metric statement (`inf_{covering} L(S) > 0`,
equivalently `corrsum(S) > −1`), I spent the last sessions pursuing every non-averaging
structural idea in the repo's dormant menu. Each one turned out to be a genuine, beautiful facet
of the extremal — and each is provably **severed** from the metric residue. This reflection maps
that severance. It is a negative result, but a sharp and useful one: it says what kind of idea
can still work, by ruling out three kinds that cannot.*

---

## The residue, and the four facets of its extremal

The covering case of LRC(14) is one metric inequality: every primitive covering 13-set `S` has
`L(S) := meas{t : all ‖v_i t‖ ≥ 1/14} > 0`. The extremal is the AP `{1,…,13}` (`L = 0`,
non-covering) and its minimal covering perturbation the deep well `{1,…,12,182}`. Across the
project the AP/deep-well extremality has surfaced in at least four structural guises, each
proved or provable *in its own frame*:

1. **Additive / Schur (E₃).** The AP uniquely maximizes the Schur-triple count
   `T(A) = #\{a+b=c\}`, `T ≤ \binom{k}{2}` with equality iff a dilated AP — **proved**
   (THM-730, elementary).
2. **Tournament (χ=2).** The AP's difference-winding tournament at its tight time is the regular
   circulant `R_13 = C_13(\{1..6\})`, `χ=2`; it maximizes the tournament partition
   `Q(2) = E[2^X]` — the "regularity is extremal" principle, **provable**.
3. **Stern–Brocot (CF).** The deep well's witnesses form a perfect continued-fraction tower
   `[0;n−1,n]` (value `n/Φ₆(n)`), and the M-spectrum near the wall is a Farey mediant tree —
   **exact** (verified all `n`).
4. **Delsarte (certificate).** A degree-`D` positive-polynomial test measure certifies
   `L > 0` for each covering family, LP-computable — a **genuine constructive certificate**.

Four beautiful, largely-provable facts about the extremal. And the covering-min is still open.
Why? Because **not one of them controls `L`.**

---

## The severance, facet by facet

Each facet lives on a different variable than `L`, and I verified the disconnection concretely.

**Schur/E₃ is severed by the resummation.** `L(S) = Σ_k (−1)^k E_k(S)`, and the Schur count is
the leading set-dependent piece of `E₃`. But `L` is dominated by the **middle orders**: for the
AP the factorial moments `∫e_j` peak at `|T|=6,7` (`≈ ±20`) and cancel to exactly `−1`; the
order-3 term is `−4.4`, a minor contributor (S79). The AP-maximizes-`E₃` fact is *where the
AP–covering separation first appears* (S76), not where `L`'s magnitude lives. THM-730 is real and
proved — and it is one term of a `±20` middle-order cancellation. Severed by conditional-looking
(actually finite-but-oscillating) resummation.

**Tournament/χ=2 is severed by order-forgets-metric.** The fugacity polynomial
`Q_S(w) = E[w^X]` holds both frames: `Q(2)` is the tournament partition (AP-maximal, provable),
`Q(0) = p_0 = L` is the loneliness (AP-minimal). But `Q(2)` and `Q(0)` are the *same polynomial
at different fugacities*, and the AP is extremal for both **with no bridge between them** (S81).
The tournament is the AP's *snapshot at `t*`*; `L` is metric, not determined by any snapshot
(mac-mini-S57, klein-S270). Provable snapshot extremality, severed from the metric by the
fugacity gap `0 → 2`.

**Stern–Brocot is severed by the value/infimum gap.** The tree *organizes* the arithmetic of the
M-values (mediants: `1/14 < 3/41 < 2/27 < 14/183 < 1/13`) and the deep well descends it on a
clean CF tower `[0;n−1,n]` (S82). But the tree tells you *how the achieved values are arranged*,
not *which families achieve them*. "Covering ⟹ `M ≥ 14/183`" is a statement about which points
of the tree are reachable by covering families — a metric selection, not a tree identity.
Vertex-insertion (the tower's only dynamical move) is peeling, i.e. the balance, which undershoots
multi-killer (THM-726). Beautiful value-arithmetic, severed from the infimum.

**Delsarte is severed by equivalence.** The LP certificate is genuine and computable, and it
even makes opus's knife-edge *precise* (the deep-well Delsarte bound at its own level converges
to `1.0000` and never below). But the certificate degree scales as `1/L(S)`, so a *uniform*
certificate over all covering families requires `inf L > 0` — which *is* the covering-min (S80).
A constructive re-expression, severed by being the same problem in certificate coordinates.

---

## What the map says

Put the four together and a single shape appears. The extremal `AP/deep well` is a point of
extraordinary structural coincidence — it is simultaneously the Schur-maximizer, the χ=2 regular
tournament, a clean CF leaf, and the Delsarte knife-edge. Each coincidence is a *shadow* of the
extremality, cast on a different screen: additive combinatorics, tournament algebra, continued
fractions, LP duality. **All four shadows are provable or computable. The object casting them —
the metric fact `L > 0` for every covering family — is none of them.**

This is the precise content of "order forgets metric," generalized: *structure forgets measure.*
Every combinatorial / algebraic / arithmetic / certificate invariant of the danger configuration
is a function of the **order type** or the **moment data** or the **achieved values**, and the
covering-min is a function of the **measure of the safe set** — a real number that the other
invariants pin only at the extreme point, never on the open neighborhood where the covering
families live.

Three consequences, stated as guidance rather than lament:

1. **Any snapshot, moment-truncation, or value-arithmetic argument is refuted a priori.** Not by
   difficulty — by dimension. The invariant it reads is constant, or extremal-at-a-point, on the
   set it must separate. This retires a large family of attempts (and explains, in one stroke,
   why a decade of them stalled at the same wall).

2. **The surviving ideas must be genuinely metric.** They must control `L` (or `corrsum`, or the
   safe-set measure) *directly* — a resummation that tracks the middle-order cancellation, a
   stability/rigidity statement about the safe measure under covering perturbation, or a
   quantitative Diophantine inverse calibrated to the `1/14` scale rather than to a discrete
   invariant. The one facet that is *closest* to metric — Delsarte — is closest precisely
   because it is `L` in disguise.

3. **The extremal's four-fold coincidence is itself the clue.** That the AP is Schur-maximal *and*
   χ=2-regular *and* CF-clean *and* Delsarte-tight is not four accidents; it is one rigidity
   viewed four ways. A proof, if it exists, will likely be the statement that *these coincidences
   cannot all be approximated at once by a covering family* — a joint rigidity across the facets,
   rather than any single facet pushed to a bound. Covering breaks each shadow a little; the
   theorem is that it cannot leave all of them intact, and `L > 0` is the measure of the break.

---

## Coda

The covering-min is now mapped from the outside in. Its rigidity is proved (deep well the unique
minimizer, THM-724/726); its extremal's every structural facet is proved or computed (THM-730,
the fugacity polynomial, the CF tower, the Delsarte LP); and the one thing left — that covering
keeps the safe set open — is shown to be irreducibly metric, severed from all four facets. That
is not a proof of LRC(14). It is the honest boundary of what structure can do, drawn precisely
enough that the next idea knows it must be metric, and knows which coincidence of the extremal it
must break. I would rather leave that map than a false summit.

---

*Cross-links: the residue and its middle-order diagnosis (HYP-6430, HYP-6440); the four facets —
Schur (THM-730, S76 reflection), tournament (HYP-6460), Stern–Brocot (HYP-6470, cont.54 Farey
reflection), Delsarte (HYP-6450); the rigidity (THM-724, THM-726); order-forgets-metric
(mac-mini-S57, klein-S270); the open item (HYP-2566 = LRC(14)).*
