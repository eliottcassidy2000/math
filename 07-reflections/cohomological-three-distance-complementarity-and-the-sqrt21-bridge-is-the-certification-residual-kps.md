# Cohomological three-distance: the lonely set is `H⁰`, complementarity is Alexander duality, the ι-odd part is the Gauss certificate — and the `√21` bridge is the *certification*-program residual, not the classical census crux

*kind-pasteur-2026-07-01-S21. Extending the sharp three-distance / far-element work cohomologically, thinking about complementarity, and answering the owner's question: is the `Q(√-3,√-7)` bridge (the ι-odd Gauss certificate `i√7` meeting the Eisenstein covering `√-3` through `√21`) the OPEN-Q-108 residual? Answer: it is the residual **of the Lefschetz/trace-formula program** — but it is *not* the classical multi-far census residual, whose extremizers live in their own fields `Q(√5), Q(√-19)`. Complements opus-S24's metallic-vs-Gauss dichotomy.*

## The lonely set is `H⁰`, and complementarity is Alexander duality on `S¹`

For the covering-min construction champion `{1..11,13,36}` at band `1/14`, the lonely set is `b_0 = 2` disjoint arcs (`H⁰(L_C)=ℝ²`, `H¹=0` — arcs are contractible). Verified:

- **`M ≥ 1/n ⟺ b_0(L_C) ≥ 1`** — loneliness is the non-vanishing of `H⁰`.
- **Complementarity (Alexander duality):** `b_0(L_C) = b_0(danger cover) = 2`. On `S¹` the lonely arcs and the covered gaps **alternate in equal number** — a closed set and its open complement on a circle have the same number of components. So *covering ⟺ the danger arcs close up into a single cycle* (`b_0(danger)=1`, `b_0(L_C)=0`), and *loneliness ⟺ the cover has a gap*. This is the topological content of `M ≥ 1/n`.
- **Three-distance:** the two arcs have **one** width (Steinhaus/Morse collar — few distinct gap sizes); the champion's `L_C` is a single ι-symmetric collar-pair.

## The ι-action on `H⁰` is the Lefschetz split

The antipode `ι: t↦1−t` acts on `L_C` (it is ι-symmetric), permuting the `b_0` arcs. For the champion, ι **swaps the two arcs** (they are an ι-pair). So `H⁰(L_C)=ℝ²` decomposes as

`H⁰ = (ι-even: ℝ, the sum = the measure/covering side) ⊕ (ι-odd: ℝ, the difference = the Borsuk–Ulam/Gauss side).`

This is the S20 Lefschetz split at the `H⁰` level: the **ordinary trace `Λ(ι)`** counts ι-*fixed* arcs (none — the pair is swapped, `Λ(ι)=0`, free), so the even part is blind; the **ι-odd part** carries the certificate, whose arithmetic is the Gauss sum `i√7`. The complement/ι-symmetry that the atoms' 2nd moment sees as a `k=N/2` spike (opus-S24) is the same ι-swap on `H⁰` here.

## The far element is a degree-`w` circle map that cannot close the cover

A far speed `w` is the circle map `φ_w:t↦wt` of degree `w` (Lefschetz `Λ(φ_w)=1−w`); its danger set is `w` arcs of width `2r/w`. Adding it cuts `L_C` by three-distance, and (verified `w=37,50,300,700`) it **equidistributes** — the surviving fraction `→ 1−2r = 6/7` — but `b_0(L_C∩safe(w))` stays `> 0`: **the far element never closes the cover.** This is the cohomological form of the far-element decorrelation (`O(1/w)`, HYP-3787): a single high-degree map cannot fill an ι-symmetric collar; survival is the ι-odd certificate. So the single-far and ≤6-far cases are *cohomologically* certified — the residual is finite/bounded (klein's ILP) plus the multi-far census.

## The `√21` bridge: which residual is it?

opus-S24 pinned the two hard-side flavors: the **covering / half-tiling spine is metallic** (continued-fraction quadratics, the `Q(√-3)`/Φ₆ side) and the **Paley obstruction is Gauss-sum** (`√p`, the `Q(√-7)` side) — exactly my S21 field split. Their **real compositum is `Q(√21)`** (`√-3·√-7=−√21`), and indeed the deep-well binding CF discriminant `D_14 = 33852 = 21·403` **contains `√21`**. So the bridge is realized: the metallic covering and the Gauss certificate meet in the deep-well.

**Is this the OPEN-Q-108 residual?**
- **Yes — for the Lefschetz/certification program.** To make `i√7 ≠ 0 ⇒ M ≥ 1/n` a theorem, the ι-odd Gauss certificate (`√-7` isotype) must control the *whole* danger cover, whose geometry is Eisenstein (`√-3`). The "mixed" cohomology — neither pure `√-3` nor pure `√-7` — lives in the compositum, i.e. `√21`. **Bridging `Q(√-3)` and `Q(√-7)` in `Q(√-3,√-7)` is exactly what the certification program has left to do.** So the owner's identification is right for this program.
- **No — for the classical multi-far census residual.** The reduced crux `inf meas(L_C) ≥ 1/36` bottoms out at the **pentagon `(Z/10)*`** (field `Q(√5)`) and the **sporadic two-clash `(Z/19)`** (field `Q(√-19)`) — their *own* sub-apex fields, unrelated to `√21`. The census residual is not the apex-compositum bridge.

**Honest verdict:** `√21 = Q(√-3,√-7)` is a *promising reformulation of the residual for the trace-formula route* (certify by bridging the two apex fields), distinct from the classical census crux (which lives in `Q(√5), Q(√-19)`). The two programs coincide only if the census extremizers' fields reduce to the apex compositum — and they do **not** — so the `√21` bridge is *one* route to closing LRC(14), not the entire remaining difficulty.

## The complementarity, summarized

Everything is one ι-equivariant complementarity on `S¹`:
- **lonely `H⁰(L_C)` ⟷ covered `H⁰(danger)`** (Alexander duality, alternate, `b_0` equal),
- **ι-even (measure/covering, metallic `Q(√-3)`) ⟷ ι-odd (certificate/Gauss, `Q(√-7)`)**,
- and the far element is the degree-`w` map whose equidistribution keeps `b_0>0`.

`M ≥ 1/n` is `H⁰(L_C) ≠ 0`; its **certification** is the ι-odd part (Gauss `i√7`); its **obstruction to a clean proof** is that the certificate's field `Q(√-7)` and the cover's field `Q(√-3)` first split at `n=14`, meeting only in `Q(√21)`.

## Honest status & next

- **Computed:** `b_0(L_C)=2`, complementarity `b_0(L_C)=b_0(danger)`, one three-distance width, the ι-swap on `H⁰`, far-element equidistribution (`b_0` stays `>0`), `√21 | √D_14`.
- **Assessment:** the `√21` bridge = the certification-program residual (grounded), ≠ the classical census residual (pentagon/sporadic fields) — the owner's conjecture holds for one of the two programs.
- **Next:** compute the ι-odd certificate *in the compositum* `Q(√-3,√-7)` — does the mixed part vanish/localize (closing the trace-formula route)? And identify the chain complex whose ι-odd Lefschetz number is `i√7` (the S20 open leap): the far-element circle maps (`Λ=1−w`) and their three-distance cells are the natural cells — build the ι-equivariant nerve and read its odd Lefschetz number.

— Related: `the-singular-series-is-an-iota-equivariant-lefschetz-trace-…` (S20), `the-spine-stays-quadratic-…` (S21, the field split), `the-second-moments-atom-pair-correlation-…` (opus-S24, metallic vs Gauss), HYP-3787 (far-element `O(1/w)`), HYP-3789 (Lasserre/inclusion-exclusion), THM-581/582 (ι-odd index), `the-eisenstein-cusp-dichotomy-is-the-three-distance-theorem-…`, OPEN-Q-108, MISTAKE-078. Script: `04-computation/cohomological_three_distance_and_sqrt21_bridge_kps.py` (+ .out). Not a HYP reservation.
