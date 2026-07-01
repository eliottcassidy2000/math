# The ι-odd certificate in `Q(√-3,√-7)` localizes to `√21 = g_3·g_7` — the bridge stays quadratic, `√21` is the product of the LRC(6) and LRC(14) apex Gauss sums; plus three external leads

*kind-pasteur-2026-07-01-S22. Computing the next step: the ι-odd Gauss certificate in the compositum `Q(√-3, √-7)`. It localizes cleanly to the real quadratic subfield `Q(√21)`, and `√21` is exhibited explicitly as the product of the covering-apex (p=3) and certification-apex (p=7) Gauss sums. Plus an honest assessment of three external leads the owner flagged (LTC/left-right Cayley complexes, Cornell CS6840 Hotelling games, an AI theorem-proving pipeline).*

## The compositum certificate localizes to `√21` (verified)

The bridge field is the biquadratic `Q(√-3, √-7)`, `Gal = (Z/2)²`, with three quadratic subfields: `Q(√-3)` (the Eisenstein **covering** field, Φ₆, the `p=3`/LRC(6) apex), `Q(√-7)` (the Gauss **certification** field, the `p=7`/LRC(14) apex), and the **real** `Q(√21)`. Computing the ι-odd (Gauss) certificate in the compositum:

> `g_3 = i√3`, `g_7 = i√7`, and **`g_3·g_7 = −√21`** — a *single real number* in `Q(√21)`. Hasse–Davenport: `g_21 = χ_3(7)·χ_7(3)·g_3·g_7 = (−1)·(−√21) = +√21`, the real quadratic Gauss sum for `D=21≡1 mod 4`. All verified numerically.

So **the "mixed" part of the certificate — neither pure covering nor pure certification — localizes to the real quadratic `Q(√21)`.** It is degree 2, *not* degree 4: **bridging the two apex fields costs a field extension but not spectral complexity.** For the trace-formula route this is the good news — the S22 open leap (make `i√7 ⇒ M≥1/n` rigorous by bridging to the covering) does not blow up; the bridged certificate is still one quadratic number.

## `√21` exhibited explicitly, three ways

1. **As a Gauss-sum product:** `√21 = g_3·g_7 = (\text{LRC(6) certificate})·(\text{LRC(14) certificate})`. The bridge certificate for the *open* case is literally the product of the *proved* case's (`p=3`) and the open case's (`p=7`) apex Gauss sums. LRC(6) and LRC(14) — the two triple-pillar cases (S18) — multiply into the LRC(14) bridge.
2. **In the deep-well CF:** `D_14 = n(n−1)(n²−n+4) = 33852 = 21·403·4`, so `√21 ∣ √D_14` — the badly-approximable binding continued fraction (Herman rigidity) carries the compositum in its discriminant.
3. **As the real subfield** of `Q(√-3,√-7)`: the covering (imaginary `√-3`) and certification (imaginary `√-7`) meet in the *real* `Q(√21)` — the certificate is real precisely because it is a product of two imaginary units.

This answers S21's question ("exhibit `√21` in the covering-min certificate"): the covering geometry (`√-3`) and the loneliness certificate (`√-7`) are conjugate imaginary quadratics whose product is the real `√21` that the deep-well continued fraction already contains.

## Three external leads (honest assessment)

**(1) `Good Locally Testable Codes` (Dinur–Evra–Livne–Lubotzky–Mozes, Annals 2026 203-2) — STRONG.** The c³-LTC breakthrough builds codes of constant rate/distance/query from **left-right Cayley complexes** (a 2-dimensional square complex from two generating sets). Two real bridges: **(a) local-testability = local-to-global = reconstruction** — the entire "does `H` close reconstruction / the moment ladder / the n=7 wall" thread (S14–S19) is *exactly* asking how much of the global tournament class a *local* (bounded-invariant) view certifies; an LTC is a structure where local views *do* certify the global. **(b)** The tournament **tiling model / half-tiling is a square complex** (rows × columns of the staircase), and Lubotzky's Ramanujan/expander arithmetic is the **`√p` Gauss-sum** flavor (S23/S20) — the Paley skew spectrum `{0,±i√p}` is a Ramanujan-type bound. *Lead: is the tiling/half-tiling square complex a (left-right) Cayley complex, and does its "local testability" quantify the reconstruction wall?*

**(2) Cornell CS6840 Algorithmic Game Theory, Hotelling games (Tardos, Sept 2025) — MODERATE-STRONG.** Sept 10 was Hotelling games (facility location on a line/circle). **The LRC covering-min *is* an adversarial facility-location game on the circle:** the runners are facilities, every point is "served" by its nearest runner, the **lonely observer is the point farthest from all facilities** (`max_t min_v ‖vt‖`), and the **covering-min is the adversary's (min over runner configs) value** — the attacker-defender framing (my earlier thread). AGT tools — potential functions, price-of-anarchy/stability — are candidate levers for the `inf meas` bound. *Lead: read the covering-min as a Hotelling/facility-location game value and import a potential-function argument.*

**(3) `Pengbinghui/pipeline-math` (AI prover-verifier, Lean 4) — MODERATE (methodological).** An agentic LLM prover-verifier pipeline that has resolved open problems including **"tiling and complement problems"** — the project's exact domain (tournament tilings + complement). *Lead: the methodology (prover-verifier + Lean formalization) is a tool for the finite/MSS residual of LRC(14) and for formalizing the reconstruction/half-tiling results; the specific tiling-complement problems it solved are worth reading for overlap.*

## Honest status & next

- **Verified:** `g_3·g_7 = −√21`, `g_21 = √21` (Hasse–Davenport), the biquadratic subfield structure, `√21 ∣ √D_14`, the localization to degree-2 real `Q(√21)`.
- **Significance:** the bridge certificate is a single real quadratic number — the trace-formula route stays low-complexity across the field split. `√21` = (proved case)×(open case) apex Gauss sums.
- **Leads (tangential, unproven):** LTC/Cayley-complex ↔ reconstruction/tiling-square-complex (strongest); Hotelling game ↔ covering-min game value; pipeline-math ↔ AI/Lean for the finite residual.
- **Next:** (a) test whether the tiling/half-tiling square complex is a left-right Cayley complex and whether "local testability" is the reconstruction threshold; (b) formalize the covering-min as a Hotelling game and look for a potential function; (c) build the ι-equivariant nerve (S22) and check its odd Lefschetz number is the `√21`-controlled bridge.

Sources for the external leads: [Good Locally Testable Codes (Annals 2026 203-2)](https://annals.math.princeton.edu/2026/203-2/p03), [Cornell CS6840 AGT Fall 2025 lectures](https://www.cs.cornell.edu/courses/cs6840/2025fa/lectures/), [Pengbinghui/pipeline-math](https://github.com/Pengbinghui/pipeline-math).

— Related: `cohomological-three-distance-…` (S22, the `√21` bridge as certification residual), `the-spine-stays-quadratic-…` (S21, field split), `the-singular-series-is-an-iota-equivariant-lefschetz-trace-…` (S20, `i√7`), `does-H-close-reconstruction-…` + `the-quarter-tiling-…` (reconstruction/local-to-global), HYP-3815/3816 (opus Lefschetz/metallic), OPEN-Q-108. Script: `04-computation/compositum_certificate_sqrt21_kps.py` (+ .out). Not a HYP reservation.
