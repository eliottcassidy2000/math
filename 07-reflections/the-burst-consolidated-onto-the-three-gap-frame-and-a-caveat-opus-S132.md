# The LRC(14) burst, consolidated onto the three-gap frame — prior art + a caveat on the proof path

**opus-2026-07-07-S132.** Owner: *look back through past work; ensure contemporary agents are not
reinventing the wheel; more structural clarifications.* A historical audit (2031 reflections) plus
re-reading the foundations shows the last two weeks (agents opus/kps/mac-mini/klein, "S30–S56") have
been **re-deriving, piecemeal, the synthesis mac-mini-S15 already wrote down** — and points to where
the genuinely-new work sits. This note consolidates the burst onto S15's frame, credits the prior
art, and records a caveat that sharpens the open crux.

## The governing frame (mac-mini-S15, HYP-4412, 2026-06-06) — already established

LRC is an **additive↔multiplicative duality mediated by continued fractions / the three-gap
theorem**:
- **Loneliness `M`** is additive-metric (orbit gaps, arcs, three-distance).
- **Covering / resonance-killing** is multiplicative (`b ∣ v_i`, "force a multiple of 14", dilation).
- The **AP is the unique fixed point** where both optimize; the spectrum is the **Ostrowski ladder**
  `k/(k·13+1)` (for 12 speeds) / `k/(14k−?)`; the gaps are **Farey cells**; the tail is the core
  **dilated** (Steinhaus scale-invariance).

S15 already stated the density floor's role and the proof path verbatim: *"the density floor /
contraction rate IS the quantitative three-gap rigidity: how far below the next rung a non-`{kα}`
family is pushed. The AP is the `g=2` fixed point; any detuning increases `g` and the value jumps to
the next rung."* Proof path: **`M(S) < rung ⇒ g(S) ≤ 3 ⇒ {kα}-orbit ⇒ M ∈ rungs`** (converse-three-gap,
Sós / Świerczkowski / van Ravenstein).

## The prior art the burst rests on (cite these; do not re-derive)

| burst result (S30–S56) | it is a facet/refinement of | earliest |
|---|---|---|
| "saturated reduction" / "counterexample ⟹ saturated" (kps-S55/56, opus-S131) | the **denominator sieve** `sieve_frac` / `counterexample_needs_all_divisors` | **oracle-S18** (2026-05-31) |
| "saturated ⟹ spread ⟹ margin" (opus-S131, kps-S55) | "**to kill resonance `b` you need spread, which raises `M`; the AP is the unique least-spread killer**" | **kps-S31p** (thread 3 of S15) |
| μ_{1/7} AP-minimality + exact constants (opus-S130) | "the density floor **IS** the quantitative three-gap rigidity" | **mac-mini-S15** |
| Farey/Ostrowski ladder, gap `(1/14,2/27)` empty (kps-S54, mac-mini-S26/39) | "spectrum = Ostrowski ladder; gaps = Farey cells" | **mac-mini-S26** / S15 |
| coarse bound `M(v) ≥ M(K) − A/L` (kps-S52/53, mac-mini-S38) | "the tail is the core dilated (Steinhaus self-similarity)" | S15 thread 5 |
| AP conjugate witness (klein-S152) | the AP's roots-of-unity `t = c/(n·d·L)`, `c ∈ (ℤ/14)*` | (new, on the frame) |

**Genuinely new in the burst** (build on these, they aren't re-derivation): opus-S130's **exact
μ_{1/7} constants** (`477/1078` etc.) + the **PZ reduction `μ_{1/7} ≥ E[U]`** — the *quantitative*
content S15 called for; klein-S152's **provable conjugate witness**; kps-S56's **dilation-invariance
fix** (primitive saturated); mac-mini-S39's **single-scale-has-no-escape** (the finiteness the
multi-scale escapes lack). The rest is convergence from different angles onto S15's frame.

## A caveat that sharpens the crux (new, opus-S132)

Extending S15's `g`-vs-`M` table from the 12-speed `(C)` gap to **LRC(14) (13 speeds, threshold
1/14)** (`lrc14_gcount_rigidity_opus_S132`), the clean "near-tight ⟺ small `g`, loose ⟺ large `g`"
**does not survive** at the *optimal* witness:

| family | M | g (at M-witness) |
|---|---|---|
| AP {1..13} (tight) | 1/14 | **1** |
| GW {1..11,13,24} (tight) | 1/14 | **2** |
| {1..12,182} (deep well) | 14/183 | 2 |
| {2..14} (LOOSE, M=1/8) | 1/8 | **2** ← loose but small g |
| primes (LOOSE, M=1/2) | 1/2 | **1** ← loose but g=1 |

So `g` at the optimal witness is small for *structured* families irrespective of `M`. The step S15's
path needs — **`g ≤ 3 ⇒ {kα}-orbit`** — is a *non-classical converse* and is **false in general**:
`{2,…,14}` at `t=1/16` has `g=2` yet is **not** an arithmetic orbit (it is a proper subset of one).
The three-gap theorem is one-directional (`{kα} ⇒ g≤3`); the reverse needs the full **Sós /
van Ravenstein three-distance *structure*** (which subsets of `{jα}` occur with `g≤3`), not the gap
*count*. So routing `(G)`/the moat through "small `g`" alone is not a clean classical reduction — the
missing ingredient is *which* `g≤3` configurations are forced by near-tightness, i.e. the actual
Sós word, not the count. This is exactly where the quantitative density floor (μ_{1/7}, `E[U]`) does
the work the `g`-count cannot: it measures *how much* margin a non-orbit configuration loses.

## Recommendation (to stop the re-derivation)

1. **Re-read before extending:** `…three-gap-quantization…macmini-S15` (the frame),
   `lrc-denominator-sieve-lean-and-two-regimes-s18` (the sieve foundation), and kps-S31p's
   resonance-killing (spread⟹raises-M). Most "new" saturated/spread/rigidity claims are there.
2. **The crux, stated once in the unified language:** *a primitive single-scale non-AP 13-family has
   `M > 1/14`*, i.e. detuning the AP's `{kα}`-orbit witness forces `M` up to the next Farey rung
   `2/27`. The tool is the **quantitative density floor** (μ_{1/7} margin / `E[U]`), NOT the `g`-count
   (which this note shows doesn't separate tight from loose). Multi-scale is handled (coarse + klein);
   far-element is peeled; the residue is this one single-scale rigidity.
3. **Do not** re-formalize the coarse bound (exists 2×: LRCCoarseReduction, LRCDecorrelation) or the
   sieve (LonelyRunner `sieve_frac`, oracle-S18).

## Net

The burst is *convergent refinement*, not waste — but it under-cited S15/S18/S31p and re-derived the
frame several times. The frame is settled; the open crux is the **quantitative** density floor on
single-scale primitive saturated families (the Sós-structure content, not the `g`-count), and the
new bricks (μ_{1/7} constants, `E[U]`, conjugate witness, primitive-saturated fix) are the pieces
that actually advance it.
