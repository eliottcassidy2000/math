# The abandoned attempts, decoded: the crux is additive dimension 2 (pinned by 169), not any scalar

*boxeph-2026-07-18-S92. Mining the cluster's history for crux-adjacent dead-ends that were recorded
but not decoded (owner's ask). A single missing structure explains almost all of them, and it connects
two threads (opus-S187, death-star-S56) that were never linked. NOT a proof — a diagnosis that
sharpens the crux and names the one right tool. Verified `lrc_inverse_freiman` + the S92 checks.*

## The recurring failure mode: order-2 scalars are blind to additive DIMENSION

Every *scalar* invariant the cluster built to characterize the tight locus **fails the same way**, and
none decoded why:

| attempt (agent) | invariant | why it failed (recorded) | the undecoded reason |
|---|---|---|---|
| `saw` (mac-mini-S113, HYP-7390) | `Σ_pairs[ρ − 4/169]` (pair-overlap excess) | "measures divisibility, not tightness" | it is **order-2**; divisor-rich sets have high pairwise coherence but high additive dimension |
| additive energy (opus-S181) | `E = #{a+b=c+d}` | "necessary not sufficient; 2-D GAPs carry max E yet loose" | `E` is "only the order-2 shadow" of dimension |
| moments (opus-0630) | `∫m²` | "no moment functional characterizes tightness" | L¹/L² averages wash out the pointwise gradient |
| `deg Q` / function-field (boxeph-S91) | the level | "discreteness erases the tightness gradient" | the value group has no fine structure |

**Decoded (verified this session).** The `saw`-beater `{1,2,3,4,6,8,12,24,36,48,72,144}` (which beat the
tight AP on `saw`, `+0.84` vs `+0.60`) has doubling `|S−S| = 103` — additive **dimension ≈ 4.5** — while
the tight AP has `|S−S| = 23` = the dim-1 minimum. So the beaters are **multiplicatively coherent but
additively high-dimensional.** Order-2 coherence (what `saw`/`E` see) and additive dimension (what
tightness needs) genuinely diverge, and *that* is why every scalar fails. The right invariant is not a
number — it is the **additive dimension**.

## What the crux actually is: additive dimension ≤ 2, pinned by 169 = 13²

death-star-S56 (THM-1028) found the correct statement and half-built it: the deep well's maximizer
**residue set** `R = 14·{1,…,12} ∪ {169}` has `|R−R| = 47 > 3·13−4 = 35`, so ordinary (dim-1) Freiman
`3k−4` **fails on the nose** — `R` is **genuinely dimension 2**. The second progression is the far
residue: `169 ≡ 1 (mod 14)`, i.e. `R ⊆ {14i + j : 1≤i≤12, j∈{0,1}}`, a `12×2` GAP. And the second
coordinate is **arithmetically pinned**: `169 = 13² = q − val = 12·val + 1` (with `q=183=13·14+1`,
`val=14`) — the far residue is forced to `13²`, which forces `val = 14`.

> **Sharpened crux (INV, dimension form).** `M < 1/13` (covering) ⟹ the maximizer residue set has
> additive **dimension ≤ 2**, with the second generator pinned by `169 = 13²`. The `12`-element core is
> the dim-1 slice (the AP); the far element is the entire second dimension.

This is strictly cleaner than "the core is an AP": it is a **dimension bound**, and dimension is what
all the scalar shadows were failing to capture.

## The connection nobody made: opus-S187 ≡ death-star-S56

opus-S187 ("the Freiman `3k−4` finish line splits at burden 12; 2-D GAPs appear at 13") found that
Freiman stability breaks at 13 elements, where **2-D GAPs of unbounded spread** appear (example
`(0,1,7,8,9,15,16)`), and stopped, calling it a "dimension obstruction." death-star-S56 found the deep
well's residue set is a **2-D GAP of 13 elements** with unbounded spread (the far residue). **These are
the same object.** opus-S187's "the burden-13 2-D-GAP obstruction" *is* death-star's "the deep well is
dimension 2." Neither reflection cites the other. Connecting them:

- opus-S187 needed the **sharp 2-D Freiman `3k−4`** (which exists — the burden-13 classification is
  finite: 25 shapes) and left it as the missing step;
- death-star supplied the **pinning** (`169 = 13²`) that collapses opus's "unbounded spread" to a single
  arithmetically-forced value.

Together they are a route: **classify the burden-13 2-D GAPs (opus's 25 shapes, finite), then rule out
all but the pinned deep-well shape via `169 = 13²` (death-star)**. Each half was parked; the composite
was never attempted.

## The two remaining walls are the same wall (union vs correlation)

- **≥6-linear cancellation** (klein-S279, "unproved on either side"): the covering residual is
  `≥6`-linear; the bilinear part is clean, the signal is the higher-order **correlation**.
- **Radius-7 Hamming ceiling** (klein-S315): the tail lemma walls at `r < 1/(2δ) = 6.25`.

Both are **union bounds saturating** — they assume the `r` danger arcs / the `6` sectors are
*independent/disjoint*, and the bound dies when there are too many. **Diagnosis of the radius wall
(verified this session):** the actual `M` of radius-`r` perturbations of `{1..12}` *increases* with `r`
(`r=1: 2/25`, `r=7: ≈1/9`) — rigidity **strengthens** at high radius, so the wall is a pure method
artifact: at `r ≥ 7` the seven thin danger arcs must heavily **overlap** (the true good set is large),
so the disjoint-union bound is wildly pessimistic. The genuinely hard case is `r=1` (`2/25`, already
proved). *Also:* the naive Hamming metric conflates distance with **dilation** — the dilated APs
`c·{1..12}` are "radius 12" yet tight; the right invariant is radius-from-**nearest dilated AP**, i.e.
the **dimension-quotient**, closing the loop with the dimension crux above.

## The one right tool (the parked bridge)

opus-S181 flagged it and no one picked it up: **BSG → Freiman `3k−4` / PFR (Polynomial Freiman–Ruzsa)**
appear in the repo "only as parked sidecars." The dimension crux is exactly a **sharp 2-D Freiman /
PFR** statement. General PFR gives *bounded* dimension with lossy constants; what this needs is the
**sharp `3k−4` in dimension 2** (opus's finite burden-13 classification) plus the **`169`-pinning**
(death-star). That composite — not any scalar functional, not a stronger tail lemma — is the crux.

## Net (honest)

LRC(14) is not closed. But the abandoned attempts, decoded, agree: the tight locus is a **dimension-2
coset-progression** condition pinned by `169 = 13²`, invisible to every order-2 scalar and to every
union/independence bound. The unattempted composite is **opus-S187's burden-13 2-D classification ∘
death-star's `169`-pinning**, and the tool is the **sharp 2-D Freiman**. That is the sharpest form of
the crux the history supports.

Cross-links: [[the-freiman-3k4-finish-line-splits-at-burden-12-2D-gaps-appear-at-13-opus-S187]],
[[the-inverse-theorem-is-a-function-field-freiman-3k-4-the-far-element-is-the-second-dimension-deathstar-S56]],
[[additive-energy-is-necessary-not-sufficient-the-tight-locus-is-resonance-geometry-opus-S181]],
[[both-lrc14-routes-bottom-on-the-same-multilinear-cancellation-not-one-mollification-klein-S279]],
[[THM-1017-ap-core-bridge-reduction]], [[the-lrc14-crux-is-sharp-freiman-additive-energy-and-a-discrete-markoff-spectrum-boxeph-S89]],
HYP-7390, HYP-7362, HYP-7370.
