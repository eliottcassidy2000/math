# The tight-locus rigidity of LRC(13): the free-rider's uncited dependency

*kind-pasteur-2026-07-05-S1, HYP-4096. A reflection on the one lemma the compressed
peel silently assumes, why it is true, and how far it can be proved elementarily.*

## The claim and why it matters

The LRC(14) proof has funnelled to a single open Lean leaf `hcomp` (compressed covering
⟹ lonely), which klein-S132's `hcomp_of_primitive` reduces to `hprim` (PRIMITIVE
compressed covering ⟹ lonely). The peel route splits `hprim` by the peeled base
`W = V∖{argmax}`:

* **tight base** (`M(W) = max_t min_{w∈W} ‖wt‖ = 1/13`): klein's `tight_ap_free_rider`
  closes it — but only for `W = c·{1,…,12}` (the dilated AP);
* **loose base** (`M(W) > 1/13`): the window peel / mac-mini's alignment bands.

The dispatch "base not a dilated AP ⟹ loose" is the **tight-locus rigidity**:

> **RIGIDITY.** For `W` a 12-element set of positive integers with `gcd(W)=1`,
> `M(W) = 1/13`  ⟺  `W = {1,2,…,12}`.
> Contrapositive: `W` primitive, `W ≠ {1..12}` ⟹ `M(W) ≥ 2/25` (the second value,
> attained by `{1,…,11,24} = 2/(2·13−1)`).

It is the one place the free-rider route leans on an **uncited** fact. The
Sungkawichai–Trakulthongchai LRC(≤13) paper (arXiv:2604.23906) does *not* characterize
extremizers (it only cites Goddyn–Wong for the word "tight"), so this is genuine open
mathematics, not a literature citation. mac-mini-S47 verified it computationally; klein-S132
flagged it MISTAKE-100-risk-class.

## Verified decisively (kps-S1)

Exact max-min via critical times `t=m/(w_i±w_j)`:
* **Exhaustive**: all 1820 primitive 12-subsets of `[1,16]` — only the AP is tight.
* **Residue-system adversaries**: all 531 376 primitive sets whose residues mod 13 are a
  permutation of `{1,…,12}` with lifts in `{r, r+13, r+26}` — the *only* sets that can
  even be optimal at `t=1/13` — only the AP is tight.
* **Height**: 400 000 random primitive sets with entries `≤ 200` — zero non-AP tight.
* Second value confirmed `2/25`.

The rigidity is as solid as computation makes it. What follows is how much of it is
*provable* elementarily — and why the full statement is genuinely hard.

## Two elementary necessary conditions (proved)

Write `W = {w₁ < w₂ < ⋯ < w₁₂}`.

**(N1) The reciprocal-window bound: `M(W) ≥ w₁/(w₁+w₁₂)`, so tight ⟹ `w₁₂ ≥ 12·w₁`.**
At `t* = 1/(w₁+w₁₂)` every `w_i t* = w_i/(w₁+w₁₂) ∈ [w₁/(w₁+w₁₂), w₁₂/(w₁+w₁₂)] ⊂ (0,1)`,
a sub-interval of `(0,1)` symmetric about `1/2` (since `w₁₂/(w₁+w₁₂) = 1 − w₁/(w₁+w₁₂)`).
Hence `‖w_i t*‖ ≥ w₁/(w₁+w₁₂)` for **all** `i`, with equality at `i=1` and `i=12`. So
`M(W) ≥ w₁/(w₁+w₁₂)`; tightness forces `w₁/(w₁+w₁₂) ≤ 1/13`, i.e. `w₁₂ ≥ 12 w₁`. This is
an **equality** for the AP (`12 = 12·1`) — the AP is extremal for (N1). Formalized as
`band_margin_reciprocal` / `spread12_lonely13` in `LRCBandMargin.lean` (kernel-pure); it
is the explicit-margin refinement of klein's `spread13_lonely`, and its contrapositive
says every *narrow* band (`w₁₂ < 12 w₁`) is loose with the explicit witness `1/(w₁+w₁₂)`.

**(N2) The gap bound: tight ⟹ `w₁₂ ≤ 12·w₁₁`.**
The 11-subset `W' = {w₁,…,w₁₁}` has `M(W') ≥ 1/12` by LRC(12), so there is `t₀` with
`min_{i≤11}‖w_i t₀‖ ≥ 1/12`. The min is `w₁₁`-Lipschitz, so it stays `> 1/13` on an
interval of half-width `(1/12−1/13)/w₁₁ = 1/(156 w₁₁)` — a gap `G` of length
`≥ 1/(78 w₁₁)` where all 11 small runners are safe. If `W` is tight, `w₁₂`'s danger arcs
must cover `G` entirely; but `w₁₂`'s danger arcs have width `2/(13 w₁₂)` separated by
safe gaps `11/(13 w₁₂) > 0`, so an interval `⊆ danger(w₁₂)` fits inside **one** arc:
`|G| ≤ 2/(13 w₁₂)`. Combining `1/(78 w₁₁) ≤ 2/(13 w₁₂)` gives `w₁₂ ≤ 12 w₁₁`.

Together (N1)+(N2) pin the top of the spectrum: `12 w₁ ≤ w₁₂ ≤ 12 w₁₁`. Both are tight-ish
for the AP. They are real rigidity content — but they bound only the top two speeds.

## Why the full rigidity is hard (and why d=12 is special)

The obstruction is **height**, not spread. (N1)/(N2) do not chain down the spectrum:
there is no elementary `w_{k+1} ≤ C·w_k` for the middle. A priori the whole configuration
could be large. Computationally it never is (only `{1..12}`), so a *height bound*
`w₁₂ ≤ B(explicit)` would upgrade the exhaustive check to a proof — but proving one needs
the decoupling estimate `M(W) → M(W') ≥ 1/12` as `w₁₂ → ∞` made **quantitative and
uniform over the 11-subset**, which is exactly the OPEN-Q-108-flavoured measure control
the whole project is stuck on. So the rigidity is not a shortcut around the hard core; it
is another face of it.

Why does d=12 have *no* sporadic tight family while d=13 has the Goddyn–Wong
`{1,…,11,13,24}`? The GW construction `{1,…,n−2, n, 2(n−1)}` is tight only for
`n` in a residue class mod 6; the 12-speed slot lands off it. So at d=12 the AP is
genuinely alone, which is *why* mac-mini-S47's "tight locus = dilated AP" holds here and
`tight_ap_free_rider` suffices — but the *proof* that no sporadic exists is the missing
piece.

## Status and handoff

* **Proved + formalized:** (N1) as `band_margin_reciprocal` + `spread12_lonely13`
  (`LRCBandMargin.lean`, kernel-pure). Narrow bands are loose, explicit witness.
* **Proved (reflection):** (N2) via LRC(12) + interval covering. Formalizing needs the
  gap/covering machinery (klein's `good_point_in_long_interval` is the seed).
* **Verified, open:** the full rigidity (height bound is the gap). Not a citation.

For klein/mac-mini: use the contrapositive form freely — `W` primitive `≠ {1..12}` ⟹
`M(W) ≥ 2/25`, so `β = 2/25` is available in `lonely_of_window_margin` for every non-AP
base (threshold drops `13B → (25/3)B`). The residual sub-threshold killers remain the
alignment-band job.

See [[lrc14-thread]]; HYP-4096; `04-computation/lrc13_tight_locus_rigidity_kps_S1.py`,
`lrc13_rigidity_hardened_kps_S1.py`.
