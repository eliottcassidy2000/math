# The density-floor theta-sum is tail-dominated — the Selberg bound must carry the width, and the harmonic route is n-specific

*mac-mini-2026-07-06-S22 (HYP-4512). Owner: work the sole open piece (the density
floor) collaboratively. This note develops opus's theta-sum (HYP-4446) and kps's
harmonic-leading-order route (HYP-4467) by explicit relation-length shells, and
reports two honest findings that sharpen the route: the short shells do NOT lead
the tiling cancellation (it is tail-dominated), and the harmonic-count route is
n-specific (the n=7 tiler tiles with a single harmonic relation). The rigorous
object is a Beurling–Selberg majorant whose band-limit must scale as `N ~ 2k²` —
i.e. it must carry opus-S113's width. Verified:
`lrc_theta_shells_selberg_macmini_S22.out`.*

## The object

opus HYP-4446: `safe(S,β) = Σ_{a ∈ L(S)} ∏ᵢ ĥ(aᵢ)`, `L(S) = {a : Σ aᵢvᵢ = 0}`,
`ĥ(0) = 1−2β`, `ĥ(m) = −sin(2πmβ)/(πm)` (`|ĥ(m)| ~ 1/m`).  The density floor is
`safe(S, 2/25) > 0` for every non-AP covering family.  kps HYP-4467: shell by
relation length; the length-3 **harmonic** relations `(1,−2,1)` carry the largest
weight after the main term, and (kps-S30, GREEN) a family has *all* harmonic
relations iff it is an AP — a route `safe=0 ⇒ harmonic present ⇒ AP`.

## Finding 1 — the cancellation is TAIL-dominated, not leading-order

Computed the theta-sum truncated to sparse short relations (`|aᵢ| ≤ 2`, ≤ 4 nonzero
coords) at `β = 2/25`:

| family | #harmonic | truncated θ | true safe |
|---|---|---|---|
| AP `{1..12}` | 10 | **0.196** | **0** |
| single-lift `{1..11,24}` | 9 | 0.176 | 0 |
| generic `{1..11,23}` | 9 | 0.177 | > 0 |

The main term is `(1−2β)¹² = 0.123`.  The AP tiles (`safe = 0`), yet its truncated
θ is **0.196** — the short shells *add* `+0.073`, they do not cancel.  The full
cancellation to `0` is carried by the **long tail** (contributing `−0.196`).  So
`|ĥ(m)| ~ 1/m` gives the *shortest relations the largest per-term weight*, but
there are combinatorially many long relations, and the **tail dominates the sum**.
The short-shell truncation is therefore *not* a proxy for `safe`; and the θ-value
does not even order correctly (the AP is highest, not lowest).  **The floor cannot
be read off a finite low-order truncation** — the Selberg tail bound is essential,
not a technicality.

## Finding 2 — the harmonic route is n-specific (the n=7 tiler refutes it)

At `n = 7` (`β = 2/13`):

| family | #harmonic | tiles (safe=0)? |
|---|---|---|
| AP `{1..6}` | 4 | yes |
| **n=7 gap member `{1,5,6,11,16,17}`** | **1** | **yes** |
| generic `{1,2,3,4,5,20}` | 3 | no |

The n=7 gap member **tiles with a single harmonic relation** — so
`safe = 0 ⇒ all harmonic ⇒ AP` is **false at n=7**.  It cancels to `0` through its
generalized-AP structure (`{1,6,11,16}` AP + edges `{5,17}`) *and the wide window*
(opus-S113: `1/91` at `k=6`), not through harmonic saturation.  The harmonic
*count* is a clean structural signal (AP is the max at every `n`), but — like every
structural lens — it is **necessary, not sufficient, and n-blind**.  Only the narrow
`k=12` window `1/325` forces harmonic saturation; the route works at `n=13` *because*
of the width, not the structure.

## The rigorous object: a width-carrying Beurling–Selberg majorant

The correct finite bound is a **majorant**, not a truncation.  Take a band-limited
`g⁺ ≥ g_β` with `ĝ⁺` supported on `|m| ≤ N` (the Beurling–Selberg majorant of the
arc, excess `∫g⁺ − 2β = 1/(N+1)`).  Since `g⁺ ≥ g`, `∏(1−g⁺) ≤ ∏(1−g)`, so

> `safe(S,β) ≥ ∫₀¹ ∏ᵢ (1 − g⁺(vᵢ t)) dt = Σ_{|aᵢ| ≤ N} ∏ᵢ ĝ⁺(aᵢ)` — a **finite**
> theta-sum; positivity for a fixed non-AP family is a bounded check.

But the excess `1/(N+1)` costs `~ n/(N+1)` in the bound, so it is only positive when
`safe(S) > n/(N+1)`.  For a non-AP family `safe ~ window = 1/(k(2k+1)) ~ 1/(2k²)`
(opus-S113), forcing

> **`N ≳ 2k²`** — the band-limit must scale with the *square* of the family size.

This is the analytic form of STRUCTURE × WIDTH: the majorant must resolve the arc to
the window scale `1/(2k²)`, i.e. it must **carry the width**.  A fixed-order Selberg
bound (n-uniform) cannot close the floor — matching that no n-blind lens can.  At
`n=13`, `N ~ 288`; the finite theta-sum over `|aᵢ| ≤ 288` relations is the honest
(infeasible-to-enumerate, but *rigorous*) object.  The route's remaining content is
purely: **prove the width-`N` majorant bound is positive for every non-AP covering
12-family** — a single Beurling–Selberg estimate, now with its scaling pinned.

## Net (collaborative)

- opus HYP-4446 (θ-sum) + kps HYP-4467 (harmonic leading order): correct framework,
  but **tail-dominated** — the short shells do not lead the cancellation, and the
  harmonic route is **n-specific** (n=7 tiler, 1 harmonic, tiles).
- The rigorous object is a **Beurling–Selberg majorant with `N ~ 2k²`** — it must
  carry opus-S113's width; a fixed-order bound cannot work (consistent with the
  n-specificity being unbeatable by any n-blind certificate).
- This pins the sole open piece to one classical estimate with an explicit scale,
  and warns against over-investing in the leading-order-only (structural) route.

## Pointers

- `lrc_theta_shells_selberg_macmini_S22.py/.out` (shell computation, both n).
- opus HYP-4446 (θ-sum), HYP-4456 (Freiman / structure × width, Farey wall q≥3k+2),
  kps HYP-4467 (harmonic leading order), HYP-4457 (Paley), S30 (harmonic ⇔ AP);
  mac-mini HYP-4452 (walls / n-specificity), HYP-4482 (n=7 formal), HYP-4432 (q≤2max).

## Appendix (S23): the identity converges, at a rate set by relation density

Validating opus's identity at n=7 (`β=2/13`, all `|aᵢ|≤N` relations vs the exact
arc-measure `safe`): the truncation converges — loose `{1,2,4,8,9,15}` (safe 0.164,
θ(N=3) 0.159) and the n=7 tiler `{1,5,6,11,16,17}` (safe 0, θ(N=3) 0.0046) are both
close — but the **AP `{1..6}`** (safe 0, θ(N=3) **0.053**) and the AP-fragment
`{1,2,3,4,5,20}` (safe 0.023, θ(N=3) 0.058) converge **slowly**. The truncation error
scales with **relation density**: AP-like sets (dense short relations) need large `N`;
few-relation sets converge fast. So the Beurling–Selberg `N` is driven by the *near-AP*
regime — the floor's hard case — confirming `N ~ 2k²` there. Full-enumeration θ over 12
runners is feasible only at `N=1` (`3¹²`); `N=2` is `5¹²`, infeasible. **The route is
rigorous but not computable at n=13** — the remaining work is the *analytic* majorant
estimate. Refinement: the gap-member candidates are generalized APs with *fewer*
relations than the AP, so their θ converges faster — the effective analytic bound for
the actual candidates may be more tractable than the worst-case AP.
(`lrc_theta_convergence_n7_macmini_S23.out`.)
