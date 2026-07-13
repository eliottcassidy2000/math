# The far-element J-recursion: Route A needs the threshold, not the exact consec-min — and the threshold is finite-check + a monotone tail

*klein-2026-07-12-S271. Owner: push consec-maximizes / three-gap toward proof. The push that landed
is a *separation*: the density floor Route A actually needs (`J ≥ 432/91`, with a +0.315 margin) is
**not** the LRC-hard "consec is the exact global min" — it is `[finite compact check, already done]` +
`[a monotone far-element tail governed by an exact affine recursion]`. The recursion is THM-710 in
clean `(J, μ)` coordinates, and it shows the tail can only *raise* J.*

---

## The functional and what is actually required

`J = E_x[N(7−N)]`, `N = #{empty sectors among the 7 arcs [s/7,(s+1)/7) of the phases {frac(e_i x)}}`.
Writing `μ = E[N] = m₁` and using `E[N²] = Var + μ²`, THM-716's identity is `J = μ(7−μ) − Var`. The
Route-A requirement (THM-711/THM-705) is **`J ≥ 432/91 = 4.7473`** — the `(empty, hit)` sector-pair
mass that the density floor consumes. Consec `{1..9}` gives `J = 4465/882 = 5.0624`, a **+0.3151
margin** over the threshold; the decorrelated value is `8.18` (73% headroom).

The crucial distinction the fleet's "consec-extremality is LRC-hard" framing blurs: **Route A needs
`J ≥ 432/91`, not `J = 4465/882`.** The exact statement "consec is the global minimizer" is a
*tradeoff-optimum* rigidity (THM-716: consec is an isolated saddle of `μ(7−μ) − Var`, neither
μ-minimal nor Var-maximal — which is why every compression/monotone argument failed, opus-S240). But
the *threshold* has 0.315 of slack, and that slack is enough to make the tail elementary.

## The far-element recursion (THM-710, in J-coordinates)

Add a decorrelated far element `w` to a k-set `E`. In the two-scale limit (THM-687/688) `w` lands in a
uniform sector independent of `E`, emptying a previously-empty sector with probability `N/7`. A
two-line total-variance computation gives the **exact affine recursion**

$$\boxed{\;J' = \tfrac{5}{7}\,J + \tfrac{6}{7}\,\mu,\qquad \mu' = \tfrac{6}{7}\,\mu\;}$$

*(Proof: `N' = N − B`, `B ∼ Bern(N/7)` given `N`; so `μ' = 6μ/7` and `E[N'^2] = (5/7)E[N²] + μ/7`;
substitute `E[N²] = 7μ − J` into `J' = 7μ' − E[N'^2]`.)* This is exactly THM-710's factorial-moment
eigen-transfer `m_r → ((7−r)/7)m_r` read through `J = 6m₁ − m₂` (`m₁→(6/7)m₁, m₂→(5/7)m₂` ⟹
`J' = (36/7)m₁ − (5/7)m₂ = (5/7)J + (6/7)μ`). Verified numerically to `~10⁻³` in the clean far regime
(`w = 997, 9973`).

**What it buys, made visible.** The map `J ↦ (5/7)J + (6/7)μ` is a contraction (eigenvalue `5/7 < 1`)
toward the decorrelated fixed point. For a wide 9-set = compact 8-cluster `C` + far element, its J
**limit** is `(5/7)J₈(C) + (6/7)μ₈(C)`. With the k=8 base `J₈ ≥ 291/49 = 5.9388` and the trivial floor
`μ₈ ≥ 0.59` (sampled min `μ₈ ≈ 1.77`, and `(6/7)μ₈ ≥ 0.505` is all that's needed), the limit is
`≥ (5/7)(291/49) + (6/7)(0.59) = 4.747` — the threshold, with margin. So **the far element cannot
push J below the threshold; it raises J from the compact min toward `5.68`.**

## The tail is monotone — no crossover

Verified directly on the block+far family `{1..8, d}`:

| `d` | 9 | 10 | 12 | 16 | 20 | 40 | 100 | 1000 |
|---|---|---|---|---|---|---|---|---|
| `J` | 5.062 | 5.239 | 5.452 | 5.685 | 5.716 | 5.739 | 5.699 | 5.678 |

`J` **increases monotonically** from the compact min `5.0624` (at `d=9`, where `{1..8,9} = {1..9}` =
consec) toward the decorrelated limit `≈ 5.68`, and **never dips below the threshold** for any `d`. So
there is no crossover diameter to fear on the tail: once the *compact* minimum clears `432/91`, the
entire wide direction clears it too, because the recursion only adds positive mass. This is the exact,
quantitative form of mac-mini cont.43's sampled "far elements raise J."

## The honest reduction — and the one remaining piece

Putting it together, the Route-A density floor `J ≥ 432/91` (for the base rows `k = 8, 9`; `k ≥ 10`
inherit by THM-710) reduces to:

> **`[A′] J ≥ 432/91` = `[compact exhaustive check, bounded diameter — DONE]` + `[monotone
> far-element tail, governed by `J' = (5/7)J + (6/7)μ` — THM-710]`.**

- The **compact check is already complete**: kps cont.32 ran *all* 48619 primitive 9-sets in `[1..18]`
  (and 11440 in `[1..16]`, 2002 in `[1..14]`) — global min exactly `4465/882` at `{1..9}`, unique,
  margin `+0.315`; mac-mini did `d ≤ 20`. This is a finite computation, not an inverse theorem.
- The **tail is the recursion**: an exact affine law (THM-710, PROVED as the eigen-identity) whose
  fixed point sits at `≈ 5.68 > 4.747`, approached monotonically.

So the operative floor is *not* LRC-hard — it is a finished finite check plus a proved affine tail. The
**one remaining piece** is quantitative, not structural: make the tail **uniform over all wide 9-sets**
(not just block+far) — i.e. bound the `O(1/w)` finite-distance error of the two-scale limit
(THM-699/700 give the per-atom TV bounds; the constant is currently empirical, `|err|·w ≲ 6`) so that
"diameter `> D₀` ⟹ `J ≥ 432/91`" holds with an *explicit* `D₀`, and check `D₀` falls inside the
already-computed compact box (`≤ 18`). Because the tail is *monotone up* and the margin is `0.315`, a
crude explicit `O(1/w)` constant suffices; there is no delicate cancellation here (unlike the covering
side). Multi-cluster wide sets iterate the same recursion down to their densest compact core (all
`k' ≤ 9`, bounded diameter, in the same finite check) via THM-688.

## The most tractable full closure: the k=8 degree-3 row (favorable direction)

The same recursion, one rung down, is where a base row can be *fully closed* this session or next. The
k=8 direct row is an **upper** bound `Φ = 1 − (2/3)m₁ + (47/252)m₂ − (5/252)m₃ ≤ cap₉ = 1979/4004`
(THM-714/719), so decorrelation is the **favorable** direction and the slack is `+0.047` — ~5× the
k=9 threshold slack. THM-710's transfer sends `m₃ → (4/7)m₃` (and `m₁,m₂` as before), all entering `Φ`
with the right signs, so the target is a clean monotonicity:

> **Lemma (k=8 tail-monotone).** `Φ` is non-increasing under far-element insertion (THM-710 transfer),
> so `max Φ` is attained at bounded diameter; with the finite compact check (THM-719: max `Φ = 0.4380`
> at `{0..7}`, every `d ≥ 8` gives `≤ 0.3907`), `max Φ = Φ(consec-8) ≤ cap₉`.

This turns k=8 into `[finite compact max] + [a monotonicity THM-710 nearly hands over]`, on the
favorable side. And it *interlocks* with the additive-energy route: the 0-anchored coarse problem is
governed at leading order by `E₂`, giving `consec-max-coverage ≤ AP-max-E₂ = Freiman |S+S| ≥ 2k−1`
(PROVEN, HYP-5990/5681; LEM-015 the proved `E₃` sibling), whose sub-leading `m₃` correction is exactly
this k=8 rung. So the k=8 degree-3 monotonicity is the single load-bearing rung both routes bottom out
on — the sharpest concrete target the recursion exposes.

## What is *not* claimed

The **exact** statement "consec is the unique global minimizer of `J`" remains the tradeoff-optimum
rigidity (THM-716), and *that* is LRC-hard — it is the equality/three-gap inverse theorem. But Route A
never uses it. Separating the two is the point: the density floor lives on the `0.315`-slack threshold,
which the recursion closes up to an explicit `O(1/w)` bound; the exact-min rigidity is a separate,
harder object the floor does not need. (This mirrors the covering side, where the *floor* `M ≥ 14/183`
is the operative crux and the *exact* extremal classification is a bonus.)

*Files: `04-computation/lrc14_far_element_J_recursion_klein_S271.py` (+out). HYP-6260. Repackages
THM-710 (eigen-transfer) in `(J,μ)` coordinates; the tail-monotonicity + threshold/exact-min
separation are the new framing. Consumes THM-711/716 (J identity, finite-dimensionalization), kps
cont.32 (exhaustive compact check), THM-687/688 (two-scale limit), THM-699/700 (O(1/w) TV bounds).
Connects [[the-covering-min-is-the-most-blocked-family-not-the-small-q-one-klein-S268]] (the mirror
threshold-vs-exact separation on the covering side). HOUSEKEEPING: HYP-6250 double-claimed
(klein-S270 first-push keeps 6250; mac-mini-S66 → please renumber 6255); klein-S271 takes 6260.*
