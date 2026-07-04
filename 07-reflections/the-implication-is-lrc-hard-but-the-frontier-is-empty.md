# The remaining implication is LRC-hard — but the hard family is empty

*kind-pasteur-2026-07-03-S38. Asked to prove `M = 1/14 ⟹ dilated AP`, I could not — it is the
LRC(14) rigidity, which contains the bound and so is not lighter than the conjecture. What I
could do is map exactly where the difficulty lives and show, three ways, that the family it
would live on is empty. A synthesis, not a proof.*

## Why the implication is not a shortcut

`M = 1/14 ⟹ dilated AP` looks like a structural characterization, but it is the equality case
of the LRC(14) bound bundled with the bound itself: to know that a family *off* the AP locus
has `M` **strictly** above `1/14` is to know the bound is met strictly away from a measure-zero
set — the same wall that keeps `n = 14` open. Every route I tried bottomed out there:

- **Runner removal.** `μ_v = μ_{v∖r} − Leb(D_r ∩ safe_{v∖r})`, and each danger set has measure
  `1/7`, so `μ_v ≥ μ_{v∖r} − 1/7`. But this crude bound is useless where it matters: for the
  deep well, `max_r μ_{v∖r} = 0.103 < 1/7`, yet `μ_v = 0.024 > 0`. The danger sets are
  **positively correlated** (they concentrate in the same safe regions), so `Leb(D_r ∩ safe)`
  is far below `1/7 · (mass)` for the small runners — and controlling that correlation is
  exactly opus's resonance term `R`.

## The one clean reduction: dominant-far is measure-independence

The removal identity does close one case cleanly. Remove the **large** runner. Its comb is
fine (`v_r` teeth) and equidistributes, so `Leb(D_r ∩ safe_{v∖r}) ≈ (1/7)·μ_{v∖r}`, giving

    μ_v ≈ (6/7)·μ_{v∖r} > 0.

Verified on the deep well: `μ_v = 0.0239`, `(6/7)·μ({1..12}) = 0.0292`, agreeing to `0.005`.
Since `μ_{v∖r} > 0` whenever the 12 remaining runners have `M > 1/14` — which LRC(13) (proved)
guarantees, `M(12 speeds) ≥ 1/13 > 1/14` — the dominant-far case gives `μ_v > 0`, i.e. LRC, by
LRC(13) plus one equidistribution step. This is the far-peel in measure form, and *sharper*:
its threshold scales with the piece count of `safe_{v∖r}`, not with `V²`. (The rigorous
discrepancy bound is still looser than the true error, so it does not yet reach `182` on the
nose — but the mechanism is the right one, and it connects the far-peel to the measure route.)

## The hard family is empty (three routes, no gap)

The obstruction is only the *all-comparable* case — no dominant runner to remove. But there the
family cannot be tight:

- **Tight ⟹ small-speed** (mac-mini THM-610/612): every `M = 1/14` family sits on a
  `(q*/14)`-dilated 14th-root config, which forces small primitive speeds. The tight locus is at
  least `{AP, GW}` — GW = `{1..11,13,24}` (mac-mini THM-612 refuted my S37 "no GW", which was a
  search artifact) — but both are small-speed and non-covering, so the frontier claim is
  unaffected: *no large tight family*.
- **Large all-comparable ⟹ loose** (this session): minimizing `M` over families with all
  speeds in `[N, kN]` (`N` up to 2000, `k = 2, 3`) bottoms out at `M ≈ 0.25–0.33` — **3.5–4.7×
  the danger radius**. No large compressed family comes near `1/14`.

So a large all-comparable *tight* family does not exist, and every covering family falls into
one of three closed routes:

| regime | route | status |
|---|---|---|
| small-speed tight (window / AP / deep well) | small-`q` census | closed per instance (finite check) |
| dominant far runner | measure-independence: `μ_v ≈ (6/7)μ_{v∖large}` + LRC(13) | mechanism clean; discrepancy bound loose |
| large all-comparable | looseness: `M ≈ 0.25 ⟹ μ` large | uniform bound `M ≥ c > 1/14` unproved (HYP-2566) |

## Honest placement

Not a proof of the implication. What this adds to the crux: the runner-removal identity and the
measure-form of dominant-far (`μ_v ≈ (6/7)μ_{v∖large}`, cleaner and sharper than the far-peel),
and the confirmation that the family the implication would have to handle — large,
all-comparable, tight — **is empty**. So LRC(14) for covering families is not one monolithic
inequality; it is a three-way partition with *no missing case*, and the residual rigor is two
computational-looking bounds: the small-`q` census completeness for small-speed families, and
the uniform looseness `M ≥ c > 1/14` for large compressed ones. Both are confirmed
computationally and neither is proved. That is the honest shape of the wall.

---
*Linked: [[the-tight-locus-is-the-arithmetic-progression]] (S37), [[the-covering-min-and-the-gcd-refinement]]
(S36). mac-mini THM-610 (deep-hiding), HYP-2566 (uniform looseness), HYP-4058/opus (measure
form `R`), HYP-4044/kps (far-peel — this is its measure twin). Scripts:
`lrc14_removal_measure`, `lrc14_compressed_frontier_kps_S38.py`. HYP-4064.*
