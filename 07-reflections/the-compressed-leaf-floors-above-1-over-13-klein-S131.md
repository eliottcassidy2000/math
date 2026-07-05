# The compressed target "M ≥ 1/13" needs primitivity — and the non-primitive case closes by scale→sieve

*klein-2026-07-04-S131 (HYP-4093). Owner: keep working toward closure. The LRC(14) formalization's
sole open leaf is `hcomp`: every COMPRESSED covering family is lonely. mac-mini-S45/S46 proposed the
clean tight target `compressed ⟹ M ≥ 1/13` (dilated deep wells attain it). This session found that
target is **false without a primitivity hypothesis**, identified the counterexample, and — the useful
part — showed the non-primitive case closes cleanly by the existing scale→sieve normalization, so
`hcomp` still holds. A scope correction that keeps the proof on track.*

## The finding: a non-primitive compressed covering family below 1/13 (and below the deep well)

`hcomp` (kps-S6 `lrc14_of_compressed`) quantifies over ALL nonzero covering compressed families —
no primitivity. The proposed proof route "`compressed ⟹ M ≥ 1/13`" hits a counterexample:

> **`{2,4,6,…,24,182}` is compressed and covering, with `M = 7/92 ≈ 0.076087`** — below `1/13 =
> 0.076923` AND below the deep well `14/183 = 0.076503`, but above `1/14 = 0.071429` (so still
> lonely). Witness `t = 7/184`. (Verified exactly.)

It is compressed (`182 ≤ 13·24 = 312`), covering (`q=2..14` all divided), and **non-primitive**
(`gcd = 2`). It is `2·{1,…,12,91}`, and `M(2·W) = M(W)` by scaling. So the minimal-tightener
enumerations that reported "covering-min `14/183`" (klein-S128/S129) and "compressed floor `1/13`"
(mac-mini-S46) were implicitly over PRIMITIVE families — the non-primitive dilated class (elements in
the `(13, 26]` band, e.g. `14,16,…,24`) sits outside the `A ⊆ {1,…,13}` decomposition and dips lower.

**Consequences.** `compressed ⟹ M ≥ 1/13` is false as stated; `covering-min = 14/183` is a
*primitive*-only statement. The correct tight target is
> **PRIMITIVE compressed covering ⟹ `M ≥ 1/13`** (tight: `{3,6,…,36,182}` is primitive, `gcd = 1`,
> compressed, covering, `M = 1/13` exactly — mac-mini's dilated deep well, correctly a *primitive*
> extremizer).

## The resolution: the non-primitive case closes by scale→sieve (so hcomp still holds)

The good news — and the point of this note — is that the non-primitive families need no new bound.
Let `v` be compressed covering, `c = gcd(v)`, `w = v/c` (primitive). `w` is compressed (same ratios).
Two cases:
- **`w` covering** ⟹ `w` is a *primitive* compressed covering family ⟹ `M(w) ≥ 1/13` (the corrected
  target) ⟹ `M(v) = M(w) ≥ 1/13 > 1/14`. Lonely.
- **`w` NOT covering** ⟹ some `q ≤ 14` divides no `wᵢ` ⟹ the denominator sieve gives `w` lonely at
  `t = 1/q` with clearance `≥ 1/q ≥ 1/14` ⟹ `v` lonely. (Verified: `{2,…,24,182}/2 = {1,…,12,91}`
  misses `q=14`; `{2,…,22,26,364}/2 = {1,…,11,13,182}` misses `q=12` — both sieve out.)

Either way `v` is lonely. So `hcomp` is TRUE, but its proof must **split on primitivity**: the sharp
`≥ 1/13` bound is needed (and correct) only for PRIMITIVE compressed covering families; every
non-primitive one reduces, via the dispatch's existing scale normalization, to a primitive family
that is *either* a smaller primitive compressed covering case *or* non-covering (sieve). This is
exactly why the naive uniform "`compressed ⟹ M ≥ 1/13`" over-reaches: it tries to prove for the
non-primitive families a bound they don't satisfy — but don't need, because they were already going
to be peeled by the scale→sieve path.

## Why this helps closure

`hcomp` is the sole open leaf. The clean statement to prove is now pinned:
`PRIMITIVE compressed covering ⟹ M ≥ 1/13` (the 12-runner LRC bound, tight at the dilated deep
wells), with the non-primitive families discharged by `gcd`-normalization + the covering
sieve — both already in the corpus (`lonely_exists_of_scale`, `sieve_one_div`). This removes a
false lead (the uniform `≥1/13`) and confirms the target and extremizers are primitive. It also
re-confirms LRC-safety for the missed class: `{2,…,24,182}` and its siblings are `≥ 1/14` (they were
never a threat, but they *are* below the "covering-min," which was the surprising part).

## Honest scope

No new theorem; a scope correction + a constructive resolution. The finding is verified exactly; the
resolution uses only existing corpus lemmas (scale, sieve) + mac-mini's (primitive) `≥1/13` target.
I did NOT prove `primitive compressed ⟹ M ≥ 1/13` (that remains mac-mini's open peel-route target,
now correctly scoped). My earlier same-session "compressed floors at 7/89" pass was the same
minimal-tightener-scope error and is retracted — the correct picture is above.

## Links

- Scripts: `04-computation/lrc14_compressed_extremizers_klein_S131.py`,
  `lrc14_compressed_floor_klein_S131.py` (+ outs). Exact. HYP-4093.
- Corrects / credits: mac-mini-S45/S46 (compressed floor 1/13, dilated deep wells, offset-forcer/
  free-rider — all correct for PRIMITIVE families); kps-S6 `LRCEndgameAssembly` (`hcomp` = sole open
  leaf); klein-S129 (deep well unique — primitive scope). Uses corpus `lonely_exists_of_scale`,
  `sieve_one_div`. Open: prove `primitive compressed covering ⟹ M ≥ 1/13`.
