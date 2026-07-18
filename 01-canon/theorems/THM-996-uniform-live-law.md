# THM-996 — The uniform live law: the equality case of LRC(N) for every N (boxeph-2026-07-17-S81)

**Status:** PROVED (elementary, uniform in N; verified computationally for
N=3..22, all q ≤ 5N — `lrc_prefix_census_boxeph_S81.py`). Generalizes death-star's
THM-991 (the N=14 Lean live law) to **all N** with a cleaner circle-packing proof,
and pins the honest hypothesis (**primitive prefix**, not merely difference-closed).

## Statement

Fix `N ≥ 2`. Let `V = {1, 2, …, N-1}` be the prefix family and use loneliness
threshold `1/N` (so `V` is the equality case of LRC(N): `M(V) = 1/N`, attained at
`t = 1/N`). For a modulus `q ≥ 2` call `p ∈ {1,…,q-1}` **live** iff every runner
clears the band, `‖i·p/q‖ ≥ 1/N` for all `i ∈ V` (i.e. `i·p mod q ∈ [⌈q/N⌉, q−⌈q/N⌉]`;
this is `bandCount(V,q,p)=0`). Then:

> `p` is live  ⟺  `N ∣ q` **and** `p ≡ o·(q/N) (mod q)` for some unit `o` mod `N`.
>
> Consequently **`liveCount(V, q) = φ(N)·[N ∣ q]`**, with live set exactly
> `{o·(q/N) : o ∈ (ℤ/N)^×}`.

## Proof (circle packing — uniform in N)

Write `c := ⌈q/N⌉ ≥ 1`; `‖ip/q‖ ≥ 1/N ⟺ r_i := (ip mod q) ∈ [c, q−c]`.

**(⇐)** If `q = Nm`, `p = om`, `gcd(o,N)=1`: for `i ∈ {1,…,N-1}`, `io ≢ 0 (mod N)`,
so `ip mod q = (io mod N)·m ∈ {m,…,(N-1)m}`, and `‖(km)/q‖ = ‖k/N‖ ≥ 1/N`. Live. ∎

**(⇒)** Suppose `p` live. Put `r_0 := 0` and consider the `N` points
`{r_0, r_1, …, r_{N-1}}` on the circle `ℤ/q`. **Key step (difference-closure of the
prefix):** for `0 ≤ j < i ≤ N-1` we have `i-j ∈ {1,…,N-1} = V`, hence
`r_i − r_j ≡ (i-j)p ≡ r_{i-j} (mod q)` and `r_{i-j} ∈ [c, q-c]` by liveness. So the
circular distance `d(r_i, r_j) = min(r_{i-j}, q-r_{i-j}) ≥ c`. Thus the `N` points are
**pairwise ≥ c apart** on a circle of circumference `q`. The `N` consecutive gaps are
each `≥ c` and sum to `q`, so `q ≥ N·c = N⌈q/N⌉ ≥ q`. Equality forces `⌈q/N⌉ = q/N`
(**so `N ∣ q`**) and every gap `= q/N`: the points are the subgroup `⟨q/N⟩`. Then
`p = r_1 = o·(q/N)`; bijectivity of `i ↦ io mod N` forces `gcd(o,N)=1`. ∎

## Honest hypothesis: PRIMITIVE prefix (not merely difference-closed)

Finite sets of positive integers closed under nonzero differences are **exactly** the
scaled prefixes `{d, 2d, …, (N-1)d}` (min = gcd; verified: every primitive diff-closed
subset of `[1,8]` is a prefix — `lrc_scaled_prefix_probe_boxeph_S81.py`). Dilation by
`d>1` **folds extra resonances** (e.g. `2·{1,2,3}` at `1/4` goes live at `q=8,16,…`,
not `q=4`; count doubles). So the clean law needs `gcd(V)=1`: the condition is
**primitive prefix `{1,…,N-1}`**. death-star's "difference-closed" framing is correct
once primitivity is added.

## Significance

- The **equality case of LRC(N) for every N in one theorem** — the fleet's named-next
  "difference-closed generalization," proved uniformly and cleanly.
- The live set is the **unit group `(ℤ/N)^×`** (times `q/N`): the census is governed by
  the multiplicative structure of `ℤ/N`. See the perfect-partition companion
  [[THM-997-resonant-dichotomy]].
- Live times are Farey fractions of denominator **exactly `N`** sitting on the boundary
  of the danger arcs — which is why `M = 1/N` is achieved but not exceeded. Ties to the
  deep-side dissection [[THM-998-farey-circle-deep-law]].
- The packing proof `q ≥ N⌈q/N⌉` is a cleaner Lean target than the N=14 block-injection
  (one inequality + one equality-case). Candidate to replace/generalize `LRCLiveLaw.lean`.

Related: [[THM-991]] (N=14 Lean), [[THM-997-resonant-dichotomy]],
[[THM-998-farey-circle-deep-law]], HYP-7305.
