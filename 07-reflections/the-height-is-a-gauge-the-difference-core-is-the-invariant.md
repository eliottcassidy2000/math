# The height is a gauge; the difference-core is the invariant

**opus-2026-07-06-S110** (HYP-4426). The sole open piece of LRC(14) is the *height bound*:
a single-cluster gap-member (M ∈ (1/13, 2/25), covering, pinned) has bounded height, which —
via the witness-denominator lever (S109: M = c/q ⟹ q ∣ (|vᵢ|+|vⱼ|) ⟹ q ≤ 2·max) — finitizes
(G). Asked to find *indirect* routes to the bound, the productive move is to notice that
**height is not a real degree of freedom.** It is a gauge coordinate. The invariant content of
a family lives in its ratios — equivalently, in its height-independent difference-core.

## Height is a gauge (now formal)

`M(v) = ⨆_t margin v t` is dilation-invariant: `M(c·v) = M(v)` for `c ≠ 0`. Scaling every
speed by `c` merely reparametrizes time by `c`, and `t ↦ c·t` is a surjection of `ℝ`, so the
supremum is unchanged. This is `iSup_margin_const_mul` (LRCDilationInvariance.lean, GREEN,
standard trio) — the formal version of the "up to dilation" the project has always used
informally. Its consequence for (G): the overall *scale* of a gap-member is pure gauge; only
the ratio pattern carries loneliness.

## The difference-core is height-independent (verified)

The consecutive block `{c, c+1, …, c+11}` has difference-core `{1, 2, …, 11}` **at every
height `c`** (verified c = 1, 5, 50, 500). For any single cluster `vᵢ = rᵢ + 13·kᵣ`, the
difference-core scale is the *internal spread* (the range of the `kᵣ`), **not** the height
`13·min kᵣ`. mac-mini's two-scale decorrelation ties the near-carrier loneliness of the
cluster to the loneliness of its difference-core — a height-*blind* object. So the
"height → ∞" limit of the cluster's loneliness is governed by a fixed, bounded structure.
The height simply does not appear.

This is the reframe: **(G)'s open piece is not a height bound — it is a *spread* bound.**

## Two indirect routes the reframe opens

**Route A — the 13-adic spread descent (well-founded).** Write `vᵢ = rᵢ + 13·kᵣ`. The `kᵣ`
are themselves a 12-family; if their spread is large they recurse (their own difference-core),
a 13-adic tower in which each level's "height" is the previous level's spread. The descent
strictly contracts (or exposes a sub-scale-gap, killed by decorrelation), terminating at a
bounded base. The AP — all `kᵣ` equal, difference-core `{1,…,11}` — is the depth-0 fixed
point. The base case is a bounded-height covering family: a **finite check** (S109 lever:
`q ≤ 2·max`). Verified: random `kᵣ` bottom out at descent depth ≤ 1 (well-founded).

**Route B — compactness / the rational–irrational dichotomy.** Suppose there is *no* height
bound: gap-members at heights `H_k → ∞`. By dilation invariance (now formal) normalize
`wᵢ = vᵢ / H_k ∈ [1, ρ]` (single cluster ⟹ bounded ratio `ρ`). Bolzano–Weierstrass gives a
subsequence `w^(k) → w*`, a *real* bounded-ratio family with `M(w*) ∈ [1/13, 2/25]` (M
continuous). But the AP — the *only* real extremal at 12 runners — has ratio 12, so
`w* ≠ AP` whenever `ρ < 12`. The integer gap-members converge to a non-AP real family in the
window. Integrality is lost in the limit: `wᵢ` become irrational. The obstruction the argument
*needs* is an integrality/rationality witness that **survives the limit** — and the S109 lever
supplies exactly one: `M(S_k) = c_k / q_k` with `q_k ∣ (pair-sum)`. The open hinge is whether
`q_k` stays bounded (⟹ rational limit ⟹ finite) or `q_k → ∞` (⟹ irrational limit). Honest
status: Route B does not yet close, but it **localizes** the entire difficulty to the
behaviour of the witness denominator along a height-blowup sequence — a single scalar to
control, not a family.

## The window is a forbidden band (finite-residual confirmation)

The descent base was searched directly: 16,511 covering 12-families across three strata —
lifted-AP, random single-cluster, and pinned (residues a permutation of `1..12 mod 13`) — with
`M` computed *exactly* (grid `b ≤ 2·max`, which the lever certifies captures the maximizer).
**Zero** land in the open window `(1/13, 2/25)`. The near-misses land *exactly at* `M = 2/25`,
never inside. Lifts jump *over* the window: tight AP (`1/13`) → loose (`≥ 2/25`), never into
it. This is the discrete Farey-ladder gap — `2/25 = mediant(1/13, 1/12)` has no integer
realizer strictly between the AP floor and the loose regime. (G) holds on the descent base,
computationally.

## The unifying picture

Height is gauge; the invariant is the ratio pattern (the difference-core). The AP is the
unique fixed point of the ratio-descent, and it is the unique real extremal. The window is
forbidden not because integer families *avoid* a region of scale, but because the achievable
loneliness values are Farey rungs and the open interval `(1/13, 2/25)` sits *between* rungs —
unreachable except by the AP's exact `1/13`. So the density floor (the remaining analytic
input) is asking a gauge-fixed question: among bounded-ratio patterns, the AP alone touches
`1/13`, and everything else has already jumped the band. The witness denominator `q` is the
one scalar that measures "which rung" — bounding it (Route B) or descending the spread tower
to its base (Route A) are the two faces of the same reduction, and both now rest on formal
keystones: dilation invariance (S110) and the witness lever (S109).
