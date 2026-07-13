---
source: opus-2026-07-11-S252
status: RESULT + relocation (working the Chebyshev-equioscillation route toward "tight => {k*alpha}"). The
  tight-side {k*alpha} is ELEMENTARY: for any family with NO multiple of 14, M=1/14 forces t*=1/14 to attain
  it, so the phases lie in {j/14} = {k*(1/14)} outright -- no equioscillation needed. Hence "tight => {k*alpha}"
  reduces to "tight => no multiple of 14", which EQUALS the clean-ruler residual ("mult-of-14 => loose",
  verified 11995/11995). So on the tight side the Ostrowski rigidity is not a separate problem -- it IS the
  residual. The equioscillation tool (mac-mini S40, 2-point binding pair) and the genuinely hard {k*alpha} both
  live on the COVERING side (deep well, alpha=14/183), which needs a dual certificate.
tags:
  - lrc14
  - k-alpha-structure
  - chebyshev-equioscillation
  - elementary-lemma
  - no-multiple-of-14
  - residual-unification
  - covering-min
---

# "tight ⟹ {kα}" is elementary and equals the residual; the hard {kα} is the covering side

**opus-2026-07-11-S252.** Owner: work the Chebyshev-equioscillation route toward "tight ⟹ {kα}". Worked it —
and the tight side turns out **elementary**, not needing equioscillation at all. Honest clarification: this
dissolves "tight ⟹ {kα}" as a separate open problem and relocates the real difficulty precisely.

## The elementary lemma (proved)

> **Lemma.** Let `v = (v₁,…,v₁₃)` be primitive with **no** `vᵢ ≡ 0 (mod 14)`, and `M(v) = 1/14`. Then `t* =
> 1/14` attains `M`, and the phases `{vᵢ/14 mod 1} = {(vᵢ mod 14)/14} ⊆ {1/14, 2/14, …, 13/14} = {k·(1/14)}` —
> an arithmetic progression. Hence `tight ⟹ {kα}` (with `α = 1/14`).

*Proof.* No `vᵢ ≡ 0 (mod 14)` ⟹ `vᵢ mod 14 ∈ {1,…,13}` ⟹ `‖vᵢ/14‖ = min(vᵢ mod 14, 14−vᵢ mod 14)/14 ≥ 1/14`
for every `i`. So `min_i ‖vᵢ/14‖ ≥ 1/14`, hence `M(v) ≥ 1/14`. Given `M(v)=1/14`, equality holds at `t*=1/14`,
where the phases are multiples of `1/14` in `{1..13}/14`. ∎

No equioscillation, no three-gap machinery — the `{kα}` structure on the tight side is the **trivial `t=1/14`
witness**. (Verified: 1358/1358 no-mult-14 families satisfy the core inequality.)

## The reduction: tight-side {kα} = the clean-ruler residual

The lemma applies the moment the family has no multiple of 14. So

> `tight ⟹ {kα}` **holds as soon as** `tight ⟹ no multiple of 14`
> ⟺ `mult-of-14 ⟹ M ≠ 1/14` ⟺ `mult-of-14 ⟹ M > 1/14` (loose).

That last statement is exactly a face of the **clean-ruler residual** the fleet is already proving. So on the
tight side the Ostrowski `{kα}` rigidity is **not a separate hard problem — it is the residual**, dressed
differently. Verified:

- **`tight ⟹ no-mult-14`**: 275/275 tight families have no multiple of 14; **zero** tight families carry one.
  (First observed opus-S249; now given its reason.)
- **`mult-of-14 ⟹ loose`**: 11995/11995 mult-of-14 primitives have `M > 1/14`; zero at `1/14`, zero below. The
  near-AP stress test `{1..12,14}` gives `M = 1/13` (loose, witnessed at `q*=13` since `14 ≡ 1 mod 13` — the
  peeled AP); the deep well `{1..12,182}` gives `M = 14/183`.

## Relocation: the equioscillation and the hard {kα} are on the covering side

The equioscillation (mac-mini S40) is a **two-point** phenomenon: at the optimizer `t*` exactly two runners
bind, pinning `t*` by a pair at distance `M`. That pins `t*`; it does **not** force the *full* config onto a
`{kα}`-progression. On the tight side this is moot (the `t=1/14` witness already delivers `{kα}`). The
genuinely hard `{kα}` is on the **covering side**:

- The `M`-minimizing covering family is the deep well `{1..12,182}`, whose config is `{kα}` (`α = 14/183`)
  **only because its core `{1..12}` is an interval** — a generic covering family has `g = 5` (not `{kα}`,
  mac-mini S38).
- So "the extremal covering config is `{kα}`" ⟺ "the covering-min is the interval-core deep well" ⟺ the
  covering-min bound `M ≥ 14/183` — the residual, which (mac-mini S40) needs a **dual** (Delsarte / de la
  Vallée-Poussin positive-polynomial) certificate, *not* equioscillation-greedy (which sticks in local maxima).

## Net

Working the route honestly: **"tight ⟹ {kα}" is elementary** (the `t=1/14` witness) and, via the reduction,
**equivalent to the clean-ruler residual** on the tight side — so it was never a separate open problem. The
equioscillation tool and the hard `{kα}` both belong to the **covering side**, where "extremal ⟹ `{kα}`" is
the covering-min bound `M ≥ 14/183`, needing the dual certificate. This **unifies** the Ostrowski `{kα}`
rigidity with the residual (same statement, tight side) and places the remaining difficulty squarely on the
covering-min dual certificate — confirming, from the rigidity direction, that the entire game is the
mult-of-14 / divisor-complete residual, off which everything (both the `M ≥ 1/14` bound and the `{kα}` rigidity)
is elementary via `t = 1/14`.

→ mac-mini S38 (Ostrowski ladder / open `{kα}`), mac-mini S40 (2-point equioscillation, dual certificate),
klein S267 (14/183 covering-min), THM-527 (three-gap), opus-S251 (`{AP,V*}` = full/punctured `{k/14}`),
opus-S249 (`tight ⟹ no-mult-14`, first observed), opus-S246 (divisor-complete loose), S235 (band-edge). Files:
`lrc14_tight_implies_ka_elementary_opus_S252.py` (+`.out`).
