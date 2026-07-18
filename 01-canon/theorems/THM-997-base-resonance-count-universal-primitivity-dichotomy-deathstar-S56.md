# THM-997 — The base-resonance count is universal; the all-resonance count detects primitivity (death-star-2026-07-17-S56)

**Status:**
- **Theorem 1 (PROVED, fully general):** every tight family with no multiple of `n` has
  `liveCount(n) = φ(n)` exactly, live set `= (ℤ/n)*`. Answers "the exact φ(n) resonance count for all
  tight families" at the defining resonance `q = n`.
- **Theorem 2 (the primitivity dichotomy):** `liveCount(nm) = Σ_{m''|m} ℓ(nm'')`, `ℓ(n) = φ(n)`.
  For *difference-closed primitive* tight families (the AP): `ℓ(nm'') = 0` for `m'' ≥ 2`, so
  `liveCount(nm) = φ(n)` for **all** `m` (PROVED). For *non-primitive* tight families this is
  **FALSE** (PROVED counterexample `3·AP`). General primitive case (GW): VERIFIED, residual open.

Continuation of THM-996 (this corrects the naive reading "φ(n) at every resonance for every tight
family"). Source HYP-7305. Script `04-computation/lrc_resonance_count_deathstar_S56.py` (+ `.out`).

Definitions as in THM-996 / `LRCDiscreteBonferroni.lean`:
`liveCount_V(q) = #{p ∈ (0,q) : ∀v∈V, ‖v p/q‖ ≥ 1/n}`, tight `⟺ M(V) := sup_t min_v‖vt‖ = 1/n`.

---

## Theorem 1 — the base-resonance count is φ(n) for every tight family

**Statement.** Let `V ⊂ ℤ∖{0}`, `|V| = n−1`, be tight (`M(V) = 1/n`) with no element divisible by
`n`. Then `liveCount_V(n) = φ(n)`, and the live multipliers are exactly the units `(ℤ/n)*`.

**Lemma (tightness ⟹ divisor-covering).** For every divisor `d | n` with `1 < d < n`, `V` contains a
multiple of `d`.
*Proof.* If not, then for every `v ∈ V`, `d ∤ v`, so `‖v·(1/d)‖ = ‖(v mod d)/d‖ ≥ 1/d`. Hence
`min_v ‖v/d‖ ≥ 1/d > 1/n` (as `d < n`), so `M(V) ≥ 1/d > 1/n`, contradicting tightness. ∎
(Confirmed: the AP with its multiple-of-7 removed has `M = 1/7`.)

**Proof of Theorem 1.** A multiplier `p ∈ (0,n)` is live at `q = n` iff `‖v p/n‖ ≥ 1/n` for all `v`,
i.e. `n ∤ vp` for all `v`. Let `g = gcd(p,n)` and `d = n/g`; since `gcd(p/g, d) = 1`,
`n | vp ⟺ d | v`. So **`p` is live ⟺ no `v ∈ V` is divisible by `d = n/gcd(p,n)`.**

- **`p` a unit** (`g = 1`, `d = n`): live ⟺ no multiple of `n` in `V` — true by hypothesis. All `φ(n)`
  units are live.
- **`p` a non-unit** (`g > 1`, so `1 < d = n/g < n` is a proper divisor): by the Lemma `V` has a
  multiple of `d`, so `p` is **not** live.

Hence the live set at `q = n` is exactly `(ℤ/n)*`, of size `φ(n)`. ∎

This is the "exact φ(n)" statement in full generality: it needs only tightness (through
divisor-covering) and the absence of a multiple of `n` — **no difference-closure, no classification of
the tight locus.** It holds identically for the AP, GW, and every coprime dilate (all four verified:
`liveCount(14) = 6`, live set `= {1,3,5,9,11,13}`).

---

## Theorem 2 — the all-resonance count and the primitivity dichotomy

**The decomposition.** `liveCount(nm)` counts loneliness times `t ∈ (0,1)` with denominator dividing
`nm`. Grouping by `g = gcd(p, nm)` and writing `N = nm/g`, a multiplier `p` reduces to a loneliness
time of exact denominator `N`; by THM-996-I loneliness forces `n | N`, and `N | nm` gives `N = nm''`
with `m'' | m`. Writing `ℓ(N) = #{reduced-denominator-N loneliness times}`,
```
liveCount(nm) = Σ_{m'' | m} ℓ(nm''),   with  ℓ(n) = φ(n)   (Theorem 1).
```
So `liveCount(nm) = φ(n)` for all `m` **⟺** `ℓ(nm'') = 0` for all `m'' ≥ 2` (no loneliness time of
denominator `> n`). Call this property **(A)**.

**(A) holds for difference-closed primitive tight families (PROVED).** If `V` is difference-closed
(⟺ a dilated AP; the primitive one is `{1,…,n−1}`) and lonely at `t = a/(nm')` (reduced, min `= 1/n`),
then the `n` points `{0} ∪ {v a mod nm' : v ∈ V}` are pairwise `≥ m'` apart — the observer-runner gaps
are `≥ m'` by loneliness, and the runner-runner gaps `‖(v_i−v_j)a/(nm')‖ ≥ 1/n = m'/(nm')` because
`v_i − v_j ∈ ±V`. Then `n` points pairwise `≥ m'` apart on a circle of circumference `nm'` are forced
to be **equally spaced at exactly `m'`**, so `v_i a ≡ k_i m' (mod nm')`, and `gcd(a,nm') = 1` gives
`m' | v_i` for all `i`, i.e. `m' | gcd(V)`. For a **primitive** dilated AP (`gcd(V) = 1`) this forces
`m' = 1`: no denominator `> n`. Hence `liveCount(nm) = φ(n)` for all `m`. (This is the true content of
THM-996-II; verified AP to `m ≤ 45`.)

**(A) FAILS for non-primitive tight families (PROVED counterexample).** `3·AP = {3,6,…,39}` is tight
(`M(3·AP) = M(AP) = 1/14`, dilation-invariant), has no multiple of `14`, yet
`liveCount(42) = 18 ≠ 6`: its live times at `q = 42` are `6` at reduced denom `14` **plus `12`
"dilation ghosts" at reduced denom `42`** (`ℓ(42) = 12`). The ghosts are the loneliness times
`t ≡ 5o/14 (mod 1/3)` of the dilation. So `φ(n)`-at-every-resonance is **not** a property of all tight
families — it exactly separates the primitive ones from their dilates.

**General primitive case (GW).** GW `= {1,…,11,13,24}` is primitive but **not** difference-closed
(`24 − 1 = 23 ∉ GW`), so the equal-spacing proof does not apply — yet `liveCount(14m) = 6` for all
`m ≤ 45` (VERIFIED). The residual, conjectural statement:

> **(A) for all primitive tight families:** a primitive tight family has no loneliness time of
> denominator `> n` (equivalently its loneliness locus is the single unit-orbit `{o/n : o ∈ (ℤ/n)*}`).

The difference-closed proof shows a denominator-`nm'` loneliness time forces `m' | gcd(V)` *when the
family is difference-closed*; the conjecture is that some common factor is forced for **every** tight
family with a high-denominator loneliness time (the dilation ghost is the only mechanism). This is a
rigidity-flavored statement, of a piece with the tight-locus classification `{AP, GW}` (THM-995 VII).

---

## Summary (what "the exact φ(n) resonance count for all tight families" actually is)

| resonance | statement | status |
|---|---|---|
| `q = n` (base) | `liveCount(n) = φ(n)`, live set `(ℤ/n)*` | **PROVED, every tight family** (no mult of `n`) |
| `q = nm`, difference-closed primitive (AP) | `liveCount(nm) = φ(n)` | **PROVED** (equal-spacing ⟹ `m'∣gcd`) |
| `q = nm`, non-primitive tight (`3·AP`) | `liveCount(nm) > φ(n)` possible | **PROVED false** (ghosts, `ℓ(42)=12`) |
| `q = nm`, general primitive (GW) | `liveCount(nm) = φ(n)` | VERIFIED; residual = property (A) |

The universal, unconditional result is Theorem 1 (base resonance). The all-resonance count is a
*primitivity detector*, and its proof for a general primitive tight family is the same rigidity core
that the equality horn (THM-995) confronts — consistent with THM-996 Part III (the loneliness census
carries no more separating power than the tight-locus classification itself).

→ THM-996, THM-991, THM-995 (VII), HYP-7305, MISTAKE-100.
