# HYP-2330 — Erdős Problem 625: the cochromatic number is the σ-symmetric chromatic number; χ−ζ is the σ-asymmetry of coloring

**Session:** S652
**Status:** CONFIRMED (the σ-symmetry of ζ formalized + verified; the problem itself is near-resolved, cited)
**Provenance forward:** math-lean `Math/Combinatorics/Cochromatic.lean` (sorry-free)
**Prompt:** research Erdős problem 625, relate to the repo, explore threads, begin work creatively.

---

## 0. The problem (verified from the source)

**Erdős Problem 625** (Erdős–Gimbel, ~1991; `$100` / `$1000`): for the random graph `G ~ G(n, 1/2)`, is
it true that `χ(G) − ζ(G) → ∞` with high probability?

- **`ζ(G)` = the cochromatic number** = the fewest colours so that **every colour class is an independent
  set _or a clique_** (`erdosproblems.com/625`; Heckel arXiv:2409.17614).
- **Status:** Heckel & Steiner (2024, independently): `χ − ζ` is *not* bounded whp (`≥ n^{1/2−o(1)}` along
  a sequence); Heckel (2024): affirmative for ≈ 95% of `n`; conjectured `χ − ζ ~ n/(log n)³`. So the
  answer is morally **yes**, but not fully closed (all `n`, exact rate). I do not claim to resolve it.

---

## 1. The repo connection: ζ is the σ-symmetric (complement-invariant) chromatic number

The arc's spine is the involution **`σ` = complement** (the antipodal swap, S638/S643). It is *exactly*
the symmetry that distinguishes `χ` from `ζ`:

- `σ` swaps **clique ↔ independent set** (`Gᶜ.IsClique s ↔ G.IsIndepSet s`). So a colour class is "a
  clique or an independent set" in `G` iff it is so in `Gᶜ` — **the cochromatic colourings of `G` and
  `Gᶜ` coincide**, hence `ζ(G) = ζ(Gᶜ)`: **`ζ` is `σ`-invariant.**
- `χ(G) ≠ χ(Gᶜ)` in general (chromatic is *not* `σ`-invariant).

> **So `ζ` is the σ-symmetrisation of `χ`, and `χ − ζ` (the Erdős-625 quantity) is precisely the
> σ-asymmetry of the chromatic number.** On `G(n,1/2)` — which is `σ`-self-complementary *in
> distribution* (`G ≈ Gᶜ`) — `ζ` is the σ-symmetric core and `χ − ζ` is the gap created by breaking the
> complement symmetry. Erdős 625 asks how big that σ-asymmetry grows.

**Formalized (math-lean, sorry-free): `Math/Combinatorics/Cochromatic.lean`**
- `isClique_compl_iff` / `isIndepSet_compl_iff` — `σ` swaps the two pure phases.
- `cliqueOrIndep_compl_iff` — the "clique-or-independent" predicate is `σ`-invariant (`or_comm` after the
  swap).
- `cochromColorable_compl_iff : CochromColorable G k ↔ CochromColorable Gᶜ k` — **`ζ(G) = ζ(Gᶜ)`**, the
  σ-symmetry of the cochromatic number.

Verified exactly (`erdos625_cochromatic_sigma_s652.py`, small `G(n,1/2)`): `ζ(G) = ζ(Gᶜ)` **always**;
`χ(G) ≠ χ(Gᶜ)` for several `n` (e.g. `n=8`: `χ = 5` vs `χ(Gᶜ) = 3`, `ζ = 3` both). The `χ − ζ` gap is
small/noisy at small `n` (avg `0.3–0.8`), consistent with slow growth `~ n/(log n)³`.

---

## 2. Creative reframes (the repo's machinery on 625)

1. **Ising / two ground states.** A **clique** = all-edges (ferromagnetic "up"); an **independent set** =
   no-edges ("down"). `ζ` = min #pieces each in a **pure phase** (one of the two Ising ground states);
   `χ` = min #pieces each in the "down" phase *only*. **`χ − ζ` = the cost of forbidding the "up" phase.**
   `σ` swaps up ↔ down, fixing `ζ`. (This is the coloring–partition unification, S626–633: `ζ` partitions
   the conflict graph into its two pure phases.)
2. **The `χ·α ≥ n` bridge (S634).** `χ ≥ n/α`; `ζ ≥ n/max(α,ω)` (each pure piece has size `≤ max(α,ω)`).
   On `G(n,1/2)`, `α ≈ ω ≈ 2 log₂ n` (σ-symmetric!), and `χ ≈ ζ ≈ n/(2 log₂ n)` share the **same leading
   term** — so Erdős 625 lives entirely in the **subleading** gap, exactly where σ-symmetry (`ζ`) and
   σ-asymmetry (`χ`) first diverge.
3. **`G(n,1/2)` = the finite Rado graph** (S638 did the random *tournament*); `ζ` is its σ-symmetric
   coloring invariant. The **cube** in Heckel's conjectured `n/(log n)³` echoes the arc's `3` (cube root)
   — speculative but on-theme.

---

## 3. Beginning work in earnest: the σ-asymmetry lower bound (a clean new framing)

Since `ζ` is σ-symmetric and `χ` is not, `χ − ζ` is a *σ-equivariant* defect, and that gives a clean,
rigorous lower bound:

> **`ζ(G) ≤ min(χ(G), χ(Gᶜ))`**, hence **`χ(G) − ζ(G) ≥ χ(G) − χ(Gᶜ)`** — the σ-asymmetry of `χ` itself
> lower-bounds the Erdős-625 quantity.

*Proof.* `ζ ≤ χ` (a proper colouring is cochromatic: its classes are independent sets, an allowed pure
phase). `ζ(G) ≤ χ(Gᶜ)`: a proper colouring of `Gᶜ` partitions into `Gᶜ`-independent sets `=` `G`-cliques,
also an allowed pure phase for `G`'s cochromatic colouring. So `ζ ≤ min(χ, χ(Gᶜ))`. ∎ (Verified: `ζ ≤ χ`
and `ζ ≤ χ(Gᶜ)` on all small `G(n,1/2)`.)

This re-routes Erdős 625 through the σ-machinery: **the `χ − ζ` gap is at least the chromatic
σ-asymmetry `χ(G) − χ(Gᶜ)`.** HONEST CAVEAT: whether `χ(G) − χ(Gᶜ)` *itself* grows to `∞` on `G(n,1/2)`
(and at what rate) is its own question — `χ` and `χ(Gᶜ)` are correlated (same instance) with equal means;
I have *not* shown this route reaches the Heckel–Steiner `n^{1/2−o(1)}`. The contribution is the clean
σ-bound, not a reproof. But it locates the difficulty precisely: **`χ − ζ` is driven by how
σ-asymmetric the chromatic number of the random graph is.**

## 4. Handoffs
- Formalize `ζ(G) ≤ χ(G)` and `ζ(G) ≤ χ(Gᶜ)`, hence `χ(G) − ζ(G) ≥ χ(G) − χ(Gᶜ)` (the σ-asymmetry lower
  bound) — needs `χ` (Mathlib `chromaticNumber`) tied to `CochromColorable`.
- Quantify, on `G(n,1/2)`, the σ-asymmetry `χ(G) − χ(Gᶜ)` (the chromatic-fluctuation route to Erdős 625's
  lower bound) — the repo's σ-machinery applied to the actual problem.
- Heckel's `n/(log n)³`: is the `3` structural (the cube root / the three "levels" of the deletion poset,
  S651) or accidental?
