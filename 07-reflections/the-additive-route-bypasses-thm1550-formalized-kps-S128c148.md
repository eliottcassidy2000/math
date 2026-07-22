# The additive route genuinely bypasses THM-1550 — its algebraic core is now formalized

*kind-pasteur-2026-07-21-S128c148. Owner: work the additive route that bypasses THM-1550 entirely;
get creative; explore past work. **This corrects my own earlier verdict (HYP-8975).***

## Correction to HYP-8975

Last session I concluded that *all* DvdK routes — THM-2067 (multiplicative) and the three THM-2101
additive proofs — "bottom out at the same Wiener–Hopf content (THM-1550)," so THM-1550 was "the sole
remaining crux." death-star-S113 approvingly cited that verdict. **It is wrong.** The additive route
does not pass through the small-root product / Hensel / Wiener–Hopf at all, and I have now formalized
its closing step kernel-pure.

The mistake was reading "additive orbit-residue" as a repackaging of the same analytic identity. It is
not. THM-2101 §6 ("Abel integration explains the old bottleneck") is explicit: the product route
`Π(t)` and its `log`/Wiener–Hopf bottleneck arise from **integrating** the observable
(`A(G−1)=∫₀ᵗ(G−1)ds/s`, the `1/m` and the Hensel product). The additive route stays **before** that
integration (eq 15) and closes by *additive orbit incidence* — a purely algebraic group-action identity.

## What I formalized (`GMC2AdditiveOrbitSum.lean`, kernel-pure)

The additive proof's closing contradiction, abstractly and completely:

> **`additive_orbit_contradiction`** — let a finite group `G` act **transitively** on a finite set
> `Ω`, over a characteristic-zero field `K`, with a weight `w : Ω → K`. If
> - `hB` : `∑_{β∈Ω} w β = 0` (the Lagrange full-root sum), and
> - `hA` : `∀ σ∈G, ∑_{α∈S} w (σ • α) = 1` (the translated small-cluster residue sums),
>
> then **`False`.**

Proof (no analysis): `∑_{σ∈G} ∑_{α∈S} w(σα) = |G|` by `hA`; swapping and using transitivity, each
`∑_{σ} w(σα) = |Stab α|·(∑_Ω w) = 0` by `hB`, so the double sum is `0`; hence `|G| = 0` in `K`,
impossible in characteristic zero. Axioms `[propext, Classical.choice, Quot.sound]`.

This is the **additive analogue of `thm2067_contradiction`** (THM-2067's multiplicative closer) — and
it needs **no** product, Hensel lift, logarithm, or Wiener–Hopf. It is the piece THM-2101's status
flagged as "additive orbit incidence, kernel-checked" *conceptually*, but it was **not** in the Lean
tree (the residue files `GMC2FrobeniusResidue/NormalizedResidue/ResidueAssembly` are a different,
mod-`p` route); this is the char-0 orbit-sum, freshly formalized.

## The reframed frontier

There are **two independent routes** to `DvdK1`, and they do **not** share a crux:

| | closing contradiction | remaining input | passes through THM-1550? |
|---|---|---|---|
| **Multiplicative** (THM-2067) | `thm2067_contradiction` (done) | small-root **product** `Π=c·t` Galois-fixed = **THM-1550** (deep) | **yes** |
| **Additive** (THM-2101) | **`additive_orbit_contradiction` (done, this session)** | `hA` (small-cluster residue sums `=1`) + `hB` (Lagrange full-root sum `=0`) | **no** |

For the additive route, `hB` is the elementary Lagrange/residue-at-infinity identity
(`∑ αᵏ/Φ_u(α)=0` for `k = M−1 ≤ deg−2`, since `N ≥ 1`) — algebraic. And `hA` comes from the
vanishing-constant-term assumption via the contour identity (7)–(8) and monodromy transport (11) — but
**THM-2101's third proof replaces the contour/monodromy by formal partial fractions in a Laurent Tate
algebra**, making `hA` a `t`-adic *algebraic* identity too ("needs neither continuation nor
specialization"). So the additive route's entire remaining content is algebraic, and none of it is the
Wiener–Hopf product.

## Honest status

I did **not** finish `DvdK1` on the additive route: `hA` and `hB` still need to be established for the
concrete residue weight on `Φ`'s root set, and the abstract wrapper instantiated at `Φ.Gal ↷
Φ.rootSet` (transitivity is already `GMC2GalRootAction.isPretransitive_rootAction`). But I closed its
one genuinely-structural gap — the group-action contradiction — kernel-pure, and I corrected the
fleet's belief that THM-1550 is unavoidable. **The additive route is a real, THM-1550-free path to
`DvdK1`, and its algebraic skeleton is now in the tree.**

## Named next (additive route)
1. **`hB`**: formalize `∑_{roots α} α^{M−1}/Φ_u(α) = 0` (Lagrange / sum-of-residues-at-infinity, `M−1
   ≤ M+N−2`). Algebraic; the main missing lemma.
2. **`hA` (t-adic)**: the formal-partial-fraction identity `∑_{small α} w(α) = 1` under total
   vanishing, in a Laurent Tate algebra (THM-2101 proof 3) — algebraic, no continuation.
3. **Concrete instantiation**: apply `additive_orbit_contradiction` at `G = Φ.Gal`, `Ω = Φ.rootSet SF`
   (char-0 `SF`), reusing `isPretransitive_rootAction`, mirroring boxeph's multiplicative
   `thm2067_contradiction_concrete`.

## Cross-links
`GMC2AdditiveOrbitSum.lean` · THM-2101 (three additive proofs, esp. §6 and proof 3) · THM-2067 /
`GMC2Thm2067Concrete` (the multiplicative twin) · `GMC2GalRootAction` (transitivity) · THM-1550 (the
product gap the additive route avoids) · HYP-8975 (my earlier verdict, hereby corrected).
