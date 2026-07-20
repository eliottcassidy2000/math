# THM-1460 — The oriented-spanning-tree switching-sum (Hamiltonian paths and arborescences are two instances)

**Status:** PROVED (verified n=3,4,5 exactly for four families; general proof below).
**Author:** death-star-2026-07-20-S61c (HYP-8315). (Subsumes the switching-class
Hamiltonian-path identity Σ H = n! that I first filed as "THM-1445" in S61b; that number was
already taken by opus/kp who pushed 11 min earlier, so it is ceded and the result lives here.)

## Statement

Fix n ≥ 2. Call S an **oriented spanning tree** of K_n if S is a set of n−1 arcs whose
underlying edges form a spanning tree of K_n (so S is a spanning tree together with a chosen
direction on each of its edges). A tournament T (an orientation of every edge of K_n)
**realizes** S, written S ⊆ T, if T's arc on each edge of S agrees with S's direction.

Let 𝓕 be **any** set of oriented spanning trees of K_n, and let 𝒞 be any switching class of
tournaments (orbit under vertex switching = coset of the cut space B = ⟨δ(v)⟩ ⊆ F₂^{E(K_n)}).
Then

  **Σ_{T ∈ 𝒞} #{ S ∈ 𝓕 : S ⊆ T } = |𝓕| .**

## Proof

A switching class is a coset D₀ + B of the cut space B, dim B = n−1. Fix S ∈ 𝓕; its edge set
E(S) is a spanning tree of K_n. Restriction to the tree edges, ρ : F₂^{E(K_n)} → F₂^{E(S)},
maps B **isomorphically** onto F₂^{E(S)} — the fundamental cut-sets of the spanning tree E(S)
form a basis of B, one per tree edge, each containing exactly its own tree edge, so ρ carries
that basis to the standard basis. Hence within the coset 𝒞 the tree-edge orientations ρ(T)
range over all 2^{n−1} patterns, **each exactly once**; in particular **exactly one** T ∈ 𝒞
agrees with S on E(S), i.e. realizes S. Therefore
  Σ_{T∈𝒞} #{S∈𝓕 : S⊆T} = Σ_{S∈𝓕} #{T∈𝒞 : S⊆T} = Σ_{S∈𝓕} 1 = |𝓕|. ∎

## Instances (all verified n = 3,4,5)

| family 𝓕 | |𝓕| | class-sum |
|---|---|---|
| directed Hamiltonian paths | n! | **n!** (the Ham-path corollary) |
| Ham paths starting at fixed r | (n−1)! | **(n−1)!** |
| out-arborescences (rooted spanning trees, any root) | n^{n−1} | **n^{n−1}** |
| out-arborescences rooted at fixed r | n^{n−2} | **n^{n−2}** (Cayley) |

So **the Ham-path identity (Σ H = n!) is the "path" instance and the total-arborescence
count (→ n^{n−1}) is the "all-trees" instance of the same theorem.** Hamiltonian paths are exactly
the *tallest* (depth n−1, path-shaped) out-arborescences; the identity is blind to the shape
of the tree — it depends only on the spanning-tree edge set, which the cut space equidistributes.

## Corollaries

- **Global sums** (sum over all 2^{C(n−1,2)} classes): Σ_T H(T) = n!·2^{C(n−1,2)} and
  Σ_T a(T) = n^{n−1}·2^{C(n−1,2)}, where a(T) = Σ_r (out-arborescences rooted at r). Both
  factor as (combinatorial family size) × (the even-graph / cycle-space count 2^{C(n−1,2)}).
- **Ratio** Σ_T a(T) / Σ_T H(T) = n^{n−1}/n! — the even-graph factor cancels; paths are a
  fraction n!/n^{n−1} of all rooted spanning trees (exponentially rare, yet both average to a
  pure combinatorial constant over each switching class).

## Cross-links
the Ham-path corollary Σ_𝒞 H = n! · odd/even reflection S61b (H is odd-valued; even-graph projection
constant) · kp THM-1430 / THM-474 (tilings = switching classes; cut space simply transitive on
spanning-tree orientations) · Matrix-Tree / Tutte (a_r = Laplacian cofactor) · Cayley n^{n−2} ·
BEST theorem (regular tournaments) · HYP-8315.
