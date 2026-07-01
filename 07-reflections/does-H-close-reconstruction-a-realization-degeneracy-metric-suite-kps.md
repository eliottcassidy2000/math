# Does H close reconstruction? No — but it cuts the degeneracy ~5× and extends node-identifiability to n=5. A realization-degeneracy metric suite.

*kind-pasteur-2026-07-01-S14. Chasing the owner's question: does adding the Hamiltonian-path count `H` to the merged-metagraph node data close reconstruction (make the blue/black multigraph the unique realization)? Answer: not fully — but the effect is sharp and quantifiable, and it exposes the first irreducible ambiguity at n=6. Here is the honest verdict plus a suite of realization-degeneracy metrics.*

## Three notions of "H closes reconstruction," kept apart

Reconstruction = recover the blue/black multigraph from node data. The answer depends on *what data*:

1. **Per-node `(category, blue-deg, black-deg)` + bare `H` label.** A degree-preserving 2-swap rewires edges without touching any node's `H`, so bare `H` gives **zero** extra constraint: the abstract degeneracy is unchanged. (n=5: `σ = 84` black + `16` blue swaps, with or without the `H` label.) **Bare `H` does not help.**
2. **The `H`-level *connection pattern`** (which `H`-values link, per colour). Requiring 2-swaps to preserve the multiset of endpoint-`H`-pairs cuts the degeneracy sharply: `σ_H = 20+2` at n=5, `15208+54` at n=6 — a **~4–6× reduction**, but **not to zero**. Knowing which `H`-levels connect is a strong partial constraint.
3. **`(category, degree, H)` as a node *fingerprint`.** It is **injective (a complete node-identifier) for n≤5** and **first fails at n=6** (6 twin-pairs). So `H` extends node-identifiability exactly one step past degree alone.

## The numbers

| n | nodes | σ baseline (bk+bl) | σ_H (H-connection) | ratio | twins `(cat,deg,H)` | twins `+c3` |
|---|---|---|---|---|---|---|
| 4 | 3 | **0** (unique) | 0 | — | 0 | 0 |
| 5 | 10 | 100 | 22 | 4.5× | **0** (injective) | 0 |
| 6 | 34 | 85099 | 15262 | 5.6× | **6** | 5 |

Two thresholds emerge, one n apart:
- **Swap-uniqueness** (degree+colour *determines* the graph): holds **only at n=4**.
- **Fingerprint-injectivity** (`(cat,deg,H)` *identifies* every node): holds **through n=5**, breaks at **n=6**.

So `H` buys exactly one more level of identifiability. And it is a genuine break-at-6 — the same n where the "sea" switches on (pure-black↔pure-black lines) and black self-loops appear.

## The first irreducible ambiguity: the n=6 twins

Six nodes collapse into three... actually **six twin-pairs** share `(category, blue-deg, black-deg, H)`:

| type | (blue,black,H) | nodes (idx, tiling-count, c₃) | c₃ separates? |
|---|---|---|---|
| pure-black | (0, 2, 3) | (1, 2, 1), (3, 2, 1) | no |
| pure-black | (0, 18, 9) | (4, 18, 3), (6, 18, 3) | no |
| **mixed** | (5, 12, 17) | (12, 17, 4), (22, 17, 4) | no |
| pure-black | (0, 46, 23) | (21, 46, **5**), (24, 46, **6**) | **yes** |
| pure-black | (0, 74, 37) | (38, 74, 7), (42, 74, 7) | no |
| **mixed** | (7, 22, 29) | (43, 29, 6), (44, 29, 6) | no |

Adding `c₃` breaks only one pair (21,24); **five twin-pairs survive `(cat, deg, H, c₃)`** — genuinely distinct tournament classes with identical low-order invariants. These are the irreducible reconstruction obstruction: no bounded local fingerprint we tried separates them. Notably two of the six are **mixed** (SC) nodes — the ambiguity is not confined to the anonymous sea.

## Why the graph is still "reconstructible from itself"

A separate fact: **1-WL colour refinement seeded by category alone is already discrete** (34 distinct colours at n=6). So the metagraph is asymmetric and 1-WL-identifiable *from its own adjacency*. The tension is exactly the point: the graph *is* determined once you know its edges (WL sees them), but it is **not** determined by node fingerprints `(cat, deg, H, c₃)` — the twins prove local data is insufficient from n=6. Reconstruction needs a *relational* invariant (the `φ`-pairing / complement-tiling map), not just richer node labels.

## The realization-degeneracy metric suite (proposed trackables)

For any decorated multigraph arising this way, track:
- **σ(n)** — baseline swap number (degree+colour): the raw realization freedom. `0, 100, 85099`.
- **σ_I(n)** — `I`-refined swap number for an edge-invariant `I` (e.g. `I=H`-pair): the freedom *given* `I`. Report the **reduction ratio `σ/σ_I`** (`H` buys ~5×).
- **τ_I(n)** — twin count under node-fingerprint `I`: `τ_{(cat,deg,H)} = 0,0,6`; `τ_{+c₃}=0,0,5`. The **identifiability threshold** = last n with `τ=0` (here 5).
- **WL-level** — coarsest seed whose 1-WL refinement is discrete (here: category, always) — measures intrinsic graph rigidity.
- **the twin list** — the explicit irreducible ambiguity (above), the objects any complete invariant must separate.

These form a clean ladder: `σ` (no invariant) ≥ `σ_H` (H-connection) ≥ … ≥ 0 (full edges), with `τ_I` measuring where node-fingerprints stop identifying. The "reconstruction defect" is `σ_I` for the best local `I`; it is `0` only at n=4, and `H` reduces but never zeroes it for n≥5.

## Verdict and next targets

**`H` does not close reconstruction.** It (a) reduces the abstract degeneracy ~5× via the `H`-connection pattern, and (b) makes `(category, degree, H)` a complete node-fingerprint through n=5, breaking at n=6. The residual is a *relational* obstruction — the five `c₃`-robust twin-pairs and the same-`H` rewiring freedom — that only the complement-tiling `φ`-map resolves.

- **(a) The twin separator.** Find the *minimal* invariant that splits the five surviving n=6 twins (candidates: score sequence, `|Aut|`, the OCF independence polynomial `I(Ω,x)`, or the *neighbour-φ-image* `H(φ(·))`). Whatever separates them is the missing reconstruction coordinate.
- **(b) A relational invariant.** Test whether `(H, H∘φ)` per node — pairing each node with its complement-tiling image's `H` — closes it (this is the edge-relation, not a node label). If `σ_{(H,Hφ)}=0`, the complement-`H`-pairing is the reconstruction key.
- **(c) Asymptotics.** `σ(n)`, `τ(n)`, ratio `σ/σ_H` — do the twins proliferate like the sea (`→` most nodes) or stay sparse?

## Honest status

- **Computed (n≤6):** σ, σ_H, the ~5× reduction, τ = 0,0,6 (identifiability breaks at n=6), the explicit twin list, c₃ splits 1 of 6, WL-discrete from category.
- **Conclusion:** `H` is a strong *partial* reconstruction invariant (one extra identifiability level + ~5× degeneracy cut), not a complete one; the obstruction is relational (`φ`).
- **Open:** the twin separator (a); the `H∘φ` relational test (b); asymptotics (c).

— Related: `buckets-and-pairs-the-merged-metagraph-as-a-constraint-conjecture-multitude-kps.md` (R2 reconstruction rigidity), `the-blueblack-line-pairing-is-a-degree-tiling-count-realization-kps.md` (the pairing process), THM-549/550 (half-tiling / `φ`), THM-584 (complement=antipodal), `merged-metagraph-invariants.md`. Script: `04-computation/merged_metagraph_H_reconstruction_kps.py` (+ .out). Not a HYP reservation.
