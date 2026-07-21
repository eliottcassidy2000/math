---
id: THM-1750
title: "THE PER-ROOT ARBORESCENCE VECTOR IS A POLY-TIME RANKING (= the stationary distribution of the dominance random walk, MCTT) THAT SPLITS THM-1580's UNBREAKABLE RESIDUE — but is complementary to H, not dominant. (1) On strongly connected tournaments {a_r}/Σa is EXACTLY the stationary distribution of the against-the-arcs 'who-dominates-me' walk (Markov-Chain-Tree Theorem, unit rates), verified to machine precision on all n=5 classes — so the arborescence vector is a genuine tournament centrality ranking, source-heavy (the transitive source takes all the weight). (2) {a_r} refines Σa and SPLITS BOTH cospectral groups THM-1580 flagged as resisting spec(A)+Σa+H: (Σa,H)=(1680,47) and (2380,143), each into two distinct sorted vectors. So the exact wall THM-1580 named unbreakable 'by the adjacency spectrum, the arborescence count, or the Hamiltonian-path count' IS broken — by the arborescence RANKING (the vector, not its sum). (3) As a combined fingerprint the finer vector drops the (spec A, ·, H) indistinguishable-pair count from 27 (with Σa) to 4 (with {a_r}). (4) HONEST CEILING: (spec A, Σa, {a_r}, H) is still NOT complete at n=7 — a NEW 4-pair residue resists it, (Σa,H) = (1080,23),(1088,29),(1092,31),(1224,43). The vector and the #P-hard path count are COMPLEMENTARY (THM-1580's theme at the vector level), neither dominant."
status: >
  (1) PROVED (MCTT is classical) + VERIFIED to machine precision: {a_r}/Σa equals the
  stationary distribution of the unit-rate continuous-time generator with rates along the
  REVERSED arcs, on every strongly connected n=5 class (0 mismatches). Direction matters:
  a_r counts OUT-trees from r, which is the against-the-arcs walk's stationary measure.
  (2) VERIFIED-EXACT (exact integer arborescence cofactors): both THM-1580 residue groups
  carry two distinct sorted {a_r} vectors — explicit vectors in §2.
  (3),(4) VERIFIED-EXACT, exhaustive over all 456 classes at n=7, exact integer arithmetic
  (Faddeev-LeVerrier char poly, Matrix-Tree cofactors, direct hp).
  SELF-CORRECTION on record: an interim summary claimed {a_r} "resolves the entire n=7
  residue" — false, it conflated THM-1580's 2 monochromatic cospectral GROUPS with the
  combined fingerprint's 27 collision PAIRS; corrected in §3–4 (the 4-pair residue).
  Advances no open problem; sharpens THM-1580 and adds the ranking identity.
source: klein-2026-07-20-S379 (owner: work another cutting-edge math session, think arborescences)
depends_on:
  - THM-1460  # arborescences = determinantal relaxation of H; Σa is Laplacian-spectral
  - THM-1580  # the poly-time-beats-#P-hard 14x finding + the 2-group residue this splits
related:
  - THM-923   # arborescence transfer / flip law
script: 04-computation/arborescence_vector_klein_S379.py (+ .out)
---

# THM-1750 — the arborescence vector is a ranking, and it splits THM-1580's residue

## 1. The vector is the dominance-walk stationary distribution (a ranking)

`a_r` = number of spanning arborescences out of root `r` (the `(r,r)` cofactor of the
in-Laplacian `L_in = D_in − A`; `Σ_r a_r` is THM-1460's Laplacian-spectral count). The **vector**
`(a_r)_r` — not its sum — is the object here.

**Markov-Chain-Tree Theorem.** For the continuous-time unit-rate walk whose rate along the
*reversed* arcs is 1 (from `v`, flow toward the vertices that **dominate** `v`), the stationary
distribution is `π_r ∝ a_r`. Verified to machine precision on every strongly connected `n=5`
class (0 mismatches). So:

> `{a_r}` normalized **is** the stationary distribution of the "who-dominates-me" random walk —
> a poly-time tournament **centrality ranking**, source-heavy: in the transitive tournament the
> source takes all the weight, `(n−1)!` at one vertex.

This ties the arborescence thread (THM-1460/1580) to the project's core theme — ranking a
tournament — and gives an alternative to the Hamiltonian/Kendall ranking that is `O(n^4)` exact.

## 2. It splits THM-1580's unbreakable residue

THM-1580 isolated **two** cospectral groups at `n=7` whose entire membership shares `(Σa, H)`,
calling them the wall for any "determinantal-plus-path" invariant — unresolvable "by the
adjacency spectrum, the arborescence count, or the Hamiltonian-path count." The arborescence
**vector** resolves both:

| `(Σa, H)` | the two `{a_r}` sorted vectors |
|---|---|
| `(1680, 47)` | `(48,72,120,240,240,240,720)` vs `(60,96,100,144,200,480,600)` |
| `(2380, 143)` | `(140,220,236,340,348,532,564)` vs `(172,196,236,340,348,476,612)` |

Same spectrum, same total `Σa`, same `H` — **different arborescence rankings.** The wall was a
wall for the *sum*; the *distribution* sees through it. (The directed-Laplacian char poly
`charpoly(L_in)` does **not** — it fails on both, so it is the vector's per-root resolution, not
the Laplacian spectrum, that does the work.)

## 3. Combined-fingerprint gain

Replacing the scalar `Σa` by the vector `{a_r}` in the poly-time fingerprint, over all 456
classes at `n=7`:

```text
(spec A, Σa,   H)  →  27 indistinguishable pairs
(spec A, {a_r}, H)  →   4 indistinguishable pairs
```

The finer vector resolves 23 of the 27, **including** both §2 groups. (`{a_r}` refines `Σa`, so
`4 ≤ 27` is forced.)

## 4. The honest ceiling — a new 4-pair residue

`{a_r}` is **not** a complete fingerprint. Four pairs resist `(spec A, Σa, {a_r}, H)` at `n=7`:

```text
(Σa, H) = (1080, 23),  (1088, 29),  (1092, 31),  (1224, 43)
```

each a pair sharing adjacency spectrum, arborescence count, arborescence **vector**, and
Hamiltonian-path count. These are the true wall for the arborescence-ranking-plus-path
fingerprint — a strictly deeper residue than THM-1580's, and the right next reconstruction
target.

**The shape of the result, honestly.** The arborescence vector is not a silver bullet that
dominates `H`; it is *complementary* to it. It breaks the specific 2-group wall THM-1580 named
(where the sum and the path count both fail), and it is the natural ranking-theoretic invariant
— but a 4-pair residue survives everything poly-time plus `H`. This is THM-1580's "nearly
complementary failure sets" theme, now sharpened to the vector level: hardness and
discriminating power remain independent axes.

*Self-correction:* an interim run's summary claimed the vector "resolves the entire `n=7`
residue." That conflated THM-1580's 2 monochromatic cospectral *groups* with the combined
fingerprint's 27 collision *pairs*; §§3–4 give the corrected counts and the surviving 4-pair
residue.

## 5. Scope

Exhaustive and exact at `n ≤ 7`. Introduces the ranking identity (MCTT) and sharpens the
reconstruction wall from THM-1580's `(Σa, H)`-monochromatic pairs to the `(spec A, Σa, {a_r},
H)`-residue. No claim beyond `n=7`; the 4-pair residue at `n=7` is the concrete open target.

*Files: `04-computation/arborescence_vector_klein_S379.py` (+ `.out`).*
