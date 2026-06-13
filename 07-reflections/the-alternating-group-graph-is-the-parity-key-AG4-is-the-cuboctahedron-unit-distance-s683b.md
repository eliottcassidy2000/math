---
source: oracle-2026-06-06-S683b
status: synthesis + verified instance (the alternating group graph as the non-abelian / parity node of the LRC=HN=unit-distance unification; AG_4 = cuboctahedron = a unit-distance graph)
tags:
  - lonely-runner
  - hadwiger-nelson
  - unit-distance
  - alternating-group
  - cayley-graph
  - chromatic-number
  - parity
  - odd-cycle
---

# The Alternating Group Graph Is the Parity Key — and AG_4 *Is* the Cuboctahedron (a Unit-Distance Graph)

The unification (HYP-2265): LRC, Hadwiger-Nelson, and unit-distance are all **distance Cayley
graphs** `Cay(A,S)`, sharing the **connection-set Fourier / Hoffman-Delsarte spectral bound**
`χ ≥ 1 − λ_max/λ_min`. That frame is for **abelian** `A` (`ℝ²`, `ℤ`, `ℤ_m`). The user asked to
bring in the **alternating group graph** — which turns out to be the unification's
**non-abelian node and its parity key**, with one *literal* shared instance.

## 1. The parity mechanism is the chromatic obstruction (the shared key)

Chromatic number = the failure of a **sign (parity) 2-coloring**, obstructed by **odd cycles**
(the repo's OCF / parity divide). On Cayley graphs this is exact and *visible in the
generators' parity* (`alternating_group_graph_unification_s683b.py`):

| graph | generators | parity | result |
|---|---|---|---|
| `Cay(S_n, transpositions)` | transpositions | **odd** | the sign character `ε:S_n→±1` 2-colors every edge ⟹ **bipartite, χ=2** (verified S_4, S_5) |
| `Cay(A_n, 3-cycles)` = **AG_n** | 3-cycles | **even** | `ε≡+1` on `A_n` ⟹ sign 2-coloring *fails*; order-3 generators force **triangles** ⟹ **χ ≥ ω = 3** |

So **even generators kill the parity 2-coloring and manufacture odd cycles** — the exact
chromatic obstruction that is, verbatim:
- **LRC**: the interval/AP "additive" beat (χ=2, sign-respecting, S582/S626) vs the Paley/QR
  even-structured beat (χ=3);
- **Hadwiger-Nelson**: the unit **equilateral triangle** (3 mutual unit distances = a triangle
  ⟹ `χ(ℝ²) ≥ 3`), and the Moser/de Grey rigid gadgets gluing such triangles.

A triangle in `AG_n` (`g — g·(012) — g·(021) — g`), a unit equilateral triangle, and an LRC
odd resonance are **the same odd cycle**.

## 2. The spectral bound is EXACTLY TIGHT on AG_n (the method works on the non-abelian node)

`AG_4`: 12 vtx, 4-regular, `χ=3`, `ω=3`, **Hoffman `1 − λ_max/λ_min = 1 − 4/(−2) = 3`** (tight).
`AG_5`: 60 vtx, 6-regular, `χ=3`, **Hoffman `1 − 6/(−3) = 3`** (tight). So HYP-2265's shared
spectral key is **exact** here, with the Cayley eigenvalues coming from the **irreducible
representations of `A_n`** (non-abelian Fourier) in place of the abelian characters — exactly
the generalization the unification needs: *connection-set Fourier transform → connection-set
image under the irreps.*

## 3. The literal shared instance: AG_4 = the cuboctahedron = a unit-distance graph

`AG_4` has spectrum `{4, 2³, 0³, (−2)⁵}`, 12 vertices, 4-regular, **8 triangles, 6 squares** —
this is **exactly the cuboctahedron graph** (verified). The cuboctahedron is a genuine
**unit-distance graph**: its 12 vertices are `(±1,±1,0)` and permutations = the **FCC / `A_3`
lattice kissing configuration** (kissing number 12), all edges equal length. Therefore

> **AG_4 is, simultaneously and literally, (i) the Cayley graph of the alternating group `A_4`,
> (ii) a unit-distance graph (the cuboctahedron = FCC kissing graph), and (iii) the χ=3 odd-cycle
> parity object.** The three problems are not merely analogous here — they share one graph.

This concretely lands the unification's geometry: HYP-2170/2265 placed the unit-distance /
de Grey construction on the **Eisenstein lattice `ℤ[ζ_6]`** (2-D triangular, kissing 6); `AG_4`
is the **3-D analogue** (FCC, kissing 12) — the same "kissing / norm-layer" geometry
(`norm-layer-dichotomy-kissing-vs-popular-norm-s702`) one dimension up, and it is the
alternating group that realizes it.

## 4. How each side's insight is the other's key (via the alternating group)

- **Parity (repo's OCF) → both χ's.** The even/odd-generator dichotomy *is* the bipartite /
  non-bipartite dichotomy, i.e. the chromatic number. The repo's odd-cycle machinery is the
  native language of `χ` for LRC, HN, and the alternating group graph at once.
- **HN's rigid gadget ↔ LRC's tight config ↔ AG_n's clique/triangle.** de Grey's finite
  4-/5-chromatic unit-distance subgraph, the LRC tight AP/V*, and `AG_n`'s triangles
  (`ω=3`) are all *finite extremal certificates* pinning the spectral bound (which is tight on
  `AG_n`).
- **Kissing/lattice ↔ resonance shells.** `AG_4`=cuboctahedron=FCC kissing (12); the LRC
  resonance shells and de Grey live on `ℤ[ζ_6]` (kissing 6). The alternating group exhibits
  the lattice kissing geometry as a *group* Cayley graph.

## Caveats / honest notes
- "χ" splits into two invariants and they must not be conflated: the **undirected chromatic
  number** of the distance graph (the unification's `χ`, matching HN; the AP gives the
  *complete* graph `K_m`, χ=m, extremal/tight — consistent with HYP-2265) vs the **dichromatic
  number** of the oriented tournament (`R_m` → 2, S582). The parity story here is about the
  undirected `χ` and is the one shared with HN.
- AG_4 = cuboctahedron is a unit-distance graph in `ℝ³` (the polytope); whether its abstract
  graph also has a *planar* (`ℝ²`) unit-distance embedding is a separate, open-flavoured
  question (not needed for the unification, which lives in any `ℝ^d` / lattice).
- AG_5 (60 vtx, 6-regular, χ=3) being a known polytope/unit-distance graph is plausible (A_5 =
  icosahedral rotation group) but not identified here — a clean next check.

## Verdict / next
- **The alternating group graph is the non-abelian node and the parity key of the
  LRC=HN=unit-distance unification; AG_4 = cuboctahedron = FCC-kissing unit-distance graph is a
  literal shared instance; the Hoffman spectral bound is exactly tight on AG_n** with the
  eigenvalues from `A_n` irreps.
- Next: (1) identify AG_5 / AG_6 as polytope/lattice-kissing unit-distance graphs (the
  alternating-group → kissing-graph functor); (2) push the **odd-cycle (OCF) = chromatic
  obstruction** dictionary as the single shared engine across all three; (3) the non-abelian
  Delsarte LP on `A_n` irreps as the analogue of HN's `J_0`-Bessel and LRC's sinc bounds.

## Artifacts
```
04-computation/alternating_group_graph_unification_s683b.py
05-knowledge/results/alternating_group_graph_unification_s683b.out
```
Related: HYP-2265 (LRC=HN=unit-distance Cayley + Fourier), HYP-2170 (LRC↔unit-distance,
ℤ[ζ_6]), `the-five-platonic-tournaments.md` (A_5=icosahedron), `norm-layer-...-kissing-s702`,
S582 (LRC χ=2 interval vs χ=3 Paley), de Grey 2018, Barajas-Serra (circular chromatic LRC@7).
