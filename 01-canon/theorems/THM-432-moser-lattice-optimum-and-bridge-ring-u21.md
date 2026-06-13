---
id: THM-432
name: the-u(21)-optimum-lives-in-the-Moser-lattice-(bridge-ring)-not-the-triangular-lattice
status: VERIFIED (cited published proof: Alexeev–Mixon–Parshall 2025) + construction exact-verified here
date: 2026-06-06
session: monad-explorer-2026-06-06-S4
depends_on:
  - THM-421   # the 3N common-neighbour floor (unit-distance lane, S709) — this CORRECTS its scope
  - HYP-2267  # triangular sqrt(7) DISK construction (suboptimal: gives 47 at n=21, not 57)
  - HYP-2262  # "bridge group between triangular and CM field" — the Moser lattice IS it (reportedly closed prematurely)
supersedes_framing_in:
  - 07-reflections/sc-complement-cosets-and-unit-distance-n21-s630.md  # "57 = 20+37 centered-hex" is the WRONG realization
---

# THM-432: the u(21) optimum lives in the MOSER LATTICE (the bridge ring), not the triangular lattice

> **Namespace note.** `u(21)=57` (the headline exact value) is recorded first-come in
> **THM-431** (monad-explorer-S710). THM-432 records the distinct STRUCTURAL content:
> the optimum is a *Moser-lattice* graph (triangular gives only 47), the exact `W₆⊕Δ`
> construction certifying `u(21) ≥ 57`, and the Moser lattice as the bridge ring that
> reopens HYP-2262. A 3-way THM-431 collision (this + two S710 files) was resolved by
> renumbering this later-arriving file to THM-432 (cf. MISTAKE-057 pattern).

## The dispatched campaign is already solved in the literature

The dispatched seed asked to "work toward a rigorous proof for N=21 of u(N) = the
maximum number of unit distances among N points in the plane." **This is now a
theorem in the published literature, and the value is exactly 57.**

`u(n) := max{ |E(G)| : G is a unit-distance graph on n vertices }`  (OEIS A186705).

### Statement [VERIFIED — published proof]

```
        u(21) = 57.
```

**Source.** B. Alexeev, D. G. Mixon, H. Parshall, *The Erdős unit distance problem
for small point sets*, arXiv:2412.11914v2 (12 Feb 2025), Theorem 1(a). They prove
the exact value of `u(n)` for **every** `n ≤ 21`:

| n    | 16 | 17 | 18 | 19 | 20 | 21 |
|------|----|----|----|----|----|----|
| u(n) | 41 | 43 | 46 | 50 | 54 | 57 |

Before this paper the best known bounds were `57 ≤ u(21) ≤ 68` (lower bound from
Schade 1993 / Engel et al. 2024; the paper closed the upper bound down to 57).

### Number of optimal graphs [VERIFIED — published enumeration]

There are **exactly 5** maximally dense unit-distance graphs on 21 vertices, up to
isomorphism (Alexeev–Mixon–Parshall, Theorem 1(c) / Table 1, last column; almost
all were previously found by Engel et al., arXiv:2406.15317). The counts of densest
graphs for the newly-settled range:

| n           | 16 | 17 | 18 | 19 | 20 | 21 |
|-------------|----|----|----|----|----|----|
| # densest   |  1 |  7 | 16 |  3 |  1 |  5 |

### Neighbouring values (only BOUNDS are known) [VERIFIED — published]

For `n ∈ {22,…,30}` only bounds are proven (Theorem 1(b)):

| n   | 22 | 23 | 24 | 25 | 26 | 27 | 28 | 29 | 30 |
|-----|----|----|----|----|----|----|----|----|----|
| ≥   | 60 | 64 | 68 | 72 | 76 | 81 | 85 | 89 | 93 |
| ≤   | 61 | 66 | 72 | 78 | 84 | 90 | 96 |103 |110 |

`n = 22` is the smallest `n` with `u(n) < ū(n)` (the F-free combinatorial bound
`ū(22)=62` overshoots; `u(22) ≤ 61`). The exact frontier of the *exact-value*
problem is therefore now at **n = 22**, not n = 21.

## The optimum is NOT in the triangular lattice — it is in the MOSER LATTICE

This is the structurally important correction for this project.

The **triangular (Eisenstein) lattice** `ℤ[ζ₆]` gives, by Harborth's theorem, only
`⌊3n − √(12n−3)⌋ = ⌊63 − √249⌋ = 47` unit distances at n=21. **47 < 57:** the
triangular lattice is suboptimal at n=21 by a full +10. (It is in fact beaten
already at n=9: Harborth gives 16, u(9)=18.) So the repo's recurring assumption
that the n=21 optimum is a triangular-lattice / √7-Eisenstein patch is **false**.

The densest small unit-distance graphs (including all 5 at n=21) live in the
**Moser lattice** (Engel et al. 2024, Def. 2.2):

```
   M_L = { a·1 + b·ω₁ + c·ω₃ + d·ω₁ω₃ : a,b,c,d ∈ ℤ } ≅ ℤ⁴,
   ω_t = exp( i · arccos(1 − 1/(2t)) ),
   ω₁ = ζ₆ = (1 + i√3)/2,          ω₁² − ω₁ + 1 = 0   (cos = 1/2, triangular generator)
   ω₃ = (5 + i√11)/6,              ω₃² − (5/3)ω₃ + 1 = 0   (cos = 5/6, the "Moser angle")
```

So `M_L = ℤ[ζ₆]` (the triangular lattice) **extended by the second unit direction
ω₃**, with `ω₃ = (5 + i√11)/6 ∈ ℚ(√−11)`. The whole lattice sits in the
biquadratic CM field `ℚ(√−3, √−11)`. The Moser lattice has **exactly 18 unit
vectors** (Engel et al., Thm 2.5: integer solutions of `|z|²=1` ∧ `ad=bc`), so a
Moser-lattice vertex has degree ≤ 18 (vs ≤ 6 for the triangular rosette). Every
known maximally dense UD graph is embeddable in `M_L`.

### Why 57: the Minkowski-sum structure

`57 = 3·12 + 7·3 = |Δ|·E(W₆) + |W₆|·E(Δ)` is the edge count of the **Minkowski sum
`W₆ ⊕ Δ`** of a **6-wheel** `W₆` (hub + hexagon, 7 vertices, 12 edges = 6 spokes +
6 rim) with a **unit triangle** `Δ` (3 vertices, 3 edges), `|W₆⊕Δ| = 7·3 = 21`.
Alexeev–Mixon–Parshall describe the n=21 optima as "a non-disjoint Minkowski sum of
three edges with a 6-wheel." The triangle must use a unit direction **off** the
triangular sublattice (i.e. an ω₃-direction) so the three wheel-copies stay
disjoint and acquire no extra coincident unit distances — which is exactly what the
extra Moser generator ω₃ provides and the pure triangular lattice cannot.

Exact-arithmetic construction + verification of a 57-edge 21-vertex Moser-lattice
graph: see `04-computation/unit_distance_moser_lattice_u21_monad_s4.py` (exact
ℚ(√3,√11) arithmetic; coefficients in basis (1,√3,√11,√33); no float decides any
adjacency). **Confirmed results** (`...s4.out`):
- M_L has **exactly 18** unit vectors (|coeff|≤4 exhaustive), all satisfying ad=bc:
  6 triangular (c=d=0) + **12 genuinely-Moser** (c or d ≠ 0). [Engel Thm 2.5 ✓]
- best triangular-lattice 21-pt patch (exact count): **U = 47** (= Harborth). ✗ optimum.
- of the 18 unit triangles Δ through 0: the **6** using only triangular directions
  collapse `W₆⊕Δ` to (|V|,U)=(12,24); the **12** using ≥1 transverse (ω₃/√11)
  direction give faithful **(|V|,U) = (21, 57)** — exactly u(21). Witness
  Δ = {0, (−2,1,2,−1), (−1,−1,1,1)} (both edges transverse); degree sequence
  **{deg 5 ×18, deg 8 ×3}**, ΣE = 57. This **certifies u(21) ≥ 57** by exact
  integer arithmetic.
- a naive greedy densifier in M_L stalls at the triangular blob (U=47): reaching 57
  *requires* the ω₃ direction and the Minkowski structure, it is not found by local
  growth. This is the concrete sense in which the Moser generator is essential.

## Consequences for the repo's unit-distance canon

1. **The seed's goal "prove u(21)" needs no new proof** — it is settled (=57). The
   open exact-value frontier is **n=22** (`60 ≤ u(22) ≤ 61`).
2. **The "20 + 37 centered-hex" decomposition** (reflection s630) had the right
   *number* 57 but the wrong *realization*: it is a Minkowski sum `W₆⊕Δ` in the
   Moser lattice, not a spine + centered-hexagonal bulk on a triangular/Eisenstein
   patch. The "37 = centered-hex H₃" was a numerology coincidence on the value 57,
   not the structure of the optimum.
3. **THM-421's "beats 3N" lane is consistent but not tight.** Since `u(n) < 3n` for
   all `n ≤ 21` (proven) and the best known constructions give `u(27) ≥ 81 = 3·27`
   and `u(28) ≥ 85 > 84 = 3·28`, the smallest N with a known config beating 3N is
   **N = 28**, not the `≤ 32` from THM-421(B) and far above its floor of 17. The
   floor 17 is a genuine *necessary* condition but very loose against the true ~28.
4. **HYP-2262 ("bridge group between triangular κ=6 and the CM field") should be
   REOPENED.** The Moser lattice `ℤ[ζ₆, ω₃]`, `ω₃=(5+i√11)/6`, IS such a bridge
   ring — triangular generator ⊕ a √−11 CM direction — and it is exactly where the
   density beyond the triangular rosette comes from. See HYP-2298.

## Honesty / scope

`u(21)=57` and the "5 densest graphs" count are VERIFIED by citing a published,
computer-assisted proof (F-free enumeration matching the lower bound + a custom
exact embeddability solver); they are not re-proved here. The Minkowski/Moser-
lattice **lower-bound** construction is exact-verified in this repo (integer
arithmetic over ℚ(√3,√11)). The matching **upper** bound u(21) ≤ 57 is the
computer-assisted part and is taken on the authority of the citation.

## References
- Alexeev, Mixon, Parshall, arXiv:2412.11914v2 (2025) — exact u(n), n≤21; 5 densest at 21.
- Engel, Hammond-Lee, Su, Varga, Zsámboki, arXiv:2406.15317 (2024) — Moser lattice; densest-known UD graphs to n≤100.
- Schade, *Exakte maximale Anzahlen gleicher Abstände*, TU Braunschweig thesis (1993) — exact to n≤14.
- Ágoston, Pálvölgyi, Studia Sci. Math. Hungar. 59 (2022) — u(15); improved constant factor.
- OEIS A186705.
