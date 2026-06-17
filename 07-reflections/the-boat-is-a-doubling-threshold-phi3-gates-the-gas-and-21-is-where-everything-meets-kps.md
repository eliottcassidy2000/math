# The boat dichotomy is a doubling threshold, Φ₃ gates both the ideal gas and the forbidden gaps, and 21 is where Fibonacci, triangular, the (2n−1)(4n−1) family, and the phantom volume all meet

**Source:** kind-pasteur-2026-06-16-S2. Dispatch: the deep structure of the Alcuin small-boat
vs big-boat dichotomy for conflict graphs; Fibonacci/triangular/square/prime and the repo's
prior bridges; the Heronian triangle paper (the 5 area=perimeter triangles ↔ 5 Platonic
solids); Heron's 4-factor and four-square recursion; equidecomposability (arXiv:1809.09936);
octal squares; Fib∩Tri; Fibonacci divisibility / Euler's criterion. All number-facts below
are VERIFIED (`04-computation/triangular_fibonacci_heron_boat_kps.py`); this is a synthesis
of mostly-existing repo threads plus a few new precise observations, scoped honestly.

## 1. The deep small/big-boat structure is a DOUBLING (dyadic) threshold

Csorba–Hurkens–Woeginger's cleanest large-boat criterion: the complete multipartite graph
`K_{n₁≤…≤nₖ}` is **large-boat (Alcuin=τ+1) ⟺ nₖ > 2·n_{k−1}** — the biggest part exceeds
**twice** the rest. So the boat `+1` is a **factor-of-2 / dyadic dominance threshold**: you
need the big boat exactly when one independent class *dominates* (more than doubles) its
nearest rival. This is the same `×2` the repo keeps meeting — the dyadic-valuation atom (the
three-twos), the pseudo-doubling `2−1/(n−2)→2=√2²`, the `4=2²` Boolean square. For the OCF
**conflict graph** the threshold degenerates to its sharpest form (THM-519): large-boat ⟺ `Ω`
edgeless ⟺ `H=3^{α₁}` (the ideal gas) — "one part dominates" becomes "every odd cycle is its
own isolated component." The boat tax is the dyadic dominance tax; tournament **dominance /
condensation** (THM-398 dodge) is the same axis.

## 2. Φ₃ gates BOTH the ideal-gas big-boat AND the forbidden gaps

The third cyclotomic polynomial `Φ₃(x)=x²+x+1` (the minimal polynomial of the 3-cycle's
eigenvalues `e^{±2πi/3}` — the Eisenstein angle) evaluated on the dyadic ladder:

> `Φ₃(1)=3` = the **big-boat / ideal-gas unit**: large-boat `H = 3^{α₁} = Φ₃(1)^{α₁}` (one
> factor of 3 per disjoint odd cycle, THM-519).
> `Φ₃(2)=7` = the **forbidden phantom** `H=7` (S599v: a scissors-congruence phantom volume,
> never realized; THM-029).
> `Φ₃(1)·Φ₃(2)=21` = the second forbidden gap `{7,21}` (S599v: `21=3·Φ₃(2)`).

So the **same cyclotomic Φ₃ that defines the 3-cycle governs both ends of the OCF lattice
gas** — its value at the ideal-gas fugacity-shadow `x=1` is the gas unit `3`, and at the
interacting `x=2` is the first phantom `7`. (Eisenstein `π/3` is also the dihedral angle of
the ideal hyperbolic tetrahedron whose volume `Cl₂(π/3)` is S599v's scissors invariant — the
3-cycle, Φ₃, the Eisenstein lattice, and the scissors density are one object.)

## 3. 21 is where four threads converge

`21` is, simultaneously and verifiably:
- the **4th Fibonacci∩Triangular** value (Fib∩Tri = {1,3,21,55}, exactly four);
- `T₆ = C(7,2)` (the 6th triangular = staircase area of the 7-vertex tournament);
- the `n=2` term of `(2n−1)(4n−1) = T_{4n−2}`, factoring as `3·7`;
- `Φ₃(1)·Φ₃(2) = 3·7` = the **forbidden H-gap** (3× the phantom 7).

And Fib∩Tri `= {1,3,21,55}` is **exactly the Fibonacci prefix of the triangular family
`(2n−1)(4n−1)=T_{4n−2}`** (`1,3,21,55,105,171,…`; Fibonacci-ness breaks at `n=4 → 105`).
The same quadruple is the **2×6-marble Burnside sequence** `[1,3,21,55,135,…]` (k=0..4 marbles
up to the grid's `Z₂×Z₂`; a Pólya count — *the same Burnside machinery as A000568*/THM-516 —
with `135` breaking both patterns at `k=4`). So `{1,3,21,55}` is a genuine *Pólya/triangular*
quadruple, and `21` is its tournament-meaningful member (forbidden H, `C(7,2)`).

## 4. Fibonacci divisibility IS the quadratic-residue structure that builds Paley

The Fibonacci rank-of-apparition obeys `α(p) | p − (5/p)` (the "p−1/p/p+1 rule"), where
`(5/p)` is the Legendre symbol — **Euler's criterion**, the QR of 5 (VERIFIED `p≤31`). This is
the *same* QR/Euler-criterion machinery the repo's **QR-resonance** (THM-305–308) and **Paley
tournaments** are built from, and THM-486 already proves the **involution modulus = Pisano
modulus**. So "Fibonacci mod p" and "Paley tournament on `𝔽_p`" are two faces of the Legendre
symbol: the Pisano period and the tournament's QR arc-set are governed by the same `(·/p)`.

## 5. Heron's 4-factor and the recursive 4 = 2²

`16·Area² = (a+b+c)(−a+b+c)(a−b+c)(a+b−c)`, i.e. with slacks `x=s−a, y=s−b, z=s−c` and
`s=x+y+z`: **`Area² = (x+y+z)·xyz`** — a perfect square written as a **product of 4 natural
numbers** (`s` the "glue" binding the three slacks). Dually, Fermat's polygonal / Lagrange:
every natural is a **sum of 4 squares**. The recurring `4 = 2²` is the same dyadic unit as the
boat threshold (§1) and the `n=4` Boolean square (THM-510). The `m=1` area=perimeter condition
is `xyz = 4(x+y+z)`, whose 5 solutions all have **inradius exactly 2** (`r = Area/s = P/(P/2)
= 2`). This `4`/`2`-fold structure is the "creative invariant" the dispatch asked for: the
square *is* a constrained 4-product, mirroring the 4-square sum — both pinned by the dyadic 2.
(Honest: this is a structural reframe, not a theorem about tournaments.)

## 6. The 5 triangles ↔ 5 Platonic solids: a count coincidence; the real 5 is elsewhere

There are exactly **5** integer triangles with Area=Perimeter — `(5,12,13),(6,8,10),(6,25,29),
(7,15,20),(9,10,17)`, all inradius 2 — and exactly **5** Platonic solids. Honest verdict:
**this is a count coincidence.** The triangles are 5 because the Diophantine `xyz=4(x+y+z)` has
5 positive solutions; the Platonic solids are 5 because `1/p+1/q>1/2` has 5 solutions. No
natural invariant of the triangles (areas `{30,24,60,42,36}`, the `r=2` inradius) maps to the
solids' `{V,E,F}` or symmetry orders. The **genuine** "5 ↔ tournament" structure is already in
the repo: `the-five-platonic-tournaments.md` (S20bh) maps the 5 solids to 5 tournament levels
(tetrahedron=K₄, cube=tiling cube, octahedron=arc space, icosahedron=G₅, dodecahedron=dual),
with `n=5` the Platonic/hyperbolic boundary and φ woven throughout. I keep the 5-triangle fact
as a verified aside and do **not** force a false correspondence.

## 7. Equidecomposability: the repo's own lens, and the source corrected

arXiv:1809.09936 is **not** about equidecomposability — it is Hirakawa–Matsumura, *A unique
pair of triangles*: there is a **unique** rational right triangle and rational isosceles
triangle with **equal perimeter and equal area** (via 2-descent + Coleman integration). That
"equal perimeter & equal area" matching rhymes with the `m=1` Heronian (area=perimeter) family
and with 2D **Bolyai–Gerwien** (equal-area polygons are scissors-congruent — no obstruction in
2D). The genuine *equidecomposability* content the dispatch points at is the repo's own
**S599v**: the forbidden `{7,21}` are **scissors-congruence phantom volumes**, with invariant
`Cl₂(π/3)` (the ideal-tetrahedron volume, the 3-cycle's Eisenstein angle) — exactly the Φ₃
object of §2. So the equidecomposability thread closes the loop back onto Φ₃.

## One-line synthesis

> The Alcuin `+1` is a dyadic *dominance* threshold (`nₖ>2n_{k−1}`); for the conflict graph it
> is the ideal-gas corner `H=3^{α₁}=Φ₃(1)^{α₁}`, and the *same* `Φ₃` that defines the 3-cycle
> hands back the forbidden phantom `Φ₃(2)=7` and the convergence number `21=Φ₃(1)Φ₃(2)` — the
> point where Fibonacci∩triangular, `T₆=C(7,2)`, the `(2n−1)(4n−1)` family, and the
> scissors-phantom gap all coincide; Fibonacci-mod-p is the Legendre symbol that also builds
> Paley; and the recurring dyadic `4=2²` ties Heron's 4-product, the four-square sum, and the
> boat threshold together.

## Status / honesty

VERIFIED: all number-facts (5 area=perimeter triangles & inradius 2; Fib∩Tri={1,3,21,55}=
Fibonacci prefix of `(2n−1)(4n−1)=T_{4n−2}`; 2×6-marble Burnside `[1,3,21,55,135,…]`;
`(2k+1)²=8T_k+1`; `α(p)|p−(5/p)`; `Φ₃(1)=3,Φ₃(2)=7,21=3·7`). SYNTHESIS (HYP-2553): the
boat=doubling-threshold and Φ₃-gates-both claims are precise re-readings of CHW + THM-519 +
S599v, not new theorems. The 5-triangle↔5-Platonic is a flagged COUNT COINCIDENCE. The
Fibonacci↔QR↔Paley bridge rests on THM-486 (proved) + Euler's criterion (classical). Cross-
links: THM-519 (boat=ideal gas), THM-517 (`H≤3^{α₁}`), THM-029/S599v (phantom `{7,21}`, Φ₃,
Cl₂(π/3)), THM-486 (Pisano), THM-516 (Burnside = the marble count), THM-510 (`4=2²` B₂ atom),
[[the-five-platonic-tournaments]], [[equidecomposability-the-scissors-congruence-of-the-cross-problem-bridge-s599]],
[[the-alcuin-boat-tax-is-the-ideal-gas-and-kuratowski-counts-overlapping-odd-cycles-kps]].
Sources: Csorba–Hurkens–Woeginger SIAM JDM 24(3) 2010; Hirakawa–Matsumura arXiv:1809.09936.
