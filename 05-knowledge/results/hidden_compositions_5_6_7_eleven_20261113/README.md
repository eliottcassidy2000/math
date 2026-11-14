# Hidden compositions: 5, 6, 7, the eleven of {2,3,11}, square-sum squares, and tournaments

mac-mini-2026-11-13-cuts (fourth prompt) · **Chan–Pak, HN and JC2 OPEN.** Owner prompt: relate 5, 6, 7 to spherical,
Euclidean and hyperbolic geometry; merge {2,3,11}, the golden ratio and the packing of 11 squares; see how all of
this relates to 4+7 = 5+6. Find which geometric "hidden composition" argument works for tournaments. Do more like
4+196 = 100+100 and 1–3–97–99, and bring in planar Jacobian themes.

Status: **FINITE-EXACT** (§2, §3, from complete censuses). **ELEMENTARY / PROVED** where stated. **ANALOGY**
labels the remaining dictionary entries.

## 1. The geometry (CITED, THM-4838 and THM-4592)

- **5, 6, 7.** The (2,3,n) triangle group is spherical, Euclidean or hyperbolic as 1/2 + 1/3 + 1/n is >, = or < 1.
  For n = 5 it is icosahedral, A5, with the golden ratio. For n = 6 it is the hexagonal wallpaper group. For n = 7 it
  is the Hurwitz group of the Klein quartic, with PSL2(7) of order 168 and the Fano plane.
- **{2,3,11}** is the owner's triangle in THM-4838. There, χ(2,3,p) = −B_{p−1} holds exactly for p ∈ {5, 7, 11},
  and the modular curves X(p) have genus 0, 3, 26.
- **11 = 4 + 7 = 5 + 6.** The reciprocal sums are 1/4 + 1/7 = 11/28 and 1/5 + 1/6 = 11/30. The numerator 11 is
  forced, and the denominators 28 and 30 differ by ((7−4)² − (6−5)²)/4 = 2.
- **Golden ratio and 11 squares.** THM-4592 reads Collatz period-5 parity words in base φ on Z[φ]/(φ^5 − 1). The
  eleven cells of the Fibonacci torus from the eleven-squares paper are 0, the 1/13 cycle and the −5 cycle. Their
  labels are the quadratic residues {1, 3, 4, 5, 9} mod 11 and their negatives, i.e. the arc set of the Paley
  tournament **P11**. P11 is the smallest 4-dichromatic tournament (Neumann-Lara), so order 11 is the first order
  that cannot be split into three transitive parts. F4 = 65537 also first occurs at order 11.

## 2. The tournament version of 4+7 = 5+6 (FINITE-EXACT)

Compose 11 vertices as strong components of sizes 4+7 or 5+6. The counts H(R_2)·H(C7) = 5·Sp(7) and H(C5)·H(C6) =
Sp(5)·Sp(6) take 67 and 48 values, and they share exactly **19**:

135, 165, 195, 225, 255, 275, 285, 325, 345, 375, 405, 435, 465, 495, 555, 585, 615, 645, 675.

**Mechanism (exact).** 5m lies in Sp(5)·Sp(6) only through the 5-multiples of the two spectra, namely 15 ∈ Sp(5) and
15, 25, 45 ∈ Sp(6). So the collisions are exactly 5m with m ∈ Sp(7) ∩ (3Sp(6) ∪ 3Sp(5) ∪ 5Sp(5) ∪ 9Sp(5)). For
example, 495 = 5·99 = 11·45 = 15·33.

**Correction (next session, [brick_moves](../brick_moves_odd_coordinates_20261113/README.md)):** the sentence that followed here was too strong. 5 of the 19 collisions also move an 11 or a 13; for example 165 = 5·33 = 11·15 trades a 5 for an 11. The exact law, checked at 11 and 12 vertices, is re-association by the primes {3, 5, 11, 13} of the order-5 strong spectrum. Many collisions do move a Busch brick, 3 = F0 or 5 = F1, from one strong component to the other. This is the
tournament form of a hidden composition: the vertex composition 4+7 = 5+6 is matched on H by re-associating
Fermat bricks between components. It is the same re-association that underlies Busch's floor and the rooted
core grammar of [ternary_core_grammar](../ternary_core_grammar_20261113/README.md).

## 3. Square-sum squares: 1–3–97–99 and its family (FINITE-EXACT, ELEMENTARY)

**1–3–97–99 is the only 4-cycle of distinct odd numbers below 100 with all four adjacent sums square.** The sums are
4, 100, 196 and 100, so 4 + 196 = 100 + 100 = 200.

Below 600 there are **168** such odd 4-cycles. The first ones are:

| 4-cycle | identity |
|---|---|
| 1–3–61–195 | 4 + 256 = 64 + 196 = 260 |
| 1–3–193–483 | 4 + 676 = 196 + 484 = 680 |
| 1–3–253–323 | 4 + 576 = 256 + 324 = 580 |
| 1–15–85–399 | 16 + 484 = 100 + 400 = 500 |
| 1–15–129–195 | 16 + 324 = 144 + 196 = 340 |

There are none among all integers below 40.

**The hidden composition (PROVED).** A square-sum 4-cycle (a, b, c, d) is exactly a pair of two-square
representations p² + q² = r² + s² of N = a + b + c + d, plus one free vertex weight. The vertices are recovered as
b = p² − a, c = r² − b and d = s² − a. In Z[i] the two representations are two ways of distributing Gaussian
conjugates over the same factorization. For N = 200 = 2³·5²:

- −2 + 14i = (1 + i)·2·(2 + i)² (norm 4 + 196)
- 10 + 10i = (1 + i)·2·(2 + i)(2 − i)

So moving one conjugate of 5 = (2 + i)(2 − i) is what turns 4 + 196 into 100 + 100.

**Tetrahedron and Klein.** View the four numbers as vertex weights on K4. Its three perfect matchings, the pairs of
opposite edges of a tetrahedron, all have the same total weight: 4 + 196 = 100 + 100 = 98 + 102, the last being the
diagonals. They are permuted by S4 → S3, whose kernel is the Klein group V4: rotations about opposite-edge axes. This
is exactly the diagonal-flip group of the four-vertex tournament atom, where a square plus two flip-able diagonals
gives H = 1, 3, 3, 5. The pattern 4+7 = 5+6 is the same three-matching identity with weights (2, 2, 3, 4).

## 4. Which hidden-composition argument transfers (summary)

| setting | visible compositions | hidden data | re-association move |
|---|---|---|---|
| Gaussian integers | p² + q² = r² + s² | Gaussian prime factorization | swap a conjugate of a split prime (5 = (2+i)(2−i)) |
| tetrahedron or K4 | three opposite-edge pairs | vertex weights | V4 = kernel of S4 → S3 |
| tournaments | SCC compositions (4,7) and (5,6) | strong-core counts | move a Busch brick 3 or 5 between components |
| planar JC (ANALOGY) | multidegree of an automorphism, deg = ∏ d_i | Jung–van der Kulk amalgam | re-bracket compositions; Abhyankar–Moh degree divisibility |

What transfers to tournaments is the multiplicative half: H is multiplicative over strong components, and the bricks
3 and 5 behave like split primes that can sit in either factor. The additive half (square sums) has no tournament
analogue beyond the OCF expansion H = Σ 2^|S|.

JC2 entries (CITED or ANALOGY):
- The repo's dictionary already matches the Jung–van der Kulk Bass–Serre tree with the Farey or Stern–Brocot tree
  (PROBLEM-LEDGER, Farey/Mandelbrot Collatz).
- Moh's theorem settles JC2 for degree ≤ 100. That is the "100" a counterexample must exceed, and the 100 of
  100 + 100.
- A Keller map that is a hidden composition of elementary maps is an automorphism. JC2 asserts every Keller map is
  one, much as every tournament is an ordinal composition of strong primes.

## 5. Next

- Test whether the 19 collisions persist as an exact rule at 12 vertices, with compositions 4+8, 5+7 and 6+6. The
  prediction is that collisions occur exactly where a brick can move.
- Recast the square-sum Hamiltonian-path lane (square-sum growth rules) on K4 vertex weights. Its 4-cycles are pairs
  of Gaussian two-square representations, which gives a factorization-based generator.
- Compare Klein-quartic (2,3,7) incidence with P_7/Fano and P11 with (2,3,11) as triangle-group colourings for the HN
  gates. The interior bound is now 10 (see the [interiors note](../hn_interiors_atoms_tilings_20261113/README.md)).
