# Odd zeta values and the Petersen graph: the ζ(5), ζ(7), ζ(9) preprint fails at an impossible arithmetic gain; what irrationality proofs are made of; Petersen inside M̄₀,₅, Brown's dinner parties, and Ziegler's question

**opus S15, 2026-10-01.**

- Script: `04-computation/experiments/odd_zeta_parallels_20261001.py` (+ `.out`, ALL CHECKS PASSED).
- Companion: the twenty-first Collatz note
  [`collatz_golden_holonomy_20261001.md`](collatz_golden_holonomy_20261001.md) (holonomy mirror).

**Status.**
- Audit of arXiv:2609.22316: its proof is INVALID. The decisive error is elementary, is checked here, and is
  confirmed by an independent numerical audit (§2).
- Structural parallels: typed (DICTIONARY / STRUCTURAL); facts about M̄₀,ₙ are standard (Keel).
- Ziegler's question (is the product of two Petersen graphs polytopal?): OPEN. Nothing is claimed beyond the
  literature except one elementary remark.
- No new irrationality result.
- Independent audit of this note: DONE (2026-10-01, blind).
  - One claim was UNSOUND and is corrected: the dinner-table counts are not Brown's convergent configurations
    from `N = 8` on, and the OEIS number was wrong.
  - The other wording corrections are applied (MISTAKE-556), and the record is in §6.

## 0. The owner's prompts

1. "Is the cartesian graph product of two Petersen graphs polytopal? If so, does the polytope have dimension 4 or
   5?"
2. "Pursue proofs related to the irrationality of ζ(2n+1) by drawing impressive surprising parallels to the deep
   structure of other fields of math, think deeply and freely about what things represent:
   arXiv:2609.22316v1."

## 1. Ziegler's question: open in dimensions 4 and 5

`P□P` has 100 vertices and 300 edges and is 6-regular. Ziegler asked whether it is polytopal (CRM 2009).
Pfeifle–Pilaud–Santos (Israel J. Math. 192 (2012)) settled what can be settled locally:
- **Dimension ≤ 6.** A d-polytope's graph has minimum degree ≥ d.
- **Not 6.** A 6-dimensional realization would be simple. Their Theorem 2.3 says a product of graphs is simply
  polytopal iff its factors are, and the Petersen graph is not polytopal. They settle dimension 6 via Corollary
  1.16(iii): a simply polytopal graph cannot contain an induced Petersen graph.
- **Not ≤ 3.** The graph is non-planar (Steinitz).
- **Dimensions 4 and 5: OPEN.** Their own words: "we have no answer for dimensions 4 and 5". The Wikipedia
  article "Graph of a polytope" says it was still unknown as of July 2025.
- **Local combinatorics cannot settle it**, as PPS remark after their Proposition 2.18. `P□P` is the graph of a
  cellular decomposition of the closed 4-manifold `RP² × RP²` (36 products of two pentagons; Petersen = the
  hemi-dodecahedron on `RP²`).

**Answer to the owner: unknown.** If `P□P` is polytopal, the polytope has dimension 4 or 5, and nobody knows
which, or whether either occurs.

*One elementary remark (not in PPS, easy).* In a 4-dimensional realization:
- The vertex figure at `(u,v)` is a 3-polytope on its 6 neighbours: 3 horizontal and 3 vertical.
- Suppose every 2-face through `(u,v)` had one horizontal and one vertical edge at `(u,v)` ("mixed at v"). Then the
  vertex figure's graph, which has minimum degree 3, would be `K₃,₃`, which is non-planar.
- So **at every vertex some 2-face contains two edges of the same Petersen fibre** (a fibre polygon, or a
  "bent" induced cycle such as a 7-cycle running through an adjacent fibre).
- **At most 8 of the 2-faces at each vertex are mixed at v**, in particular at most 8 of the 9 rectangles (4-cycles)
  through it. A planar bipartite graph on 6 vertices has ≤ 8 edges, and the bound is sharp among the 7 polyhedral
  graphs on 6 vertices.

This narrows the search but decides nothing.

## 2. arXiv:2609.22316 ("At least one of the three numbers ζ(5), ζ(7), ζ(9) is irrational"): the proof fails

The preprint (6 pages, math.GM) builds Zudilin-type very-well-poised linear forms with `h₀ = 3n+2`,
`h₁ = … = h₁₁ = n+1`:

`R(t) = n!² (3n+2+2t)(t+1)_n³(t+2n+2)_n³ / (t+n+1)_{n+1}⁸`, `F_n = ½ Σ_{t≥−n} R''(t)`,

and claims `C₀' = 5.7535 > C₂' = 5.2769`.

**The construction of the linear forms is standard.**
- `R = O(t^{−2n−7})`, so the residue sum `Σ_k B_{1,k} = 0` and the ζ(3) coefficient vanishes.
- The symmetry `R(−t−h₀) = −R(t)` kills the even zeta values.
- So `F_n ∈ Q + Qζ(5) + Qζ(7) + Qζ(9)`.

**The arithmetic is not.** Its eq. (15) exponent, for these parameters, works out as follows.

1. **The exponent is exact.** `v_p = 3` if `3(n mod p) ≥ 2p − 1`, and `v_p = 0` otherwise (Check A, all 4341 pairs
   `(n, p)` tested with `√h₀ < p ≤ n`).
   - These are Kummer carries: a carry-free split of the binomial-type coefficients exists iff some residue `r`
     has `2m − p < r < p − m`, where `m = n mod p`.
   - So `φ(x) = 3·1_{[2/3,1)}(frac x)`.
2. **The gain is 0.72306**, confirmed by the direct prime sum at `n = 3·10⁶` (0.7207):

   `lim (1/n) log Φ_n = ∫₀¹ φ(y) ψ'(1+y) dy = ∫φ dψ − ∫φ dy/y² = 3(ψ(1) − ψ(2/3)) − 3/2 = 2.22306 − 1.5`.

3. **The preprint subtracts the second integral instead of adding it back.** Its `C₂' = 9 − 2.22306 − 1.5 =
   5.27694374` reproduces its stated value to 8 digits. It evaluates its own (41) without the parentheses that
   appear in its (38). The correct value, keeping its own denominator exponent 9,
   is `C₂' = 9 − 0.72306 = 8.27694`. With the exponent that Zudilin's lemma actually supports (10), it is 9.27694.
4. **Even the trivial bound kills it.** `v_p ≤ 3` and `Σ_{p≤n} log p ~ n`, so the gain is at most 3 and
   `C₂' ≥ 6 > C₀' = 5.7535`.

   The preprint's implied gain of 3.7231 exceeds what any exponent pattern bounded by 3 can give.

**Verdict: INVALID.** The criterion's bound `e^{(C₂'−C₀')n}` has a positive exponent: 2.52 with the preprint's
denominator exponent 9, and 3.52 with the correct exponent 10. So it cannot force the integer linear forms to 0.
Zudilin's theorem (one of ζ(5), ζ(7), ζ(9), ζ(11) is irrational, 2001) is not improved.

**Independent numerical audit (subagent, own code; record in the session scratchpad `s21/zaudit/AUDIT.md`).**

*Sources.*
- The preprint's construction is Zudilin's *Arithmetic of linear forms involving odd zeta values* (J. Théor.
  Nombres Bordeaux 16 (2004); arXiv:math/0206176), Section 8, with `q = 13` replaced by `q = 11`. The citations'
  page numbers do not exist in the works cited.
- Zudilin's Lemma 11 gives `lim (1/n) log Φ_n = ∫φ dψ − ∫₀^{1/m} φ dx/x²`. His Proposition 5 gives
  `C₂ = r m₁ + m₂ + … − (∫φ dψ − ∫φ dx/x²)`.

*The sign.* Read with the subtraction, the formula reproduces Zudilin's published `C₂ = 226.24944266…` for his own
example (`q = 13`, `η₀ = 91`). Read with the preprint's flip, it gives 127.25, which is wrong.

*Confirmed.*
- `C₀' = 5.75349395301` (claim (c) is true). The saddle point `τ₀ = 2.86852453 + 0.11960091 i` is reproduced, and
  regressions on the exact `F_n` up to `n = 450` give 5.7535.
- The linear forms are exact: `A₃ = A₄ = A₆ = A₈ = A₁₀ = 0` for all `n = 3..300`, and they agree with direct
  summation to 45–85 digits (claim (a) is true).
- `(1/n) log Φ_n = 0.6282, 0.7178, 0.7199, 0.7202, 0.7216, 0.7226` at `n = 10³…10⁸`, converging to 0.72306.

*Further errors.*
- The integrality lemma (claim (b): `2D_n⁹Φ_n⁻¹ A_i ∈ Z`) is false. It fails for 280 of the 298 values
  `n = 3..300`, and the first counterexample is `n = 13, p = 7`. The step that creates the extra power silently
  changes an exponent `−(11−j)` into `−(9−j)`. The correct form, `D_n^{10} Φ_n⁻¹`, is Zudilin's Lemma 19.
- The denominator exponent is 10, not 9: the paper drops `m₈`.
- `φ` is never defined in the preprint.

*Consequence.* With the true least common denominators `d_n`, `ln(d_n|F_n|)` grows at slope `+2.73` per unit `n`.
The integer linear forms explode instead of tending to 0: `ln(d_n|F_n|) = 1231.7` at `n = 450`.

## 3. What the objects represent

### 3.1 An irrationality proof is an adelic inequality

For a nonzero rational `q`, `∏_v |q|_v = 1`. Apéry-type proofs build `q_n = d_n L_n` that is
- **integral at every prime** (the denominators `d_n`, rate `C₂`), yet
- **tiny at the real place** (`|L_n| ≈ e^{−C₀ n}`).

`C₀ > C₂` forces a nonzero integer below 1.

The Collatz problem in this repo has the same two-place shape:
- `7·Φ_T(n) ∈ Z` is a 2-adic integrality;
- whether `n` reaches 1 is an archimedean statement about size;
- Baker's method for cycles is the same inequality in disguise.

Irrationality proofs win the two-place tug-of-war with an explicit construction. Collatz has no construction to
offer. (STRUCTURAL ANALOGY.)

### 3.2 The arithmetic gain is carry statistics

`Φ_n` collects the primes that must divide every coefficient. By Kummer's theorem, `v_p` of a binomial counts the
carries when adding in base `p`. Zudilin's `φ(x)` is the minimal number of carries over all admissible splits, and
here carries are forced exactly when `frac(n/p) ≳ 2/3` (§2).

Collatz is carry statistics in base 2: `3n = n + 2n`, "11 as a bit shift with memory". The twenty-first note
shows that reading the parity word in base φ absorbs the carry into the golden normal form `11 → 100`.
- In irrationality proofs, forced carries are the resource.
- In Collatz, unforced carries are the obstacle.
(STRUCTURAL ANALOGY.)

### 3.3 Reflections with a centre at a half-integer

The very-well-poised symmetry `t ↦ −t − h₀` kills ζ(even), the zeta values whose irrationality is free (`π^{2k}`).

Collatz's sheet reflection `x ↦ −1 − x` (nineteenth note) pairs the `3n+1` and `3n−1` sheets.

Both are involutions that split a problem into a free part and a hard part. The Collatz centre is `−1/2`; the
well-poised centre `−h₀/2 = −(3n+2)/2` is half-integral when `n` is odd.
(DICTIONARY.)

### 3.4 The Petersen graph is the dual graph of the boundary of the moduli space of ζ(2)

Brown ("Irrationality proofs for zeta values, moduli spaces and dinner parties", 2016) rewrites Apéry-type
integrals as **cellular integrals** on `M₀,N`. A cell is a dihedral order `δ` of `N` points on `P¹`. A form is
another dihedral order `δ'`. The integral converges iff no set of `k` elements (`2 ≤ k ≤ N−2`) is consecutive in
both orders. For `N ≤ 7` this is the dinner-table condition: nobody sits next to a former neighbour.
- `N = 5` gives ζ(2), `N = 6` gives Apéry's ζ(3), and ζ(5) lives at `N = 8`.
- **The boundary of `M̄₀,₅` is the Petersen graph.** Its 10 boundary lines `D_ij` meet iff `{i,j}` and `{k,l}`
  are disjoint, so their intersection graph is `KG(5,2)` = Petersen (the 10 lines on the quintic del Pezzo
  surface).
- **Brown's ζ(2) configuration is the Petersen graph split in two** (Check C). For `δ = (12345)` and its unique
  convergent partner `δ' = (13524)`:
  - the consecutive pairs of `δ` (the boundary of the integration cell) form the outer pentagon;
  - those of `δ'` (the poles of the form) form the inner pentagram;
  - the two are joined by the perfect matching of spokes.

  Convergence means exactly that these two Petersen pentagons are disjoint.
- **ζ(5)'s moduli space contains `M̄₀,₅ × M̄₀,₅`.** By Keel, each of the 35 boundary divisors `D_S ⊂ M̄₀,₈` with
  `|S| = 4` is isomorphic to `M̄₀,₅ × M̄₀,₅` (Check D). Its 100 corner surfaces `D_ij × D_kl` meet along curves
  exactly when one coordinate agrees and the other pair of lines meets. So **the curve-adjacency graph of these
  corners is the Cartesian product of two Petersen graphs**, Ziegler's graph.

  (DICTIONARY: correct, but it says nothing about polytopality.)

### 3.5 Dinner parties and Brown's configurations

**Dinner-table rearrangements** are the Hamiltonian cycles of the complement of the `N`-cycle: `1, 3, 23, 177, 1553`
for `N = 5..9` (OEIS A002816), forming `1, 1, 5, 19, 112` classes up to the dihedral symmetry of `δ`.

**Brown's convergent configurations** satisfy his stronger condition. There are `1, 3, 23, 169, 1463` of them,
forming `1, 1, 5, 17, 105` classes (Check B, reproducing Brown's `C_N`). The two counts agree up to `N = 7` and
differ from `N = 8` on.
- Brown's 17 configurations at `N = 8` give 13 families of linear forms in 1, ζ(3), ζ(5) (his §10.2.4).
- The tournament notes of this session study the same objects from the other side: Hamiltonian paths that do or
  do not close, and the cycles they avoid. (DICTIONARY.)

### 3.6 Holonomy: Apéry and Collatz are mirror images

| | holonomic? | integral? | where the work is |
|---|---|---|---|
| Apéry | given (a Picard–Fuchs equation; Beukers' modular parametrisation) | needs denominators `d_n` | arithmetic and decay |
| Collatz | is the question: no divergent orbit ⟺ every parity series `F_n` is D-finite (Szegő; twenty-first note, Theorem G3) | free (coefficients 0/1) | holonomy |

Calegari–Dimitrov–Tang's 2024 arithmetic holonomy method bounds the holonomy rank (the dimension over `Q(x)`) of
integral power series with prescribed analytic continuation. Rationality is its rank-one (Borel–Pólya) case. For
Collatz, the analytic continuation of `F_n` is what is missing.

## 4. What a proof of "one of ζ(5), ζ(7), ζ(9) is irrational" would need

In this family the gap is `C₂' − C₀' = 8.277 − 5.753 = 2.52`.
- Raising `η₀` improves the decay but enlarges the denominators.
- Zudilin's optimisation over very-well-poised families stops at four zeta values. The best known result remains
  his 2001 theorem.
- A genuinely new input is needed: a larger arithmetic gain (more forced carries) or a better analytic
  construction (for example Brown's `N = 8` cellular integrals, or the arithmetic holonomy method).

**Open.**

## 5. Directions

- **D78.** Start from Brown's 17 convergent configurations at `N = 8` (13 families of linear forms, his §10.2.4)
  and compute `C₀` and `C₂` for each family.
- **D79.** Ziegler's question in dimension 4: enumerate the possible vertex figures (3-polytopes on 6 vertices
  with at least one fibre angle) and the compatible 2-face systems, as a SAT problem. A first step toward a
  non-realizability proof, or toward a construction.

## 6. Audit record (2026-10-01)

A blind auditor subagent used its own code and the texts of Brown (arXiv:1412.6508), PPS (arXiv:1009.1499) and the
preprint. Record: session scratchpad `audit21/AUDIT.md`.
- **Z1, Z2** (the exponent and the gain): SOUND, checked for all 14,097 pairs with `n ≤ 400`. Exact expansions
  confirm the gain is genuine.
- **P1, P2**: SOUND. The audit derived the 100 corners and their 300 curve adjacencies from stable trees and
  found them isomorphic to `P□P`.
- **P3**: UNSOUND as first written (dinner-table counts presented as Brown's configurations, wrong OEIS number).
  Corrected in §3.5 and Check B.
- **P4**: SOUND WITH CORRECTION ("mixed at v" instead of "rectangle"; manifold; PPS's Theorem 2.3 quoted as
  stated).

## 7. Reproduction

```bash
python 04-computation/experiments/odd_zeta_parallels_20261001.py
```
