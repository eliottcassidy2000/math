# The harmonic boundary p=1 recurs — three independent lenses land on the same edge

*kind-pasteur-2026-07-21-S128c144. Owner: reciprocal of a sequence = a subset of the harmonic
numbers; `Σ 1/T_n = 2` while `1+½+⅓+¼+⅕ > 2`; study our sequences extensively and extend. Companion
to THM-1990 (the figurate reciprocal ladder) and THM-1980 (the 2-adic edge of H).*

> **Correction (MISTAKE-209; THM-2000/2005).** The simplex ladder below is
> exact, but the scalar reciprocal sum is not a fingerprint and does not
> predict computational complexity.  Literal supports deduplicate repeated
> values, so labeled-tournament and switching rows have equal support mass.
> The arithmetic “gradient” and the alignment with formula/`#P` boundaries
> are heuristic only; the rigorous analytic edge is the full
> Abel--Dirichlet/Bertrand profile.

The original session noticed a useful analogy among three different
boundaries.  Only the first column below is literally an analytic `p=1`
statement; the other two are complexity/length filtrations and must not be
identified with it as theorems.

| lens | the `p<1` / easy side | the `p=1` edge | the `p>1` / hard side |
|---|---|---|---|
| **simplex reciprocal convergence** (THM-1990) | — | vertices `n` (harmonic, **diverges**) | binomial faces `C(n,k)`, `k>=2`, converge |
| **2-adic formula depth** (THM-1980) | bits below parity: `#P` | H `mod 2` (Rédei): the **one** poly bit | full spectral ladder: poly |
| **cycle length** (THM-1870) | Hamiltonian `k=n`: `#P` | — | sub-Hamiltonian `k≤n−1`: spectral |

The shared language is mnemonic, not predictive.  Outside the simplex ladder,
derived objects may lie on either side of the Abel--Dini boundary, and
reciprocal convergence says nothing by itself about formula complexity.

## Why `2`

`Σ 1/T_n = 2` is the `k=2` rung of the figurate ladder `k/(k-1)`, and `2` is where the ladder *first
lands* after the harmonic divergence at `k=1`. The owner's observation — that the harmonic partial
sum passes `2` by term 5 while the *triangular* reciprocals only *reach* `2` in the limit — is the
statement that **the arc-count row is the first convergent rung of this
particular simplex ladder.** Arcs are `K_n`'s
2-dimensional faces; the telescoping `1/C(n,2)=2(1/(n−1)−1/n)` is literally the map from arcs back to
vertices (dim 2 → dim 1), and it deposits its whole mass `2` at the first vertex. The staircase
triangle — the project's foundation — is the dim-2 figurate object, and `2` is its harmonic signature.

## The reciprocal signature as a lossy coordinate

Reading `D_A(1)=Σ_(m in support(a))1/m` across the corpus turns each support
into one real coordinate.  Some useful named rows are:
- **figurate / staircase:** `2, 3/2, 4/3, …` (rational, the dimensional ladder), and `var(λ²)=2C(n,3) → ¾`.
- **analytic constants:** `π²/6` (squares), `e` (factorials), `4/3+2π√3/27` (central binomials), the
  golden reciprocal-Fibonacci constant, Erdős–Borwein (Mersenne) — the corpus quietly contains the
  headline reciprocal-series constants.
- **2-adic / theta:** labeled tournaments and switching classes have the same
  support and hence the same theta mass; their **indexed** sums differ by the
  collision tax `1`.  The signed pentagonal version is Euler's
  `∏(1−2⁻ⁿ)`, linking to THM-488's pentagonal hub.

## The profile hierarchy

Growth controls convergence regions, not arithmetic type.  The faithful
hierarchy is support `A`, counting function `A(x)`, block occupancies, full
Dirichlet profile `D_A(s)`, and only then the scalar slice `D_A(1)`.  Distinct
supports can share both mass and abscissa, so no arithmetic or complexity
classification should be inferred from one decimal.

## What to carry forward

- For every new sequence, first fix the support/indexing convention, then ask
  for its counting law, abscissa, boundary behavior, and tail certificate.
- **Signed reciprocal sums** (`Σ (−1)ⁿ/a_n`) are the unexplored half — the pentagonal case shows they
  can hit partition/modular constants; worth computing for tournaments and even graphs.
- Separate the support mass `Σ_(h achieved)1/h` from multiplicity-weighted
  sums over labeled tournaments or isomorphism classes when studying Rédei
  counts; these are three different objects.

## Cross-links
THM-1990 (figurate ladder, signatures) · THM-1980 (2-adic edge of H) · THM-1870 (cycle-length edge) ·
THM-1885 (the `a`-monoid; telescoping = Mode-A) · THM-1880 (triangular numbers as `a/b` coefficients) ·
THM-488 (pentagonal / η²⁴ hub) · CLAUDE.md (the staircase triangle, Mode A/B).
