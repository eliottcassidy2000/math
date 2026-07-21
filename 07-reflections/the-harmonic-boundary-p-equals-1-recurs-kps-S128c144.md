# The harmonic boundary p=1 recurs — three independent lenses land on the same edge

*kind-pasteur-2026-07-21-S128c144. Owner: reciprocal of a sequence = a subset of the harmonic
numbers; `Σ 1/T_n = 2` while `1+½+⅓+¼+⅕ > 2`; study our sequences extensively and extend. Companion
to THM-1990 (the figurate reciprocal ladder) and THM-1980 (the 2-adic edge of H).*

The owner has now pointed at the harmonic series from a second direction, and it closes a loop. The
same boundary — `p = 1`, the marginal case of the `p`-series `Σ 1/nᵖ` — shows up in **three**
independent parts of this project, each time as the exact dividing line between "expressible" and
"not":

| lens | the `p<1` / easy side | the `p=1` edge | the `p>1` / hard side |
|---|---|---|---|
| **reciprocal convergence** (THM-1990) | — | vertices `n` (harmonic, **diverges**) | arcs & up: converge (2, 3/2, …) |
| **2-adic formula depth** (THM-1980) | bits below parity: `#P` | H `mod 2` (Rédei): the **one** poly bit | full spectral ladder: poly |
| **cycle length** (THM-1870) | Hamiltonian `k=n`: `#P` | — | sub-Hamiltonian `k≤n−1`: spectral |

The alignment is not a pun. In all three, **dimension 1 — the ground set, the parity bit, the
marginal length — is where a formula runs out.** Everything the project builds *on top of* the
vertices (arcs, tournaments, cycles, even graphs) sits on the convergent / formula-rich side; the
bare ground set sits on the edge.

## Why `2`

`Σ 1/T_n = 2` is the `k=2` rung of the figurate ladder `k/(k-1)`, and `2` is where the ladder *first
lands* after the harmonic divergence at `k=1`. The owner's observation — that the harmonic partial
sum passes `2` by term 5 while the *triangular* reciprocals only *reach* `2` in the limit — is the
statement that **the arc-count sequence is the first sequence sparse enough to converge, and it
converges to exactly the boundary value the harmonic series races past.** Arcs are `K_n`'s
2-dimensional faces; the telescoping `1/C(n,2)=2(1/(n−1)−1/n)` is literally the map from arcs back to
vertices (dim 2 → dim 1), and it deposits its whole mass `2` at the first vertex. The staircase
triangle — the project's foundation — is the dim-2 figurate object, and `2` is its harmonic signature.

## The reciprocal signature as a fingerprint

Reading `σ(a) = Σ 1/a_n` across the corpus (THM-1990) turns each sequence into one real number, and
they sort into families that mirror the project's own structure:
- **figurate / staircase:** `2, 3/2, 4/3, …` (rational, the dimensional ladder), and `var(λ²)=2C(n,3) → ¾`.
- **analytic constants:** `π²/6` (squares), `e` (factorials), `4/3+2π√3/27` (central binomials), the
  golden reciprocal-Fibonacci constant, Erdős–Borwein (Mersenne) — the corpus quietly contains the
  headline reciprocal-series constants.
- **2-adic / theta:** labeled tournaments and switching classes give partial-theta values differing
  by exactly `1`, and the *signed* pentagonal version is Euler's `∏(1−2⁻ⁿ)` — the partition/modular
  side, linking to THM-488's pentagonal hub.

## The transcendence gradient

There is even a gradient of arithmetic difficulty *along the ladger*, matching the growth of the
sequence: polynomial growth → rational or classical-constant signatures (`2`, `π²/6`); exponential →
irrational analytic constants (golden, Erdős–Borwein); lacunary/2-adic → theta values (no elementary
closed form). The faster the sequence, the "harder" its harmonic signature — the reciprocal sum is a
thermometer for how far a sequence sits from the harmonic edge.

## What to carry forward

- The `p=1` boundary is the project's **invariant edge**: whenever a new object is introduced, ask
  which side of `p=1` it sits on (does `Σ 1/(its count)` converge? is it a ground-set or a built-on
  object?). The answer predicts whether its invariants are formula-expressible.
- **Signed reciprocal sums** (`Σ (−1)ⁿ/a_n`) are the unexplored half — the pentagonal case shows they
  can hit partition/modular constants; worth computing for tournaments and even graphs.
- **`Σ 1/H(T)`** (Rédei-count reciprocal over all tournaments) is the bridge from this lens directly
  onto THM-1980's `#P` object.

## Cross-links
THM-1990 (figurate ladder, signatures) · THM-1980 (2-adic edge of H) · THM-1870 (cycle-length edge) ·
THM-1885 (the `a`-monoid; telescoping = Mode-A) · THM-1880 (triangular numbers as `a/b` coefficients) ·
THM-488 (pentagonal / η²⁴ hub) · CLAUDE.md (the staircase triangle, Mode A/B).
