---
id: THM-1990
title: "THE FIGURATE RECIPROCAL LADDER AND THE HARMONIC EDGE — reading every integer sequence as a SUB-SERIES of the harmonic series (sum 1/a_n) sorts the project's sequences by a single boundary at p=1, and the project's own dimensional ladder vertices -> arcs -> tetrahedra IS the figurate reciprocal ladder crossing that boundary. PROVEN (telescoping): sum_{n=k}^{N} 1/C(n,k) = k/(k-1)*(1 - 1/C(N,k-1)), so sum_{n>=k} 1/C(n,k) = k/(k-1) exactly, via 1/C(n,k) = k/(k-1)[1/C(n-1,k-1) - 1/C(n,k-1)]. THE LADDER: k=1 (vertices n, dim 1) = the HARMONIC series, DIVERGES = the edge; k=2 (arcs C(n,2) = triangular, dim 2) = 2 (the owner's value, sum 1/T_n = 2); k=3 (tetrahedral) = 3/2; k=4 = 4/3; ... -> 1 as k->inf. THE HARMONIC EDGE IS THE SAME p=1 BOUNDARY AS THM-1980: dimension 1 (the ground set / vertices) is exactly where reciprocals diverge AND where H's formula content collapses to one bit -- both are 'p=1', the marginal case, and the project's objects (arcs and up) all live on the convergent/formula-poor side. THE RECIPROCAL SIGNATURE sigma(a)=sum 1/a_n is a real-number fingerprint per sequence, VERIFIED to named constants: arcs 2, tetrahedral 3/2, var(lambda^2)=2C(n,3) -> 3/4, squares pi^2/6, factorial e, central binomial 4/3+2 pi sqrt3/27, Catalan 2+4 pi sqrt3/27, Fibonacci = reciprocal-Fibonacci constant, 2^n-1 = Erdos-Borwein, 2^n = 1; and the 2-ADIC / LACUNARY family gives THETA values: labeled tournaments sum 1/2^{C(n,2)} = 1.6416325606551539... (partial theta at q=1/2), switching classes sum 1/2^{C(n-1,2)} = 1 + that = 2.6416325606551539... (the +1 = the extra n=1 term, PROVEN exact), and the SIGNED pentagonal version is Euler's product prod(1-2^{-n})=0.288788...=(1/2;1/2)_inf (pentagonal number theorem, the repo's pentagonal thread). THE CONVERGENCE DICHOTOMY: sum 1/a_n diverges iff a_n grows at most linearly; EVERY combinatorial sequence in the corpus (degree >=2: arcs, tournaments, even graphs, tangent, ...) converges, only the linear ground-set sequences (vertices, odd numbers, the H-spectrum odds\\{7,21}) sit on the harmonic edge."
status: >
  The figurate ladder is PROVEN: the telescoping identity is exact (verified by exact rational
  partial sums, match=True all N,k), giving sum_{n>=k}1/C(n,k) = k/(k-1) for k>=2 and divergence at
  k=1.  The reciprocal-signature values are VERIFIED to 18-40 digits vs their closed forms (central
  binomial = 4/3+2pi sqrt3/27 EXACT to 18 digits; squares = pi^2/6; factorial = e; Mersenne =
  Erdos-Borwein; switching = 1 + tournament EXACT).  The theta/lacunary constants (tournament,
  switching) are genuine partial-theta values (no elementary closed form expected; a lacunary series).
  The census sums (A000568=3.8535..., A002854=3.0618..., score=3.9325..., tangent=1.5663...,
  secant=2.2171...) are convergent numerical fingerprints.  The HARMONIC-EDGE UNIFICATION with
  THM-1980 (both boundaries are p=1 = dimension 1 = the ground set) is a REFRAMING of two proven
  facts (figurate divergence at k=1; H's 2-adic depth -> 1).  Classical background: figurate
  reciprocal sums k/(k-1), Basel, Erdos-Borwein, reciprocal-Fibonacci (Andre-Jeannin irrational),
  pentagonal number theorem are known; the contribution is the UNIFIED classification of the corpus'
  sequences by their harmonic-subset signature and the ladder=dimension=p=1 synthesis.
source: kind-pasteur-2026-07-21-S128c144 (owner: reciprocal of a sequence = a subset of the harmonic numbers; sum 1/T_n = 2 while 1+1/2+..+1/5 > 2; study our sequences' reciprocal sums extensively and extend)
depends_on:
  - THM-1980    # the 2-adic edge of H (the p=1 / harmonic boundary this unifies with)
  - THM-1870    # the cycle-length edge (Hamiltonian length = the other p=1 boundary)
related: [THM-1885, THM-1880, THM-488]
external:
  - "Figurate reciprocal sums sum 1/C(n,k)=k/(k-1); Basel pi^2/6; Erdos-Borwein constant (sum 1/(2^n-1)); reciprocal-Fibonacci constant (Andre-Jeannin, irrational); central-binomial/Catalan sums (4/3+2pi sqrt3/27, 2+4pi sqrt3/27); Euler pentagonal number theorem / q-Pochhammer (1/2;1/2)_inf; Jacobi partial theta."
  - "OEIS: A000217 (triangular), A000292 (tetrahedral), A000568, A002854, A000571, A000182, A000364, A000108, A000045, A000225."
script: 04-computation/reciprocal_sums_of_our_sequences_kps_S128c144.py (+ .out)
---

# THM-1990 — the figurate reciprocal ladder and the harmonic edge

**The owner's frame.** A sequence `a_n ⊆ ℕ` turns its reciprocals `Σ 1/a_n` into a *sub-series of the
harmonic series* `Σ 1/m`. The full harmonic series diverges; a sub-series may converge, and its sum
is a fingerprint. The anchor: `1 + 1/2 + 1/3 + 1/4 + 1/5 > 2` already, yet `Σ 1/T_n = 2` *exactly*.

## The figurate ladder (proven)

**Telescoping identity:** `1/C(n,k) = k/(k-1) · [1/C(n-1,k-1) − 1/C(n,k-1)]`. Summing,

> **`Σ_{n=k}^{N} 1/C(n,k) = k/(k-1) · (1 − 1/C(N,k-1))`** (exact), hence
> **`Σ_{n≥k} 1/C(n,k) = k/(k-1)`** for `k ≥ 2`.

| dim `k` | sequence | `Σ 1/a_n` |
|---|---|---|
| **1** | vertices `n` | **∞ (harmonic — the edge)** |
| 2 | arcs `C(n,2)` = triangular | **2** ← the owner's value |
| 3 | tetrahedral `C(n,3)` | 3/2 |
| 4 | `C(n,4)` | 4/3 |
| k | `C(n,k)` | `k/(k-1)` |
| ∞ | — | → 1 |

The **arc-count** of `K_n` is `C(n,2)`, so the corpus' central "triangle" sequence (arcs, the
staircase cells) has reciprocal signature exactly **2**. Its telescoping `1/C(n,2) = 2(1/(n-1) − 1/n)`
*unfolds the arc-reciprocal mass into vertex-reciprocal differences* — the arc↔vertex (dim-2↔dim-1)
reduction **is** the telescoping, and it lands on `2 = 2·(1/1)`.

## The harmonic edge = p = 1 = dimension 1 (unifies with THM-1980)

The ladder **diverges only at `k = 1`** — the ground set (vertices) itself. Every object the project
builds *on top of* the vertices (arcs, tournaments, cycles, even graphs) has a **convergent**
reciprocal sum. So:

> **The dimension-1 ground set is exactly the harmonic boundary `p = 1`** — and it is *the same*
> boundary as THM-1980, where H's formula-expressible content collapses to a single bit ("H at
> `p=1`"), and as THM-1870, where cycle counts turn `#P` at the Hamiltonian length. Three independent
> lenses — reciprocal convergence, 2-adic depth, cycle length — all place the marginal case at the
> same `p=1` corner. The project's dimensional ladder `n → C(n,2) → C(n,3) → …` is the figurate
> reciprocal ladder *crossing* that edge, from divergence (2, the harmonic origin) into convergence.

## The reciprocal-signature table (verified)

| sequence | `Σ 1/a_n` | closed form |
|---|---|---|
| arcs / triangular | 2 | figurate `k=2` |
| tetrahedral | 3/2 | figurate `k=3` |
| `var(λ²)=2·C(n,3)` (THM-1880/1930) | 3/4 | `½·(3/2)` |
| squares `n²` | 1.6449… | `π²/6` (Basel) |
| factorial `n!` | 2.71828… | `e` |
| central binomial `C(2n,n)` | 1.73639… | `4/3 + 2π√3/27` (exact) |
| Catalan | 2.80613… | `2 + 4π√3/27` |
| Fibonacci | 3.3599… | reciprocal-Fibonacci const (irrational) |
| Mersenne `2ⁿ−1` | 1.60669… | Erdős–Borwein const |
| powers of 2 `2ⁿ` (Cayley–Dickson) | 1 | geometric |
| **labeled tournaments `2^{C(n,2)}`** | **1.6416325607…** | partial theta at `q=½` |
| **switching classes `2^{C(n−1,2)}`** | **2.6416325607…** | `= 1 + tournaments` (proven exact) |
| A000568 tournaments (unlabeled) | 3.8535… | super-exp fingerprint |
| A002854 even graphs `V(Eₙ)` | 3.0618… | super-exp fingerprint |
| A000571 score sequences | 3.9325… | fingerprint |
| A000182 tangent | 1.5663… | fingerprint |
| A000364 secant/Euler | 2.2171… | fingerprint |
| **vertices / odds / H-spectrum** | **∞** | harmonic edge |

## The 2-adic / theta family and the pentagonal thread

The lacunary `2^{C(n,2)}` sequences give **partial-theta** values at `q=½`. Two exact structural facts:

- **`Σ 1/2^{C(n−1,2)} = 1 + Σ 1/2^{C(n,2)}`** — the switching-class signature is the tournament
  signature plus 1 (the extra `n=1` term; verified exact). The `+1` is the lone singleton tournament.
- The **signed** pentagonal version is **Euler's product** `∏_{n≥1}(1−2^{−n}) = 0.2887880951… =
  (½;½)_∞` — the pentagonal number theorem `∏(1−qⁿ) = Σ(−1)^k q^{k(3k−1)/2}` evaluated at `q=½`. This
  ties the 2-adic / switching-class thread (and THM-488's pentagonal/η²⁴ hub) to partition theory and
  modular forms via the reciprocal-sum lens.

## Extensions (the "extend" ask)

1. **The reciprocal signature `σ(a) = Σ 1/a_n` as a sequence invariant.** It fingerprints each corpus
   sequence by a single real number; the figurate ones give the rational ladder `k/(k-1)`, the
   analytic ones give `π`/`e`/golden constants, the 2-adic ones give theta values.
2. **Telescoping = the `a`-monoid action (THM-1885).** `1/C(n,k)` unfolds into dimension-`(k−1)`
   reciprocal differences; the generator `a: n ↦ n+1` shifts the ladder, and the reduction
   Mode-A (`n → n−1`) is one telescoping step. Reciprocal sums are a representation of the `a/b`
   monoid on figurate dimensions.
3. **The convergence dichotomy.** `Σ 1/a_n` diverges **iff** `a_n` grows at most linearly. The corpus
   lives almost entirely on the convergent side (everything degree ≥2); only the linear ground-set
   sequences (vertices, odds, the H-spectrum) sit on the harmonic edge — a clean separation of
   "ground set" from "structure built on it."

## Named next

- **Identify the theta constant `1.6416325607…`** (`Σ 2^{−C(n,2)}`) in closed form / as a specific
  Jacobi `θ` value; likewise the census fingerprints (are A000568/A002854 sums related by the
  Burnside `÷n!`?).
- **The signed reciprocal sums** `Σ (−1)^n /a_n` for the corpus sequences — do the tournament / even-
  graph alternating sums hit partition/modular constants (as the pentagonal one does)?
- **`Σ 1/H(T)` over all tournaments** — the Rédei-count reciprocal sum, a bridge from this lens to
  THM-1980's `#P` object.
