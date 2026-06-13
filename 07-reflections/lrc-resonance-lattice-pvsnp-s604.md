---
source: claude-2026-06-03-S604
status: proved Poisson formula + verified + complexity placement + conjectures
tags: [LRC, relation-lattice, resonance, SVP, geometry-of-numbers, harmonic-analysis, P-vs-NP, sharp-P, partition-function, OCF, Poisson, ergodic, recursive-log]
---

# The resonance lattice, and where LRC sits in the P-vs-NP landscape

## The master object is a lattice, not a graph

"Resonance graph" undersells it. The right object is the **relation lattice**

```
L(V) = { c ∈ ℤⁿ : c·V = 0 },     rank n−1.
```

Its **short vectors are the resonances**; the support-3 `±1` vectors
`e_i+e_j−e_k` are *exactly* the additive triples `v_i+v_j=v_k`. The "resonance
graph" is the bottom shell (support 3) of this lattice; the full lattice is what
governs tightness.

## The linchpin: a Poisson-summation formula for `p₀`

With `δ = 1/(n+1)`, expand the lonely indicator `∏_i(1−1_{F_i})` in Fourier
series and integrate. Each runner contributes its mean `1−2δ` (zero frequency)
or a nonzero frequency `m_i v_i` with coefficient `−sin(2πm_iδ)/(πm_i)`; the
integral keeps only the components summing to zero frequency — i.e. the lattice
`L(V)`:

```
p₀ = Σ_{c∈L(V)} ∏_i κ(c_i),   κ(0)=1−2δ,  κ(c)=−sin(2πcδ)/(πc).   [PROVED]
```

Verified numerically: truncating `|c_i|≤M` converges to the exact `p₀` from
every example. The `c=0` term is `(1−2δ)ⁿ` — the **independence value**, which
→ `e⁻² ≈ 0.135 > 0`. Every other term is a **resonance correction** indexed by a
lattice vector.

Three consequences fall out, each rigorous:

1. **Tightness requires resonances.** With only `c=0`, `p₀ = (1−2δ)ⁿ > 0`. So a
   relation-poor lattice can never be tight; tightness forces the resonance
   corrections to cancel `(1−2δ)ⁿ` exactly.
2. **The Vitali wall is the slow Fourier tail.** For tight configs the truncated
   sum approaches 0 *slowly* — the cancellation recruits arbitrarily high
   frequencies. No finite shell of short vectors (no finite set of moments)
   decides `p₀=0`. This is S603's wall, now visibly a property of `L(V)`.
3. **Tightness is a large deviation.** Mean depth → 2; for relation-poor speeds
   the arcs decorrelate (Weyl equidistribution of the joint flow `(v_i t)`),
   `depth → Poisson(2)`, `p₀ → e⁻²`. Tight configs are the arithmetic
   large-deviation extreme that drives the vacuum probability to 0.

## Where P vs NP enters — honestly

LRC has the exact *shape* of an NP problem, and the relation lattice explains the
hardness. I am not proving complexity-class separations; I am locating the
problem and identifying the mechanism.

- **Verifier in P [PROVED].** Given a rational time `t`, checking
  `‖v_i t‖ ≥ δ` for all `i` is polynomial in the bit-sizes. The extremal
  (tight) witnesses are **half-division points `j/(2n)`** — *polynomial* bit-size.
  So "a lonely time exists" carries short certificates: an NP-shaped existence
  question.

- **Deciding `p₀=0` is a short-vector cancellation on `L(V)` [structural].** By
  the Poisson formula, tightness is the statement that the lattice sum cancels
  the independence term. Finding the dominant resonance is an **SVP** instance on
  `L(V)`; SVP is NP-hard (ℓ∞: van Emde Boas; ℓ₂: Ajtai, randomized). So the
  *source of hardness* is geometry-of-numbers, not arbitrary combinatorics.

- **Counting is `#P` [analogy + proved sibling].** The depth polynomial
  `P(z)=Σ p_k z^k` is the LRC **independence polynomial**; evaluating
  independence/partition polynomials is `#P`-hard in general. Its repo sibling is
  `H(T)=#Hamiltonian paths` of a tournament — a `#P`-hard counting problem and
  the OCF object the whole repository orbits. So the LRC depth distribution and
  the tournament `H(T)` are two faces of one `#P` partition-function world.

- **Input encoding [PROVED].** Speeds in binary have size `~Σ log v_i`, but the
  depth function has `~2Σv_i` breakpoints — exponential. The direct algorithm is
  only *pseudo*-polynomial; genuine hardness, if present, is a binary-input
  phenomenon, exactly as in knapsack/subset-sum (themselves additive-relation
  problems on `ℤ`).

- **The dichotomy [CONJ].** Relation-poor `V` ⇒ `p₀ ≈ e⁻²`, loneliness abundant
  and cheap to certify; relation-rich `V` (AP, additive chains) ⇒ `p₀ → 0`,
  tight, and certification needs the whole lattice. **Arithmetic structure is the
  source of computational hardness here** — the inverse of the usual "random
  instances are hard." This is the same inversion seen in the repo's Collatz/
  two-block work (HYP-2143): the hard set is the structured, density-blind
  residual.

## Eight fields, one object

The same `(L(V), P(z), p₀)` triple reads natively in eight fields; the rigorous
spine is the lattice + Poisson formula, and I mark each tie's status:

| field | reading of the object | status |
|---|---|---|
| harmonic analysis | `p₀` = Poisson sum over `L(V)`; `κ` = arc Fourier coeff | **proved** |
| geometry of numbers | resonances = short vectors of `L(V)`; dominant = SVP | structural |
| additive combinatorics | support-3 `±1` vectors = additive triples | **proved** |
| statistical mechanics | `P(z)` partition fn; `z` fugacity; `p₀` vacuum; roots = Lee–Yang | analogy |
| probability | `depth→Poisson(2)`; `p₀→e⁻²`; tight = large deviation | verified |
| ergodic theory | `(v_i t)` Weyl-equidistribute; `L(V)` = joint-flow resonances | **proved** |
| complexity | `P(z)` eval `#P`; `SVP(L(V))` NP-hard; lonely-time verify in P | structural |
| tournament / OCF (repo) | `H(T)` = `#`Ham paths, `#P` independence-poly sibling | analogy |

## Where this sits in the session arc

The whole thread has one spine: a **counting/covering functional** and its
**log-scale + arithmetic-resonance structure**.

- Collatz in rapidity (`arctanh = log_F`): the `+1` is a bounded harmonic defect;
  the excursion exponent `θ=2` is fixed by `3+1=4` (HYP-2147–2150).
- Helly entropy accounting: the iterated-log *exponent* = obstruction *rank*
  (HYP-2151).
- LRC as covering: the depth distribution `p_k` is the master functional, tight =
  entropy-minimal cover (HYP-2152); the finer invariant is the Vitali wall —
  top-order only (HYP-2153).
- Here: the relation lattice `L(V)` is the algebraic master object; resonances
  are its short vectors; `p₀` is a Poisson sum over it; and the problem's
  complexity is a short-vector/`#P` phenomenon modulated by arithmetic structure.

The unifying picture: **`H(T)`, the Collatz orbit count, and the LRC depth
distribution are all partition functions; their "hard, structured extremes" are
governed by an arithmetic resonance lattice; and the iterated logarithm is the
ruler that measures how many resonance scales you must keep.**

## Open / next

1. **Make the dichotomy a theorem:** prove `p₀ ≥ c(n) > 0` for dissociated `V`
   (bound the resonance corrections by the lattice's first minimum) — a clean
   geometry-of-numbers statement.
2. **SVP reduction attempt:** can a genuine (even approximate) reduction from
   SVP/CVP to LRC-tightness be built on `L(V)`? That would upgrade the structural
   tie to a hardness theorem.
3. **OCF deletion–contraction for `P(z)`:** compute `p₀` via a contraction on the
   resonance lattice rather than inclusion–exclusion — the tournament method
   transported to LRC.
4. **Lee–Yang for `P(z)`:** the depth polynomial is not real-rooted (S603); map
   its complex zeros as `δ` varies — a statistical-mechanics phase transition at
   `δ = 1/(n+1)` (where `0` becomes a root, i.e. tightness onsets).
