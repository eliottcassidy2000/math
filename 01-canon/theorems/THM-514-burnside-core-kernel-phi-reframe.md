---
id: THM-514
title: The Burnside core-kernel φ-reframe — peeling the 1-tail isolates the difficulty of A000568 in a positive-definite GCD/Euler-φ quadratic form e(μ)=C(t,2)+½Σφ(d)M_d², with an exact add-a-part recurrence; the same compression is one member of a metagraph-enumerator family (tournaments, graphs, self-complementary, edge-colored)
status: PROVED (the φ/Smith-GCD identity for the Burnside edge-orbit exponent; the add-a-part recurrence). VERIFIED computationally: φ-reframe on all 1113 cores of mass ≤40 (0 failures), the recurrence on all cores of mass ≤35 (0 failures), full Burnside == (m,t)-compressed == A000568 through n=16 (integer, matches all known terms), the enumerator family T_n(2)=A000568 / T_n(1)·n!=A000246 / G_n(2)=A000088 through n=8.
source: kind-pasteur-2026-06-15-S7
depends_on:
  - THM-505   # OCF non-spectral defect dim = A000009(n)-3 = partitions into odd parts ≥3 (the SAME core family)
related:
  - THM-511   # converse-parity / the score (level-1) layer = the 1-tail (fixed points)
  - HYP-2538  # the θ-function / generating-function form of the core kernel
  - HYP-2539  # the P_n(x) edge-orbit triangle is a new sequence; the enumerator family is unified
  - OPEN-Q-103
  - T825
  - reflection: the-burnside-core-kernel-is-a-gcd-quadratic-form-and-a-metagraph-enumerator-family-kps
---

# THM-514 — the Burnside core-kernel φ-reframe

**Setup (the user's 1-tail peeling).** `A000568(n)` = #tournaments on `n` nodes
`= Σ_{λ⊢n, all parts odd} 2^{e(λ)}/z_λ`, where `e(λ)=Σ_i (λ_i−1)/2 + Σ_{i<j} gcd(λ_i,λ_j)`
(the number of free edge-orbits of a cycle-type-`λ` permutation; even-part types fix no
tournament) and `z_λ = ∏_k k^{m_k} m_k!`. Split `λ = μ ∪ 1^r` with `μ` the **core** (odd
parts `≥3`) and `1^r` the **tail**. Then `e(λ)=e(μ)+C(r,2)+r·ℓ(μ)`, `z_λ=z_μ·r!`, so

> `a(n) = Σ_{m,t} B[m,t] · 2^{C(n−m,2)+(n−m)t}/(n−m)!`,  `B[m,t]=Σ_{μ:|μ|=m,ℓ(μ)=t} 2^{e(μ)}/z_μ`.

The kernel `B[m,t]` is **`n`-independent**; the whole `n`-dependence is the explicit
factor. (Verified: full == compressed == A000568, `n≤16`, exact integers.)

## 1. The φ / GCD-matrix reframe (PROVED)

The residual difficulty — the cross-term `Σ_{i<j} gcd` that does not factor over `(m,t)` —
collapses via Smith's identity `gcd(p,q)=Σ_{d|p, d|q} φ(d)`. Writing `M_d :=` #parts of `μ`
divisible by `d`:

> **`e(μ) = C(t,2) + ½ · Σ_{d odd, d≥3} φ(d) · M_d²`**   (`t = ℓ(μ)` = #parts).

This is a **positive-definite quadratic form on the divisor-multiplicity lattice**
`(M_3, M_5, M_7, …)`, weighted by Euler's `φ` — a Smith-determinant / GCD-matrix structure,
hence a **theta-function exponent**: `2^{e(μ)} = 2^{C(t,2)} ∏_{d≥3} (2^{φ(d)/2})^{M_d²}`
(`φ(d)` even for `d≥3`, so the exponents are integers). VERIFIED: equals `e(μ)` on all
**1113 cores of mass ≤40, 0 failures**.

## 2. The exact add-a-part recurrence (PROVED)

Adding one odd part `p≥3` to a core `μ'`:

> **`Δe = (p−1)/2 + Σ_{d | p} φ(d) · M'_d`**   (sum over **all** odd divisors `d≥1`; the
> `d=1` term `φ(1)M'_1 = ℓ(μ')` = current #parts).

VERIFIED: reproduces `e` on all cores of mass `≤35`, 0 failures. **This pinpoints the
obstruction:** the increment needs the full **divisor profile** `{M'_d}`, not just `(m,t)` —
so a closed `(m,t)`-only recurrence cannot exist; the difficulty is *exactly* the GCD
cross-term, and the recurrence closes on the divisor-profile state.

## 3. The 1-tail = the trivial/score layer; the cores = the OCF non-spectral family

The peeled `1^r` are the **fixed points** of the automorphism — the analogue of the
score/level-1 "ranking" layer (THM-511): the leading term `m=t=0` gives the rigid
asymptotic `a(n) ~ 2^{C(n,2)}/n!`. The cores `μ` (**odd parts ≥3**) are exactly the family
that indexes the **OCF non-spectral defect** `dim = A000009(n)−3` (THM-505, disjoint-odd-cycle
packings). Same family, two roles: *automorphism core of the iso-count* and *cycle-packing
core of the OCF* — both "the irreducible part after removing the trivial layer."

## 4. The metagraph-enumerator family (the missing metric)

The compression is one member of a `S_n`-cycle-index family; the 1-tail peeling and the
φ-reframe apply to each (verified `n≤8`):

| enumerator | parts | base | `x=2` value | other anchor |
|---|---|---|---|---|
| `T_n(x)=Σ_{odd λ} x^{e}/z` (orientations) | odd | 2 | **A000568** (tournaments) | `T_n(1)·n!`=**A000246** (odd-cycle perms) |
| `G_n(x)=Σ_{all λ} x^{e_graph}/z`, `e_graph=Σ⌊λ_i/2⌋+Σgcd` | all | 2 | **A000088** (graphs) | colored `K_n` |
| `SC(n)` (self-complementary) | — | **4** | spine count | base-4 Burnside + 2^{#parts} fixed-vertex tax (survey) |

The integer polynomial **`P_n(x) = n!·T_n(x) = Σ_{odd-cycle σ} x^{#edge-orbits(σ)}`** is the
edge-orbit enumerator of odd-cycle permutations: `P_n(1)=A000246(n)`, `P_n(2)=n!·A000568(n)`.
Its coefficient triangle (`P_7: 720,504,280,70,1`; rows sum to A000246) appears to be a
**new sequence** (not located in OEIS this session).

## Scope / honesty

PROVED & VERIFIED: the φ-reframe, the add-a-part recurrence, the `(m,t)` compression
(== A000568, `n≤16`), the enumerator-family `x=2`/`x=1` anchors (`n≤8`). The φ-reframe does
**not** by itself yield a closed form for `a(n)` — A000568 remains open; the contribution is
that the difficulty is now an explicit GCD-quadratic-form theta sum over a divisor lattice
(HYP-2538) and that the same compression unifies a family of metagraph counts (HYP-2539).
SOURCED (survey + WebSearch): Burnside fixes `2^{c(g_A)/2}` tournaments for odd-order `g`;
the SC base-4 formula; A000009−3 = THM-505. Cross-links: THM-505, THM-511, the Γ-species
quotient view (arXiv:1204.1402), tournaments↔even-graphs equinumerous (arXiv:2204.01947).
