---
id: THM-1990
title: "THE FIGURATE RECIPROCAL LADDER AND HARMONIC EDGE (CORRECTED SUPPORT SEMANTICS)"
status: >
  PARTIALLY RETRACTED by MISTAKE-209 and repaired by THM-2000.  The exact
  simplex telescoping ladder remains PROVED.  The old linear-growth iff,
  switching-support +1, unnamed arithmetic-type claims, and global H-spectrum
  divergence are retracted.  Support semantics require deduplication and the
  exact near-linear boundary is Abel--Stieltjes/Bertrand.
source: kind-pasteur-2026-07-21-S128c144 (owner: reciprocal of a sequence = a subset of the harmonic numbers; sum 1/T_n = 2 while 1+1/2+..+1/5 > 2; study our sequences' reciprocal sums extensively and extend)
depends_on:
  - THM-1980    # the 2-adic edge of H (the p=1 / harmonic boundary this unifies with)
  - THM-1870    # the cycle-length edge (Hamiltonian length = the other p=1 boundary)
related: [THM-1885, THM-1880, THM-488, THM-2000, MISTAKE-209, THM-1370-h-spectrum-omits-7-21-all-n.md]
external:
  - "Figurate reciprocal sums sum 1/C(n,k)=k/(k-1); Basel pi^2/6; Erdos-Borwein constant (sum 1/(2^n-1)); reciprocal-Fibonacci constant (Andre-Jeannin, irrational); central-binomial/Catalan sums (4/3+2pi sqrt3/27, 2+4pi sqrt3/27); Euler pentagonal number theorem / q-Pochhammer (1/2;1/2)_inf; Jacobi partial theta."
  - "OEIS: A000217 (triangular), A000292 (tetrahedral), A000568, A002854, A000571, A000182, A000364, A000108, A000045, A000225."
script: 04-computation/reciprocal_sums_of_our_sequences_kps_S128c144.py (+ .out)
---

# THM-1990 — the figurate reciprocal ladder and the harmonic edge

> **CORRECTION (MISTAKE-209; THM-2000).**  The simplex telescoping theorem in
> this file is correct.  The broader atlas originally conflated a termwise
> **multiset** sum with the user's literal subset of harmonic numbers, asserted
> a false linear-growth iff, promoted the conjectural all-odd H-spectrum in
> `THM-1370-h-spectrum-omits-7-21-all-n.md` to a
> theorem, and missed Gauss's product for the triangular partial theta.  Under
> support semantics, repetitions are removed: labeled-tournament and switching
> values have the same support and the same mass; factorial, Fibonacci, and
> Catalan each lose one duplicated `1`.  Near-linear convergence is governed by
> Bertrand's iterated-log boundary, and H-spectrum divergence remains open.
> THM-2000 supplies the corrected support/Abel--Dini framework.

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

The **simplex ladder** diverges only at `k = 1` — the ground set (vertices)
itself — and every higher pure binomial column converges.  This is not a claim
about every sequence built on vertices: THM-1127's affine toothpick support,
for example, contains an arithmetic progression and diverges.  For the
simplex columns:

> **The dimension-1 ground set is exactly the harmonic boundary `p = 1`** — and it is *the same*
> boundary as THM-1980, where H's formula-expressible content collapses to a single bit ("H at
> `p=1`"), and as THM-1870, where cycle counts turn `#P` at the Hamiltonian length. Three independent
> lenses — reciprocal convergence, 2-adic depth, cycle length — all place the marginal case at the
> same `p=1` corner. The project's dimensional ladder `n → C(n,2) → C(n,3) → …` is the figurate
> reciprocal ladder *crossing* that edge, from divergence (`infinity` at dimension 1) into the
> convergent triangular mass `2` and beyond.

## The term-multiset signature table (verified numerically; support correction above)

| sequence | `Σ 1/a_n` | closed form |
|---|---|---|
| arcs / triangular | 2 | figurate `k=2` |
| tetrahedral | 3/2 | figurate `k=3` |
| `var(λ²)=2·C(n,3)` (THM-1880/1930) | 3/4 | `½·(3/2)` |
| squares `n²` | 1.6449… | `π²/6` (Basel) |
| factorial `n!`, indexed from `0` | 2.71828… termwise; **support `e-1`** | duplicate `0!=1!` |
| central binomial `C(2n,n)`, `n>=0` | 1.73639… | `4/3 + 2π√3/27` (exact) |
| Catalan, `n>=0` | 2.80613… termwise; **support is one less** | duplicated initial `1` |
| Fibonacci, `F_1=F_2=1` | 3.3599… termwise; **support is one less** | duplicated initial `1` |
| Mersenne `2ⁿ−1`, `n>=1` | 1.60669… | Erdős–Borwein const |
| powers of 2 `2ⁿ`, `n>=1` (Cayley–Dickson) | 1 | geometric |
| **labeled tournaments `2^{C(n,2)}`, `n>=1`** | **1.6416325607…** | partial theta at `q=½` |
| **switching classes `2^{C(n−1,2)}`, `n>=1`** | **2.6416325607… termwise** | support equals labeled-tournament support; `+1` is collision tax |
| A000568 tournaments (unlabeled) | historical prefix retracted | offset/error audit in MISTAKE-219 |
| A002854 even graphs `V(Eₙ)` | historical prefix only | no canonical tail value here |
| A000571 score sequences | historical prefix only | no canonical tail value here |
| A000182 tangent | historical prefix only | no canonical tail value here |
| A000364 secant/Euler | historical prefix only | no canonical tail value here |
| **vertices / odds** | **∞** | harmonic edge |
| **global H-spectrum** | **OPEN** | THM-1370 proves only two all-`n` omissions plus finite coverage |

## The 2-adic / theta family and the pentagonal thread

The lacunary `2^{C(n,2)}` sequences give Gauss's triangular-number theta value at `q=½`:

`sum_(r>=0)q^(r(r+1)/2)=(q^2;q^2)_infinity^2/(q;q)_infinity
=theta_2(0,sqrt(q))/(2q^(1/8))`.

Two exact structural facts:

- **Termwise**, `Σ 1/2^{C(n−1,2)} = 1 + Σ 1/2^{C(n,2)}`.  As harmonic
  **subsets**, the two ranges are identical and their support masses are equal;
  the termwise `+1` is precisely the repeated value `1`.
- The **signed** pentagonal version is **Euler's product** `∏_{n≥1}(1−2^{−n}) = 0.2887880951… =
  (½;½)_∞` — the pentagonal number theorem
  `∏(1−qⁿ) = Σ_(k in Z)(−1)^k q^{k(3k−1)/2}` evaluated at `q=½`. This
  ties the 2-adic / switching-class thread (and THM-488's pentagonal/η²⁴ hub) to partition theory and
  modular forms via the reciprocal-sum lens.

## Extensions (the "extend" ask)

1. **The support reciprocal feature
   `σ(A)=Σ_(m in support(a))1/m`.**  It is invariant under reindexing and
   repetition, unlike `Σ_n1/a_n`, but is not a complete fingerprint: THM-2000
   gives equal-mass supports with radically different block and achievement-set
   geometry.  The figurate examples give the rational ladder `k/(k-1)` and
   the 2-adic examples give theta values.
2. **Telescoping = the `a`-monoid action (THM-1885).** `1/C(n,k)` unfolds into dimension-`(k−1)`
   reciprocal differences; the generator `a: n ↦ n+1` shifts the ladder, and the reduction
   Mode-A (`n → n−1`) is one telescoping step. Reciprocal sums are a representation of the `a/b`
   monoid on figurate dimensions.
3. **The convergence boundary (corrected).**  For the strictly increasing
   support enumeration `b_n`, `b_n=O(n)` forces divergence and
   `b_n>=c n^(1+epsilon)` forces convergence, but neither is an iff.  Indexed
   repetitions invalidate these statements for `a_n` itself.  For example
   `b_n~n log n` is superlinear and still divergent.  THM-2000 gives the exact
   counting-function criterion and the Bertrand iterated-log scale.

## Named next

- The theta constant `1.6416325607…` is now identified by Gauss's product above.
  The arithmetic nature of the census support masses remains open; no
  transcendence claim is justified by fast convergence alone.
- **The signed reciprocal sums** `Σ (−1)^n /a_n` for the corpus sequences — do the tournament / even-
  graph alternating sums hit partition/modular constants (as the pentagonal one does)?
- **`Σ_(h in global H-value support)1/h`** — its convergence is genuinely open.
  By contrast, summing `1/H(T)` with one term per labeled tournament trivially
  diverges: `H(T)<=n!` while there are `2^C(n,2)` tournaments at order `n`.
