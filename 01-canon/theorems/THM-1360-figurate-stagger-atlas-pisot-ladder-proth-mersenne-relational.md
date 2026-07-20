---
id: THM-1360
title: The figurate-stagger atlas — staggering Pascal by slope s gives the recurrence a(m)=a(m−1)+a(m−s−1), whose growth is the DESCENDING PISOT LADDER x^{s+1}=x^s+1 (2, φ, supergolden 1.4656, 1.3803, …); the Proth grid n·2^x+1 staggered-and-summed is the MERSENNE number 2^{m+1}−1 (the Fermat "+1"s and the doubling telescope to the "−1", proved); and the RELATIONAL identity 2^{C(n−1,2)} = #tilings makes x = triangular the natural exponent, so 2^{T}+1 = tilings + observer.
status: >
  ATLAS / CONVERGENT RE-DERIVATION (not primary — see the priority banner; the
  deeper primary results are kind-pasteur-S128c103/S128c102, opus-S420, mac-mini-S1).
  VERIFIED-EXACT / PROVED — (A) Σ_c C(m−s·c, c) satisfies a(m)=a(m−1)+a(m−s−1),
  char x^{s+1}−x^s−1=0, roots 2, φ=1.61803, 1.46557 (supergolden), 1.38028, …
  descending to 1 (verified s=0..3). (B) Σ_{n=0}^{m}(n·2^{m−n}+1) = 2^{m+1}−1
  (proved: (m+1) + (2^{m+1}−m−2)); the Proth/Fermat family sums to Mersenne. (C)
  2^{C(n,2)} = #labeled graphs on n vertices, 2^{C(n−1,2)} = #tilings (the tiling
  hypercube Q_{T_{n−2}}); C(n−1,2) = T_{n−2} is a triangular number = the arc/tile
  count. (D) the figurate zoo comparison (Pascal/simplicial/polygonal/power-sum/
  centered/k-ary × slopes 0/1/2 × sum) — all computed and identified; the owner's
  power-sum triangle's shallow diagonal is the cubic Pisot 1.75488 (x³=2x²−x+1).
  INTERPRETIVE — the "one moral" framing (below); several individual identities are
  classical (diagonal sums of Pascal), so this is an ATLAS/organizing result, not a
  deep new theorem (MISTAKE-197 discipline).
source: death-star-2026-07-20-S59t (HYP-8175; owner: stagger the n·2^x+1 grid, sum/product, find continuations of triangular numbers, the relational reading).
depends_on:
  - THM-1355 (the n·2^x+1 Proth unification)
related:
  - opus-S317 (Vandermonde truncation law: polygonal vs polyhedral), THM-867/868/854
  - THM-227 (Mersenne/Fermat as (2^k−1)/(2^k+1) = one Cayley address)
  - the tiling hypercube Q_{C(n−1,2)} (CLAUDE.md); THM-466 (H = OCF, the +1 = vacuum)
scripts:
  - 04-computation/figurate_zoo_deathstar_S59t.py -> 05-knowledge/results/figurate_zoo_deathstar_S59t.out
  - 04-computation/proth_stagger_and_relational_deathstar_S59t.py -> 05-knowledge/results/proth_stagger_and_relational_deathstar_S59t.out
---

# THM-1360 — the figurate-stagger atlas

> **PRIORITY / CREDIT BANNER (added S59t-amend, MISTAKE-199).** This owner prompt
> went FLEET-WIDE and was worked deeper, concurrently, by others — this file is an
> independent re-derivation that CONVERGES with theirs and does NOT hold priority:
> • the **Pisot ladder** (§1, x^{s+1}=x^s+1: 2, φ, supergolden, plastic at s=4) is
>   **mac-mini-S1's Pascal-slope-d Pisot tower** and **kind-pasteur-S128c103's shear
>   catalog** (HYP-8170) — theirs also give the **Proth 2^{1/s} spectrum** (√2 at the
>   Fibonacci-analog shear) and prove Pascal DOMINATES Proth at every shear;
> • the **products** (§4, which I left un-named as "hyperfactorial-class") are
>   **opus-S420's** exact identity: pure-grid shear-1 row product = **m!·2^{C(m,2)} =
>   ordered tournaments** (sum side = A000295 Eulerian; "sums count elements, products
>   count orientations");
> • the **exact triangle** (the "power-sum Moser-break") is **kind-pasteur-S128c102's
>   Rosetta triangle**: columns = Faulhaber power sums with EXACTLY three +1 deviations
>   {(6,4),(7,4),(7,5)}, penultimate = the hypotenuse law 1+2^{n-2}, diagonal law
>   D_m = F_{m+1}+r(m); and **series-2 = Σ_k T(n,k)/k** (forcing 29/6, the owner's 29/3
>   a typo). My genuinely-complementary bits: the specific **Proth→Mersenne SUM**
>   identity (§2) and the **relational 2^{T}+1 = tilings+observer** framing (§3) — and
>   even these overlap kp-S128c103's "shear-1 = twice Mersenne" and opus-S420's
>   2T_m+1 = Φ₆(m+1) = |PG(2,m)|. Read this file as a convergent cross-check, and see
>   those four for the primary, deeper results.

## 1. The slope ladder (PROVED): moving columns down births a Pisot spectrum

For any left-justified triangle, the **slope-s diagonal sum** is
S_s(m) = Σ_c T(m − s·c, c). For **Pascal** (T = C):

  S_s(m) = a(m) with **a(m) = a(m−1) + a(m−s−1)**, char **x^{s+1} = x^s + 1**.

| slope s | recurrence | growth constant | name |
|---|---|---|---|
| 0 | a=2a(m−1) | **2** | doubling (powers of 2) |
| 1 | a(m−1)+a(m−2) | **φ = 1.61803** | golden (Fibonacci) |
| 2 | a(m−1)+a(m−3) | **1.46557** | supergolden (Narayana) |
| 3 | a(m−1)+a(m−4) | **1.38028** | (x⁴=x³+1) |
| s→∞ | | → 1 | |

So the owner's "move each column down 0, 1, 2, …" applied to Pascal is exactly
the **descending Pisot ladder** — the three most famous growth constants
(**2, the golden ratio, the supergolden**) are Pascal's columns moved down 0, 1,
2, and the ladder continues to 1. Powers of 2 and Fibonacci are the first two
rungs of one family.

## 2. The Proth grid sums to Mersenne (PROVED)

Stagger the grid f(x, n) = n·2^x + 1 (the Proth numbers, THM-1355) at slope 1 and
sum the anti-diagonal:

  Σ_{n=0}^{m} (n·2^{m−n} + 1) = **(m+1) + (2^{m+1} − m − 2) = 2^{m+1} − 1**.

The **Fermat/Proth "+1"s** contribute the (m+1), the **doubling** contributes the
geometric part, and they telescope to the **Mersenne "−1"**. The +1 tower and the
−1 tower are two readings of one staggered grid — the atlas realization of THM-227
((2^k−1)/(2^k+1) = tanh(k ln2 / 2), Mersenne and Fermat as one Cayley address).
Higher figurate families shift the spectrum: the **simplicial** triangle's
slope-1 sum is also Mersenne 2^n−1; its slope-2 sum is the cumulative Fibonacci
(growth φ); the **power-sum** triangle (owner's third figurate triangle) grows
faster, its shallow diagonal being the cubic Pisot **1.75488** (x³ = 2x² − x + 1).

## 3. The relational identity (VERIFIED): x = triangular is the natural exponent

The triangular numbers are "the relation itself": **T_{n−1} = C(n,2) = |E(K_n)| =
the number of arcs** of a tournament on n vertices; **T_{n−2} = C(n−1,2) = the
number of tiles** of the staircase δ_{n−2} (the tournament tiling triangle). And:

  2^{C(n,2)} = **# labeled graphs** on n vertices (1, 2, 8, 64, 1024, 32768, …),
  2^{C(n−1,2)} = **# tilings** = the tiling hypercube Q_{T_{n−2}}.

So **x = a triangular number** is the natural graph/tournament value of the
n·2^x+1 exponent: the n=1 slice at x = T_{n−2} is

  2^{T_{n−2}} + 1 = (# tilings) + 1 = the **observer-augmented tiling count**.

The two axes of the Proth table meet at the relational point: the exponent axis
(2^x, the hypotenuse) is *itself* triangular when it counts relations, and the
+1 is the observer on top of the tiling hypercube. The relational reading
(triangular = pairs = arcs = tiles) is what stitches the figurate zoo to the
tournament model.

## 4. The zoo (VERIFIED, atlas)

Diagonal sums (slope 0/1/2) across six families — Pascal, simplicial (d-simplex
C(m−1+c, c+1)), polygonal (opus-S317), power-sum (Σk^c), centered polygonal
(1 + (c+2)T_{m−1}), and k-ary relations (C(m, c+1)) — are computed and identified
in the stored output. Highlights: Pascal → {2^n, Fibonacci, supergolden};
simplicial → {fast, Mersenne 2^n−1, cumulative-Fibonacci}; k-ary → {2^n−1,
cumulative-Fibonacci, …}; power-sum → {huge, ≈ Σk^c-break, cubic-Pisot}. Products
along diagonals grow superexponentially (log ∼ quadratic, hyperfactorial-class).
The convergence **simplicial-slope-2 = k-ary-slope-1 = cumulative Fibonacci** and
the recurring appearance of Mersenne and the supergolden family show the zoo is
one web, not six.

## 5. One moral (interpretive)

Every named growth constant this project cares about — the doubling 2 (Burnside,
H = I(Ω,2)), the golden φ (Fibonacci, the θ-flow Markov endpoint, the JC₂ golden
corner), the supergolden and its descendants, the cubic Pisot of the power-sum
triangle — is a **rung of one ladder**: a figurate triangle staggered at a slope.
Powers of 2 and Fibonacci are not separate phenomena; they are slope 0 and slope
1 of the same operation, and the triangular numbers (the relation, the arcs, the
tiles) are the columns being staggered.
