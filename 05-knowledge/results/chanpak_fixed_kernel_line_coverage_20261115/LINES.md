# Does the pair spectrum cover every line? Fixed kernels at orders 9 to 12

Session klein-2026-11-15-lines. **Chan–Pak remains OPEN.**

## The question

Let V_n = {(H(A), P₂(A)) : A a tournament of order n}, where P₂ counts ordered covers by two nonempty disjoint paths.
[Wave 35](../chanpak_waves_20261113/README.md#wave-35-positivity-confinement--no-kernel-gives-cancellation-so-prime-targets-have-finitely-many-port-pairs)
reduces prime targets to *pair-spectrum fullness*: for some fixed small kernel, V should meet every line
αx + βy = p.

The cyclic join with a reader B (THM-5222) is such a kernel:

    H(B ⇒ ... ) :  N = (X_B + 2Y_B)·H(A) + X_B·P₂(A),    X_B = H(B) + P₂(B),  Y_B = H(B).

| reader B | kernel (α, β) | statistic | tournament realizing N |
|---|---|---|---|
| one vertex K₁ | (3, 1) | 3H + P₂ = H + F₂ | a → A → c → a, with A plus 2 vertices |
| TT₂ | (5, 3) | 5H + 3P₂ | the cyclic join with a 2-vertex reader |
| C₃ | (15, 9) | 15H + 9P₂ | the cyclic join with a triangle reader |

**Each pair lies on exactly one line of a fixed kernel.** So "V_n meets every line N in a window" says exactly
that the single statistic S = αH + βP₂ is onto the admissible N of that window. For one fixed kernel, line coverage
is therefore a **one-dimensional** statement. What it needs is enough pairs per line (σ > 1 in wave 35's notation),
not occupancy of the whole two-dimensional band (σ > 2).

## 1. Exact results through order 10 (FINITE-EXACT)

The census `paircensus.c` covers all 9,733,056 classes of order 10. It reproduces wave 25's counts of distinct
pairs, 2, 3, 8, 24, 156, 1,727, 50,173 and 1,550,357 for orders 3–10.

| kernel | order | distinct values | longest run with every admissible N covered | ratio |
|---|---:|---:|---|---:|
| K₁ (3, 1) | 9 | 7,087 | [10,198, 41,588], N ≡ 3 (mod 6) | 4.08 |
| K₁ (3, 1) | 10 | 125,869 | **[36,834, 241,522], all odd N** | **6.56** |
| TT₂ (5, 3) | 10 | 329,444 | [188,268, 623,232], all odd N | 3.31 |
| C₃ (15, 9) | 10 | 329,444 | [564,802, 1,869,698], N ≡ 3 (mod 6) | 3.31 |
| TT₃ (9, 7) | 10 | 664,721 | [840,710, 891,008] | 1.06 |
| order-4 readers (17, 15), (27, 21) | 10 | ≥ 664,721 | runs of ratio ≤ 1.06 | — |

So at order 10 a single fixed kernel already meets **every** line in a window of ratio 6.6. The cyclic joins
K₁ → A → c → K₁ with |A| = 10 realize every odd N in [36,835, 241,521].

This window starts exactly where wave 25's coverage started (36,833). That run used every pair of census cores; here
one fixed reader suffices on that stretch.

**Multiplicity per line in the bulk** (the lines carrying the central 80% of pairs):

| kernel, order | mean pairs per line | minimum | 1% quantile | uncovered lines |
|---|---:|---:|---:|---:|
| K₁, 8 | 0.53 | 0 | 0 | 1,490 of 2,588 |
| K₁, 9 | 9.5 | 2 | 4 | 0 |
| K₁, 10 | 17.9 | 4 | 9 | 0 |
| TT₂, 9 | 1.2 | 0 | 0 | 7,007 of 34,525 |
| TT₂, 10 | 6.5 | 1 | 3 | 0 |

Coverage switches on as the mean multiplicity passes a few pairs per line. The kernels with β/α near 1 spread the
same pairs over more lines and switch on later.

## 2. An order-parity obstruction for odd readers (PROVED, from THM-5172 (C))

THM-5172 (C) gives N ≡ 0 (mod 3) for a cyclic join whose two cores both have odd order. Hence:
- with any odd-order reader (K₁, C₃, TT₃, ...), the classes N ≡ 1, 2 (mod 3) come **only from even-order cores**;
- a single odd reader therefore needs windows two orders apart to overlap. Their scales differ by about
  (n+1)(n+2)/4 · (2n+1)/(2n−3), which is about 37 at n = 10;
- an even-order reader such as TT₂ has N ≡ 2H(A) (mod 3) for cores of every order, so every order contributes all
  three classes, and consecutive windows are only about (n+1)/2 apart.

## 3. The junctions at orders 11 and 12 (WITNESSED; sampled, not exhaustive)

Order 11 has 2.2·10⁹ classes, too many for an exact census of P₂ here. We therefore sample the **tilted ensemble**
of THM-5127 (B): orient i → j with probability q for i < j.
- Tilting moves the statistic from the bulk (q = 1/2) down to the transitive end (q → 1). So the junction regions,
  in the low tails of the larger order, are reached directly.
- Every sampled tournament is an explicit witness. The sampler is deterministic in (n, q, seed, index).
- `verify_witnesses.py` replays its xorshift128+ stream in Python and recounts (H, P₂) by a separate DP. It agreed
  on all spot checks: 3 order-12 witnesses and 6 random order-11 witnesses from the certified window.

Runs: order 11 used 8,000,000 samples at q ∈ {0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.88, 0.9}, giving 5,344,956 distinct
pairs. Order 12 used 4,800,000 samples at q ∈ {0.85, 0.88, 0.9, 0.92}, giving 2,579,405 distinct pairs, plus a
targeted junction run (§3.3).

**3.1. K₁ at odd orders closes the junction with room to spare.**
- The sampled order-11 cores meet every K₁ line N ≡ 3 (mod 6) in [47,385, 1,294,989], a window of ratio 27.
- Together with order 10, every odd N ≡ 0 (mod 3) in [36,835, 1,294,993] is 3H(A) + P₂(A) for an explicit A of
  order 10 or 11.

**3.2. TT₂ at order 11, and the alternating chain, are sampling-limited.**
- TT₂ at orders 10–11 covers [188,269, 674,517] before its first gap.
- K₁ at order 10 followed by TT₂ at order 11 stalls at 242,207, just above the order-10 top.
- Here the statistic spreads the same samples over more lines. These gaps show the sampling density, not empty
  lines.

**3.3. K₁ at order 12 for N ≢ 0 (mod 3): the decisive junction.**
- In [241,521, 300,000] the base order-12 run put 287,200 samples on the 19,493 lines, 8.9 per line. It left 197
  lines empty, against about 3 expected under Poisson.
- The empty lines spread evenly over residues mod 4, 8, 16, 32, 64, 9, 27, 5 and 7. So the excess is
  over-dispersion (near-transitive samples repeat pairs), not a local obstruction.
- The saturation run below tests whether the empty lines vanish.

**Saturation (targeted run).** The order-12 runs used tilts q ∈ [0.91, 0.945], with H filtered to [1,500, 6,500];
4.78M samples landed in the window. The empty order-12 lines in [241,521, 400,000] with N odd and N ≢ 0 (mod 3)
fell as follows:

| targeted samples added | 0 | 1.19M | 2.39M | 3.58M | 4.78M |
|---|---:|---:|---:|---:|---:|
| empty lines (of 52,826) | 505 | 126 | 47 | 15 | 7 |

The decay is steady. **Order 12 by itself fills the junction** up to a handful of rare lines, and §3.4 covers those
too.

**3.4. Junctions are reducible, so they are partly exact.**
- *Why reducible tournaments dominate.* The low tail of a larger order (near-transitive tournaments) consists
  largely of *reducible* tournaments.
- *Their pairs multiply.* For an ordinal sum A ⇒ B every q-slot cover splits slot by slot, so
  P_(A⇒B)(q) = P_A(q)·P_B(q). In particular (H, F₂) = (H_A·H_B, F_A·F_B).
- *So they are computable exactly.* The reducible part of V_11, ..., V_14 follows from the census of orders ≤ 10
  (`reducible.py`). That gives 1.55M, 1.60M, 3.16M and 4.91M exact pairs, with no sampling.

| kernel | data | every odd N covered in | misses up to 4·10⁶ (N ≡ 0, 1, 2 mod 3) |
|---|---|---|---|
| K₁, N = H + F₂ | exact (census ≤ 10 + reducible 11–14) | [36,835, 243,197] | 18,113 / 60,261 / 59,672 |
| K₁ | exact + witnessed (samples, and products with blocks of order ≤ 4) | **[36,835, 1,238,097]** | 0 / 2,616 / 2,583 |
| TT₂, N = 3F₂ − H | exact | [188,269, 818,949] | 3,279 / 46,682 / 46,684 |
| TT₂ | exact + witnessed | **[188,269, 1,597,157]** | 0 / 104 / 88 |

- *A cleaner fact.* With exact data alone, K₁ misses **no** odd multiple of 3 up to 1.4·10⁶.
- *Where the misses are.* All remaining misses are in the classes N ≢ 0 (mod 3), which §2 assigns to even orders.
  They lie beyond the reach of the samples.

## 4. Residues in the tilted ensemble (NUMERICAL)

At order 12 the statistics 3H + P₂ and 5H + 3P₂ realize every residue class modulo 3, 5, 7, 9, 16 (odd classes)
and 25.
- At q = 0.85 the class frequencies are within 1–2% of uniform modulo 5 and 7, and within 14% modulo 3.
- Toward the transitive end (q = 0.92) multiples of 3 become over-represented: the max/min frequency ratio is 1.65
  modulo 3 and 2.3 modulo 9. This is the odd-cycle crowding seen in THM-5192 (E).

No class is missing, so a local limit law would carry a nonuniform singular series rather than an obstruction.

## 5. What this says about a proof

For one fixed kernel, Chan–Pak for large N needs two things at every order n ≥ n₀.

- **(a) Bulk coverage.** A unit-scale local limit theorem for the single statistic S = αH + βP₂ under the tilted
  measure μ_q, which is THM-5127 (B)'s ensemble. This is one-dimensional, strictly weaker than the two-dimensional
  band fullness of HYP-9692 and wave 35. The data show it at order 10 (ratio 6.6, minimum multiplicity 4) and at
  order 11.
- **(b) Junction coverage.** The low tail of the next contributing order must reach the top of the previous window.
  - With an odd reader that next order is n + 2 for two of the three residue classes (§2).
  - Under tilting, (b) is again a bulk statement for μ_q at a larger q, so (a) uniformly in q would give (b).
  - The 10 → 11 junction holds with a factor-5 overlap. The 10 → 12 junction also holds: the order-12 run
    saturates it (§3.3), and the union in §3.4 covers every odd N above the order-10 window up to 1,238,097.

**The derived-kernel chain (PROVED reduction).**
- *Adding a transitive block moves the kernel.* Append the transitive block TT_k to a core A of order n. Then
  A ⇒ TT_k has (H, F₂) = (H_A, 2^k F_A). So the K₁ statistic of the order-(n + k) core A ⇒ TT_k is the **derived
  statistic** H + 2^k·F₂ of A.
- *Each step doubles the scale.* The derived statistics slide the order-n window up by factors of about 2 per step
  of k. Wave 25's padding and readers do the same from the other side.
- *Hence the reduction.* K₁ coverage of every odd N ≥ N₀ follows from:
  - bulk coverage of the finitely many statistics H + 2^k F₂, k = 0, ..., 4, at every order n ≥ n₀ (each a 1D
    statement for μ_q);
  - window ratios above 2, so that consecutive windows overlap;
  - for the classes N ≢ 0 (mod 3), only the combinations with n + k even count.
- *What is already true.* At order 10 the k = 0 window alone has ratio 6.6. The exact reducible data close the
  TT₂ junction to 818,949 with no sampling.

**Conjecture (single kernel, K₁; HYP-9772).** Every odd N ≥ 36,835 equals 3H(A) + P₂(A) = H(A) + F₂(A) for some
tournament A. Equivalently, N is the Hamiltonian-path count of the double cone a → A → c → a.
- *Record.* Certified for every odd N in [36,835, 1,238,097] (exact + witnessed), and for every odd multiple of 3
  up to 1.4·10⁶ by exact data alone.
- *Relation to Chan–Pak.* It implies Chan–Pak above 36,835, which other constructions already verify to 2³².
- *Why state it.* It is cleaner than HYP-9692: one statistic, one family of tournaments, one dimension. By the
  derived-kernel chain it reduces to finitely many one-dimensional tilted local-limit statements per order.
