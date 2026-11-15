# A local limit law for αH + βP₂: what the exact data say, what rules out local proofs, and host operators

Session klein-2026-11-15-llt. **Chan–Pak remains OPEN; no local limit theorem is proved here.**

[HYP-9772](../../hypotheses/HYP-9772-single-kernel-conjecture-every-odd-n-from-36835-is-h-plus-f2-of-some-tournament-the-double-cone-count-and-fixed-kernel-line-coverage-is-one-dimensional.md)
reduces single-kernel coverage to a one-dimensional local limit law for S = αH + βP₂, where P₂ counts ordered covers
by two nonempty disjoint paths. This note measures that law exactly, rules out the standard local proof technique,
and records an exact operator calculus for the hosts.

## 1. Exact unit-scale data (FINITE-EXACT)

We count the classes on each line S = 3H + P₂ = N for every tournament of order n, using `paircensus.c` on
`gentourng`. The bulk is the central 60% of classes, and the smooth density is a moving average over ±25 lines.

| order | lattice of S | bulk lines | classes per line μ | relative fluctuation | dispersion D | worst line / smooth | correlation, lags 1–6 |
|---:|---|---:|---:|---:|---:|---:|---|
| 8 | 2ℤ + 1 | 1,832 | 2.3 | 1.54 | 5.6 | 0 | −0.31, 0.10, 0.00, −0.28, 0.44, −0.34 |
| 9 | 6ℤ + 3 | 2,941 | 39.1 | 0.23 | 2.0 | 0.25 | 0.00, 0.01, 0.04, 0.00, −0.01, 0.01 |
| 10 | 2ℤ + 1 | 46,892 | 124.5 | 0.17 | 3.7 | 0.37 | −0.38, 0.17, 0.08, −0.32, 0.44, −0.43 |

- *Ratio form.* The fluctuation equals √(D/μ) to two digits. So the counts behave like a smooth density times
  Poisson-like noise with bounded over-dispersion D.
- *Why it should improve.* μ grows superexponentially with the order (about 2.2·10⁹ classes over about 10⁵ lines at
  order 11). If D stays bounded, P(S = N) = smooth(N)·(1 + O(μ^(−1/2))) uniformly on the bulk, which is a ratio local
  limit law.
- *Even orders.* The oscillating correlations are fine structure, not an obstruction: the worst bulk line still
  carries 37% of the smooth density.

**The singular series is flat (NUMERICAL).** 10.1M uniform order-10 samples (`llt_test.py`) give class densities
equal to within 1% modulo 8, 16, 32, 3, 9, 27, 5, 7 and 24. There are no local factors at the uniform measure.

## 2. Single local moves cannot carry a unit-scale proof

**The lattice of S (PROVED, from THM-5172 (C)).** At odd order every S lies in 6ℤ + 3, so every change is a multiple
of 6. At even order the lattice is 2ℤ + 1.

**Local moves are coarse (NUMERICAL).** This is the input a switching proof would need: a move changing S by one
lattice step, available for most tournaments.

| order | moves examined | ΔS = 0 | ΔS = one lattice step | mean \|ΔS\|/S |
|---|---|---:|---:|---:|
| 9, single arc flips | 10,800 | 96 | 18 (\|ΔS\| = 6) | 0.16 |
| 10, single arc flips | 5,400 | 31 | 1 (\|ΔS\| = 2) | 0.13 |
| 10, cyclic-triangle reversals | 3,603 | 81 | 2 (\|ΔS\| = 2) | — |

- So a tournament has a unit move with probability about 10⁻².
- The moves number polynomially in n, while S grows like n!/2ⁿ, so this rate tends to 0. A **switching proof
  through single local moves is impossible**.
- A unit-scale law has to come from a global mechanism: the cancellation of many independent choices, as in a
  bilinear or circle-method argument.

## 3. The local limit conjecture, and what would prove it

**Conjecture (ratio local limit law).** Fix a kernel (α, β), and let μ_q be the tilted ensemble (THM-5127 (B)).
There are densities g_(n,q) such that, uniformly for N in the bulk of μ_q on the lattice of S,

    #{T of order n : αH(T) + βP₂(T) = N} = (1 + o(1)) · 𝔖_n(N) · g_(n,q)(N),

with a singular series 𝔖_n that is constant on the lattice for the uniform measure. Together with the
derived-kernel chain of HYP-9772, this implies HYP-9772, and hence Chan–Pak beyond the census.

**Proof routes, honestly assessed.**
- *L² version (ARGUMENT).* "Most lines are near the smooth density" is equivalent to the collision energy of S
  being (1 + o(1)) times its smooth value. This is the same kind of statement as HYP-9752, for the statistic S
  over all tournaments instead of a carrier row.
- *Bilinear version (ARGUMENT).*
  - For a reducible core A₁ ⇒ A₂ with independent random blocks, the double-cone count is the bilinear form
    H₁H₂ + F₁F₂ = ⟨(H₁, F₁), (H₂, F₂)⟩.
  - A Cauchy–Schwarz circle method, as in Bourgain–Kontorovich for Zaremba, would reduce minor arcs to 2D spread
    of the pair (H, F₂) of each block.
  - Major arcs would need equidistribution modulo q, which the flat singular series supports.
  - Such a method naturally yields *almost all* N, not all N.
- *Pointwise version.* This needs concentration for the class counts per line. Its data are the dispersion D ≈ 2–4
  and the worst-line ratios 0.25 and 0.37 above.

## 4. Host operators: the double cone acts linearly on the path-cover polynomial (PROVED)

Let D(B) be the double cone a → B → c → a. Then for every q ≥ 1

    P_(D(B))(q) = q · [ P_B(q+1) + q·P_B(q) − P_B(q−1) ],      P_B(0) = 0.

This follows from the cyclic-join slot decomposition with a one-vertex reader (P_(K₁)(q) = q). Brute force agrees on
60 random cores of order 1–6 for q = 1, 2, 3 (`doublecone.py`).

Consequences:
- *The kernel itself.* q = 1 gives H(D(B)) = H(B) + F₂(B), the single-kernel statistic of HYP-9772. q = 2 gives
  F₂(D(B)) = 2F₃(B) + 4F₂(B) − 2H(B).
- *Nested double cones.* D^m is a fixed linear operator 𝒟^m on the path-cover polynomial. So H(D^m(B)) is an
  explicit integer combination of P_B(1), …, P_B(m+1); for example H(D²(B)) = 2F₃ + 5F₂ − H.
- *Closed forms.* For B = TT_k (P(q) = q^k) the counts are exponential polynomials:
  - H(D(TT_k)) = 2^k + 1;
  - H(D²(TT_k)) = 2·3^k + 5·2^k − 1;
  - in general H(D^m(TT_k)) = Σ_i λ_(m,i) (i+1)^k.
- *The general framework.* Every host is a linear functional on the path-cover polynomial of its module, and hosts
  compose, so they form a semigroup of functionals. Chan–Pak follows if **one** fixed functional of that semigroup
  is onto the odd numbers above a bound; HYP-9772 is the case 𝒟 read at q = 1.
- *A trade-off for later kernels.* Functionals that weight higher P(q) mix in more independent path structure, but
  they also stretch the range of values. At a fixed order, the first kernel is the best-covered (§1).

## 5. Exact repairs: where they can and cannot come from

The swarm's repairs are compared in [README](README.md). This session adds two constraints:
- *Local moves cannot repair at unit scale.* §2: unit changes of S are vanishingly rare among local moves.
- *Bilinear repairs are coarse.* Swapping a block inside a product or a join moves the count in steps of the other
  block's coordinates. These are huge unless that block is tiny, and then the supply is tiny too.

So an exact repair wider than a steering error must come either from equidistribution (§3) or from a family with
explicit additive digits: carrier rows, the near-join H = h + 2^k B, and fibre swaps in steps of 24 (wave 37). Such
digit families are where the exact work now lives.
