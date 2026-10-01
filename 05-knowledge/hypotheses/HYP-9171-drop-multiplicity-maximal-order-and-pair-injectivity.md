---
id: HYP-9171
title: "(a) The maximal multiplicity of a Syracuse drop satisfies max_{d <= X} m(d) = (1+o(1)) sqrt(2 log2 X), where m(d) = #{v >= 2 : (2^v - 3) | 6d+1}. (b) For every odd q, two consecutive drops determine the point: A -> (K(A), K(T(A))) is injective for T(A) = oddpart(qA+1) on odd A >= 1"
status: >
  OPEN.
  (a) PROVED bounds (THM-4530): sqrt(2 log2 X) - 2.5 <= max <= 1 + log2(6X+1)/2. FINITE-EXACT: exact records up to 13 copies
  (exhaustive below 2^100); log2 of the k-th record fits k^2/2 + 2k. A Zsigmondy-type lower bound for lcm(2^v - 3) would
  settle it (lane direction D76).
  (b) PROVED for q = 3 on both sheets (3x+1 and 3x-1; THM-4530, Theorem 4.5). EMPIRICAL for q = 1, 5, 7, 9, 11, 13 and 5x-1
  (A < 4e5).
source: collatz-procgen-20260922 session, cdiff lane (2026-10-01), Conjecture 2.5' and direction D78; promoted with THM-4530
related:
  - 01-canon/theorems/THM-4530-syracuse-drop-multiplicities-two-copies-over-z-densities-and-pair-injectivity.md
  - 01-canon/theorems/THM-4527-syracuse-in-owner-labels-microcosm-and-difference-spectrum.md
---

# HYP-9171 — drop multiplicities: maximal order, and two drops determine the point

**(a)** `m(d) - 1` is the number of moduli `2^v - 3` (`v >= 3`) dividing `6d + 1`.

- **Large values are easy to produce.** To make `m(d)` large, `6d + 1` must be divisible by many of the moduli. Their lcm grows like `2^(sum of the v's)`, apart from the small shared primes, and the cheapest choice uses the smallest `v`. That gives `log2 X ~ k^2/2`, i.e. `k ~ sqrt(2 log2 X)`.
- **The difficulty is the other direction.** One must exclude a large lcm-collapse among the moduli. Known Zsigmondy-type results for `2^v - 1` do not apply directly to `2^v - 3`.

**(b)** For `q = 3`, the proof eliminates the second drop and reduces to a single residue contradiction (`65 != 1 mod 6`). For general odd `q` the same elimination gives
`n_1 E = M_v M_v' (2^delta - 1)`-type identities with `M_v = 2^v - q`. A uniform argument is missing.
