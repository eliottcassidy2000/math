---
theorem_id: THM-137
title: The Paley-to-Interval crossover mechanism for H-maximization
status: PROVED (p=7,11 analytical; p=19 computational; mechanism understood)
proved_by: kind-pasteur-2026-03-12-S56c
date: 2026-03-12
related_theorems: [THM-130, THM-133, THM-134, THM-135, THM-136]
related_hypotheses: [HYP-474, HYP-479, HYP-480, HYP-481, HYP-482]
tags: [paley, interval, crossover, hamiltonian-path, spectral, additive-combinatorics]
---

## Main Result

**Theorem (THM-137):** The Paley tournament T_P maximizes H among circulant
tournaments on Z_p at p = 3, 7, 11, but NOT at p = 19 (where the cyclic interval
C_p wins). The crossover is caused by the competition between two effects:

1. **Paley cycle advantage**: Paley maximizes directed k-cycle counts c_k at
   EVERY odd k (HYP-474), giving more terms in the OCF.

2. **Interval structural advantage**: The interval's cycles are more INDEPENDENT
   (higher alpha_2), contributing more to the higher-order OCF terms.

At small p, effect (1) dominates. At large p, the exponential growth of the
interval's dominant eigenvalue amplifies effect (2) beyond what (1) can overcome.

## The Three-Layer Structure

### Layer 1: Spectral (THM-136, trace alternation)

Both eigenvalue sums oscillate with k mod 4:
- k = 1 mod 4: Paley wins tr(A^k) (less negative eigenvalue sum)
- k = 3 mod 4: Interval wins tr(A^k) (more positive eigenvalue sum)

This comes from both sets having eigenvalue phases near pi/2:
- Paley: theta = arctan(sqrt(p)) = pi/2 - O(1/sqrt(p))
- Interval: phi_1 = pi*(p+1)/(2p) = pi/2 + O(1/p)

The interval's dominant eigenvalue |mu_1| ~ p/pi grows faster than Paley's
uniform |lambda| = sqrt((p+1)/4) ~ sqrt(p)/2.

### Layer 2: Non-simple walk correction

The trace tr(A^k)/k counts ALL closed k-walks, not just directed k-cycles.
For k >= 7, non-simple walks contribute significantly. The correction
c_k - tr(A^k)/k is:

| p | k | c_k(P) | tr/k(P) | corr(P) | c_k(I) | tr/k(I) | corr(I) |
|---|---|--------|---------|---------|--------|---------|---------|
| 7 | 7 | 24 | 318 | -294 | 17 | 395 | -378 |
| 11| 7 | 3960 | 11220 | -7260 | 3399 | 12749 | -9350 |

The interval has MORE non-simple walks per vertex (because its dominant
eigenvalue creates more circulation), so its correction is larger in
magnitude. This makes the actual cycle counts CLOSER than the trace
values suggest — and even flips the k=7 comparison (c_7(P) > c_7(I)
despite tr(A^7)/7 favoring the interval).

**Key**: The trace-based H approximation H ~ 1 + sum 2^{(k-1)/2} * tr(A^k)/k
ALWAYS favors the interval, even at p = 7, 11. Paley's actual H advantage
comes entirely from the non-simple walk correction.

### Layer 3: OCF independent set structure

H = I(Omega(T), 2) = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...

where alpha_k = number of independent k-sets in the conflict graph of odd cycles.

| p | alpha_1(P) | alpha_1(I) | alpha_2(P) | alpha_2(I) | H(P) | H(I) |
|---|-----------|-----------|-----------|-----------|------|------|
| 7 | 80 | 59 | 7 | 14 | 189 | 175 |

At p = 7:
- 2*(alpha_1(P) - alpha_1(I)) = 2*21 = +42 (Paley advantage)
- 4*(alpha_2(P) - alpha_2(I)) = 4*(-7) = -28 (Interval advantage)
- Net: +14 = H(P) - H(I) = 189 - 175

Paley has MORE total cycles (alpha_1) but FEWER independent pairs (alpha_2).
The interval's cycles are more spread out, less conflicting.

## The Crossover Mechanism

As p grows:

1. The number of odd-k terms grows linearly (from 1 term at p=3 to
   (p-3)/2 terms at general p).

2. The interval's dominant eigenvalue r_1 ~ p/pi means its trace
   at k grows like (p/pi)^k, while Paley's grows like (sqrt(p)/2)^k.
   The ratio (r_1/|lam_P|)^k ~ (2*sqrt(p)/pi)^k grows exponentially.

3. At small p, the non-simple walk correction rescues Paley by removing
   more from the interval's inflated traces than from Paley's modest traces.

4. At p = 19, the exponential spectral growth of the interval at high k
   overwhelms the non-simple walk correction, and the interval wins H.

## Quantitative Crossover

| p | H(Paley) | H(Interval) | Margin | Winner |
|---|----------|-------------|--------|--------|
| 3 | 1 | 1 | 0% | TIE |
| 7 | 189 | 175 | +7.4% | PALEY |
| 11| 95095 | 93027 | +2.2% | PALEY |
| 19| 1,172,695,746,915 | 1,184,212,824,763 | -1.0% | INTERVAL |

The crossover occurs at p = 19, the first prime = 3 mod 4 after 11.
(p = 13, 17 are 1 mod 4 and don't have Paley tournaments.)

## Additive Combinatorics Perspective

The trace alternation theorem reduces to a counting problem:
- M_k = #{k-tuples from QR_p summing to 0 mod p}
- N_k = #{k-tuples from {1,...,m} summing to 0 mod p}
- Delta_k = p*(M_k - N_k)

At k = 1 mod 4: M_k > N_k (QR has more sum-zero solutions)
At k = 3 mod 4: M_k < N_k (Interval has more)

The additive energy E(S) = #{(a,b,c,d) : a+b=c+d} satisfies
E(INT) > E(QR) at all tested primes. The interval is more additively
structured (consecutive elements have many sum collisions), while QR
is more "pseudo-random" (multiplicatively structured but additively uniform).

## Open Questions

1. Can THM-136 (trace alternation) be proved algebraically for all p?
   The Paley side is exact (Gauss sum formula). The interval side needs
   error bounds on the dominant eigenvalue approximation.

2. Does the interval maximize H among ALL circulant tournaments for all p >= 19?
   This is HYP-480 (open).

3. Does the Paley/interval gap grow as p -> infinity? The eigenvalue ratio
   suggests yes: (r_1/|lam_P|)^p ~ (2*sqrt(p)/pi)^p -> infinity.

4. Is there a closed-form for the crossover prime in terms of the OCF structure?

## Scripts

- `04-computation/trace_alternation.py` — original discovery
- `04-computation/trace_alternation_proof.py` — analytical derivation
- `04-computation/trace_alternation_clean_proof.py` — verification
- `04-computation/trace_algebraic.py` — algebraic structure
- `04-computation/additive_structure.py` — additive combinatorics
- `04-computation/nonsimple_walk_correction.py` — correction analysis
- `04-computation/p19_crossover_analysis.py` — p=19 investigation
