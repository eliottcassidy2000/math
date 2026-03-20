# THM-258: Paley Sub-Tournament Uniformity Theorem

**Status:** PROVED (computational, p=7 and p=11)
**Filed by:** kind-pasteur-2026-03-20-S1

## Statement

For Paley tournaments T_p (p = 3 mod 4 prime), the sub-tournament structure is
dramatically more uniform than for the interval tournament I_p. Specifically:

### p=11 Sub-Tournament Uniformity Comparison

| Subset size k | Paley CV(H_k) | Interval CV(H_k) | Paley # distinct H | Interval # distinct H |
|---------------|---------------|-------------------|---------------------|----------------------|
| 5 | 34.6% | 39.0% | 5 | 4 |
| 7 | 10.7% | 31.2% | 4 | 8 |
| 9 | **0.0%** | 15.6% | **1** | 5 |

### p=7 Sub-Tournament Uniformity (from THM-257)

| Subset size k | Paley CV(H_k) | Type II CV(H_k) | Type III CV(H_k) |
|---------------|---------------|-----------------|-----------------|
| 5 | **0.0%** | 20.8% | 14.9% |

## Key Discovery: Perfect Uniformity Cascade

**Paley T_7:** ALL 21 five-vertex sub-tournaments have H=13 (zero variance).
**Paley T_11:** ALL 55 nine-vertex sub-tournaments have H=2879 (zero variance).

The pattern: perfect uniformity occurs at subset size k = p-2 = n-2.

**Conjecture:** For all Paley primes p, every (p-2)-vertex sub-tournament of T_p
has the same H value. (Verified: p=7 with H_5=13, p=11 with H_9=2879.)

## Mechanism: Uniformity Maximizes Total Cycles

The total directed odd cycle count is the sum of directed Hamiltonian cycles
over all odd-size vertex subsets:

    alpha_1 = sum_{k odd, S in C(p,k)} c_k(T_p[S])

For uniform sub-tournaments (all H_k equal), every subset contributes equally.
For non-uniform sub-tournaments, some subsets contribute few cycles, dragging
down the total.

### Quantitative comparison at p=11:

| Metric | Paley | Interval | Paley/Interval |
|--------|-------|---------|----------------|
| 9-cycles: all H_9=2879 | 11055 | 9350 | 1.182 |
| 7-cycles: CV=10.7% | 3960 | 3399 | 1.165 |
| 5-cycles: CV=34.6% | 594 | 484 | 1.227 |
| 3-cycles: both c3=55 | 55 | 55 | 1.000 |

The advantage ratio is LARGEST where uniformity is weakest (5-cycles) and
approaches 1 where uniformity is strongest (3-cycles, identical for both).

## Spectral Interpretation

Paley's uniformity comes from its FLAT spectrum. The eigenvalues of T_p's
adjacency matrix are:
- Lambda_0 = (p-1)/2
- Lambda_k = (-1 +/- sqrt(-p))/2 for k != 0

All non-trivial eigenvalues have |Lambda_k| = sqrt(p)/2. This spectral flatness
creates uniform cycle distribution across vertex subsets.

The interval tournament has a CONCENTRATED spectrum with dominant eigenvalue
|Lambda_1| ~ p/pi, creating highly variable sub-tournament structure.

## Proof

Computational verification:
1. p=7: Exhaustive (all C(7,5)=21 five-vertex sub-tournaments). THM-257.
2. p=11: Exhaustive for all C(11,k) subsets at k=5,7,9.
   - 5-vertex: 462 subsets, H values {5:110, 9:110, 11:55, 13:110, 15:77}
   - 7-vertex: 330 subsets, H values {115:110, 127:55, 143:110, 151:55}
   - 9-vertex: 55 subsets, H values {2879:55} (PERFECT UNIFORMITY)

Scripts: `04-computation/paley11_subtournament_uniformity.py`
Results: `05-knowledge/results/paley11_subtournament_uniformity.out`

## Corollary: The Spectral Advantage Has Finite Range

As p grows, the interval's concentrated spectrum eventually wins because:
1. The dominant eigenvalue p/pi creates exponentially more long cycles
2. Paley's flat sqrt(p)/2 can't keep up at large p

The crossover occurs between p=11 and p=19 (THM-256).

## Open Questions

1. **Prove the perfect uniformity conjecture:** H_{p-2} is constant for all
   (p-2)-vertex sub-tournaments of T_p.
2. **Characterize the k-dependent uniformity:** What determines the number of
   distinct H_k values at each subset size?
3. **Jensen's inequality formalization:** Does the uniformity → higher alpha_1
   relationship follow from a convexity argument?

## Related

- THM-255: SC Maximizer Dichotomy at n=6
- THM-256: Paley vs Interval comparison
- THM-257: BIBD Uniformity at n=7
- THM-162: Spectral flatness and H-maximization
- THM-126: Paley uniquely maximizes among circulants at p=7
