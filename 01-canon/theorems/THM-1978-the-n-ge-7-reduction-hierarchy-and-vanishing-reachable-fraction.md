---
id: THM-1978
title: "THE n≥7 REDUCTION HIERARCHY AND THE VANISHING-REACHABLE-FRACTION LAW. Tournament iso classes stratify by three nested reduction principles of increasing fineness — order-join/SCC atoms = STRONG (1,1,6,35,353 for n=3..7), modular/substitution atoms = MODULAR-PRIME seeds (1,1,1,0,3,15,197 for n=1..7, the n=7 entry NEW, completing opus THM-1960's open census), and circulant character-generated (1,0,1,0,2). The reduction-reachable (reducible) fraction FALLS monotonically — strong-fraction 0.25,0.50,0.625,0.774,0.873 (n=4..8), prime-fraction jumps 0.268→0.432 at n=7 — so reduction principles reach an asymptotically null set. This is why n=7 is the wall: the transitive/stability laws (srange≤tr THM-1862, srange≤β THM-1845, GIT-unstable⇔transitive THM-1825) all break at n=7 driven by the SAME THM-1830 witness (3-cycle atom + (n−3) singletons, impossible below n=7), and the spectral completeness facts (skew-Seidel THM-1440, odd-cycle count THM-500) both first fail at n=7 while char-poly-tie collapse goes 89%→99.1% (n=7→8)."
status: >
  VERIFIED (boxeph-2026-07-21-S197). Modular-prime seed census computed exactly over all iso classes
  n≤7: 1,1,1,0,3,15,197 — the n=7 value 197 is new (completes opus THM-1960's open census); for n≥3
  every reducible tournament has a nontrivial module (its dominated part), so modular-primes ⊆ strong
  (197 ≤ 353). Strong counts 1,0,1,1,6,35,353 (n=1..7) and reducible fractions verified against
  A000568 (…,456,6880 at n=7,8) and A051337. Circulant census 1,0,1,0,2 from THM-1955. The "breaks
  at n=7" catalog and the common THM-1830 driver are a documented synthesis of the cited theorems
  (THM-1862/1845/1825/1440/500/925/HYP-3732), each independently verified in its own source; this
  theorem is the census + the vanishing-fraction law + the unifying observation, not a re-proof.
source: boxeph-2026-07-21-S197 (owner: look back through all tournament work at size n≥7 for patterns in the hard-to-enumerate larger sizes)
depends_on: []
related:
  - THM-1960  # opus: modular seeds + spectral substitution (the n=7 seed census completed here)
  - THM-1955  # my reduction DAG (reducible/circulant/prime census n≤7)
  - THM-1926  # my zeta (factors over strong components)
  - THM-1830  # the 3-cycle-atom witness — the common driver of the n=7 breaks
  - THM-1440  # skew-Seidel completeness fails n=7 (spectral); THM-500 odd-cycle count too
  - "07-reflections/the-n-ge-7-regime-what-breaks-what-survives-boxeph-S197.md"
script: 04-computation/modular_prime_census_n7_boxeph_S197.py, large_n_circulant_patterns_boxeph_S197.py (+ .out)
---

# THM-1970 — the n≥7 reduction hierarchy and the vanishing-reachable-fraction law

## The stratification (three nested reductions, finest last)

| reduction | atoms | count n=3..7 | invariants that factor |
|---|---|---|---|
| order-join / SCC | **strong** | 1, 6, 35, 353 (and 1 at n=3) | char_A, H, signed-R, ζ (THM-1862/1926/1936) |
| modular / substitution | **modular-prime seeds** | (n=1..7) `1,1,1,0,3,15,197` | regular-seed spectral split (opus THM-1960) |
| circulant / character | **circulant** | 1, 0, 1, 0, 2 | Gauss/Chebyshev (THM-1955) |

Modular-primes ⊆ strong for n≥3 (a reducible tournament's dominated part is a nontrivial module).
The **n=7 seed count 197** is computed here, completing opus THM-1960's open census; the full seed
sequence is `1,1,1,0,3,15,197`.

## The vanishing-reachable-fraction law

```
   n            : 4     5     6     7      8
   strong-frac  : .25   .50   .625  .774   .873
   modprime-frac: 0     .25   .268  .432
   asymmetric   :  —     —     —    .875           (|Aut|=1: 399/456 at n=7)
```

The **reducible** fraction (what the order-join recursion collapses to products) is only 12.7% at
n=7 and 8.7% at n=8; the strong / modular-prime / asymmetric interior grows to full measure. So the
reduction principles — powerful *where they apply* (they factor char_A, H, R, ζ over atoms and pin
the circulant thread to Gauss/Chebyshev characters) — reach an **asymptotically null set**. That is
the precise sense in which n=7 forces enumeration.

## Why n=7 (the common driver)

The transitive/stability laws all break at n=7 on the **same** witness — THM-1830's *one 3-cycle atom
+ (n−3) transitive singletons* — which cannot exist below n=7 (instability needs #singletons > n/2
while a nontrivial strong component needs ≥3 vertices): `srange ≤ tr` (THM-1862), `srange ≤ β`
(THM-1845), "GIT-unstable ⇔ transitive" (THM-1825). Independently, spectral completeness dies at the
same n: skew-Seidel spectrum (THM-1440, 11 vs 12 classes, one cospectral pair) and odd-cycle count
(THM-500) both stop determining at n=7, and char-poly-tie collapse then runs 89%→99.1% (n=7→8). n=7
is simply the first size where a small cyclic atom fits inside an otherwise-transitive frame — and
that is exactly the point where the reducible island stops covering the tournaments.

## The surviving islands (large-n structure that stays clean)

Regular/symmetric families remain enumerable and exact at every n: Paley `T_p` (p≡3 mod 4) doubly
regular, self-complementary, H-maximal, `c₃ = p(p²−1)/24`, eigenvalues `(−1±i√p)/2`, `|λ|²=(p+1)/4`;
every regular tournament on m shares `c₃ = m(m²−1)/24`; circulant iso counts `2,4,4,6,16,16,30`
(n=7..19) all on the `Re λ = −1/2` line, Paley the flat (zero-spread Gauss) member. These are the
theorem-friendly islands; the generic sea is the n≥7 wall.
