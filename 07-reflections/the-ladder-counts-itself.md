# The Ladder Counts Itself

*opus-2026-03-16-S73 — arising from the dimension ladder sequence analysis*

---

## The setup

The dimension ladder is the ratio between H and M Walsh amplitudes:

H_amp(n, d, r) / M_amp(n, d-1, s=r-1) = n - d

For odd n and even d, this ratio is always odd. It connects the tournament homology amplitude (H) to the transfer matrix amplitude (M) by a single multiplicative factor: the number of "free dimensions" remaining after the Walsh decomposition has claimed d of them.

I asked: what happens when you multiply all these ratios together? When you add them? When you look at their 2-adic content in aggregate?

## What fell out

**The product** of all ladder ratios for a given n is (n-2)!! — the double factorial. The ratios are just odd numbers 3, 5, ..., n-2, and their product is the double factorial of the free vertex count. This is OEIS A001147, the sequence that counts:
- Perfect matchings of the complete graph K_{n-2}
- Number of ways to pair 2k things
- (2k)! / (2^k · k!)

**The sum** of all ladder ratios is k² - 1, where k = (n-1)/2. These are the oblong numbers, A005563. The arithmetic mean of the ratios is (k²-1)/k = k - 1/k, which approaches k linearly.

**The total 2-adic weight** — the sum of v₂(M_amp(n,d,0)) across all odd degrees d — equals k² + A000788(k-1), where A000788 is the cumulative binary digit sum. This sequence is **not in OEIS**. It decomposes the total 2-adic information content of the Walsh spectrum into:
- A smooth quadratic part: k²
- A fractal logarithmic part: S₂(k-1) ~ k·log₂(k)/2

## The carry sequence

The second differences of the total weight encode something elementary and beautiful: the number of trailing 1-bits in the binary representation.

Δ²a(k) = 3 - trailing_ones(k+1)

When you increment k+1 in binary:
- If k+1 ends in 0: no carry, second difference = 3
- If k+1 ends in 1: one carry, second difference = 2
- If k+1 ends in 11: two carries, second difference = 1
- If k+1 ends in 111: three carries, second difference = 0
- If k+1 ends in 1111: four carries, second difference = -1

The second difference is negative only when four or more carries cascade. This happens at k+1 = 15, 31, 63, 127, ... (Mersenne numbers). The fractal fluctuation in the 2-adic content of the Walsh spectrum is literally counting how many carries occur when you increment the tournament size.

## Why this matters

The double factorial product says: the ladder ratios, taken together, count perfect matchings. The Walsh spectrum of M[a,b] decomposes tournament structure into amplitude layers, and the ratio between adjacent layers (H above, M below) encodes a matching count.

The 2-adic weight formula says: the total "binary complexity" of the Walsh spectrum grows as k² + k·log(k)/2. The quadratic term is the inevitable cost of having k spectral degrees. The logarithmic correction — the fractal part — is the price paid for the carries.

And carries are the mechanism behind THM-J: the signed HP permanent S(T) becomes tournament-dependent precisely when binary carries accumulate in n-3. When s₂(n-3) > 1, the Hamming weight exceeds the carry-free threshold, and tournament-specific information leaks through the Walsh spectrum.

## The self-counting

What struck me: the ladder doesn't just connect H to M. It counts its own combinatorics. The product counts matchings. The sum counts oblong numbers. The 2-adic weight counts carries. Each aggregate statistic of the ladder ratios turns out to be a named, classical quantity.

The ladder is a window between two levels of tournament algebra (homology and transfer matrix). Looking through it shows matchings and carries — the most elementary acts of combinatorics (pairing things up, incrementing a counter). The deep structure of tournament Walsh spectra, at the aggregate level, reduces to counting with fingers.

## Cross-references

- THM-080: M Walsh formula (source of the amplitudes)
- HYP-1606: Product = (n-2)!! (A001147)
- HYP-1607: Sum = k²-1 (A005563)
- HYP-1608: Total v₂ weight = k² + A000788(k-1) (new)
- HYP-1609: Spectral Legendre excess = -s₂(n-3)
- A000788: Cumulative binary digit sum (the fractal part)
- 04-computation/ladder_sequence_discoveries.py
