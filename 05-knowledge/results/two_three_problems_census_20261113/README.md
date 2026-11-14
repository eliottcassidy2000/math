# Problems of 2 and 3 touched by this repository, and the 3-smooth decoder lemma

mac-mini-2026-11-13-cuts (overnight session). File counts are the number of navigation, canon and hypothesis files
whose text matches the problem name, a rough index of how much each lane has been touched.

## Census (twenty lanes)

| # | problem | the 2 and the 3 | files |
|---|---|---|---:|
| 1 | Collatz (3x+1) | ×3 + 1 and ÷2 | 323 |
| 2 | Chan–Pak spectrum of H(T) | 3 = H(C3) and 5 = H(R_2) bricks; two-hub counts are sums of 2^α·3^β frames | 209 |
| 3 | Lonely runner | HYP-9174: the ×2×3 lonely spectrum is the Artin coordinate of (2, 3) | 180 |
| 4 | Mersenne primes 2^p − 1 | 7 = 2³ − 1 is the first Chan–Pak exception | 121 |
| 5 | Pillai / Catalan, 3^a − 2^b | Mihăilescu; Pillai's conjecture | ~98 (some hits are Catalan numbers) |
| 6 | Berggren tree of primitive Pythagorean triples | a ternary tree; Berggren non-backtracking and legal paths (THM-4807/4808) | 91 |
| 7 | Zeckendorf / golden numeration | golden, not 2 and 3 strictly; the eleven-squares torus links it to Collatz period 5 | 66 |
| 8 | Mahler's 3/2 problem (Z-numbers) | the MAHLER-THREE-HALVES-FRONTIER navigation file | 57 |
| 9 | Fractional parts of (3/2)^k | Mahler and ideal Waring | 40 |
| 10 | Wieferich primes, 2^(p−1) ≡ 1 mod p² | Wieferich orbit explosion (THM-4644) | 26 |
| 11 | Hadwiger–Nelson | 7 colours = Φ3(2) = N(2 − ω); interior bound now ≥ 10 | 25 |
| 12 | 3-smooth representations | this note's lemma; Chan–Pak frames | 23 |
| 13 | Ideal Waring problem g(k) | ⌊(3/2)^k⌋ | 14 |
| 14 | Thue / S-unit equations in {2, 3} | | 9 |
| 15 | Furstenberg ×2 ×3 | | 5 |
| 16 | abc for 2 and 3 | smooth triples 2^a + 3^b = c | 3 |
| 17 | Gelfond's digit-sum problems in bases 2 and 3 | | 2 |
| 18 | Erdős: ternary digits of 2^n (a 2 appears for n > 8) | | 1 |
| 19 | Fermat primes F_r = 2^(2^r) + 1 | Busch bricks 3 and 5 are F0 and F1 (THM-4847, THM-4644) | many |
| 20 | Primes with ord_p(2) odd (density 7/24, Hasse) | 7 is the first; P_7 = Cay(Z_7, ⟨2⟩) | this session |

## The 3-smooth decoder lemma (PROVED, folklore-style; not found in the repository)

**Lemma.** Every positive integer n is a sum of numbers 2^a·3^b, none of which divides another.

*Proof (induction, which is a decoder).* For n = 1, use 3^0. If n is even, double a representation of n/2;
doubling preserves the antichain. If n is odd, let 3^k ≤ n < 3^(k+1). Then n − 3^k is even and less than 2·3^k.
Represent (n − 3^k)/2 and double it. Each doubled term is even, so it does not divide 3^k. Each is less than
2·3^k, so its power of 3 is at most 3^(k−1), and 3^k does not divide it. ∎

The recursion "subtract the largest 3^k, halve" is a balanced-ternary-like decoder. It matches the two-hub horizon
X(L) ≈ 3^L/2.

**Relevance to Chan–Pak.** Two-hub counts are frame sums Σ 2^α·3^β (waves 6 and 8 of the
[waves note](../chanpak_waves_20261113/README.md)). The lemma shows the arithmetic target set is complete. For
example 7 = 3 + 4 and 21 = 9 + 12 are antichain sums, so the obstruction to 7 and 21 is that the frames are not
freely choosable. A proof of the Two-Hub Conjecture would show that the admissible frame multisets contain the
antichains produced by this decoder, possibly up to bounded corrections, for every odd n ≥ 23 other than 21.

## Pythagorean triples and ternary digits

Primitive Pythagorean triples form the Berggren ternary tree: every triple has a unique address in {A, B, C}*, and
the hypotenuses are an N-rational series over 3 letters. Two-hub counts are an N-rational series over 4 letters (36
states) whose "ternary digits" are the triply covered spine vertices.

Primitive hypotenuses are the odd numbers whose prime factors are all ≡ 1 mod 4. They are exactly the totals of the
square-sum 4-cycles (two two-square representations; see the hidden-compositions note). THM-4807/4808 already link
the Berggren non-backtracking tree to legal Hamiltonian-path reconstruction.

**Prediction P23.** For a primitive triple of Berggren depth d, its hypotenuse c is a two-hub count with words of
length ≤ d + O(1). *Test:* the decoder on all primitive hypotenuses below 10^6.
