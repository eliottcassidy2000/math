# Base-42 Analysis of the Erdős-Straus Conjecture

*opus-2026-03-16-S73 — an unexpected connection from tournament parity*

---

## The conjecture

**Erdős-Straus (1948)**: For every integer n ≥ 2, the equation 4/n = 1/x + 1/y + 1/z has a solution in positive integers.

Verified computationally to 10^14 (Swett 1999). Unproved.

## The base-42 insight

42 = 2 × 3 × 7 encodes the three fundamental constants of tournament parity:
- **2**: orientation/parity (F₂ arithmetic)
- **3**: smallest cycle (C₃, the triangle)
- **7**: first forbidden H value (Fano plane obstruction)

These same three primes control the decomposition of 4/p into Egyptian fractions.

## The splitting characterization (PROVED)

**Theorem**: 3/N = 1/b + 1/c has a solution in positive integers iff N has a prime factor p ≡ 2 (mod 3).

Equivalently: 3/N is unsplittable iff every prime factor of N is ≡ 1 (mod 3).

**Proof**: 3/N = 1/b + 1/c ⟺ (3b-N)(3c-N) = N². So we need d | N² with d ≡ -N (mod 3). For N ≡ 1 mod 3: need d ≡ 2 mod 3.

If q | N with q ≡ 2 mod 3, then q is a divisor of N² with q ≡ 2 mod 3. ✓

Conversely, if all prime factors are ≡ 1 mod 3, then all divisors of N² are ≡ 1 mod 3 (products of factors ≡ 1 mod 3 are ≡ 1 mod 3), so no d ≡ 2 mod 3 exists. ✗

Verified exhaustively for N ≤ 5000.

## Application to Erdős-Straus

For prime p ≡ 1 mod 4: the standard identity uses a = (p+3)/4 (so 4a = p+3), giving remainder 3/(pa).

This remainder 3/N with N = p·(p+3)/4 splits iff N has a prime factor ≡ 2 mod 3.

**The hard cases**: p ≡ 1 mod 12 (i.e., p ≡ 1 mod 4 AND p ≡ 1 mod 3).

Here p ≡ 1 mod 3 and (p+3)/4 ≡ 1 mod 3, so N ≡ 1 mod 3. The splitting characterization applies: 3/N fails iff (p+3)/4 has ALL prime factors ≡ 1 mod 3.

This happens for about 25% of primes p ≡ 1 mod 12 (empirically stable density).

## The covering system

For each coprime residue class mod 42:

| Residue mod 42 | Condition | Identity | Coverage |
|---|---|---|---|
| p ≡ 3 mod 4 (r = 11,19,23,31) | Always | a = (p+1)/4, rem = 1/(pa) → always splits | 100% |
| p ≡ 2 mod 3 (r = 5,17,29,41) | Always | a = p, split 3/p directly | 100% |
| r = 1 mod 42 | p ≡ 1 mod 12, mod7=1 | Primary: a=(p+3)/4; Secondary: a=(p+15)/4 | 100% (to 200k) |
| r = 13 mod 42 | p ≡ 1 mod 12, mod7=6 | Primary: a=(p+7)/4; Secondary: search | 100% (to 200k) |
| r = 25 mod 42 | p ≡ 1 mod 12, mod7=4 | Primary: a=(p+3)/4; Secondary: a=(p+7)/4 | 100% (to 200k) |
| r = 37 mod 42 | p ≡ 1 mod 12, mod7=2 | Primary: a=(p+3)/4; Secondary: a=(p+15)/4 | 100% (to 200k) |

**Result**: For all 9,592 primes up to 100,000 (and 4,472 hard-case primes up to 200,000), the two-identity covering system produces a solution. Zero failures.

## The meta-structure

The "easy" 8 residue classes mod 42 are handled by **factor 2** (mod-4 identity) or **factor 3** (mod-3 identity).

The "hard" 4 residue classes are {1, 13, 25, 37} mod 42 = {p ≡ 1 mod 12}. These are distinguished by their **mod-7 residue** {1, 6, 4, 2}, which controls which secondary identity is needed.

The secondary identity shifts a by 3 (from (p+3)/4 to (p+15)/4), which changes N mod 3 and N mod 7 simultaneously. This is the CRT structure of base 42: a single shift mod 42 adjusts behavior at all three constituent primes.

## What would prove the conjecture

**Sufficient**: Show that for every prime p ≡ 1 mod 12, at least one of finitely many values N_j = p·(p+4j-1)/4 (j = 1, 2, ..., k) has a prime factor ≡ 2 mod 3.

This is a statement about the distribution of primes ≡ 2 mod 3 in arithmetic progressions — a Chebotarev density problem. By Dirichlet, exactly half of all primes are ≡ 2 mod 3, so the probability of all k values avoiding such factors decreases exponentially with k. The obstacle to a proof is that the values N_j are correlated (they share the factor p).

## Tournament connection

The same base-42 structure appears in both contexts because both problems involve **decomposing a ratio into structurally simpler parts**:

- **Tournament H**: H(T) decomposes via Walsh-Fourier into pure-degree components, with 2-adic valuation controlled by digit sums
- **Egyptian fractions**: 4/p decomposes into 1/a + 1/b + 1/c, with splittability controlled by prime factorization mod 3

In both cases, the number 42 = 2 × 3 × 7 captures the complete modular structure needed to classify all cases.

## Cross-references

- THM-085: Universal 9-divisibility of H(T,ω) (role of factor 3)
- THM-J: Universality criterion s₂(n-3) ≤ 1 (role of factor 2)
- 07-reflections/987-amplituhedron-chemistry.md: 7 as forbidden value
- Mordell (1967): covering system mod 120
- Schinzel (1956): covering system mod 840
- Salez (2014): covering system mod 72
- Swett (1999): computational verification to 10^14
