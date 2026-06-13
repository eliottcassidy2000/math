# THM-094: F_k(T) mod 2 is Tournament-Independent

**Type:** Theorem (proved modulo universal Taylor zeros mod 2)
**Certainty:** 4.5 — PROOF SKETCH COMPLETE (verified exhaustively n<=6, sampled n=7,8)
**Status:** PROOF SKETCH + VERIFIED
**Added by:** kind-pasteur-2026-03-07-S38
**Tags:** #f-polynomial #modular-arithmetic #redei #eulerian-numbers

---

## Statement

**Theorem:** For any tournament T on n >= 1 vertices,

  F_k(T) = A(n,k) (mod 2)  for all k = 0, 1, ..., n-1,

where A(n,k) is the Eulerian number.

Equivalently: F_k(T) = C(n-1, k) (mod 2) for all k, since A(n,k) = C(n-1,k) mod 2.

Equivalently: F(T,x) = (1+x)^{n-1} (mod 2) for all tournaments T.

**Corollary (Rédei's theorem):** F_{n-1}(T) is odd for all tournaments T. This is exactly
Rédei's theorem: every tournament has an odd number of Hamiltonian paths.

---

## Proof Sketch

### Step 1: Universal Taylor zeros mod 2

**Claim:** c_j(T) = 0 mod 2 for all j < n-1 and all tournaments T on n vertices.

Verified computationally:
- n=2,...,6: exhaustive (0 failures)
- n=7: 5000 random samples (0 failures)
- n=8: 3000 random samples (0 failures)

The (x-1)-adic valuation of A_n(x) mod 2 is exactly n-1 (since A_n(x) = (1+x)^{n-1} mod 2
and (1+x) = (x-1) + 2 = (x-1) mod 2). The claim says this maximal valuation carries over
to F(T,x) for ALL tournaments.

*The algebraic proof of this step follows the THM-086 pattern (DC induction + palindrome),
generalized from p=3 to p=2. The key ingredient is that (n-j-1)! is always even for the
relevant range of j, so tournament-dependent terms in c_j vanish mod 2.*

### Step 2: Rédei's theorem

c_{n-1}(T) = F_{n-1}(T) = #{permutations P where P[i]->P[i+1] for all i in T}.

This is the number of Hamiltonian paths following all tournament arcs. By Rédei's theorem
(1934), this count is always odd. Therefore c_{n-1}(T) = 1 mod 2 for all T.

### Step 3: Combining

From Steps 1 and 2:
  F(T,x) = 1 * (x-1)^{n-1} mod 2 = (1+x)^{n-1} mod 2 = A_n(x) mod 2.

Therefore F_k(T) = [x^k](1+x)^{n-1} = C(n-1,k) = A(n,k) mod 2. QED.

---

## Consequences

### F_k parity from Lucas' theorem

Since F_k(T) = C(n-1, k) mod 2, by Lucas' theorem:
  F_k(T) is odd  iff  every binary digit of k is <= the corresponding digit of n-1.

Examples:
- n=3 (n-1=10_2): F_k odd at k=0,2 (00,10). Even at k=1 (01).
- n=5 (n-1=100_2): F_k odd at k=0,4 only.
- n=7 (n-1=110_2): F_k odd at k=0,2,4,6.
- n=8 (n-1=111_2): ALL F_k odd.

### Strongest individual F_k result for any prime

For p=3, we proved c_j(T) = 0 mod 3 for j < val(n), but this gives LINEAR COMBINATIONS
of F_k vanishing, not individual F_k. The mod-2 result is stronger: every individual F_k
is determined mod 2.

### Connection to THM-086

THM-086 proves the mod-3 analogue: c_j(T) = 0 mod 3 for j < 2*floor((n-1)/2).
THM-094 is the mod-2 analogue with val = n-1 (maximal).
Combined with Rédei (the mod-2 "free parameter" is always 1), this pins down F mod 2 exactly.

For mod 3, the "free parameter" c_{val(n)} is NOT always the same mod 3 (~30-35% are 0),
so we cannot determine F mod 3 completely.

---

## Generalization to other primes

Computational investigation shows:

| p | Universal zeros match Eulerian val for n >= | Notes |
|---|----------------------------------------------|-------|
| 2 | ALL n >= 1                                   | F(T,x) = A_n(x) mod 2 |
| 3 | ALL n >= 3 (THM-086)                         | 1 free parameter for odd n |
| 5 | n >= 7 (= p+2)                               | Fails at n=5,6 |
| 7 | n >= 9 (= p+2, conjectured)                  | Fails at n=7 |

The threshold n >= p+2 corresponds to p | (n-2)!, which is needed for the
tournament-dependent part of c_2 to vanish mod p.

---

## Verification Record

| n  | Method     | Samples | F mod 2 patterns | Status |
|----|------------|---------|------------------|--------|
| 2  | exhaustive | 2       | 1                | PASS   |
| 3  | exhaustive | 8       | 1                | PASS   |
| 4  | exhaustive | 64      | 1                | PASS   |
| 5  | exhaustive | 1024    | 1                | PASS   |
| 6  | exhaustive | 32768   | 1                | PASS   |
| 7  | sampled    | 5000    | 1                | PASS   |

---

## Scripts

- `04-computation/fk_mod2_proof.py` — verification and DC step check
- `04-computation/taylor_zeros_mod_p.py` — multi-prime analysis
