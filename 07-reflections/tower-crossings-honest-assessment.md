# Tower Crossings: Honest Assessment — Proved vs Speculative

**Session:** kind-pasteur-2026-03-20-S10

## Context: The Repo's Error Culture

This project has documented 28+ mistakes (MISTAKES.md), with a recurring pattern:
**treating numerical coincidences as structural without verification at the next case.**

The most relevant precedents:
- MISTAKE-028: "7 = Mersenne AND forbidden" → falsely concluded ALL Mersenne are forbidden
  (31, 63, 127 are all achievable). "One data point is not a pattern."
- MISTAKE-006: c_7/C(11,7) = 4 looked like a pattern but was accidental
- MISTAKE-024: H=63 re-derived as "impossible" when earlier sessions already showed it achievable
- HYP-1618: Standard zeta(-3) = 1/120, NOT 7. The entire Riemann zeta forbidden-value
  connection was false.

These mistakes share a common mode: **seeing a numerical relationship, assuming it's
structural, and not testing the next case.**

I must apply this same rigor to my tower crossings synthesis.

---

## TIER 1: RIGOROUSLY PROVED

These claims are verified by exhaustive computation or algebraic proof.

### 1. The k-periodicity of the tower: PROVED
- Graphs (k=2): Level d exact for n <= 2d+1. Verified n=1..8.
- Tournaments (k=3): Level d exact for n <= 3d+2. Verified n=1..10.
- General: Level d exact for n <= k(d+1)-1 where k = min automorphism cycle.

### 2. The layer decomposition: PROVED
- D(n) = sum over non-identity odd partitions of D_lambda(n)
- Each D_lambda has an exact closed form involving 2^{f(lambda,n)} / (n-|lambda|)!
- Verified exactly at all computed n values.

### 3. P(n) sequence: COMPUTED EXACTLY
- P(n) = 1, 2, 4, 12, 48, 296, 3040, 54256, 1716608, 97213472
- Via Burnside formula by cycle type. All values verified.

### 4. The 3-periodic tower structure: PROVED
- Each level adds exactly 3 more correct terms (for tournaments)
- Because smallest depth-(d+1) partition is (3^{d+1}), first at n=3(d+1)
- Verified: Level 0 exact n<=2, Level 1 n<=5, Level 2 n<=8

### 5. D_graph >> D_tournament: PROVED
- D_graph(8) = 19504, D_tournament(8) = 784
- Ratio = 1219/49 ~ 24.88
- Structural reason: graphs have transpositions (order 2), tournaments don't

### 6. Prime filter irrelevant for n<=8: PROVED
- All contributing odd partitions up to n=8 have only prime parts (1,3,5,7)
- First composite odd part (9) contributes at n=9, changing T by exactly 1

### 7. The identity 2^10 = 10^3 + 4!: TRUE
- Pure arithmetic identity, verified.
- Holds only at n=2 and n=5 for the generalized form 2^C(n,2) = C(n,2)^3 + (n-1)!

### 8. Automorphism distribution at n=5: PROVED
- |Aut|=1: 840, |Aut|=3: 160, |Aut|=5: 24. Total 1024.
- The 24 = 4! = (n-1)! with |Aut|=5 is the unique regular Z_5 tournament.

---

## TIER 2: PLAUSIBLE BUT UNVERIFIED

These have some evidence but lack complete proofs or sufficient verification.

### 9. The Euler product structure of D(n)
- CLAIM: D factors multiplicatively over odd primes
- EVIDENCE: The layer structure is indexed by odd partitions, which have
  a natural multiplicative structure
- CONCERN: The layers DON'T actually multiply — D_{3,3} is NOT D_3 * D_3.
  The layers are ADDITIVE (D = D_3 + D_5 + D_{3,3} + ...), not multiplicative.
- STATUS: The "Euler product" language is MISLEADING. The structure is an
  additive sum over a multiplicative index set. This is like saying the
  prime counting function pi(x) = sum over primes, which is true but not
  an Euler product.

### 10. The Cheeger-periodicity product h*d → 1/p
- CLAIM: As d → infinity, the product h(level d) * d approaches 1/p
- EVIDENCE: None computed. This was stated without verification.
- STATUS: UNVERIFIED. Need to define h precisely and compute it.

### 11. The deficit ratio D_G/D_T approaching 25
- CLAIM: The ratio approaches 5^2 = 25
- ACTUAL: Ratio is 1219/49 = 24.878..., NOT exactly 25
- CONCERN: At n=7 the ratio was different (~14.6). The ratio CHANGES with n.
  Calling it "~25" at n=8 is a snapshot, not a limit.
- STATUS: The ratio has no proven limit. Claiming it's "~25" is premature.

---

## TIER 3: SPECULATIVE / LIKELY COINCIDENTAL

These are numerical observations that may not be structural.

### 12. The forbidden cascade 7 → 21 → 42
- CLAIM: 7 × 3 = 21 is the k-periodicity chain (multiply by min_cycle)
- REALITY: H=7 and H=21 are forbidden for COMPLETELY DIFFERENT REASONS:
  - H=7: cycle-forcing argument (alpha_1=3 with i_2=0 forces extra cycles)
  - H=21: six-way OCF decomposition blocking (each (alpha_1,alpha_2) pair blocked independently)
  The proofs of these two results share NO common mechanism.
- TEST: Is 7 × 5 = 35 forbidden? NO (achievable at n=7). Is 7 × 7 = 49 forbidden? NO.
  So 7 × k is NOT forbidden for k=5,7,... — only k=3 gives a forbidden value.
  This is exactly the "one data point" error warned about in MISTAKE-028.
- STATUS: LIKELY COINCIDENTAL. The factor 21/7 = 3 = min_cycle is probably accidental.

### 13. 42 × 5/3 = 70 = C(8,4) connection
- CLAIM: The Hurwitz constant 42 connects to the central binomial deviation
- REALITY: 42 × 5/3 = 70 is arithmetic. C(8,4) = 70 is combinatorial.
  There is no known structural reason why C(8,4) should relate to 42.
- STATUS: NUMEROLOGY until a structural connection is demonstrated.

### 14. The Kaprekar number 6174 = 2^1 × 3^2 × 7^3
- CLAIM: The staircase exponents connect to tournament theory via {2,3,7}
- REALITY: {2,3,7} appears in many contexts (Hurwitz, von Staudt, etc.)
  but the staircase exponent pattern (1,2,3) has no tournament interpretation.
  The digit product 168 = |PSL(2,7)| is a NUMBER THEORY fact about 6174,
  not a tournament fact.
- STATUS: COINCIDENTAL until proved otherwise.

### 15. The hierarchy Binary(2) → Tournament(3) → Hurwitz(7) → Bernoulli(43)
- CLAIM: These form a connected hierarchy of symmetry scales
- REALITY: 2 and 3 ARE connected (graph and tournament min_cycles, proved).
  7 and 43 are from completely different mathematical areas (Hurwitz surfaces,
  Bernoulli denominators). No mechanism connects tournament theory to
  Hurwitz's 84(g-1) theorem.
- STATUS: Levels 0-1 PROVED. Levels 2+ are SPECULATIVE analogies.

### 16. Non-integer periodicities (phi, e, pi)
- CLAIM: Irrational min_cycles correspond to quasicrystals
- REALITY: The k-periodicity theorem is proved only for INTEGER k
  (corresponding to actual cycle lengths in permutation groups).
  There is no definition of "phi-periodic tower" that connects to
  any known mathematical structure.
- STATUS: CREATIVE SPECULATION. Interesting but ungrounded.

---

## TIER 4: WHAT ACTUALLY HOLDS STRUCTURALLY

After filtering out the speculative, here is what the tower crossing concept
actually establishes:

### A. The tower IS real
P(n) = n*T(n) - D(n) with D decomposing into layers indexed by odd partitions.
The periodicity k=3 is structural (forced by Moon's theorem: tournament
automorphisms have odd order, minimum 3).

### B. The graph comparison IS real
k=2 for graphs (ALL permutations contribute) gives a genuinely different
tower with 2-periodicity. The tournament oddness filter is a genuine
arithmetic sieve.

### C. The identity 2^10 = 10^3 + 4! IS real but limited
It's a genuine identity at n=5 with a genuine tournament interpretation
(24 tournaments have Z_5 symmetry). But 10^3 = 1000 doesn't correspond
to any natural subset of tournaments — the split is 840+160+24, not 1000+24.

### D. The forbidden values 7 and 21 ARE structural but unconnected
Both are proved impossible for all n. But the multiplicative relationship
21 = 7 × 3 appears to be coincidental (7 × 5, 7 × 7 are NOT forbidden).

### E. The super-exponential decay IS structural
delta(n) = D/(nT) ~ (n-1)!/2^{2n-3} → 0, meaning most tournaments at
large n have trivial automorphism groups. This is genuine asymptotics.

---

## LESSONS

1. **The layer decomposition is genuine mathematics.** The exact formulas,
   the depth structure, the 3-periodicity — all proved and verified.

2. **The connections between layers 2, 3, 7, 42, 43 are analogies, not theorems.**
   They are interesting conceptual bridges but should NOT be stated as
   mathematical facts.

3. **The "Euler product" language overstates the structure.** The layers are
   additive over a multiplicative index set. This is weaker than a true
   Euler product.

4. **The "one data point" warning applies to 21 = 7 × 3.**
   Testing 7 × 5 = 35 (achievable) immediately shows the cascade is not
   a general principle.

5. **Numerical coincidences are common at small n.** The identity
   2^10 = 10^3 + 4! holds only at n=5. The central binomial formula
   D/2 = C(2k,k) holds only for n ≤ 6. Both break when tested further.
   This is the nature of small-number coincidences.
