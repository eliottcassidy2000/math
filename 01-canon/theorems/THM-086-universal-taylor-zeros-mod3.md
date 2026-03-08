# THM-086: Universal Taylor Zeros of F(T,x) Modulo 3

**Type:** Theorem (proved for j < 3; inductive proof structure identified for all j)
**Certainty:** 4.5 — PROOF SKETCH COMPLETE (algebraic proof for c_0,c_1,c_2; DC induction + palindrome for higher j; verified n<=10)
**Status:** PROOF SKETCH + VERIFIED
**Added by:** kind-pasteur-2026-03-07-S37
**Tags:** #f-polynomial #modular-arithmetic #taylor-expansion #eulerian-numbers

---

## Statement

Let T be any tournament on n >= 3 vertices, F(T,x) = sum_k F_k x^k the forward-edge
polynomial, and c_j(T) = sum_k C(k,j) F_k(T) the j-th Taylor coefficient at x=1.

Define val(n) = 2 * floor((n-1)/2) = n-1 if n is odd, n-2 if n is even.

**Theorem A (Universal Taylor Zeros):**
  3 | c_j(T) for all tournaments T on n vertices and all j < val(n).

Equivalently: (x-1)^{val(n)} divides F(T,x) modulo 3.

**Corollary (Eulerian Conjecture):**
  If 3 | A(n,k) (Eulerian number), then 3 | F_k(T) for all tournaments T on n vertices.

**Theorem B (Sharpness):**
  c_{val(n)}(T) is NOT always divisible by 3. The fraction of tournaments with
  3 | c_{val(n)} is approximately 30-35% at small n.

---

## Connection to THM-085

THM-085 proves the j < 3 case: c_0, c_1, c_2 are all 0 mod 3 for n >= 5.
THM-086 extends this to j < val(n), which is much stronger at large n:
- n=5,6: val=4, so adds c_3=0 (from palindrome)
- n=7,8: val=6, so adds c_3,c_4,c_5=0 (c_3 from palindrome, c_4,c_5 NEW)
- n=9,10: val=8, so adds c_3,...,c_7=0
- n=2k+1: val=2k, so c_0,...,c_{2k-1} all 0 mod 3

---

## Proof of the Corollary from Theorem A

The (x-1)-adic valuation of A_n(x) mod 3 is exactly val(n). This follows from
the Desarmenien factorization A_n(x) = prod_i A_{d_i+1}(x^{3^i}) mod 3.

For n odd: F(T,x) mod 3 = c_{n-1}(T) * (x-1)^{n-1} mod 3.
The coefficient [x^k] of (x-1)^{n-1} is C(n-1,k)*(-1)^{n-1-k} mod 3.
This vanishes at k iff C(n-1,k) = 0 mod 3. By Lucas' theorem, C(n-1,k) = 0 mod 3
iff some base-3 digit of k exceeds the corresponding digit of n-1.
These are exactly the positions where A(n,k) = 0 mod 3 (by Desarmenien).

For n even: F(T,x) mod 3 = c_{n-2}(T)*(x-1)^{n-2} + c_{n-1}(T)*(x-1)^{n-1} mod 3.
Palindrome forces c_{n-1} = f(c_{n-2}) mod 3, leaving one free parameter.
The resulting zero pattern in F_k mod 3 again matches the Eulerian zeros.

**Why p=3 is special (S38 discovery):** The Eulerian conjecture holds for p=3 because
val_3(n) = n-1 at odd n (single free parameter) and palindrome forces one parameter at
even n. For p >= 5, val_p(n) < n-1, giving MULTIPLE free parameters in F(T,x) mod p.
Different tournaments can have different alpha, beta, gamma, producing F_k != 0 mod p
even where A(n,k) = 0 mod p. VERIFIED: the Eulerian conjecture FAILS for p=5 at n=7
(F_1 and F_5 are NOT always 0 mod 5 despite A(7,1) = A(7,5) = 0 mod 5).

---

## Proof of Theorem A (partial)

### Proved cases: j = 0, 1, 2

See THM-085 for the complete algebraic proof. Summary:
- c_0 = n! (tournament-independent, 3|n! for n>=3)
- c_1 = n!(n-1)/2 (tournament-INDEPENDENT, 3|c_1 for n>=3)
- c_2 = A_non + (n-2)!*dp(T), both terms 0 mod 3 for n>=5

### Proved case: j = 3 (via palindrome)

At n >= 5, palindrome F_k = F_{n-1-k} combined with c_0=c_1=c_2=0 mod 3
forces c_3 = 0 mod 3. The proof is by linear algebra over F_3:
the palindrome constraints on c_3,...,c_{n-1} with c_0=c_1=c_2=0
always produce c_3=0 as a consequence.

### Inductive proof via Deletion-Contraction (j >= 4)

**Key identity (THM-083 in Taylor form):** c_j(T) = c_j(T\e) + c_{j-1}(T/e),
where e is any arc of T, T\e is the deletion, T/e is the contraction.

**Induction scheme (verified computationally n=5..8, 0 failures):**

1. **T/e induction:** c_{j-1}(T/e) = 0 mod 3 for j-1 < val(n-1) by induction
   (T/e is an (n-1)-vertex tournament).

2. **Almost-tournament claim:** c_j(T\e) = 0 mod 3 for j < val(n)-1.
   T\e is a tournament with one arc deleted. Verified:
   - n=5: val(T\e) = 3 = val(5)-1. HOLDS (exhaustive)
   - n=6: val(T\e) = 4 >= val(6)-1=3. HOLDS (5000 samples)
   - n=7: val(T\e) = 5 = val(7)-1. HOLDS (5000 samples)
   - n=8: val(T\e) = 6 >= val(8)-1=5. HOLDS (5000 samples)

3. **Combined:** c_j(T) = 0 mod 3 for j < val(n)-1.

4. **Palindrome upgrade:** c_{val(n)-1}(T) = 0 mod 3 follows from palindrome
   F_k = F_{n-1-k} combined with c_j = 0 for j < val(n)-1.

**What remains for a full proof:** The almost-tournament claim (Step 2) itself
needs proof. It can likely be proved by a similar DC argument applied to T\e,
giving a nested induction. The base cases (n=3,4) are verified exhaustively.

### General pattern

For each j, 3 | c_j(T) for all tournaments T on n >= N(j) vertices:
- N(0) = N(1) = 3
- N(2) = N(3) = 5
- N(4) = N(5) = 7
- N(6) = N(7) = 9
- General: N(2k) = N(2k+1) = 2k+3

---

## Structural insight: F(T,x) mod 3 has few free parameters

For n odd: F(T,x) mod 3 = alpha * (x-1)^{n-1} mod 3, with alpha = c_{n-1}(T) mod 3.
This is a SINGLE free parameter! All F_k mod 3 are determined by one number.

For n even: F(T,x) mod 3 = alpha * (x-1)^{n-2} + beta * (x-1)^{n-1} mod 3,
with palindrome constraining beta = f(alpha). So again ONE free parameter.

This is a remarkably strong rigidity result: modulo 3, the entire n!-term polynomial
F(T,x) is determined by a single parameter, regardless of the tournament.

---

## Verification Record

| n  | val | Method    | Samples | THM-086 | Eulerian |
|----|-----|-----------|---------|---------|----------|
| 5  | 4   | exhaustive| 1024    | PASS    | PASS     |
| 6  | 4   | exhaustive| 32768   | PASS    | PASS     |
| 7  | 6   | sampled   | 10000   | PASS    | PASS     |
| 8  | 6   | sampled   | 5000    | PASS    | PASS     |
| 9  | 8   | sampled   | 3000    | PASS    | PASS     |
| 10 | 8   | sampled   | 500     | PASS    | PASS     |

---

## Scripts

- `04-computation/eulerian_zeros_from_palindrome.py` — palindrome analysis
- `04-computation/taylor_cj_mod3_analysis.py` — c_j divisibility patterns
- `04-computation/thm086_verify.py` — combined verification
- `04-computation/fk_mod3_conjecture.py` — original Eulerian conjecture check (S36)
