# THM-227: k-nacci Mersenne Identity

**Status:** PROVED
**Session:** opus-2026-03-16-S90ap
**Verification:** Exact symbolic computation (sympy) for k = 2, 3, 4, 5, 7, 11, 13, 17, 19

---

## Statement

**Theorem.** For the k-nacci companion matrix M_k (k x k, k >= 2) with characteristic polynomial p(x) = x^k - x^{k-1} - ... - x - 1, and for any 1 <= n <= k:

Tr(M_k^n) = 2^n - 1.

In particular, Tr(M_k^k) = 2^k - 1 (the k-th Mersenne number).

**Corollary 1 (k-independence).** The trace Tr(M_k^n) is independent of k for 1 <= n <= k. All k-nacci matrices agree on Tr(M^n) = 2^n - 1 within the "Newton range" n <= k.

**Corollary 2 (first departure).** At n = k+1, the Mersenne formula breaks by exactly (k+1):

Tr(M_k^{k+1}) = 2^{k+1} - k - 2.

---

## Proof

The trace sequence t(n) = Tr(M_k^n) satisfies Newton's identities. The k-nacci characteristic polynomial p(x) = x^k - x^{k-1} - ... - x - 1 has all coefficients equal to -1 (i.e., c_j = -1 for j = 0, ..., k-1).

Newton's identities, applied with -c_j = 1 for all j, give:

t(n) = sum_{j=1}^{n-1} t(n-j) + n, for 1 <= n <= k.

Equivalently: t(n) = sum_{m=1}^{n-1} t(m) + n.

**Proof by induction on n:**

Base case: t(1) = (empty sum) + 1 = 1 = 2^1 - 1.

Inductive step: Assume t(m) = 2^m - 1 for 1 <= m <= n-1 (where n <= k). Then:

t(n) = sum_{m=1}^{n-1} t(m) + n
     = sum_{m=1}^{n-1} (2^m - 1) + n
     = (2 + 4 + ... + 2^{n-1}) - (n-1) + n
     = (2^n - 2) - n + 1 + n
     = 2^n - 1.  QED.

Note: The value t(0) = k never enters the induction. The formula t(n) = 2^n - 1 is independent of k.

For Corollary 2: at n = k+1, the k-nacci recurrence (without the +n correction) applies:

t(k+1) = sum_{m=1}^{k} t(m) = sum_{m=1}^{k} (2^m - 1) = (2^{k+1} - 2) - k = 2^{k+1} - k - 2.

---

## Connections

1. **Forbidden H values:** Tr(M_3^3) = 7, the first forbidden Hamiltonian path count. The forbidden value of the k-nacci theory at power k is the k-th Mersenne number.

2. **Mersenne primes:** When 2^k - 1 is prime (k = 2, 3, 5, 7, 13, 17, 19, ...), the forbidden value is a Mersenne prime.

3. **Cayley line:** tanh(k*ln(2)/2) = (2^k - 1)/(2^k + 1). The numerator of the Cayley address of 2^k is Tr(M_k^k).

4. **Why all coefficients are -1:** The k-nacci polynomial p(x) = x^k - x^{k-1} - ... - 1 satisfies (x-1)*p(x) = x^{k+1} - 2x^k + 1. The uniform -1 coefficients are what make Newton's identities produce the geometric series sum, which is what yields the Mersenne numbers.

---

## Verification

| k | Tr(M_k^k) | 2^k - 1 | Match | Mersenne prime? |
|---|-----------|---------|-------|-----------------|
| 2 | 3 | 3 | YES | YES |
| 3 | 7 | 7 | YES | YES |
| 4 | 15 | 15 | YES | no (3*5) |
| 5 | 31 | 31 | YES | YES |
| 7 | 127 | 127 | YES | YES |
| 11 | 2047 | 2047 | YES | no (23*89) |
| 13 | 8191 | 8191 | YES | YES |
| 17 | 131071 | 131071 | YES | YES |
| 19 | 524287 | 524287 | YES | YES |

## Files

- `04-computation/prove_knacci_mersenne_s90ap.py` — full proof with verification
- `04-computation/tanh_powers_of_2_s90ao.py` — discovery script (tanh connection)
- `05-knowledge/results/prove_knacci_mersenne_s90ap.out` — computation output
