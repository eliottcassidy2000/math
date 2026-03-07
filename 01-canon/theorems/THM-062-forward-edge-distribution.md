# THM-062: Forward-Edge Distribution and Derivatives of W(r)

**Type:** Theorem (PROVED)
**Certainty:** 5 -- PROVED (algebraic + verified exhaustive/sampled)
**Status:** PROVED
**Added by:** opus-2026-03-07-S32
**Tags:** #W-polynomial #Eulerian #forward-edges #derivatives #master-polynomial

---

## Definitions

For a tournament T on n vertices and permutation P in S_n, define:
- **Forward edge at step i:** A[P(i), P(i+1)] = 1
- **a_k(T):** number of permutations with exactly k forward edges

The forward-edge distribution {a_0, a_1, ..., a_{n-1}} encodes how W(r) "spreads" its weight.

---

## Statement

### (i) Palindromy (reversal symmetry)

**a_k(T) = a_{n-1-k}(T)** for all tournaments T and all k.

*Proof:* Reversing P maps k forward edges to (n-1-k) forward edges (since reversed edges flip direction). Reversal is an involution on S_n. QED.

This is equivalent to the r-parity of W(r): W(-r) = (-1)^{n-1} W(r).

### (ii) Endpoint values

- a_{n-1}(T) = H(T) (all edges forward = Hamiltonian path)
- a_0(T) = H(T) (by palindromy)

### (iii) First derivative: a_{n-2}

**W'(1/2) = (n-1) H(T) + a_{n-2}(T)**

*Proof:* W'(r) = sum_P sum_j prod_{i != j} (r + s_i). At r = 1/2, the product is nonzero only when all non-j positions have forward edges. If position j is the unique backward edge, the product is (-1). If no backward edges (Hamiltonian path), each of the (n-1) choices of j contributes 1. So W'(1/2) = (n-1) H(T) + sum over perms with exactly 1 backward edge of (-1), but actually: for perms with k forward edges, the contribution is k * 1^{k-1} * 0^{n-2-k} (from choosing j among the forward positions) ... more carefully: only perms with n-1 or n-2 forward edges contribute. Those with n-1 (H paths) contribute (n-1) each; those with n-2 (one backward) contribute 1 * (-1) ... wait, no.

*Cleaner proof via Eulerian expansion:* Since W(r) = sum_k a_k (r+1/2)^k (r-1/2)^{n-1-k}, differentiating and setting r = 1/2 (so (r-1/2) = 0): only terms where at most one factor of (r-1/2) survives contribute. The k = n-1 term gives (n-1) a_{n-1} = (n-1) H(T). The k = n-2 term gives a_{n-2} * 1 * (n-2 choose 0 from powers) ... Precisely, by the product rule applied to (r+1/2)^k (r-1/2)^{n-1-k}:

At r = 1/2: d/dr [(r+1/2)^{n-1}] = (n-1), and d/dr [(r+1/2)^{n-2}(r-1/2)] = 1.

So **W'(1/2) = (n-1) a_{n-1} + a_{n-2} = (n-1) H(T) + a_{n-2}(T).** QED.

### (iv) OCF formula for a_{n-2}

From THM-059, W'(1/2) = F'_{n-1}(1/2) + sum_I 2^{parts(I)} F'_{f_I}(1/2) I(T), where **F'_f(1/2) = 2^{f+1} - 2**.

Therefore: **a_{n-2}(T) = W'(1/2) - (n-1) H(T)**, which is an explicit OCF polynomial.

Concrete formulas:

| n | a_{n-2}(T) |
|---|-----------|
| 5 | 26 + 4 t3 - 8 t5 |
| 7 | 120 + 48 t3 - 12 t7 |

Note: at n=7, the t5 and bc terms cancel in the subtraction W'(1/2) - 6 H(T).

### (v) Second derivative: a_{n-3}

**W''(1/2) = 2 [C(n-1,2) H(T) + (n-2) a_{n-2}(T) + a_{n-3}(T)]**

The factor of 2 comes from d^2/dr^2 of a product: each pair (j, k) of differentiated positions contributes twice (from the symmetric second derivative).

From THM-059: W''(1/2) = F''_{n-1}(1/2) + sum_I 2^{parts} F''_{f_I}(1/2) I(T), where

**F''_f(1/2) = f(f-1) + 2(f-1)(2^{f+1} - f - 2) + 2 A(f+1, 2)**

with A(n, k) the Eulerian number.

At n=7: **a_{n-3}(T) = 1191 + 30 t3 - 18 t5 + 30 t7 - 36 bc**

The constant 1191 = A(7, 2) is the Eulerian number.

### (vi) COMPLETE CLOSED FORM: Deformed Eulerian Numbers

**MAIN THEOREM.** For any tournament T on n vertices:

**a_k(T) = A(n, k) + sum_I 2^{parts(I)} * c_k^{(f_I, n-1)} * I(T)**

where the **inflated Eulerian coefficient** is:

**c_k^{(f, d)} = sum_{j=max(0,k-d+f)}^{min(f,k)} A(f+1, j) * C(d-f, k-j) * (-1)^{d-f-k+j}**

*Proof:* The forward-edge distribution gives W(r) = sum_k a_k p^k q^{n-1-k} where p = r+1/2, q = r-1/2. By THM-059, W(r) = F_{n-1}(r) + sum_I 2^{parts} F_{f_I}(r) I(T). The transitive tournament has all invariants = 0, and for it a_k = A(n,k) (forward edges = ascents). Each F_f(r) = sum_j A(f+1,j) p^j q^{f-j} must be "inflated" to the degree-(n-1) basis via (p-q)^{n-1-f} = 1, giving the c_k formula. QED.

**Special cases:**
- c_k^{(0, d)} = (-1)^{d-k} C(d, k) [signed Pascal's triangle row]
- c_k^{(d, d)} = A(d+1, k) [standard Eulerian numbers]

**Properties:**
- c_k^{(f, d)} = c_{d-k}^{(f, d)} [palindromic — matches a_k palindromy]
- sum_k c_k^{(f, d)} = 0 for f < d [total a_k = n! is invariant-independent]

**Explicit inflated coefficients at n=7 (d=6):**

| k | A(7,k) | c_k^{(4,6)} [t3] | c_k^{(2,6)} [t5,bc] | c_k^{(0,6)} [t7] |
|---|--------|-------------------|----------------------|-------------------|
| 0 | 1 | 1 | 1 | 1 |
| 1 | 120 | 24 | 0 | -6 |
| 2 | 1191 | 15 | -9 | 15 |
| 3 | 2416 | -80 | 16 | -20 |
| 4 | 1191 | 15 | -9 | 15 |
| 5 | 120 | 24 | 0 | -6 |
| 6 | 1 | 1 | 1 | 1 |

(t3 gets 2*c_k^{(4,6)}, t5 gets 2*c_k^{(2,6)}, t7 gets 2*c_k^{(0,6)}, bc gets 4*c_k^{(2,6)})

---

## Key Formula: F^{(k)}_f(1/2)

| k | F^{(k)}_f(1/2) |
|---|-----------------|
| 0 | 1 (for all f) |
| 1 | 2^{f+1} - 2 |
| 2 | f(f-1) + 2(f-1)(2^{f+1}-f-2) + 2 A(f+1, 2) |

---

## Verification

| Result | n | Tested | Status |
|--------|---|--------|--------|
| palindromy | 5, 7 | all / 30 samples | PASS |
| W'(1/2) formula | 5, 7, 9 | 5 each | PASS |
| a_{n-2} at n=5 | 5 | 10 | PASS |
| a_{n-2} at n=7 | 7 | 30 | PASS |
| a_{n-3} at n=7 | 7 | 30 | PASS |
| **full a_k formula at n=5** | **5** | **75 coefficients** | **PASS** |
| **full a_k formula at n=7** | **7** | **140 coefficients** | **PASS** |

---

## Connection to Other Results

- **THM-059:** The master polynomial F_f(r) provides the building blocks. The inflated Eulerian coefficient c_k^{(f,d)} is the change-of-basis from F_f's natural degree-f expansion to the degree-d forward-edge basis.
- **THM-061:** W(-1/2) = (-1)^{n-1} H(T) follows from evaluating both sides of the deformed Eulerian formula at r = -1/2.
- **Palindromy = r-parity:** The reversal bijection proving a_k = a_{n-1-k} is the simplest proof of W(-r) = (-1)^{n-1} W(r).
- **Transitive tournament:** For the transitive tournament (all invariants = 0), a_k = A(n, k), since forward edges correspond to ascents in the permutation.

---

## Scripts

- `04-computation/forward_edge_distribution.py` -- Part 1-4: W'(1/2), palindromy, a_{n-2} regression
- `04-computation/forward_edge_higher_derivs.py` -- F''_f(1/2) formula, a_{n-3} at n=7, general structure
- `04-computation/deformed_eulerian_verify.py` -- Complete closed-form verification at n=5 and n=7
- `04-computation/transitive_eulerian_check.py` -- Transitive tournament = Eulerian numbers
- `04-computation/eulerian_deformation.py` -- Delta analysis and regression
