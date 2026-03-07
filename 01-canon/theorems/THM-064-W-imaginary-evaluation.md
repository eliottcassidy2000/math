# THM-064: W(i/2) Evaluation and H-Maximizer Characterization

**Type:** Theorem + Conjecture
**Certainty:** 4 -- VERIFIED (exhaustive n=3,4,5; sampled n=6,7; partially at n=11)
**Status:** PARTIALLY PROVED
**Added by:** kind-pasteur-2026-03-07-S27
**Tags:** #W-polynomial #forward-edge #Paley #maximizer #complex-evaluation

---

## Statement

### (i) W(i/2) is real at odd n, pure imaginary at even n

For any tournament T on n vertices:
- If n is odd: W(i/2) is real (W has only even powers of r)
- If n is even: W(i/2) is pure imaginary (W has only odd powers of r)

This follows immediately from the r-parity of W (THM-061).

### (ii) Polar decomposition

Each permutation P contributes to W(i/2) a term of modulus (1/sqrt(2))^{n-1} and argument pi/4 * (3(n-1) - 2k_P), where k_P = number of forward edges.

Therefore:

**W(i/2) = exp(i*pi*3(n-1)/4) / 2^{(n-1)/2} * A(-i)**

where A(x) = sum_k a_k x^k is the forward-edge generating polynomial.

### (iii) W(i/2) = 0 iff (x^2+1) | A(x)

Since A(x) has real coefficients and is palindromic (a_k = a_{n-1-k}), the condition A(-i) = 0 is equivalent to:

**(x^2 + 1) divides A(x)**

which decomposes into two real conditions:
- sum_k (-1)^k a_{2k} = 0 (alternating sum of even-indexed coefficients)
- sum_k (-1)^k a_{2k+1} = 0 (alternating sum of odd-indexed coefficients)

At odd n (where only even powers of r appear), the imaginary condition is automatic, and only one real condition remains.

### (iv) Explicit formulas via Fourier coefficients

Using W(i/2) = sum_j w_{n-1-2j} * (-1)^{(n-1)/2-j} / 2^{n-1-2j} at odd n:

**n=3:** W(i/2) = -2 + 2*t3

**n=5:** W(i/2) = 16 - 4*t3 + 2*t5

**n=7:** W(i/2) = -272 + 32*t3 - 4*t5 - 8*a33 + 2*t7

### (v) Vanishing characterizes H-maximizers at small odd n (VERIFIED)

| n | W(i/2) = 0 tournaments | Total | All max-H? | Max H |
|---|------------------------|-------|------------|-------|
| 3 | 2 (both 3-cycles) | 8 | YES | 3 |
| 5 | 24 (all regular) | 1024 | YES | 15 |
| 7 | 240 (Paley/BIBD class) | 32768 | YES | 189 |

At n=3,5,7: **W(i/2) = 0 if and only if T is an H-maximizer.**

### (vi) Failure at n=9 and n=11

At n=9: The H-maximizer (circulant {1,2,3,5}, H=3357) has W(i/2) = 342, NOT zero.
However, |W(i/2)| is minimized at the maximizer class. Corr(H, |W(i/2)|) = -0.845 across all n=9 tournaments.

At n=11: Paley T_11 (the H-maximizer with H = 95095) has W(i/2) = -10010, NOT zero.

So W(i/2) = 0 characterizes H-maximizers ONLY at n=3,5,7. At larger odd n, |W(i/2)| is minimized (but not zero) at maximizers.

### (vii) Even n: W(i/2) never vanishes

At n=4: W(i/2) is pure imaginary and nonzero for all 64 tournaments.
At n=6: W(i/2) is pure imaginary and nonzero for 5000 sampled tournaments.

---

## Invariant Condition

W(i/2) = 0 imposes a linear constraint on the OCF invariants:

**n=3:** t3 = 1

**n=5:** t5 = 2*t3 - 8 (satisfied only by regular tournaments with t3=5, t5=2)

**n=7:** 32*t3 - 4*t5 - 8*a33 + 2*t7 = 272 (Paley: 448-168-56+48 = 272)

These are EXTREMAL conditions: only the H-maximizing tournament class satisfies them.

---

## Connection to Forward-Edge Polynomial

The condition (x^2+1) | A(x) means the forward-edge polynomial A(x) has 4th roots of unity (+-i) as roots. For palindromic polynomials of degree d, this is equivalent to B(0) = 0 where A(x) = x^{d/2} B(x + 1/x).

At n=3: A(x) = 3(1 + x^2) = 3(x - i)(x + i). Complete factorization.

---

## Connection to Reduced Polynomial P(u, x) (THM-063)

**THEOREM (kind-pasteur-S28):** W(i/2) = (-1)^m / 2^m * p_0(2), where m = (n-1)/2 and p_0(x) is the constant coefficient of the reduced polynomial P(u, x) = G_T(t, x) / t^m evaluated at u = t + 1/t.

**Proof:** The palindromic forward-edge polynomial satisfies A(-i) = (-i)^{n-1} A(i). Since W(i/2) = ((-1+i)/2)^{n-1} A(-i) and P(0, 2) = A(i)/i^m, the ratio W(i/2)/P(0,2) = ((1+i)/2)^{n-1} * i^m = (-1)^m / 2^m.

**Consequences:**
- W(i/2) = 0 iff p_0(2) = 0 iff u divides P(u, 2) iff (t^2+1) divides E_T(t)
- For Paley T_7: p_0(x) = 112(x-2)(x-68/7), so x=2 is a root and P(u,2) = u*(189u^2 + 504u + 756)
- The leading coefficient p_m(x) = I(Omega(T), x) is ALWAYS the independence polynomial
- The constant coefficient p_0(x) controls the "imaginary evaluation"
- P(2, x) = n! for all x (universal constraint)

This connects THM-064 (W(i/2) vanishing) with the trivariate GF framework of THM-063.

---

## Open Questions

1. **Why does W(i/2) = 0 characterize maximizers at n=3,5,7 but not n=11?** Is there a structural reason related to the number of independent OCF invariants?

2. **Does W(i/2) = 0 imply H = max for ALL odd n?** NO — fails at n=9 (maximizer has W(i/2) = 342). The vanishing characterization is specific to n=3,5,7.

3. **Is |W(i/2)| monotonically related to H?** Strong evidence: Corr(H, |W(i/2)|) = -0.845 at n=9. At n=7, W(i/2) is always negative, and more negative for lower H. The anti-correlation persists even when vanishing fails.

4. **Other roots of unity:** Do higher roots of unity (cube roots, 5th roots) at r = e^{2*pi*i*k/m}/2 give additional constraints that help characterize maximizers?

---

## Scripts

- `04-computation/W_half_i_vanishing.py` -- Main verification
- `04-computation/W_complex_evals.py` -- Complex evaluations and Parseval analysis
- `04-computation/W_half_i_n9.py` -- n=9 evaluation showing failure of vanishing
