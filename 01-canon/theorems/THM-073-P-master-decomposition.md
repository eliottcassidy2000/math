# THM-073: Master Decomposition of P(u,x)

**Status:** PROVED (algebraic identity, verified computationally at n=5 exhaustive, n=7 all cycle types)
**Proved by:** kind-pasteur-2026-03-07-S29

## Statement

Let T be a tournament on n vertices (n odd), with G_T(t,x) = t^m * P(u,x) where
m = (n-1)/2 and u = t + 1/t.

Define the **base polynomials** P_k(u,0) for odd k by:
  P_k(u,0) = A_k(t) / t^{(k-1)/2}  (expressed in u = t + 1/t)

where A_k(t) is the Eulerian polynomial of degree k-1.

Then the correction polynomial for an independent set I = {C_1, ..., C_p} of odd directed
cycles with lengths l_1, ..., l_p and S_I = sum(l_i - 1) satisfies:

**g_I(u) = P_{n-S_I}(u, 0) * (u - 2)^{S_I/2}**

and the full master decomposition is:

**P(u,x) = P_n(u,0) + sum_I c_I * x^{|I|} * P_{n-S_I}(u,0) * (u-2)^{S_I/2}**

where c_I = number of independent sets of type I in Omega(T).

## Key Properties

1. **P(2,x) = n! for all tournaments:** All corrections vanish at u=2 due to the (u-2) factor.

2. **The correction depends only on S = sum(l_i - 1):** Different cycle compositions
   with the same S have the SAME correction polynomial. This explains the null space
   of the typed independence polynomial.

3. **Factorization structure:** Each correction is a product of:
   - P_{n-S}(u,0): the Eulerian polynomial of the "free" vertices (those not consumed by cycles)
   - (u-2)^{S/2}: a zero of multiplicity S/2 at u=2, reflecting descent degrees of freedom consumed

4. **Verified values of P_k(u,0):**
   - P_1(u,0) = 1
   - P_3(u,0) = u + 4
   - P_5(u,0) = u^2 + 26u + 64
   - P_7(u,0) = u^3 + 120u^2 + 1188u + 2176
   - P_9(u,0) = u^4 + 502u^3 + 14604u^2 + 86728u + 126976

   All satisfy P_k(2,0) = k! (evaluating at u=2, i.e., t=1).

## Proof Sketch

The correction function in t-space for an independent set with f = (n-1) - S free positions is:

  g(t) = A_{f+1}(t) * (t-1)^{n-1-f}

Converting to u-space by dividing by t^m:
- A_{f+1}(t) is palindromic of degree f, so A_{f+1}(t)/t^{f/2} = P_{f+1}(u,0)
- (t-1)^{n-1-f} = (t-1)^S has S even (since each l_i-1 is even), so (t-1)^S/t^{S/2} = (u-2)^{S/2}

since (t-1)^2/t = t - 2 + 1/t = u - 2.

The total: g(t)/t^m = P_{f+1}(u,0) * (u-2)^{S/2} = P_{n-S}(u,0) * (u-2)^{S/2}.

## Examples

### n=5 (m=2)
P(u,x) = P_5(u,0) + t3 * x * 2 * (u-2) * P_3(u,0) + t5 * x * 2 * (u-2)^2

= (u^2+26u+64) + 2*t3*x*(u-2)(u+4) + 2*t5*x*(u-2)^2

### n=7 (m=3)
P(u,x) = P_7(u,0) + 2*t3*x*(u-2)*P_5(u,0) + [2*t5+4*bc]*x^{1,2}*(u-2)^2*P_3(u,0) + 2*t7*x*(u-2)^3

More precisely, grouping by S:
- S=0: P_7(u,0) [base]
- S=2: 2*t3 * x * (u-2) * P_5(u,0) [3-cycles]
- S=4: [2*t5*x + 4*bc*x^2] * (u-2)^2 * P_3(u,0) [5-cycles and (3,3)-pairs share same u-factor!]
- S=6: 2*t7 * x * (u-2)^3 [7-cycles]

## Connection to THM-067 (Mersenne Vanishing)

The c_1 coefficient of P_{m-1}(x) involves evaluating g_I(u) at the u^{m-1} level.
The coefficient c_1^{(f,d)} = 2^{f+1} - d - 2 vanishes when n = 2^{f+1} - 1 (Mersenne).
At n=7 (f=2, d=6): c_1^{(2,6)} = 0, making P_2(x) linear in x (the bc/alpha_2 term vanishes).

## Connection to Null Space

The null space of the map (cycle counts) -> (P-coefficients) arises because g_I(u) depends
only on S_I, not on the individual cycle lengths. Any two independent set types with the same
(|I|, S_I) produce the same contribution to P(u,x). The null space dimension equals:
  #{distinct cycle types} - #{distinct (|I|, S_I) pairs}
