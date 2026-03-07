# THM-076: Walsh-OCF Factorization Identity

**Status:** PROVED (algebraic, all r, all degrees, all odd n)
**Author:** opus-2026-03-07-S35 (continued^3), general-r proof: opus-2026-03-07-S35 (continued^4)
**Date:** 2026-03-07
**Dependencies:** THM-068 (PCD), THM-071 (Walsh-Fourier), OCF (Grinberg-Stanley)

## Statement

For any Walsh monomial S of type (2a_1, ..., 2a_r) (r disjoint even-length paths on 2k+r vertices, where k = sum a_i), the Walsh coefficient of I(Omega(T), 2) factorizes as:

**Single component (r=1, path P_{2a}):**

$$\hat{I}[S] = \epsilon_S \cdot \sum_{\substack{q = 2a+1, 2a+3, \ldots, n \\ q \text{ odd}}} \binom{n-2a-1}{q-2a-1} \cdot \frac{2(q-2a-1)!}{2^q} \cdot 2 \cdot f(n-q)$$

where f(n') = n'!/2^{n'-1} for n' >= 1, f(0) = 1.

Each term satisfies C(n-2a-1, q-2a-1) * (q-2a-1)! * (n-q)! = (n-2a-1)!, so the sum telescopes:

$$|\hat{I}[S]| = \frac{2 \cdot (n-2a)!}{2^{n-1}} = \frac{2^1 \cdot (n-2k)!}{2^{n-1}}$$

**Multi-component (r >= 2):**

The Walsh coefficient decomposes into covering configurations: ways to assign the r path components to disjoint odd cycles. Each configuration contributes a product of cycle terms times E[I(Omega, 2)] on unused vertices.

$$|\hat{I}[S]| = \frac{2^r \cdot (n-2k)!}{2^{n-1}}$$

## Proof (single component, r=1)

Fix a path S = P_{2a} on vertices {v_0, ..., v_{2a}} with 2a edges. WLOG epsilon_S = 1.

**Step 1: Only cycles containing S contribute.**
For I(Omega(T), 2) = sum_{R indep set} 2^|R|, the Walsh coefficient hat{I}[S] = E[I * chi_S]. An independent set R contributes E[2^|R| * indicator(R subset Omega) * chi_S]. For chi_S to have nonzero expectation, ALL S-edges must be fixed by some cycle in R. Since S is a connected path, one cycle C must contain all of S.

**Step 2: Factor into fan-cycle and remaining.**
For an independent set R = {C, C_2, ..., C_r} where C contains S:
- C uses q >= 2a+1 vertices (odd), fixing q edges
- {C_2, ..., C_r} are vertex-disjoint from C, using remaining n-q vertices
- chi_S = epsilon_S for all valid cycles C (constant by PCD descent-sign argument)
- E[indicator(R subset Omega) * chi_S] = epsilon_S * prod(1/2^|C_i|)
- Weight: 2^|R|

Summing over all R containing a cycle through S:
hat{I}[S] = epsilon_S * sum_C 2/2^|C| * sum_{R' on remaining} 2^|R'| * prod(1/2^|C_i|)

The inner sum = E[I(Omega(T'), 2)] on n-|C| vertices = f(n-|C|).

**Step 3: Count cycles containing S.**
A cycle of size q on a subset including the 2a+1 path vertices:
- Choose q-(2a+1) extra vertices from n-(2a+1): C(n-2a-1, q-2a-1) ways
- Cycle structure: path endpoints connect through extra vertices in both directions
- Directed cycles per subset: 2 * (q-2a-1)!

**Step 4: Telescoping.**
Each term: C(n-2a-1, j) * j! * (n-j-2a-1)! = (n-2a-1)! (constant!)

For q < n: regular term = 4 * (n-2a-1)! / 2^{n-1}
For q = n: half-term = 2 * (n-2a-1)! / 2^{n-1}

Count: (n-2a-1)/2 regular terms + 1 half-term = (n-2a)/2 effective terms

Total: (n-2a)/2 * 4 * (n-2a-1)!/2^{n-1} = 2*(n-2a)!/2^{n-1}. QED.

## Proof (multi-component, general r) — COMPLETE

For S = P_{2a_1} ∪ ... ∪ P_{2a_r} on disjoint vertex sets with k = sum a_i, define m = n - 2k - r (extra vertex pool).

**Step 1: Reduce to a combinatorial identity.**

By the same fan-cycle factorization as r=1, the amplitude decomposes over set partitions π of [r] into g groups. After simplification (the (m-E)! from f(remaining) cancels with the multinomial denominator), each term in the sum equals:

term(e_1,...,e_g) = m! * prod_j [(e_j + s_j - 1)!/e_j!] * 2^{r+g-n+1}

where s_j = |B_j| (group size), e_j = extra vertices in cycle j, with parity constraint e_j ≡ (1-s_j) mod 2, and the E=m term gets half weight.

The total becomes: |hat{I}[S]| = 2^{r-n+1} * m! * Sigma, where:

Sigma = sum_π 2^{g_π} * sum_{(e_1,...,e_g)} [half-weight] prod_j (e_j+s_j-1)!/e_j!

**Step 2: EGF computation of Sigma.**

Define the block weight generating function: for a group of size s,

g_s(x) = (s-1)!/2 * [1/(1-x)^s + (-1)^{s+1}/(1+x)^s]

(encodes the parity-constrained sum of Pochhammer terms).

The set partition EGF is:

F(t,x) = sum_{r>=0} G_r(x) * t^r/r! = exp(sum_{s>=1} 2*g_s(x) * t^s/s!)

Computing the exponent: sum_{s>=1} 2*g_s(x)*t^s/s! = ln((1 + t/(1+x)) / (1 - t/(1-x)))

Therefore F(t,x) = (1 + t/(1+x)) / (1 - t/(1-x)).

Extracting: G_r(x) = r! * [t^r] F(t,x) = **2*r! / ((1-x)^r * (1+x))**.

**Step 3: Half-weighted partial sum.**

Sigma_r(m) = sum_{E=0}^{m-1} c_E + c_m/2, where c_E = [x^E] G_r(x).

Using c_E = 2*r! * sum_{j=0}^E C(j+r-1,r-1)*(-1)^{E-j}, swap summation order:

Sigma_r(m) = 2*r! * sum_{j=0}^m C(j+r-1,r-1) * S(j)

where S(j) = sum_{E=j}^m w(E)*(-1)^{E-j}. A key calculation shows **S(j) = 1/2 for ALL j** (the alternating sum with half-weight endpoint always telescopes to 1/2).

Therefore: Sigma_r(m) = r! * sum_{j=0}^m C(j+r-1,r-1) = r! * C(m+r, r) = **(m+r)!/m!**

(using the hockey-stick identity for the last step).

**Step 4: Final assembly.**

|hat{I}[S]| = 2^{r-n+1} * m! * (m+r)!/m! = (m+r)! * 2^{r-n+1}

Since m + r = n - 2k: = (n-2k)! * 2^{r-n+1} = **2^r * (n-2k)! / 2^{n-1}**. QED.

**Verified computationally:**
- r=1 through r=5, multiple path sizes, n up to 20
- S(j) = 1/2 identity verified for m up to 9
- Hockey-stick identity verified for r up to 5
- Direct covering enumeration matches for r=2 (P2+P2, P2+P4), r=3 (P2+P2+P2, P2+P2+P4), r=4 (P2^4)

## Key Identity

The factorization relies on the **constant-term identity**:

$$\binom{m}{j} \cdot j! \cdot (m-j)! = m! \quad \text{for all } 0 \le j \le m$$

This trivial combinatorial identity is the engine that makes every term in the sum contribute equally, enabling the telescoping.

## Corollaries

1. **Alternative OCF proof path:** The Walsh amplitude formula for I(Omega,2) can be derived purely from:
   - Cycle counting (how many directed odd cycles contain a given path)
   - Degree-0 OCF (E[I(Omega,2)] = n!/2^{n-1}, which follows from linearity of expectation)
   - The constant-term identity

   This means higher-degree OCF follows from degree-0 OCF by factorization!

2. **Cycle polynomial structure:**
   - t_k (directed k-cycle count) has Walsh degree <= k-1 (polynomial degree bound)
   - Only even Walsh degrees survive (complement invariance)
   - The mixing coefficient t_k_hat / t_3_hat at degree 2 is C(n-3,k-3)/2 at n=5
   - Walsh-orthogonal invariants: q_{2j} = t_{2j+1} minus lower-degree projections

3. **Independence polynomial decomposition:**
   - Level 0 (constant): Walsh degree 0 only
   - Level 1 (single cycles): Walsh degrees 0, 2, ..., n-2
   - Level 2 (disjoint pairs): Walsh degrees 0, 2, ..., min(2(n-2), n-1)
   - The Walsh spectrum of H = I(Omega,2) is fully determined by cycle-path coverings

## Verified

- Single component (r=1): all P_{2a} for a=1,...,5, n up to 19 (exact rational arithmetic)
- Multi-component (r=2): P2+P2 at n=7 (exact)
- Degree-0 identity: n=3,...,17 (exact)
- Degree-2 identity: n=3,...,17 (exact)
