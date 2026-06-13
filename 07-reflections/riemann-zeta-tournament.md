# The Zeta Function Knows About Tournaments

*Written: opus-2026-03-17-S74. A reflection on the six genuine connections between the Riemann zeta function (and its relatives) and tournament parity theory.*

---

## The Question

When we found that 7 and 21 are the forbidden values in the H-spectrum — the only odd numbers that can never be the number of consistent rankings for any tournament — we noticed they looked like evaluations of a zeta function at negative odd integers. The standard Riemann ζ(-3) = 1/120, which is not 7. So the naive connection fails. But the correction is more interesting than the original claim.

## Six Genuine Connections

### 1. Von Staudt-Clausen: Why 42 is Inevitable

The Von Staudt-Clausen theorem says:

$$\text{denom}(B_{2k}) = \prod_{(p-1) | 2k} p$$

For k=3: the primes p with (p-1)|6 are p=2 (1|6), p=3 (2|6), p=7 (6|6). So denom(B_6) = 2·3·7 = 42.

These are *exactly* the three primes of tournament parity:
- **2**: orientation (a tournament is a directed graph)
- **3**: the smallest odd cycle (Condorcet paradox)
- **7**: the Fano obstruction (smallest forbidden H)

The number 42 is not a coincidence. It is the sixth Bernoulli denominator, the number that encodes exactly those primes whose "period" divides the dimension where tournament structure first becomes interesting.

### 2. Kummer's Carries = THM-J's Carries

THM-J says the signed Hamiltonian permanent S(T) is universal mod 2^{n-1} if and only if s_2(n-3) ≤ 1 — the binary digit sum of n-3 is at most 1. Kummer's theorem says v_p(C(m+n,n)) equals the number of carries when adding m and n in base p. Both are carry-counting conditions on binomial coefficients.

More precisely, Legendre's formula v_2(k!) = k - s_2(k) appears *directly* in our Walsh amplitude formula:

$$v_2(\text{amp}) = s - d - s_2(n-2-d)$$

The s_2 function — itself a fractal — controls tournament universality in the same way it controls p-adic valuations of factorials. The mechanism is the same: carries in binary addition.

### 3. Odd Cycles as Primes: The Finite Euler Product

The OCF says:

$$H(T) = I(\Omega(T), 2) = 1 + 2\alpha_1 + 4\alpha_2 + \cdots$$

where α_k counts collections of k non-overlapping odd cycles. This is formally identical to expanding a finite Euler product:

$$\prod_{\text{odd cycles } C} (1 + 2 \cdot \mathbf{1}_C)$$

In number theory, the Euler product ζ(s) = ∏_p (1-p^{-s})^{-1} represents the multiplicative structure of the integers, with primes as generators. In tournament theory, odd cycles are the "primes" — they generate all contributions to the Hamiltonian path count multiplicatively, through disjoint collections.

The OCF is a *finite* Euler product. A tournament on n vertices has finitely many odd cycles, so the product terminates. This is why tournament invariants are more tractable than number-theoretic ones: the "zeta function" of a tournament is a polynomial, not an infinite series.

### 4. Gauss Sums Bridge Paley Eigenvalues and L-Functions

The Paley tournament T_p has adjacency eigenvalues (-1 ± √(p*))/2 where p* = (-1)^{(p-1)/2} · p. The quantity √(p*) is exactly the Gauss sum g(χ_p) = Σ_{a=0}^{p-1} χ(a) ω^a, where χ is the Legendre symbol and ω = e^{2πi/p}.

This same Gauss sum appears in the functional equation of the Dirichlet L-function:

$$L(s, \chi) \leftrightarrow g(\chi) \cdot L(1-s, \bar{\chi})$$

So Paley tournament eigenvalues are *direct functions of the Gauss sum* that controls the L-function at p. When we compute Betti numbers from these eigenvalues (via the eigenspace decomposition THM-125), we are extracting topological information from the same algebraic object that controls prime distribution in arithmetic progressions.

The Betti ratio β_{m+1}/β_m = (m+1)/(m-3) for Paley tournaments may have a deeper interpretation through this L-function lens.

### 5. k-nacci Traces: The Real Source of 7 and 21

The corrected story: the forbidden values 7 and 21 come not from Riemann's ζ but from Newton's identities applied to k-nacci companion matrices.

For the tribonacci matrix M_3, Newton's identity gives:

$$p_3 = e_1^3 - 3e_1 e_2 + 3e_3 = 1 + 3 + 3 = 7$$

This uses e_1 = 1, e_2 = -1, e_3 = -1 (the first three k-nacci coefficients). Crucially, these coefficients are the SAME for ALL k ≥ 3, so Tr(M_k^3) = 7 universally. The value 7 also happens to be a Mersenne number (2^3 - 1) and appears in the tribonacci trace sequence.

For p_5 = 21: this is specific to the tribonacci (k=3) case, where p_5 = p_2 · p_3 = 3 × 7.

**Important correction (opus-S74 audit):** While the k-nacci traces *hit* both forbidden values, they do not *cause* forbiddenness. The tribonacci trace sequence [1, 3, 7, 11, 21, 39, 71, ...] contains both forbidden values (7 and 21) but also many achievable values (1, 3, 11, 39, 71, ...). The Mersenne numbers 31, 63, and 127 are all achievable H values (at n=6, 8, and 7 respectively). The actual mechanism is combinatorial: cycle-forcing for H=7 (THM-029) and OCF decomposition blocking for H=21 (THM-079). The best characterization: {7, 21} = {7·3⁰, 7·3¹} with the 7-obstruction having nilpotency 2 (HYP-1231), or equivalently {Φ₃(2), Φ₃(4)} (HYP-1317).

### 6. The Critical Strip Analogy

In the Riemann zeta function:
- **Trivial zeros** at s = -2, -4, -6, ... (negative even integers)
- **Nontrivial zeros** on the critical line Re(s) = 1/2 (conjectured)
- **The functional equation** connects s and 1-s

In tournament Walsh analysis:
- **Trivial vanishing** at odd Walsh degrees (parity cancellation — every monomial with an odd-degree component vanishes identically)
- **Nontrivial structure** concentrated at even Walsh degrees 0, 2, 4, ...
- **Complement duality** H(T) = H(T^op), β(T) = β(T^op) connects a tournament to its "reflection"

The parallel is structural, not numerical. The mechanism is the same: a symmetry (functional equation / complement involution) forces half the spectrum to vanish, concentrating all information on the surviving half.

## What Compels

The most striking aspect is not any single connection but their convergence. Five different aspects of the Riemann zeta function (Bernoulli denominators, carry-counting, Euler products, Gauss sums, functional equations) all find natural analogues in tournament theory. The tournament versions are *finite* — they terminate, they are computable, they are provable — but they share the same algebraic DNA.

This suggests that tournaments sit at a crossroads where combinatorics, algebra, and analytic number theory can see each other. The OCF is a finite Euler product. The Paley eigenvalues are Gauss sums. The signed permanent has Mertens-like cancellation bounded by carry conditions. And 42 = denom(B_6) is the number that organizes it all.

Perhaps the deepest lesson: the zeta function doesn't know about tournaments *directly*. Rather, both the zeta function and tournaments are downstream of the same arithmetic — the multiplicative structure of the integers, the carry pattern in binary addition, the quadratic reciprocity that builds Paley tournaments from Legendre symbols. They are siblings, not parent and child.

---

*Cross-references: THM-227, THM-J, THM-125, THM-130, HYP-1618 (corrected), HYP-1622-1628. See `04-computation/zeta_deep_dive.py` and `04-computation/riemann_zeta_tournament_connections.py` for all computations.*
