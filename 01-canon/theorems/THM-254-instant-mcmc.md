# THM-254: The Instant MCMC Theorem (Polynomial Mixing)

**Status:** PROVED
**Session:** kind-pasteur-2026-03-17-S116n33

## Statement

For the random flip Markov chain on the m-dimensional hypercube {0,1}^m, and any observable f: {0,1}^m → R with max Walsh degree D:

**E[f(X_t) | X_0 = x] is a polynomial of degree D in z = exp(-2t/m)**

with rational coefficients depending on x.

Specifically:

**E[f(X_t) | X_0 = x] = Σ_{k=0}^{D} B_k(x) · z^k**

where B_k(x) = Σ_{|S|=k} ĥ_f(S) · χ_S(x) and ĥ_f(S) = (1/2^m) Σ_y f(y) · χ_S(y) are the Walsh-Fourier coefficients of f.

## Consequences

1. **At log-rational times** t = (m/2)·ln(q) for rational q: z = 1/q is rational, and E[f] is an **exact rational number**. No floating-point approximation needed.

2. **The polynomial has only D+1 coefficients**, regardless of the size 2^m of the state space. For tournament H at n=6: D=4, so 5 coefficients encode the ENTIRE mixing dynamics.

3. **After O(2^m) preprocessing** (computing the Walsh transform), any query E[f(X_t) | X_0 = x] at any time t is answered in **O(D) = O(1) time**.

4. **P(0) = B_0 = mean(f)** always. P(1) = f(x) always. P(z) interpolates between the starting value and the equilibrium mean.

## Proof

The heat kernel of the hypercube flip chain is:

K_t(x,y) = (1/2^m) Σ_S exp(-μ_S · t) · χ_S(x) · χ_S(y)

where μ_S = 2|S|/m. Therefore:

E[f(X_t) | X_0 = x] = Σ_y K_t(x,y) · f(y)
= Σ_S ĥ_f(S) · exp(-2|S|t/m) · χ_S(x)
= Σ_S ĥ_f(S) · z^{|S|} · χ_S(x)

where z = exp(-2t/m). Grouping by degree k = |S|:

= Σ_{k=0}^{m} [Σ_{|S|=k} ĥ_f(S) · χ_S(x)] · z^k
= Σ_{k=0}^{D} B_k(x) · z^k

since ĥ_f(S) = 0 for |S| > D (f has Walsh degree D). ∎

## Corollary (Heat Kernel at Log-Rationals)

At t = (m/2)·ln(q) with q rational: z = 1/q, and:

K_t(x,y) = ((q+1)/q)^{m-d} · ((q-1)/q)^d / 2^m

where d = Hamming distance d(x,y). This is an **exact rational number** depending only on d.

## Corollary (The Backward Trick)

For the transitive tournament at n=6: P(0) = P(2) = 29 = mean H. Evaluating at z=2 (corresponding to "backward time" t = -5·ln(2)) gives the EXACT equilibrium mean. Running backward for 3.47 flips achieves the same as running forward forever.

## Implementation

The `InstantMCMC` class in `04-computation/instant_mcmc.py` implements this theorem:
- `__init__(f_table, m)`: O(2^m) Walsh transform preprocessing
- `.polynomial(x)`: returns [B_0, ..., B_D] for starting point x
- `.expected_f_at_z(x, z)`: O(D) polynomial evaluation, exact rational
- `.robustness(x, k)`: E[f] after k random flips
- `.mixing_time(x, eps)`: time to reach eps-close to mean
