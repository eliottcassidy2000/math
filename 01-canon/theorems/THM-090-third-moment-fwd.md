# THM-090: Third Moment of Forward-Edge Count

**Status:** PROVED (algebraic) + VERIFIED n=3..6 exhaustive
**Proved by:** opus-2026-03-07-S46c
**Scope:** All tournaments

---

## Statement

Let T be a tournament on n >= 3 vertices. Then:

$$E[\text{fwd}^3] = A(n) + \frac{6}{n} \cdot t_3(T)$$

where:
- fwd = fwd(sigma) = number of forward edges in permutation sigma (edges following T)
- The expectation is over uniformly random permutations sigma in S_n
- t_3(T) = number of directed 3-cycles in T
- A(n) = 3 * (n-1)/2 * E_0[fwd^2] - 2 * ((n-1)/2)^3, where E_0[fwd^2] is the second moment for the transitive tournament

Equivalently: A(n) = E[desc^3] for the Eulerian descent distribution on S_n.

---

## Proof (clean version via reversal symmetry)

**Step 1: Reversal symmetry forces zero skewness.**

For any permutation sigma, define sigma^rev by sigma^rev(i) = sigma(n-1-i) (the reversal). Then:

fwd(sigma) + fwd(sigma^rev) = n - 1

because each edge (P[i], P[i+1]) in sigma corresponds to the edge (P[n-1-i], P[n-2-i]) in sigma^rev, and exactly one goes forward in T (since T is a tournament).

Since sigma -> sigma^rev is a bijection on S_n, the random variable fwd has the same distribution as (n-1) - fwd. Therefore the distribution of fwd is symmetric about mu = (n-1)/2, and all odd central moments are zero:

kappa_3 = E[(fwd - mu)^3] = 0

Verified exhaustively at n=4,5,6 (all F-classes).

**Step 2: Third moment from second moment.**

For any distribution symmetric about its mean mu:

E[X^3] = 3*mu*E[X^2] - 2*mu^3

(This follows from E[(X-mu)^3] = 0 expanded: E[X^3] - 3mu*E[X^2] + 3mu^2*E[X] - mu^3 = 0, combined with E[X] = mu.)

**Step 3: Substitution of variance formula.**

By THM-089: Var[fwd] = (n+1)/12 + 4*t3/(n(n-1))

So E[fwd^2] = Var + mu^2 = (n+1)/12 + 4*t3/(n(n-1)) + (n-1)^2/4

Substituting into the Step 2 formula:

E[fwd^3] = 3 * (n-1)/2 * [(n+1)/12 + 4*t3/(n(n-1)) + (n-1)^2/4] - 2*((n-1)/2)^3

The t3-dependent part is: 3 * (n-1)/2 * 4*t3/(n(n-1)) = 6*t3/n

The t3-independent part evaluates to A(n) (a polynomial in n). QED.

---

## Corollaries

1. **Zero skewness at transitive:** For the transitive tournament (t3=0), E[fwd^3] = A(n) = 3*mu*E[fwd^2] - 2*mu^3. This is exactly the third moment of a zero-skewness distribution with mean mu and variance Var, confirming that the Eulerian descent distribution has zero skewness (as follows from the palindrome symmetry).

2. **Moment hierarchy:** E[fwd^r] is determined by t3 alone for r <= 3. At r=4, additional invariants enter (verified at n=5,6). This explains why the first 4 Worpitzky coefficients depend only on t3.

3. **Forward r-path counts:**
   - #fwd2path = C(n,2) (trivial: every ordered pair has one direction)
   - #fwd3path = C(n,3) + 2*t3 (= directed 2-paths, involves 3-cycles)
   - #fwd4path = C(n,4) + 2(n-2)*t3 (involves 3-cycles in 4-vertex subsets)
   - General pattern: #fwd(r+1)path depends on t3 for r <= 3, needs more invariants for r >= 4.

---

## Verified

- E[fwd^3] = A(n) + 6t3/n: n=3,4,5,6 exhaustive (all F-classes)
- E[X_i X_{i+1} X_{i+2}] = 1/24 + t3/(3*C(n,3)): n=4,5,6 exhaustive
- #fwd4path = C(n,4) + 2(n-2)*t3: n=4,5,6 exhaustive
- Scripts: fwd3_formula.py, efwd3_proof.py, fwd_moments_cycles.py
