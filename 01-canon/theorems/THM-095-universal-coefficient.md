# THM-095: Universal Coefficient Theorem for Cumulant-Cycle Hierarchy

**Status:** PROVED (algebraic + numerical verification)
**Proved by:** opus-2026-03-07-S46e
**Resolves:** OPEN-Q-023

## Statement

For any tournament T on n >= 2k+1 vertices, the coefficient of t_{2k+1} (the count of directed (2k+1)-cycles) in the even cumulant kappa_{2k} of the forward-edge distribution is exactly:

    coeff(t_{2k+1} in kappa_{2k}) = 2 / C(n, 2k)

where C(n,2k) is the binomial coefficient.

## Examples

- k=1: coeff(t_3 in kappa_2) = 2/C(n,2) = 4/(n(n-1))
  (This is the t_3 coefficient in Var = (n+1)/12 + 4t_3/(n(n-1)), from THM-089.)

- k=2: coeff(t_5 in kappa_4) = 2/C(n,4)
  (Verified in THM-093.)

- k=3: coeff(t_7 in kappa_6) = 2/C(n,6)
  (Verified at n=7 in S46d.)

## Proof

### Step 1: Forward path formula
The number of directed r-edge paths through r+1 distinct vertices equals:

    #fwd(r)path = sum_{S in C(V,r+1)} H(T[S])

where H(T[S]) is the Hamiltonian path count of the subtournament on S.

### Step 2: OCF contribution of new cycles
By OCF (H = I(Omega, 2)), when r+1 = 2k+1:

    H(T[S]) = 1 + 2*t_3(S) + ... + 2*t_{2k+1}(S) + [higher alpha terms]

Each global directed (2k+1)-cycle appears in exactly ONE (2k+1)-element subset
(its own vertex set). Therefore:

    sum_S t_{2k+1}(S) = t_{2k+1}

giving a contribution of 2*t_{2k+1} to #fwd(2k+1)path.

### Step 3: Moment expansion
The centered moment mu_{2k} = E[(fwd - mean)^{2k}] expands as:

    mu_{2k} = sum over multisets of 2k position indices

A "maximal cluster" is a set of 2k consecutive X_i indicators involving
2k+1 vertices. The contribution is:

    (2k)! * (n - 2k) * E[X_i * X_{i+1} * ... * X_{i+2k-1}]

where the (2k)! comes from the multinomial coefficient (choosing each
position once), and (n-2k) is the number of valid starting positions.

### Step 4: Combining
Each cluster expectation decomposes as:

    E[X_i...X_{i+2k-1}] = #fwd(2k+1)path / P(n, 2k+1)

The t_{2k+1} part is 2*t_{2k+1} / P(n, 2k+1). Multiplying:

    coeff(t_{2k+1} in mu_{2k}) = (2k)! * (n-2k) * 2 / P(n, 2k+1)
                                = (2k)! * 2 / P(n, 2k)
                                = 2 * (2k)! * (n-2k)! / n!
                                = 2 / C(n, 2k)

### Step 5: Moment-to-cumulant preservation
The cumulant kappa_{2k} = mu_{2k} - [products of lower cumulants/moments].
Since mu_{2j} (for j < k) depends only on cycle counts t_{2j+1} with 2j+1 < 2k+1,
the correction terms cannot introduce t_{2k+1}. Therefore:

    coeff(t_{2k+1} in kappa_{2k}) = coeff(t_{2k+1} in mu_{2k}) = 2/C(n, 2k). QED.

## Numerical verification

- k=1: Exact for all n >= 3 (proved algebraically in THM-089)
- k=2: Verified exhaustively at n=5,6; sampled at n=7,8 (THM-093)
- k=3: Verified at n=7 (149 F-classes, S46d)
- k=4: Algebraic formula verified for n=9,10,11,12
- Forward path formula verified at n=5,6,7

## Significance

This theorem establishes a beautiful structural result: the even cumulants
of the forward-edge distribution form a GRADED hierarchy where each level
captures one new family of odd cycles with a universal binomial coefficient.

The factor of 2 comes from OCF (each directed odd cycle contributes 2 to
the independence polynomial at x=2), and the 1/C(n,2k) comes from the
probability that a specific (2k+1)-vertex subset appears as a maximal
cluster in a random permutation position.

## Scripts

- `04-computation/universal_coeff_proof.py` — Full proof + verification
- `04-computation/kappa4_coeff_pattern.py` — k=2 discovery
- `04-computation/kappa6_test.py` — k=3 verification
