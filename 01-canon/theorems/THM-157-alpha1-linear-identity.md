# THM-157: alpha_1 = linear(S_4,...,S_{p-1}) Identity

**Status:** PROVED (algebraic identity + exact rational verification)
**Author:** kind-pasteur-2026-03-13-S60
**Dependencies:** Eigenvalue expansion of c_k (standard)

## Statement

For any prime p and any circulant tournament T on Z_p with connection set S of size m = (p-1)/2, the total directed odd cycle count

  alpha_1(T) = sum_{k=3,5,...,p} c_k(T)

is an **exact linear function** of the eigenvalue moments S_{2r} = sum_{t=1}^m D_t^{2r}:

  alpha_1 = sum_{r=2}^{m} a_r * S_{2r} + C(p)

where the coefficients are:

  a_r = sum_{k=2r+1, k odd, k<=p} (2/k) * C(k,2r) * (-1/2)^{k-2r} * (-1)^r

and C(p) is a constant depending only on p (not the orientation).

## Key Properties

1. **S_{p-1} coefficient**: a_m = (-1)^{(p+1)/2}
   - p = 3 mod 4: a_m = +1
   - p = 1 mod 4: a_m = -1

2. **Number of free parameters**: m-1 = (p-3)/2 (since S_2 = mp/4 is constant)

3. **Overconstrained verification**: At p=23, 90 orbit types with 11 parameters (79 excess constraints), all satisfied exactly with rational arithmetic.

## Proof

Direct algebraic consequence of the eigenvalue expansion:

  c_k = (1/k)[m^k + 2 * sum_{r=0}^{floor(k/2)} C(k,2r) * (-1/2)^{k-2r} * (-1)^r * S_{2r}]

Summing over all odd k = 3, 5, ..., p:

  alpha_1 = sum_k c_k = sum_k (1/k)[m^k + 2 * sum_r coeff(k,r) * S_{2r}]

Interchanging the order of summation:

  alpha_1 = const(p) + sum_{r=2}^m [sum_{k=2r+1,odd,<=p} (2/k)*C(k,2r)*(-1/2)^{k-2r}*(-1)^r] * S_{2r}

The inner sum gives the coefficient a_r. The constant absorbs the m^k terms and the S_0, S_2 contributions (which are tournament-independent).

## Coefficient Table

| p  | S_4 coeff | S_6 | S_8 | S_{p-1} | n_types | n_params | Excess |
|----|-----------|-----|-----|---------|---------|----------|--------|
| 5  | -1        | -   | -   | -1      | 1       | 2        | -1     |
| 7  | -9/4      | 1   | -   | 1       | 2       | 3        | -1     |
| 11 | -115/32   | 143/24 | -19/4 | 1   | 4       | 5        | -1     |
| 13 | -975/256  | 385/48 | -175/16 | -1 | 6      | 6        | 0      |
| 17 | -2013/512 | 2569/256 | -2973/128 | -1 | 16 | 8       | **8**  |
| 19 | -16155/4096 | 10633/1024 | -13881/512 | 1 | 30 | 9    | **21** |
| 23 | -1035445/262144 | 695237/65536 | ... | 1 | 90 | 11    | **79** |

## Significance

1. **First overconstrained algebraic identity** in this project (THM-155, THM-156 were verified but not overconstrained)
2. The S_{p-1} sign change with p mod 4 connects to the Gauss sum structure (g^2 = (-1)^{(p-1)/2} * p)
3. Combined with THM-156 (alpha_2 = (p/4)*E + C), gives: H = 1 + 2*linear(S4,...,S_{p-1}) + 4*linear(S4) + 8*alpha_3 + ...
4. The question of whether H is itself a linear function of moments reduces to whether alpha_3 (and higher) are linear in moments

## Computational Evidence

- `04-computation/alpha1_identity_all_primes.py` — full verification
- `04-computation/p17_exact_arithmetic.py` — first overconstrained verification
- `05-knowledge/results/alpha1_identity_all_primes.out` — output
