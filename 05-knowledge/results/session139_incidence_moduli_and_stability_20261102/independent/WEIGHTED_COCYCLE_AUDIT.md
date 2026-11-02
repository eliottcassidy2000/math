# Weighted second-moment cocycle audit and a necessary sidecar

Status: VERIFIED independently, with symbolic identities and direct rational
rounding. The proof under audit is `../cocycle/PROOF.md`, section4.

For one four-set, the three signed perfect-matching products in Pf(S) use
pairwise disjoint sets of independent edge-sign variables. Each product has
square1 and mean equal to its corresponding product of K entries. Therefore

    E Pf(S)^2 = Pf(K)^2 + 3 - sum_matchings K_e^2 K_f^2,
    E eta(S)^2 = [Pf(K)^2 + 2 - sum_matchings K_e^2 K_f^2]/8.

Here `E eta(S)^2` means the expectation of the square, not `(E eta(S))^2`.
The mean is the cubic cocycle delta(t(K))/4 and may vanish by cancellation
while the second moment is positive. The uniform K=0 case has mean0 but
bad-face probability1/4 on each four-set.

Write q_e=(1-K_e^2)/4. Summing the products in the fixed-face formula gives

    sum_(foursets,matchings) K_e^2 K_f^2
      =3 binom(n,4)-4 binom(n-2,2) V+16 M22,

where V=sum q_e and M22 sums q_e q_f over unordered disjoint edge pairs.
Each fixed edge appears with exactly binom(n-2,2) disjoint partners. Since
e4(K) is the sum of squared four-principal Pfaffians, this gives

    E |support eta(S)|
      =[e4(K)-binom(n,4)+4 binom(n-2,2)V-16M22]/8.

No independence is asserted between different four-faces; linearity of
expectation is sufficient for this global step.

## The same cospectral family proves the sidecar is genuinely missing

Use K14=a,K25=b, the skew opposites, all other entries0, and the root-fixing
Householder Q from the weighted five-vertex spectral audit. Let K'=QKQ^T.
Custom exact bivariate polynomial arithmetic verifies

    V(K')=V(K),          e4(K')=e4(K),
    M22(K')-M22(K)=-a^2 b^2/54,
    Ebad4(K')-Ebad4(K)=a^2 b^2/27.

Both matrices have bipartite support (first three versus last two vertices),
so all triangle-product means vanish and every mean eta-coordinate is zero.
Thus even the complete mean charge, V, and the full skew spectrum together
do not determine the expected support. M22 is not only a convenient formula
coordinate: it distinguishes an explicit fiber of all those retained data.

For a=1/4,b=1/2, direct enumeration of all1024 binary roundings of each matrix
returns

| quantity | K | K' |
|---|---:|---:|
| V | 155/64 | 155/64 |
| e4 | 1/64 | 1/64 |
| all five mean charges | zero | zero |
| M22 | 901/1024 | 24319/27648 |
| expected bad-four count | 5/4 | 2161/1728 |

The M22 difference is -1/3456 and expected-support difference is1/1728.
This is distinct from the M23 disjoint edge/triangle sidecar in the weighted
Hamiltonian-zeta audit, though the same family exposes both losses.

## Reproduction and scope

`audit_weighted_cocycle.py` imports only the independently written standard-
library audit helpers in `audit_bridges.py` and `audit_integral_cocycle.py`.
It checks the full symbolic M22 difference and3136 exact binary roundings:
1024 for each cospectral member,1024 for uniform order5, and64 for an
unstructured rational order4 example. The latter has nonzero mean charge
-61/3360 and second moment251/1080, guarding against a zero-mean-only test.

Normal and -O outputs agree byte-for-byte. Replacing the required +2 in the
fixed-face formula by +1 is rejected by explicit RuntimeError in both modes;
see `weighted_cocycle_negative_controls.out`. The first local implementation
accidentally allowed Python float zero entries in a numeric matrix and its
mass-one guard rejected that run; numeric entries are now explicitly coerced
to Fraction before rounding. All reported receipts come from the repaired
exact-arithmetic implementation.

The binary support gap is pointwise in each rounded realization. It does not
become a gap for support of the mean cocycle. No new weighted-stability
constant is inferred from these identities alone.
