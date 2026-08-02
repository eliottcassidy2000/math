---
id: THM-3115
title: "Low-degree monomial-fibre Newton refinement transport"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT / UNDER INDEPENDENT AUDIT.  In degrees
  five through eight, every chamber-Newton coefficient of both normalized
  THM-3110 product-Gamma response banks is a nonnegative transport from a
  Young-subgroup fibre type to coarser fibre types.  Consequently the
  row-normalized central operator is positive semidefinite for every integer
  anchored support 0<a<b.  This is new normalized row-minimum positivity,
  stronger than THM-3110's unnormalized low-degree Schur positivity, but not
  an all-degree product-Gamma or Gaussian-moment proof.
source: root/multiscale-newton-flag-2026-08-02
depends_on:
  - THM-3110-arbitrary-anchored-product-gamma-dominant-tail-and-low-histogram-reduction
  - THM-3112-cycle-weighted-young-subgroup-gap-and-uniform-octopus-boundary
script: 04-computation/gmc_monomial_fibre_newton_refinement_thm3115.py
output: 05-knowledge/results/gmc_monomial_fibre_newton_refinement_thm3115.out
script_sha256: 81c5a6381a66cddacd902c536c67062ab88466c858f1c4fba96bab2212d2fee5
output_sha256: fc2597d3cf94479ca7c11b1ff4096159cfbb3e0afe6d4f70ff7e09fabf96b2cf
secondary_script: 04-computation/gmc_monomial_fibre_refinement_n8_scout.py
secondary_output: 05-knowledge/results/gmc_monomial_fibre_refinement_n8_scout.out
secondary_script_sha256: 283c1d99b7a36514e97eed711e4dca2b2c2c4e15abd17e4f63732d48c7a9ada6
secondary_output_sha256: 30ce74dc4077b5738861f5247ce963c68e869538f3e616a72dec099d2bffef3b
hash_basis: LF-normalized bytes
---

# THM-3115 -- low-degree monomial-fibre Newton refinement transport

**PROVED CANDIDATE + VERIFIED-EXACT / UNDER INDEPENDENT AUDIT.**

THM-3112 decomposes each individual cycle-weighted gap into positive labelled
Young-subgroup gaps, but its degree-five hostile proves that the signed
product-Gamma bank is not positive term by term in that decomposition.  The
missing operation is to forget labels only after grouping equal fibre sizes,
and then transport negative mass upward in the refinement order.  That
operation succeeds globally in the first four nonzero degrees.

This theorem strengthens the scalar Schur checks in THM-3110: those prove
`Phi_j(s_lambda)>=0`, whereas the operator below proves that the normalized
ratio is minimized by the row partition.  Its degree range is deliberately
smaller.

## 1. Uniform Young-fibre gaps

Fix `N>=1`.  For a set partition `pi` of `[N]`, put

```text
H_pi=product_(B in pi) Sym(B),
P_pi=|H_pi|^(-1) sum_(sigma in H_pi) sigma,
L_pi=identity-P_pi.                                             (1)
```

If `mu` is an integer partition of `N`, let `Part_mu` be the set of set
partitions whose block-size multiset is `mu`, and define

```text
Lbar_mu=|Part_mu|^(-1) sum_(pi in Part_mu) L_pi.                 (2)
```

For every finite labelled nonnegative alphabet `S`, including zero and
repeated letters, the THM-3112 gap has the exact monomial-fibre expansion

```text
beta_N(S)=sum_(mu partition N) m_mu(S) Lbar_mu.                  (3)
```

### Proof of `(3)`

A colouring `f:[N]->S` induces the set partition of `[N]` into its nonempty
colour fibres.  Fix `pi in Part_mu`.  Summing the THM-3112 coefficient `a_f`
over the injective assignments of colours to the blocks of `pi` gives

```text
[product_i mu_i!/N!] [product_r mult_r(mu)!] m_mu(S).            (4)
```

There are

```text
N!/[product_i mu_i! product_r mult_r(mu)!]                      (5)
```

set partitions of type `mu`.  Multiplying `(4)` and `(5)` gives exactly
`m_mu(S)`, proving `(3)`.  Zero letters merely kill the corresponding
injective assignments.

## 2. The refinement order survives averaging

Say that `nu` coarsens `mu` if the parts of `mu` can be grouped into sums
whose multiset is `nu`.  Then

```text
nu coarsens mu  ==>  Lbar_mu <= Lbar_nu.                         (6)
```

Indeed, form the bipartite incidence graph of pairs `(pi,tau)` with
`pi in Part_mu`, `tau in Part_nu`, and `pi` refining `tau`.  The symmetric
group is transitive on each vertex class and preserves incidence, so this
graph is biregular.  A uniformly chosen edge therefore has uniform marginals.
For every edge, `H_pi subset H_tau`; hence THM-3112 gives
`L_pi<=L_tau`.  Averaging this pointwise comparison over the uniform edge
coupling proves `(6)`.

The boundary type `(1^N)` is special:

```text
Lbar_(1^N)=0,                                                    (7)
```

because its Young subgroup is trivial.

## 3. The normalized product-Gamma coefficient vector

Use either residual signed-alphabet bank `Phi_j` (`j=1,2`) from THM-3110,
after its exact common-root deletion.  Let `Q_j` be the distinguished
Ferrers-dominant residual alphabet, `c_Q>0` its collected coefficient, and
let

```text
base_1=a^4 b^2 (b-a)^3,
base_2=a^3 b^2 (b-a)^4.                                        (8)
```

For `mu partition N`, define

```text
G_(j,N,mu)=
 [Phi_j(h_N) m_mu(Q_j)-Phi_j(m_mu) h_N(Q_j)]/base_j.             (9)
```

If `D_(j,N)` is the row-normalized central operator of THM-3112, then

```text
[h_N(Q_j)/base_j] D_(j,N)
   =sum_(mu partition N) G_(j,N,mu) Lbar_mu,                    (10)

sum_mu G_(j,N,mu)=0.                                           (11)
```

Equation `(10)` follows by substituting `(3)` into THM-3112 equation `(13)`
and clearing `h_N(Q_j)`.  Equation `(11)` is the identity
`h_N=sum_mu m_mu`.  Since `base_j>0` and `h_N(Q_j)>0`, positivity of the
right side of `(10)` is equivalent to `D_(j,N)>=0`.

## 4. Sharp chamber degree

In either THM-3110 chamber, every residual root multiset is a finite union and
difference of intervals with affine endpoints in the chamber coordinates.
Its power sum of order `r` therefore has degree at most `r+1`.  The standard
set-partition expansion of `m_mu` in power sums gives

```text
degree m_mu(S_R) <= N+length(mu).                                (12)
```

The collision divisor `(8)` has degree nine.  THM-3110 proves it for every
homogeneous response coefficient; equivalently, one may expand an arbitrary
homogeneous symmetric function in the Schur basis and use linearity.  Applying
the same estimate inside the signed response yields

```text
degree [Phi_j(m_mu)/base_j] <= N+length(mu)-9,
degree [Phi_j(h_N)/base_j] <= 2N-9.                              (13)
```

Thus the coefficient in `(9)` has degree at most
`3N+length(mu)-9`.  Every nonzero gap in `(10)` has
`length(mu)<=N-1`, so its degree is at most `4N-10`.  The only omitted type
is `(1^N)`, whose gap is zero by `(7)`; moreover `(11)` expresses its
coefficient as the negative sum of all the others.  Therefore the complete
coefficient vector has the sharp safe bound

```text
degree G_(j,N,mu) <= 4N-10.                                     (14)
```

This proves the degree bound before interpolation; it is not inferred from a
vanishing numerical shell.

## 5. Exact Newton-refinement certificate

Use the nonnegative integer chamber coordinates

```text
wide:   a=A+1,       b=2a+D,
tight:  a=A+D+2,     b=A+2D+3.                                 (15)
```

They cover, respectively, every integer `b>=2a` and every integer
`a<b<2a`.  For `N=5,6,7,8`, expand each coefficient in `(9)` as

```text
G_(j,N,mu)(A,D)
 =sum_(r+s<=4N-10) c_(mu;r,s) binom(A,r)binom(D,s).              (16)
```

For every fixed Newton slot `(r,s)`, the exact coefficient vector
`c_mu=c_(mu;r,s)` has total mass zero.  The companions construct
nonnegative rational numbers `x_(nu,mu)` only on refinement edges
`nu coarsens mu`, such that

```text
sum_(nu coarsens mu) x_(nu,mu)=-c_mu        when c_mu<0,
sum_mu x_(nu,mu)=c_nu                       when c_nu>0.           (17)
```

The initially unused positive supply is routed to the singleton root
`(1^N)`.  This is possible because total mass is zero and every passing
nonzero slot has `c_(1^N)<=0`.  Equations `(6)--(7)` then give

```text
sum_mu c_mu Lbar_mu
 =sum_(nu,mu) x_(nu,mu)(Lbar_nu-Lbar_mu)
 >=0.                                                           (18)
```

Every transitive coarsening edge can be replaced by any saturated chain of
single block merges; the vector differences telescope.  Thus `c` is
literally the boundary of a nonnegative one-chain on the partition-lattice
Hasse diagram.  Equation `(18)` is a positive refinement holotopy, not only a
feasible max-flow certificate.

The exact certificate universe is

```text
N       degree bound      Newton slots per bank/chamber
5           10                         66
6           14                        120
7           18                        190
8           22                        276.                       (19)
```

Thus `(18)` is checked on `2,608` exact Newton coefficient vectors across
the two banks and two chambers.  Every rational transport is supported on
an actual coarsening edge.  One excess interpolation shell and three
off-grid points per case verify the polynomial reconstruction independently
of the transport pass.  No floating-point comparison is used.

Because every binomial factor in `(16)` is nonnegative at integer
`A,D>=0`, summing `(18)` proves

```text
D_(j,N)>=0
for j=1,2, N=5,6,7,8, and every integer 0<a<b.                  (20)
```

Equivalently, in every irreducible representation of `Sym_N`, the normalized
THM-3110 response is minimized by the row partition `(N)`.

## 6. Failure boundary and scope

The transport in `(17)` is load-bearing.  THM-3112 gives an exact negative
termwise coefficient at `N=5`, `(a,b)=(1,2)`, and fibre type
`(2,1,1,1)`.  Thus `(20)` cannot be obtained by declaring each fibre
coefficient nonnegative.  It is the averaged merge coupling `(6)` that pays
this debt.

This theorem does **not** prove any of the following:

1. the refinement transport for `N>=9`;
2. strict positivity on every nontrivial representation;
3. membership in the smaller conjugacy-averaged uniform-octopus template
   cone ruled out by THM-3112;
4. positivity of the full finite histogram bank in all degrees;
5. the product-Gamma conjecture or the Gaussian Moment Conjecture.

What is new is a quotient-stable positive mechanism: labelled fibre gaps may
cancel, but their uniform type averages admit a monotone stochastic merge.
That is the exact sidecar missing from the failed termwise proof.

## 7. Exact companion

Reproduce with

```bash
python 04-computation/gmc_monomial_fibre_newton_refinement_thm3115.py
python -O 04-computation/gmc_monomial_fibre_newton_refinement_thm3115.py
```

The primary companion independently rebuilds the `N=5,6,7` monomial
symmetric evaluations, coarsening relation, sharp Newton triangles, excess
shells, off-grid evaluations, and rational flows.  It hashes the complete
coefficient/edge certificate rather than relying on printed decimal signs.
The secondary companion independently rebuilds `N=8` by a power-sum and
augmented-monomial recurrence.  It additionally proves computationally that
hooks alone, two-row types alone, and their union all miss valid slots; the
full coarsening poset is genuinely used.

QED.
