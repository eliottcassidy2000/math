---
id: THM-3127
title: "Partition-refinement Strassen upset dual and filter response"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On every
  finite partition-refinement poset, a zero-mass current is the boundary of
  a nonnegative fine-to-coarse Hasse flow iff it is nonnegative on every
  coarsening upset, equivalently on every isotone test function, equivalently
  it admits a refinement-supported Strassen coupling.  For the normalized
  product-Gamma current these cuts are exactly cross-multiplied monomial-filter
  response inequalities.  Principal upsets alone are insufficient already
  in degree four.
source: root/multiscale-newton-flag/2026-08-02
audit: >
  Two independent hostile audits rederived the max-flow cut orientation and
  capacity M+c(U), Strassen path decomposition, rational/integral typing,
  isotone layer-cake, Young-carrier implication and nonconverse, normalized
  filter identity including the zero-denominator boundary, coefficientwise
  Newton equivalence, and the nonprincipal-upset hostile.  Fresh normal and
  optimized runs byte-match the stored transcript; all load-bearing checks
  remain active under optimization and the declared LF hashes match.
depends_on:
  - THM-3115-low-degree-monomial-fibre-newton-refinement-transport
related:
  - THM-3119-factorial-normalized-labelled-deletion-and-young-carrier-order
  - THM-3120-row-pole-prefix-newton-flag-positivity
  - THM-3122-labelled-deletion-positive-kernel-ghost-and-no-upward-induction
script: 04-computation/gmc_partition_refinement_strassen_thm3127.py
output: 05-knowledge/results/gmc_partition_refinement_strassen_thm3127.out
script_sha256: 355324fd3289d2226d746e059c25968d6f8c0145eab157b7a6cb87873d660796
output_sha256: 08d4dbc2db6366695c9e89b46c4ad5f3f223d4825b94e6230d3f21e1433b3759
hash_basis: LF-normalized bytes
---

# THM-3127 -- partition-refinement Strassen upset dual and filter response

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3115 certifies thousands of normalized product-Gamma coefficient vectors
by max flow on the integer-partition refinement order.  The max flow is not
merely an implementation detail.  Its exact dual is a family of symmetric
function filters, and a failed flow has a canonical coarsening-upset witness.
This theorem identifies that dual and separates three notions which had been
easy to conflate: positive refinement transport, positivity against the
particular Young carriers, and pointwise positivity of the final scalar row.

## 1. Refinement currents

Let `P` be a finite poset.  Write `u lessdot v` when `v` covers `u`, and
orient every Hasse edge upward.  A vector

```text
c=sum_(v in P)c_v e_v,                  sum_v c_v=0,           (1)
```

is a **nonnegative upward Hasse boundary** if there are `x_(u,v)>=0` with

```text
c=sum_(u lessdot v)x_(u,v)(e_v-e_u).                            (2)
```

For the application, `P=P_N` is the set of integer partitions of `N`, with

```text
mu <= nu  iff  nu is obtained by merging parts of mu.          (3)
```

Thus upward means fine-to-coarse.  A subset `U subset P` is an upset if
`u in U` and `u<=v` imply `v in U`.  A function `f:P->R` is isotone if
`u<=v` implies `f(u)<=f(v)`.

## 2. Finite Strassen--Hasse duality

For every zero-mass vector `(1)`, the following are equivalent.

1. `c` is a nonnegative upward Hasse boundary.
2. Every isotone function has nonnegative pairing:

   ```text
   sum_v c_v f(v)>=0.                                          (4)
   ```

3. Every upset has nonnegative mass:

   ```text
   c(U):=sum_(v in U)c_v>=0.                                   (5)
   ```

4. If `c=c^+-c^-` and `M=sum c^+=sum c^->0`, there is a coupling
   `gamma(u,v)>=0`, supported on `u<=v`, whose first marginal is `c^-`
   and whose second marginal is `c^+`.

When `c=0`, item 4 is interpreted using the zero coupling.  If `c` is
rational, `x` and `gamma` may be chosen rational.  If `c` is integral, they
may be chosen integral.

### Proof

Equation `(2)` implies `(5)` edge by edge.  An upward edge contributes zero
to `c(U)` unless it crosses into `U`, in which case it contributes its
nonnegative weight.

Conversely, put `M=sum c^-` and form a network with source `s`, sink `t`, and
the vertices of `P`.  Give the edges

```text
s -> v       capacity -c_v       when c_v<0,
u -> v       capacity M          when u lessdot v,
v -> t       capacity  c_v       when c_v>0.                  (6)
```

A source-side cut which crosses a Hasse edge already has capacity at least
`M`.  Otherwise its poset vertices form an upset `U`, and its capacity is

```text
M+c(U).                                                        (7)
```

Condition `(5)` therefore says every cut has capacity at least `M`.  Max
flow/min cut gives a flow of value `M`; its restriction to the Hasse edges
has divergence `c`, proving `(2)`.  Integral capacities give an integral
flow, and clearing denominators gives the rational assertion.

Path-decomposing a Hasse flow gives the coupling in item 4.  Conversely,
route each coupled pair along any saturated chain from `u` to `v`; internal
vertices cancel.  This proves `1 iff 4`.

Finally, an upset indicator is isotone, so `(4)` implies `(5)`.  If an
isotone `f` has distinct values `a_1<...<a_k`, then

```text
f=a_1 1_P+sum_(j=2)^k(a_j-a_(j-1))1_{f>=a_j}.                 (8)
```

Every superlevel set in `(8)` is an upset.  The constant term pairs to zero
by `(1)`, while the other coefficients are positive.  Hence `(5)` implies
`(4)`, completing the proof.

## 3. Consequence for monotone Young carriers

Let `T_v` be self-adjoint operators satisfying

```text
u<=v  ==>  T_u<=T_v.                                          (9)
```

Then any of the equivalent conditions above gives

```text
sum_v c_v T_v
 =sum_(u lessdot v)x_(u,v)(T_v-T_u)>=0.                       (10)
```

Taking `T_mu=Lbar_mu`, the uniform Young-subgroup gaps of THM-3115, recovers
its positive refinement holotopy.  Taking the factorial-rescaled carriers
from THM-3119 gives the same implication in that gauge.

The converse to `(10)` is not asserted.  A particular carrier family is a
quotient of the full isotone dual and can erase a nonzero current.  Refinement
transport is therefore stronger than positivity of one operator or one
irreducible scalar.

## 4. The normalized monomial-filter criterion

Fix `N`.  Let `Phi` be any real signed response functional on homogeneous
symmetric functions of degree `N`, let `Q` be a distinguished nonnegative
alphabet with `h_N(Q)>0`, and let `b>0`.  Put

```text
q_mu=m_mu(Q),                 p_mu=Phi(m_mu),
H_Q=h_N(Q)=sum_mu q_mu,       H_Phi=Phi(h_N)=sum_mu p_mu,
G_mu=[H_Phi q_mu-p_mu H_Q]/b.                                (11)
```

Then `sum_mu G_mu=0`.  For an upset `U subset P_N`, define the monomial
filter

```text
F_U=sum_(mu in U)m_mu.                                        (12)
```

Summing `(11)` over `U` gives the exact identity

```text
G(U)=
[Phi(h_N)F_U(Q)-Phi(F_U)h_N(Q)]/b.                            (13)
```

Consequently `G` has a nonnegative fine-to-coarse Hasse transport if and
only if every coarsening filter obeys the cross-multiplied inequality

```text
Phi(h_N)F_U(Q)-Phi(F_U)h_N(Q)>=0.                             (14)
```

When `F_U(Q)>0`, this is the normalized response comparison

```text
Phi(F_U)/F_U(Q) <= Phi(h_N)/h_N(Q).                           (15)
```

The cross-multiplied form `(14)` is the canonical statement.  It remains
valid when `F_U(Q)=0`, where division in `(15)` would be illegal.

For the two THM-3115 product-Gamma banks, substitute `Phi=Phi_j`, `Q=Q_j`,
and `b=base_j`.  Equation `(13)` is exactly the sum over `U` of the
coefficient vector in THM-3115 equation `(9)`.  Thus a failed refinement
flow is not an anonymous max-flow failure: a minimum cut returns an explicit
symmetric filter `F_U` violating `(14)`.

## 5. Coefficientwise Newton form

Suppose the current `(11)` depends polynomially on chamber variables `z` and
is expanded in a basis

```text
G_mu(z)=sum_alpha c_(mu,alpha) B_alpha(z).                    (16)
```

By linearity, the coefficient vector `c_(.,alpha)` is a nonnegative Hasse
boundary if and only if, for every upset `U`,

```text
[B_alpha]
 { [Phi(h_N)F_U(Q)-Phi(F_U)h_N(Q)]/b } >=0.                   (17)
```

This is an exact equivalence for each slot.  In the wide and tight chambers
of THM-3115,

```text
B_(r,s)(A,D)=binom(A,r)binom(D,s)>=0                           (18)
```

on every integer chamber point.  Hence coefficientwise filter inequalities
imply a refinement flow at every support after summing the slots.  The
`4,120` exact flows in THM-3115 are equivalently `4,120` simultaneous families
of inequalities `(17)`.

This reformulation is useful beyond verification: an all-degree proof may
target the normalized responses of the upper filters directly, without
constructing a separate transport path for every Newton slot.  It does not
reduce the filter family to principal upsets.

## 6. Principal filters are not enough

Degree four is the first branching partition poset.  In the order

```text
(1,1,1,1) < (2,1,1) < {(3,1),(2,2)} < (4),                  (19)
```

consider

```text
c=e_(4)-e_(3,1)-e_(2,2)+e_(2,1,1).                          (20)
```

Every principal upset has mass zero or one.  Nevertheless

```text
U={(4),(3,1),(2,2)}                                         (21)
```

is an upset and `c(U)=-1`.  Therefore `(20)` is not a nonnegative Hasse
boundary.  Tests only on filters `F_(up mu)` miss the obstruction.

The degree-three vector

```text
-e_(3)+3e_(2,1)-2e_(1,1,1)                                 (22)
```

gives a second sharp orientation control: its coarsest upset has mass `-1`.

## 7. Exact companion and scope

Run

```text
python 04-computation/gmc_partition_refinement_strassen_thm3127.py
python -O 04-computation/gmc_partition_refinement_strassen_thm3127.py
```

and compare byte-for-byte with

```text
05-knowledge/results/gmc_partition_refinement_strassen_thm3127.out.
```

The companion exhausts all bounded zero-mass integer vectors in the first
four partition-poset sizes used there, comparing integral max flow with every
upset inequality.  It verifies `(20)--(22)` and checks `(13)` on a genuine
signed alphabet bank by exact integer monomial-symmetric evaluation.  Normal
and optimized runs must match the stored transcript.

This theorem is a duality and separation theorem, not a new sign proof for
the product-Gamma bank.  It does not prove `(14)` in degree ten or above,
positivity of arbitrary pole flags, arbitrary anchored width-three
product-Gamma, SFC, the Gaussian Moment Conjecture, LRC(14), JC(2), or DC(2).
Its exact contribution is to identify the missing object: a violating degree
has a coarsening filter, while a passing degree admits a monotone Strassen
coupling and a positive refinement holotopy.

QED.
