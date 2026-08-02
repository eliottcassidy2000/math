---
id: THM-3119
title: "Factorial-normalized labelled deletion and Young-carrier order"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Factorial
  block weights uniquely conjugate the algebraic unweighted occupancy
  deletion to genuine size-biased same-label deletion.  Same-label deletion
  preserves partition refinement stochastically, while the simultaneously
  rescaled Young-subgroup carriers remain monotone under coarsening.  Raw
  unweighted deletion does not preserve the positive Young-gap cone; the
  first nontrivial failures occur in degree four.  This repairs a vertical
  typing error but does not propagate the product-Gamma theorem to all
  degrees.
source: root/multiscale-newton-flag/product-gamma-width3-2026-08-02
audit: >
  An independent hostile audit checked the factorial conjugacy coefficient by
  coefficient, global-scalar uniqueness, same-uniform-label stochastic
  coupling and boundary orientation, rescaled carrier order, monomial
  multiplicities, zero/repeated-letter typing, both degree-four sign
  hostiles, the exact cover/max-flow census, hashes, and no-induction scope.
  Fresh normal and optimized runs line-match the stored transcript.
depends_on:
  - THM-3112-cycle-weighted-young-subgroup-gap-and-uniform-octopus-boundary
related:
  - THM-3115-low-degree-monomial-fibre-newton-refinement-transport
script: 04-computation/gmc_factorial_labelled_deletion_thm3119.py
output: 05-knowledge/results/gmc_factorial_labelled_deletion_thm3119.out
script_sha256: 85ce009278d42e5fda50a10d4218a3f6b2d0a8f11d4b9747c92e017eca0578d2
output_sha256: 814d74dcf5be7324bf987b487554b81fc423f9e20df9a0db460fd49d3b822aaa
hash_basis: LF-normalized bytes
---

# THM-3119 -- factorial-normalized labelled deletion and Young-carrier order

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3115 realizes each low-degree normalized product-Gamma coefficient vector
as the boundary of a nonnegative one-chain on the integer-partition Hasse
diagram.  A tempting vertical induction in the degree deletes one unit from
each occupied block.  The algebraic recurrence is exact, but that deletion is
not the same as deleting one uniformly labelled atom and does not preserve the
positive Hasse cone.

This theorem gives the exact repair.  The required diagonal gauge is the
factorial block weight.  It is unique, it converts block deletion into genuine
same-label deletion, and it preserves the monotone Young carrier after the
dual rescaling.  A degree-four hostile shows that the gauge is load-bearing.

## 1. Two deletion operators

For an integer partition `lambda`, write `r_k(lambda)` for the multiplicity
of its part `k`, and let `lambda down k` be the partition obtained by replacing
one part `k` by `k-1` and deleting a resulting zero.  On the free type module
define

```text
B e_lambda=sum_(k>=1) r_k(lambda) e_(lambda down k),
A e_lambda=sum_(k>=1) k r_k(lambda) e_(lambda down k).           (1)
```

`B` chooses an occupied block with weight one.  `A` chooses one of the
labelled atoms, so its total outgoing mass in degree `N` is `N`.

Put

```text
w_lambda=product_i lambda_i!,       W e_lambda=w_lambda e_lambda. (2)
```

Then, between consecutive degrees,

```text
W_N^(-1) B W_(N+1)=A.                                      (3)
```

Indeed,

```text
w_lambda/w_(lambda down k)=k,                              (4)
```

so the coefficient of `e_(lambda down k)` on the left of `(3)` is exactly
`k r_k(lambda)`.

### Uniqueness

Suppose positive diagonal weights `v_lambda` satisfy `(3)` with `W` replaced
by `V` in every degree.  Comparing the `lambda -> lambda down k` entry gives

```text
v_lambda/v_(lambda down k)=k.                               (5)
```

Lowering all parts to zero, in any order, yields

```text
v_lambda=v_emptyset product_i lambda_i!.                    (6)
```

Thus `(2)` is unique up to one global positive scalar.  Allowing an
independent scalar in each degree merely inserts the corresponding scalar
ratio in `(3)`.

## 2. Same-label deletion preserves refinement

Normalize `A` by `N` on partitions of `N`.  It is a stochastic kernel: choose
a set partition uniformly from its block-size orbit, choose a label uniformly
from `[N]`, and delete that same label.  The resulting type is
`lambda down k` with probability `k r_k(lambda)/N`.

If `nu` coarsens `mu`, use the biregular incidence coupling of uniform set
partitions `pi in Part_mu` and `tau in Part_nu` with `pi` refining `tau`.
Independently choose one uniform label `x` and delete `x` from both.  Then

```text
pi\{x} refines tau\{x}                                      (7)
```

pointwise, while both marginals are the kernel just described.  Therefore

```text
nu coarsens mu
  ==> A e_nu-A e_mu is a nonnegative fine-to-coarse Hasse boundary. (8)
```

The statement holds for every comparable pair, not only Hasse covers.

## 3. The factorial Young carrier remains ordered

Let `Lbar_lambda` be the uniform Young-subgroup gap from THM-3115 and set

```text
Kbar_lambda=w_lambda Lbar_lambda.                             (9)
```

If `nu` properly coarsens `mu`, every block merge multiplies `w` by a binomial
coefficient at least two.  In particular `w_nu>=w_mu`.  THM-3115 gives
`Lbar_mu<=Lbar_nu`, so

```text
Kbar_nu-Kbar_mu
 =w_mu(Lbar_nu-Lbar_mu)+(w_nu-w_mu)Lbar_nu
 >=0.                                                        (10)
```

Thus the same diagonal gauge that repairs deletion does not destroy the
positive carrier order.

For a finite nonnegative alphabet `S`, define the occupancy vector

```text
M_N(S)=sum_(lambda partition N) m_lambda(S)e_lambda.           (11)
```

The monomial product recurrence is exactly

```text
B M_(N+1)(S)=p_1(S) M_N(S).                                  (12)
```

For a fixed target type `mu`, the terms on the left are obtained either by
adjoining a new singleton part or by increasing one existing part of `mu` by
one.  Their multiplicities are exactly the standard coefficients in
`p_1 m_mu`, proving `(12)` coefficient by coefficient.

Conjugating by `(3)` gives the genuine labelled-deletion law

```text
A[W_(N+1)^(-1)M_(N+1)(S)]
 =p_1(S)W_N^(-1)M_N(S).                                      (13)
```

At the same time the THM-3112 gap becomes

```text
beta_N(S)=sum_lambda [m_lambda(S)/w_lambda] Kbar_lambda.       (14)
```

Equations `(10)`, `(13)`, and `(14)` are the factorial-normalized vertical
carrier.  They are exact for zero and repeated alphabet letters as well.

## 4. Why the raw deletion is false

Raw block deletion has outgoing mass `length(lambda)`, not `N`.  Across any
proper Hasse cover `mu -> nu`, the total mass of

```text
B e_nu-B e_mu                                                   (15)
```

is `-1`, so it cannot be an ordinary zero-mass Hasse boundary.  More strongly,
it need not even define a positive Young-gap operator.

The first nontrivial failures occur in source degree four.  For the covers

```text
(3,1) -> (4),              (2,2) -> (4),                       (16)
```

one gets in target degree three

```text
B e_4-B e_(3,1)=-e_(2,1),
B e_4-B e_(2,2)= e_3-2e_(2,1).                                (17)
```

On the sign representation of `Sym_3`, both `Lbar_3` and `Lbar_(2,1)` have
scalar one: each relevant Young subgroup contains a transposition and its
signed average vanishes.  Both lines of `(17)` therefore have scalar `-1`.
This is an exact operator counterexample, not only a coefficient-mass warning.

In degrees two and three the negative mass defect is absorbed by the
singleton gap `Lbar_(1^n)=0`, leaving a zero or positive operator.  Degree
four is the first visible obstruction.

## 5. Exact cover census

The companion exhausts every type and every Hasse cover in source degrees
two through nine.  It checks `(3)` entry by entry and verifies `(8)` by exact
integer max flow on the target refinement poset.  The universe is

```text
source degrees                         2,...,9
integer partitions                    95
distinct lowering entries            186
distinct Hasse covers                 182
labelled-deletion equality covers       1
minimum proper-merge weight ratio       2.                     (18)
```

All `182` same-label deletion comparisons pass.  The sole equality is the
degree-two cover `(1,1)->(2)`, whose deletion marginals both equal the unique
degree-one type.  The sign-representation raw census has `161` negative,
`14` zero, and `7` positive cover differences; hence the degree-four hostile
is the beginning of a pervasive obstruction, not an isolated accident.

## 6. Scope and remaining sidecar

This theorem repairs the vertical deletion typing, but it does not by itself
induct THM-3115 from degree `N` to `N+1`.  Diagonal normalization changes a
coefficient current `g_lambda` into `g_lambda/w_lambda`; an unscaled positive
Hasse divergence need not remain a positive divergence after that change.
The missing sidecar is factorial-normalized current compatibility, or an
explicit rooted/size-biased chain homotopy for the product-Gamma response
bank.

The theorem does not prove the row-normalized inequality in degree nine, the
all-degree product-Gamma conjecture, the Gaussian Moment Conjecture, SFC in
width three, NC2, LRC(14), JC(2), or DC(2).  It proves a universal deletion
and carrier lemma plus the sharp raw-deletion no-go.

## 7. Reproduction

Run

```text
python 04-computation/gmc_factorial_labelled_deletion_thm3119.py
python -O 04-computation/gmc_factorial_labelled_deletion_thm3119.py
```

and compare byte-for-byte with

```text
05-knowledge/results/gmc_factorial_labelled_deletion_thm3119.out.
```

The companion uses exact integers and rational numbers only.

QED.
