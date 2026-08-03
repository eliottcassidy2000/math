---
id: THM-3209
title: "Depth-eight complete quotient reset and negative-singleton tangent"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For support (1,3), bank I2, the dominant quotient alphabet is exactly
  sigma0=(1,3,3,4,5,6,7,8).  Removing this complete alphabet makes every
  positive-degree partition response vanish, so delta_sigma0 is an
  all-degree boundary selector law.  It is the unique complete degree-five
  response zero through physical depth eight.  Every legal one-extra-pole
  tangent state has negative degree-five coarsest-upset response, so the
  reset has no feasible first-order ray inside that four-state star.
source: root/multiscale-newton-flag/product-gamma-width3/2026-08-02
audit: >
  The exact companion pins both transitive response helpers and THM-3184,
  reconstructs all 678 depth-eight states and the ten-row cocircuit, checks
  its 15/662/1 sign census and unique zero, separately checks the complete
  degree-five response on all 2,498 states through depth eight, verifies all
  25 one-exchange neighbours, and replays the negative-singleton formula on
  all partitions through degree fourteen.  Normal, optimized, and stored
  replay agree exactly.  An independent immutable audit rederived the
  all-degree empty-alphabet reset, uniqueness implication, one-exchange
  hostile, negative-singleton specialization, and local-star separation;
  replayed both modes; matched both hashes; and accepted the local-versus-
  global depth-nine boundary.
depends_on:
  - THM-3115-low-degree-monomial-fibre-newton-refinement-transport
  - THM-3127-partition-refinement-strassen-upset-dual-and-filter-response
  - THM-3184-depth-seven-degree-fourteen-farkas-death
related:
  - THM-3160-complete-pluecker-pole-holotopy-and-selector-projection-no-go
  - THM-3182-factorial-gauss-manin-rank-one-reset-and-two-transverse-smith-bands
  - THM-3185-iterated-factorial-frobenius-descent-and-witt-carry-reset-hierarchy
script: 04-computation/gmc_depth_eight_complete_quotient_reset_thm3209.py
output: 05-knowledge/results/gmc_depth_eight_complete_quotient_reset_thm3209.out
script_sha256: c04b80a5e4a9cbfc4d7a59fac338147f29582a36751853e0da98752a83288fe8
output_sha256: cb8c20f1a8fffd1b594cf09283ee87269b71638cb8df089db8c74c147717345f
hash_basis: LF-normalized bytes
---

# THM-3209 -- depth-eight complete quotient reset and negative-singleton tangent

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3184 proves that the complete support-`(1,3)`, bank-`I2` selector cone is
empty through physical prefix depth seven once degree fourteen is imposed.
At depth eight, its separating cocircuit meets one state exactly.  This is not
an accidental finite-degree zero: the state deletes the entire dominant
quotient alphabet and therefore annihilates every positive-degree response.

The resulting resurrection is a boundary point, not a strict selector.  Its
first local tangent is also rigid: adding any one still-available pole creates
a negative-singleton virtual alphabet whose degree-five coarsest upset points
strictly out of the selector cone.

## 1. Complete quotient-alphabet deletion

Retain THM-3184's reduced pole multiset

```text
P=(8,7,6,5,5,4,4,3,3,2,2,2,1,1,1,1).                 (1)
```

For the dominant support-`(1,3)`, bank-`I2` row, direct cancellation in the
reduced product gives the distinguished quotient alphabet

```text
Q=(1,3,3,4,5,6,7,8).                                  (2)
```

Write a physical prefix state as a submultiset `sigma<=P` and use lambda-ring
difference for the virtual quotient

```text
Q^sigma=Q-sigma.                                       (3)
```

The fixed-`Q` partition response inherited by THM-3184 has determinant form

```text
G_N^sigma(lambda)
 =Phi^sigma(h_N)m_lambda[Q^sigma]
  -Phi^sigma(m_lambda)h_N[Q^sigma],                    (4)
```

for every partition `lambda` of `N`.  Set

```text
sigma0=Q=(1,3,3,4,5,6,7,8).                            (5)
```

Then `Q^sigma0=0`.  Every homogeneous symmetric function of positive degree
vanishes at the empty alphabet, so for every `N>=1` and every partition
`lambda` of `N`,

```text
m_lambda[Q^sigma0]=h_N[Q^sigma0]=0,
G_N^sigma0(lambda)=0.                                  (6)
```

No property of the possibly nonzero numerator functionals
`Phi^sigma0(h_N)` and `Phi^sigma0(m_lambda)` is needed.  The vanishing is an
exact quotient-side reset.

Consequently, for every finite `D>=5`,

```text
delta_sigma0 in C_D^(<=8),                             (7)
```

and every partition-coarsening upset inequality holds with equality.  Thus
depth eight resurrects weak feasibility after THM-3184's depth-seven death,
simultaneously in all degrees.  It does not give strict feasibility.

## 2. Uniqueness through depth eight

The ten-row normalized cocircuit `H` from THM-3184 is strict on all `1,820`
states of depths one through seven.  Exact evaluation on the `678` states of
depth eight gives

```text
#(H>0)=15,       #(H<0)=662,       #(H=0)=1,            (8)
```

with

```text
{sigma: |sigma|=8 and H(sigma)=0}={sigma0}.             (9)
```

Thus the old cocircuit ceases to separate at depth eight but has a unique
contact point.  Independently, exact evaluation of the complete degree-five
response on all

```text
1,820+678=2,498
```

states through depth eight gives

```text
sigma0 is the unique state in S_<=8
whose complete degree-five response is zero.            (10)
```

In particular `sigma0` is the unique all-degree reset state through depth
eight.  This is uniqueness inside the fixed physical prefix bank, not among
arbitrary virtual alphabets.

The reset is already isolated under one exchange.  There are exactly `25`
legal depth-eight states obtained by deleting one entry of `sigma0` and
inserting a different available pole.  Every one has nonzero `(5)` response.
The smallest absolute value is attained at

```text
sigma1=(2,3,3,4,5,6,7,8).                              (11)
```

Here `Q^sigma1={1}-{2}`, and exact evaluation gives

```text
Phi^sigma1(h_5)=1440,     h_5[{1}-{2}]=-1,
Phi^sigma1(p_5)=0,        p_5[{1}-{2}]=1^5-2^5=-31,
G_5^sigma1((5))=-44640.                                (12)
```

This is a sharp local hostile, not a claim that one exchange is the only
possible route away from the reset.

## 3. The negative-singleton tangent star

After deleting `sigma0`, the unused pole multiset is

```text
P-sigma0=(1,1,1,2,2,2,4,5).                            (13)
```

Hence the four distinct legal one-extra-pole states are

```text
sigma_r=sigma0+{r},             r in {1,2,4,5}.         (14)
```

Their quotient virtual alphabet is `Q^sigma_r=-{r}`.  If a partition `lambda`
of `N` has length `ell` and `a_j` parts equal to `j`, the monomial recurrence,
or equivalently expansion in power sums `p_k[-r]=-r^k`, gives

```text
m_lambda[-r]
 =(-1)^ell ell!/(product_j a_j!) r^N.                  (15)
```

Also

```text
h_N[-r]=0                    for every N>=2.            (16)
```

Substitution in `(4)` therefore gives the universal tangent direction

```text
G_N^sigma_r(lambda)
 =Phi^sigma_r(h_N)
  (-1)^ell ell!/(product_j a_j!) r^N,       N>=2.       (17)
```

At degree five the numerator scalar is independent of the legal root:

```text
Phi^sigma_r(h_5)=1440,             r in {1,2,4,5}.      (18)
```

The principal coarsest upset `{(5)}` consequently has responses

```text
r=1: -1440,   r=2: -46080,
r=4: -1474560, r=5: -4500000.                           (19)
```

For any probability law `nu` supported on these four tangent states,

```text
R_(5,{(5)})(nu)=-1440 sum_r nu_r r^5<0.                (20)
```

Because `delta_sigma0` has zero response, `(20)` also says that every
nontrivial ray

```text
(1-epsilon)delta_sigma0+epsilon nu,       0<epsilon<=1, (21)
```

leaves the selector cone already at degree five.  Thus no mixture supported
on the legal one-extra-pole star can produce strict resurrection through
degree fourteen; it is not even weakly feasible at degree five.

Equivalently, if

```text
Delta_star=conv({delta_sigma0} union
                {delta_sigma_r:r in {1,2,4,5}}),
```

then for every `D>=5`,

```text
C_D^(<=9) intersect Delta_star={delta_sigma0}.           (22)
```

This is a local tangent statement.  It does not exclude a depth-nine law
using states outside the star, whose responses can compensate `(20)`.

## 4. Holotopy interpretation and cross-frontier analogy

The depth-degree bifiltration now has two different kinds of event:

```text
depth 7, degree 14:     strict cocircuit death;
depth 8, every degree:  one zero-stratum reset.          (23)
```

The new point is not an interior chamber.  It is an isolated contact of the
old cocircuit and all partition-upset walls, and its visible outgoing star is
oriented outward by one common degree-five facet.  Any strict continuation
must therefore arrive nonlocally in the state graph.

There is a mechanism-level analogy with THM-3182 and THM-3185: a complete
factorial block resets a transfer state there, while deletion of the complete
quotient alphabet resets the partition response here.  The preserved idea is
full-block annihilation.  The objects, gauges, and carriers are different;
there is no map from the factorial Gauss--Manin lattice to this selector bank,
and no consequence is transported across that analogy.

## 5. Exact evidence and scope

Run

```text
python 04-computation/gmc_depth_eight_complete_quotient_reset_thm3209.py
python -O 04-computation/gmc_depth_eight_complete_quotient_reset_thm3209.py
```

and compare LF-normalized bytes with the declared output.  The companion uses
integer and `Fraction` arithmetic only.  It pins the two transitive response
helpers and THM-3184's immutable script/output, independently rebuilds the
ten upsets, enumerates all states through depth eight, checks `(8)`--`(12)`,
verifies the empty-alphabet identities through degree fourteen, and exhausts every
partition in the negative-singleton formula through degree fourteen.  The
all-degree quantifier in `(6)` and `(15)`--`(17)` is algebraic; the finite
replay is a regression witness, not the logical source of that quantifier.

This theorem concerns the fixed support-`(1,3)`, bank-`I2` averaged virtual
prefix model.  The atom `delta_sigma0` is a legal finite-bank law, but no
response-compatible stopping rule, nonlinear or Markov carrier, or
decomposition of the original product-Gamma current is supplied.  The result
does not prove strict selector feasibility at depth eight or nine, does not
exclude compensating unrelated depth-nine states, and proves no Gaussian
moment conjecture, `NC(2)`, `GMC(2)`, or `LRC(14)` conclusion.

QED.
