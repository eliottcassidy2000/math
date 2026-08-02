---
id: THM-3086
title: "Arbitrary cluster-composition chambers and alternant clutch holotopy"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT REQUESTED.  Every
  ordered composition of the remote slots into fixed-gap clusters gives an
  explicit nonempty simultaneous scale chamber above an arbitrary physical
  three-slot base.  The final carrier has a closed factorial ledger, every
  individual support prefix is positive, and each generalized alternant
  survives with total determinant exponent K!.  A separate multiscale
  alternant-clutch theorem identifies the symbolic gluing law, but its
  growing-gap physical realization is not asserted here.
source: root-gmc-cluster-composition-2026-08-01
depends_on:
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
  - THM-3069-one-normal-remote-terminal-suspension-and-physical-tropical-flag
  - THM-3073-upper-triangular-resultant-norm-and-torsion-character-death-barcode
  - THM-3082-admissible-suspension-word-simultaneous-chambers-and-scale-tree-holotopy
  - THM-3085-multi-normal-fixed-gap-cluster-and-unconditional-all-width-tail
related:
  - THM-2985-multiparameter-normal-map-and-arc-factor-separation
  - THM-3047-formal-corner-width-product-gamma-moment-and-strict-hankel-positivity
  - THM-3077-weighted-norm-augmentation-isogeny-and-kummer-carry-reconstruction
script: 04-computation/gmc_cluster_composition_chambers_thm3086.py
output: 05-knowledge/results/gmc_cluster_composition_chambers_thm3086.out
script_sha256: 2e560ec516423163be94126ff9a009ff21ce49bc415649f9abf905f39f7931ff
output_sha256: e2f1f47afb4f9383e0a334cb08befd3b92af870267842a1df02cfd09c8fbb793
hash_basis: LF-normalized bytes
---

# THM-3086 -- arbitrary cluster-composition chambers

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT REQUESTED.**

THM-3082 allows one-point and two-point suspension words.  THM-3085 supplies
one cluster of arbitrary fixed size.  Combining the exact quotient before
the entropy estimate shows that **every** ordered composition of the remote
slots is physical: each part is one fixed-gap cluster and each cut is one
genuine separation of scales.  The tropical exponent ledger forgets the
composition, but the generalized alternants do not.  They are the clutch
sidecar of the resulting scale-tree holotopy.

## 1. Cluster-composition supports

Fix

```text
n>=1,                  0<=a<b<c,                         (1)
```

and let `S_3=R_3(a,b,c)>0`, as supplied by THM-2824.  Fix a final width
`K>=4` and an ordered composition

```text
q_1+...+q_s=K-3,       q_i>=1.                           (2)
```

Put

```text
m_0=3,                 m_i=m_(i-1)+q_i,                 (3)
```

so node `i` has child width `m_(i-1)` and parent width `m_i`.  At that node
fix distinct integer gaps

```text
0=h_(i,1)<...<h_(i,q_i).                                (4)
```

Choose scales

```text
0<sigma_1<...<sigma_s,
L_i(T)=sigma_i T+O(1) in Z,                              (5)
```

and append the cluster

```text
L_i(T)+h_(i,1),...,L_i(T)+h_(i,q_i).                    (6)
```

All fixed data stay fixed as `T->infinity`.  Set `sigma_0=0` and

```text
delta_i=sigma_(i-1)/sigma_i.                             (7)
```

The clusters and the base are strictly ordered for large `T`.

For a parent width `k`, retain THM-3082's full factorial bill

```text
B_k=(k-1)k!,
J_k(delta)=B_k delta[log(k/delta)+1],
J_k(0)=0.                                                (8)
```

The local scalar gaps are

```text
Gamma_i=
  log[m_i^m_i/(m_i-1)^(m_i-1)]             if q_i=1,
  m_(i-1)log[(m_(i-1)+1)/m_(i-1)]           if q_i>=2.   (9)
```

The first line is THM-3069's sharp one-normal resultant-Newton gap.  The
second is `-log rho_m` for THM-3085's coefficientwise multi-normal gap
`rho_m=(m/(m+1))^m`.

Assume at every node

```text
0<=delta_i<1/e,
J_(m_i)(delta_i)<Gamma_i.                               (10)
```

Every gate has a nonempty interval at zero because `J_k(0)=0<Gamma_i`.
Thus `(10)` defines a nonempty relative-open simultaneous scale chamber for
every ordered composition `(2)`.

## 2. Exact node quotient and the local gap

At node `i`, pivot the determinant-one normal coordinates on the original
fixed offset `c`, never on a moving cluster point.  Put

```text
p_i=m_(i-1)+1.                                          (11)
```

For `p_i<=r<=m_i`, let `U_(r,i)>0` be THM-3085's exact physical high
carrier at `L_i(T)`.  Let `E_i(T)` be the exact intrinsic `q_i`-variable
all-high normal resultant.  Its limiting row forms are

```text
(-r^h_(i,q_i),
 r^h_(i,1)-r^h_(i,q_i),...,
 r^h_(i,q_i-1)-r^h_(i,q_i)),             p_i<=r<=m_i.  (12)
```

Writing

```text
M_i=(r^h_(i,j))_(p_i<=r<=m_i,1<=j<=q_i),               (13)
```

the generalized-Vandermonde theorem gives `det(M_i)>0` and

```text
E_i(T)=(-1)^[q_i m_i!/m_(i-1)!]
       det(M_i)^[m_i!/m_(i-1)!][1+O(T^-1)].             (14)
```

For `q_i=1`, the normalized all-high form is exactly `(-z)^m_i`, so its
powered contribution below is one.  For `q_i>=2`, the determinant exponent
in `(14)` is even and `E_i(T)>0` for large `T`.

Let `I_i` be the coefficient ideal generated by every lower form term that
contains a new normal variable.  THM-3073 and THM-3085 give the exact
associated-graded node quotient

```text
R_(m_i) mod I_i
 =R_(m_(i-1))^[m_i!/m_(i-1)!]
  E_i(T)^[m_(i-1)!]
  product_(r=p_i)^m_i U_(r,i)^[m_i!/r].                 (15)
```

At a one-point node, THM-3069 sharpens the complementary resultant bank to
the first gap in `(9)`.  At a multi-point node, THM-3085 bounds every
coefficient layer outside `(15)` by the second gap.  Turning on all child
scales can raise any outer resultant atom by at most `J_(m_i)(delta_i)` per
unit parent scale: this is exactly THM-3082's multislow entropy lemma and the
bill `(8)`.

The quotient is taken **before** this estimate.  Its grouped carrier has
nonnegative child-scale rate, so there is no factor two in the relative
loss.  The entire complementary bank is therefore

```text
O(poly(T) exp(-sigma_i[Gamma_i-J_(m_i)(delta_i)]T))      (16)
```

relative to `(15)`.  Strictness in `(10)` is load-bearing.

## 3. Closed simultaneous carrier

Define

```text
eta=min_i sigma_i[Gamma_i-J_(m_i)(delta_i)]>0.           (17)
```

Induction from the innermost node through `(15)--(16)` gives

```text
R_K(T)=S_3^[K!/6]
       product_(r=4)^K U_(r,node(r))^[K!/r]
       product_(i=1)^s E_i(T)^[K! m_(i-1)!/m_i!]
       [1+O(poly(T)e^(-eta T))].                         (18)
```

The unique node containing degree `r` determines `node(r)`.  Every factor
in `(18)` is positive after its displayed exponent for large `T`; hence the
physical first-window resultant is positive throughout every chamber
`(10)` for all sufficiently large `T`.

The factorial recursion is exact:

```text
S_3: K!/6,          U_r: K!/r,
E_i: K! m_(i-1)!/m_i!.                                  (19)
```

Moreover `(14)` shows that the determinant `det(M_i)` at **every** node has
total exponent

```text
[m_i!/m_(i-1)!][K!m_(i-1)!/m_i!]=K!.                   (20)
```

Thus the tropical `S_3,U_r` ledger forgets the ordered composition, while
the exact carrier remembers every cluster through one generalized
alternant, all with the same final exponent `K!`.

There are

```text
2^(K-4)                                                   (21)
```

ordered compositions of `K-3`, hence that many sufficient
cluster-composition sectors.  They are not claimed disjoint or maximal.

## 4. Every internal support prefix is positive

Suppose a node has child width `m` and cluster size `q>=2`.  Its endpoint
condition is

```text
J_(m+q)(delta)<m log[(m+1)/m].                           (22)
```

For a partial cluster of size `2<=j<=q`, the THM-3085 gap is the same right
side, while

```text
J_(m+j)(delta)<=J_(m+q)(delta).                          (23)
```

For the first point, THM-3069 has the strictly larger gap

```text
log[(m+1)^(m+1)/m^m]
 >m log[(m+1)/m].                                        (24)
```

Equations `(22)--(24)` therefore control every internal prefix without an
extra chamber inequality.  The theorem proves one nested physical flag of
positive resultants through **all** widths `4,...,K`, not merely positivity
at cluster endpoints.

## 5. A two-cluster eight-slot chamber

Take the composition `(2,3)`: append a fixed-gap pair at scale `D`, then a
fixed-gap triple at scale `C`.  The only nonzero ratio is `delta=D/C`.  At
the second node the aperture is bracketed rigorously by

```text
0.00000017320136922 < epsilon_(5,8)
                    < 0.00000017320136924.               (25)
```

Thus `D/C=1/10000000` is safe.  Formula `(18)` is

```text
R_8=S_3^6720 E_(45,D)^2016 E_(678,C)^120
    U_(4,D)^10080 U_(5,D)^8064
    U_(6,C)^6720 U_(7,C)^5760 U_(8,C)^5040
    [1+O(poly(T)e^(-eta T))]>0.                          (26)
```

Both alternant determinants in `(26)` occur to total exponent `8!=40320`,
as predicted by `(20)`.

## 6. Alternant clutch theorem

The sidecar in `(20)` has a general algebraic explanation.  Let

```text
0<x_1<...<x_q,
0<=A_1<...<A_s in Z,
0<=b_(t,1)<...<b_(t,q_t),       sum q_t=q,               (27)
```

and, for positive integer `H`, form the `q` exponent columns

```text
h_(t,j)(H)=A_t H+b_(t,j).                                (28)
```

Let `I_t` be the consecutive block of `q_t` rows paired with block `t`, in
increasing order.  Then

```text
det[x_r^h_(t,j)(H)]
 =product_t [(product_(r in I_t)x_r)^(A_t H)
              det(x_r^b_(t,j))_(r in I_t,1<=j<=q_t)]
  [1+O(e^(-epsilon H))],                                (29)
```

for some `epsilon>0`.  Every local determinant in `(29)` is a positive
generalized alternant.

Indeed, a determinant permutation has `H`-rate

```text
sum_(t,j) A_t log x_(pi(t,j)).                           (30)
```

Strict rearrangement uniquely assigns the consecutive largest row block to
the largest `A_t`, and so on.  Permutations inside each assigned block are
tied.  Summing their signed coefficients is exactly the product of the local
determinants in `(29)`; row and column blocks have the same order, so no
extra sign appears.  The finite next-rate gap gives `epsilon`.

Formula `(29)` is the **alternant clutch law**.  It is valid for any number
of blocks and arbitrary strictly increasing within-block integer gaps.  It
explains how one large alternant breaks into the node sidecars of a scale
tree.  Here it is an algebraic symbol theorem only: if the physical gaps in
`(4)` themselves grow with `T`, uniform control of the factorial
coefficients and the condition number of the alternant is additionally
required.  No such growing-gap physical chamber is claimed.

## 7. Composition-cube holotopy and boundaries

An ordered composition is equivalently a choice of cuts in the `K-4`
internal remote-slot gaps.  Hence the labels form a Boolean composition
cube.  Within one label, the edge ratios `(7)` satisfying `(10)` form a
star-shaped rate chamber; sending an edge ratio to zero recovers the ordered
parenthesization.  The clutch law `(29)` records the sidecar carried across
refinements.  This is a rate-space scale-tree holotopy, not a continuous
path through integer supports and not a wall-crossing theorem.

The proof uses arbitrary physical three-slot bases, but fixed final width,
fixed distinct gaps in each cluster, and sufficiently separated positive
scales.  Repeated gaps kill a generalized alternant.  A zero child
resultant kills the carrier.  Equality in `(10)` admits exact carrier
cancellation hostiles from THM-3082.  Equal-scale arbitrary supports,
growing cluster diameter, a width-uniform aperture, maximality of the
composition fan, arbitrary-radial GMC(2), NC2, LRC(14), JC(2), and DC(2)
remain outside the theorem.

The relation/magnitude/carry distinction is exact here.  The alternants and
composition give the symbolic relation state; `(10)` supplies the strict
magnitude margin; positivity of the physical `U,E` carrier trivializes the
phase carry.  In an abstract phased norm problem, THM-3077's augmentation
and continuation sidecar would still be necessary.

## 8. Exact evidence

The companion verifies:

1. `7,768` nonsurviving physical coefficient cells and `288` declared
   survivors for child widths three through ten and cluster sizes through
   four;
2. all `511` ordered-composition ledgers through width twelve, including all
   `256` width-twelve sectors, every exponent in `(19)`, and the determinant
   identity `(20)`;
3. the exact eight-slot carrier `(26)`, the safe ratio `10^-7`, and the
   rational threshold bracket `(25)`;
4. fifty-six two-block unique-row-set and local-sign controls; and
5. `360` full multiblock tied-permutation sums: every nontrivial ordered
   composition through rank seven and three within-block gap patterns gives
   exactly the positive product of its local generalized alternants.

Run

```text
python 04-computation/gmc_cluster_composition_chambers_thm3086.py
python -O 04-computation/gmc_cluster_composition_chambers_thm3086.py
```

Both modes must equal the stored transcript after LF normalization.

**QED, pending independent audit and status promotion.**
