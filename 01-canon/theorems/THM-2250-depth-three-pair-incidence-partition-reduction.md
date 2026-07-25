---
id: THM-2250
title: "Depth-three pair-incidence partition reduction"
status: >
  PROVED + VERIFIED-EXACT + HOSTILE-AUDITED. In the remaining scalar
  valuation profile (3,4,5), a centered charge built from the whole
  three-danger OR clause reduces the signed residual to one nonnegative
  guard-weighted cylinder. If at least two normalized blocker cores are
  distinct, THM-1166's same-time danger-pair cap supplies a Lagrange tax.
  Exact co-adapted Bellman certificates close the all-distinct partition
  and each of the three one-equality partitions. Thus every (3,4,5)
  survivor must have all three normalized blocker cores equal. The
  all-equal branch remains open, and LRC(14) remains open.
source: codex-klein-2026-07-25-depth-three-pair-incidence
depends_on:
  - THM-1166-seven-wall-fano-gcd-discrepancy
  - THM-2198-scalar-five-plus-three-image-pump-and-first-depth-exclusion
  - THM-2222-scalar-transfer-parity-tower-and-four-checkpoint-survivor-reduction
related:
  - THM-2232-same-core-signed-eigen-markov-dual-exclusion
  - THM-2233-guard-danger-hidden-state-bellman-profile-exclusion
  - THM-2239-unrestricted-multicore-signed-dual-profile-exclusion
script: 04-computation/lrc14_depth_three_pair_incidence_thm2250.py
output: 05-knowledge/results/lrc14_depth_three_pair_incidence_thm2250.out
script_sha256: 55d8517d660dd4117e47fc181d3275c01593b752072f6eeff2e2ba077d098916
output_sha256: de1370371262a84500f071198e1d3a54dbc666f6ba0038c8f2317baf50b1afda
hash_basis: working-tree bytes (LF)
---

# THM-2250 -- pair incidence leaves only the all-equal core

THM-2239 reduces the scalar ledger to the first-depth-one rows and the
single depth-three profile `(3,4,5)`. This theorem resolves the equality
partition of the latter profile. It does not eliminate the profile itself:
the all-equal normalized-core branch is the exact residual.

## 1. Universe and the diagonal clause word

On `T=R/Z`, put

```text
D_a={x:||ax||<1/14},       C_H={x:||Hx||>1/7}.       (1)
```

Work in THM-2198's scalar five-unit/three-blocker branch with

```text
c_0=13^3 u_0,       c_1=13^4 u_1,       c_2=13^5 u_2, (2)
```

where every `u_j` is a positive thirteen-unit. Let

```text
R=1_(C_H)-sum_(i=1)^5 1_(D_(q_i)),

p=integral R_+=integral R_- >=delta_5:=961/6930.     (3)
```

As in THM-2239, set `Sx=13x mod 1` and

```text
G_t=1_(C_H) o S^t,
X_(j,t)=1_(D_(u_j)) o S^t.                           (4)
```

The profile `(3,4,5)` has a useful diagonal alignment. Define

```text
W_t=X_(0,t) OR X_(1,t+1) OR X_(2,t+2).              (5)
```

The checkpoint union `U_k` from THM-2222 is exactly

```text
U_k=W_(3-k),                         0<=k<=3.         (6)
```

Consequently the positive and negative supports of the same residual obey

```text
support(R_+) subset P:=W_1 intersection W_3,
support(R_-) subset N:=W_0 intersection W_2.         (7)
```

Both odd checkpoints are retained in `N`; no single-odd-clause projection
is made.

## 2. Center the whole OR clause

Put `rho=-1/13`. Since `W_2=W_0 o S^2`, transfer duality for the exact
eigenlaw `L^2R=R` gives

```text
integral R W_2=rho^2 integral R W_0.                 (8)
```

Define the nonnegative whole-clause charge and its complementary dual by

```text
C=W_2+rho^2(1-W_0)
 =W_2+(1/169)(1-W_0),

q=1-C.                                                (9)
```

Because `integral R=0`, equation (8) implies

```text
integral R C=integral R q=0.                         (10)
```

On `N`, both `W_0` and `W_2` equal one, hence

```text
C=1 and q=0                         on N.             (11)
```

Using (7), (10), and (11) in the signed identity gives

```text
p
 =integral [R_+(1-q)+R_-q]
 =integral R_+ C
 <=integral F,                                       (12)

F
 =1_(G_0)1_(W_1)1_(W_3)
    [W_2+(1/169)(1-W_0)].                            (13)
```

The last inequality uses the exact facts that `R_+` is the indicator of
the positive residual set, that this set lies in `G_0 intersection P`,
and that `C>=0`. Centering the composite OR clause is essential: centering
its three atoms separately double-counts overlaps.

## 3. The distinct-core incidence cap

If `u_a!=u_b`, THM-1166's exact danger-pair law gives, for every time `t`,

```text
integral X_(a,t)X_(b,t)
 =measure(D_(13^t u_a) intersection D_(13^t u_b))
 =measure(D_(u_a) intersection D_(u_b))
 <=1/14.                                             (14)
```

The middle equality is common-dilation invariance of Haar measure. The
last inequality can also be read directly from THM-1166's folded formula.
For a reduced coprime pair `a<b`,

```text
measure(D_a intersection D_b)
 =1/49+[fold_14(a+b)-fold_14(b-a)]/(196ab),

0<=fold_14(r)<=49.
```

If `ab>=5`, this is at most `1/49+1/(4ab)<1/14`. The only smaller
distinct reduced pairs `(1,2),(1,3),(1,4)` have exact intersections
`1/14,1/21,1/28`. Thus (14) includes its equality boundary and does not
depend on a numerical scout. The upper bound applies only to distinct
normalized cores. If the cores are equal, the intersection has mass
`1/7`, so substituting `1/14` there would be false.

There are exactly five equality partitions of three labelled cores:

```text
distinct,
u_0=u_1!=u_2,
u_0=u_2!=u_1,
u_1=u_2!=u_0,
u_0=u_1=u_2.                                         (15)
```

The theorem closes the first four.

## 4. Lagrange form of the causal bound

For one of the first four partitions, replace equal danger coordinates by
one group bit. A hidden state consists of the guard bit and the resulting
one, two, or three danger-group bits. Conditional on the exact future point,
root counting gives the coordinate marginals

```text
P(G_t=1 | G_(t+1)=g)=(10-g)/13,

P(X_(a,t)=1 | X_(a,t+1)=x)=(2-x)/13.                 (16)
```

The bits use the same root digit. Their joint current law is therefore
enlarged to every coupling with the marginals in (16); it is never replaced
by an independent product. At terminal time, the joint law is enlarged to
every coupling with marginals

```text
5/7, 1/7, ..., 1/7.                                  (17)
```

The actual root process is feasible in these coupling polytopes pointwise.

Choose nonnegative rational penalties `theta_(t,ab)` only for pairs of
distinct danger groups. Equation (14) gives

```text
integral F
 <=sum_(t,ab) theta_(t,ab)/14
   +sup E[
       F-sum_(t,ab)theta_(t,ab)X_(a,t)X_(b,t)
     ],                                              (18)
```

where the supremum ranges over the co-adapted coupling relaxation in
(16)--(17). The second term is an exact backward Bellman problem. Its state
retains

```text
time t in {5,...,0};
the two-bit mask for W_1,W_3;
the two-bit mask for W_0,W_2;
the full labelled future hidden bit vector.          (19)
```

Thus the same root digit, every equality relation, both parity carriers,
and the time-zero guard bit are retained. Pair incidence enters as a
same-time stage cost before the current root digit is forgotten.

## 5. Four exact certificates

For `distinct`, group labels are `0,1,2`. In a one-equality row, labels
`0,1` denote the two distinct groups. All omitted penalties are zero.

```text
partition   nonzero (time,pair):theta

distinct    (5,12):1/128,
            (4,01):1/256, (4,12):1/14,
            (3,01):7/64,  (3,12):1/2

u_0=u_1     (5,01):1/128, (4,01):1/8,  (3,01):1/2

u_0=u_2     (5,01):1/128, (4,01):3/32, (3,01):2/3

u_1=u_2     (4,01):1/256, (3,01):2/3.                (20)
```

The exact Bellman term, cap charge, total bound, and strict target gap are

```text
partition   Bellman term                    cap charge
            total bound                     delta_5-total

distinct    40924182149/787117397248         1241/25088
            159719273895/1574234794496
            28998933529879/779246223275520

u_0=u_1     15834701/270301304               81/1792
            897681961/8649641728
            149383550777/4281572655360

u_0=u_2     18763571/439239619                295/5376
            32921235043/337336027392
            2286565797041/55660444519680

u_1=u_2     7571445773/112445342464          515/10752
            11106308699/96381722112
            2609352610927/111320889039360.           (21)
```

Every final entry in (21) is positive. Equations (3), (12), (18), and
(21) contradict one another in each of the four partitions. Therefore a
surviving `(3,4,5)` profile must satisfy

```text
u_0=u_1=u_2.                                         (22)
```

## 6. Exact finite verification

For `d` hidden bits, including the guard, every local coupling LP has
`2^d` columns

```text
(1,b_1,...,b_d)^T
```

and rank `d+1`. The companion enumerates every invertible column basis,
keeps every nonnegative basic distribution, and deduplicates coincident
vertices. For the all-distinct case, `d=4`: there are `3008` invertible
bases, and the vertex counts for the sixteen future bit vectors followed
by the terminal marginal vector are

```text
(94,73,73,41,73,41,41,30,
 34,71,71,36,71,36,36,30,34).                       (23)
```

For every one-equality case, `d=3`: there are `58` invertible bases and

```text
(9,8,8,6,6,8,8,6,6)                                 (24)
```

vertices. Every Bellman maximum is first obtained by exact rational primal
vertex enumeration. For an optimal basis `B`, the companion independently
solves

```text
A_B^T y=c_B
```

and checks `A^T y>=c` on every column and equality of the primal and dual
values. The four certified runs make respectively

```text
309, 169, 137, 135
```

local requests and verify

```text
159, 79, 83, 77
```

distinct exact primal-dual optima. Any failed assertion aborts the script.

Reproduce with

```bash
python3 04-computation/lrc14_depth_three_pair_incidence_thm2250.py
```

## 7. Hostile controls and scope

With every pair penalty removed, the same exact relaxation returns

```text
distinct       8963699/62748517,
u_0=u_1        62694053/439239619,
u_0=u_2        62692901/439239619,
u_1=u_2        62691380/439239619,                   (25)
```

all strictly above `delta_5`. Thus the incidence caps, not equality
substitution alone, are load-bearing.

For the all-equal partition there is no legal distinct-pair cap. Its exact
hostile value is

```text
52940/371293
 =delta_5+10061627/2573060490
 >delta_5.                                           (26)
```

The script asserts this failure rather than promoting it. It also treats
the three one-equality partitions separately: the profile is labelled and
not invariant under every permutation of the cores. No joint Markov law,
independence assumption, false cap on an equal pair, or stationary joint
coupling is used.

What is proved is the partition reduction (22). The all-equal normalized
core remains an open branch of `(3,4,5)`; the 165 first-depth-one profiles
also remain. Hence this theorem does not prove LRC(14). QED.
