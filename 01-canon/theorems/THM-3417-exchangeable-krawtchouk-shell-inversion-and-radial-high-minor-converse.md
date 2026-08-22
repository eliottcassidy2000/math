---
id: THM-3417
title: "Exchangeable Krawtchouk shell inversion and radial high-minor converse"
status: PROVED + VERIFIED-EXACT + INDEPENDENTLY PROOF/REPLAY-AUDITED
source: hadamard-krawtchouk-shell-2026-08-15
audit: independent Krawtchouk-normalization, radial-inverse, exchangeable-detector, convolution, Hermite, hostile, replay, hash, dependency, routing, and scope audit clean
depends_on:
  - THM-3413-strength-k-orthogonal-array-toggle-filtration-and-high-minor-converse
related:
  - THM-3407-hadamard-core-multitoggle-response-plaquette-shells-and-trade-distance
  - THM-3411-pairwise-independent-toggle-filter-and-sharp-high-minor-norm
  - THM-3396-four-bit-pairwise-independent-fourier-cone
verified_by:
  - 04-computation/hadamard_exchangeable_krawtchouk_shell_thm3417.py
  - 05-knowledge/results/hadamard_exchangeable_krawtchouk_shell_thm3417.out
script: 04-computation/hadamard_exchangeable_krawtchouk_shell_thm3417.py
output: 05-knowledge/results/hadamard_exchangeable_krawtchouk_shell_thm3417.out
script_sha256: f5c2a4f4c3917d25bbcca078ccbb28348c8e512d7724a626d8c61099db90ba4c
output_sha256: d8753c02e75d4987cc74dc7cb28f0838d0d1715760a6da30f88d9342e0c4f2d5
semantic_sha256: d2471a54b365e186d605bdc5640434219013efe1cff787a18fbf66989c4e355b
hash_basis: LF-normalized bytes
---

# THM-3417 -- exchangeable Krawtchouk shell inversion and radial high-minor converse

## 1. Setting and radial coordinates

Let `H` be a normalized real Hadamard matrix of order `4m`, with sign core
`K` and binary maxdet core `B=(J-K)/2`.  Choose `t>=1` **distinct core
positions**

```text
e_a=(i_a,j_a),              delta_a=K_(i_a,j_a),       a in [t].    (1)
```

Rows or columns may repeat.  For `S subseteq [t]`, in inherited event order,
put

```text
M_S=(product_(a in S) delta_a)
    det K[(i_a)_(a in S),(j_a)_(a in S)],       M_empty=1.          (2)
```

For `x in {+-1}^t`, toggle event `a` when `x_a=-1`, and let

```text
f(x)=det B_(toggle(x))/det B.                                     (3)
```

THM-3413 proves that the Walsh coefficient at nonempty `T` is

```text
C_T=sum_(S superseteq T) (-1)^(|S|+|T|) M_S/(4m)^|S|.             (4)
```

Extend this notation to the empty set by the same formula.  Thus

```text
C_empty=E_u f,                                                     (5)
```

where `u` is uniform on the cube.  Define the radial Walsh and signed-minor
shells

```text
A_s=sum_(|T|=s) C_T,                 0<=s<=t,
P_j=sum_(|S|=j) M_S,                 0<=j<=t.                      (6)
```

In particular `P_0=1`.  These are sums over labelled subsets, not absolute
values or multisets.

Use the binary Krawtchouk normalization

```text
K_s(w;t)=sum_r (-1)^r binom(w,r) binom(t-w,s-r),                  (7)
```

where an out-of-range binomial coefficient is zero, and put

```text
kappa_s(w)=K_s(w;t)/binom(t,s).                                  (8)
```

## 2. Fixed-weight shell transform and both exact inverses

Let `nu_w` be uniform on the Hamming shell containing exactly `w` minus
signs, and write

```text
g_w=E_(nu_w) f,                         0<=w<=t.                   (9)
```

Then the complete fixed-shell response is

```text
g_w=sum_(s=0)^t kappa_s(w) A_s.                                  (10)
```

All `t+1` shell means invert exactly:

```text
A_s=2^(-t) sum_(w=0)^t binom(t,w) K_s(w;t) g_w.                  (11)
```

There is a second, triangular expression directly in signed-minor shells:

```text
g_w=sum_(j=0)^w (-1/(2m))^j
               [binom(w,j)/binom(t,j)] P_j,                     (12)

P_j=(2m)^j binom(t,j)
    sum_(w=0)^j (-1)^w binom(j,w) g_w.                           (13)
```

The radial Boolean--Walsh change of coordinates is itself exact:

```text
A_s=sum_(j=s)^t (-1)^(j+s) binom(j,s) P_j/(4m)^j,                (14)

P_j=(4m)^j sum_(s=j)^t binom(s,j) A_s.                           (15)
```

Thus fixed Hamming weights lose no information **inside the permutation-
invariant quotient**.  They recover every aggregate shell `A_s` and `P_j`,
but not the individual labelled coordinates being summed.

## 3. Proof of the shell transform

Fix `T` of size `s`.  A uniformly random `w`-set of toggled coordinates
meets `T` in `r` coordinates with probability

```text
binom(s,r) binom(t-s,w-r)/binom(t,w).
```

Consequently

```text
E_(nu_w) chi_T
 =sum_r (-1)^r binom(s,r)binom(t-s,w-r)/binom(t,w)
 =K_w(s;t)/binom(t,w)
 =K_s(w;t)/binom(t,s)
 =kappa_s(w).                                                     (16)
```

The middle equality follows by counting pairs of an `s`-set and a `w`-set
with a prescribed intersection.  Fourier pairing of `f` with `(16)` gives
`(10)`.

For completeness, the Krawtchouk orthogonality needed for inversion is

```text
sum_(w=0)^t binom(t,w) K_s(w;t)K_r(w;t)
  =2^t binom(t,s) 1_(s=r).                                      (17)
```

Indeed `K_s(w;t)=sum_(|T|=s)chi_T(x)` for any `x` of weight `w`.  Summing
the product over the whole cube leaves only equal characters.  Multiply
`(10)` by `binom(t,w)K_s(w;t)`, sum in `w`, and use `(17)` to obtain `(11)`.

THM-3407's Boolean event polynomial is

```text
F(z)=sum_(S subseteq [t])(-1/(2m))^|S| M_S product_(a in S)z_a.  (18)
```

On the `w`-shell, a fixed `j`-set lies inside the toggled set with probability

```text
binom(w,j)/binom(t,j).
```

Averaging `(18)` therefore gives `(12)`.  Ordinary binomial inversion gives
`(13)`.  Finally, summing `(4)` over all `T` of size `s` gives `(14)`, because
each source set of size `j` contains `binom(j,s)` such targets.  Summing
THM-3413's individual inverse over all source sets of size `j` gives `(15)`.
This proves every displayed transform.  QED.

## 4. Every exchangeable law and the radial high-minor converse

Let `nu` be exchangeable under all permutations of the `t` event labels, and
let

```text
p_w=nu(number of minus signs=w).
```

Every character moment then depends only on its degree:

```text
lambda_s=E_nu chi_T
        =sum_(w=0)^t p_w kappa_s(w),          |T|=s.              (19)
```

Hence

```text
E_nu f=sum_(s=0)^t lambda_s A_s.                                (20)
```

The shell laws are the extreme exchangeable laws, so `(10)` and `(20)` are
equivalent descriptions of the entire exchangeable observer class.

Fix `0<=k<=t`.  The following are equivalent:

1. `E_nu f=E_u f` for every exchangeable strength-`k` Rademacher law `nu`;
2. equality holds for the finite symmetrized parity bank defined in Section 5
   at every degree `s>k` and both signs;
3. `A_s=0` for every `s>k`;
4. `P_j=0` for every `j>k`.

The equivalence of 3 and 4 is the triangular pair `(14)--(15)`.  Conditions
1--3 follow from the exact finite detectors below.  This is the radial
high-minor converse.  It is strictly weaker than THM-3413's labelled converse:
same-degree signed minors may cancel in `P_j`, and same-degree responses may
cancel in `A_s`.

## 5. Symmetrized parity halfcubes are exact finite detectors

For `1<=s<=t` and `sigma in {+-1}`, average THM-3413's parity-halfcube laws
over all characters of degree `s`:

```text
bar(nu)_(s,sigma)
 =binom(t,s)^(-1) sum_(|T|=s) nu_(T,sigma).                       (21)
```

This is an exchangeable probability law.  At a point of Hamming weight `w`,
its density relative to the uniform cube is

```text
d bar(nu)_(s,sigma)/du
 =1+sigma K_s(w;t)/binom(t,s).                                  (22)
```

It is nonnegative because `(21)` is an average of probability laws.  Its only
nonconstant Walsh moments are

```text
E_bar(nu) chi_T=sigma/binom(t,s)       when |T|=s.                (23)
```

Thus it has strength `s-1` and

```text
E_bar(nu) f-E_u f=sigma A_s/binom(t,s).                          (24)
```

Equation `(24)` proves the missing detector direction in Section 4.  It also
has a literal finite orthogonal-array compiler: concatenate the
`binom(t,s)` parity halfcubes.  This gives `binom(t,s)2^(t-1)` rows, with row
multiplicity at a point of weight `w` equal to

```text
[binom(t,s)+sigma K_s(w;t)]/2.                                  (25)
```

No Hadamard completion is asserted.

## 6. Fixed-shell convolution approaches the total-parity boundary

Let `nu_w^(star r)` be the coordinatewise product of `r` independent samples
from `nu_w`.  Characters diagonalize this convolution, so

```text
E_(nu_w^(star r))f-E_u f
 =sum_(s=1)^t A_s kappa_s(w)^r.                                 (26)
```

Suppose `0<w<t`.  Then

```text
kappa_t(w)=(-1)^w,
|kappa_s(w)|<1                    for 1<=s<t.                    (27)
```

For the strict inequality, choose `a in T` and `b notin T`.  There is a
`w`-set containing `a` but not `b`; swapping them changes the sign of
`chi_T`.  Hence that character is not constant on the shell, so its average
cannot have magnitude one.

Put `rho_w=max_(1<=s<t)|kappa_s(w)|<1`.  Then

```text
|E_(nu_w^(star r))f-E_u f-(-1)^(wr)A_t|
 <=rho_w^r sum_(s=1)^(t-1)|A_s|.                               (28)
```

More strongly, Fourier inversion shows that `nu_w^(star r)` approaches the
parity halfcube

```text
chi_[t]=(-1)^(wr).                                               (29)
```

For odd `w` the target alternates with `r`; for even `w` it is fixed.  The
endpoint shells `w=0,t` are point masses and are the sharp exceptions: many
characters then have moment of magnitude one.  This identifies a concrete
exchangeable route to THM-3413's fixed/period-two parity boundary.

## 7. Central Krawtchouk shells discretize Hermite chaos

For fixed `s`, take integers `w_t` such that

```text
x_t=(t-2w_t)/sqrt(t) -> x.                                      (30)
```

Then, with `He_s` denoting the probabilists' Hermite polynomial,

```text
t^(s/2) kappa_s(w_t) -> He_s(x).                                (31)
```

The exact first three formulas, with `y=t-2w`, are

```text
kappa_1=y/t,
kappa_2=(y^2-t)/[t(t-1)],
kappa_3=[y^3-(3t-2)y]/[t(t-1)(t-2)].                            (32)
```

To prove `(31)`, use the generating function

```text
sum_s K_s(w;t) z^s=(1+z)^(t-w)(1-z)^w.                          (33)
```

Substitute `z=u/sqrt(t)` and `(30)`.  The logarithm of the right side tends
coefficientwise to

```text
xu-u^2/2,
```

whose exponential is `sum_s He_s(x)u^s/s!`.  Since
`binom(t,s)~t^s/s!`, equation `(31)` follows.

This is a precise but limited bridge to Saha--Li--Xue et al.,
[*New Lower and Upper Bounds for the Grothendieck Constant*](https://arxiv.org/abs/2608.11158v2),
a **CITED VERY RECENT PREPRINT** whose rounding analysis uses Gaussian Hermite
projections and an explicit cubic Hermite coordinate.  The map here preserves
only chaos degree and the aggregate coefficient `A_s`.  It destroys labelled
minor coordinates, Gaussian sign geometry, cross-function inner products,
and inverse-majorant data.  Therefore `(31)` supplies no Grothendieck bound;
the preprint is context, not a dependency.

## 8. Sharp information-loss hostiles

### 8.1 Exchangeable strength two is strictly weaker than labelled strength two

In the Paley order-eight core, take the zero-based events

```text
((0,0),(1,2),(2,1),(3,4)).                                      (34)
```

Direct determinants give

```text
(P_0,P_1,P_2,P_3,P_4)=(1,4,4,0,0),
(A_0,A_1,A_2,A_3,A_4)=(9/16,3/8,1/16,0,0),
(g_0,g_1,g_2,g_3,g_4)=(1,3/4,13/24,3/8,1/4).                   (35)
```

Nevertheless the labelled high packet, listed at masks
`7,11,13,14,15`, is

```text
(1/128,0,-1/128,0,0).                                          (36)
```

Every exchangeable strength-two law is blind by `(35)`, while individual
parity halfcubes see the opposite nonzero responses in `(36)`.  The missing
coordinate is within-level localization, not high-degree mass.

### 8.2 Equal shell means need not even determine within-level multisets

In the Paley order-four core, the two physical packets

```text
((0,0),(0,1),(0,2),(2,0)),
((0,0),(0,1),(1,0),(1,1))                                      (37)
```

have the same shell means

```text
(1,1/2,1/6,0,0).                                                (38)
```

Their sorted degree-one Walsh profiles are respectively

```text
(0,1/8,1/8,1/4),
(1/8,1/8,1/8,1/8).                                              (39)
```

Because the multisets differ, no event relabelling repairs the loss.  Thus the
kernel is larger than a cosmetic choice of coordinate order.

## 9. Repeated rows and columns survive; duplicate positions do not

If selected events repeat a row or column, every affected signed minor in
`(2)` vanishes and all formulas above remain valid.  The order-four exhaustive
control includes `478` such nonmatching position sets among its `511`
nonempty sets.

Exact duplicate **positions** are excluded by `(1)`.  Two nominal toggles at
one entry act by XOR: activating both returns the original entry.  The
multilinear additive event compiler instead adds two copies of the same
rank-one update.  At the order-four position `(0,0)`, the naive additive
two-event values are

```text
(1,1/2,1/2,0),
```

whereas actual XOR toggling gives

```text
(1,1/2,1/2,1).                                                   (40)
```

Repeated positions must first be consolidated by event parity; they cannot be
silently admitted to this theorem.

## 10. Connection and loss ledger

| field | content |
|---|---|
| source | the full labelled THM-3413 toggle response `(C_T)` or signed-minor packet `(M_S)` |
| target | determinant means under all exchangeable laws on the event cube |
| map | average under `S_t`, then use the Krawtchouk transform `(10)` |
| preserved | every exchangeable expectation, radial strength filtration, aggregate signed-minor shells, total parity, exact finite inversion |
| destroyed | event localization, anisotropic same-degree contrasts, row/column addresses, Hadamard equivalence, Gaussian sign geometry |
| restoration sidecar | the labelled parity-halfcube bank, equivalently the individual `C_T` |
| cheapest decisive tests | H4 fibre collision `(37)--(39)`, H8 strict-loss packet `(34)--(36)`, repeated-row sets, duplicate-position XOR hostile `(40)` |

The quotient has only `t+1` coordinates rather than `2^t`; above strength
`k` it retains `t-k` aggregate coordinates rather than
`sum_(s>k)binom(t,s)` labelled ones.  This theorem proves no Hadamard
existence or completion result, determinant improvement, Grothendieck or
Crouzeix bound, LRC(14), or JC(2).

## 11. Exact companion and status

Run

```text
python3 04-computation/hadamard_exchangeable_krawtchouk_shell_thm3417.py
python3 -O 04-computation/hadamard_exchangeable_krawtchouk_shell_thm3417.py
```

The standard-library companion checks `819` Krawtchouk orthogonality
identities through `t=12`, both radial inversions on abstract packets through
`t=9`, and `90` symmetrized parity laws with `16,298` nonconstant moment
checks.  It directly checks all `511` nonempty distinct-position sets in the
Paley order-four core, including `478` repeated-row or repeated-column sets,
the H4 fibre collision, the H8 strict-loss hostile, an H8 top-minor
convolution through twelve powers at weights one and two, and an H12
seven-toggle control.  It also checks `506` internal-shell strict mixing
inequalities through `t=12`, the duplicate-position failure, all three exact
formulas in `(32)` for `2,565` inputs, and central Hermite scaling controls.

Normal and optimized runs are byte-identical to the frozen output.  The
LF-normalized script/output hashes are respectively

```text
f5c2a4f4c3917d25bbcca078ccbb28348c8e512d7724a626d8c61099db90ba4c
d8753c02e75d4987cc74dc7cb28f0838d0d1715760a6da30f88d9342e0c4f2d5
```

and the semantic digest is
`d2471a54b365e186d605bdc5640434219013efe1cff787a18fbf66989c4e355b`.
An independent immutable-file proof/replay audit rederived every transform,
checked the detector, convolution, Hermite, and hostile boundaries, and
reproduced the frozen artifact under normal and optimized Python.
