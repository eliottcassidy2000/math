---
id: THM-3392
title: "Bipartite sign lift, optimal factor-two synchronization, and the ternary support sidecar"
status: >
  PROVED + VERIFIED-EXACT + HOSTILE-AUDITED. For every real symmetric
  zero-diagonal matrix A, the bipartite sign norm
  B(A)=max_(x,y in {+-1}^n)|x^T A y| and synchronized sign norm
  Q(A)=max_(z in {+-1}^n)|z^T A z| satisfy Q(A)<=B(A)<=2Q(A).
  The factor two is asymptotically optimal even when every off-diagonal
  entry is a sign, by an explicit Sylvester-Hadamard two-block family.
  A 6-by-6 sign matrix has Q=10 and B=12; an exhaustive switching-class
  audit proves that order six is the first strict sign-matrix loss.
  Zero diagonal and symmetry are essential. This transfers a real
  Grothendieck bilinear relaxation to symmetric switching with an optimal
  universal synchronization tax, but does not solve the plus-minus
  extremum or make a symmetric sign matrix into a tournament.
source: codex-2026-08-14-grothendieck-transfer
related:
  - THM-2501-switching-fourth-moment-signed-c4-and-gram-energy
  - THM-3315-tournament-cut-switching-centered-coronal-walk-compiler
  - THM-3322-tournament-switching-second-moment-deletion-gram-and-order-join-law
script: 04-computation/switching_bipartite_synchronization_thm3392.py
output: 05-knowledge/results/switching_bipartite_synchronization_thm3392.out
script_sha256: fad03eca02ec8da0a4bd51ba1ff54e3830636b6c49b4fa563299e25d5865fae9
output_sha256: d2ef8e250c05eeaae2930c72737f67de6315faeb0406a081302286de9913f2f3
hash_basis: working-tree bytes (LF)
---

# THM-3392 -- bipartite sign lift and optimal synchronization loss

**PROVED + VERIFIED-EXACT + HOSTILE-AUDITED.**

Let `A` be a real symmetric `n x n` matrix with zero diagonal. Define

```text
Q(A) = max_(z in {+-1}^n) |z^T A z|,

B(A) = max_(x,y in {+-1}^n) |x^T A y|.                       (1)
```

Then

```text
Q(A) <= B(A) <= 2 Q(A).                                      (2)
```

Both the hypotheses and the constant have real content:

- the factor `2` is asymptotically optimal even among symmetric
  zero-diagonal matrices with every off-diagonal entry in `{+-1}`;
- there is already strict loss `B/Q=6/5` at order six; and
- if either symmetry or the zero-diagonal condition is removed, one can
  have `Q=0<B`, so no finite comparison remains.

The exact coordinate lost by the bipartite lift is synchronization of the
two sign vectors. Its minimal repair is a pair of disjointly supported
ternary sign vectors.

## 1. The ternary polarization identity

Fix `x,y in {+-1}^n` and put

```text
u=(x+y)/2,                 v=(x-y)/2.                          (3)
```

Thus `u,v in {-1,0,1}^n`, their supports are disjoint, and
`x=u+v`, `y=u-v`. Symmetry of `A` cancels the mixed terms:

```text
x^T A y
 =(u+v)^T A(u-v)
 =u^T A u-v^T A v.                                            (4)
```

Every ternary quadratic value is controlled by `Q(A)`. Indeed, given
`w in {-1,0,1}^n`, extend its nonzero coordinates to a random
`Z in {+-1}^n` by choosing independent unbiased signs on the zero
coordinates. Since `A` has zero diagonal, every term involving a newly
random coordinate has expectation zero, and hence

```text
E[Z^T A Z]=w^T A w.                                            (5)
```

Therefore

```text
|w^T A w| <= E|Z^T A Z| <= Q(A).                              (6)
```

Applying (6) to `u` and `v` in (4) gives

```text
|x^T A y| <= |u^T A u|+|v^T A v| <= 2Q(A).
```

Taking the maximum proves the upper bound in (2). The lower bound follows
by taking `x=y=z` in the bilinear maximum.

This proof identifies the sidecar exactly:

```text
bipartite pair (x,y)
   <-> disjoint ternary packet (u,v)
   -> two synchronized quadratic values u^T A u and v^T A v.  (7)
```

Keeping only the scalar bilinear optimum forgets which two disjoint
supports carry the competing quadratic signs. That support partition is
the synchronization data which a transfer must restore.

## 2. The factor two is sharp even for sign matrices

Let `m=2^r`, and let `H_m` be the Sylvester Hadamard matrix, so

```text
H_m H_m^T=mI_m.                                                (8)
```

On two blocks `L,R` of size `m`, form the symmetric sign matrix

```text
        [ -(J-I)    H_m ]
A_m  = [                  ].                                  (9)
        [  H_m^T    J-I  ]
```

It has zero diagonal and all other entries are signs. Take `x` to be all
ones and take `y` to be `-1` on `L` and `+1` on `R`. The two cross-block
contributions cancel, while each internal block contributes `m(m-1)`:

```text
B(A_m) >= x^T A_m y=2m(m-1).                                  (10)
```

For any sign vectors `s,t` on `L,R`, respectively,

```text
(s,t)^T A_m(s,t)
  =(sum t_i)^2-(sum s_i)^2+2s^T H_m t.                         (11)
```

The first difference has absolute value at most `m^2`. From (8) and
Cauchy--Schwarz,

```text
|s^T H_m t| <= ||s||_2 ||H_m t||_2=m^(3/2).                   (12)
```

Consequently

```text
Q(A_m) <= m^2+2m^(3/2),

B(A_m)/Q(A_m)
 >= 2m(m-1)/(m^2+2m^(3/2)) -> 2.                              (13)
```

Thus no universal constant smaller than `2` can replace the upper constant
in (2), even on the exact signed-complete-graph class of THM-2501.

## 3. The first strict finite sign witness

The following first-row-positive sign matrix has

```text
Q(A)=10,                 B(A)=12:

     0  1  1  1  1  1
     1  0 -1 -1  1  1
     1 -1  0  1 -1  1
A =  1 -1  1  0  1 -1 .                                      (14)
     1  1 -1  1  0 -1
     1  1  1 -1 -1  0
```

One bilinear maximizer is

```text
x=(-1,-1,-1,-1,-1, 1),
y=(-1, 1, 1,-1,-1,-1).
```

Its ternary packet is

```text
u=(-1, 0, 0,-1,-1, 0),       u^T A u= 6,
v=( 0,-1,-1, 0, 0, 1),       v^T A v=-6,                     (15)
```

so (4) gives `x^T A y=12`. Exhaustion gives `Q=10`.

The companion checks every switching class of sign matrices through order
six. Switching by `D=diag(d_i)` preserves both norms, and every class has a
unique representative with first row positive after fixing `d_1=1`.
Therefore the exact universe at order `n` has
`2^binom(n-1,2)` representatives. The ratio histograms are

```text
n=2:  {1:1}
n=3:  {1:2}
n=4:  {1:8}
n=5:  {1:64}
n=6:  {1:1012, 6/5:12}.                                      (16)
```

Hence `B=Q` through order five and order six is the first strict sign-matrix
loss. This minimal-order statement is **FINITE-EXACT**; no uncomputed
classification is inferred from it.

## 4. Equality and hypothesis boundaries

If a simultaneous sign switch makes all off-diagonal entries of `A`
nonnegative (or all nonpositive), then the all-one sign vector attains the
sum of all absolute entries. Thus `B(A)=Q(A)` on this balanced switching
stratum. A complete equality classification is not claimed.

The two structural hypotheses in the theorem are indispensable:

- without zero diagonal, `A=diag(1,-1)` has `z^T A z=0` for every sign
  vector, while `B(A)=2`;
- without symmetry, the skew matrix `[[0,1],[-1,0]]` likewise has
  `z^T A z=0` identically and `B(A)=2`.

The second hostile is the same representation boundary already emphasized
by THM-2501: a skew tournament arc matrix is not the symmetric switching
object.

## 5. Exact Grothendieck transfer and its stopping point

For a real matrix define the bipartite vector relaxation

```text
Gamma(A)=sup_(unit p_i,q_j) |sum_(i,j) a_ij <p_i,q_j>|.        (17)
```

The real Grothendieck inequality gives the cited implication

```text
Gamma(A) <= K_G B(A) <= 2K_G Q(A).                             (18)
```

The synchronized vector relaxation obtained by requiring `p_i=q_i` is a
subproblem of `Gamma(A)`, so it obeys the same bound. Combining (18) with
the **CITED preprint** upper bound

```text
K_G <= 1.7818666069360661
```

from Saha--Li--Xue et al., *New Lower and Upper Bounds for the
Grothendieck Constant*, [arXiv:2608.11158v2](https://arxiv.org/abs/2608.11158),
gives the numerical corollary

```text
Gamma(A) <= 3.5637332138721322 Q(A).                           (19)
```

The external preprint bound is not a dependency of (2), and (19) is only a
rounding/relaxation estimate. The typed transfer is:

```text
source:     bipartite Gaussian sign rounding;
target:     THM-2501 symmetric switching objective;
map:        use A on independent left and right sign copies;
preserved:  bilinear objective and its vector relaxation;
destroyed:  the constraint x=y;
sidecar:    the disjoint ternary support packet (u,v);
cost:       an optimal universal factor two;
next test:  seek a coupling-aware rounding on a restricted Gram/C4 class.
```

The Hadamard family proves that no argument which sees only `B(A)` and
`Q(A)` can improve that universal cost. Any sharper plus-minus theorem must
use additional structure, such as THM-2501's signed-four-cycle/Gram packet,
not merely a better value of `K_G`.

## 6. Exact companion

Run

```text
python3 04-computation/switching_bipartite_synchronization_thm3392.py
python3 -O 04-computation/switching_bipartite_synchronization_thm3392.py
```

Both commands reproduce

```text
05-knowledge/results/switching_bipartite_synchronization_thm3392.out
```

byte-for-byte. The dependency-free referee checks:

- all `1,099` switching-normalized sign matrices through order six;
- all `67,732` sign vectors used to compute both exact norms;
- the ratio histograms (16) and witness (14)--(15);
- the diagonal and skew hostiles; and
- the Sylvester Gram identity through order `64` for the sharp family.

The factor-two proof and the Hadamard asymptotic are symbolic; the finite
companion is an audit and the minimal-order certificate, not their logical
source. QED.
