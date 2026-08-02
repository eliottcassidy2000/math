---
id: THM-3093
title: "Arbitrary-gap remote-cluster Monge flag compactification"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT REQUESTED.  Above any
  fixed nonzero physical child resultant, an arbitrary fixed-rank remote
  cluster is eventually positive uniformly over every distinct integer gap
  vector, with no bound on its diameter.  Strict factorial-response Monge
  inequalities expose a staircase flag; multivariate Jensen contracts every
  physical error atom, while divergent adjacent gaps split the normal system
  into positive generalized-alternant blocks.  The threshold in the first
  remote offset is independent of all internal gaps.  This is fixed child
  and fixed width, not arbitrary-support GMC(2).
source: root-gmc-arbitrary-cluster-flag-2026-08-01
depends_on:
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
  - THM-3073-upper-triangular-resultant-norm-and-torsion-character-death-barcode
  - THM-3085-multi-normal-fixed-gap-cluster-and-unconditional-all-width-tail
  - THM-3089-logarithmic-moving-gap-cluster-cone-and-condition-number-boundary
related:
  - THM-3086-arbitrary-cluster-composition-chambers-and-alternant-clutch-holotopy
  - THM-3091-arbitrary-gap-remote-pair-desuspension-and-exact-Jensen-contraction
script: 04-computation/gmc_arbitrary_gap_remote_cluster_flag_thm3093.py
output: 05-knowledge/results/gmc_arbitrary_gap_remote_cluster_flag_thm3093.out
script_sha256: 808b69c71533d576f41bbc38f12cec70907665e5b2a90371f5572c335c02750d
output_sha256: d120b193190a83ea3af4b0f2b1e15a845aa9338cfffae9a1d4c6c26ad8806b81
hash_basis: LF-normalized bytes
---

# THM-3093 -- arbitrary-gap remote-cluster Monge flag compactification

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT REQUESTED.**

THM-3085 permits any fixed internal gaps.  THM-3089 moves them only through
a square-root cone because it compares the whole normal system with one line
power.  THM-3091 removes every restriction for a pair by a second triangular
normalization.  The correct higher-rank object is not one line: it is a
staircase flag exposed by the strictly Monge factorial-response matrix.

## 1. Uniform arbitrary-gap theorem

Fix a physical lower support of width `m>=1` with normalized first-window
resultant

```text
S_m!=0.                                                     (1)
```

For `m=1`, this means the empty resultant `S_1=1`.  For `m=2`, it is the
positive one-variable quadratic coefficient

```text
S_2=L((f_(n+a)-f_(n+b))^2)/L(f_n^2)>0.
```

Fix `q>=2`, and put

```text
p=m+1,                   k=m+q,
r_i=p+i-1                (1<=i<=q),
D_z=product_i r_i=k!/m!,
rho_m=(m/(m+1))^m<1.                                      (2)
```

For every positive integer `C`, choose an arbitrary distinct gap vector

```text
0=h_1<h_2<...<h_q,                                        (3)
```

with no upper bound and no prescribed dependence on `C`, and append the
physical offsets

```text
C+h_1,...,C+h_q.                                          (4)
```

There is a threshold `C_0`, depending on the fixed child, `m`, and `q`, but
**independent of the entire gap vector**, such that for every `C>=C_0` the
width-`k` physical first-window resultant satisfies

```text
R_k(C;h_1,...,h_q)>0.                                     (5)
```

More precisely, Section 3 defines positive carriers and a normalized normal
resultant `Etilde_(C,h)`.  There is `eta_(m,q)>0` such that, uniformly in
`(3)`,

```text
Etilde_(C,h)>=eta_(m,q),                                  (6)
```

and

```text
R_k(C;h)=S_m^D_z Etilde_(C,h)^m!
 product_(i=1)^q U_(r_i,C+h_i)^[k!/r_i]
 [1+O_(child,m,q)(poly(C)rho_m^C)].                      (7)
```

The product `D_z` contains two consecutive integers and is even.  Every `U`
and, by Section 5, `Etilde` are positive, so `(7)` proves `(5)` even if the
child resultant has negative sign.  This also handles `m=1`, where `m!=1`.

The cases `m=1,2` above are unconditional, and for `m=3`, THM-2824 supplies
`S_3>0` for every arbitrary physical three-slot base.  In particular, at
fixed minimum exponent and fixed width `k>=3`, every support whose first gap
is sufficiently large is good, regardless of every later gap.  At widths
above three the same statement also holds over any fixed two- or three-slot
child.  The threshold is uniform over the internal tail geometry.

## 2. Exact response matrix and its Monge flag

Set `N=n+C` and define, for an integer `h>=0`,

```text
V_r(h)=L(f_(N+h)^r)/L(f_N^r)
      =(rN+1)_(rh)/(N+1)_h^r.                            (8)
```

For an adjacent gap interval put

```text
phi_t(r)=(1/r)log[V_r(h_(t+1))/V_r(h_t)].                (9)
```

Extending the endpoint temporarily to a real variable gives the exact
digamma integral

```text
phi_t(r)=integral_(h_t)^(h_(t+1))
 [psi(r(N+x)+1)-psi(N+x+1)] dx.                         (10)
```

Its derivative in `r` is the integral of
`(N+x)psi'(r(N+x)+1)>0`.  Therefore every `phi_t(r)` is strictly increasing
in `r`.  Equivalently,

```text
[V_a(h_(t+1))/V_a(h_t)]^b
 <[V_b(h_(t+1))/V_b(h_t)]^a       (0<a<b),              (11)
```

the exact rational strict-Monge inequality after roots are cleared.

Choose column and row potentials by

```text
lambda_1=1,
lambda_(t+1)/lambda_t
 =[V_(r_t)(h_t)/V_(r_t)(h_(t+1))]^(1/r_t),

a_i=V_(r_i)(h_i)lambda_i^r_i,
b_i=a_i^(1/r_i).                                         (12)
```

The pure normalized root weight in row `i`, column `t`, is

```text
W_(i,t)=V_(r_i)(h_t)^(1/r_i)lambda_t/b_i.                (13)
```

By `(11)`, as `t` increases the row weights rise strictly up to `i`, tie at
`i,i+1`, and then fall strictly; the last row has its unique maximum at `q`:

```text
0<W_(i,t)<=1,
W_(i,i)=W_(i,i+1)=1       (i<q),
W_(q,q)=1.                                                (14)
```

Indeed, if `g_t(r)=exp(phi_t(r))`, then

```text
lambda_i=product_(t<i)g_t(r_t)^(-1),
b_i=product_(t<i)g_t(r_i)/g_t(r_t)>=1.                  (15)
```

Thus `r_i^(-1)log V_(r_i)(h_t)` is strictly Monge, while the logarithms of
`lambda_t,b_i` are assignment-dual potentials exposing the staircase graph

```text
{(i,i),(i,i+1):i<q} union {(q,q)}.                       (16)
```

This is the finite symbolic state of the proof.

## 3. Normal forms and exact carrier covariance

In direct high coordinates `v_1,...,v_q`, the normalized all-high form in
degree `r` is

```text
A_r(v)=sum_(|alpha|=r) multinom(r;alpha) Q_(r,alpha)v^alpha,

Q_(r,alpha)
 =(rN+1)_(sum_t alpha_t h_t)
   /product_t (N+1)_(h_t)^alpha_t.                       (17)
```

Its pure coefficient at column `t` is exactly `V_r(h_t)`.  Define

```text
B_i(v)=a_i^(-1)A_(r_i)(lambda_1v_1,...,lambda_qv_q),
Etilde_(C,h)=Res(B_1,...,B_q).                           (18)
```

Apply `(18)` to the entire THM-3085 outer-transformed physical system, not
only to its normal face.  Scaling the variables contributes
`product_t lambda_t^k!`; scaling equation `i` contributes
`a_i^[-k!/r_i]`.  Since `a_i=V_(r_i)(h_i)lambda_i^r_i`, every lambda power
cancels:

```text
product_t lambda_t^k! product_i a_i^[-k!/r_i]
 =product_i V_(r_i)(h_i)^[-k!/r_i].                     (19)
```

Combining this with the common-base outer covariance gives

```text
Rhat_k(C,h)=R_k(C,h)/
 product_i [U_(r_i,C)V_(r_i)(h_i)]^[k!/r_i]

            =R_k(C,h)/product_i U_(r_i,C+h_i)^[k!/r_i]. (20)
```

The direct high coordinates and THM-3085 fixed-pivot coordinates are related
by THM-3089's unimodular map.  Its determinant sign is raised to even `D_z`,
so no hidden orientation sign enters `(20)`.  The absolute carrier forgotten
by the Monge dual is restored exactly; no root choice remains.

## 4. Multivariate Jensen controls every physical atom

Let `F_r(s)=log(rN+1)_s`.  Its discrete increments increase, so it is convex
and `F_r(0)=0`.  For a multi-index `alpha` with `j=sum alpha_t<=r`, add
`r-j` copies of the zero gap and apply Jensen at the barycentre.  This gives

```text
Q_(r,alpha)^r<=product_t V_r(h_t)^alpha_t.               (21)
```

Take an unsigned physical expansion atom in row `r`.  In the decomposition

```text
Y=Y_low+sum_t v_t(f_(N+h_t)-f_c),
```

let `alpha_t` be its actual remote-factor counts, `beta_t` its direct normal
exponents, and `A` the total offset of its lower factors.  Each actual remote
factor consumes its displayed coordinate; fixed-pivot subtraction may add
normal degree but cannot remove that occurrence.  Thus

```text
beta_t>=alpha_t,
j=sum alpha_t,
D=(r-j)n+A,
S=sum alpha_t h_t.                                       (22)
```

Relative to collapsing all internal gaps to zero, its exact multiplier is

```text
Q_atom=(jN+D+1)_S/product_t(N+1)_(h_t)^alpha_t.          (23)
```

After one child-dependent threshold, `jN+D<=rN`; hence `(21)` also bounds
`Q_atom`.

The potentials decrease, `0<lambda_t<=1`.  For a lower row `r<=m`,

```text
P_(r,t)=V_r(h_t)^(1/r)lambda_t
       =product_(s<t)g_s(r)/g_s(r_s)<=1,                (24)
```

because `r<r_s`.  Therefore

```text
Q_atom product_t lambda_t^beta_t
 <=product_t P_(r,t)^alpha_t<=1.                        (25)
```

For upper row `r_i`, equations `(13)--(15)` give

```text
a_i^(-1)Q_atom product_t lambda_t^beta_t
 <=b_i^[j-r_i] product_t W_(i,t)^alpha_t<=1.            (26)
```

These estimates are unsigned and precede all inclusion cancellation.  The
flag normalization never enlarges a nonsurviving physical layer.  For
`m=2`, the same atom inequalities apply with its one-variable lower block
and the strictly positive factorial variance `S_2`.  For `m=1` there are no
lower variables or lower forms: set `H_p=0`, use the convention
`Res_empty=S_1=1`, and apply THM-3073 only to the `q`-dimensional normal
block.  With these explicit low-child conventions, the base-one bank still
has only

```text
H_2,...,H_m,H_p+B_1,B_2,...,B_q                         (27)
```

at base one; every other coefficient remains

```text
O_(child,m,q)(poly(C)rho_m^C)                            (28)
```

uniformly over gap vectors of arbitrary diameter.  This is the magnitude
sidecar of the proof.

## 5. Composition-cube compactification and uniform nonvanishing

All coefficients of every `B_i` are bounded by their multinomial
coefficients, by `(21)` and `(26)`.  Suppose uniform nonvanishing failed.
Choose a countersequence `C->infinity`; after passing to a subsequence, every
adjacent integer difference

```text
Delta_t=h_(t+1)-h_t                                     (29)
```

is either fixed or tends to infinity, and every bounded coefficient
converges.  The fixed edges partition `{1,...,q}` into consecutive blocks.

For fixed degrees `a<b`, blocking the rising factorial as in THM-3091 gives,
on an interval of length `Delta`,

```text
g_Delta(a)/g_Delta(b)<=c_(a,b)^Delta,
c_(a,b)<1.                                               (30)
```

Hence a divergent edge kills every backward pure weight exponentially;
mixed coefficients containing a killed column vanish by `(21)`.  Equations
strictly to the right of a block lose all earlier variables.  Only the last
row of an earlier block may retain bounded forward variables.  The limit is
block upper triangular in the exact sense of THM-3073.

For a block `I=[a,b]`, put `d_j=h_j-h_a` and retain its fixed internal gaps.
Its diagonal restriction consists of powers of linear forms.  Up to positive
column multipliers `Lambda_j`, the line matrix is

```text
L_(i,j)=r_i^d_j Lambda_j/(r_i^d_i Lambda_i),      i,j in I.
```

The multipliers cancel in the determinant:

```text
det L_I=det(r_i^d_j)_(i,j in I)/product_(i in I)r_i^d_i

 >=Vandermonde(r_a,...,r_b)
    /product_(i=a)^b r_i^(i-a)>0.                        (31)
```

This is the generalized-Vandermonde/Schur monomial floor from THM-3089 and
is independent of the fixed gaps' values.

Iterating THM-3073 removes every arbitrary forward coefficient.  If `P` is
the consecutive-block partition, the limiting normal resultant is

```text
lim Etilde_(C,h)=product_(I in P)(det L_I)^D_z>0.         (32)
```

There are only `2^(q-1)` partitions, and `(31)` gives a positive row-only
floor on each.  Compactness contradicts the countersequence, proving `(6)`.

Finally, THM-3073 applied to `(27)` gives

```text
S_m^D_z Etilde_(C,h)^m!.                                 (33)
```

Equations `(19)--(20)`, `(28)`, and `(6)` prove `(7)`.

## 6. Holotopy and the composition extension

The exposed staircase `(16)` is a path.  Sending one adjacent gap to
infinity deletes an edge; sending several gaps to infinity produces a face
of the `(q-1)`-dimensional composition cube.  Formula `(32)` is its clutch
law.  The bounded cluster, every partially desuspended cluster, and the fully
iterated one-normal flag are one compactified physical holotopy with no sign
wall.

This physicalizes the growing-gap alternant clutch that THM-3086 retained
only symbolically.  At each fixed-width node of one of its already admissible
scale compositions, internal fixed gaps may be replaced by arbitrary
distinct integer gaps.  Equations `(25)--(26)` do not worsen that node's
child-scale entropy bill, while `(32)` supplies its positive normal sidecar.
The pre-existing moving-child entropy chamber remains necessary; this
theorem removes only the **internal cluster-diameter** restriction.

In the puzzle language of the project, `(16)` is the finite relation state,
`(31)--(32)` are the algebraic defect/carry carrier, and the strict
`rho_m^C` bound in `(28)` is the magnitude sidecar.  None can replace the
other two.

## 7. Boundaries and exact evidence

The child support, `m`, and `q` are fixed while `C->infinity`.  Repeated gaps
give repeated columns and kill a diagonal alternant.  A moving child, growing
width, `S_m=0`, unordered/nonconsecutive degree systems without a fresh
orientation/parity audit, and arbitrary equal-scale supports are outside the
theorem.  The `q=1` case is THM-3069 and does not inherit the multi-normal
even exponent automatically.

The result does not prove arbitrary-radial GMC(2), NC2, LRC(14), JC(2), or
DC(2).  Its unconditional low-child corollary is a remote-tail theorem, not
a claim that the first gap of every support is already beyond its
minimum- and width-dependent threshold.  Thus any counterexample at fixed
minimum and width is confined to a bounded first-gap tube, but that tube is
not emptied here.

The exact companion verifies:

1. `960` cleared rational response inequalities `(11)`;
2. `12,532` multivariate Jensen cells `(21)`;
3. `870` root-free flag weights and every staircase tie `(14)`;
4. `4,880` raw-start lower/upper physical atom contractions `(25)--(26)`;
5. `945` diagonal block determinants, Schur floors, and arbitrary-forward
   triangular controls over all composition faces through rank five;
6. all `60` displayed rank-three boundary faces;
7. `112` exact lambda-covariance cancellations and `1,344` carrier
   telescopes;
8. five explicit low-child controls: the empty resultant/product at `m=1`
   and four exact positive two-slot factorial variances at `m=2`, together
   with the repeated-gap and one-normal parity boundaries.

Run

```text
python 04-computation/gmc_arbitrary_gap_remote_cluster_flag_thm3093.py
python -O 04-computation/gmc_arbitrary_gap_remote_cluster_flag_thm3093.py
```

Both modes must equal the stored transcript after LF normalization.

**QED, pending independent audit and status promotion.**
