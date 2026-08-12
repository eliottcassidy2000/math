---
id: THM-3347
title: "U-spine signed prime-clock Gram kernel and projective-fold boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Signed Hensel-branch
  features make the U-spine log-content polarization an exact positive
  semidefinite Gram kernel.  At fixed admissible grade, conjugation is the
  antipode: the naive absolute Gram and folded content distance have a sharp
  rank-four failure, while the squared/projector and cosh kernels descend
  positively and retain the parent state.  The projective realization embeds
  the cubical quotient but needs its tautological-line sidecar to retain the
  conjugation cover.  No LRC, Berggren, tournament, or JC transfer follows.
audit: >
  Independent arithmetic and topology audits rederived the unrestricted
  Hensel-layer Gram law, both Hilbert distance invoices, normalized equality
  scope, fixed-grade Walsh spectra, the universal rank-three boundary, the
  all-rank Catalan level-four obstruction, and the literal N=32045 hostile.
  They separately checked the projector, finite even-tensor rank filtration,
  cosh repair, tautological unit-sphere cover, and quotient-loss typing.
  Normal, optimized, and stored transcripts byte-match; both recorded hashes
  match.
source: codex-2026-08-12-u-spine-prime-clock-gram
depends_on:
  - THM-3346-u-spine-prime-toggle-root-atlas-and-conjugation-monodromy
related:
  - THM-3336-primitive-gaussian-multiplication-content-curved-farey-triangulation
  - THM-2216-residual-capacity-hinge-gram-law
  - THM-2221-tournament-context-cut-metric-and-pinned-transport-response
  - THM-584-complement-is-antipodal-map-level-parity-spectrum
script: 04-computation/u_spine_signed_clock_gram_thm3347.py
output: 05-knowledge/results/u_spine_signed_clock_gram_thm3347.out
script_sha256: c0048e61640598f88e3992579215030c3087c1664d1d07986455546eabe8c01d
output_sha256: 9dc5417ae300384f837f983f80f0bf656fb5ea16ed1b06fa339c53d19495632e
hash_basis: working-tree bytes (LF)
---

# THM-3347 -- Signed prime clocks are Hilbertian before the antipodal fold

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The proof below is exact in its stated scopes.  No literature-priority claim
is made.

## 1. Inheritance and the connection contract

THM-3346 proves the U-spine content splitter.  For

```text
z_t=(t+1)+it,                   C_t=2t^2+2t+1,           (1)
d_-(r,s)=gcd(C_r,s-r),          d_+(r,s)=gcd(C_r,r+s+1), (2)
```

one has

```text
gcd(C_r,C_s)=d_-(r,s)d_+(r,s),
gcd(d_-(r,s),d_+(r,s))=1.                              (3)
```

Here `z_t` is the consecutive Euclid spinor and

```text
Phi(z_t)=(2t+1,2t(t+1),C_t)                            (3a)
```

is a primitive Pythagorean triple with inradius `t`.  Thus the prime clocks
below live simultaneously on a sum-of-two-squares norm and on the U-spine of
primitive right triangles.

The difference channel records equal local roots of `-1`; the reflected-sum
channel records opposite roots.  THM-3336 supplies the underlying Gaussian
content/charge law.  THM-2221 is the nearest proved signed-Gram/negative-type
mechanism, and THM-584 is the nearest antipodal Fourier mechanism.  Neither
identifies the U-spine kernel or the projective-fold boundary below.
The unsigned log-gcd layer-cake kernel is standard and is not claimed as new.

The connection is typed as follows.

| Item | Content |
|---|---|
| source | integral U-spine indices, then the fixed-grade CRT root cube |
| target | a prime-clock Hilbert space, then real projective space |
| map | Hensel branch signs, followed by rank-one projectors |
| preserved | shared valuations, both content channels, folded weight and parent state |
| first loss | conjugation forgets the signed Gaussian lift |
| second loss | a scalar vertex Gram forgets cell attachment and the tautological line |
| needed sidecar | raw sign section or tautological line bundle, plus any consumer owner/phase/height |
| cheapest hostile | `N=32045=5*13*17*29` |

## 2. The unrestricted signed log-content kernel is Gram

For every split prime `p=1 mod 4` and `e>=1`, independently choose one of the
two roots `rho_(p,e)` of `C_t=0 mod p^e`.  No compatibility between successive
Hensel layers is needed.  For nonnegative integral indices define

```text
sigma_(p,e)(t)=
  +1,  t=rho_(p,e) mod p^e,
  -1,  t=-1-rho_(p,e) mod p^e,
   0,  p^e does not divide C_t.                          (4)
```

Changing the choice at one layer flips that coordinate globally and therefore
does not change the Gram matrix.  Put

```text
Psi(t)=(sqrt(log p) sigma_(p,e)(t))_(p,e) in ell^2.      (5)
```

Every `C_t` is odd and all its prime divisors split in `Z[i]`, so (5) has
finite support.  Define the signed log-content polarization

```text
Lambda(r,s)=log d_-(r,s)-log d_+(r,s).                  (6)
```

Then

```text
Lambda(r,s)=<Psi(r),Psi(s)>,
Lambda(t,t)=||Psi(t)||^2=log C_t.                       (7)
```

Indeed a common `p^e` layer contributes `+log p` precisely when the two roots
lie on the same Hensel branch and `-log p` precisely when they lie on opposite
branches.  Summing over `e<=min(v_p(C_r),v_p(C_s))` gives (6).  Thus every
finite `Lambda` matrix is positive semidefinite.

Let

```text
g=gcd(C_r,C_s),           A=C_r/g,          B=C_s/g.     (8)
```

The two Hilbert distances have exact arithmetic invoices

```text
||Psi(r)-Psi(s)||^2=log(A B d_+(r,s)^4),
||Psi(r)+Psi(s)||^2=log(A B d_-(r,s)^4).                (9)
```

This is an unrestricted statement: it separates the coprime height tails
`A,B` from the shared same/opposite Gaussian contents.

For positive indices `r,s`, the normalized kernel

```text
kappa(r,s)=Lambda(r,s)/sqrt(log C_r log C_s)             (10)
```

is positive semidefinite and lies in `[-1,1]`.  Equality in absolute value
means collinearity of the clock-feature vectors.  For `r,s>=1`, this gives
`|kappa(r,s)|=1` if and only if `r=s` (and then `kappa=1`).  Indeed collinearity
forces identical nonzero layer support, hence identical valuations and
`C_r=C_s`; strict growth on the nonnegative chart gives `r=s`.  The value at
`t=0` is not normalized because `log C_0=0`.  Allowing negative indices would
reintroduce conjugate equality, so the chart restriction is load-bearing.

## 3. A fixed grade is a weighted sign cube

Let

```text
N=product_(j=1)^k q_j,
q_j=p_j^e_j,                         w_j=log q_j,
W=sum_j w_j=log N.                                      (11)
```

Here `N>1`, the `p_j` are distinct primes congruent to `1 mod 4`, and every
`e_j>=1`, exactly as in THM-3346.

For the modular root atlas `R_N` of THM-3346, choose CRT signs
`epsilon_j(t) in {+-1}` and define

```text
Psi_N(t)=(sqrt(w_j) epsilon_j(t))_(j=1)^k.               (12)
```

Writing `delta_-,delta_+` for the `N`-primary channels of THM-3346 gives

```text
Lambda_N(r,s)=<Psi_N(r),Psi_N(s)>
 =log delta_-(r,s)-log delta_+(r,s),                    (13)
||Psi_N(t)||^2=W.                                       (14)
```

On the raw `2^k` root cube, the Gram matrix has rank `k`; its nonzero Walsh
eigenvalues are `2^k w_j` on the coordinate characters.  Writing
`bar(t)=N-1-t`, conjugating one lift negates its feature vector and its pairing
with an unchanged lift, while simultaneous conjugation preserves the pairing:

```text
Psi_N(bar(t))=-Psi_N(t),
Lambda_N(bar(r),s)=-Lambda_N(r,s),
Lambda_N(bar(r),bar(s))=Lambda_N(r,s).                  (15)
```

Hence signed `Lambda_N` does not descend to the fixed-parent quotient.

Let `h_w(r,s)` be the weighted Hamming distance of the two sign words.  Then

```text
Lambda_N=W-2h_w,
d_N([r],[s])=min(h_w,W-h_w)=(W-|Lambda_N|)/2.            (16)
```

Thus the folded Gaussian content is an absolute-correlation invoice.  The
temptation to call `|Lambda_N|` a Gram kernel is false.

## 4. Sharp low-rank survival and the rank-four hostile

For every positive weight vector, the entrywise absolute kernel
`|Lambda_N|` is positive semidefinite when `k<=3`; equivalently the folded
metric `d_N` is conditionally negative definite there.

Only rank three needs proof.  Order the weights `a,b,c` with `c=max`.  On the
four parent classes, the four Walsh eigenvalues of the absolute kernel are:

```text
if c>=a+b:     (4c,4a,4b,0),
if c<a+b:      2(a+b+c, a-b+c, -a+b+c, a+b-c).          (17)
```

They are nonnegative.  Ranks one and two are immediate.  This universal
statement is sharp.  The least squarefree admissible rank-four grade is

```text
N=5*13*17*29=32045.                                     (18)
```

Choose the eight quotient representatives whose `5`-adic sign is positive,
with roots

```text
(1081,3546,4851,7316,5501,7966,9271,11736),             (19)
```

and coefficients

```text
c_x=(-1)^popcount(x)=(1,-1,-1,1,-1,1,1,-1).             (20)
```

Then `sum c_x=0` and exact prime-log collection gives

```text
c^T [d_N] c=16 log 5>0,
c^T [|Lambda_N|] c=-32 log 5<0.                         (21)
```

Therefore `d_N` is not of conditional negative type and `|Lambda_N|` is not
positive semidefinite.  The first failure is Fourier level four: equal
weights already give positive level-four coefficients in every rank `k>=4`.
For `k=2m` that unnormalized quotient coefficient is
`2 Catalan_(m-2)`; for `k=2m+1` it is `Catalan_(m-1)`.  Expanding the level-four
Krawtchouk sum and applying Pascal twice reduces it to these Catalan values.
The theorem does **not** claim failure for every rank-four
weight vector.  If one weight satisfies `w_j>=sum_(i!=j) w_i`, orient each
parent by that sign;
then (16) unfolds to ordinary weighted Hamming distance and is CND.

One tempting geometric reformulation fails with the same witness.  The naive
antipodal quotient-sphere chord satisfies

```text
q_sph([r],[s])^2=min(||u_r-u_s||^2,||u_r+u_s||^2)
                  =2-2|Lambda_N|/W=4d_N/W.             (21a)
```

Consequently its square is not CND in rank four.  This quotient-sphere chord
is not the projector chord used next.

## 5. Two exact positive repairs

### 5.1 Projector / Veronese repair

Normalize `u_t=Psi_N(t)/sqrt(W)` and map a parent to the rank-one projector

```text
P_[t]=u_t u_t^T.                                        (22)
```

This is independent of the sign lift.  Its kernel and squared chordal
distance are

```text
K_2([r],[s])=<P_[r],P_[s]>_F=(Lambda_N/W)^2,            (23)
delta_ch^2=1-K_2
 =4(d_N/W)(1-d_N/W),                                    (24)
||P_[r]-P_[s]||_F^2=2 delta_ch^2.                       (25)
```

Thus `K_2` is positive semidefinite and `delta_ch^2` is CND for every rank
and every positive weight vector.  The transform in (24) is strictly
increasing for `0<=d_N<=W/2`, so it loses no parent-distance information.

On the quotient cube, `Lambda_N^2` has Fourier support exactly at levels zero
and two, directly from

```text
Lambda_N(x,y)^2=sum_j w_j^2
 +2 sum_(i<j) w_i w_j chi_{ij}(x-y).                   (25a)
```

Its rank is

```text
1+binom(k,2),                                            (26)
```

More precisely, its unnormalized Walsh eigenvalue is

```text
2^(k-1) sum_j w_j^2                 at level zero,
2^k w_i w_j                          on the pair {i,j},
0                                    at every other even level. (26a)
```

Hence the centered Euclidean projector embedding has dimension `binom(k,2)`.
The pair characters are the positive replacement for the level-four and
higher curvature introduced by the absolute fold.

### 5.2 Strictly positive cosh repair

For every `tau>0`, define

```text
K_tau([r],[s])=cosh(tau Lambda_N(r,s)).                 (27)
```

It descends because cosh is even.  On an even quotient character
`A subset {1,...,k}`, the exact eigenvalue is

```text
lambda_A=2^(k-1)
 product_(j in A) sinh(tau w_j)
 product_(j notin A) cosh(tau w_j)>0.                  (28)
```

So `K_tau` is strictly positive definite on all `2^(k-1)` parents.  It is
lossless for the folded content:

```text
|Lambda_N|=tau^(-1) arcosh K_tau,
d_N=(W-|Lambda_N|)/2,                                  (29)
{delta_-,delta_+}={exp(d_N),exp(W-d_N)}.                (30)
```

Here `arcosh` is applied entrywise.  Formula (28) follows by product-expanding
`exp(tau sum_j w_j chi_j)`, averaging it with its antipodal sign reversal so
only even characters remain, and halving the raw-cube Fourier sum.  Every
displayed factor is positive.

At `tau=1`, all entries and eigenvalues are rational functions of the integer
prime powers `q_j`, since `cosh(log q)=(q+q^(-1))/2` and
`sinh(log q)=(q-q^(-1))/2`.

### 5.3 The finite even-tensor hierarchy

The projector and cosh repairs are the first and the generating-function ends
of an exact hierarchy.  For `m>=1`, put

```text
K_(2m)([x],[y])=Lambda_N(x,y)^(2m)
 =<Psi_N(x)^(tensor 2m),Psi_N(y)^(tensor 2m)>.          (30a)
```

It is PSD and descends to parents.  Index quotient characters by even subsets
`A subset {1,...,k}` and define

```text
P_(2m,A)(w)=sum_(alpha_1+...+alpha_k=2m,
                      alpha_j=1_A(j) mod 2)
 (2m)!/(product_j alpha_j!) product_j w_j^alpha_j.      (30b)
```

The exact unnormalized eigenvalue is

```text
lambda_A=2^(k-1) P_(2m,A)(w).                          (30c)
```

Because every `w_j>0`, this is positive exactly when `|A|<=2m` and is zero
otherwise.  Consequently

```text
rank K_(2m)=sum_(j=0)^min(m,floor(k/2)) binom(k,2j),    (30d)
```

and `K_(2m)` is strictly positive definite first at
`m=max(1,floor(k/2))` among `m>=1`.  The proof expands the `2m` tensor words:
their odd-occurrence set is precisely `A`; such a word exists exactly in the
stated range.  Thus degree `2m` sees exactly the even prime-toggle interactions
through size `2m`.  It is already lossless for the folded arithmetic since

```text
d_N=(W-K_(2m)([x],[y])^(1/(2m)))/2.                    (30e)
```

It still loses the oriented content channel.  Finally,

```text
K_tau=sum_(m>=0) tau^(2m) K_(2m)/(2m)!                 (30f)
```

explains why cosh activates every even character at once.

Do not confuse algebraic rank with point separation.  The even Veronese map
`[u]->u^(tensor 2m)` embeds `RP^(k-1)` for every `m>=1`, so already `K_2`
separates every parent even when its Gram matrix is not strictly PD (the first
case is rank `7<8` at `k=4`).  Higher `m` adds higher even interactions and
changes the extrinsic Veronese geometry, not the underlying projective
topology or cubical skeleton dimension.

## 6. The projective topology and the sidecar it still needs

Assume throughout this section that `k>=3`; ranks one and two have the fixed-
cell boundaries isolated in THM-3346.

Center the cube as `I=[-1,1]^k`, put `D=diag(sqrt(w_j))`, and define

```text
F:partial I -> S^(k-1),                F(y)=Dy/||Dy||.   (31)
```

Every ray through the origin meets the boundary of `D(I)` once, so `F` is an
odd homeomorphism.  It sends the sign vertices to `u_t` and descends to

```text
partial I/{+-1}  isomorphic to RP^(k-1).                (32)
```

Restricted to the cubical two-skeleton, this realizes THM-3346's `K_N` as a
weighted projective two-skeleton.  The Veronese map `[u]->uu^T` embeds the
whole `RP^(k-1)` into the symmetric matrices.

There are two different forgetting operations:

1. projectors retain the parent point but forget which raw sign lift was
   chosen;
2. a scalar vertex kernel alone also forgets the cubical cell attachment.

Pulling back the tautological line over the rank-one projector locus and then
taking its unit-sphere bundle recovers the unit-vector double cover.  On
`K_N` this is

```text
Q_k^(2) -> K_N.                                         (33)
```

Its first Stiefel--Whitney class is exactly THM-3346's unique conjugation
class in `H^1(K_N;F_2)`.  If the projector locus and tautological-line sidecar
are discarded inside the ambient Euclidean matrix space, that monodromy is
lost.  A positive scalar kernel is therefore not an owner, orientation, or
ancestry section.

## 7. Consequence and stopping boundary

The theorem proves:

1. an unrestricted PSD U-spine kernel whose entries are signed Gaussian
   content ratios and whose norm is `log C_t`;
2. exact Hilbert distance invoices separating height tails from two content
   channels;
3. the sharp universal `k<=3` survival and a literal rank-four arithmetic
   hostile for the naive projective fold;
4. a graded even-tensor hierarchy from lossless projector to strictly
   positive cosh kernel at every grade;
5. a weighted projective realization of the prime-toggle complex with the
   tautological line as its conjugation sidecar.

It does not prove that the folded metric is CND in arbitrary rank, produce a
Gaussian lift from an unpointed scalar Gram matrix, preserve integral height,
or supply an LRC owner, phase, saturated column map, current, clock packet, or
global exit.  The cosh/projector kernels are diagnostic interfaces, not LRC
certificates.  No Berggren, tournament, Jacobian, or design conclusion follows.

## 8. Provisional exact evidence

The companion audits the unrestricted signed-content Gram law and both
Hilbert distance invoices on all pairs `0<=r,s<=300`; exhausts positive
integer weight rows through eight in ranks at most three; checks balanced
level-four obstructions through rank ten, the Catalan formula through rank 30,
and `294,080` dominant-section Fourier rows; verifies the raw Gram, quotient-
sphere, exact projector and cosh spectra through rank seven; checks the even-
tensor support/rank hierarchy through rank seven across eight cubes and tensor
half-degree four; and reproduces (19)--(21) using formal prime-log coefficient
vectors.
Reproduce with

```bash
python3 04-computation/u_spine_signed_clock_gram_thm3347.py
python3 -O 04-computation/u_spine_signed_clock_gram_thm3347.py
```

The frozen normal, optimized, and stored transcripts byte-match.  The hashes
recorded in frontmatter are for the exact working-tree bytes audited here.
