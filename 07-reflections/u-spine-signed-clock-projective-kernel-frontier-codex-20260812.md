# Signed U-spine clocks: Hilbert structure, projective curvature, and prime-toggle sidecars

**Research synthesis, 2026-08-12.**  THM-3347 is the proof source for the
exact statements summarized here.  This note records the concept board and
new questions.  It is not an LRC certificate or an additional proof source.

## Outcome first

The U-spine now has a useful analytic coordinate that is genuinely arithmetic.
For the consecutive Gaussian spinor

```text
z_t=(t+1)+it,       C_t=N(z_t)=2t^2+2t+1,
Phi(z_t)=(2t+1,2t(t+1),C_t),       inradius=t,
```

attach one signed coordinate for every split-prime Hensel layer dividing
`C_t`.  Equal local roots have the same sign; reflected roots have opposite
sign.  The resulting Hilbert inner product is exactly

```text
Lambda(r,s)=log gcd(C_r,s-r)-log gcd(C_r,r+s+1).
```

So a single object simultaneously records:

- a primitive Pythagorean triple and its inradius;
- a sum-of-two-squares norm;
- same/opposite Gaussian prime allocation; and
- a positive semidefinite Gram coordinate.

The subtlety begins at the fixed-hypotenuse parent quotient.  Conjugation is
the antipode, so the signed Gram does not descend.  Taking absolute values
looks natural but is wrong in rank four: at the least four-prime grade
`N=5*13*17*29`, the absolute kernel is indefinite and the folded content
metric is not of negative type.  Two repairs survive:

1. the rank-one projector kernel `Lambda^2`, which is low-rank, Hilbertian,
   and lossless for pairwise folded distance; and
2. the full-rank kernel `cosh(tau Lambda)`, which is strictly positive
   definite and reconstructs the complete unordered divisor pair.

The projector does not abolish the conjugation topology.  Its tautological
line has the original raw sign cube as unit-sphere cover, and its first
Stiefel--Whitney class is THM-3346's surviving `Z/2` monodromy.

## Updated Anchor / Niche / Wildcard portfolio

| Lane | Exact object | Gain | Boundary |
|---|---|---|---|
| Anchor -- LRC(14) | labelled determinant/Kelvin decks | a new PSD arithmetic feature family and explicit quotient hostile | no map preserving saturation, owner, phase, or first exit |
| Niche -- primitive triples | U-spine spinors with Hensel-layer signs | exact content Gram and two Hilbert distance invoices | signed lift is lost by conjugation |
| Wildcard -- projective topology | weighted sign cube in `RP^(k-1)` | Veronese and cosh repairs; tautological cover recovers `w_1` | scalar vertex matrices alone forget cell attachment |

The niche again overtook the anchor.  Its leverage for LRC is a design rule:
Hilbertization can make a carrier easier to compare without making it a
certificate, and entrywise absolute value can destroy precisely the positivity
one hoped to exploit.

## Concept board after the Gram pull

| Object | Representation | Invariant | Operation | Lost coordinate / next test |
|---|---|---|---|---|
| U-spine triple | `Phi(t+1,t)` | `C_t=4T_t+1`, radius `t` | Gaussian product | unordered legs lose spinor gauge |
| prime clock | `sigma_(p,e)(t)` | Hensel branch at every valuation layer | compare two depths | unsigned valuation loses branch |
| raw clock vector | `Psi(t)` | norm `log C_t` | Hilbert sum/difference | no fixed consumer owner |
| fixed-grade cube | weighted signs `sqrt(log q_j) epsilon_j` | signed content ratio | prime toggle | conjugation negates the vector |
| parent metric | `min(h,W-h)` | smaller complementary divisor | antipodal fold | not CND in rank four |
| projector | `uu^T` | squared correlation, levels `0,2` | Veronese | lift sign and cell structure absent from scalar kernel |
| cosh kernel | `cosh(tau Lambda)` | every even Fourier character | Schur/exponential repair | full rank but still no physical ownership |
| tautological line | line over projector locus | cover class `w_1` | choose a section | no global section on the prime-toggle complex |

## Why signed content is PSD

This positivity is not a numerical accident.  Each prime-power layer is a
rank-one feature.  A layer common to `C_r` and `C_s` contributes `+log p` when
the roots agree and `-log p` when they are reflected.  Summing the layers gives
the exact two-channel content ratio.  Consequently

```text
||Psi(r)-Psi(s)||^2
 =log((C_r/g)(C_s/g) gcd(C_r,r+s+1)^4),
||Psi(r)+Psi(s)||^2
 =log((C_r/g)(C_s/g) gcd(C_r,s-r)^4),
```

where `g=gcd(C_r,C_s)`.  The two distances separate unshared height from
shared orientation.  This is stronger than an unsigned log-gcd kernel: the
sign tells which Gaussian composition absorbs a shared factor.

## Why the fold first fails at four prime directions

At fixed grade, the raw inner product is `W-2h`, where `h` is weighted Hamming
distance and `W=log N`.  Parent distance replaces it by

```text
d=min(h,W-h)=(W-|Lambda|)/2.
```

Ranks at most three have no dangerous even Fourier layer beyond pairs; the
absolute kernel remains PSD for all weights.  At rank four a level-four mode
appears.  The arithmetic witness at `N=32045` gives a centered vector `c` with

```text
c^T d c=16 log 5>0,
c^T |Lambda| c=-32 log 5<0.
```

This is the exact first failed implication:

```text
signed Gram PSD  does not imply  entrywise-absolute Gram PSD,
antipodal quotient of a Hilbert cut metric need not have negative type.
```

The strongest simple survivor is a hemisphere condition.  If one log-prime
weight dominates all the others combined, orient every parent by that sign;
the fold becomes ordinary weighted Hamming distance.

## Two repairs encode different amounts of information

The projector kernel has the expansion

```text
Lambda(x,y)^2=sum w_j^2+2 sum_(i<j) w_i w_j chi_ij(x-y).
```

It retains only Fourier levels zero and two, yet the nonlinear chord transform
is strictly increasing in the parent distance.  It is therefore pairwise
lossless while occupying only `binom(k,2)` centered Euclidean dimensions.

The cosh kernel instead retains every even character:

```text
lambda_A=2^(k-1) product_(j in A)sinh(tau w_j)
                    product_(j notin A)cosh(tau w_j)>0.
```

It is full rank on the parent torsor and reconstructs both complementary
content divisors.  This gives a useful hierarchy:

```text
signed lift -> raw linear Gram
parent point -> low-rank projector Gram
full parent function algebra -> cosh Gram.
```

The interpolation is exact.  The degree-`2m` symmetric-tensor kernel
`Lambda^(2m)` has precisely the even character levels `0,2,...,2m`; its rank is

```text
sum_(j=0)^min(m,floor(k/2)) binom(k,2j).
```

It becomes strictly positive definite first at `m=floor(k/2)` for `k>=2`.
Each Fourier coefficient counts weighted length-`2m` prime-clock words by the
parity set of coordinates used oddly.  Cosh is the positive generating
function of this filtration.  Thus higher even degree has a concrete meaning:
it pays for larger interacting subsets of prime toggles rather than merely
adding opaque nonlinear features.

No member of the hierarchy restores an owner or physical phase by itself.

## Frontier experiments

### 1. Kernelized LRC hostile before any positive claim

Source: a lawful saturated 13-column coefficient plane.  Map each column or
parameter direction to U-spine clock features only when an exact integral
spinor correspondence exists.  Preserved predicate: full determinant and
content data.  Destroyed information: owner, endpoint, and row height.  The
cheapest decisive hostile is to compare a plane with its separately normalized
Gaussian image from THM-3336.  A Gram improvement must not be confused with
the sufficient determinant gate or the actual lonely-runner margin.

### 2. Tautological-section obstruction as a consumer test

Source: `K_N` with projector coordinates.  Target: any proposed sign/owner
bundle supplied by a frontier problem.  Ask whether the consumer pulls back a
section of the tautological line.  At rank three the test is the cubical
`RP^2`; a claimed global orientation must fail on its unique `w_1` loop unless
the consumer provides new side data.  This is sharper than saying merely that
conjugation was forgotten.

### 3. Compare low-rank and full-rank repairs on Pell bridge grades

THM-3346 supplies Pell grades with unbounded prime-toggle rank.  Evaluate the
projector and cosh kernels on their two distinguished parents and adjacent
U-spine outputs.  The projector should retain their folded distance; the cosh
kernel should reconstruct `{M_-,M_+}`.  The open question is whether any
natural recurrence acts sparsely on even Fourier characters rather than on
the source-dependent Berggren tree.

### 4. Curvature lives at Fourier level four

The failure mode suggests a controlled higher-order invariant.  Level two is
the projector/Hilbert part; level four is the first projective-fold curvature.
Test whether the level-four coefficient has a direct count of interacting
prime-toggle squares or a signed two-chain interpretation.  A positive answer
would connect THM-3346's growing `H_2` to the analytic obstruction without
pretending that ambient Euclidean homology sees it.

## Stopping boundary

The session found a precise positive structure, a minimal negative boundary,
and two exact repairs.  It did not solve LRC(14), produce a source-independent
Berggren action, or turn prime allocations into a tournament.  The durable
architecture is

```text
signed prime clock + projective parent + tautological lift sidecar.
```

Dropping the sign too early breaks positivity; retaining only a positive
kernel drops topology; retaining the projective locus and its line bundle
keeps the conjugation obstruction honestly typed.
