---
id: THM-3147
title: "Length-singleton endpoint-jet facet observer and finite portability wall"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  The length endpoint jet gives exact rank-tail Hasse facets and evolves under
  virtual pole subtraction by a triangular Gregory recurrence.  On the full
  THM-3136 active-prefix bank through degree nine it detects exactly 27 of the
  43 failures missed by every principal upset.  The remaining 16 all use one
  N=7 within-length antichain; adjoining only the number of singleton parts
  detects all 43.  This is a finite obstruction observer, not a positive
  current or original-response reconstruction theorem.
source: root/multiscale-newton-flag/low-child-flag-extension-2026-08-02
audit: >
  Two independent hostile audits rederived both marked generating factors,
  the virtual-pole kernels, Gregory recurrence, endpoint/binomial inversion,
  rank-tail and singleton facet geometry, information-loss hostile, and the
  related-only heptic S7/D14 classifier.  Fresh normal and optimized runs
  both reproduce the stored transcript, exact 27+16 census, witnesses,
  counts, and declared LF hashes.
depends_on:
  - THM-3120-row-pole-prefix-newton-flag-positivity
  - THM-3127-partition-refinement-strassen-upset-dual-and-filter-response
  - THM-3129-bounded-poset-upset-facet-irredundancy
  - THM-3136-one-sided-fixed-reference-elementary-tail-hasse-no-go
related:
  - THM-3123-heptic-e3-remaining-accessory-classification-and-s7-monodromy
  - THM-3134-tournament-endpoint-jet-and-c3-newton-profile-transform
  - THM-3137-finite-stochastic-pole-selector-polytope-and-portability-wall
  - THM-3144-mixed-depth-selector-persistence-death-barcode
script: 04-computation/gmc_length_singleton_endpoint_jet_observer_thm3147.py
output: 05-knowledge/results/gmc_length_singleton_endpoint_jet_observer_thm3147.out
script_sha256: 8b1c9572a5b03423a60fa352aecd59cc393152e5991a0ba9a16b210e6a063fb0
output_sha256: 3ee7108bad50bf3182fa54bfa5269af2c6b1b6254bec367604595ae620fb0acb
hash_basis: LF-normalized bytes
---

# THM-3147 -- length-singleton endpoint-jet facet observer and finite portability wall

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3134 shows that a scalar endpoint can hide a complete tournament
path-cover profile.  The analogous sidecar for the THM-3136 Young current is
the profile by number of parts.  Here the sidecar is a detector rather than a
positivity repair: it exposes 27 genuinely nonprincipal Hasse facets, misses
16 shape-within-rank facets, and becomes complete for that finite residual
after one minimal singleton-count refinement.

## 1. The group-like length and singleton profiles

For an alphabet `X`, define

```text
Z_X(t,u)=sum_(N,lambda partition N)
         u^ell(lambda) m_lambda[X] t^N
        =product_(x in X) (1+(u-1)xt)/(1-xt).                (1)
```

The product identity follows by either leaving a letter `x` unused or
assigning it one positive exponent.  It extends group-likely to virtual
alphabets.  At `u=1`,

```text
Z_X(t,1)=H_X(t).                                             (2)
```

For one virtual pole `M`, division by the one-letter factor gives the exact
profile transform

```text
Z_(X-M)(t,u)=Z_X(t,u) (1-Mt)/(1-(1-u)Mt).                   (3)
```

To retain the smallest statistic needed below, mark singleton parts as well:

```text
Z_X(t,u,w)=sum_lambda
 u^ell(lambda) w^m_1(lambda) m_lambda[X] t^|lambda|

=product_(x in X)
 [1+(uw-1)xt+u(1-w)x^2t^2]/(1-xt).                         (4)
```

Consequently

```text
Z_(X-M)(t,u,w)
=Z_X(t,u,w) (1-Mt)
 /[1+(uw-1)Mt+u(1-w)M^2t^2].                              (5)
```

These are identities in formal power series.  No positivity assumption is
being made.

## 2. Endpoint jets and rank-tail facets

Put `s=u-1` and define the endpoint Taylor array

```text
A_(N,j)(X)=[t^N s^j]Z_X(t,1+s)
           =e_j[X]h_(N-j)[X].                               (6)
```

Multiplying `(3)` by its denominator gives a triangular Gregory recurrence:

```text
A'_(N,j)=A_(N,j)-M A_(N-1,j)-M A'_(N-1,j-1),               (7)
```

where primes mean `X-M` and out-of-range entries are zero.  Iterating `(7)`
transports the full length jet through any pole prefix.  Equation `(5)` gives
the analogous triangular two-variable transport.  Around `u=w=1`, writing
`s=u-1,r=w-1`, its denominator is

```text
1+sMt+rMt(1-Mt)+srMt(1-Mt).                                 (8)
```

Now let `C_N(lambda)` be any zero-mass current on partitions of `N`, and set

```text
C_N(u)=sum_lambda C_N(lambda)u^ell(lambda),
J_j=C_N^(j)(1)/j!.                                          (9)
```

For `1<=m<N`, the rank tail

```text
U_m={lambda:ell(lambda)<=m}                                 (10)
```

is a coarsening upset.  Binomial inversion at the endpoint gives its exact
mass:

```text
C_N(U_m)=sum_(j=m+1)^N
 (-1)^(j-m) binom(j-1,m) J_j.                              (11)
```

This is not merely a necessary inequality chosen for convenience.  The
induced Hasse graph on `U_m` is connected through the top partition `(N)`,
and its complement is connected through `(1^N)`.  THM-3129 therefore makes
every nontrivial `(10)` an extreme upset ray, equivalently a facet of the
Hasse boundary cone.

For the THM-3136 fixed-reference current, the endpoint entries in `(9)` are

```text
J_j=Phi^R(h_N)e_j[Q]h_(N-j)[Q]
    -h_N[Q]Phi^R(e_j h_(N-j)).                              (12)
```

Thus `(11)` is an exact low-dimensional observer of a genuine Young-poset
facet.

## 3. Exact finite census: length catches 27 of 43

The companion scans all `8,241` active prefixes, both banks on all `115`
THM-3120 supports, and degrees `5<=N<=9`.  It checks `(7)` and `(5)` at every
one of the `8,011` successive pole removals.  The exact principal-upset and
rank-tail counts are

```text
N       5       6       7       8       9
principal bad
        0    5641    5386    6510    6767
rank-tail bad
        0    5641    5278    6451    6647.                 (13)
```

Among currents passing every principal upset, rank tails catch

```text
N=7: 19,                 N=9: 8,                 total: 27. (14)
```

THM-3136 proves that exactly `43` Hasse failures pass every principal upset,
so length alone is a strict but incomplete compression.

The first length-jet witness is the THM-3136 case

```text
N=7, (a,b)=(1,2), bank I1, R=(5), m=5.                     (15)
```

Its rank-tail mass is

```text
-1151498720.                                                 (16)
```

In endpoint coordinates, `(11)` is the two-term identity

```text
-J_6+6J_7=-1151498720.                                     (17)
```

The apparently complicated nonprincipal antichain in THM-3136 is exactly
the rank-two tail `ell(lambda)<=N-2`.

## 4. One singleton coordinate detects the remaining 16

Every one of the `16` principal-pass, rank-tail-pass failures occurs at
`N=7` and has the same minimal antichain

```text
{(4,1,1,1),(3,2,1,1)}.                                     (18)
```

Its upset is

```text
V={lambda:ell(lambda)<=3}
  union {lambda:ell(lambda)=4 and m_1(lambda)>=2}.           (19)
```

This is visible coefficientwise in `(4)`.  It is also a genuine facet:
`V` and its complement are induced-Hasse connected, and `(18)` is its exact
minimal antichain.  The exact census gives

```text
length observer catches                 27,
singleton-refined observer catches      16,
total                                   43.                 (20)
```

Thus the two-variable profile detects every finite nonprincipal-required
failure from THM-3136 through degree nine.

The first new witness is

```text
N=7, (a,b)=(2,10), bank I2, prefix length 24,               (21)

C_7(V)=-2816353101792484700160000.                          (22)
```

Every rank tail and every principal upset passes there.  The need for one
extra coordinate is information-theoretic: the zero-mass current

```text
delta_(4,1,1,1)-delta_(2,2,2,1)                            (23)
```

has identically zero length profile, but `V` assigns it mass `1`.  Marking
`m_1(lambda)` separates the two shapes.  This proves minimality only relative
to statistics that depend on length alone; it does not assert a unique or
globally minimal refinement of the full partition profile.

### 4.1 A heptic monodromy interpretation of the new coordinate

The singleton statistic is not an arbitrary separator on the degree-seven,
length-four slice.  Independently, THM-3123 classifies the complete heptic
`e=3` accessory atlas indexed by exactly the same three partitions:

```text
passport          m_1        unmarked monodromy
(4,1,1,1)          3                 S_7
(3,2,1,1)          2                 S_7
(2,2,2,1)          1                 D_14.                  (23a)
```

Therefore, on that proved passport atlas,

```text
m_1>=2   iff   the unmarked cover has full monodromy S_7.   (23b)
```

The `w`-coordinate in `(4)` is exactly the symmetry-breaking coordinate
separating the two full-symmetric strata from the Chebyshev/dihedral wall.
This is a shared partition-indexed classifier, not a map from the GMC Young
current to a response cover: no monodromy, accessory parameter, Keller-chart
entry, or JC consequence is transported between the two problems.

## 5. Positivity and partial-reference boundaries

The full jet is not automatically positive when its scalar endpoint is.  At
the selected `(a,b)=(1,2),I1,R=(5),N=6` row, the signed length profile is

```text
20u^3+368u^4+264u^5-40u^6,                                 (24)
```

while its endpoint Taylor coefficients are

```text
(612,2612,4308,3332,1088,24,-40).                           (25)
```

The scalar endpoint `612` is positive, but the top jet coefficient is
negative.  Hence `(7)` transports information, not a coefficientwise
positive cone.

For the normalized current, that same top coefficient is the singleton
partition coordinate.  THM-3136 strengthens the obstruction along the full
partial-reference interval `Q_z=Q-[5z]`: its degree-six bottom coordinate has
positive Bernstein coefficients and stays at least

```text
1760080                                                     (26)
```

for every `0<=z<=1`.  Thus no point of that reference holotopy crosses the
rank-one tail facet into the Hasse cone.

## 6. Scope and portability wall

The theorem is exact for the complete THM-3136 bank through degree nine.
It says that `(u,w)` is complete for the *43 already identified finite
nonprincipal-required failures*.  It does not say that length plus singleton
count determines every partition upset, remains complete in degree ten or
beyond, or proves positivity when all displayed observers pass.

The profiles are derived signed-bank sidecars.  Neither `(3)` nor `(5)`
reconstructs the original nonrow normalized product-Gamma response, and no
probability law over prefixes is supplied.  THM-3144's mixed-depth selector
barcode is a distinct convex persistence question.  No Gaussian moment
conjecture, NC2, LRC(14), JC(2), or DC(2) is proved here.

## 7. Exact companion

Run

```text
python 04-computation/gmc_length_singleton_endpoint_jet_observer_thm3147.py
python -O 04-computation/gmc_length_singleton_endpoint_jet_observer_thm3147.py
```

and compare byte-for-byte with

```text
05-knowledge/results/gmc_length_singleton_endpoint_jet_observer_thm3147.out.
```

The companion uses only integer/rational arithmetic.  It reconstructs every
profile from monomial symmetric functions, checks both pole kernels on every
successive prefix, verifies the facet connectivity, and performs the exact
`27+16` census.

QED.
