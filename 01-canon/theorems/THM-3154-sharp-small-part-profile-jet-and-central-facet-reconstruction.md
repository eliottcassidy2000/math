---
id: THM-3154
title: "Sharp small-part profile jet and central-class/facet reconstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For partitions of N, length together with the multiplicities of parts
  1 through floor((N-2)/2) is a faithful code.  The cutoff is sharp among
  initial small-part flags.  Its multivariate endpoint factorial jet recovers
  every partition coordinate and hence every Hasse facet, while an exact
  rational kernel transports the whole jet through virtual pole subtraction.
  The same code reconstructs every symmetric-group conjugacy-class index and
  hence evaluates any specified central class function.  This is a complete
  observer, not a positivity proof or a source of unspecified coefficient
  values.
source: root/multiscale-newton-flag/low-child-flag-extension-2026-08-02
audit: >
  Two independent hostile audits rederived the sharp faithful-code bound and
  witnesses, marked one-letter and virtual-pole kernels, multivariate upper
  binomial inversion, full Hasse-facet typing, conjugacy-class and specified
  central-function scope, and the labelled-charge/Frobenius-shadow boundary.
  Fresh normal and optimized executions both reproduce the stored transcript
  and declared LF hashes exactly.
depends_on:
  - THM-3127-partition-refinement-strassen-upset-dual-and-filter-response
  - THM-3129-bounded-poset-upset-facet-irredundancy
  - THM-3147-length-singleton-endpoint-jet-facet-observer
related:
  - THM-2022-gmc2-frobenius-lowest-balanced-face
  - THM-3112-cycle-weighted-young-subgroup-gap-and-uniform-octopus-boundary
  - THM-3123-heptic-e3-remaining-accessory-classification-and-s7-monodromy
  - THM-3134-tournament-endpoint-jet-and-c3-newton-profile-transform
  - THM-3149-depth-three-selector-persistence-and-cross-support-wall
script: 04-computation/gmc_sharp_small_part_profile_jet_thm3154.py
output: 05-knowledge/results/gmc_sharp_small_part_profile_jet_thm3154.out
script_sha256: dc6516dd7bf666daf2b44516b396382fa8094f61c3e47fd5fb258faadb93928c
output_sha256: 6833ba109704e4029f073122c3467e86b2c7e60da88016f16c4c46f5f408cad2
hash_basis: LF-normalized bytes
---

# THM-3154 -- sharp small-part profile jet and central-class/facet reconstruction

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3147 finds that length and one singleton coordinate expose every
nonprincipal Hasse failure in its finite active-prefix bank.  That empirical
compression has a sharp universal completion.  In fixed degree one need not
mark every part size: length plus an initial half-range of small-part counts
already determines the entire partition.  The resulting endpoint jet is a
complete coordinate chart for Hasse facets, virtual-pole holotopies, and
specified central symmetric-group class functions.

Completeness is deliberately separated from sign.  THM-3129 makes every
nontrivial upset a genuine facet, and the inversion below has alternating
coefficients.  Nothing here turns those reconstructed facet masses positive.

## 1. The marked profile and its exact pole kernel

Fix `k>=0`.  For a partition `lambda`, let `ell(lambda)` be its number of
parts and `m_r(lambda)` its number of parts equal to `r`.  Define

```text
Z_X^[k](t,u,w_1,...,w_k)
 =sum_lambda u^ell(lambda) product_(r=1)^k w_r^m_r(lambda)
              m_lambda[X] t^|lambda|.                         (1)
```

For one letter `x`, put `z=xt`.  The choice of leaving that letter unused or
assigning it one positive exponent gives

```text
L_k(z)
 =1+u sum_(r=1)^k w_r z^r+u sum_(r>k)z^r

 =[1+(u-1)z+u(1-z)sum_(r=1)^k(w_r-1)z^r]/(1-z).              (2)
```

Multiplying `(2)` over the letters proves `(1)` and shows that the series is
group-like in the alphabet.  Therefore one plethystic pole subtraction has
the exact formal kernel

```text
Z_(X-M)^[k]
 =Z_X^[k] (1-Mt)
  /[1+(u-1)Mt
     +u(1-Mt)sum_(r=1)^k(w_r-1)(Mt)^r].                      (3)
```

At `k=0`, `(3)` is the THM-3147 length kernel.  At `k=1`, it specializes to

```text
(1-Mt)/[1+(uw_1-1)Mt+u(1-w_1)M^2t^2],                       (4)
```

exactly the length-singleton kernel of that theorem.

Put `s_0=u-1` and `s_r=w_r-1`.  The denominator in `(3)` becomes

```text
1+s_0Mt+(1+s_0)(1-Mt)sum_(r=1)^k s_r(Mt)^r.                  (5)
```

Multiplying by `(5)` gives a triangular recurrence for every coefficient
`[t^n s_0^j product s_r^alpha_r]`: apart from the leading coefficient, each
term has smaller `n`, smaller total marker order, or both.  Thus the complete
jet can be transported through any finite pole prefix without reconstructing
individual monomial symmetric functions at every step.

The same formula supplies an exact plethystic holotopy chart.  Replacing `M`
by `zM` in `(3)`, with a scalar parameter `z`, transports every marked
coordinate along `X-[zM]`.  This charts a proposed deformation; it does not
assert that the path enters or remains in a positive cone.

## 2. The sharp faithful code

Fix a degree `N>=2` and put

```text
k_N=floor((N-2)/2).                                          (6)
```

Then the code

```text
kappa_N(lambda)
 =(ell(lambda),m_1(lambda),...,m_(k_N)(lambda))              (7)
```

is injective on the partitions of `N`.

### Proof

Suppose `(7)` is known.  Remove the marked parts of sizes at most `k_N`.
Both the number `q` and total size `S` of the remaining parts are known.
Every remaining part is at least `k_N+1`.

If `q=0` or `q=1`, the residual partition is forced.  If `q>=2`, subtract
`k_N+1` from each residual part.  Total surplus zero has only the all-zero
partition, and total surplus one has only one partition.  Two distinct
residual partitions of the same length and sum therefore require surplus at
least two, hence

```text
S>=q(k_N+1)+2>=2k_N+4>N,                                    (8)
```

contrary to `S<=N`.  Thus the residual is unique and `(7)` is faithful.

The cutoff is sharp among initial small-part flags.  If `N=2d>=4`, the two
same-length partitions

```text
(d+1,d-1)                   and                   (d,d)       (9)
```

agree in every `m_r` for `r<=d-2` and first separate at `m_(d-1)`.  If
`N=2d+1>=5`, the corresponding pair is

```text
(d+2,d-1)                   and                 (d+1,d).      (10)
```

Degrees two and three need no marked count: length alone is faithful.  This
minimality concerns the nested family of initial part-size markers.  It does
not rule out a different nonlinear encoding with fewer displayed scalars.

For example, `N=7` universally needs at most

```text
(ell,m_1,m_2),                                               (11)
```

while the special THM-3147 residual bank happens to need only `(ell,m_1)`.
That finite saving is structure of the bank, not a universal coding theorem.

## 3. Endpoint factorial jets recover every Hasse facet

Let `C(lambda)` be any scalar current on the partitions of `N`.  Write

```text
C(u,w)=sum_(lambda partition N)
       C(lambda)u^ell(lambda)product_(r=1)^k w_r^m_r(lambda), (12)
```

where `k=k_N`.  Its endpoint factorial jets are

```text
J_(j,alpha)
 =[s_0^j product_r s_r^alpha_r]
   C(1+s_0,1+s_1,...,1+s_k)

 =sum_lambda C(lambda) binom(ell(lambda),j)
                   product_r binom(m_r(lambda),alpha_r).     (13)
```

Extend `C(v)` by zero from the code image to the finite ambient integer box.
The multivariate binomial inversion formula is

```text
C(v)=sum_(beta>=v)(-1)^(|beta-v|)
      product_i binom(beta_i,v_i) J_beta.                    (14)
```

Indeed `(13)` is the product of the ordinary one-dimensional upper binomial
transforms, and applying binomial inversion in each coordinate proves `(14)`.
By the injectivity of `(7)`, `(14)` recovers every individual `C(lambda)`.

Consequently, for every coarsening upset `U`, its mass is the explicit jet
functional

```text
C(U)=sum_(lambda in U) C(kappa_N(lambda)),                   (15)
```

with each summand expanded by `(14)`.  If `C` has zero total mass, THM-3127
identifies the inequalities `C(U)>=0` with nonnegative upward Hasse transport,
and THM-3129 says every nontrivial one is a distinct facet.  Thus the jet is
a complete observer of the full facet system.

This is a coordinate reconstruction, not a facet compression.  There are
exactly `p(N)` admissible codes because `(7)` is injective, and the alternating
signs in `(14)` provide no automatic positivity.

## 4. Symmetric-group central-class reconstruction

Apply `(7)` to the cycle partition of a permutation `sigma in S_N`.  Then

```text
ell(lambda)=number of cycles of sigma,
m_r(lambda)=number of r-cycles of sigma.                     (16)
```

The faithful code therefore determines the complete conjugacy-class index.
Consequently, once a central class function is specified by a formula or
class table, its value can be evaluated from the code.  The code does not
invent or recover an arbitrary table of coefficient values which was never
supplied.  For the specified cycle-weighted operator of THM-3112,

```text
p_lambda(S)=product_(r>=1) p_r(S)^m_r(lambda)                (17)
```

is exactly evaluable from the code even though only the initial `m_r` are
displayed: the missing large cycle, if present, is reconstructed by `(8)`.
Together with the given alphabets and signed-bank formula, this evaluates the
central coefficients of `Omega_N(S)` and its product-Gamma combinations.

It supplies no positive-semidefinite sign.  THM-3112's equal-support-size
hostile cycle types `(4,1,1)` and `(2,2,1,1)` are separated by `(7)`, but
knowing their distinct coefficients does not build the required positive
Young-subgroup or rooted-octopus transport.

## 5. What is shared, and what is not transported

Two proved finite results occupy very low coordinates of the universal jet:

1. THM-3147's `43` nonprincipal-required active-prefix failures through
   degree nine are all seen by `(ell,m_1)`.
2. THM-3149's six exact depth-at-most-three selector-wall facets also factor
   through `(ell,m_1)`.

These are special low-coordinate compressions.  They do not make
`(ell,m_1)` universally faithful, and the depth-three selector statement is
not a depth-uniform obstruction.

The heptic passport atlas of THM-3123 uses the same partition coordinate and
has the special rule `m_1>=2 iff` the listed cover has monodromy `S_7`.
Outside that proved three-passport atlas, a passport does not generally
determine monodromy.  No cover, accessory parameter, or JC trajectory is
constructed from a GMC current here.

Likewise, reconstruction of a partition or conjugacy class does not recover
labels, incidence, endpoint phase, or a root-sheet coupling.  The theorem
cannot be used as an LRC or tournament bridge when the target varies inside
a profile fibre.

There is a particularly sharp boundary against replacing the successful
GMC2 face argument by a still larger partition jet.  In characteristic `p`,

```text
m_lambda[X]^p=m_(p lambda)[X],                              (18)
```

so the Kummer/Frobenius layer of THM-2022 has an exact unlabeled shadow on
the dilated partition types `p lambda`.  The code recognizes this shadow:

```text
m_r(p lambda)=0 unless p divides r,
m_(pr)(p lambda)=m_r(lambda).                               (19)
```

But sorting a labelled multiplicity vector into a partition erases its
charge assignment.  The minimal hostile is

```text
charges q=(1,-1,2),
s=(1,1,0),              q dot s=0,
t=(1,0,1),              q dot t=3.                          (20)
```

The positive entries of both `s` and `t` have the identical partition
`(1,1)`, hence identical complete profile code, while only `s` is balanced.
Thus no refinement by further unlabeled `m_r` coordinates can recover the
labelled balanced face or its coefficient phase `Q^p`.  The indispensable
THM-2022 sidecar is the support/charge assignment together with the whole-face
coefficient sum, not another partition statistic.

## 6. Exact companion and scope

Run

```text
python 04-computation/gmc_sharp_small_part_profile_jet_thm3154.py
python -O 04-computation/gmc_sharp_small_part_profile_jet_thm3154.py
```

and compare byte-for-byte with

```text
05-knowledge/results/gmc_sharp_small_part_profile_jet_thm3154.out.
```

The companion uses only integer arithmetic.  It checks all `28,627`
partition codes through degree thirty and all `27` sharp cutoff witnesses;
checks the exact marked virtual-pole kernel in degrees through twelve for
cutoffs zero through five; reconstructs `913` arbitrary signed coordinates
by the multivariate endpoint jet through degree sixteen; checks all `98`
partition upsets through degree seven; and independently verifies `913`
conjugacy classes and the symmetric-group class-size mass formula in every
degree through sixteen.

This theorem is a faithful observer and transport kernel.  It does not prove
that a product-Gamma current is Hasse-positive, that any pole-prefix selector
exists, that a central operator is positive semidefinite, or that a positive
holotopy exists.  It proves no new instance of NC2, GMC, LRC(14), JC(2), or
DC(2).

QED.
