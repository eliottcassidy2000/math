---
id: THM-3713
title: "LRC orbitwise deep-offset detector and successor quotient hostile"
status: >
  PROVED + FINITE-EXACT + VERIFIED.  Fiberizing the THM-2365 tensor by its
  deep offset u=r-t turns target drift into the average of thirteen ordinary
  two-torus variances.  The two coordinate-edge banks and, separately, the
  THM-3665 three-site mask on every offset vanish exactly when target drift
  vanishes, with sharp frame bounds and explicit one-defect tariffs.  For a
  rational lawful tensor, one nonzero orbitwise defect forces a nonzero
  target coefficient in every one of the twelve nonzero deep residues.  The
  successor marginal is their sum over u and can cancel: an explicit
  nonnegative diagonal-zero tensor has constant successor marginal but
  D_H=E_dt=4/169 and 78 nonzero orbitwise defects.  On the typed THM-3672
  non-cover control, every chart has common active offsets {1,2,10,11,12}
  with the same ++--- sign pattern.  This proves no cover-row nonvanishing
  and no LRC(14) conclusion.
source: codex-lrc14-20260822 / deep-offset fiberization and typed control
audit: >
  PASS -- independent proof audits verified the offset isometry, both sharp
  frames, one-defect tariffs, successor quotient, anchored-zero cyclotomic
  saturation, exact calH/B target typing, and THM-2365 extraction scope.  An
  alternate whole-period interval reconstruction reproduced all 390 typed
  defects, their digest, five-colour signs, two-site signs, strongest cases,
  tariffs, all thirteen mod-79 transforms, and the balanced-moment range.
  Normal and optimized transcripts equal the stored output; both raw-LF
  hashes match this metadata.  No cover-row or LRC(14) overclaim was found.
depends_on:
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-3665-lrc-support-minimal-three-twist-target-detector
  - THM-3672-lrc-successor-mass-all-packet-positive-control
  - THM-3674-sharp-successor-variance-drift-and-target-energy-tariff
related:
  - THM-2350-owner-pivot-dual-dipole-normal-form
  - THM-2536-deep-comb-selector-flow-target-drift-and-centered-toothpick-boundary
  - THM-2580-hasse-bockstein-carry-tower-and-salem-local-unit-boundary
  - THM-3662-lrc-eleven-cell-exceptional-flux-and-high-digit-variation-gate
  - THM-3701-lrc-radial-successor-mass-gate-and-star-frame
  - THM-3710-lrc-successor-endpoint-182-skeleton-bad-prime-collapse
  - THM-3718-lrc-complete-atom-orbit-defect-saturation-and-semantic-boundary
script: 04-computation/lrc_orbitwise_deep_offset_detector_thm3713.py
output: 05-knowledge/results/lrc_orbitwise_deep_offset_detector_thm3713.out
script_sha256: 178c434f3d37b186d105495250f836c234668524d8dbb5275af5f12b6c2c4af2
output_sha256: af366cc70a5ecf7cd9bc97737f9c0da46590a09984619c783632d020a2f372a7
hash_basis: raw LF bytes
---

# THM-3713 -- retaining the deep offset makes the target detector complete

**PROVED + FINITE-EXACT + VERIFIED + INDEPENDENTLY AUDITED.**  THM-3701
reduces the successor screen to thirteen scalar masses, but that screen has
already summed away the deep colour.  This theorem restores precisely that
lost coordinate.  The result is a complete local detector for THM-2365
target drift and an exact account of how the successor screen can miss it.

## 1. Deep-offset fiberization is an isometry

Fix `p=13` and let

```text
H:F_p^3 -> C
```

be a THM-2365 tensor.  For each deep offset `u` put

```text
h_u(s,t)=H(t+u,s,t),
bar(h_u)=p^-2 sum_(s,t) h_u(s,t).                   (1)
```

THM-2365's target-action projection is

```text
(P H)(r,s,t)
 =p^-2 sum_(a,b) H(r+b,s+a,t+b).
```

At `r=t+u`, the substitutions `s'=s+a,t'=t+b` give the exact identity

```text
(P H)(t+u,s,t)=bar(h_u).                            (2)
```

Consequently the target-drift energy fiberizes as

```text
D_H
 =p^-3 sum_(r,s,t)|H(r,s,t)-(P H)(r,s,t)|^2
 =p^-3 sum_(u,s,t)|h_u(s,t)-bar(h_u)|^2.           (3)
```

Equivalently, for the normalized two-torus variance
`Var(h_u)=p^-2 sum_(s,t)|h_u-bar(h_u)|^2`, equation (3) reads
`D_H=p^-1 sum_u Var(h_u)`.

Thus `D_H=0` exactly when every offset slice is constant, equivalently

```text
H(r,s,t)=G(r-t).                                   (4)
```

No positivity or diagonal-zero hypothesis is needed for (1)--(4).

## 2. Two edge banks give a complete rational detector

All indices below are cyclic modulo thirteen.  Define

```text
E_orb=p^-3 sum_(u,s,t) (
 |h_u(s+1,t)-h_u(s,t)|^2
 +|h_u(s,t+1)-h_u(s,t)|^2).                        (5)
```

This is the ordinary square-torus Dirichlet energy on every offset slice.
It is the ungauged version of the magnetic Cayley-square mechanism in
THM-2350.  Its character `(a,b)` has eigenvalue
`4sin(pi*a/13)^2+4sin(pi*b/13)^2`; the least nonzero and largest eigenvalues
are respectively `4sin(pi/13)^2` and `8cos(pi/26)^2`.  Hence the sharp frame

```text
4 sin(pi/13)^2 D_H
 <=E_orb
 <=8 cos(pi/26)^2 D_H.                             (6)
```

In particular, the two full coordinate-edge banks vanish exactly when
`D_H=0`.  A single nonzero edge

```text
delta=h_u(s+1,t)-h_u(s,t)
```

or its `t`-direction analogue already gives, by the norm-two Cauchy bound,

```text
D_H>=|delta|^2/(2p^3)=|delta|^2/4394.              (7)
```

If `H` also has THM-2365's lawful diagonal zero `H(t,s,t)=0`, its
nonzero-deep, nonzero-target energy satisfies `E_dt>=D_H/p`; hence

```text
E_dt>=|delta|^2/(2p^4)=|delta|^2/57122.            (8)
```

A single two-site difference is not a universal convolution detector: it
has a one-dimensional character kernel.  There is no conflict with
THM-3665's three-site support minimum.  Completeness in (6) comes from using
both coordinate directions as a joint bank.

## 3. The support-minimal detector works slice by slice

Apply the THM-3665 mask to every retained offset:

```text
M_u(s,t)
 =h_u(s,t)+h_u(s+1,t)-2h_u(s,t+1).                 (9)
```

Its sharp frame, summed over `u`, is

```text
16 sin(pi/13)^4 D_H
 <=p^-3 sum_(u,s,t)|M_u(s,t)|^2
 <=16 cos(pi/26)^2 D_H.                            (10)
```

Therefore

```text
all M_u(s,t)=0
 iff every h_u is constant
 iff D_H=0.                                        (11)
```

One defect `Delta=M_u(s,t)` has coefficient norm squared six, so

```text
D_H>=|Delta|^2/(6p^3)=|Delta|^2/13182,             (12)

E_dt>=|Delta|^2/(6p^4)=|Delta|^2/171366            (13)
```

under the lawful diagonal-zero hypothesis.  Equations (9)--(13) are a
direct transposition of THM-3665 from one target profile to the thirteen
offset fibers; no new Fourier nonvanishing assumption is inserted.

## 4. The successor screen is a cancellation quotient

The successor marginal is

```text
S(s,t)=sum_r H(r,s,t)=sum_u h_u(s,t).               (14)
```

Its coordinate differences and three-site defect obey

```text
S(s+1,t)-S(s,t)
 =sum_u (h_u(s+1,t)-h_u(s,t)),

S(s,t+1)-S(s,t)
 =sum_u (h_u(s,t+1)-h_u(s,t)),

S(s,t)+S(s+1,t)-2S(s,t+1)
 =sum_u M_u(s,t).                                  (15)
```

Thus successor variation is a sufficient drift certificate, but the
converse fails because (15) forgets `u`.  The source is the collection of
thirteen offset gradients; the target is their sum.  Target nonconstancy is
preserved when the sum is nonzero, while cancellation through the pointwise
twelve-dimensional offset-sum kernel `ker(sum:C^13->C)` destroys the
converse.  THM-3701's thirteen sampled successor masses are a still coarser
shadow of this quotient.

### 4a. One rational orbitwise defect saturates the deep colours

Assume now that `H` is rational-valued and has the lawful diagonal zero.
Then `h_0(s,t)=0` identically, so also

```text
M_0(s,t)=0.                                        (15a)
```

Fix one centre `(s0,t0)` and put `d_u=M_u(s0,t0)`.  If `d` is nonzero, then

```text
d_hat(a)=sum_u d_u zeta^(a u) !=0
                  for every a!=0.                 (15b)
```

Indeed, if one value in (15b) vanished, the rational polynomial
`sum_(u=0)^12 d_u X^u` would be divisible by
`Phi_13(X)=1+X+...+X^12`.  Since both degrees are at most twelve, all
thirteen coefficients `d_u` would be equal.  The anchored zero `d_0=0`
would force `d=0`, a contradiction.  This is the anchored-zero cyclotomic
mechanism of THM-2536, now applied to the deep-offset defect.

The character typing is exact.  Define

```text
calH(a,b,q)
 =p^-3 sum_(u,s,t)h_u(s,t)zeta^(a u+b s+q t)
 =B(a,b,q-a),                                      (15c)

mu(b,q)=1+zeta^(-b)-2zeta^(-q).
```

Fourier inversion in the two target variables gives

```text
d_hat(a)
 =p sum_(b,q) mu(b,q)calH(a,b,q)
      zeta^(-b s0-q t0).                           (15d)
```

THM-3665 says `mu(b,q)` vanishes exactly at the zero target `(b,q)=(0,0)`.
Therefore (15b)--(15d) prove the simultaneous landing

```text
for every a!=0, some (b,q)!=(0,0) has
B(a,b,q-a)!=0.                                     (15e)
```

Thus one rational orbitwise defect supplies a nonzero finite target
coefficient in every nonzero deep residue, not merely in one residue.  When
`H` is an actual THM-2365 interval tensor, that theorem's absolutely
convergent extraction further gives, for each `a!=0`, some exact ordinary
frequency and some integer deep frequency congruent to `a mod 13` and coprime
to `91`.  The target, ordinary frequency, and integer representative may vary
with `a`; no common triangle or all-colour magnitude floor is asserted.

## 5. A bounded lawful-array hostile

The failure is not caused by signs or large coefficients.  Put

```text
f(0)=1, f(1)=-1, f(s)=0 otherwise,
g(1)=1, g(2)=-1, g(u)=0 otherwise,

H(r,s,t)=0                         if r=t,
H(r,s,t)=2+f(s)g(r-t)              if r!=t.         (16)
```

Then `H` is nonnegative, takes values in `{0,1,2,3}`, and has the exact
lawful diagonal zero.  Since `sum_u g(u)=0`,

```text
S(s,t)=24
```

for all 169 pairs, so every successor edge and every successor three-site
defect vanishes.  Nevertheless only the `u=1,2` slices vary and exact
counting gives

```text
D_H=4/169.                                         (17)
```

THM-3674's split `D_H=Var(S)/p^2+E_dt` and `Var(S)=0` give

```text
E_dt=D_H=4/169.                                    (18)
```

Here

```text
M_u(s,t)=g(u)(f(s+1)-f(s)).                         (19)
```

There are exactly 78 nonzero orbitwise defects, with values
`{-2,-1,1,2}`, and

```text
E_orb=p^-3 sum|M_u|^2=12/169.                      (20)
```

They form 39 nonzero offset profiles, and (15e) makes every one of those
profiles live in all twelve nonzero deep residues despite its zero successor
sum.  The companion verifies the same saturation in an exact order-13
subgroup of `F_79^*`.

This is a minimal-mechanism refutation of

```text
"nonnegative + diagonal zero + constant successor marginal
 implies zero target drift."                       (21)
```

The strongest survivor is the exact opposite energy attribution: when the
successor is constant, all drift lies in `E_dt`.  The tensor is an abstract
lawful-array hostile, not an interval factorization, a scalar cover, or a
counterexample to LRC(14).

## 6. The typed THM-3672 control reveals five offset colours

The companion next reconstructs THM-3672's pinned typed non-cover row on its
common numerator grid `N=50334435734703120`.  For each ordered legal chart
`(k,l)`, put `H^num_(k,l)=N H_(k,l)` and write

```text
h^0_u       =H^num_(k,l)(u,0,0),
h^A_(k,u)   =H^num_(k,l)(u,1,0),
h^B_(l,u)   =H^num_(k,l)(u+1,0,1).                 (22)
```

The numerator-valued offset defect at the three sampled twists is

```text
Delta_(k,l)(u)
 =h^0_u+h^A_(k,u)-2h^B_(l,u).                      (23)
```

The shift `u+1` in the last term is forced by the definition `u=r-t`, not a
choice of gauge.  Summing (23) over `u` reproduces THM-3672's exact successor
defect numerator.

All thirty charts have precisely the common nonzero offset support

```text
{1,2,10,11,12}={1,2,-3,-2,-1} in F_13.            (24)
```

Moreover, every chart has the same oriented sign pattern:

```text
Delta_(k,l)(1), Delta_(k,l)(2) >0,
Delta_(k,l)(10),Delta_(k,l)(11),Delta_(k,l)(12)<0. (25)
```

Let `u_tilde` be the balanced lift in `{-6,...,6}` and define the oriented
first moment

```text
L_(k,l)=sum_u u_tilde Delta_(k,l)(u).               (25a)
```

Every term on the five-point support is positive by (24)--(25), so this one
scalar cannot cancel under the cross-offset cancellation that mixes opposing
signs in the successor sum.  Exact reconstruction gives

```text
430999196213378
 <=L_(k,l)<=441421418113351                         (25b)
```

on all thirty charts.  Modulo thirteen, (25a) is the same first binomial
coefficient functional that appears in THM-2580; its Cayley-filling
interpretation there additionally requires augmentation zero, which these
typed defects do not have.  Over the integers it is an oriented dual
functional, analogous to THM-3662's reversal-odd flux.  The map from the
five-colour vector to `L_(k,l)` preserves strict positivity on the cone
(24)--(25) but destroys the individual colours and magnitudes.  Its sign uses
the physical origin `u=0` and the balanced lift, so it is not licensed after
an arbitrary rotation or orientation-forgetting quotient.

The finer two-site edges are uniform over graft labels as well:

```text
h^A_(k,u)-h^0_u: negative exactly at {1,2,11,12},
h^B_(l,u)-h^0_u: negative at {1,2}, positive at {10,11,12}. (26)
```

The three largest orbitwise defects occur at

```text
(k,l,u)=(5,0,11),(5,1,11),(5,2,11)
```

and have common numerator `-66828200140260` over `N`, reducing to

```text
Delta=-28559059889/21510442621668.                 (27)
```

Using the normalized defect in (27), equations (12)--(13) give the exact
positive-control tariffs

```text
D_H >=815619901743488692321
      /6099300086944899889559253516768
     =1.3372352403...*10^-10,

E_dt>=815619901743488692321
      /79290901130283698564270295717984
     =1.0286424925...*10^-11.                      (28)
```

These are about `124.04` times the strongest successor-only tariffs in
THM-3672 on the same three charts.  The gain is exactly what (15) predicts:
keeping one offset before summation avoids cross-offset cancellation.

Statements (22)--(28) are **FINITE-EXACT** on one typed positive non-cover
row.  In particular, (15e) holds for all twelve deep residues on each of its
thirty charts.  The common five-colour support and first-moment sign are not
asserted for a cover.

## 7. Preserved data and the narrowed cover target

Fiberization preserves the deep offset, both target-twist coordinates,
amplitudes, signs, and the exact drift norm.  Successor summation preserves
only the total over offsets and destroys which deep colour paid or cancelled
the defect.  Neither representation retains by itself the endpoint owner,
the delayed-word ancestry address, visible height, terminal phase, or the
all-nine-coordinate `91`-unit projector.

THM-3718 subsequently composes THM-2445/2452 with this theorem.  For any
hypothetical scalar-cover realization among the current 165 rows, it forces
both detector banks to fire on one adaptively selected complete matched
Boolean atom, and transfers a defect to a literal
occupied-to-excluded-target selector class.  That atom is not identified
with the canonical exclusive owner, semantic word/repair, physical root
section, or prescribed first-expiration clock.

The resulting cover-specific target is therefore precise:

```text
for a genuine scalar cover with an eligible low target,
transport one THM-3718 adaptive selector defect, target character and
deep residue to the canonical owner/root/word packet, or independently
force a defect in that canonical table (equivalently close its
circulant branch) with the missing semantic sidecars.              (29)
```

The five-colour typed control suggests that a signed offset-residue sidecar
can be much smaller than the full `13^3` tensor.  THM-3718 now supplies an
adaptive cover-derived sidecar, but not the canonical semantic transport in
(29), which remains open.  LRC(14) is not proved.

The standard-library-only companion pins THM-3672 and its interval engine,
checks the bounded hostile exactly, reconstructs all 390 typed coloured
defects, and stores a digest of the complete ledger.

```bash
python -B 04-computation/lrc_orbitwise_deep_offset_detector_thm3713.py
python -B -O 04-computation/lrc_orbitwise_deep_offset_detector_thm3713.py
```

**QED.**
