---
id: THM-2524
title: "Translated chi7 Hamilton polarization and lossless centred-collision inversion"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  Retaining all thirteen
  relative translations of one chi_7 Hamilton polarization repairs the
  isotropic-diagonal failure of THM-2523.  If B is THM-2519's rational
  collision profile and
  b=B-average(B), then for every nonzero K_14 slope tau the translated bank
  R_tau(t) is exactly 13 A_tau b.  The operator A_tau is invertible on the
  twelve-dimensional centred space, with inverse
  (A_tau^5-39A_tau^3+299A_tau)/325.  Hence R_tau recovers b and the drift
  b_0, and positive rational drift gives all twelve nontrivial root modes
  even when the diagonal contrast R_tau(0) vanishes.  Self-correlation is
  still antipodally even.  For a supplied ordered cross-cospan the same map
  preserves its odd component if and only if that component was already
  present.  Applied to THM-2522, every live THM-2349 word event has this bank
  at its uniform first collision L=1 while its shallow and old-deep banks
  remain separately retained.  The bank is signed and recovers an
  autocorrelation, not the predecessor phase, a Boolean owner/source event,
  or an emitted scalar row; it does not prove LRC(14).
source: codex-2026-07-27-translated-chi7-polarization
depends_on:
  - THM-2519-last-digit-collision-drift-and-k13-dirichlet-boundary
  - THM-2522-intrinsic-collision-depth-toothpick-descent-and-late-owner-decoupling
  - THM-2523-chi7-hamilton-energy-split-form-and-fano-boundary
related:
  - THM-2514-cyclic-k14-factor-chart-and-six-phase-ordinary-degree-reconstruction
  - THM-2521-k13-drift-k14-potential-module-bridge
script: 04-computation/lrc14_translated_chi7_polarization_thm2524.py
output: 05-knowledge/results/lrc14_translated_chi7_polarization_thm2524.out
script_sha256: e8c9ba91eda0d05ccdae17939b82982332c0cc210a67a9e5528e62a4e90c1b49
output_sha256: 8226cd9a9d2045a4cf2c82d74cc789c245d4de9b9c9e02478ad156ecbb027d4e
hash_basis: working-tree bytes (LF)
---

# THM-2524 -- the full translated `chi_7` polarization is lossless

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2523 proves that the multiplicative `chi_7` contrast of the six
Hamilton energies is a nondegenerate quadratic form, but its value on one
profile can vanish.  That is a diagonal-sampling failure, not a kernel of
the underlying operator.  Retain the thirteen cyclic polarizations instead:

```text
R_tau(t)=<p,A_tau T_t p>,                    t in F_13.        (1)
```

The entry `t=0` is the old scalar contrast.  The whole vector (1) is exactly
the image of the centred collision profile under the invertible THM-2523
operator.  Thus one fixed slope recovers every collision gap modulo the
irrelevant constant profile, even when its diagonal entry is zero.

This is a genuine scalar strengthening of the collision chart: for positive
rational drift all twelve nontrivial Fourier coefficients of (1) are
nonzero scalar integrals.  It is not phase retrieval.  The collision profile
is an autocorrelation, so the uncontracted predecessor vector and every
source/arrival orientation erased by self-correlation remain absent.

## 1. The centred collision profile

Use the notation of THM-2519.  At one collision depth let

```text
B=(B_u)_(u in F_13),

B_bar=1/13 sum_u B_u,

b_u=B_u-B_bar.                                                (2)
```

Then `sum_u b_u=0`.  Let

```text
p_r(z)=q_r(z)-1/13 sum_j q_j(z),              r in F_13       (3)
```

and regard `p=(p_r)` as a vector in the real Hilbert space
`L^2(G(z)dz)^13`.  For

```text
(T_t p)_r=p_(r+t),                                            (4)
```

the exact disintegration in THM-2519 gives

```text
b_t=1/13 <p,T_t p>,

D_13=b_0=1/13 ||p||^2.                                       (5)
```

In particular `b_(-t)=b_t`.  The constant part `B_bar` contains no
last-digit drift and is deliberately removed; it can be retained as one
extra scalar if reconstruction of the full `B` rather than `b` is desired.

Fix `tau in F_13^*`.  Put

```text
chi_7(s)=+1,       s in {1,2,4},
chi_7(s)=-1,       s in {3,5,6},                              (6)

(A_tau x)_r
 =-sum_(s=1)^6 chi_7(s)
      [x_(r+2tau s)+x_(r-2tau s)].                            (7)
```

This is THM-2523's symmetric signed circulant.  Its quadratic form is

```text
<p,A_tau p>
 =sum_(s=1)^6 chi_7(s)
    sum_h ||p_(h-tau s)-p_(h+tau s)||^2.                      (8)
```

## 2. The thirteen translated Hamilton polarizations

Define the full translated bank by

```text
R_tau(t)=<p,A_tau T_t p>.                                    (9)
```

Equivalently, polarize the six cycle gradients:

```text
R_tau(t)
 =sum_(s=1)^6 chi_7(s) sum_h
    <p_(h-tau s)-p_(h+tau s),
     p_(h+t-tau s)-p_(h+t+tau s)>.                           (10)
```

At `t=0`, equation (10) is exactly (8), the multiplicative Hamilton-energy
contrast of THM-2523.  For general `t`, expand (9) and use (5):

```text
R_tau(t)
 =-13 sum_(s=1)^6 chi_7(s)
       [b_(t+2tau s)+b_(t-2tau s)]

 =-13 sum_(s=1)^6 chi_7(s)
       [B_(t+2tau s)+B_(t-2tau s)]

 =13(A_tau b)_t.                                             (11)
```

The second equality uses `sum_s chi_7(s)=0`, so the constant `B_bar`
cancels.  Formula (11) is an exact finite identity.  The translate `t` is a
relative predecessor-root displacement, not a new physical time or a
future-owner delay.

## 3. Fourier diagonalization and explicit inversion

Let `zeta=zeta_13` and use normalized transforms

```text
f_hat(k)=1/13 sum_(t in F_13) f_t zeta^(-kt).                 (12)
```

The THM-2523 eigenvalue is

```text
lambda_(tau,k)
 =-sum_(s=1)^6 chi_7(s)
    [zeta^(2ktau s)+zeta^(-2ktau s)].                         (13)
```

It is nonzero for every `tau,k!=0`.  Fourier transform of (11) gives

```text
R_hat_tau(k)=13 lambda_(tau,k) B_hat(k),       k!=0.          (14)
```

Here `b_hat(k)=B_hat(k)` off the zero character.  At `k=0`, both sides of
(14) vanish because `A_tau` kills constants.

More is true than Fourier injectivity.  THM-2523 proves on the centred
space

```text
A_tau^6-39A_tau^4+299A_tau^2-325I=0.                         (15)
```

Consequently

```text
A_tau^(-1)
 =(A_tau^5-39A_tau^3+299A_tau)/325,                          (16)

b
 =1/4225
   (A_tau^5-39A_tau^3+299A_tau)R_tau.                        (17)
```

Thus every one of the six unoriented slope classes gives a lossless linear
transform of the twelve-dimensional centred collision profile.  In
particular the drift itself is recovered by the zeroth coordinate of (17):

```text
D_13
 =1/4225 [(A_tau^5-39A_tau^3+299A_tau)R_tau]_0.              (18)
```

The mean `B_bar` is the exact one-dimensional kernel on the uncentred
thirteen-space.  No stronger reconstruction of `B` without that scalar is
possible.

## 4. All root modes survive even when the diagonal vanishes

Assume the THM-2519 data are rational step functions and `D_13>0`.
THM-2519 gives

```text
B_hat(k)>0                         for every k!=0.             (19)
```

Combining (13), (14), and (19) yields

```text
R_hat_tau(k)!=0
       for every tau,k in F_13^*.                             (20)
```

This is a nonzero scalar-integral statement, not merely pointwise or
Hilbert-valued support.  Rationality can also be read directly: any nonzero
centred `b in Q^13` has all twelve primitive characters by irreducibility of
`Phi_13`, and the invertible rational operator `A_tau` preserves that law.

The diagonal can nevertheless vanish.  Take one point mass

```text
q_0=1,                 q_r=0 for r!=0.                        (21)
```

Then

```text
B_0=1/13,              B_u=0 for u!=0,

D_13=12/169.                                                 (22)
```

For every `tau!=0`,

```text
R_tau(0)=0,

{R_tau(t):t!=0}={six copies of -1, six copies of +1}.        (23)
```

Thus (23) lies on THM-2523's isotropic cone while its translated bank has
all twelve root modes.  The exact binary census in Section 8 contains
`9,828` such zero-diagonal slope/profile controls.  The repair is genuinely
the retained translation coordinate, not a claim that the old scalar was
nonzero.

## 5. Ordered cross-cospans preserve, but do not create, orientation

The algebra has a useful polarized form.  Let `p,q` be two supplied real
centred predecessor profiles on one common weighted ancestry space and put

```text
b^(p,q)_t=1/13 <p,T_t q>,

R_tau^(p,q)(t)=<p,A_tau T_t q>.                              (24)
```

The profile `b^(p,q)` is centred because `q` is centred.

The same expansion and inversion give

```text
R_tau^(p,q)=13A_tau b^(p,q),                                 (25)

b^(p,q)
 =1/4225
  (A_tau^5-39A_tau^3+299A_tau)R_tau^(p,q).                   (26)
```

Let `Jf(t)=f(-t)` and split `f=f^++f^-` into its even and odd parts.  The
kernel in (7) is even, so `A_tau J=J A_tau`.  Since `A_tau` is invertible on
centred vectors,

```text
(R_tau^(p,q))^-=13A_tau (b^(p,q))^-,

(R_tau^(p,q))^-=0       iff       (b^(p,q))^-=0.              (27)
```

For a self-correlation `p=q`, the left profile in (24) is even, and so is
the translated bank.  For an ordered cross-cospan,

```text
b^(q,p)=J b^(p,q),               R_tau^(q,p)=J R_tau^(p,q).  (28)
```

Therefore the bank retains an oriented odd component exactly when the
input cospan already has one.  It cannot manufacture source/arrival order
from THM-2519's self-product.  In complex Hilbert space (28) includes the
usual conjugation; the real LRC application above needs none.

## 6. Direct live depth-one application

THM-2522 proves that every one of the `165` live THM-2349 shallow
owner--word events has

```text
D_0!=0,                    m_*=0,                    L_*=1.   (29)
```

The events are rational Boolean step functions.  Hence their unweighted
depth-one collision profiles satisfy the hypotheses of Section 4.  For
every live event and every `tau!=0`, its bank (11) therefore has all twelve
nontrivial scalar root modes and reconstructs its centred collision profile
and positive drift by (17)--(18).

At `L=1`, THM-2522 also retains by label permutation both the shallow
septimal bank with `nu_13(c_j)=1` and the old deep bank with
`nu_13(c_3)>1`.  These are three simultaneous but differently typed
sidecars:

```text
translated collision bank R_tau:    lossless centred F_13 autocorrelation;
shallow septimal bank:              retained old source labels;
old deep bank:                      retained old c_3 labels.             (30)
```

No identity between the coordinates in (30) is asserted.  A sufficiently
late positive rational BV owner block can be added by THM-2522 without
moving the collision from `L=1`; its resulting positive weighted collision
profile admits the same translated inversion.  This remains a typed signed
packet, not a scalar-cover row or a Boolean emission theorem.

## 7. Exact gain and stopping boundary

The proved chain is

```text
positive rational collision drift
  -> rational nonconstant centred autocorrelation b
  -> one full translated chi_7 polarization R_tau=13A_tau b
  -> explicit recovery of b and D_13
  -> all twelve nontrivial scalar root modes of R_tau.        (31)
```

This closes the linear-algebraic information loss in taking only the
diagonal `chi_7` energy contrast.  It does not close the LRC seam.

1. `R_tau` is a signed linear combination of lawful collision
   intersections.  Neither (16) nor (17) is positivity-preserving or
   Boolean.
2. The recovered `b` is an autocorrelation.  It does not determine the
   phase of the uncontracted predecessor vector `p`; THM-2521's potential
   module still requires that vector as a sidecar.
3. The self-profile is even, `A_(-tau)=A_tau`, and no tournament/source
   orientation follows.
4. The slope chart and predecessor torsor are chosen affine gauges.  No
   theorem here identifies them with a runner, target, source, arrival,
   owner, old-deep, or future-deep label.
5. A cross-cospan can retain an existing odd part, but this theorem does not
   supply two semantically lawful ordered legs on one common ancestry sheet.
6. Nonzero scalar root modes of a typed collision packet are not yet an
   emitted Boolean owner-loop current or a scalar-cover row exclusion.

Accordingly LRC(14) remains open.  The next bridge must realize or pair the
signed translated bank on the live owner/source/deep ancestry object without
discarding its root phase.

## 8. Exact dependency-free referee

Run

```bash
python3 04-computation/lrc14_translated_chi7_polarization_thm2524.py
python3 -O 04-computation/lrc14_translated_chi7_polarization_thm2524.py
```

Both executions reproduce

```text
05-knowledge/results/lrc14_translated_chi7_polarization_thm2524.out
```

byte-for-byte.  The companion verifies:

- for all twelve nonzero slopes, the integer matrix identity

  ```text
  A^6-39A^4+299A^2-325I=-25J,
  ```

  and both sides of the inverse identity on all `1,872` centred-basis
  entries;
- all `144` nonzero slope/frequency multipliers and `1,728` basiswise
  Fourier factorizations in `F_547`;
- all `8,190` nonconstant Boolean profiles, all `49,140` unoriented
  slope/profile pairs, and all `589,680` rational primitive modes using the
  exact `Phi_13` remainder test;
- the `9,828` positive-drift controls with `R_tau(0)=0`;
- the centred-delta values (22)--(23); and
- direct gradient polarization, reversal, a nonzero odd cross component,
  and the even self-correlation boundary.

Normal and optimized executions reproduce the stored transcript exactly.
The finite-field checks audit the intertwining, while nonvanishing over the
cyclotomic field is proved by THM-2523 and rational irreducibility, not
inferred from reduction modulo one prime.  An independent root audit
rederived (5), the translated polarization identity (11), the normalized
factor `13` in (14), the polynomial inverse and denominator `4,225`, the
centred-delta hostile, the cross-cospan parity equivalence, every finite
count, and all Boolean/phase/orientation stopping claims; it also reproduced
normal, optimized, and stored output byte-for-byte. **QED.**
