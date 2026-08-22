---
id: THM-3665
title: "LRC support-minimal three-twist target detector"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT; PENDING INDEPENDENT HOSTILE AUDIT.
  A nonzero mean-zero convolution mask on F_p^d whose Fourier transform is
  nonzero at every nontrivial character needs at least d+1 support sites,
  and an explicit d+1-site mask attains the bound.  For THM-2334's
  two-dimensional twist group, the optimal mask is delta_0+delta_a-2delta_b
  for any ordered basis (a,b).  Thus nonzero target survival is equivalent
  to one three-twist imbalance H(s)+H(s-a)-2H(s-b) being nonzero.  Exact
  spectral extrema give sharp l2 frame constants.  This does not prove that
  a covering-row imbalance is nonzero or identify twists with physical
  ancestry addresses.
source: kps-s192 / support-capacity sharpening of THM-3664, 2026-08-21
depends_on:
  - THM-2334-relation-residue-current-and-character-twist-pushforward
related:
  - THM-3661-lrc-exceptional-detector-simple-spectrum-convolution-rigidity
  - THM-3664-lrc-sparse-eight-twist-target-detector-frame
script: 04-computation/lrc_minimal_three_twist_target_detector_thm3665.py
output: 05-knowledge/results/lrc_minimal_three_twist_target_detector_thm3665.out
script_sha256: a5ef75f038d80c5d91308bb5379303970b44e9e538323eb49cf8779386356938
output_sha256: 172fe3e32fc27bb2abb21f4c7a7af59e71cdfa4c604586b2d4f8725a5ac6211a
semantic_sha256: 7130760feb9c252c6f7d1b2dc5fdd2b224f3469718aa4363380511de88c586d3
hash_basis: raw LF bytes
---

# THM-3665 -- three twists are necessary and sufficient

**PROVED + FINITE-EXACT + VERIFIED-EXACT; PENDING INDEPENDENT HOSTILE
AUDIT.**  The eight-point detector of THM-3664 is useful because it descends
from the exceptional two-current geometry.  Once the task is posed purely on
the freely based target-twist group, however, a universal rank-plus-one law
gives the optimal detector.

## 1. The rank-plus-one support law

Let

```text
V=F_p^d, d>=1,                                      (1)
```

and let `f:V->C` be nonzero with mean zero.  Suppose

```text
f_hat(q)!=0 for every q!=0.                         (2)
```

Then

```text
|supp(f)|>=d+1.                                     (3)
```

Indeed, if the support has `r<=d` points, choose one point `x_0`.  The
differences

```text
{x-x_0:x in supp(f)}                               (4)
```

span a subspace of dimension at most `r-1<d`.  A nonzero character `q`
annihilates that span, so its phase is constant on the support.  Therefore

```text
f_hat(q)=chi_q(-x_0) sum_x f(x)=0,                  (5)
```

contradicting (2).

The bound is attained for every prime and dimension.  Given a basis
`e_1,...,e_d`, put

```text
m_d=delta_0+sum_(i=1)^(d-1)delta_(e_i)-d delta_(e_d).
                                                               (6)
```

It has `d+1` support sites and mean zero.  If `m_d_hat(q)=0`, then a sum of
`d` unit complex numbers equals `d` times another unit complex number.
Equality in the triangle inequality forces all `d` summands to have the
same phase.  One summand is the phase at zero, hence all phases are one.
The character annihilates the whole basis and `q=0`.  Thus (6) satisfies
(2), proving

```text
minimum support=d+1.                                (7)
```

This is a support theorem for linear, translation-covariant, mean-zero
detectors.  It does not claim optimality among arbitrary nonlinear tests.

## 2. The optimal THM-2334 detector

THM-2334's twist group is a two-dimensional vector space

```text
T=G^ isomorphic to F13^2.                           (8)
```

Choose any ordered basis `(a,b)` of `T`.  Specializing (6) gives

```text
m_(a,b)=delta_0+delta_a-2delta_b.                   (9)
```

For the THM-2334 boundary twist profile `H:T->C`, define

```text
M_(a,b)(s)=(m_(a,b)*H)(s)
           =H(s)+H(s-a)-2H(s-b).                   (10)
```

With the same Fourier convention as THM-2334,

```text
M_hat(q)=m_hat(q)H_hat(q)=169m_hat(q)A(q).          (11)
```

The rank-plus-one law says `m_hat(q)` vanishes exactly at `q=0`.  Hence the
following are equivalent:

1. `A(q)!=0` for some nonzero target vector;
2. the twist profile `H` is nonconstant;
3. `M_(a,b)` is not identically zero;
4. `H(s)+H(s-a)-2H(s-b)!=0` for some centre `s`.

Thus three distinct twist values are sufficient, and (7) proves that no
two-site mean-zero convolution detector can have the same iff property.
More concretely, a two-site mask is a scalar translate of
`delta_0-delta_c`; all twelve nontrivial characters in `c^perp` kill it.

As before, the complete nonconstant profile is reconstructed by the inverse
multiplier

```text
h_hat(0)=0,
h_hat(q)=1/m_hat(q), q!=0,
h*M_(a,b)=H-H_bar.                                  (12)
```

## 3. Sharp spectrum at level thirteen

Use coordinates `a=(1,0)`, `b=(0,1)` and

```text
zeta=exp(2*pi*i/13).
```

The multiplier is

```text
m_hat(u,v)=1+zeta^(-u)-2zeta^(-v).                  (13)
```

Choose the representative `u` in `[-6,6]`, set

```text
c=cos(pi*u/13),
delta=pi*(2v-u)/13 mod 2*pi, |delta|<=pi.           (14)
```

Then

```text
|m_hat(u,v)|=2|c-exp(-i delta)|.                    (15)
```

For the minimum, if `delta=0`, parity makes `u` even.  Outside the trivial
frequency, `|u|>=2`, and

```text
2(1-c)>=2(1-cos(2*pi/13))=4sin(pi/13)^2.            (16)
```

If `0<|delta|<pi`, its imaginary part gives

```text
2|c-exp(-i delta)|>=2sin(pi/13)>4sin(pi/13)^2.      (17)
```

The strict inequality follows from `sin(pi/13)<1/2`.  At `|delta|=pi`,
the left side is at least two.  Equality in (16) occurs exactly at

```text
(u,v)=(2,1),(11,12).                                (18)
```

Therefore

```text
min_(q!=0)|m_hat(q)|=4sin(pi/13)^2.                 (19)
```

For the maximum, first take `|delta|<pi`.  If `cos(delta)<=0`, the distance
in (15) increases with `c in [0,1]`, and its maximum is bounded by

```text
2|1-exp(-i delta)|<=4cos(pi/26).                    (20)
```

If `cos(delta)>0`, convexity in `c` puts the maximum at `c=0` or `c=1`,
and the same bound holds.  If `|delta|=pi`, parity forces `u` odd, so

```text
2(1+c)<=2(1+cos(pi/13))
       =4cos(pi/26)^2<4cos(pi/26).                  (21)
```

Equality in (20) occurs exactly at

```text
(u,v)=(0,6),(0,7).                                  (22)
```

Consequently

```text
max_q |m_hat(q)|=4cos(pi/26).                       (23)
```

## 4. The sharp frame inequality

Parseval, (19), and (23) yield

```text
16sin(pi/13)^4 sum_s|H(s)-H_bar|^2
 <=sum_s|M_(a,b)(s)|^2
 <=16cos(pi/26)^2 sum_s|H(s)-H_bar|^2.              (24)
```

Both constants are attained by character profiles at the frequencies in
(18) and (22), so (24) is sharp.  Compared with the inherited eight-point
frame, the three-site mask improves both support and conditioning.

## 5. Exact cyclotomic ledger

For projective representatives `(1,t)`, `0<=t<=12`, followed by `(0,1)`,
the fourteen exact norms of (13) are

```text
(13,13,35503,20969,18603,6773,26039,
 169,26039,6773,18603,20969,35503,53248).           (25)
```

Their product, the characteristic-zero augmentation determinant, is

```text
9072750758202804116713220630638041462433831593234432
 =13*26417870930056774319865792^2.                  (26)
```

Modulo the pinned exact field used by the companion this is

```text
408889968250631621108841.                           (27)
```

The norm ledger is not needed for the triangle-inequality proof, but it is
an independent exact nonvanishing certificate and a reproducible arithmetic
fingerprint.

## 6. Basis freedom and the new concrete target

There are exactly

```text
|GL(2,13)|=26208                                    (28)
```

ordered bases `(a,b)`, and all give distinct based masks (9).  The decisive
covering-row target is therefore:

```text
choose a basis adapted to the interval-coordinate action and prove
H(s)+H(s-a)-2H(s-b)!=0 for one s.                   (29)
```

Equation (10) has a useful geometric reading: at one centre, the value on
one basis edge cannot equal twice the value on the other edge minus the
centre.  Possible proof mechanisms include a unique interval jump, strict
convexity along one translated coordinate, unequal first nonzero boundary
jets, or a valuation that singles out one of the three terms.

The theorem does not establish (29) for any hypothetical covering row,
does not apply directly to the all-`91`-unit projector, and does not align
`a,b` with the two-current ancestry digits of THM-3657.  LRC(14) remains
open.

## 7. Exact companion

Reproduce with

```bash
python3 -B 04-computation/lrc_minimal_three_twist_target_detector_thm3665.py
python3 -B -O 04-computation/lrc_minimal_three_twist_target_detector_thm3665.py
```

Both streams match the stored transcript.  The assertion-free companion
checks all `169^2` transfer identities, the complete 168-difference
two-site obstruction census, the fourteen cyclotomic norms and determinant,
the exact extremizer case ledger, and all `26208` ordered bases.  **QED.**
