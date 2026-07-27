---
id: THM-2510
title: "Truncated-Radon tight frame and integral dipole boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. For every row-zero
  array on F_p x {0,...,q-1}, 2<=q<p, the complete nonzero-slope
  truncated-Radon bank is a tight frame on the horizontal-mean-zero
  quotient: its total counting energy is exactly p times the squared
  distance from the h-independent kernel. The kernel has dimension q-1.
  For a nonvertical integral array the sharp energy floor is 2(p-1),
  with equality exactly at one horizontal unit dipole plus an arbitrary
  vertical kernel vector. At (p,q)=(13,7) the floor is 24; the 42-cut
  THM-2508 bank has energy 42 times the base bank and sharp floor 1008.
  THM-2436 essential parents therefore have pointwise energy at least 24
  and inherited integrated floors 48/7 and 72/7. Energy is quadratic and
  phase-forgetting; no signed physical current, live-row transplant, or
  LRC(14) closure follows.
source: codex-2026-07-27-truncated-radon-tight-frame
depends_on:
  - THM-2435-top-blocker-essential-parent-and-punctured-stalk-carrier
  - THM-2436-punctured-ninety-one-stalk-repeated-step-spectrum
  - THM-2507-truncated-radon-toothpick-tomography-and-nonaffine-root-boundary
  - THM-2508-affine-cut-bundle-covariance-and-carry-permutation
related:
  - THM-2506-punctured-stalk-primitive-module-saturation-and-thirteen-primary-pushforward-no-go
  - THM-2509-antipodal-radon-cospan-and-lossless-septimal-chart
  - THM-2511-affine-cut-quadratic-root-service-and-pair-space-boundary
script: 04-computation/lrc14_truncated_radon_tight_frame_thm2510.py
output: 05-knowledge/results/lrc14_truncated_radon_tight_frame_thm2510.out
script_sha256: 8e5bcfe9b23f00648cb4d55ccd39975ff88a1d7752c848c31431e7f1731225fe
output_sha256: c2064f4cc59a6e585c3c514bd009aef014053097473353d880798b176c9c9061
hash_basis: working-tree bytes (LF)
---

# THM-2510 -- the toothpick bank is a tight frame

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2507 determines the rank and worst-case number of surviving toothpick
slopes.  The complete bank has a stronger metric structure: after removing
the exact vertical kernel, it is a tight frame.  There is no hidden condition
number and no direction-dependent loss:

```text
total Radon energy
  = prime times horizontal fluctuation energy.                    (1)
```

For integral defects this identity has a sharp first quantum.  It is one
unit dipole between two septimal columns at one horizontal site.

## 1. General truncated-Radon identity

Let `p` be prime and

```text
2<=q<p.
```

Write `I_q={0,...,q-1} subset F_p`.  On the complex vector space

```text
V_(p,q)
 ={d:F_p x I_q -> C:
     sum_(r in I_q)d(h,r)=0 for every h},                        (2)
```

use counting-measure inner products and define

```text
R_tau d(v)=sum_(r in I_q)d(v-tau r,r),
                  tau in F_p^*.                                 (3)
```

Let `E_h` be horizontal averaging within every column:

```text
(E_h d)(h,r)=dbar_r=(1/p)sum_(x in F_p)d(x,r).                  (4)
```

Then

```text
sum_(tau in F_p^*) ||R_tau d||_2^2
 =p||d-E_h d||_2^2                                             (5)

 =sum_(r in I_q)
   [p sum_h |d(h,r)|^2-|sum_h d(h,r)|^2].                       (6)
```

By polarization, if

```text
T:d -> (R_tau d)_(tau!=0),                                     (7)
```

then as an operator on `V_(p,q)`,

```text
T^*T=p(I-E_h).                                                  (8)
```

Thus the full nonzero-slope bank is a tight frame with frame constant `p`
on the orthogonal complement of its kernel.

### Proof

Fix a primitive `p`-th root `zeta`.  For the horizontal Fourier transform

```text
D_alpha(r)=sum_h d(h,r)zeta^(-alpha h),
P_alpha(X)=sum_(r=0)^(q-1)D_alpha(r)X^r,                        (9)
```

THM-2507's slice calculation gives

```text
(R_tau d)_hat(alpha)=P_alpha(zeta^(-alpha tau)).                (10)
```

The row-zero law says `P_alpha(1)=0` for every `alpha`.  If `alpha!=0`,
the numbers `zeta^(-alpha tau)`, `tau!=0`, are all nontrivial `p`-th roots.
Root orthogonality, with the value at `1` restored and then removed, gives

```text
sum_(tau!=0)|P_alpha(zeta^(-alpha tau))|^2
 =p sum_r |D_alpha(r)|^2-|P_alpha(1)|^2
 =p sum_r |D_alpha(r)|^2.                                      (11)
```

The `alpha=0` output is zero because the total mass is zero.  Apply
unnormalized Parseval first to every output and then to every column:

```text
sum_tau ||R_tau d||_2^2
 =sum_(alpha!=0,r)|D_alpha(r)|^2
 =p sum_(h,r)|d(h,r)|^2-sum_r|D_0(r)|^2.                       (12)
```

Since `D_0(r)=sum_h d(h,r)`, this is (5)--(6).  Polarization proves (8).
QED.

The same identity has a direct pair-count form.  For each column `x_h`,

```text
p sum_h |x_h|^2-|sum_h x_h|^2
 =sum_(h<h')|x_h-x_(h')|^2.                                   (13)
```

Equation (13) makes both the kernel and the integral boundary transparent.

## 2. Exact kernel and normalization

Every term on the right of (13) vanishes exactly when the column is constant
in `h`.  Hence

```text
ker T
 ={d(h,r)=b_r:sum_r b_r=0},

dim ker T=q-1.                                                  (14)
```

This is exactly the vertical kernel already identified by the rank theorem
in THM-2507.  Equation (8) adds the metric statement: all nonzero singular
values on the quotient are `sqrt(p)`.

If both spaces use probability norms, (5) becomes

```text
sum_(tau!=0)||R_tau d||_(F_p,prob)^2
 =p q ||d-E_h d||_(F_p x I_q,prob)^2.                          (15)
```

At `(p,q)=(13,7)`, the coefficient in (15) is `91`.  No normalization is
implicit in (5)--(6); those equations use literal counting sums.

## 3. Sharp integral energy quantum

Suppose `d` is integral and not in the kernel (14).  Then

```text
sum_(tau!=0)||R_tau d||_2^2>=2(p-1).                           (16)
```

This bound is sharp.  Equality holds exactly for

```text
d(h,r)=b_r+epsilon 1_(h=h_*)
                    (1_(r=r_0)-1_(r=r_1)),                    (17)

r_0!=r_1,       epsilon in {+1,-1},       b_r in Z,
sum_r b_r=0.
```

To prove the floor, let `x in Z^p` be a nonconstant column.  If `m` entries
attain its minimum, the `m(p-m)` pairs crossing from a minimum entry to a
larger entry contribute at least one each to (13).  Therefore

```text
p sum_h x_h^2-(sum_h x_h)^2>=m(p-m)>=p-1.                      (18)
```

Equality in (18) requires `m in {1,p-1}` and all other entries to be the
adjacent integer: exactly one entry differs from the other `p-1` by one.

A pointwise row-zero array cannot have exactly one nonconstant column,
because the sum of all its constant columns would force the last column to
be constant.  Thus at least two columns contribute (18), proving (16).

For equality, precisely two columns vary, each by one exceptional unit.  At
some horizontal site outside their at most two exceptional sites, the
constant baselines already sum to zero.  The two exceptional terms must then
cancel pointwise, forcing a common site and opposite signs.  This is exactly
(17).  Conversely substitution of (17) into (13) gives `2(p-1)`.

For arrays whose denominators divide `M`, scaling gives the sharp floor

```text
2(p-1)/M^2.                                                     (19)
```

There is no positive scale-free floor over unrestricted rational arrays.

## 4. The `13 x 7` punctured stalk

For `(p,q)=(13,7)`, equations (5) and (16) read

```text
sum_(tau=1)^12 ||R_tau d||_2^2
 =13||d-E_h d||_2^2>=24                                      (20)
```

for every nonvertical integral `d`.  THM-2436 proves that its integral
punctured-stalk defect is vertical exactly when it is zero, and that zero is
exactly the flat locus.  Hence every essential parent has pointwise Radon
energy at least `24`.

Let `mathcal E` be the flat locus in the parent `P`.  THM-2435 gives

```text
mu(P minus mathcal E)>=(1+k)/7.                                (21)
```

Integrating (20) yields the exact inherited floors

```text
k=1: integral_P sum_tau ||R_(Y,tau)||_2^2 dY >=48/7,

k=2: integral_P sum_tau ||R_(Y,tau)||_2^2 dY >=72/7.           (22)
```

These bounds do not pay THM-2507's slope-selection fraction.  They retain
all slopes quadratically.  They also remain entirely inside THM-2436's
already-empty high-septimal branch.

## 5. The complete affine-cut bank

Use THM-2508's `42` ordered affine cuts

```text
R_(tau,a,c)d(v)
 =sum_(r in F_7)d(v-tau rep(ar+c),r),

a in F_7^*,       c in F_7.                                   (23)
```

For fixed `(a,c)`, reindexing `s=ar+c` only permutes the seven columns of
`d`.  Both `||d-E_h d||_2` and integrality are unchanged.  Applying (20) to
each cut and summing gives

```text
sum_(a,c,tau)||R_(tau,a,c)d||_2^2
 =42*13||d-E_h d||_2^2.                                       (24)
```

For a nonvertical integral array the sharp full-bank floor is therefore

```text
42*24=1008.                                                     (25)
```

The THM-2506 two-row hostile has

```text
13||d-E_h d||_2^2=46,
```

so (24) gives

```text
42*46=1932,                                                     (26)
```

exactly explaining the quadratic energy printed by THM-2508's independent
cut-bundle probe.

Equation (24) is a total-bank identity.  It does not prejudge the more
refined pointwise pair-space observable reserved under THM-2511.

## 6. What energy preserves and destroys

The tight frame proves that the complete Radon bank loses no **quadratic
size** beyond the known vertical kernel.  It deliberately sums over slopes,
cuts, output sites, signs, and phases.  Consequently:

- it is invariant under the antipodal leg swap of THM-2509;
- it cannot orient a tournament on slopes;
- it does not choose one fixed signed marginal;
- it does not defeat cancellation after integration over a moving parent;
- it does not identify the cut character with a THM-2471 collision colour;
  and
- it does not provide an owner, terminal word, deep sheet, or old address.

THM-2508 proves the corresponding linear boundary exactly: the zero cut
character and the only affine-invariant linear scalar both vanish.  Passing
to the positive energy in (24) restores a gauge-invariant number by forgetting
the phase that a physical current needs.  Thus the energy theorem is a sharp
invoice and equality classifier, not a transplant from the closed
high-septimal branch to the `165` live rows.

## 7. Exact companion and independent audit

Run

```bash
python3 04-computation/lrc14_truncated_radon_tight_frame_thm2510.py
python3 -O 04-computation/lrc14_truncated_radon_tight_frame_thm2510.py
```

Both executions reproduce the stored transcript byte-for-byte.  The
dependency-free referee checks:

- the complete `78`-dimensional row-zero chart and all `3,081` entries of
  the symmetric operator Gram matrix;
- operator rank `72` and the six explicit vertical kernel vectors;
- all `728` nonzero one-row controls with free entries in `{-1,0,1}`, whose
  minimum is `24` and whose `42` equality cases are exactly unit dipoles;
- all `546` unit dipoles after adding a nonzero vertical baseline;
- `200` deterministic integral arrays and `20` direct full-`42`-cut checks;
  and
- the THM-2506 energies `46` and `1932`.

The independent audit separately derived (5) by root orthogonality and by
physical pair counting, checked the probability normalization, kernel,
integral minimum, complete equality classification, cut multiplier, and
THM-2436 mass invoices.  It explicitly found no linear or live-row
consequence beyond the stated quadratic boundary. QED.
