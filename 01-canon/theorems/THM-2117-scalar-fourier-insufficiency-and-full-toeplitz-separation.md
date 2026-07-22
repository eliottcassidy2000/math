---
id: THM-2117
title: "Scalar Fourier insufficiency and full Toeplitz separation on the half-fiber"
status: >
  PROVED. There is an explicit rank-eight binary-shifted half-fiber family
  which passes every THM-2105 clock m=2,...,7, has maximum Hunter-tree mass
  strictly below 1/7, and satisfies |hat H(r)|<1/7 at every nonzero scalar
  frequency, but nevertheless has an explicit open safe interval. Thus all
  2-by-2 Toeplitz minors can pass while a larger finite Toeplitz section must
  fail. The scalar all-frequency claim is reduced to an exact 181-state LCM
  closure census; the safe phase and pair tree are exact rational audits.
  This separates proof carriers and does not prove LRC(14).
source: codex-2026-07-22-LRC-scalar-toeplitz-separation
depends_on:
  - THM-641
  - THM-2115
related:
  - THM-1185
  - THM-2105
  - THM-2114
  - THM-2121
  - HYP-2974
---

# THM-2117 -- scalar Fourier insufficiency and full Toeplitz separation

On the guard half-fiber of THM-2115, write a terminal character as

```text
c_i=a_i p+n_i q,
e_i=a_i mod 2 in {0,1},        v_i=|n_i|,
B_i(beta)=1_{||v_i beta+e_i/2||<1/14},
N(beta)=sum_i B_i(beta),       H(beta)=N(beta)-1.       (1)
```

Consider the following eight `(v_i,e_i)` pairs:

```text
(840,0), (4228,1), (4984,1), (5691,1),
(7070,1),(4256,1), (6195,1), (4445,1).                (2)
```

We prove four simultaneous facts:

1. every affine clock `m=2,...,7` from THM-2105 is covered;
2. the maximum half-fiber Hunter tree has mass strictly less than `1/7`;
3. every nonzero scalar coefficient obeys `|hat H(r)|<hat H(0)=1/7`;
4. `H=-1` on an explicit open interval, so some finite Toeplitz quadratic is
   negative.

Consequently the clock, pair-tree, and complete scalar-frequency gates do not
imply full Toeplitz positive semidefiniteness.

## 1. Every small clock is blind

The clock at `m` asks whether

```text
Z/(2m)=union_i {j:v_i j+m e_i=0 mod 2m}.               (3)
```

For every `2<=m<=7`, one has `2m|840`. Since the first shift in (2) is zero,
the first terminal alone contains every residue in (3). This deliberately
hostile device makes all six clocks pass for the same reason; no accidental
union of partial carriers is involved.

## 2. The exact pair-tree test also passes

Let `w_ij=measure(B_i intersect B_j)`. For two labels `(v,e),(w,f)`, reduce
`q=v/g,r=w/g`, where `g=gcd(v,w)`, and set

```text
a=e/2-1/14,       b=f/2-1/14,       phi=q b-r a.       (4)
```

THM-641, with both widths equal to `1/7`, gives

```text
w((v,e),(w,f))=1/49+
 [-B2(phi+q/7)+B2(phi)+B2(phi+(q-r)/7)-B2(phi-r/7)]/(2qr),   (5)
```

where `B2(x)={x}^2-{x}+1/6`. Exact Kruskal maximization over the 28 values in
(5) gives the tree

```text
edge       mass
(4,6)      2554/125139
(1,4)      10894/533785
(4,5)      5483/268660
(3,6)      11421/559615
(2,7)      6459/316484
(2,6)      15003/735140
(0,3)      387/18970,                                  (6)
```

with total

```text
tau=3142001236410829/21994608830427660
   =1/7-85739364551/21994608830427660 <1/7.            (7)
```

Hunter says `measure(union_i B_i)<=8/7-tau`. A contradiction to coverage
would require `tau>1/7`; hence this strongest pair-tree invoice passes, by a
small but strict rational margin.

## 3. Every scalar Toeplitz minor passes

THM-2115 gives, for `r!=0`,

```text
hat H(r)=sum_(i:v_i|r)
 (-1)^(e_i r/v_i) sin(pi r/(7v_i))/(pi r/v_i),
hat H(0)=1/7.                                           (8)
```

The `2 by 2` principal Toeplitz minor on frequencies `{0,r}` is PSD exactly
when `|hat H(r)|<=1/7`. We verify this for every integer `r`, not merely over a
finite scan.

If no `v_i` divides `r`, (8) is zero. If exactly one divides it, put
`|r|=k v_i`. The elementary induction

```text
|sin(kx)|<=k sin(x),       0<x<pi, k in N,              (9)
```

at `x=pi/7`, followed by `sin(x)<x`, gives

```text
|hat H(r)|<=sin(pi/7)/pi<1/7.                           (10)
```

It remains to control simultaneous divisors. For a frequency `r`, let

```text
D(r)={v_i:v_i|r},       L=lcm(D(r)).                    (11)
```

If `|D(r)|>=2`, then `L|r`, and the set of speeds dividing `L` is exactly
`D(r)`. Therefore

```text
sum_(v in D(r)) v/|r| <= sum_(v in D(r)) v/L.           (12)
```

Every possible right side arises by taking a nonempty subset of the eight
speeds, forming its LCM, and closing under divisibility. There are only 181
distinct `(LCM,closure)` records with closure size at least two. Exact
enumeration gives

```text
max sum_(v|L) v/L =67/472<1/7,                          (13)
```

uniquely at the maximal record used here

```text
L=49560,       D={840,6195}.                            (14)
```

Using `|sin y|<=1` termwise in (8), `pi>1`, and (12)--(13),

```text
|hat H(r)|<=sum_(v_i|r) v_i/(pi|r|)<1/7.               (15)
```

Thus every scalar inequality, equivalently every `2 by 2` Toeplitz principal
minor, is strictly positive.

## 4. An exact open safe interval

Take

```text
beta_0=223163/6285230.                                  (16)
```

The exact circle distances and their margins over `1/14` are

```text
(v,e)       ||v beta_0+e/2||       margin
(840,0)       15714/89789           18601/179578
(4228,1)      341993/897890         19847/64135
(4984,1)      414471/897890         25024/64135
(5691,1)      195603/448945         327071/897890
(7070,1)      841/1778              51/127
(4256,1)      347231/897890         141548/448945
(6195,1)      41240/89789           69653/179578
(4445,1)      229/707               51/202.              (17)
```

All margins are positive. Direct rational endpoint union goes further: the
largest safe-gap length is `51/448945`, and one such component is

```text
I=(1103/31115,176/4949),       beta_0=midpoint(I).       (18)
```

There are 12,824 safe components in all, with total safe mass

```text
1717225315810035/5865229021447376>0.                    (19)
```

In particular `H=-1` throughout `I`.

## 5. Why a finite higher Toeplitz section must fail

For a finite coefficient vector `(c_j)`, the Toeplitz quadratic is

```text
sum_(j,k) conjugate(c_j)c_k hat H(j-k)
=integral_T H(beta)|sum_j c_j exp(2 pi i j beta)|^2 d beta.   (20)
```

Center the normalized Fejer kernels at `beta_0`. Each is the squared modulus
of a finite trigonometric polynomial, so its integral against `H` is a
quadratic of the form (20). Since `H` is constant and equal to `-1` on a
neighborhood of `beta_0`, the approximate-identity theorem gives

```text
integral_T H(beta) F_K(beta-beta_0) d beta -> -1.        (21)
```

Hence (20) is negative for some finite `K`. This proves that a larger finite
Toeplitz section fails even though every one of its `2 by 2` principal minors
passes. No numerical eigenvalue or floating-point sign is used. THM-2121 later
makes the section size effective: the displayed cell is detected by every
integer `K>8/(51/448945)`, in particular by `K=70423`.

## 6. What this changes in the LRC route

The example kills the proposed implication

```text
small clocks + Hunter tree + all scalar Fourier bounds
   => half-fiber cover.                                 (22)
```

The surviving harmonic carrier must retain joint frequency phase through a
Toeplitz minor, a Fejer polynomial, or an equivalent positive trigonometric
certificate. A scalar packet bank cannot do this. Conversely, THM-1185's
symbol argument warns that an unrestricted full-Toeplitz search merely
re-encodes open-set coverage and cannot see boundary-only lonely points. The
useful target is therefore a **bounded or structurally forced** negative
minor after the finite clocks/content/toothpick gates, with endpoint owners
kept as a sidecar.

The family (2) is not an LRC(14) counterexample. It is a hostile half-fiber
row proving strict separation among necessary cover tests. The binary shifts
`0,1/2` are forced by character parity; no arbitrary-shift conjecture is used.

## 7. Assumption challenge and Tournament Analysis

The challenged assumption was that scalar frequency packets exhaust the new
information after pairwise overlap. They do not: positivity couples several
Fourier differences inside one quadratic form.

Candidate tournament vertices were runners, gaps, pair edges, clock residues,
LCM closures, frequencies, Toeplitz minors, and proof obligations. The faithful
finite objects are the divisibility hypergraph of LCM-closed packets and the
Toeplitz moment matrix. Orienting packets by their scalar envelope, with
increasing LCM as the tie Hamiltonian path, preserves search priority but
destroys joint phase and therefore PSD. Such a tournament can schedule the 181
records; it cannot be the certificate.

## 8. Exact verification

The companion script uses `Fraction` throughout the clocks, THM-641 pair
matrix, Kruskal tree, LCM closure census, phase margins, and interval union.
It uses floating point only to print the illustrative value
`sin(pi/7)/pi`; inequality (10) is analytic. Runtime checks use `require`, so
normal and `python -O` executions are byte-identical and end in `PASS`.

```text
04-computation/lrc14_signed_toeplitz_scalar_separation_codex_20260721.py
05-knowledge/results/lrc14_signed_toeplitz_scalar_separation_codex_20260721.out
source SHA-256 ade4e9f14f37714adcdb88a6137b272a7de9f4ad75dc55f604592266ff439110
output SHA-256 9516fd53a1b50d8bc17b56f9f51d0712e28dc21cd03f04a85174edbcd60737f3
```

The hash basis is the working-tree files with LF line endings. QED.
