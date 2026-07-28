---
id: THM-2645
title: "Eleven-sheet multiplicity full-character spectrum and energy split"
status: >
  PROVED + VERIFIED-EXACT + TWO INDEPENDENT HOSTILE AUDITS.  Let p be
  odd, let S=F_p\A and T=F_p\B with A,B two-point sets, and let
  r=1_S*1_T be the two-edge representation multiplicity.  Every nonzero
  normalized Fourier coefficient satisfies
  rhat(k)=p Ahat(k)Bhat(k)!=0, with the sharp uniform amplitude floor
  4 sin^2(pi/(2p))/p.  Its centered energy is
  2(3p-8)/p^2 when the two complement steps agree up to sign and
  4(p-4)/p^2 otherwise.  At p=13 all twelve charged carry characters
  survive; the energies are 62/169 and 36/169 over exactly 1,014 and 5,070
  ordered relation pairs, so some charged mode has square at least 3/169.
  Thus THM-2642's support saturation loses carry selection but its exact
  multiplicity retains full charged spectrum.  The profile remains blind to
  a free thirteen-fold common-origin gauge and has no LRC consequence absent
  a physical same-base positive transition-count table.
source: deep-energy-audit-2026-07-28-multiplicity-character
depends_on:
  - THM-2642-cyclic-difference-relation-saturation-and-thick-holotopy-no-go
related:
  - THM-2380-cross-word-charged-target-correlation-and-pair-twist-gate
  - THM-2634-endpoint-pair-two-carry-cospan-and-single-carry-no-go
  - THM-2637-derangement-character-fixed-branch-holotopy-principle
script: 04-computation/lrc14_eleven_sheet_multiplicity_spectrum_thm2645.py
output: 05-knowledge/results/lrc14_eleven_sheet_multiplicity_spectrum_thm2645.out
script_sha256: c3cd514ab3a6ce392d07bf33ce802709e038ea182e5cb34e0e6828bc2f092772
output_sha256: 569a140ab0b731e526cdff3148577d79278cae5e9bea7e83c84c16e25b2a2554
hash_basis: LF-normalized bytes
---

# THM-2645 -- multiplicity restores every charged carry colour

**PROVED + VERIFIED-EXACT + TWO INDEPENDENT HOSTILE AUDITS.**

THM-2642 proves that two eleven-sheet relations on `C_13` have positive
sections for every affine clutch.  Support has therefore become useless for
selecting a carry.  The exact multiplicity profile behaves very differently:
its small two-point complement forces every charged Fourier colour to
survive.  This is the cheapest positive sidecar left by the saturation
no-go.

## 1. Complement factorization

Let `p` be an odd prime, `G=F_p`, and choose two-point sets

```text
A={a_0,a_1},                 B={b_0,b_1}.                (1)
```

Put

```text
S=G\A,                       T=G\B,                      (2)
```

and let

```text
r(c)=(1_S*1_T)(c)
    =#{(s,t) in S x T:s+t=c}.                            (3)
```

As in THM-2642, this is the number of two-edge increment decompositions per
base sheet; the total number of `c`-clutched sections is `p r(c)`.

In the integral group algebra,

```text
1_S=1_G-1_A,                  1_T=1_G-1_B.               (4)
```

For a fixed `c`, inclusion--exclusion gives

```text
r(c)=p-4+m(c),                m=1_A*1_B.                 (5)
```

At `p=13`, this is the exact THM-2642 formula `r=9+m`.

## 2. Every charged character survives

Fix a primitive `p`th root `zeta` and use the normalized transform

```text
fhat(k)=p^(-1) sum_(x in G) f(x) zeta^(-kx).             (6)
```

For unnormalized convolution,

```text
(f*g)^hat(k)=p fhat(k)ghat(k).                           (7)
```

When `k!=0`, `1_G^hat(k)=0`, so (4) and (7) give

```text
rhat(k)=p Ahat(k)Bhat(k)
       =p^(-1) zeta^(-k(a_0+b_0))
          (1+zeta^(-k d_A))(1+zeta^(-k d_B)),            (8)

d_A=a_1-a_0,                 d_B=b_1-b_0.                (9)
```

Both `k d_A` and `k d_B` are nonzero.  A two-term factor in (8) could vanish
only if an odd-order root of unity equalled `-1`.  Therefore

```text
rhat(k)!=0                    for every k!=0.             (10)
```

The same argument gives a sharp uniform amplitude floor.  Among nonzero
`m mod p`,

```text
|1+zeta^m|>=2 sin(pi/(2p)),                               (11)
```

with equality at `m=(p-1)/2` or `(p+1)/2`.  Hence

```text
|rhat(k)|>=4 sin^2(pi/(2p))/p.                           (12)
```

Equality is attainable by choosing both complement steps so that their
`k`-multiples are closest to the half turn.  Thus (12), unlike a mere Galois
nonvanishing statement, is a sharp per-colour floor.

## 3. Exact energy and its equality classes

The mean is

```text
rhat(0)=|S||T|/p=(p-2)^2/p=p-4+4/p.                     (13)
```

Define the charged energy

```text
E(S,T)=sum_(k!=0)|rhat(k)|^2
      =p^(-1) sum_c |r(c)-rhat(0)|^2.                   (14)
```

There are exactly two cases.  Call the complement steps **matched** when

```text
d_B=+-d_A.                                               (15)
```

Then the four sums in `A+B` have multiplicities `2,1,1`, so

```text
sum_c m(c)^2=6.                                          (16)
```

If the steps are distinct up to sign, all four sums are distinct and

```text
sum_c m(c)^2=4.                                          (17)
```

Unnormalized Parseval for `m` says

```text
sum_(k!=0) |sum_c m(c)zeta^(-kc)|^2
   =p sum_c m(c)^2-16.                                   (18)
```

Equations (8) and (18) give the complete classification

```text
E(S,T)=
  2(3p-8)/p^2,       if d_B=+-d_A,
  4(p-4)/p^2,        otherwise.                          (19)
```

At `p=13`, the `78` two-point complements split into six undirected step
classes of size `13`.  Thus among the `78^2=6,084` ordered pairs, exactly

```text
6*13^2=1,014                                           (20)
```

are step-matched and `5,070` are step-distinct.  Formula (19) becomes

```text
E=62/169             on 1,014 pairs,
E=36/169             on 5,070 pairs.                    (21)
```

In particular every pair satisfies

```text
max_(k!=0)|rhat(k)|^2
 >=E/12
 >=3/169.                                                 (22)
```

The equality classification in (21) is exact.  Equation (22) is the derived
universal maximum-mode floor; no claim that its final Cauchy bound is attained
is needed.

## 4. What multiplicity still forgets

The positive result does not reconstruct the two endpoint relations.  For
every `u in G`, the common-origin gauge

```text
(A,B) -> (A+u,B-u)                                       (23)
```

leaves `m=1_A*1_B`, hence `r`, exactly unchanged.  On proper two-point sets
the action is free, so every profile has at least thirteen labelled
realizations.  Swapping `A` and `B` is another invariance.

For example

```text
({0,1},{0,2})       and       ({3,4},{10,12})             (24)
```

have identical multiplicity profiles.  Thus (10) recovers charged carry
**spectrum**, not an allocated left/right carry pair, endpoint chronology, or
one private increment decomposition.  This is the multiplicity analogue of
the allocation-gauge boundary in THM-2634.

## 5. Full charged spectrum is quantitatively not purity

THM-2644 supplies a nonlinear fixed-branch gate for a nonnegative torsor
weight `mu`: with

```text
M=sum_c mu(c),        E_raw=sum_c mu(c)^2,
delta=M^2-E_raw,      R=(mu*mu)(0),                       (25)
```

the strict inequality `R>delta` forces identity mass on an odd torsor.  The
fully charged profile `mu=r` lies decisively on the opposite side.

From the two multiplicity shapes in Section 3,

```text
M=121,
E_raw=1129           on the 5,070 step-distinct pairs,
E_raw=1131           on the 1,014 step-matched pairs,     (26)
```

so `delta` is respectively `13,512` or `13,510`.  Reflection preserves the
`l^2` norm, hence Cauchy--Schwarz gives

```text
R=sum_c r(c)r(-c)<=E_raw<=1131<13510<=delta.             (27)
```

The complete finite bank sharpens the separation: `max R=1131`, and the
smallest deficit `delta-R` is `12,379`.  Therefore THM-2644's robust gate
cannot fire on any THM-2645 multiplicity profile, despite all twelve charged
characters surviving with the floors (12), (21), and (22).

This proves that full character spectrum and fixed-branch purity are
transverse invariants.  Another Fourier measurement cannot by itself turn
this same dense profile into a branch.  Conditioning or thinning the physical
transition is necessary.

## 6. LRC interface and stopping boundary

The theorem identifies a precise conditional consumer.  If a physical LRC
construction supplies a same-base positive two-edge transition table whose
allowed increments are two eleven-sheet `C_13` relations, then:

1. THM-2642 says support services every carry and cannot select one;
2. THM-2645 says the exact representation multiplicity has all twelve
   charged carry colours and the energy floor (21); but
3. the origin gauge (23) still prevents recovery of the two endpoint carries
   without a pair reference.

No current LRC theorem supplies that physical table.  The existing
eleven-sheet rows are static/coefficient fibres, and the physical predecessor
digit is target-neutral rather than a lawfully transported carry coordinate.
Signed Fourier or tomographic rows are not positive transition counts.
Therefore (10)--(22) cannot yet be inserted into a canonical THM-2334 current,
THM-2634 endpoint pair, or THM-2637 local system.

No row is excluded and the scalar ledger remains `165`.

## 7. Exact companion

Run

```bash
python3 04-computation/lrc14_eleven_sheet_multiplicity_spectrum_thm2645.py
python3 -O 04-computation/lrc14_eleven_sheet_multiplicity_spectrum_thm2645.py
```

Both executions byte-match

```text
05-knowledge/results/lrc14_eleven_sheet_multiplicity_spectrum_thm2645.out.
```

The dependency-free exact referee uses the integral group-ring criterion

```text
sum_(j=0)^12 n_j zeta^j=0  iff  n_0=...=n_12            (28)
```

for integer coefficient vectors.  It exhausts all `6,084` ordered complement
pairs and checks:

1. all `73,008=6,084*12` charged character numerators are nonzero;
2. the exact step-class and energy censuses (20)--(21);
3. normalized physical-space and group-ring Parseval energies agree;
4. all `79,092=6,084*13` common-origin gauge transforms preserve the profile
   and form free thirteen-element orbits; and
5. all `6,084` dense multiplicity profiles obey (26)--(27), with exact
   `max R=1,131` and minimum return-gate deficit `12,379`.

Every logical guard survives optimized Python.  The LF-normalized SHA-256
hashes are declared in the frontmatter.

An independent immutable audit of the original theorem rederived the
complement factorization, the normalized convolution sign and factor `p` in
(8), all-mode nonvanishing and
the sharp odd-root amplitude floor, both Parseval energy classes and their
`5,070+1,014` census, the maximum-mode deduction, and the freeness and example
for the common-origin gauge.  It also confirmed that the same-base positive
transition-count table is genuinely absent from current LRC canon.  Normal
and optimized executions byte-match the stored transcript, and the declared
LF-normalized hashes were independently reproduced.

An independent re-audit of the post-promotion transversality corollary
recomputed `E_raw=1,129/1,131`, `delta=13,512/13,510`, the reversal
Cauchy bound, the attained `max R=1,131`, and the sharp minimum deficit
`12,379` across all `6,084` pairs.  It replayed normal and optimized modes,
matched the stored transcript, and independently reproduced the new hashes.

QED.
