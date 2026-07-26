---
id: THM-2431
title: "Repeated-step rounding exclusion of guard-top zero-blocker types"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT + THREE INDEPENDENT
  HOSTILE AUDITS ACCEPT; AWAITING THM-2430 PROMOTION. In a THM-2427
  residual with t=5 and b=0, the top guard and five top ordinary words
  give the THM-2430 exact common-91-root tiling on an image set of
  parent phases of Haar mass at least 4/7=52/91. The five normalized
  ordinary steps are fixed speed-residue data. Every one of the 62
  exact tilings repeats a step, so one fixed labelled pair repeats on
  every parent. Its centre difference, together with the
  nearest-integer addition defect, confines the parent to at most 50
  of the 91 equal-mass phase residues. The resulting 50/91<52/91
  contradiction excludes all t=5,b=0 valuation types. The deep-c_3
  residual becomes two types at M=0 and two at M>0. This is a
  regime-typed valuation reduction, not a scalar-profile decrement or
  a proof of LRC(14).
source: codex-2026-07-26-repeated-step-rounding-exclusion
depends_on:
  - THM-2427-guard-top-thirteen-root-capacity-and-residual-types
related:
  - THM-2430-guard-top-common-ninety-one-root-tiling-spectrum
script: 04-computation/lrc14_repeated_step_rounding_exclusion_thm2431.py
output: 05-knowledge/results/lrc14_repeated_step_rounding_exclusion_thm2431.out
script_sha256: 8fc14647db14516c02d961a2b75a14002fdd376b8bbe9f4ac6439f4b872ab4eb
output_sha256: 6fa4ad91ce4c5a99ea5722f488fd89bfeb0e19eb950c68b6b0310e7bf19b44cb
hash_basis: working-tree bytes (LF)
---

# THM-2431 -- a repeated step cannot persist on the required parent mass

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT + THREE INDEPENDENT
HOSTILE AUDITS ACCEPT; AWAITING THM-2430 PROMOTION.**

The corrected THM-2430 atlas has a feature which its local
realizability does not expose: every normalized tiling repeats an
ordinary step. For fixed physical speeds, that repeated step belongs
to one fixed labelled pair, independent of the parent phase. Its
centre difference converts the finite tiling atlas into a strict Haar
mass cap:

```text
fixed repeated pair
  -> nearest-integer residue bank of size at most 50
  -> parent mass at most 50/91
  -> contradiction with the blocker-safe image mass 4/7=52/91.     (1)
```

The word `fixed` is load-bearing. Choosing a different repeated pair
on different parents would give only a `66/91` union bound and no
contradiction.

## 1. Typed setting and the large parent set

Retain THM-2427's primitive scalar cover and suppose

```text
nu_7(c_3)>M,                    nu_7(H)=M,

t=#{i:nu_7(q_i)=M}=5,          b=#{j<=2:nu_7(c_j)=M}=0.           (2)
```

Put

```text
L=91,                          N=7^(M+1),

u_0=H/7^M,                     u_i=q_i/7^M.                       (3)
```

The six `u_i` are positive units modulo `91`, and the five ordinary
speeds are pairwise distinct.

Let

```text
A=D_(C_1)^c intersection D_(C_2)^c intersection D_(C_3)^c.
```

THM-2427 gives

```text
mu(A)>=4/7.                                                         (4)
```

Remove from `A` the finite union of pulled-back scalar-cover,
endpoint, and inherited orbit-identity null sets used in THM-2427,
and call the result `A_good`. Define

```text
P=T_N(A_good).                                                      (5)
```

The circle endomorphism `T_N` preserves Haar measure under inverse
image, while

```text
A_good subset T_N^(-1)(P).
```

Therefore

```text
mu(P)>=mu(A_good)>=4/7=52/91.                                      (6)
```

The image of every removed null set is again null because `T_N` is
finite-to-one.

For `y in A_good`, write `Y={Ny}`. THM-2430 composes the thirteen
inverse roots and the seven top-bin representatives. If

```text
z=7^M(y+h)/13+r/7,
```

then

```text
Lz
 =N(y+h)+13r
 =Y+floor(Ny)+Nh+13r.                                              (7)
```

As `(h,r)` ranges over `F_13 x F_7`, the integer part on the right
ranges bijectively over `Z/LZ`. The omitted `floor(Ny)` in an
unlabelled set description is only a common translation. Thus (7) is
the full common root fibre

```text
{(Y+s)/L:s in Z/LZ}.                                                (8)
```

THM-2427 has already proved that the guard and five top ordinary words
cover every top bin. Their sizes on (8) are `26,13,...,13`, with total
`L`; hence the top-word cover is an exact one-fold tiling for every
`Y in P`.

Low blockers can reappear away from the original `r=0` slice. The
proof does not call the whole `91`-stalk physically blocker-free and
does not use low-blocker absence there.

## 2. Intrinsic centres on the common root fibre

Define nearest-integer rounding by

```text
R(x)=floor(x+1/2).                                                  (9)
```

Away from the already removed endpoints, the ordinary mask `D_u` on
(8) is the length-thirteen progression

```text
{-u^(-1)R(uY)+j u^(-1): -6<=j<=6}          in Z/LZ.                (10)
```

Indeed its defining inequality is

```text
|uY+us-Lk|<13/2.
```

Normalize the guard by one fixed orientation and a
`Y`-dependent translation:

```text
s -> epsilon u_0 s+kappa(Y),                epsilon in {+-1}.     (11)
```

For ordinary label `i`, put

```text
a_i=epsilon u_0 u_i^(-1) mod L.                                  (12)
```

Its signed normalized step is `a_i`, its intrinsic centre is

```text
c_i(Y)=kappa(Y)-a_i R(u_iY),                                      (13)
```

and its canonical unsigned step

```text
d_i=min(a_i,L-a_i) in {1,...,45}                                  (14)
```

depends only on the fixed speed residues, not on `Y`.

## 3. The fixed repeated pair

THM-2430's complete `62`-tiling atlas has eighteen unsigned
five-step multisets, and every one contains a repetition. Because the
five values in (14) are fixed before the parent varies, there are
fixed labels `i!=j` and a fixed

```text
d=d_i=d_j.
```

No parentwise partition is needed. The only repeatable steps are

```text
d in {1,2,3,4,5,44,45};                                          (15)
```

step `30` never repeats.

There is a fixed `sigma in {+-1}` with

```text
a_j=sigma a_i.
```

Equation (12) then gives

```text
u_j=sigma u_i+Ln                                                   (16)
```

for an integer `n`. If `sigma=1`, pairwise distinctness gives
`n!=0`; if `sigma=-1`, positivity gives `n>0`.

Nearest-integer rounding satisfies

```text
R(x)+R(y)=R(x+y)+e_+,          R(x)-R(y)=R(x-y)+e_-,

e_+,e_- in {-1,0,1}.                                             (17)
```

Using (13), (16), and (17), the directed centre difference is

```text
c_i(Y)-c_j(Y)
 =tau d [R(LnY)+e(Y)] mod L,            tau in {+-1}.             (18)
```

## 4. The exact centre-difference banks

For each repeatable `d`, let `S_d` be the set of all directed
differences between intrinsic centres of two same-`d` ordinary blocks
among the `62` normalized THM-2430 tilings. Each `S_d` is symmetric.
The exact atlas gives:

| `d` | `|S_d|` | `|d^(-1)S_d+{-1,0,1}|` |
|---:|---:|---:|
| 1 | 28 | 50 |
| 2 | 24 | 44 |
| 3 | 4 | 8 |
| 4 | 6 | 12 |
| 5 | 8 | 16 |
| 44 | 8 | 12 |
| 45 | 4 | 12 |

The symmetry absorbs both `tau` in (18) and the reflected guard
normalization. Consequently every `Y in P` satisfies

```text
R(LnY) mod L
 in B_d:=d^(-1)S_d+{-1,0,1},              |B_d|<=50.              (19)
```

For every nonzero integer `n`, multiplication by `n` preserves Haar
measure. On a Haar-uniform `x`, the variable

```text
R(Lx) mod L
```

has exactly `L` equal-mass fibres: residues `1,...,90` each occupy
one interval of length `1/L`, while residue zero is the union of two
endpoint half-intervals of total length `1/L`. Hence (19) implies

```text
mu(P)<=50/91.                                                      (20)
```

Equations (6) and (20) contradict one another with the exact slack

```text
52/91-50/91=2/91.                                                 (21)
```

This excludes every THM-2427 residual satisfying `t=5,b=0`.

## 5. Residual valuation types

At `M=0`, THM-2427's three types reduce from

```text
(0,5,0,7), (1,5,1,8), (2,5,2,9)
```

to

```text
(1,5,1,8), (2,5,2,9).                                            (22)
```

At positive `M`, its four primitive types reduce from

```text
(1,5,0,7), (2,0,0,2), (2,5,0,7), (2,5,1,8)
```

to

```text
(2,0,0,2), (2,5,1,8).                                            (23)
```

This is a strict regime-typed valuation reduction. It does not remove
one of the `165` scalar profiles, identify an owner or terminal
current, or prove LRC(14).

## 6. Exact verification and stopping boundaries

Run

```text
python3 04-computation/lrc14_repeated_step_rounding_exclusion_thm2431.py
python3 -O 04-computation/lrc14_repeated_step_rounding_exclusion_thm2431.py
```

The dependency-free companion:

- reconstructs the `3,276` unit progressions, `182` eligible blocks,
  and all `62` covers by safe point-pivot Algorithm X;
- independently checks the `620=10*62` pair/triple invoice;
- directly checks the physical centre formula (10);
- extracts the seven directed `S_d` banks and every thickened size in
  the table;
- checks the nearest-integer sum and difference defects exactly on a
  large rational control bank;
- checks the equal-mass nearest-residue fibres and the exact `2/91`
  contradiction; and
- verifies the residual lists (22)--(23).

The union of all seven `B_d` has size `66`, so a proof which chooses
the repeated pair after seeing the parent would give only `66/91`.
This is the sharp quantifier boundary of the argument. Local physical
THM-2430 packets remain possible on smaller phase sets; the theorem
uses the scalar cover only through the large image set (6).

Normal and optimized runs must match

```text
05-knowledge/results/lrc14_repeated_step_rounding_exclusion_thm2431.out
```

byte-for-byte.

## 7. Independent audit

Three independent hostile audits reconstructed the `62` covers and
all seven banks, checked the fixed-pair quantifier, the signed centre
law, the common-root parent map including `floor(Ny)`, Haar image
measure, and the low-blocker scope. All accepted the proof conditional
only on promotion of THM-2430's atlas. QED conditional on that
dependency promotion.
