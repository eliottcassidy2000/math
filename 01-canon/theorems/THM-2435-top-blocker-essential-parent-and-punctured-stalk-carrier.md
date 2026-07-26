---
id: THM-2435
title: "Top-blocker essential parent and punctured-stalk carrier"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; UNDER INDEPENDENT
  HOSTILE AUDIT. On each of the three surviving deep-c_3 valuation
  shapes, THM-2427 supplies a common-91-root parent image of mass at
  least 4/7. THM-2431 caps the locus where the guard and five ordinary
  unit words cover the whole stalk at 3/7. Hence parents of mass at
  least 1/7 have a literal gap in the six-unit cover which must be
  filled by one of the one or two exact-depth top blockers. The
  punctured physical packet has mass at least 1/637; one fixed actual
  top-blocker label has mass at least 1/1274 and parent mass at least
  1/14. Its seven-root profile is at most one-sheeted, so every
  septimal Fourier residue class has exact energy rho/7. In particular
  every six nonzero classes have energy at least 1/8918 each, total at
  least 3/4459. This does not force a thirteen-unit frequency, a unique
  owner, a terminal current, a valuation-type exclusion, or LRC(14).
source: codex-2026-07-26-top-blocker-essential-parent
depends_on:
  - THM-2427-guard-top-thirteen-root-capacity-and-residual-types
  - THM-2431-repeated-step-rounding-exclusion-of-guard-top-zero-blocker-types
  - THM-2432-guard-top-pair-cage-and-low-blocker-residual-exclusion
related:
  - THM-2398-prime-cyclic-rational-restoration-dichotomy
  - THM-2409-unfiltered-septimal-source-orbit-and-real-word-obstruction
  - THM-2424-coprime-common-root-crt-and-unit-residue-spectrum
script: 04-computation/lrc14_top_blocker_essential_parent_thm2435.py
output: 05-knowledge/results/lrc14_top_blocker_essential_parent_thm2435.out
script_sha256: 99a442b35d21f6abd44b45d96cbc738a9e8d93bab5d9170baf810038c3ece9c3
output_sha256: e6cf58a236e5ff6c48fcc5d1d0fb75627d388758a014b6ef5383cb8e8da0b867
hash_basis: working-tree bytes (LF)
---

# THM-2435 -- the surviving shapes have a physical top-blocker carrier

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; UNDER INDEPENDENT
HOSTILE AUDIT.**

THM-2431 eliminates the zero-top-blocker branch because its required
parent image is too large for an exact six-unit tiling. On the three
blocker-bearing survivors, the same mass mismatch cannot be a
contradiction: it becomes a positive, physically labelled carrier:

```text
parent image at least 4/7
  - six-unit exact-tiling locus at most 3/7
  = essential top-blocker parent mass at least 1/7.                 (1)
```

Disintegrating one such gap through the common `91`-root stalk gives
mass `1/637`. A fixed actual top blocker retains at least half of it.
Because that blocker occupies one sheet of every septimal root fibre,
the carrier has exact, positive energy in all six nonzero septimal
Fourier residue classes.

## 1. The three typed survivors

Retain THM-2427's primitive scalar cover on the deep-`c_3`,
top-guard side:

```text
nu_7(H)=M<nu_7(c_3),                    t=5.                       (2)
```

After THM-2431 and THM-2432, the exact residual is

```text
M=0:   (k,t,b,W)=(1,5,1,8), (2,5,2,9),

M>0:   (k,t,b,W)=(2,5,1,8).                                    (3)
```

Thus

```text
J={j in {1,2}:nu_7(c_j)=M}
```

has cardinality

```text
b=1,2,1
```

in the three cases. Put

```text
N=7^(M+1),                    L=91,

u_0=H/7^M,                   u_i=q_i/7^M,

v_j=c_j/7^M                 for j in J.                           (4)
```

The six `u` speeds are units modulo `91`. Each top blocker has

```text
v_j=13w_j,                  7 not|w_j.                            (5)
```

## 2. The large parent image and the six-unit locus

As in THM-2427, let

```text
A=D_(C_1)^c intersection D_(C_2)^c intersection D_(C_3)^c,

P=T_N(A).                                                         (6)
```

Endpoint and inherited exceptional images are null and will always be
ignored. Since

```text
A subset T_N^(-1)(P)
```

and `T_N` preserves Haar measure under inverse image,

```text
mu(P)>=mu(A)>=4/7.                                                (7)
```

For `Y in P`, use the physical common-root coordinate

```text
z_s=(Y+s)/L,                       s in Z/LZ.                      (8)
```

The common translation `floor(Ny)` in the labelled CRT index is
absorbed into `s`; no guard normalization is used in the packet below.
Let

```text
U=E_(u_0) union D_(u_1) union ... union D_(u_5).                  (9)
```

The guard has `26` points and the five ordinary words have `13` each,
so their total incidence on (8) is exactly `91`.

Define the global exact-unit locus

```text
mathcal E={Y: U contains every z_s in (8)}.                       (10)
```

Whenever (10) holds, the six masks are an exact one-fold tiling.
THM-2430 puts it in the corrected `62`-tiling atlas, and THM-2431's
fixed-pair sharp phase bound applies to these same fixed speeds.
Therefore

```text
mu(mathcal E)<=3/7.                                               (11)
```

Equations (7) and (11) give

```text
mu(P minus mathcal E)>=1/7.                                      (12)
```

## 3. Puncture the common stalk

Define the physical Boolean packet

```text
F={z:T_L z in P,                    z notin U}.                    (13)
```

The map

```text
(Y,s) -> (Y+s)/L
```

is an almost-everywhere measure-preserving bijection from Haar measure
times normalized counting measure on `Z/LZ` to physical Haar measure.
Every parent in `P\mathcal E` has at least one gap in (9). Hence

```text
mu(F)
 =integral_P #{s:z_s notin U}/L dY
 >=1/(7*91)
 =1/637.                                                         (14)
```

THM-2427's inequality `W>k` has already proved that the top words cover
every septimal bin. Thus every point of `F` is covered by at least one
of the actual top blockers:

```text
F subset union_(j in J) D_(v_j)                     a.e.          (15)
```

Low blockers can also be present and the two top blockers can overlap.
Equation (15) is not a unique-owner assertion.

For `j in J`, put

```text
F_j=F intersection D_(v_j),

P_j={Y in P: some z_s belongs to F_j}.                            (16)
```

The sets `P_j` cover `P\mathcal E`. Since `b<=2`, one fixed actual
blocker label satisfies

```text
mu(P_j)>=1/(7b)>=1/14.                                           (17)
```

Each parent in `P_j` supplies at least one physical root, so

```text
rho:=mu(F_j)>=mu(P_j)/91>=1/1274.                                (18)
```

When `b=1`, the sharper bounds are

```text
mu(P_j)>=1/7,                     rho>=1/637.                     (19)
```

The label in (17)--(19) is fixed before any Fourier colour is chosen.

## 4. Exact septimal polyphase energy

Let `f=1_(F_j)` for the fixed label in Section 3 and define, for
`ell in Z/7Z`,

```text
C_ell(t)
 =1/7 sum_(r in Z/7Z)
       f((t+r)/7) zeta_7^(-ell r).                               (20)
```

Because `f<=1_(D_(v_j))` and (5) holds, the seven-root profile in
(20) has at most one occupied sheet for almost every `t`. Therefore
every colour has the same pointwise squared magnitude:

```text
|C_ell(t)|^2
 =1/49 sum_r f((t+r)/7).                                        (21)
```

Integration and the standard polyphase identity give the exact law

```text
sum_(n congruent ell mod 7)|fhat(n)|^2
 =integral |C_ell(t)|^2dt
 =rho/7                                                       (22)
```

for every `ell`, including zero. In particular, for each of the six
nonzero septimal colours,

```text
sum_(n congruent ell mod 7)|fhat(n)|^2
 >=1/8918,                                                       (23)
```

and the total charged energy is at least

```text
6/8918=3/4459.                                                   (24)
```

In the two one-blocker shapes, (19) sharpens (23)--(24) to

```text
per colour >=1/4459,                  total >=6/4459.             (25)
```

Thus every nonzero residue class modulo seven contains a genuine
physical Fourier frequency of the same fixed labelled packet.

## 5. Bounded lifts and the exact stopping boundary

The set `A`, its image `P`, all combs in (9), and `F_j` are rational
finite-step sets. For a nonzero class `ell mod 7`, group the oriented
jumps of `f` along the progression `ell+7n` into `L_ell` effective
nodes. THM-2424's endpoint-Prony lift gives a surviving frequency with

```text
n_physical congruent ell mod 7,

0<|n_physical|<=7 ceil(L_ell/2)-1.                               (26)
```

This is a physical frequency, not a relation-current address.

No conclusion modulo thirteen follows from (22). The sharp hostile is

```text
f=1_(D_13).
```

It is `1/13`-periodic, so every Fourier frequency is divisible by
thirteen. Nevertheless, for each `ell=1,...,6`, the frequency

```text
n=13(7-ell)
```

lies in class `ell mod 7` and has nonzero central-interval
coefficient. Hence all six septimal classes can survive while every
frequency is target-dark modulo thirteen. A target-root phase,
common endpoint reference, or genuine CRT intertwiner is still
required for a `91`-unit current.

## 6. Exact companion

Run

```text
python3 04-computation/lrc14_top_blocker_essential_parent_thm2435.py
python3 -O 04-computation/lrc14_top_blocker_essential_parent_thm2435.py
```

The dependency-free companion:

- exhausts all `352,947` labelled placements of a two-bin guard and
  five ordinary singleton words on `F_7`;
- checks the complete unit-gap histogram and all one-/two-blocker
  covering and essential counts;
- verifies the exact `1/7`, `1/637`, `1/14`, and `1/1274` mass
  invoices;
- verifies the singleton-sheet DFT norm and every quantitative energy
  floor in (23)--(25); and
- checks the six explicit target-dark frequencies of Section 5.

The finite word bank contains one-gap and two-gap configurations, so
the essential-blocker alternative is a genuine sharp branch, not a
counting artefact. The universal factor `1/2` in (17) and the
one-root factor `1/91` in (18) are optimal from the present data:
two blocker labels may split the essential parents evenly, with only
one puncture on each parent.

Normal and optimized runs must reproduce

```text
05-knowledge/results/lrc14_top_blocker_essential_parent_thm2435.out
```

byte-for-byte.

## 7. Consequence and scope

Every one of the three live deep-`c_3` valuation shapes now carries:

```text
one fixed actual top-blocker label,
positive punctured physical mass,
all six septimal Fourier residue energies,
and a finite physical lift in every charged septimal class.       (27)
```

This does not identify a unique scalar-cover owner, preserve a
preselected terminal word, force a thirteen-unit frequency, exclude
one of the three shapes, remove a scalar profile row, or prove
LRC(14). The next decisive sidecar must couple (27) to a nonflat
target/root coordinate without losing the fixed blocker label.

## 8. Independent audit

Independent hostile audit is in progress. QED conditional on that
audit.
