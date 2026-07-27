---
id: THM-2592
title: "Fallback-rail digit-diagonal pullback and primitive Bockstein"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.
  On the canonical typed row, the literal route-two pullback of the 84
  THM-2586 theta-zero rail cells against all thirteen normalized THM-2585
  target sections is positive in exactly 39 of 1092 cases: precisely the
  three fallback cells (s4,ell4,v,t)=(7,4/5/6,6,12), for every q.  After one
  global primitive-content reduction, 38 of the 39 positive slices have a
  nonzero first Bockstein in all six owner colours and 37 have unit septimal
  slice polynomial.  The construction is an actual common-x fibre product,
  but is not uniform across the other 81 rails and does not identify the two
  clock labels, a named old head, or a canonical semantic endpoint.
source: codex-2026-07-28-joint-rail-projector-pullback
depends_on:
  - THM-2584-b-word-depth-five-absolute-deep-root-tensor
  - THM-2585-saturated-normalized-target-projector-and-bockstein-noncommutation
  - THM-2586-depth-five-arrival-to-future-root-diagonal
related:
  - THM-2315-marked-target-gain-corolla-and-pairwise-composition-boundary
  - THM-2452-indicator-idempotent-aggregate-endpoint-restoration
  - THM-2459-four-atom-drift-and-root-service-coarsening
  - THM-2591-theta-zero-selector-cech-coboundary-and-c91-holonomy-no-go
script: 04-computation/lrc14_fallback_rail_digit_diagonal_pullback_thm2592.py
output: 05-knowledge/results/lrc14_fallback_rail_digit_diagonal_pullback_thm2592.out
script_sha256: e75123c09aafb3a2e500e3a68f37eed2dbed5466a24d181463d1e2b9e01f5edc
output_sha256: 33443b6c0cb901d4c10cfb79a7f9555d200d6a5bc9ffefdabb218b6c2ec11827
hash_basis: working-tree bytes (LF)
---

# THM-2592 -- fallback-rail digit-diagonal pullback and primitive Bockstein

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.**

## 1. The common physical carrier

Work on the canonical typed row

```text
(H,q1,...,q5,c1,c2,c3)
 =(1,14,27,40,53,66,13,13^3,2*13^5),                    (1)
```

with

```text
T=297836897838480,   d=13^5,   R=13^6.                 (2)
```

THM-2584 writes its `sigma={b}`, depth-five packet in the route-two
physical coordinate `x` as

```text
13 U(x)V(x-s4/13),                                     (3)
```

where

```text
U=P_d(1_Qb P_169 1_E),   V=P_d 1_E.                    (4)
```

Before integration retain

```text
ell4 : frac(169x) lies in the ell4-th owner window,
v    = floor(13x),
t    = floor(26x) mod 13.                              (5)
```

For `s4!=0`, THM-2586's deterministic theta-zero choice is

```text
(v,t)=(0,0),
```

except at

```text
(s4,ell4)=(7,4),(7,5),(7,6),                           (6)
```

where `(v,t)=(6,12)`.  In every case `7t=v mod 13`.

THM-2585 has a different word and a different target action.  Let

```text
F_(ell5,s5)(x)
```

be its present `sigma={a}` source/target packet, let

```text
Delta_r(x)=d(c3 x-r/13),                               (7)
```

and let `Q^+_(ell5,h)` denote its delayed word restricted to root digit
`h`.  Thus

```text
Q^+_(ell5,h)(Rx)
 =Q^+_ell5(Rx) 1_(floor(13 frac(Rx))=h).                (8)
```

For `q in F_13` put `s5=-q`.  The joint tensor is

```text
J_(s4,ell4,q;ell5,r)
 =13 int U(x)V(x-s4/13)
       1_(ell4,v,t)(x)
       F_(ell5,-q)(x) Delta_r(x)
       Q^+_(ell5,v)(Rx) dx,                            (9)
```

where `(v,t)` is the selected cell above and `1_(ell4,v,t)` denotes all
three restrictions in (5).

Although `U,V` are Perron averages, (9) is still a finite Boolean
fibre-product bank.  With THM-2586's preimages

```text
X_a=(x+a)/d,
Z_(a,beta)=(X_a+beta)/169,
Y_e=({x-s4/13}+e)/d,
```

its exact expansion is

```text
J=13/(169d^2) sum_(a,beta,e)
    mu(S_(a,beta,e) intersection
       F_(ell5,-q) intersection Delta_r intersection
       {Q^+_(ell5,v)(Rx)=1}).                          (9a)
```

Here `S_(a,beta,e)` is THM-2586's literal Boolean ancestry sheet with the
three root/owner indicators already inserted.  Thus positivity of any entry
of (9) really supplies a positive Boolean pullback sheet; no probabilistic
coupling is being chosen after the fact.

Equation (9) is the required common carrier.  It is formed **before** the
THM-2584 rail marginal and before the THM-2585 target/deep transforms.
Moreover the two occurrences of `v` give the exact time identity

```text
floor(13x)=v=floor(13 frac(13^6 x)).                    (10)
```

Thus (9) is an actual fibre product over the same physical `x`, not a
composition inferred from equal marginal labels.  This is precisely the
distinction required by THM-2315.

## 2. One global integer normalization

Write the integer numerator profiles of (4) as

```text
U=N_U/(169d),   V=N_V/d.                                (11)
```

The exact prefix formula for the last factor in (9) gives integers
`A_(s4,ell4,q;ell5,r)` such that

```text
J_(s4,ell4,q;ell5,r)
 =13 A_(s4,ell4,q;ell5,r)/(169 d^2 R T).                (12)
```

No cell is normalized separately.  Across the complete
`84*13*7*13` bank, with the identically zero `r=0` cells retained, the
nonzero `A` have global content

```text
gcd(A)=25465440
      =2^5*3*5*7*11*13*53.                             (13)
```

Including the route-two factor `13`, the physical raw content is

```text
g=331050720
 =2^5*3*5*7*11*13^2*53.                                (14)
```

Define the one globally primitive tensor

```text
a=A/25465440.                                           (15)
```

Then every joint mass has the common scalar form

```text
J=(2/202345849311155871624914247) a.                    (16)
```

The serialization digest of the full primitive bank, ordered by
`(s4,ell4,q,ell5,r)` in the companion's deterministic rail order, is

```text
1286b1f2baa1f05299d93a1074db19afc3f99eb3d44500fd1de39c97c62c300c. (17)
```

This single normalization is load-bearing.  Re-primitivizing the 1092
pair slices separately would define different mod-13 objects and is not
used here.

## 3. Exact support: only the fallback rail attaches

The delayed word has the exact digit gate

```text
Q^+_(ell5,0)=empty                                      (18)
```

for every `ell5=0,...,6`, while every `Q^+_(ell5,6)` has positive
measure.  This is checked directly as a half-open interval identity on the
grid `T`; it is stronger than an integral cancellation.

Every one of the 81 primary THM-2586 cells has `v=0`.  Equations (8) and
(18) therefore make its complete pullback (9) identically zero for all
thirteen `q`.  The three fallback cells have `v=6`, and exact integration
gives

```text
sum_(ell5,r) J_(s4,ell4,q;ell5,r)>0

iff

s4=7, ell4 in {4,5,6}, (v,t)=(6,12), q arbitrary.       (19)
```

Consequently exactly

```text
39/1092                                                  (20)
```

rail-by-target-section pullbacks are positive.  Their unscaled `A` totals
range from

```text
163706661702634491201720960
```

to

```text
315647624785807202076027840.                             (21)
```

At the fine `(ell5,r)` level there are exactly

```text
1404/91728                                               (22)
```

positive entries.  Among the 39 positive pair slices, the number of
positive fine entries has histogram

```text
24:7, 36:25, 48:7.                                      (23)
```

Thus the attachment is genuine but maximally localized: the same digit
gate that proves positivity on the fallback cell annihilates the 81-cell
primary rail.

## 4. The globally primitive first Bockstein

For a positive pair slice define

```text
Y_ell5=sum_(r=1)^12 a_(s4,ell4,q;ell5,r) r^(-1) mod 13,

Y(z)=sum_(ell5=0)^6 Y_ell5 z^ell5
       in F_13[z]/(Phi_7).                              (24)
```

As in THM-2585, the deep augmentation identity follows from the literal
`r=0` zero.  If `kappa in F_7^*`, the first Bockstein of the corresponding
primitive deep cycle is

```text
beta(D^(kappa))=Omega Y(zeta_7^kappa),

Omega=(1,2,...,12) mod 13.                              (25)
```

The computation of (24) is made anew from the primitive tensor (15); no
THM-2585 marginal Bockstein is inherited through conditioning.

Exactly one positive slice has zero Bockstein:

```text
(s4,ell4,v,t,q)=(7,6,6,12,10),
(Y_0,...,Y_6)=(0,0,0,0,0,0,0).                         (26)
```

Every other positive slice has nonzero (25) for every one of the six
nonzero owner colours.  Hence

```text
38/39 positive slices,
228/234 positive-slice owner profiles,                  (27)
```

survive the first Bockstein.

Multiplication by `Y` in `F_13[z]/(Phi_7)` is invertible on 37 of the 39
positive slices.  Besides the zero slice (26), the unique nonzero nonunit
is

```text
(s4,ell4,v,t,q)=(7,4,6,12,11),
(Y_0,...,Y_6)=(9,11,0,10,0,0,0),                       (28)
```

whose multiplication determinant is zero mod 13.  Its six individual
owner-colour Bocksteins are nevertheless all nonzero.

## 5. Independent exact controls

The companion uses only integer and `Fraction` arithmetic and contains the
following separate checks.

1. THM-2584's tensor is rebuilt by its root-chart and physical-`x` routes;
   all entries agree before this pullback is formed.
2. The selected physical rail pieces re-sum exactly to the stored
   `K_(ell4)(s4,v,t)` numerators.
3. Three deep-comb restrictions are computed both by direct nonwrapping
   window intersection and by the inherited complement construction.
4. Eight positive delayed overlaps are computed both by one weighted
   prefix sweep and by grouping equal integer profile levels and invoking
   THM-2585's original unweighted sweep term by term.
5. Three unconditioned cells recover the original THM-2585 delayed-overlap
   formula before the rail factor is inserted.
6. Normal and optimized execution reproduce the stored transcript exactly.

The positive support census, both contents, primitive scalar, primitive
digest, support histogram, Bockstein census, and the two exceptional slices
(26)--(28) are executable assertions rather than printed approximations.

Reproduce with

```text
python 04-computation/lrc14_fallback_rail_digit_diagonal_pullback_thm2592.py
python -O 04-computation/lrc14_fallback_rail_digit_diagonal_pullback_thm2592.py
```

## 6. Exact scope

The theorem supplies a literal common-`x` span joining three THM-2586
fallback rail cells to every normalized THM-2585 target section, and a
globally primitive nonzero Bockstein on 38 of those 39 sections.  It is the
first actual attaching carrier between these two tensors; THM-2315's
pair-marginal hostile does not apply after (9) has been formed.

It does **not** prove any of the following.

- The other 81 theta-zero rail cells attach.  They are exactly zero by
  (18), so no uniform 84-cell or Cech filler follows.
- The owner clocks `ell4` and `ell5` are the same coordinate.  Both are
  retained, not identified.
- The future digit in (10) is the preselected named `k_a` boundary, the old
  semantic head, or a THM-2305 endpoint phase.
- The globally primitive Bockstein selects the same hidden Perron Boolean
  sheet as any separately chosen positive mass witness.
- A relation residue, all-coordinate unit address, all-165 conclusion, row
  decrement, or LRC(14) closure follows.

The sharp next problem is therefore relative rather than marginal: either
transport the three fallback attachments across the digit-zero hole (18),
or prove that the hole represents a genuine semantic/clock obstruction.

QED for the stated candidate.
