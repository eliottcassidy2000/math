---
id: THM-2585
title: "Saturated normalized target projector and Bockstein noncommutation"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.
  On the one fixed globally primitive THM-2571 digit-diagonal carrier, the
  normalized inverse target Fourier projector is coefficientwise integral:
  its numerator is exactly 13 times one literal target-shift slice.  Every
  one of the thirteen slice polynomials is a unit in the septimal algebra,
  so all 78 nonzero-owner/shift Bocksteins survive, with product-basis
  support histogram 48:8, 60:32, 72:38.  Target translation permutes these
  charged sections.  The operation is not integral on arbitrary independent
  representatives, even when all thirteen lie in one nonzero Cayley-cokernel
  coset and every unnormalized Fourier numerator fills.  This is one
  coefficient-level common-carrier theorem, not a positive semantic
  endpoint, a THM-2334 current, an all-165 statement, a row exclusion, or an
  LRC(14) conclusion.
source: root-holotopy-2026-07-28-saturated-target-projector
depends_on:
  - THM-2532-cyclic-tournament-cayley-algebra-and-chi7-even-quotient
  - THM-2571-deep-colour-cayley-filling-bockstein-and-norm-curvature-split
related:
  - THM-2550-canonical-typed-row-double-nondegeneracy
  - THM-2573-logarithmic-abel-normal-and-common-endpoint-jump-pairing
  - THM-2579-socle-flat-target-torsor-and-integral-difference-filling
  - THM-2580-hasse-bockstein-carry-tower-and-salem-local-unit-boundary
script: 04-computation/lrc14_saturated_target_projector_bockstein_thm2585.py
output: 05-knowledge/results/lrc14_saturated_target_projector_bockstein_thm2585.out
script_sha256: 616f4b1ee46c3fcb3eea14f137568eaf7e2638d2227bc7cdd1572385ca7cdfc4
output_sha256: a8d7ca53e77021eedea4f1ed42d213f7985f6d62d499bbc0ce93ff9394d0820a
hash_basis: LF-normalized bytes
---

# THM-2585 -- saturation restores every target-shift section

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT
PENDING.**

THM-2571 proves that the thirteen target-colour cycles on its canonical
primitive carrier have one target-flat nonzero Cayley-cokernel class.
THM-2579 records the resulting warning: every integral target difference and
every unnormalized target Fourier numerator fills.  The normalization by
`13`, however, is not innocuous.  On the actual common Fourier carrier it is
an exact coordinate projector:

```text
one primitive target-shift tensor
  -> thirteen target-character cycles
  -> unnormalized inverse Fourier numerator = 13 times one shift slice
  -> primitive division by 13 restores a nonzero Bockstein.              (1)
```

Thus the target-flat class is not the absence of target structure.  It is the
image of thirteen charged sections after the 13-cyclotomic socle has collapsed
their differences.  The division is lawful only because all thirteen cycles
come from one fixed integral tensor.  A sharp same-coset hostile below shows
that no integral group-ring operation on unrelated representatives is being
constructed.

## 1. The restricted normalized projector is integral

Put

```text
M=Z[zeta_7,zeta_13].                                      (2)
```

Use exactly the one globally primitive THM-2571 tensor

```text
x_(ell,s,r) in Z,       ell in F_7, s,r in F_13,           (3)
```

obtained by clearing the complete digit-diagonal overlap once and then
dividing its global content `13` once.  For `kappa in F_7^*` and
`b in F_13`, THM-2571 defines the integral deep cycle

```text
c_m^(kappa,b)
 =sum_(ell,s,r) x_(ell,s,r)
    zeta_7^(kappa ell) zeta_13^(b s+m r),

sum_m c_m^(kappa,b)=0.                                    (4)
```

For `q in F_13`, form the unnormalized inverse target numerator and its
normalized section

```text
N_m^(kappa,q)=sum_(b in F_13) zeta_13^(q b)c_m^(kappa,b),

D_m^(kappa,q)=N_m^(kappa,q)/13.                            (5)
```

The quotient in (5) is coefficientwise integral on this carrier.  Indeed,
exact root orthogonality gives

```text
sum_b zeta_13^[b(q+s)] = 13  if s=-q,
                         0   otherwise.                   (6)
```

Substitution into (4), before any quotient or reduction, gives the identity

```text
N_m^(kappa,q)
 =13 sum_(ell,r) x_(ell,-q,r)
       zeta_7^(kappa ell) zeta_13^(m r),

D_m^(kappa,q)
 =sum_(ell,r) x_(ell,-q,r)
       zeta_7^(kappa ell) zeta_13^(m r) in M.              (7)
```

Thus `D^(kappa,q)` is not a formal rational average: it is the literal
`s=-q` target-shift slice.  Since THM-2571's carrier has
`x_(ell,s,0)=0`, summing (7) over `m` again gives zero.  Hence every section
lies in the same integral deep augmentation lattice as the original cycles.

Equation (7) is the whole integrality mechanism.  It would be invalid to
divide an arbitrary family of thirteen integral cycles by `13`; Section 5
exhibits a family where that quotient fails for every `q`.

## 2. Primitive division and the Bockstein do not commute

On the deep augmentation lattice write

```text
beta(a)=sum_(m=0)^12 m a_m                 in M/13M,       (8)

Omega=sum_(m=0)^12 m zeta_13^m.
```

In the power basis `1,zeta_13,...,zeta_13^11`,

```text
Omega=(1,2,3,4,5,6,7,8,9,10,11,12) mod 13,               (9)
```

and THM-2571 proves

```text
sum_m m zeta_13^(m r)=r^(-1)Omega,
zeta_13^a Omega=Omega                         mod 13       (10)
```

for every `r!=0` and every `a`.  Define the thirteen septimal slice
polynomials

```text
Y_(ell,q)=sum_(r=1)^12 x_(ell,-q,r) r^(-1)       mod 13,

Y_q(z)=sum_(ell=0)^6 Y_(ell,q)z^ell
       in R_7=F_13[z]/(Phi_7(z)).                          (11)
```

Equations (7) and (10) give the exact factorization

```text
beta(D^(kappa,q))=Omega Y_q(zeta_7^kappa).                 (12)
```

By contrast,

```text
beta(N^(kappa,q))=beta(13D^(kappa,q))=0.                  (13)
```

Thus the integral numerator is a Cayley boundary, while its lawful primitive
content reduction by `13` is obstructed.  There is no operation
`beta(N)/13` in `M/13M`; (12)--(13) say precisely that primitive saturation
and passage to the first Bockstein do not commute.

The target-flat formula of THM-2571 is the sum of the sections.  Since the
socle identity in (10) kills the target phase before the `s`-sum,

```text
beta(c^(kappa,b))
 =Omega [sum_q Y_q](zeta_7^kappa),                         (14)
```

independently of `b`.  Equations (12) and (14) are compatible: flatness is
socle contraction, not slice vanishing.

## 3. All thirteen canonical slice factors are units

The carrier in (3) is the fixed THM-2550(B)/THM-2309 (25) typed-row carrier,
not THM-2550(A)'s distinct `k=2` drift packet.  Its source row and clock are

```text
(H,q_1,...,q_5,c_1,c_2,c_3)
 =(1,14,27,40,53,66,13,2197,742586),

R=13^6,   T=297836897838480.                              (15)
```

The complete `1183`-cell interval reconstruction has `1092` positive cells,
global raw content `13`, primitive denominator

```text
110584755309142754640,                                    (16)
```

and primitive serialization digest

```text
a66ba96d31a33354468392b1dabc19865e6e925158efdab059fad9a98d4390f4.
                                                                    (17)
```

The exact slice vectors, determinants of multiplication by their reduced
classes in `R_7`, and owner colours with full 72-coordinate product-basis
support are:

| `q` | `(Y_(0,q),...,Y_(6,q))` | `det` mod 13 | full-support `kappa` |
|---:|:---|---:|:---|
| 0 | `(0,9,9,0,0,4,4)` | 7 | none |
| 1 | `(5,9,7,11,11,4,9)` | 6 | `3,4` |
| 2 | `(3,11,7,10,10,11,1)` | 5 | `1,3` |
| 3 | `(9,11,4,2,10,5,10)` | 2 | `2,3,4,6` |
| 4 | `(11,5,4,10,12,11,8)` | 7 | `1,2,3,5,6` |
| 5 | `(11,1,11,2,8,1,1)` | 10 | `2,5` |
| 6 | `(6,9,12,8,4,4,7)` | 7 | `1,2,3,6` |
| 7 | `(7,6,9,9,5,1,4)` | 7 | `1,4,5,6` |
| 8 | `(2,12,12,5,11,2,12)` | 10 | `2,5` |
| 9 | `(2,5,2,1,3,9,8)` | 7 | `1,2,4,5,6` |
| 10 | `(4,3,8,3,11,9,2)` | 2 | `1,3,4,5` |
| 11 | `(10,12,2,3,3,6,2)` | 5 | `4,6` |
| 12 | `(8,4,9,2,2,6,4)` | 6 | `3,4` |

Every determinant is nonzero.  Hence every `Y_q` is a unit in the etale
septimal algebra `R_7`, and every Galois owner evaluation
`Y_q(zeta_7^kappa)` is nonzero.  Consequently

```text
beta(D^(kappa,q))!=0
   for all (kappa,q) in F_7^* x F_13:       78/78.         (18)
```

In the product power basis of
`F_13[zeta_7,zeta_13]`, the exact support histogram of the 78 classes is

```text
support 48: 8 profiles,
support 60: 32 profiles,
support 72: 38 profiles.                                  (19)
```

There are no zero profiles.  Summing the thirteen slice vectors gives

```text
sum_q Y_q=(0,6,5,1,12,8,7),                              (20)
```

exactly THM-2571's target-flat global factor.

## 4. Translation gives a section atlas, not a canonical section

For `a in F_13`, translate the target-shift coordinate by

```text
(T_a x)_(ell,s,r)=x_(ell,s-a,r).                           (21)
```

Then

```text
c_m^(kappa,b)(T_a x)=zeta_13^(a b)c_m^(kappa,b)(x),

D_m^(kappa,q)(T_a x)=D_m^(kappa,q+a)(x).                  (22)
```

Thus the thirteen nonzero sections form a translation-permuted atlas.  The
operation preserves the owner colour `kappa`, the deep colour coordinate
`m`, and the common primitive coefficient lattice, but it chooses a
target-shift position `s=-q`.  Translation supplies no invariant preferred
`q` and no semantic orientation of that position.

This is the exact source/target/loss contract:

```text
source:     one primitive x_(ell,s,r) and its 13 target characters c^b;
map:        restricted inverse Fourier section (1/13)sum_b zeta^(qb)c^b;
target:     the integral target-shift slice s=-q;
preserved:  owner colour, deep cycle, primitive lattice, absolute shift slot;
lost/not supplied:
            target-character charge b, positivity, semantic endpoint,
            relation address, and a translation-invariant choice of q.       (23)
```

## 5. Sharp positive and hostile controls

### 5.1 A singleton common carrier restores the class

On the one-cell carrier at target shift `s=1` and deep root `r=1`, put

```text
a_m^(b)=zeta_13^(b+m).                                    (24)
```

Every absolute profile has the same nonzero class because
`zeta_13^b Omega=Omega`.  Its inverse target numerators are

```text
sum_b zeta_13^(q b)a^(b)=0,             q!=-1,
                              =13a^(0), q=-1.              (25)
```

Hence the normalized `q=-1=12` section is exactly `a^(0)` and has Bockstein
`Omega!=0`.  This is the smallest common-carrier positive signal for
(12)--(13).

### 5.2 Even one common cokernel coset is insufficient

Let `C=C_1` be THM-2532's integral Cayley derivative on the deep colour
lattice and let

```text
u=e_1-e_0,

w=Cu=(1,1,-2,2,-2,2,-2,2,-2,2,-2,2,-2).                 (26)
```

The vector `w` has coefficient content `1`, and `beta(u)=1`, `beta(w)=0`.
Define thirteen integral augmentation cycles by

```text
a^(0)=u+w,                  a^(b)=u for b!=0.              (27)
```

All thirteen represent the same nonzero Cayley-cokernel coset `[u]`.  Yet
their unnormalized inverse target numerators are

```text
sum_b a^(b)=13u+w,                              q=0,

sum_b zeta_13^(q b)a^(b)=w,                    q!=0.       (28)
```

Every numerator in (28) has Bockstein zero and hence an integral Cayley
primitive, but none has coefficientwise content `13`.  Therefore none of the
thirteen normalized quotients is integral.

This hostile is stronger than taking thirteen arbitrary classes: equality of
all cokernel classes and fillability of all integral Fourier numerators still
do not reconstruct a common Fourier carrier.  The literal common tensor in
(3)--(4), not merely its torsor of classes, is load-bearing.

## 6. Exact companion

Run

```bash
python3 04-computation/lrc14_saturated_target_projector_bockstein_thm2585.py
python3 -O 04-computation/lrc14_saturated_target_projector_bockstein_thm2585.py
```

Both executions must reproduce, after LF normalization,

```text
05-knowledge/results/lrc14_saturated_target_projector_bockstein_thm2585.out
```

byte-for-byte.  The companion imports only the fixed THM-2550(B) source
interval module, reconstructs the entire rational carrier independently of
THM-2571's stored output, performs the one global primitive reduction, and
checks:

- all `1183` cells, the `1092`-cell positive locus, raw content, denominator,
  and primitive digest;
- all `169` target-root orthogonality identities and all `169` translation
  section identities;
- the thirteen slice vectors and thirteen nonzero multiplication
  determinants;
- all `78` Bocksteins coefficientwise, their support histogram, and the exact
  full-support owner sets;
- the global slice sum and primitive-reduction noncommutation;
- all thirteen singleton sections; and
- the thirteen same-coset hostile numerators, their fillability, and the
  failure of integral normalization.

There are `4603` explicit checks, none implemented with `assert`.

## 7. Stopping boundary

The theorem is confined to one fixed THM-2571 primitive carrier, whose row
and word come from the THM-2550(B)/THM-2309 (25) large-clock construction.
MISTAKE-283 forbids identifying it with THM-2550(A)'s `k=2` drift packet merely
because their row and word labels agree.

The selected `q` is a target-shift **position**, not a surviving nonzero
target-character `b`.  After owner contraction, `D^(kappa,q)` is a signed
cyclotomic coefficient cycle, not a nonnegative physical endpoint.  No
argument here identifies it with THM-2573's lawful logarithmic Abel normal,
a THM-2305 semantic word endpoint, a THM-2334 relation current, an owner/root
intertwiner, or a scalar-cover row.  Nothing is asserted for the other `164`
rows.  No ledger entry decreases, and LRC(14) remains open. **QED.**
