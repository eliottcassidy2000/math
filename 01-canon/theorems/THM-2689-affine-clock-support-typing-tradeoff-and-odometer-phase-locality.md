---
id: THM-2689
title: "Affine clock support/typing tradeoff and odometer phase locality"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY REPLAYED.  On the complete
  THM-2584 three-tooth rail envelope, let K be the states whose intrinsic
  shallow-to-owner seven-clock edge is nonconstant.  For either affine sign,
  K intersect T^(-1)K intersect T^(-2)K has positive length exactly when the
  shift beta is neither 0 nor 1/2 modulo one.  A global seven-clock label
  permutation exists exactly when 91 beta is integral.  Consequently every
  nontrivial THM-2657 carry/root lift has positive intrinsic three-event
  support, but none is globally seven-clock-covariant.  The support no-go for
  unshifted dilation is therefore a phase resonance, while the odometer escape
  necessarily needs a local phase sidecar.  Present factors, delayed words,
  primitive units, semantic endpoints, and LRC(14) row exclusion do not follow.
source: root/lrc-clock-referee-2026-07-28
depends_on:
  - THM-2584-b-word-depth-five-absolute-deep-root-tensor
  - THM-2657-odometer-carry-root-lift-nonsplit-extension-and-cech-cocycle
  - THM-2684-three-tooth-rail-envelope-diagonal-arrival-law-and-full-dilation-nilpotence
related:
  - THM-2640-predecessor-carry-private-root-atlas-and-target-action-clutching-no-go
  - THM-2680-dilation-reversed-two-edge-clock-fibre-products-and-source-drift-boundary
script: 04-computation/lrc14_affine_three_tooth_full_clock_classification.py
output: 05-knowledge/results/lrc14_affine_three_tooth_full_clock_classification.out
script_sha256: b68fd1c73bcddacded54bbe3579e28f45057a071105990c49509751c3469a77c
output_sha256: 860f80537567d246d32da7eb356432f8d9e24605bea152f0b8532f664c04db58
secondary_script: 04-computation/lrc14_half_clock_affine_handoff_positive_triple_probe.py
secondary_output: 05-knowledge/results/lrc14_half_clock_affine_handoff_positive_triple_probe.out
secondary_script_sha256: 3740dacd58739f006191e6f74dfbf3eefc624c18f0f0b911882003caac485cf7
secondary_output_sha256: 16315041702ab591abaa805196c27930a53640323468ea5e4bbf2277b47b20d1
hash_basis: LF-normalized bytes
---

# THM-2689 -- every physical odometer lift escapes support and loses global clock labels

**PROVED + VERIFIED-EXACT + INDEPENDENTLY REPLAYED.**

THM-2684 proves that ordinary dilation has no clock-legal third event on the
complete inherited rail bank.  That zero is sharp but nongeneric: translating
the dilation restores positive intrinsic support almost everywhere.  The same
translation usually destroys the global label transport needed to interpret
that support as one inherited chronological packet.  For the physical
THM-2657 odometer lifts, the two phenomena are mutually exclusive in the
strongest possible way.

## 1. The intrinsic carrier

Put

```text
p=13,                 q=7,                 R=13^6,
D(x)={13x},            T_(s,beta)(x)={s 13x+beta},
s in {+1,-1}.                                             (1)
```

THM-2684 identifies the exact union of the `324` positive THM-2584 rail
supports as

```text
B=[0,1/28) union [13/28,15/28) union [27/28,1).           (2)
```

Use the nearest seven-clock

```text
c_7(u)=floor(7{u}+1/2) mod 7.                             (3)
```

An inherited rail event at `u` stores the intrinsic edge

```text
c_7(Du) -> c_7(D^2u).                                    (4)
```

Let

```text
K={u in B : c_7(Du)!=c_7(D^2u)}.                         (5)
```

Cutting `(2)` at every boundary of the two clocks in `(4)` gives exactly
`148` open rational components and

```text
measure(K)=146/1183.                                     (6)
```

Endpoints are irrelevant to every positive-length statement below.

## 2. Exact affine-parameter classification

For both signs independently,

```text
measure(K intersect T_(s,beta)^(-1)K
          intersect T_(s,beta)^(-2)K)>0

iff beta is not 0 or 1/2 modulo 1.                       (7)
```

Here is a finite exact proof rather than a parameter sample.  Fix three
components `I_0,I_1,I_2` of `K` and the integer wrap branches `m_1,m_2`.
The corresponding points `(x,beta)` satisfy six strict rational affine
inequalities:

```text
x in I_0,
s 13x+beta-m_1 in I_1,
13^2x+(s13+1)beta-m_2 in I_2.                            (8)
```

Every lower bound on `x` is affine in `beta`, as is every upper bound.
Fourier--Motzkin elimination says that `(8)` has positive `x`-width exactly
on one open rational interval in `beta`.  The exact companion constructs
these intervals from interior witnesses and verifies their inequalities.
For `s=+1`, `264` greedy corridors plus `23` seam corridors cover

```text
(0,1/2) union (1/2,1);                                   (9)
```

for `s=-1`, the corresponding counts are `252+14`.  Exact interval union,
not decimal comparison, gives `(9)` in both cases.

At `beta=0` and `beta=1/2`, raw three-tooth returns still have total mass
`1/1183`; thus the exception is not geometric emptiness.  On every such raw
return the first edge `(4)` is diagonal.  This proves the reverse direction
of `(7)` and supplies hostile controls for confusing rail support with legal
clock support.

Because `B` is a finite union of actual positive rail supports, distributivity
shows that positive support in `(7)` contains some positive raw rail-labelled
triple.  It does not say that the three rail labels share the remaining
present, word, unit, or endpoint data.

## 3. The global clock-normalizer lemma

The boundaries of `(3)` are the coset

```text
{(2a+1)/14 : a in F_7}.                                  (10)
```

Therefore the circle map `u -> {s u+t}` permutes the seven clock cells if
and only if

```text
7t is integral.                                          (11)
```

Necessity follows because a cell permutation must preserve the boundary set
`(10)`; sufficiency follows by direct translation or reflection.  If
`y=T_(s,beta)x`, then

```text
c_7(Dy)=c_7(s D^2x+13 beta).                             (12)
```

Consequently the following shallow label is a fixed permutation of the
current owner label exactly when

```text
7*13*beta is integral.                                   (13)
```

Writing `13 beta=k/7`, the forward rule is `a -> s a+k`; the inverse
transport is `b -> s(b-k)`.  Formula `(13)` is a global statement.  A local
component may have a consistent clock word without satisfying it.

## 4. Universal odometer support/typing tradeoff

THM-2657 classifies every physical carry/root translation as

```text
beta=k/R,                   k=7 delta mod 13.             (14)
```

For a nonzero quotient label `delta`, `k` is a thirteenth-unit.  There are
exactly

```text
12*13^5=4,455,516                                           (15)
```

such residues modulo `R`.  None gives either exceptional value in `(7)`:
`k/R` is nonzero, and the odd denominator `R` cannot represent `1/2`.
Thus **every nontrivial physical lift has positive intrinsic three-event
support for both affine signs**.

On the other hand `(13)` becomes

```text
91 k/R is integral  iff  13^5 divides k.                 (16)
```

The clock-covariant residues are exactly

```text
{j 13^5 : 0<=j<13}.                                      (17)
```

They are the unique order-thirteen subgroup of the translation odometer and
all have `delta=0` in `(14)`.  No nontrivial carry/root lift is therefore
globally seven-clock-covariant.

Equations `(7)` and `(16)` give the promised exact tradeoff:

```text
nonzero physical root motion  =>  positive intrinsic support
                              +  no global clock-label action.  (18)
```

The companion also enumerates all residues modulo `R` and finds no violation,
but the divisibility proof above is uniform and explains the census.

## 5. Typed kernel controls and a half-clock boundary witness

The generator `beta=1/13` belongs to `(17)` and moves no quotient root.  It
has explicit positive labelled triples.  For `s=+1`, the rail keys are

```text
(12,1,3,12), (0,1,2,0), (6,1,3,12),                     (19)
```

on the component

```text
(398445/399854,56921/57122),       length 1/199927,       (20)
```

with stored edges `0->3,3->2,2->3`.  Clock transport is the identity and
the digit law is `j_next=h`.  For `s=-1`, the keys

```text
(12,1,3,12), (6,1,5,0), (6,1,3,12)                      (21)
```

give

```text
(386615/399854,55231/57122),       length 1/199927,       (22)
```

with stored edges `4->3,4->5,2->3`, reflected clock transport, and
`j_next=12-h`.  These witnesses show that the clock-normalizing kernel is
not itself support-empty.

There is a complementary boundary-scale construction.  For every `q>=3`,
put

```text
p=2q-1,  delta=1/(4q),  beta=1/(2q),
F(x)={px+beta},
x_*=delta-(p-3)delta/p^3.                                (23)
```

Then `(p+1)beta=1`, so `F^2=D_p^2`, and exact algebra gives

```text
F(x_*)=1/2+delta-(p-3)delta/p^2,
F^2(x_*)=3delta/p.                                       (24)
```

The orbit has tooth word low--central--low and intrinsic clock prefix
`floor(q/2)->0->1`.  At `q=7`, the literal THM-2584 keys

```text
(0,1,0,0), (6,1,1,0), (0,1,3,0)                         (25)
```

have one positive common component

```text
(14215/399854,2031/57122),         length 1/199927.       (26)
```

The weighted common-grid pullback and an independent local inverse-branch
intersection agree exactly.  But `beta=1/14` violates `(13)`: it shifts a
half clock cell and requires a phase sidecar.  It is not a THM-2657 lift.

## 6. Lawful local cycle and coefficient-only signal

Let `S=13^5` and use the lawful alternating lifts `+/-(S+1)`.  They exchange

```text
x_+=1/2+(S+1)/(14R),       x_-=1/2-(S+1)/(14R).          (27)
```

Their intrinsic edges are `4->3` and `3->4`; owner-to-next-shallow labels
glue, carries are `7,5`, future digits are both `6`, and quotient root steps
are `11,2`.  By `(16)` this two-cycle is necessarily phase-local, not the
restriction of a global seven-label action.

Four THM-2640 coefficient rows on the surrounding rails were rebuilt from
the weighted rail, present, delayed, carry, edge, and root data.  Clock
reflection `rho(Y)_ell=Y_(-ell)` gives

```text
(0,5,2,2,0,0,0)=10 rho(0,0,0,0,8,8,7),
(0,6,5,5,0,0,0)= 4 rho(0,0,0,0,11,11,8)          mod 13. (28)
```

All four are unit rows.  They have different carry/rail/edge/root
metadata, however, and `(28)` is only a coefficient-line covariance.  It
does not identify their physical supports or turn `(27)` into a packet map.

## 7. Exact scope and reproduction

The theorem changes the LRC carrier diagnosis:

* the THM-2684 zero at `beta=0` is an exceptional phase resonance, so adding
  a lawful odometer twist genuinely restores intrinsic three-event support;
* no nonzero physical root lift can transport the global seven-clock labels,
  so a phase-local chart or a different clock grammar is mandatory;
* the positive sets above stop before the present factor, delayed word,
  predecessor-carry/private-unit compatibility, semantic endpoint, and one
  common component across those coordinates.

Hence there is no row exclusion and no LRC(14) conclusion.

Run

```bash
python3 04-computation/lrc14_affine_three_tooth_full_clock_classification.py
python3 -O 04-computation/lrc14_affine_three_tooth_full_clock_classification.py
python3 04-computation/lrc14_half_clock_affine_handoff_positive_triple_probe.py
python3 -O 04-computation/lrc14_half_clock_affine_handoff_positive_triple_probe.py
```

Normal and optimized executions match their stored outputs byte for byte.
The first companion uses exact rational interval algebra throughout; the
second independently reconstructs `(26)` by two routes.  A separate parent
replay reproduced the frozen affine transcript.
