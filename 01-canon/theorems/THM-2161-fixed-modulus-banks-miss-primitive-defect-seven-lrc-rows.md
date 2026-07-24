---
id: THM-2161
title: "Fixed modulus banks miss primitive defect-seven LRC rows"
status: >
  PROVED. For every Q>=14 there is a primitive 13-speed row of AP defect
  exactly seven whose labelled residues agree with {1,...,13} modulo every
  2<=q<=Q. Consequently every modulus certificate m_q/q in that bank is at most
  1/14 and cannot prove the wider 3/41 gap. Nevertheless an explicit rational
  phase has minimum distance (L+8)/(8(L+14))>3/41, where
  L=lcm(2,...,Q), and its reduced denominator exceeds Q. Thus the row is
  loose: the obstruction is blindness of every fixed finite residue-prefix
  bank, not LRC hardness. An adaptive scale, relation, or phase sidecar is
  necessary.
source: codex-2026-07-24-critical-prefix-lrc-transfer
depends_on: []
related:
  - THM-566
  - THM-1002
  - THM-2072
  - THM-2163
---

# THM-2161 -- fixed modulus banks miss defect-seven rows

For a finite set `V` of positive speeds, write

```text
Gap(V)=max_(t in R/Z) min_(v in V)||vt||,              (1)

m_q(V)=max_(1<=a<q) min_(v in V)|av|_q,                (2)
```

where `|r|_q` is the least absolute residue modulo `q`. The elementary
modulus certificate is

```text
Gap(V)>=m_q(V)/q.                                      (3)
```

We prove that no bank with `q` bounded independently of the row can close the
current defect-at-least-seven route, even on primitive rows and even at the
wider threshold `3/41`.

## 1. A seven-defect AP mimic

Fix `Q>=14` and put

```text
L=lcm(2,3,...,Q),
V_Q={1,2,3,4,5,6,L+7,L+8,...,L+13}.                  (4)
```

The set has thirteen distinct positive elements, is primitive because it
contains `1`, and differs from the tight arithmetic progression

```text
A={1,2,...,13}                                        (5)
```

by replacing exactly the seven entries `7,...,13`. Thus its AP replacement
defect is exactly seven.

For every `2<=q<=Q`, one has `q|L`, and hence label by label

```text
L+j=j mod q,                         7<=j<=13.        (6)
```

It follows that

```text
m_q(V_Q)=m_q(A).                                      (7)
```

The classical tight-AP pigeonhole argument gives

```text
Gap(A)=1/14.                                          (8)
```

Indeed `t=1/14` attains `1/14`, while among
`0,t,2t,...,13t` the fourteen circular gaps sum to one. A gap of length at
most `1/14` joins two indexed points whose nonzero index difference is at
most thirteen, giving the matching upper bound.
Combining (3), (7), and the definition of the global maximum in (8) gives

```text
m_q(V_Q)/q=m_q(A)/q<=Gap(A)=1/14
                                         for every 2<=q<=Q. (9)
```

Therefore no modulus in the prescribed bank can certify a gap larger than
`1/14`, and in particular none can certify `Gap(V_Q)>3/41`.

## 2. The row has an explicit adaptive escape

Set

```text
D=L+14,
x=(L+8)/(8(L+14))=(D-6)/(8D).                        (10)
```

Since `8|L`, the number

```text
Dx=(L+8)/8
```

is an integer. Also

```text
0<x<1/8.                                              (11)
```

For a small speed `i=1,...,6`, both `ix>=x` and

```text
1-ix>1-6/8=1/4>x,
```

so

```text
||ix||>=x.                                            (12)
```

For a large speed `L+j`, `7<=j<=13`, put `k=14-j`, so
`1<=k<=7` and `L+j=D-k`. Because `Dx` is integral,

```text
||(L+j)x||=||kx||.                                    (13)
```

Now `kx>=x`, while

```text
1-kx>1-7/8=1/8>x.
```

Thus every distance in (13) is at least `x`. Equality already occurs for
the speed `1` and for `L+13`, so

```text
min_(v in V_Q)||vx||=x.                               (14)
```

Finally,

```text
x>3/41
 iff 41(L+8)>24(L+14)
 iff 17L>8,                                           (15)
```

which holds. Consequently

```text
Gap(V_Q)>=x>3/41>1/14.                                (16)
```

The family is therefore not a candidate counterexample or a near-tight row.
It is an exact hostile control showing that the finite modulus bank is blind
to an easy adaptive phase.

For the smallest instance `Q=14`, `L=360360` and

```text
x=22523/180187>3/41.                                  (17)
```

The successful phase is genuinely outside the bank. Write `L=8M`. Then

```text
x=(M+1)/(8M+14),
gcd(M+1,8M+14)=gcd(M+1,6)<=6.                        (18)
```

Its reduced denominator is at least `(L+14)/6`. Since the coprime integers
`Q` and `Q-1` both divide `L`,

```text
L>=Q(Q-1),
(L+14)/6>Q                         for Q>=14.          (19)
```

No floating-point estimate is used.

## 3. Consequence for the current LRC route

The source row (4) and the tight AP have identical labelled residue words at
every prescribed modulus, but their geometric gaps differ. Hence the map

```text
V -> (V mod q)_(q<=Q)                                 (18)
```

preserves every fixed-bank modulus certificate and destroys the adaptive
scale at which the witness (10) lives. Primitivity, defect seven, the
divisibility pins `7,...,13`, and the mod-seven odd-unit SDR do not repair
that loss: (4) inherits the AP
residue pattern for all fixed moduli and in particular contains a multiple of
each `q=7,...,13`.

This strengthens the fixed-bank warnings in THM-566 and THM-2072 in the exact
shape of the current defect-seven proposal. It does **not** show the
defect-seven LRC branch is hard or false. It shows that a proof must choose a
modulus from the row, use a relation/scale sidecar, or work directly with
phase geometry.

The row (4) is not a THM-523 Cover14 row: at modulus fourteen it has the
same residues as (5), and no zero residue. Thus the theorem is a no-go for a
fixed bank as a **uniform near-tight separator**, including the newly proposed
defect-seven route. It is not a no-go after an additional covering-row
restriction.

THM-2163 strengthens this loss statement in a different hostile family. Its
two primitive defect-seven rows preserve the complete labelled residue bank
through modulus thirteen, the **entire** coefficient-height-29 relation box,
and the scalar maximum, yet have different mod-17 margins. Thus neither
bounded relations nor scalar magnitude restore the adaptive residue lost by
a fixed bank. THM-2163 identifies a necessary richer state: a radix carry
together with its quotient-owner/termination sidecar. Those coordinates
alone do not yet give a finite complete state space, because the owner mask
records termination but does not bound digit depth.

## 4. The exact residue analogue only halves enumeration

The dyadic checksum of THM-2160 does have one limited residue analogue. Let
`q` be an odd prime, let `g` be a primitive root, and put `H=q-1`. Represent
a multiset of exactly `j` nonzero residues by multiplicities

```text
c_e at g^e,                         e in Z/HZ,
S(c)=sum_e e c_e mod H.                              (20)
```

If `a=nu_2(j)` and `2^(a+1)|H`, put

```text
delta=H/2^(a+1).                                      (21)
```

Global dilation by `g^delta` cyclically shifts the exponent multiplicities,
preserves `j`, and changes the checksum by

```text
j delta=H/2 mod H.                                    (22)
```

It also preserves `m_q`, because multiplication of every speed by a unit
merely permutes the numerators in (2). Therefore dilation by `g^delta` on
the lower checksum half and by `g^(-delta)` on the upper half is a
fixed-point-free involution pairing the two halves of this residue census.
For `q=17`, its divisibility hypothesis holds for every `1<=j<=13`.

This can halve a single-modulus enumeration, but it cannot be a loneliness
certificate by itself: `m_q` is constant on each pair, so certificate
success or failure is preserved rather than reversed. Cross-modulus CRT
compatibility and the adaptive phase remain required sidecars.

There is no intrinsic tournament in this argument. The faithful carrier is
the labelled residue bank together with the omitted adaptive scale; orienting
moduli or runners would discard the very sidecar exposed by (10).

QED.
