---
id: THM-2161
title: "Fixed modulus banks miss primitive defect-seven LRC rows"
status: >
  PROVED. For every Q>=14, a primitive 13-speed row of AP defect seven agrees
  labelwise with {1,...,13} modulo the whole bank while an explicit phase above
  3/41 has reduced denominator greater than Q. Thus the full labelled bank is
  information-theoretically blind on unrestricted defect-seven rows. A second
  defect-seven row is THM-523 covering and has m_q=0 throughout the same bank,
  yet the same phase succeeds. This second statement defeats the scalar m_q/q
  certificate family under covering; it does not identify everything visible
  to an arbitrary full-residue algorithm. Adaptive scale, relation, phase, or
  richer residue information is required.
source: codex-2026-07-24-critical-prefix-lrc-transfer
depends_on: []
related:
  - THM-523
  - THM-566
  - THM-1002
  - THM-2072
  - THM-2163
script: 04-computation/lrc14_fixed_modulus_bank_escape_thm2161.py
output: 05-knowledge/results/lrc14_fixed_modulus_bank_escape_thm2161.out
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

The exact regression

```text
python3 04-computation/lrc14_fixed_modulus_bank_escape_thm2161.py
python3 -O 04-computation/lrc14_fixed_modulus_bank_escape_thm2161.py
```

checks the consequence object on nine hostile prefix banks through `Q=97`,
including exact defect, primitivity, labelled and multiset residue agreement,
every bank certificate, the phase minimum, and its reduced denominator for
the AP-mimic row. It also checks covering, vanishing bank certificates, and
the exact phase minimum for the covering row below. This finite replay is a
regression for the symbolic proof, not the source of the all-`Q` quantifier.

## 3. Scalar certificate failure survives the covering restriction

The AP-mimic row (4) intentionally has no zero residue modulo fourteen. A
second seven-defect row shows that imposing
[THM-523, the constructive q-witness covering reduction](THM-523-lrc14-constructive-q-witness-reduction-to-covering-sets.md),
does not rescue the scalar certificates. Put

```text
W_Q={1,2,3,4,5,6,L,L+8,L+9,...,L+13}.                (20)
```

This again has thirteen distinct positive speeds, is primitive, and replaces
exactly the seven AP entries `7,...,13`. Since `L` is a speed and every
`2<=q<=Q` divides `L`,

```text
m_q(W_Q)=0                         for every 2<=q<=Q. (21)
```

In particular `W_Q` is divisor-complete through fourteen, hence lies in the
stated THM-523 covering class.

The same phase `x` from (10) is an explicit escape. The speeds `1,...,6` are
handled by (12), while `L+j`, `8<=j<=13`, reduce as in (13) to
`||kx||` with `1<=k<=6`. For the remaining speed `L=D-14`, use

```text
14x=7/4-21/(2D).                                      (22)
```

Here `D=L+14>42`, so `3/2<14x<7/4` and therefore

```text
||Lx||=||14x||=2-14x=1/4+21/(2D)>1/4>x.              (23)
```

The speed `1` again attains `x`. Thus

```text
min_(w in W_Q)||wx||=x>3/41,                          (24)
```

and the successful denominator is greater than `Q` by (18)--(19).

## 4. Consequence for the current LRC route

The source row (4) and the tight AP have identical labelled residue words at
every prescribed modulus, but their geometric gaps differ. Hence the map

```text
V -> (V mod q)_(q<=Q)                                 (25)
```

preserves every fixed-bank modulus certificate and destroys the adaptive
scale at which the witness (10) lives. Primitivity, defect seven, the
divisibility pins `7,...,13`, and the mod-seven odd-unit SDR do not repair
that loss: (4) inherits the AP
residue pattern for all fixed moduli and in particular contains a multiple of
each `q=7,...,13`.

The covering row (20) strengthens the fixed-bank warnings in THM-566 and
THM-2072 in the exact shape of the current defect-seven proposal: adding the
covering restriction makes every scalar `m_q/q` certificate in the bank
identically zero and still does not detect the easy phase (10). Unlike (4),
the row (20) is not compared with a bank-identical row of different gap, so
this does **not** prove information-theoretic blindness of every algorithm
using all residues. Nor does it show the defect-seven branch is hard or false.
It shows that the scalar-certificate method must choose a modulus from the
row, use a relation/scale sidecar, or be replaced by richer residue or phase
geometry.

THM-2163 strengthens this loss statement in a different hostile family. Its
two primitive defect-seven rows preserve the complete labelled residue bank
through modulus thirteen, the **entire** coefficient-height-29 relation box,
and the scalar maximum, yet have different mod-17 margins. Thus neither
bounded relations nor scalar magnitude restore the adaptive residue lost by
a fixed bank. THM-2163 identifies a necessary richer state: a radix carry
together with its quotient-owner/termination sidecar. Those coordinates
alone do not yet give a finite complete state space, because the owner mask
records termination but does not bound digit depth.

## 5. The exact residue analogue only halves enumeration

The dyadic checksum of THM-2160 does have one limited residue analogue. Let
`q` be an odd prime, let `g` be a primitive root, and put `H=q-1`. Represent
a multiset of exactly `j` nonzero residues by multiplicities

```text
c_e at g^e,                         e in Z/HZ,
S(c)=sum_e e c_e mod H.                              (26)
```

If `a=nu_2(j)` and `2^(a+1)|H`, put

```text
delta=H/2^(a+1).                                      (27)
```

Global dilation by `g^delta` cyclically shifts the exponent multiplicities,
preserves `j`, and changes the checksum by

```text
j delta=H/2 mod H.                                    (28)
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
