---
id: THM-3659
title: "LRC mod-169 carry-response cochain obstruction"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT; PENDING INDEPENDENT HOSTILE AUDIT.
  In THM-3657's pinned rank-two correction quotient, changing the address law
  from split mod-13 addition to cyclic mod-169 addition has an exact
  carry-response defect.  It vanishes off carries and, on a carry, factors
  through a 156-cell split-sum table with 139 nonzero cells, 55 distinct
  values, and rank two.  Hence a single carry bit is insufficient.  The
  generic 124-address scalar has 16 values but full cyclic Fourier support;
  neither a constant carry vector nor a sparse-frequency repair survives.
  This is a static finite-field obstruction, not a current law or LRC(14).
source: kps-s191 / THM-3657 carry-response continuation, 2026-08-21
depends_on:
  - THM-3657-lrc-two-current-quotient-address-atlas-and-reversal-gate
related:
  - THM-3658-lrc-mod169-carry-fourier-block-intertwiner
  - THM-2334-relation-residue-current-and-character-twist-pushforward
script: 04-computation/lrc_mod169_carry_response_cochain_obstruction_thm3659.py
output: 05-knowledge/results/lrc_mod169_carry_response_cochain_obstruction_thm3659.out
script_sha256: 7617ebcf84508513b71efe3009c84ba88b46d1db41123f1bf1b17530a48aa186
output_sha256: 251d7b82659c70e9c55a44b75b2e394d4bde32865da9be184fb4043504e07023
semantic_sha256: 31a8f50ea2d5a9dbbc52403a50dd49594addc45baa3de9c0043752db64a1fac8
hash_basis: raw LF bytes
---

# THM-3659 -- carry needs an address-dependent rank-two response

**PROVED + FINITE-EXACT + VERIFIED-EXACT; PENDING INDEPENDENT HOSTILE
AUDIT.**  This theorem computes the first nonlinear obstruction behind the
linear change of basis in THM-3658.  The obstruction is much smaller than an
arbitrary function on two addresses, but strictly larger than one carry bit.

## 1. The two coboundaries

Work over the pinned field

```text
F=F_p,                  p=755373809845391722745761.    (1)
```

Let `E` be THM-3657's rank-two correction plane, with its pinned RREF
coordinates, and let

```text
e:F13^2 -> E                                                (2)
```

send a two-current address to its exact correction row.  For addresses
`x=(x0,x1)` and `y=(y0,y1)`, write

```text
x (+) y = ((x0+y0) mod 13, (x1+y1) mod 13),
x (*) y = digits of (x0+13x1+y0+13y1) mod 169.          (3)
```

Thus `(+)` is split addition and `(*)` is cyclic addition under the assembly
`a=x0+13x1`.  Define the two defects and their difference by

```text
Delta_+(x,y)=e(x (+) y)-e(x)-e(y),
Delta_*(x,y)=e(x (*) y)-e(x)-e(y),
Gamma(x,y)=Delta_*(x,y)-Delta_+(x,y).                  (4)
```

On all `169^2=28,561` ordered pairs, the exact records
`(rank, nonzero rows, distinct rows)` are

```text
Delta_+: (2,26474,549),
Delta_*: (2,26079,542),
Gamma:   (2,11856, 55).                                (5)
```

The ranks and distinct-value counts are taken in the pinned coordinate copy
of `E`; they are invariant under any invertible change of basis of `E`.

## 2. Exact factorization through the split sum

Put

```text
t=x (+) y,
kappa(x,y)=1_(x0+y0>=13).                              (6)
```

Elementary base-13 arithmetic gives

```text
x (*) y=(t0,t1+kappa mod 13).                          (7)
```

Substitution in `(4)` therefore proves the identity

```text
Gamma(x,y)
 =e(t0,t1+kappa)-e(t0,t1).                            (8)
```

In particular `Gamma=0` when `kappa=0`.  When `kappa=1`, necessarily
`0<=t0<=11`, and the response factors through the 156-cell table

```text
R(t0,t1)=e(t0,t1+1)-e(t0,t1),
0<=t0<=11, 0<=t1<=12.                                 (9)
```

Exact enumeration of `(9)` gives

```text
rank(R)=2,
139 nonzero cells,
55 distinct values.                                  (10)
```

For a fixed table cell `(t0,t1)`, exactly `13(12-t0)` ordered pairs have
that split sum and a carry.  Consequently

```text
sum_(t0,t1) 13(12-t0) 1_(R(t0,t1)!=0)=11856,          (11)
```

which independently invoices every nonzero entry in the `Gamma` ledger.

Equations `(9)`--`(10)` refute the proposed one-bit model

```text
Gamma(x,y)=kappa(x,y)c                                (12)
```

for any fixed `c in E`: on carry pairs the response has 55 values and spans
both dimensions of `E`.  The strongest surviving compression is the
address-dependent response table `(9)`, not the bit alone.

## 3. The generic line is finite-state but spectrally dense

Assemble addresses by `a=r0+13r1`.  THM-3657 partitions the labels into a
37-address kernel, eight exceptional singleton directions, and a generic
124-address projective line.  In the pinned RREF coordinates, the first
coordinate parametrizes that generic line.  Across the eight generic
intervals

```text
[13,24], [26,47], [49,70], [72,77],
[91,96], [98,119], [121,142], [144,155],              (13)
```

this scalar takes exactly 16 distinct values and is reversal-even under
`a -> 168-a`.  This is a genuine finite-state compression.  Two tempting
stronger simplifications nevertheless fail exactly:

```text
minimal interpolation degrees on the eight intervals
  =(11,21,21,5,5,21,21,11),

cyclic DFT support sizes of the two full coordinate sequences
  =(169,169).                                         (14)
```

Thus the generic response is neither uniformly low-degree on its natural
intervals nor sparse in mod-169 cyclic frequency.  The 16-value automaton is
real, but any useful transport must retain its address partition or an
equivalent sidecar.

## 4. What survives for the LRC program

The failed bit model and full Fourier support rule out two cheap bridges from
THM-3657 to THM-2334.  They do not rule out a weighted or filtered bridge.
The exact residual object is now finite and typed:

```text
the response table R on 12 x 13 cells,
the eight exceptional addresses X,
and the 16-state generic scalar partition.            (15)
```

A valid next gate must show that a physically admissible cancellation cannot
average `(15)` entirely inside the generic reversal pairs.  Arbitrary
reversal-symmetric weights do not have this property, so chronology or a
current-derived positivity/variance law is indispensable.

## 5. Exact verification and boundary

Reproduce with

```bash
python3 -B 04-computation/lrc_mod169_carry_response_cochain_obstruction_thm3659.py
python3 -B -O 04-computation/lrc_mod169_carry_response_cochain_obstruction_thm3659.py
```

The assertion-free companion source-pins THM-3657, reconstructs the raw
two-current tensor, checks all 28,561 ordered pairs, the 156 response cells,
the multiplicity invoice `(11)`, the generic intervals, both full cyclic
DFTs, and every claimed digest.  Normal and optimized runs are byte-identical.
The principal exact digests are

```text
response table: 0e48e8384bcb5d98f52eba69f3e3bf1a01d739558d5bf53d90c94767aa455939,
response image: fab617c98cced507ae0963863b7d2fa3cdb35fb675fb1df6ecc9bad6ea74b309,
generic groups: d7df097a3d6be10d04b00462805c40759c553255ffe9b13c8bc5c7cdf93bdee0,
semantic:       31a8f50ea2d5a9dbbc52403a50dd49594addc45baa3de9c0043752db64a1fac8. (16)
```

This is a static finite-field response-cochain theorem.  It proves no
physical current law, chronology, admissible coefficient constraint,
exceptional-address entry, characteristic-zero transport, or LRC(14).
**QED.**
