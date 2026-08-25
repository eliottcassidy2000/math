---
id: THM-2066
title: "Dyadic-seam owner-word CRT atlas"
status: >
  PROVED. For every
  safe quotient-core clock, an odd tail determines a binary word recording
  which of the two dyadic lifts it kills. A strict counterexample requires two
  everywhere-eligible words that are bitwise complementary. This gives an
  exact finite residue-pair obstruction modulo 2N, composes by generalized
  CRT over clock banks, and closes every divisor-complete primitive quotient
  eleven-core with max(C)<=24. This is not LRC(14).
source: codex-2026-07-21-LRC-dyadic-owner-words
depends_on:
  - THM-2061
related:
  - THM-2059
  - THM-2060
  - THM-2064
  - HYP-8871
script: 04-computation/lrc14_dyadic_owner_word_crt_referee_codex_20260721.py
output: 05-knowledge/results/lrc14_dyadic_owner_word_crt_referee_codex_20260721.out
script_sha256: 0678aaf536c60279c9989fd86cb6cc5d02ad0fd983429dc83daefa2157ee98ad
output_sha256: 5d6fba079c85c532184762bcd626f78915df219f4f2b5e70a49e6ec5447c3b85
hash_basis: normalized repository blobs (LF)
---

# THM-2066 -- dyadic-seam owner-word CRT atlas

Let `C` be an eleven-speed quotient core and retain the THM-2061 seam

```text
S=2C union {x,y},       x,y distinct positive odd.       (1)
```

Write `|a|_N=min_(k in Z)|a-kN|` and label each class modulo `N` by its
canonical representative `0<=r<N`. For a clock `N`, define its weak-safe
labelled quotient-core packet

```text
A_N(C)={0<=r<N:14|cr|_N>=N for every c in C}.            (2)
```

## 1. Exact owner words

For an odd residue `z mod 2N`, call it **eligible** when

```text
7|zr|_N<N       for every r in A_N(C).                   (3)
```

The inequality makes the nearest integer to `zr/N` unique. Define the binary
owner word on the labelled packet (2) by

```text
omega_(N,C,z)(r)=nint(zr/N) mod 2.                       (4)
```

Here `nint` is the unique nearest integer. The canonical representative is a
label gauge: replacing `r` by `r+N` flips every odd tail's owner bit, while
the complementary-word relation below is unchanged by that common flip.

Both eligibility and (4) depend only on `z mod 2N`: adding `2N` changes the
nearest integer by the even number `2r`.

Let `E_N(C)` be the eligible odd residue classes and put

```text
R_N(C)={(u,v) in E_N(C)^2:
          omega_v(r)=1-omega_u(r) for every r in A_N(C)}. (5)
```

Then `(u,v)` belongs to (5) if and only if, over every safe quotient phase
`r/N`, the two odd tails together kill both dyadic lifts

```text
t_j=(r/N+j)/2,       j=0,1,                              (6)
```

strictly below `1/14`.

### Proof

THM-2061 gives the exact lift rule: an odd tail `z` is dangerous on one lift
over `r/N` precisely when (3) holds at `r`, and the killed lift is
`j=nint(zr/N) mod 2`. It kills at most one of the two lifts. Hence two tails
kill both exactly when both are eligible and their owner bits are opposite at
every packet point. This is (5). QED.

In particular,

```text
R_N(C)=empty                                               (7)
```

is an exact certificate that no strict counterexample (1) exists over `C`.
The packet must be retained with its labels. The number of eligible tails
alone loses which words are complementary.

## 2. Clock banks and divisor transport

For a finite clock bank `F`, put

```text
L=lcm_(N in F)(2N),
R_F(C)={(u,v) mod L:
          (u mod 2N,v mod 2N) in R_N(C) for every N in F}. (8)
```

This is a finite generalized-CRT join of residue-pair constraints. Every
strict counterexample has `(x mod L,y mod L)` in (8). Thus an empty bank atlas
closes the seam over `C`; a nonempty atlas is only a necessary periodic
sidecar and must still meet THM-2061's metric box.

If `N|M`, the reduction map sends

```text
R_M(C) -> R_N(C).                                        (9)
```

Indeed `r/N` is the same phase as `(rM/N)/M`, and
`r in A_N(C)` implies `rM/N in A_M(C)`. Tail eligibility and nearest-integer
parity at that phase are unchanged by reduction modulo `2N`. Consequently a
divisibility-maximal clock dominates its divisors in a bank. This is a useful
compression, not a claim that unrelated clocks are independent.

## 3. Exact closure through quotient-core maximum 24

There is no strict LRC(14) counterexample (1) with

```text
C subset {1,...,24},       |C|=11,       gcd(C)=1.        (10)
```

THM-2061 already forces every surviving quotient core to be divisor-complete:

```text
for each n=2,...,14, some c in C is divisible by n.       (11)
```

There are exactly

```text
2,496,144 total eleven-subsets of {1,...,24},
59,880 primitive divisor-complete candidates.             (12)
```

For each candidate, the exact clock bank

```text
F={15,16,...,34}                                         (13)
```

contains an `N` with `R_N(C)=empty`. The certificate-clock histogram is

```text
16:3502, 17:10779, 18:1131, 19:16079, 20:1465,
21:7885, 22:2947, 23:10082, 24:422, 25:4597,
26:733, 27:196, 28:28, 29:20, 31:10, 32:3, 34:1.         (14)
```

The counts sum to `59,880`; clocks `15,30,33` are unused. Equation (7) then
closes every candidate, proving (10). This strictly extends THM-2061's exact
`max(C)<=19` census while using only finite clock words, not interval
subdivision or an a priori tail bound.

The qualifier `gcd(C)=1` is part of the proved normalized quotient-core lane,
not cosmetic. No claim is made here for a seam presentation with a
nonprimitive quotient core.

## 4. Relation to the CRT packet carrier

THM-2059 joins one core packet to one tail by a nonnegative histogram dot
product. On the dyadic two-tail seam each tail saturates exactly one of two
sheets, so nonemptiness is no longer the right statistic: the obstruction is
**complementary ownership over every core-safe residue**. Formula (5) is the
smallest lossless sidecar. It is the signed two-sheet refinement anticipated
by THM-2060/2061 and HYP-8871; it is not a Fourier or modular-cusp identity.

## 5. Assumption challenge and tournament analysis

The vertices are owner words, not runners or tail magnitudes. The intrinsic
binary relation is symmetric complementarity `w~wbar`; it is a matching with
ties from multiple residue classes carrying the same word. Artificially
orienting it would lose exactly the pair needed in (5), so a tournament is
not the proof carrier. For search scheduling one may order words by eligible
class multiplicity and use lexicographic word order as a tie Hamiltonian path;
score histograms then rank common words but do not decide whether complements
exist. The faithful object is the complement graph together with its residue
classes modulo `2N` and the CRT sidecar (8).

## 6. Exact referee

The companion uses integer bitsets. It checks the owner-word criterion against
direct enumeration of both lifts, verifies divisor transport, exhausts all
`2,496,144` cores in (12), independently applies the divisor pins, and checks
the histogram (14). Runtime checks survive optimized Python and the frozen
output ends in `PASS`. QED.
