---
id: THM-2261
title: "Expiration-image surjectivity and one-core carrier no-go"
status: >
  PROVED + VERIFIED-EXACT. For T(x)=13x mod 1, every thirteen-unit guard H,
  every thirteen-unit owner core u, and every lambda>=1, the raw
  guard-plus-owner set has completely unrestricted expiration image:
  T^(lambda+1)(C_H intersect D_(13^lambda u))=R/Z. The proof is an exact
  fibre-restriction identity plus the facts that every 13-root fibre meets
  C_H and D_u. Moreover one rational point simultaneously realizes the
  failure on the local data of every one of the 150 strict remaining
  profiles: it is exclusively owned by the depth-one blocker, avoids five
  distinct unit masks, and its post-expiration one-core word is 00000...
  The exact one-core two-window carrier has mass 6055/28561, just below
  THM-2255's strict expiration floor, but it is not a valid target inclusion
  from raw guard/owner data. No global scalar cover, profile exclusion, or
  LRC(14) proof is asserted.
source: codex-2026-07-25-expiration-surjectivity
depends_on:
  - THM-2198-scalar-five-plus-three-image-pump-and-first-depth-exclusion
related:
  - THM-2234-first-depth-one-private-two-owner-mass-and-two-step-expansion
  - THM-2255-valuation-separated-pair-cap-and-exclusive-owner-mass
script: 04-computation/lrc14_expiration_surjectivity_one_core_no_go_thm2261.py
output: 05-knowledge/results/lrc14_expiration_surjectivity_one_core_no_go_thm2261.out
script_sha256: 788cda64a993ecc9a4c67a97674ee21b0d1213ee93848d5b1786868b95d69d99
output_sha256: 93411a4183d2f3f8f9ca7ab6ebdd10ddc0db4f49f9a09edfddf11f12dadbd8d4
hash_basis: working-tree bytes (LF)
---

# THM-2261 -- expiration forgets the raw owner completely

Put

```text
T(x)=13x mod 1,
D_s={x in R/Z:||sx||<1/14},
C_H={x in R/Z:||Hx||>1/7}.                           (1)
```

Assume throughout that `13` divides neither `H` nor `u`. In the scalar
branch `H` is also odd, but parity is not needed here.

For every integer `lambda>=1`,

```text
T^lambda(C_H intersection D_(13^lambda u))=D_u,      (2)

T^(lambda+1)(C_H intersection D_(13^lambda u))
 =R/Z.                                               (3)
```

Thus a point known only to be guard-safe and owned by a depth-`lambda`
blocker can have an arbitrary post-expiration image. This identifies the
first failed implication behind a tempting near-equality after THM-2255.

## 1. Every relevant root fibre is nonempty

Fix `y in R/Z`. The thirteen roots of `y` under `T` form a translate of the
uniform thirteen-grid. Multiplication by a thirteen-unit permutes that grid.

The danger arc has length

```text
1/7>1/13.                                            (4)
```

Every translate of the thirteen-grid therefore meets `D_u`. In fact its
danger occupancy is one or two away from the finitely many boundary
phases. In particular,

```text
T(D_u)=R/Z.                                          (5)
```

The complement of `C_H` is the closed central arc of length `2/7`. Since

```text
2/7<4/13,                                            (6)
```

that arc contains at most four points of any translated thirteen-grid.
Every fibre consequently contains at least nine points of `C_H`, including
at boundary phases. Hence

```text
T(C_H)=R/Z,
T^lambda(C_H)=R/Z                 for lambda>=1.     (7)
```

This root argument is pointwise; no null-set convention is needed for
(5)--(7).

## 2. Restrict the fibre before taking its image

The owner identity is exact:

```text
D_(13^lambda u)=T^(-lambda)(D_u).                    (8)
```

For any map `f` and sets `A,B`,

```text
f(A intersection f^(-1)(B))=f(A) intersection B.    (9)
```

Apply (9) with `f=T^lambda`, `A=C_H`, and `B=D_u`.
Equations (7)--(8) give

```text
T^lambda(C_H intersection D_(13^lambda u))
 =T^lambda(C_H) intersection D_u
 =D_u,                                               (10)
```

which is (2). Applying `T` and using (5) proves (3).

The order in (9) is the essential point. The owner bit does persist exactly
through `lambda` images, but its next image is unrestricted because every
future point has a danger-root predecessor. “Owner expiration” is literal
loss of this coordinate, not merely a weakening of its density.

## 3. Exact local witness for every strict remaining row

Surjectivity by itself leaves open whether the five unit masks and two other
blocker labels might create a pointwise restriction. They do not do so from
their local membership bits alone.

For every pair

```text
2<=b<c,                                              (11)
```

take

```text
H=1,
x=79/338,

q_i=1+338i,                         1<=i<=5,

c_1=13,             c_2=13^b,      c_3=13^c.        (12)
```

The five `q_i` are distinct thirteen-units and the three blockers are
distinct. Exact arithmetic gives

```text
1/7<x<1/2,

q_i x=x+79i                         mod 1,

13x=3+1/26,

13^b x=79*13^(b-2)/2=1/2            mod 1,
13^c x=79*13^(c-2)/2=1/2            mod 1.           (13)
```

Consequently

```text
x in C_1,
x notin union_(i=1)^5 D_(q_i),
x in D_(c_1),
x notin D_(c_2) union D_(c_3).                       (14)
```

In the notation of THM-2255, `x` lies in the locally defined exclusive
owner stratum `E_1`.

But the depth-one owner expires to

```text
y=T^2(x)=1/2.                                        (15)
```

Since `13` is odd,

```text
T^k(y)=1/2,
1_(D_1)(T^k(y))=0                   for every k>=0.  (16)
```

Thus the post-expiration owner-core word is the persistent fixed word

```text
00000... .                                           (17)
```

Restricting (11) to

```text
5<=c<=19,             2<=b<c                         (18)
```

gives all

```text
sum_(c=5)^19(c-2)=150                                (19)
```

strict profiles in the current scalar ledger. The witness is deliberately a
**local hostile control**. It does not assert that the packet in (12)
satisfies the global scalar cover, and hence it is not a counterexample to
LRC(14). It proves that a future-word restriction cannot be deduced from
the displayed guard, unit-mask, and exclusive-owner membership bits alone.

## 4. The near-equality and why it does not close a row

For one fixed thirteen-unit core `u`, define the stationary danger process

```text
X_k(y)=1_(D_u)(T^k y),               k>=0,            (20)
```

with `y` Haar-uniform. The exact thirteen-root law is

```text
P=[[11,2],
   [12,1]]/13,                                         (21)
```

with stationary distribution

```text
pi=(6/7,1/7).                                         (22)
```

The chain is reversible. Define the two-window carrier

```text
K_u={y:
  (X_0 or X_1 or X_2)
  and
  (X_2 or X_3 or X_4)}.                              (23)
```

Its complement is the union of the two events `X_0X_1X_2=000` and
`X_2X_3X_4=000`. Their intersection is the all-zero word of length five.
Therefore

```text
measure(K_u)
 =1-2(6/7)(11/13)^2+(6/7)(11/13)^4
 =6055/28561
 =0.212002380869... .                                (24)
```

THM-2255 gives some labelled exclusive-owner image in every strict row the
slightly larger floor

```text
88159/415800
 =0.212022607023... ,

88159/415800-6055/28561
 =240199/11875663800
 =0.0000202261535...>0.                              (25)
```

If that same labelled image were contained in `K_u`, (24)--(25) would close
all strict rows. Equations (3) and (15)--(17) show why no such containment
comes from raw expiration transport: the unrestricted image can hit the
permanent all-zero fixed point, which lies outside `K_1`.

The strongest surviving target is narrower and genuinely global. A valid
consumer would have to use the scalar cover across phases to manufacture a
named post-expiration set from the five unit masks and the two other blocker
flows, while retaining which roots came from the marked exclusive owner.
Neither the one-core Markov law nor the comparison of the two masses supplies
that map.

## 5. Scope and exact verification

The theorem proves:

```text
raw guard-owner expiration image = whole circle;
one strict-row local point has future word 00000...;
the 6055/28561 carrier is numerically insufficient as an unproved target.
```

It does not prove:

```text
failure of a target inclusion under the full scalar-cover hypothesis;
existence of a scalar counterexample;
exclusion of any of the 165 profiles;
LRC(14).
```

The companion checks the Markov law and stationarity, enumerates all
thirty-two length-five words, verifies (24)--(25), and verifies the rational
witness for all `150` strict rows. Reproduce with

```bash
python3 04-computation/lrc14_expiration_surjectivity_one_core_no_go_thm2261.py
python3 -O 04-computation/lrc14_expiration_surjectivity_one_core_no_go_thm2261.py
```

Both modes produce the stored transcript byte for byte. QED.
