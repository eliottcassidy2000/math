---
id: THM-3210
title: "Reflected midratio primitive floor, one-sided tooth drift, and gcd-six closure"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT.  Inside THM-2941's 561-body
  reflected sufficient-family residual, every physical extreme pair with
  3<Q/m<6, Q/m not in {4,5}, and gcd(m,Q)>=6 has a positive pair-overlap
  certificate.  The ratio-restricted primitive fibre has sharp floor 1/77,
  with equality only at (P,R)=(2,11),(3,11).  A one-sided boundary-pairing
  lemma replaces the inherited transport loss 4(a+b)/(gL) by
  a(P+1)/(gLP-a)+b(R+1)/(gLR-b).  Fourier monotonicity closes P>=13 from
  one 16,830-row base.  The remaining 135 primitive channels give 2,272,050
  exact analytic invoices; independent full-tooth calculations close all
  150 nonpositive analytic rows on located body-safe cells.  Thus the
  reflected certificate-failure residual has gcd(m,Q)<=5.  This is not a
  physical-survivor census, a proof of LRC(14), or an independent theorem
  about arbitrary reflected packets outside the inherited sufficient family.
audit: >
  Normal and optimized Python outputs byte-match.  The verifier checks the
  complete fourteen-grid primitive breakpoints, the Fourier handoff, all
  7,644 hostile boundary-displacement rows, the complete finite channel/body/
  orientation universe, and an integer two-pointer engine against the
  promoted Fraction interval engine on all 14,304 cells used by the 150
  located repairs.
source: root/lrc-midratio-2026-08-02
depends_on:
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
related:
  - THM-3135-directed-cycle-weak-order-lane-cover-and-reflected-h-boundary
script: 04-computation/lrc14_reflected_midratio_primitive_floor_g6_closure_thm3210.py
output: 05-knowledge/results/lrc14_reflected_midratio_primitive_floor_g6_closure_thm3210.out
script_sha256: 207cdd0194967c2941b2d1eb35de123c05b9619edf0e25213f369bb41b9dcd08
output_sha256: a462be53c8c7278b3d2687fe348a70a7aa6b7a33ebf35575669dae126e8c4779
hash_basis: LF-normalized bytes
---

# THM-3210 -- reflected midratio primitive floor and gcd-six closure

**PROVED + FINITE-EXACT + VERIFIED-EXACT.**

This theorem improves one precise sufficient-family frontier in THM-2941.  It
does not assert that a row left by that sufficient certificate is a physical
Lonely Runner survivor.

## 1. Inherited universe and physical coordinates

Let

```text
E subset {1,...,14},  |E|=6,       L=14 lcm(E).
```

THM-2941 leaves `561` body sets after its arbitrary-level reflected atlas.
The two bodies on which its same-level good-pair graph is not complete are
disjoint from these `561`.  Hence any word in this residual which is not
already certified has six distinct positive levels.

Fix such a word.  Let `m` and `Q` be its global minimum and maximum levels,
and let `a,b in E` be the labels occupying those physical slots.  Write

```text
m=gP,        Q=gR,        gcd(P,R)=1.                       (1)
```

The live midratio wedge is

```text
3P<R<6P,       R/P not in {4,5}.                            (2)
```

The reduced rays of ratio four and five are `(P,R)=(1,4),(1,5)`, already
closed by THM-2941.  Thus every channel in `(2)` has

```text
P>=2,       R>=7.                                          (3)
```

This typing is important: the earlier generic tail used only `P>=1,R>=7`
and therefore charged both transport and singleton debt at a nonexistent
midratio corner.

On a body-safe integer cell, each reflected clause at level `q` has exactly
`q` full teeth.  THM-2941 gives the exact singleton excess

```text
epsilon(E,q_vec)=sum_(e in E) e/[7(q_e L-e)].              (4)
```

With the physical extremes fixed, the largest possible `(4)` among
six-distinct-level words assigns the four remaining levels

```text
m+1,m+2,m+3,m+4
```

to the remaining labels in decreasing-label order.  Indeed the summand
`e/(qL-e)` increases with `e`, decreases with `q`, and its two-by-two
rearrangement gap is strictly positive when both inequalities are strict.
Call this exact worst debt

```text
D_(E,a,b)(m,Q).                                            (5)
```

Therefore a physical extreme-pair overlap larger than `(5)` closes the cell
by the inherited Bonferroni certificate.

## 2. The ratio-restricted primitive floor

For coprime `P<R`, put

```text
T_s(z)=sum_(n in Z) (s-|z+n|)_+,

F_(P,R)(z)=
 [T_((P+R)/14)(z)-T_((R-P)/14)(z)]/(PR).                  (6)
```

This is the exact fully periodized primitive overlap fibre.  The Fourier
estimate already proved in THM-2941 is

```text
F_(P,R)(z) >= 1/49-1/(2PR).                               (7)
```

If `PR>=68`, the right side of `(7)` is at least

```text
87/6664 > 1/77.                                           (8)
```

It remains to inspect `PR<68`.  Every breakpoint of `(6)` lies at `k/14`
modulo one, so the fourteen-grid is complete, not a sampling proxy.  Under
`P>=2` and `3P<R<6P`, there are exactly eleven small-product channels:

| `(P,R)` | `min_z F_(P,R)` |
|:---:|:---:|
| `(2,11)`, `(3,11)` | `1/77` |
| `(3,10)` | `1/70` |
| `(2,9)` | `1/63` |
| `(3,17)` | `2/119` |
| `(3,16)` | `1/56` |
| `(3,13)` | `5/273` |
| `(4,15)` | `2/105` |
| `(4,13)` | `1/52` |
| `(2,7)`, `(3,14)` | `1/49` |

The table contains precisely eleven channels.  Combining that complete bank
with `(8)` proves the sharp restricted statement

```text
min_z F_(P,R)(z) >= 1/77,                                 (9)
```

with equality exactly for

```text
(P,R)=(2,11),(3,11).                                      (10)
```

The global THM-2941 floor `1/105` was sharp at `(3,5)`.  That channel is not
a physical extreme pair in `(2)`, which is why `(9)` is a genuine gain rather
than a new estimate on the old full channel universe.

## 3. One-sided boundary displacement

The second gain keeps the sign of the reflected phase perturbation.  Let
`N` be a positive integer, `0<c<gLN`, `delta=1/14`, and

```text
epsilon=c/(gL),
A_epsilon={x in [0,1]: ||(N-epsilon)x-theta||<delta},
A_0      ={x in [0,1]: ||Nx-theta||<delta}.               (11)
```

The two boundary families have numerators

```text
beta=theta+k+delta,       beta=theta+k-delta.
```

A boundary relevant to `[0,1]` moves monotonically right from

```text
beta/N       to       beta/(N-epsilon).                    (12)
```

If the second root passes one, clipping only shortens its contribution.
For either residue-one arithmetic progression, the sum of its members in
`[0,N]` is at most `N(N+1)/2`.  Pairing all boundaries and summing `(12)`
therefore gives

```text
mu(A_epsilon symmetric_difference A_0)
 <= epsilon/[N(N-epsilon)] * N(N+1)
 = c(N+1)/(gLN-c).                                        (13)
```

This is an all-real-phase proof.  The verifier's `7,644` fourteen-grid rows
are hostile controls for the denominator and endpoint clipping, not the
source of its quantifier.

Now subdivide one body-safe cell by `u=(r+x)/g`.  The two exact phases are

```text
P x-a(gj+r+x)/(gL),       R x-b(gj+r+x)/(gL).             (14)
```

Dropping only the final `x` in the two numerators produces primitive clauses
of slopes `P,R`.  On each subcell their intersection is at least the global
primitive floor.  Apply `(13)` to the two clauses, use the union bound for
the lost intersection, and average over the `g` subcells.  The physical
overlap is at least

```text
min_z F_(P,R)(z)
 - a(P+1)/(gLP-a)
 - b(R+1)/(gLR-b).                                       (15)
```

Every denominator in `(15)` is positive by `(3)`, `g>=1`, `L>=168`, and
`a,b<=14`.  The inherited `4(a+b)/(gL)` bound treated the drift as two-sided;
`(13)` is the missing one-sided coordinate.

## 4. The infinite primitive tail

Use `(7)` in `(15)` and subtract `(5)`.  For fixed body and orientation:

- the Fourier floor increases with `PR`;
- `(N+1)/(gLN-c)` strictly decreases with `N` and with `g`;
- every singleton-debt denominator strictly increases with its level.

If `P>=13`, then `R>=3P+1>=40`.  Consequently every such invoice for
`g>=6` is bounded below by the single corner

```text
(g,P,R)=(6,13,40).                                        (16)
```

The complete `561*30=16,830` exact body/orientation census at `(16)` is
strictly positive.  Its weakest margin is

```text
26392028086719269675209
--------------------------------------------  > 0,        (17)
104017741322142148178386800
```

on `E=(1,2,3,4,6,12)`, orientation `(5,4)`.  At this row the Fourier floor,
one-sided loss, and exact debt are respectively

```text
991/50960,
138797/7330429,
18470013214620220016/71440756402570156715925.
```

Thus `(16)` and monotonicity close every `P>=13` channel.

## 5. Complete finite primitive bank

For `2<=P<=12`, enumerate every coprime integer `R` with `3P<R<6P`.  The
per-`P` channel counts are

```text
3,6,6,12,6,18,12,18,12,30,12,
```

whose sum is `135`.  Crossing with `561` bodies and `30` ordered physical
extreme slots gives exactly

```text
135*561*30 = 2,272,050                                   (18)
```

analytic invoices.  Each uses the exact fourteen-grid minimum of `(6)`, the
loss `(15)`, and the debt `(5)` at `g=6`.  Exactly `150` invoices are
nonpositive.  This is failure of the uniform analytic lower bound, not of the
physical pair.

For each of those `150`, the verifier enumerates every body-safe cell and
computes the exact physical overlap by an independent integer two-pointer
intersection of the `gP` and `gR` full teeth.  It compares that result with
the promoted rational interval engine on all `14,304` tested cells.  Every
comparison agrees and every row has a positive located cell.  The weakest
located margin is

```text
12808546596795938/662437112164402113 > 0                 (19)
```

on `E=(1,2,3,4,6,12)`, orientation `(3,5)`, channel
`(P,R)=(2,7)`, cell `152`.  There

```text
overlap = 6084/295261,
debt    = 841300116478834/662437112164402113.
```

This completes every finite channel at `g=6`; the same monotonicity in `g`
closes all larger scales.

## 6. Assembly and the directed-cycle boundary

Combining Sections 4 and 5 proves

```text
PROVED: inside THM-2941's reflected sufficient family, if
        3<Q/m<6, Q/m not in {4,5}, and gcd(m,Q)>=6,
        then the physical extreme pair closes the packet.                 (20)
```

Thus the current reflected sufficient-certificate obligation contracts to

```text
561 bodies;  m>=2;  3<Q/m<6;  Q/m not in {4,5};
gcd(m,Q)<=5.                                               (21)
```

THM-3135 proves that the older standard uniform single-pair comparison lanes
on the hostile body `H` form a DAG.  There is no conflict: `(20)` selects the
intrinsic physical minimum/maximum slots, retains their actual primitive
ratio, and, only where a uniform invoice fails, retains the body-safe cell.
It is a located overlap certificate, not a union of weak-order lanes which
must cover arbitrary assignments by a directed cycle.  The successful
sidecar is therefore `primitive ratio + physical extremes + cell`, precisely
the information discarded by the DAG envelope.

## 7. Verification and scope

Run

```bash
python3 04-computation/lrc14_reflected_midratio_primitive_floor_g6_closure_thm3210.py
python3 -O 04-computation/lrc14_reflected_midratio_primitive_floor_g6_closure_thm3210.py
```

Both outputs byte-match the stored transcript.  All truth-bearing gates use
explicit exceptions rather than `assert`, so optimized mode checks the same
universe.

Statement `(20)` is a theorem only inside the reflected sufficient family
and inherited `561` bodies of THM-2941.  Statement `(21)` is a residual of
that proof mechanism.  Neither statement enumerates physical LRC survivors,
closes other six-drift sectors, or proves LRC(14).
