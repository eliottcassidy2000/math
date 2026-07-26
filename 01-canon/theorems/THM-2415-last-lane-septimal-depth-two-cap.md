---
id: THM-2415
title: "Last-lane septimal depth-two cap"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. Under
  THM-2391's primitive final-lane scalar-cover hypotheses, if q_* has
  septimal depth M>=1 and c_3 has larger depth, then M<=2. The proof is
  uniform in the clean-hole mass: after the quotient-safe transition
  q_* safe -> 13q_* dangerous outside c_3, a positive-measure
  7^(M-1)-invariant set remains. For M>=3 every such orbit reaches the
  thirteenth phase digit on which the H and 13H guard words coincide.
  THM-2391's physical blocker pullback then doubles both guard roots,
  contradicting the unique-double W8 word. This removes every M>=3
  packet, including the positive-clean packets not addressed by the
  no-clean cage chain. The phase mechanism is sharp at M=2 and does not
  exclude M=1 or M=2, close a ledger row by itself, or prove LRC(14).
source: codex-2026-07-26-last-lane-depth-cap
depends_on:
  - THM-2391-blocker-caged-septimal-single-layer-address-reduction
related:
  - THM-2390-septimal-layer-kraft-peeling-and-heavy-word-reduction
  - THM-2392-clean-toothpick-or-bounded-cross-ancestor-cage
  - THM-2393-c3-safe-double-fibre-capacity-and-common-core-residual
  - THM-2405-two-level-septimal-sheet-independence-and-middle-depth-two-cage-elimination
  - MISTAKE-264
script: 04-computation/lrc14_last_lane_septimal_depth_cap_thm2415.py
output: 05-knowledge/results/lrc14_last_lane_septimal_depth_cap_thm2415.out
script_sha256: 63b9494526ea324793225b81c18f26a9002695edd39e7d51e47a5307d477a4e1
output_sha256: ad8634d146be16a53068d5828397418bd9c9e98bddfed212397495932635216b
hash_basis: working-tree bytes (LF)
---

# THM-2415 -- the final septimal lane has depth at most two

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2391 leaves a binary address at depth one and four possible lower
slopes at larger depth. The no-clean descendants THM-2392--2396 and
THM-2405 go much further inside their cage, but they do not make the
depth bound uniform in positive clean-hole mass. The missing operation is
to compare one packet before and after multiplication by thirteen.

That comparison gives a short global cap:

```text
primitive last-lane scalar cover
  + q_* safe but 13q_* dangerous
  + c_3 safe
  + a complete 7^(M-1)-phase orbit
  -> W_H=W_(13H) on one thirteenth digit
  -> both guard roots are doubled
  -> contradiction to the unique-double W8 word.                    (1)
```

The point of the proof is not merely that the transition set is
nonempty. It has a uniform positive measure, so all strict endpoints and
all almost-everywhere exceptional sets can be deleted before choosing
the phase.

## 1. Statement and divided base variables

Retain all hypotheses and notation of THM-2391. In particular there is
an almost-everywhere scalar cover

```text
C_H subset
  union_(i=1)^5 D_(q_i)
  union D_(c_1) union D_(c_2) union D_(c_3),                         (2)

D_v={x in R/Z:||vx||<1/14},
E_v={x in R/Z:||vx||<1/7},
```

with `c_j=13C_j`. There is one top label `q_*` such that

```text
nu_7(q_*)=M>=1,                 nu_7(c_3)>M.                          (3)
```

THM-2391 proves, using primitivity, that

```text
nu_7(H)=nu_7(q_i)=nu_7(c_1)=nu_7(c_2)=0
                                      for every q_i!=q_*.             (4)
```

We prove:

> **Uniform depth cap.** Every packet satisfying (2)--(4) has
> `M<=2`.

Set

```text
Q=q_*/7,                         D=c_3/7.                             (5)
```

For the seven roots

```text
x_j=(y+j)/7,                     j in Z/7Z,                           (6)
```

the correct divided base identities are

```text
1_(D_(q_*))(x_j)    =1_(D_Q)(y),
1_(D_(13q_*))(x_j)  =1_(D_(13Q))(y),
1_(D_(c_3))(x_j)    =1_(D_D)(y).                                    (7)
```

This is the base/root typing required by MISTAKE-264. In particular the
relevant base transition set is

```text
S=(D_(13Q) minus D_Q) intersection D_D^c.                            (8)
```

Every `y in S` makes `q_*` and `c_3` safe at all seven roots while
making `13q_*` dangerous at all seven roots.

## 2. The transition set has positive measure

Assume toward a contradiction that `M>=3`. Write

```text
Q=7^a u,
D=7^(a+d) w,

a=M-1>=2,             d>=1,             7 does not divide u*w.       (9)
```

Multiplication by `7^a` shows that (8) is the pullback of

```text
A intersection B,

A=D_(13u) minus D_u,
B=D_(7^d w)^c.                                                       (10)
```

First,

```text
mu(D_u intersection D_(13u))=1/91,
mu(A)=1/7-1/91=12/91.                                                (11)
```

Indeed multiplication by the seven-unit `u` preserves Haar measure.
Inside `D_1=(-1/14,1/14)`, the two neighboring components of `D_13`
begin exactly at the excluded endpoints `+/-1/14`; only the central
component of length `1/91` remains.

Disintegrate (10) over

```text
t_j=(z+j)/7.                                                         (12)
```

Away from finitely many endpoints, each of `D_u` and `D_(13u)` occupies
one root. Hence `A` occupies either zero roots or one root. By Haar
invariance,

```text
integral sum_j 1_A(t_j) dz=7mu(A)=12/13.                             (13)
```

Therefore the set of bases on which `A` occupies one root has measure
`12/13`.

The high mask `B` is constant on every fibre (12), and its safe base set
has measure `6/7`. Inclusion--exclusion now leaves base mass at least

```text
12/13+6/7-1=71/91.                                                   (14)
```

On each base counted in (14), exactly one of the seven roots belongs to
`A intersection B`. A second use of (13) gives the uniform bound

```text
mu(S)=mu(A intersection B)>=71/(91*7)=71/637>0.                      (15)
```

No relation between the units `u` and `w` was used.

## 3. A complete orbit reaches the collision digit

The set `S` is invariant under translation by `1/7^a`. Put

```text
N=7^a>=49,
I={y:6/13<{Hy}<7/13}.                                                (16)
```

Since `H` is a seven-unit, the `H`-phases on every `N`-point orbit of
`S` form a complete shifted `N`-grid. Every interval of length `1/13`
contains at least `floor(N/13)` points of such a grid. Consequently

```text
mu(S intersection I)
  >=floor(N/13)/N * mu(S)
  >0.                                                               (17)
```

Delete from (17):

- the finite comb endpoint sets;
- the pullbacks at all seven roots of the null exceptional set in the
  scalar cover (2); and
- the pullbacks at all seven roots of the null exceptional set in
  THM-2391's physical identity.

The remainder is still nonempty. Fix a base `y` in it.

## 4. The two guard words coincide

Let `r={Hy}`. Translate and relabel the seven roots by the `H`-gauge,
so that

```text
H x_j=(r+k)/7 mod 1,               k in Z/7Z.                         (18)
```

For generic `r`, the two guard roots are

```text
W_H={k:||(r+k)/7||<1/7}={0,6}.                                     (19)
```

Put `b=floor(13r)`. Since `13=-1 mod 7`, the drifted guard word is

```text
W_(13H)
 ={k:||13(r+k)/7||<1/7}
 ={b,b+1} mod 7.                                                     (20)
```

Our phase choice `r in (6/13,7/13)` gives `b=6`, and hence

```text
W_(13H)={6,0}=W_H.                                                   (21)
```

The integer digit `b` is essential. A bare reflection without this
affine digit is not the correct thirteen transport.

## 5. Two doubles contradict W8

Because `y notin D_Q union D_D`, (7) makes the top `q_*` and high
blocker `c_3` safe at every root. By (4), the lower layer on this
seven-root fibre consists of:

```text
the width-two guard E_H,
four singleton q_i words (q_i!=q_*),
two singleton physical blocker words c_1,c_2.                        (22)
```

Its total incidence is eight. The scalar cover (2) covers all seven
roots. Thus, exactly as in THM-2390--2391's W8 conclusion, precisely one
root can be doubled.

But `y in D_(13Q) intersection D_D^c`, so (7) puts every root in
`D_(13q_*) intersection D_(c_3)^c`. THM-2391's physical pullback
identity applies pointwise:

```text
1_(D_(c_1))+1_(D_(c_2))=1_(E_(13H)).                                (23)
```

Equations (21) and (23) say that the two blockers partition the same
two roots already occupied by the guard `E_H`. Each of those two
distinct roots therefore has multiplicity at least two. This gives at
least two doubled roots, contradicting (22). The assumption `M>=3` is
impossible, proving `M<=2`.

## 6. Sharp boundary of this mechanism

The orbit-hitting step really begins at `M=3`. At `M=2`, take

```text
Q=7,                 D=49,                 H=1,
y=1/91.                                                            (24)
```

Then

```text
y in (D_(13Q) minus D_Q) intersection D_D^c,                         (25)
```

but the seven-point orbit `y+j/7` misses the open phase chamber
`(6/13,7/13)`. Thus no version of the same unconditional grid argument
can remove `M=2`. This is a stopping control for the proof mechanism,
not a scalar-cover example.

## 7. Scope and interaction with the current cage chain

THM-2392, THM-2393, and THM-2405 classify or eliminate deep packets
inside the no-clean cage. The present theorem is orthogonal: (15) and
(17) never assume zero clean-hole mass. It therefore removes all
positive-clean `M>=3` packets as well.

The conclusion is only the uniform depth bound

```text
M in {1,2}.                                                          (26)
```

It does not exclude either remaining depth, produce an all-coordinate
target address, decrement the 165-row ledger by itself, or prove
LRC(14).

## 8. Exact companion and independent audit

Run:

```text
python3 04-computation/lrc14_last_lane_septimal_depth_cap_thm2415.py
```

The standard-library `Fraction` companion:

- verifies (11)--(15) exactly;
- tests (15) on `288` unrelated unit/high-mask triples;
- checks all thirteen affine guard digits in (20);
- exhausts the phase-grid endpoint chambers for `N=7,49,343`; and
- verifies the exact hostile control (24)--(25).

Stored output:

```text
05-knowledge/results/lrc14_last_lane_septimal_depth_cap_thm2415.out
```

An independent hostile audit reconstructed the divided variables,
transition mass, orbit disintegration, affine word transport, null-set
removal, and W8 contradiction. It also supplied (24) as the sharp
`M=2` control. No defect was found.
