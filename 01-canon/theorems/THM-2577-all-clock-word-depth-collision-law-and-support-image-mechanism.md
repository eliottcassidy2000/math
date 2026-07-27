---
id: THM-2577
title: "All-clock word-depth collision law and support-image mechanism"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.  On the
  canonical typed row, for every packet clock K whose prescribed word has
  positive return, the owner-normalized first-collision depth is exactly 3
  for sigma={a} and exactly 5 for sigma={b} or {a,b}.  The result is uniform
  in K, not a finite-clock extrapolation.  Before that depth the essential
  images of the owner atom and word atom are disjoint.  At the next push,
  the {a} images coincide in a 26-interval set of mass 6/7, while both
  {b}-word images are the whole circle.  Hence every nonzero returned
  subpacket collides at the stated depth.  This proves an all-clock r=5
  sidecar window on the two b-words, but it constructs no temporal
  intertwiner, physical target current, scalar-row exclusion, or LRC(14)
  conclusion.
source: root-frontier-final-2026-07-28-all-clock-collision
depends_on:
  - THM-2306-owner-normalized-disjoint-support-and-first-collision-shell
  - THM-2471-owner-first-collision-weighted-root-service-and-temporal-atom-boundary
related:
  - THM-2575-word-depth-collision-law-and-owner-clock-host-array
script: 04-computation/lrc14_all_clock_word_depth_support_thm2577.py
output: 05-knowledge/results/lrc14_all_clock_word_depth_support_thm2577.out
script_sha256: 6868063601639144bf029c84afdc5875c2df4784c1591b32df2cfc13d9113d1d
output_sha256: f1f5069feb90ed9a1bee356f13265112a1f58ec2949d9d62eb31ec0579b37b94
hash_basis: working-tree bytes (LF)
---

# THM-2577 -- the word-depth collision law is independent of the packet clock

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.**

THM-2575 computed the collision depth at packet clocks `K=2,3,4` and found
the same pattern at all three clocks:

```text
sigma={a}:                    r=3;

sigma={b} or {a,b}:           r=5.                         (1)
```

The repetition is not accidental.  The collision depth is already decided
by a finite atlas of **support images** which does not mention `K`.  Once that
atlas is exposed, (1) holds at every positive-return packet clock.

## 1. The all-clock statement

Retain the canonical typed row

```text
w=(1,14,27,40,53,66,13,13^3,2*13^5)                       (2)
```

and its unique licensed owner `j=1`.  Write `E=E_1` and, in the fixed
THM-2305 word order, write

```text
Q_a   =A_0 intersect D_(c_2) minus (D_(c_1) union D_(c_3)),

Q_b   =A_0 intersect D_(c_3) minus (D_(c_1) union D_(c_2)),

Q_ab  =A_0 intersect D_(c_2) intersect D_(c_3) minus D_(c_1). (3)
```

Here `A_0` contains the guard-safe and five unit-safe factors.  Let `K` be
**any** positive-return packet clock licensed by THM-2306/2471, and set

```text
f_(sigma,K)=1_(Q_sigma) P^K 1_E,

rho_(sigma,K)=integral_T f_(sigma,K)>0.                     (4)
```

With the owner speed `c_1=13`, put

```text
A=P_13 f_(sigma,K),              B=P_13 1_E,

I_s=integral_T (P^s A)(P^s B),

r=min{s>=1:I_s>0}.                                        (5)
```

Then, uniformly over every `K` satisfying (4),

```text
r({a},K)=3,

r({b},K)=r({a,b},K)=5.                                    (6)
```

Thus the `theta=t-2u` branch of THM-2471 is available at every
positive-return clock on both words containing `b`, not only at the three
clocks censused in THM-2575.

## 2. Support images commute with the collision question

For a measurable set `S` on the circle and an integer `m>=1`, write

```text
mathcal S_m(S)={m x mod 1:x in S}                            (7)
```

for its essential image.  If `g>=0`, the transfer formula

```text
(P_m g)(y)=1/m sum_(j=0)^(m-1)g((y+j)/m)                    (8)
```

gives the exact support identity

```text
ess supp(P_m g)=mathcal S_m({g>0}).                          (9)
```

Consequently:

- disjoint essential images in (9) force the product of the two transfers
  to vanish almost everywhere;
- if `g` is nonzero, `{g>0} subset Q`, and
  `mathcal S_m(Q) subset mathcal S_m(E)`, then `P_m g` and `P_m1_E` are
  simultaneously positive on the positive-measure image of `{g>0}`.

Transfer composition rewrites (5) as

```text
I_s=integral_T
      (P_(13^(s+1)) f_(sigma,K))
      (P_(13^(s+1)) 1_E).                                  (10)
```

The clock occurs only inside the nonzero function `f_(sigma,K)`.  Its
support is always contained in `Q_sigma`.  Thus a support-image atlas for
`E,Q_a,Q_b,Q_ab` decides (10) for **every** clock at once.

## 3. The exact image atlas

All sets in (3) are rational finite unions of intervals on the exact grid

```text
T_DEN=182 lcm(w)=297836897838480.                            (11)
```

Integer interval folding gives the following identities, modulo null
endpoints.

For the `a` word,

```text
mathcal S_(13^e)(E) intersect mathcal S_(13^e)(Q_a)=empty
                                      for e=1,2,3,            (12)

mathcal S_(13^4)(E)=mathcal S_(13^4)(Q_a).                   (13)
```

The common set in (13) consists of `26` intervals and has mass `6/7`.

For both words containing `b`,

```text
mathcal S_(13^e)(E) intersect mathcal S_(13^e)(Q_sigma)=empty
                                      for e=1,...,5,          (14)

mathcal S_(13^6)(E)=mathcal S_(13^6)(Q_sigma)=T,
                         sigma in {{b},{a,b}}.                (15)
```

The first part of each line reflects the opposite danger/safe status at the
deepest blocker named by the word: `nu_13(c_2)=3` for `{a}`, and
`nu_13(c_3)=5` for the two `b` words.  The terminal equalities (13), (15)
are the additional load-bearing fact; blocker valuation alone would give
only the lower bound.

## 4. Proof of the word-depth law

Fix a word and a clock satisfying (4).  Since

```text
{f_(sigma,K)>0} subset Q_sigma,                              (16)
```

equations (9), (10), and either (12) or (14) give

```text
I_s=0                for 0<=s<3,   sigma={a},

I_s=0                for 0<=s<5,   sigma={b},{a,b}.          (17)
```

At the next exponent, (13) or (15) contains the entire positive-measure
image of `{f_(sigma,K)>0}` inside the essential image of `E`.  Both transfer
functions in (10) are then positive on that image.  Therefore

```text
I_3>0                for sigma={a},

I_5>0                for sigma={b},{a,b}.                    (18)
```

Equations (17)--(18) prove (6).  Notice that no mixing estimate, clock
limit, or interpolation in `K` occurs.  Positive return is used exactly once:
it says that the returned subpacket in (16) is nonzero.

## 5. Boundaries, equality cases, and scope

The statement is sharp at three boundaries.

1. If the return in (4) is zero, `f_(sigma,K)=0` and no collision depth is
   defined.  The theorem does not promote a zero-return clock.
2. The empty-word residual is not a canonical THM-2305 word on this typed
   non-cover.  Its diagnostic collision depth is `1`, so including it would
   make (6) false.
3. The terminal image inclusion is row-specific.  Opposite blocker status
   by itself only proves the zeros in (17); THM-2306 has explicit delayed
   collision hostiles.  No uniform statement over the other `164` live rows
   is inferred.

The result strengthens the *availability* of the `r=5` sidecar but does not
construct its missing map.  In particular, THM-2575's `13 x 7` owner-clock
host was built on the separate `r=3`, `{a}` packet.  Multiplying it by the
present `r=5` statement would mix word strata and would repeat the
common-base error of MISTAKE-281.  The next lawful experiment is to build the
owner-clock array directly inside the `{b}`, `r=5` fibre and then test the
`theta=t-2u` toothpick contraction there.

No temporal atom intertwiner, physical endpoint current, semantic arrival,
scalar-row exclusion, or LRC(14) conclusion follows.  The live ledger remains
`165`.

## 6. Exact companion

Run

```bash
python3 04-computation/lrc14_all_clock_word_depth_support_thm2577.py
python3 -O 04-computation/lrc14_all_clock_word_depth_support_thm2577.py
```

Both executions byte-match

```text
05-knowledge/results/lrc14_all_clock_word_depth_support_thm2577.out.
```

The dependency-free referee reconstructs the four atoms in (3), checks every
zero image intersection in (12), (14), proves the terminal equality/inclusion
in (13), (15), and tests first, middle, and last positive subintervals of
each word as hostile returned subpackets.  The support lemma (9) and the
all-clock quantifier are symbolic; they are not finite extrapolations.

**QED.**
