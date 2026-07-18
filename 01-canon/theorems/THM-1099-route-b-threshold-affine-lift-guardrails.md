---
id: THM-1099
title: Route-B threshold, affine-carry, and integer-lift guardrails
status: PROVED arithmetic guardrails plus an exact cross-modulus witness; the AP-core inverse theorem and LRC(14) remain open
date: 2026-07-18
source: codex-2026-07-18-S67 Route-B audit after death-star-S57 and boxeph-S101
depends_on:
  - THM-366
  - THM-668-pair-sum-ruler-witness-structure
  - THM-769
  - THM-1008
  - THM-1017
  - THM-1091
related:
  - THM-1002
  - THM-1105
  - MISTAKE-164
results:
  - 05-knowledge/results/lrc14_routeb_affine_lift_guardrails_codex_S67.out
source_sha256: f1c5b1722cf1f837b2ecb65dcf654bdc0739371407f130d699442084e63858b5
output_sha256: aec003e07edeb898e3d1d41b65257301b222e13118aefa92ece8c28a63728408
---

# THM-1099 — Route-B threshold, affine-carry, and integer-lift guardrails

The correction in death-star-S57 is essential: at the `1/13` threshold,
covering through `13` and covering through `14` are different predicates.
This theorem records the further guardrails needed before the residue-gap
picture of boxeph-S101 can be used as a proof step.

The main point is a separation of three objects which had been conflated:

1. the **safe critical band** `1/14 < val/q < 1/13`;
2. a canonical-representative **Euclidean offset chart** modulo `val`; and
3. the **integer lifts** of a residue modulo `q`, on which covering at other
   moduli actually acts.

The first is disjoint from the hypothetical LRC-counterexample band, the
second is not a cyclic quotient when `q = 13 val + 1`, and the third is
invisible at one maximizing residue picture.  These are proved facts, not a
new inverse theorem.  LRC(14) is not claimed.

## 1. The exact threshold trichotomy

Let `s,q` be positive integers and suppose

```text
s/q < 1/13.
```

Write

```text
q = 13s + e,       e = q-13s > 0.                         (1)
```

Clearing positive denominators gives the exact trichotomy

```text
1/14 < s/q < 1/13    iff    0 < e < s,
s/q = 1/14           iff    e = s,
s/q < 1/14           iff    e > s.                         (2)
```

In particular, `s/q<1/13` proves only `13s<q`.  The further inequality
`q<14s` is *equivalent* to `s/q>1/14`; it is not a consequence of the
upper `1/13` bound.  Thus an argument which assumes

```text
13s < q < 14s
```

is analyzing the already-safe open band in (2), not a hypothetical LRC(14)
counterexample.  A counterexample would lie on the opposite side `q>14s`.

Nor does the safe band imply `e=1`.  The specialization
`q=13s+1` is the first possible denominator in that band, but the step

```text
q=13s+e  --->  e=1                                      (3)
```

is an additional spectrum theorem.  It is observed on the deep-well ladder;
it is not supplied by the inequalities in (2).

### The signed fourteen-gap identity

Let thirteen distinct canonical residues be ordered as

```text
s <= r_1 < ... < r_13 <= q-s.
```

Adjoin the observer and one wrap endpoint,

```text
r_0=0,       r_14=q,       g_i=r_(i+1)-r_i  (0<=i<=13).
```

Then telescoping gives the exact defect identity

```text
sum_(i=0)^13 (g_i-s) = q-14s = e-s.                       (4)
```

The boundary gaps `g_0,g_13` are at least `s`.  Hence in the safe critical
band `e<s`, (4) forces an **interior** gap below `s`.  This is the rigorous
content of the one-small-gap pigeonhole.  In the counterexample band `e>s`,
the total defect has the opposite sign and the pigeonhole supplies no such
conclusion.  Formula (4), rather than the number of small gaps alone, is the
threshold-correct gap object.

## 2. The `e=1` packing still has a carry case

Even after imposing `q=13s+1`, twelve points in
`[s,12s+1]` separated by at least `s` are not forced to equal
`s*{1,...,12}`.  For every `s>=2`, the set

```text
{s,2s,...,6s, 7s+1,8s+1,...,12s+1}                       (5)
```

lies in the required interval, contains `s`, and has adjacent gaps `s`
except for one gap `s+1`.  This is exactly the carry subcase already left
open in the proof discussion of the S90 packing lemma.

Therefore the implication

```text
twelve s-separated residues  --->  exactly s*{1,...,12}   (6)
```

needs an extra no-carry hypothesis.  More generally, when `q=13s+e`, the
total excess budget for twelve separated points is **at most** `e`, with
equality only when the packet saturates both endpoints.

There are two further logically separate steps after (6): the exceptional
residue must be shown to belong to `v_max`, and a modular arithmetic
progression must be lifted to an arithmetic progression in the actual
integer speeds.  Neither follows from the sorted residue gaps alone.

## 3. The offset coordinate is affine, not a ramified quotient

Suppose `q=13s+e`.  For a chosen integer representative `x`, its Euclidean
offset `x mod s` changes by `e mod s` when the representative is replaced by
`x+q`:

```text
(x+q) mod s = (x+e) mod s.                                (7)
```

Thus reduction modulo `s` is not a well-defined map on `Z/qZ` unless
`s|q`, equivalently `s|e`.  In the empirically important case `e=1`, it has
monodromy one.

The same fact appears group-theoretically.  If `gcd(s,q)=1`, multiplication
by `s` is a permutation of `Z/qZ`, so

```text
s*(Z/qZ) = Z/qZ,                                          (8)
```

and every additive translate of this set is again the whole group.  Since
`gcd(s,13s+1)=1`, phrases such as "two cosets of `s Z/qZ`" are vacuous at
`q=13s+1`: there is only the whole group, not a proper subgroup quotient.

What does exist is an affine carry chart.  For `0<=x,y<q`, put

```text
epsilon = 0 if x+y<q, and 1 otherwise,
z = x+y-epsilon*q = (x+y) mod q.
```

Then

```text
z = x+y-epsilon*(13s+e),
z mod s = (x+y-epsilon*e) mod s.                           (9)
```

Every wrap subtracts the cocycle `e`; at `e=1` it subtracts one.  Ordinary
coset/Freiman language which drops this carry loses precisely the
archimedean coordinate it is supposed to control.

This separates the present setting from the honest ramified quotients in
THM-1091 and `LRCRamifiedCosetCover`: those use maps `Z/cZ -> Z/DZ` with
`D|c`, so complete fibres and annihilator Fourier support really exist.
Route B instead needs an **affine-carry extension** of that formalism.

## 4. Residues do not determine integer speed lifts

Let `a` be a unit modulo `q`, and let `r` be the residue of a speed at
`t=a/q`.  If `u` is the least nonnegative solution of

```text
a*u = r (mod q),
```

then every positive speed with that residue is an integer lift

```text
v = u+kq.                                                  (10)
```

Consequently an AP in the residues yields only congruences among the speeds.
For example, at `q=183,a=14`, speeds `2` and `185=2+183` both have residue
`28`, although `185` is not twice the speed `1`.  Turning a residue AP into
an actual AP requires a no-wrap bound or a theorem controlling the lift
coordinates `k`; the safe residue band bounds `r`, not `v`.

Covering lives exactly on these lift coordinates.  For positive `d`, the
linear congruence theorem gives

```text
there exists k with d | (u+qk)    iff    gcd(q,d) | u.      (11)
```

When solutions exist, they form one residue class modulo
`d/gcd(q,d)`.  Proof: necessity follows because `gcd(q,d)` divides both
`qk` and `u+qk`; for sufficiency divide by the gcd and invert
`q/gcd(q,d)` modulo `d/gcd(q,d)`.  An integer solution can be shifted by
that positive period to make `k>=0`.

Equation (11) is the exact cross-modulus bridge missing from a one-point
residue picture.  The cover obligations `d=2,...,14` are owner-local
congruences on one **shared lift word**.  Keeping only the assertion that each
owner has some local carrier discards their compatibility, just as the
pre-nerve guardrail warns.

## 5. An exact lift which keeps the local maximum and changes the global value

The death-star-S57 false alarm is

```text
V_0={1,2,3,5,7,8,9,10,11,12,17,19,104}.
```

Its exact value and a maximizer are

```text
M(V_0)=8/105,       t_*=8/105.                             (12)
```

It is primitive, covers every modulus `2,...,13`, misses `14`, and has a
non-AP twelve-speed core.  At `t_*` its sorted residues and internal gaps are

```text
residues = 8,16,24,31,40,47,56,64,72,80,88,96,97,
gaps     = 8, 8, 7, 9, 7, 9, 8, 8, 8, 8, 8, 1.             (13)
```

Now replace the **interior** speed `7` by its same-residue lift

```text
112 = 7+105,
V_1={1,2,3,5,8,9,10,11,12,17,19,104,112}.                 (14)
```

Because `8*112 = 8*7 (mod 105)`, the entire residue multiset (13) at `t_*`
is unchanged.  The active pair `1,104` is unchanged as well:

```text
8*1 = 8,       8*104 = 97 = -8 (mod 105).
```

All other residues have clearance at least `9/105`; the active distance
branches have slopes `+1` and `-104`.  Thus the same opposite-slope active
pair makes `t_*` a strict local maximum of height `8/105` for both rows.
The local maximum, its active pair, its offset word, and all
of its gap data survive the lift.

Globally the rows are completely different.  The lift `112` is divisible by
`14`, so `V_1` covers every modulus `2,...,14`.  At

```text
t=3/20
```

the thirteen numerator distances are

```text
3,6,9,5,4,7,10,7,4,9,3,8,4,
```

so

```text
M(V_1) >= 3/20 > 1/13.                                    (15)
```

The complete pair-sum-ruler evaluator in the deterministic audit gives the
exact equality `M(V_1)=3/20`.  Completeness here is the proved THM-668/THM-1002
candidate theorem: at a global maximum below `1/2`, oppositely sloped active
runners `x,y` give a reduced denominator dividing `x+y`.  The evaluator includes
all pairs (also `x=y`), every divisor of every pair sum, and every reduced
numerator at that denominator; it is not a bounded-denominator heuristic.

This is a literal realization of S101's valid qualitative diagnosis:
interior residue gaps are invisible to the local germ.  It also identifies
the missing positive object.  Adding covering at modulus `14` acts on an
integer lift, and that lift creates a new witness at a different denominator
(`20`), far from the old maximizing denominator (`105`).  The global
cross-modulus comparison, not one more statistic of (13), is load-bearing.

## 6. The exact logical frontier

Write `Cover14(V)` for covering every modulus `2,...,14` and
`rho=v_max/v_2nd`.  Given the proved non-covering sieve and THM-1008, the
LRC-level residual is

```text
Cover14(V) and rho<13    implies    M(V)>=1/14.             (16)
```

The commonly stated compact floor

```text
Cover14(V) and rho<13    implies    M(V)>=1/13              (17)
```

is stronger.  The AP-core inverse theorem

```text
Cover14(V) and M(V)<1/13
  implies V minus {v_max} is a dilated twelve-term AP       (18)
```

is stronger again: with the proved AP-core/lcm half it implies (17), which
implies (16).  No reverse implication in this chain has been proved.  In
particular, (18) is a sufficient route to LRC(14), not a logically established
equivalent reformulation of it.

Moreover, the gap proof in the safe band proves neither (16) nor (18) on the
counterexample side of (2).  A correct completion must either prove (16)
directly or prove a genuinely global, lift-sensitive version of (18) which
also handles `q>=14s`.

## 7. Carrier and tournament audit

The proof-bearing finite object suggested by (10)--(11) has

```text
vertices = cover obligations d=2,...,14 and integer lift sheets,
hyperedges = one speed lift jointly carrying a set of obligations,
global word = the thirteen compatible lift coordinates.
```

It preserves cross-modulus covering and the actual integer speeds.  It loses
the continuum maximum unless rational witness obligations are added as a
second owner family.  Conversely, a winding tournament on the thirteen
runners preserves cyclic order at one time but destroys the lift coordinates
and the value of `M`; the pair `V_0,V_1` proves that loss concretely.

For methodology, the deterministic audit also forms the owner-cost
tournament: orient `d -> e` when the least available lift carrying `d` is
smaller than the least available lift carrying `e`, breaking equal costs by
the modulus label.  Its Hamiltonian path is the sorted owner-cost list; by
construction this tournament is always transitive (zero directed cycles,
singleton SCCs, one Hamiltonian path).  This fingerprint is useful telemetry but deliberately
not proof-bearing: pairwise cost order forgets whether several owners use the
same lift.  The owner--carrier hypergraph, or the pre-nerve shared-assignment
model, is the faithful quotient.

This challenges the default choice of runners as tournament vertices.  Here
the natural vertices are proof obligations and lift sheets; tournament data
is a diagnostic shadow of a ramified/affine cover problem.

## 8. Deterministic audit

`04-computation/lrc14_routeb_affine_lift_guardrails_codex_S67.py` uses exact
integer and rational arithmetic to check:

- the threshold trichotomy and signed gap defect;
- the `s+1` packing carry family;
- cyclic multiplication surjectivity and the affine carry law;
- the complete pair-sum-ruler values in (12) and (15);
- equality of the local residue data, strict active-pair slack, and covering;
- the CRT lift `7 -> 112`; and
- the owner-cost tournament fingerprints.

The audit is evidence for the concrete finite claims; the general identities
above have the displayed elementary proofs.

Frozen artifact hashes:

```text
source  41cf86eb9baeabe78a8cd77e72fab01bc608f004b3d15693a417f7ad3df8bcce
output  aec003e07edeb898e3d1d41b65257301b222e13118aefa92ece8c28a63728408
```
