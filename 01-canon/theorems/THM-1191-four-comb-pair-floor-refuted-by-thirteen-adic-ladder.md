---
id: THM-1191
title: The proposed four-comb pair floor is refuted by the thirteen-adic ladder
status: PROVED exact infinite-family refutation and exact four-subset-averaging guardrail.  For every u>=1 the four distinct combs at speeds u,13u,169u,2197u have total global pair overlap 210/2197<13/135.  Consequently the proposed R4 floor and its black-box derivation of R6>=13/54 are invalid.  No corrected universal R4 floor and no refutation of a direct R6>=13/54 theorem are claimed
source: codex-2026-07-18 pair-mass four-floor audit
depends_on: [THM-1166]
related: [THM-1176, THM-1179]
script: 04-computation/lrc14_four_comb_pair_floor_refutation_codex_20260718.py
output: 05-knowledge/results/lrc14_four_comb_pair_floor_refutation_codex_20260718.out
---

# THM-1191 -- four-comb pair-floor refutation

Put

```text
D_s={t in R/Z: ||st||<1/14},
rho(a,b)=|D_a intersect D_b|,
R(S)=sum_{{a,b} subset S} rho(a,b).                    (1)
```

The proposed universal inequality

```text
R(S)>=13/135 whenever |S|=4                              (2)
```

is false.  The obstruction is not an isolated bounded tuple.  It is an
infinite projective family carried by four consecutive levels of the
`13`-adic scale ladder.

## 1. Exact overlap word on the ladder

For coprime positive integers `A<B`, the radius-`1/14` pair-overlap formula
of THM-1166 is

```text
rho(A,B)
 =[4AB+f(A+B)-f(B-A)]/(196AB),                          (3)
f(x)=r(14-r),  r=x mod 14 in {0,...,13}.                (4)
```

Pair overlap is unchanged when both speeds are multiplied by the same
positive integer.  Hence, for `i<j`, the pair
`(u13^i,u13^j)` reduces to `(1,13^d)`, where `d=j-i`.
Since `13^d` is congruent to `-1` modulo `14` for odd `d` and to `1` for
even `d`, (3) gives the exact alternating word

```text
rho(1,13^d)=1/49-(6/49)13^(-d),  d odd,                (5)
rho(1,13^d)=1/49+(6/49)13^(-d),  d even.               (6)
```

Indeed, in the odd case

```text
f(13^d+1)=f(0)=0,   f(13^d-1)=f(12)=24,               (7)
```

and in the even case the two values are reversed.  This proves (5)--(6)
directly from (3), for every exponent gap and with no bounded search.

Now let

```text
S_u={u,13u,169u,2197u}.                                (8)
```

There are three exponent gaps of length `1`, two of length `2`, and one of
length `3`.  Therefore

```text
R(S_u)
 =1/49 [3(1-6/13)+2(1+6/13^2)+(1-6/13^3)]
 =210/2197.                                             (9)
```

Equivalently, the six exact pair values are

```text
gap 1, three times:  1/91,
gap 2, two times:   25/1183,
gap 3, one time:   313/15379.                          (10)
```

Finally,

```text
13/135-210/2197=211/296595>0.                          (11)
```

Thus every `S_u`, `u>=1`, refutes (2).  In particular the earlier finite
candidate `(1,12,27,40)`, whose value is `13/135`, is not extremal.

## 2. What failed: the scale stalk, not a finite sporadic row

The four speeds in (8) are a toothpick-self-similar carrier.  Dividing a
pair by its gcd moves its lower endpoint to exponent level zero; multiplying
the whole packet by `13` shifts every level without changing any overlap.
The correction to the bulk value `1/49` alternates with exponent-gap parity
and decays exactly as `13^{-d}`.  The object retained by the pair functional
is therefore

```text
ordered exponent stalk + parity colour + exponential gap weight,          (12)
```

not merely a four-element set of runner labels.  Any all-range reduction
which bounds the largest speed before quotienting this stalk can miss the
family at arbitrarily large height.  This is the precise all-scale reason a
finite scan suggesting `13/135` did not constitute a floor proof.

The exact referee scans the much narrower family
`(1,q,q^2,q^3)`, `2<=q<=200`, only as telemetry: `q=13` is the unique member
of that finite bank below `13/135`.  That scan is not used in (5)--(11) and
is not promoted to a corrected universal lower bound.

## 3. Six-wall four-subset averaging guardrail

Let `T` be any six-speed set and, for every four-subset `A` of `T`, form
`R(A)`.  Each of the fifteen pairs of `T` lies in exactly

```text
C(4,2)=6                                                (13)
```

of the fifteen four-subsets.  Hence the exact incidence identity is

```text
sum_(A subset T, |A|=4) R(A)=6R(T).                    (14)
```

Consequently a valid **uniform** four-comb floor `R(A)>=L` would imply

```text
R(T)>=(15/6)L=(5/2)L.                                  (15)
```

But (9) forces every such universal `L` to satisfy

```text
L<=210/2197.                                           (16)
```

Thus the black-box route (14)--(15), using only one uniform four-comb floor,
can certify at most

```text
R(T)>=(5/2)L <= 525/2197,                              (17)
13/54-525/2197=211/118638>0.                           (18)
```

The direction of (17) is a **capacity statement about the method**: the
derived lower-bound constant `(5/2)L` cannot exceed `525/2197`.  It is not an
upper bound on the actual six-comb pair mass.  In particular THM-1191 does
not refute a direct, compatibility-sensitive theorem `R(T)>=13/54`; it
refutes only the claimed derivation of that theorem from `R_4>=13/135`.

## 4. Consequence for the post-chi slow-gap inequality

In the six-on-one-slow-gap argument, a global six-pair floor `R(T)>=B`
enters the interval lower bound through the mean term

```text
(6/(7a))B.                                             (19)
```

The rejected four-floor would have supplied `B=13/54`, hence normalized
mean constant

```text
(6/7)(13/54)=13/63.                                    (20)
```

A black-box uniform-four-floor argument is capped instead by

```text
(6/7)(525/2197)=450/2197,
13/63-450/2197=211/138411.                             (21)
```

After the standard `rho(1-rho)<=13/196` conversion from pair endpoint debt
to reciprocal-gcd debt, these constant terms are multiplied by `196/13`.
Thus the proposed post-chi constant `28/9` is beyond the capacity of this
route; even if the best universal four-floor equalled the ladder value, the
route would supply only

```text
88200/28561,
28/9-88200/28561=5908/257049.                          (22)
```

The `H`-drift coefficient `72/13` is unaffected.  In the notation
`delta=aH-1`, the proposed black-box RHS

```text
28/9-(72/13)delta                                      (23)
```

must therefore be discarded.  The current rigorous `R_6>=5/24` input and
its RHS `35/13-(72/13)delta` remain valid.  Any improvement past it must use
genuine six-way compatibility, a nonuniform quartet ledger, chi classes,
Fano incidence, or another carrier which sees how the fifteen four-subsets
coexist.

## 5. Tournament and alternate-carrier audit

Tournament vertices were deliberately challenged: take the four exponent
levels `0,1,2,3`, not the four runners or their danger arcs.  On a pair
`i<j`, observe the signed folded correction

```text
rho(13^i,13^j)-1/49.                                   (24)
```

Gauge it into an orientation by pointing upward in exponent when (24) is
negative and reversing when it is positive; natural exponent order is the
tie Hamiltonian path, though this packet has no ties.  The resulting edges
are

```text
0->1, 2->0, 0->3, 1->2, 3->1, 2->3.                  (25)
```

Its fingerprint is

```text
score multiset {1,1,2,2}; directed triangles 012 and 123;
one SCC of size 4; five Hamiltonian paths; two flips from natural order.   (26)
```

This tournament preserves the alternating parity sign and exposes the
strongly connected obstruction, but destroys the magnitudes `13^{-d}` and
therefore does not determine (9).  The faithful carrier is (12).  Alternate
vertices considered in this audit were runners, slow gaps, fixed sections,
wall events, residues, valuation levels, and proof obligations.  Exponent
levels are the minimal quotient preserving the sign mechanism; attaching
gap weights restores the exact pair functional.

## 6. Exact replay

`lrc14_four_comb_pair_floor_refutation_codex_20260718.py` checks with exact
`Fraction` arithmetic:

1. the folded formula against an independent capped-trapezoid evaluator on
   `4064` pairs;
2. (5)--(6) through exponent gap `8` and (9) at five common scales;
3. every pair in (10) and all constants (11), (18), (21), and (22);
4. the six-vertex incidence multiplicity in (14); and
5. every tournament fingerprint in (26).

The script succeeds normally and under `python3 -O`; the checked-in output
is byte-identical to both runs.  SHA-256:

```text
script eafec292eb65d5e15368f53d06e8ae317029f2679ae3de151d335b70f8738728
output d473c67a82d2ec35f93f637c0c61c40649d5e1e3c7f3e576c418875bb159748d
```
