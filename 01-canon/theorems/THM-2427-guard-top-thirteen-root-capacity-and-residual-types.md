---
id: THM-2427
title: "Guard-top thirteen-root capacity and residual type reduction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. In a
  primitive THM-2367 scalar cover with nu_7(c_3)>M and nu_7(H)=M, put
  k=#{j<=2:nu_7(c_j)<=M}, t=#{i:nu_7(q_i)=M},
  b=#{j<=2:nu_7(c_j)=M}, and W=2+t+b. If W>k, THM-2367 makes the top
  words cover every septimal top bin. On a generic thirteen-root fibre
  over which all three quotient blockers are safe, every physical
  blocker disappears; the guard and t top ordinary words have capacity
  at most 4+2t. Hence t=5. Exactly three (k,t,b,W) types remain at M=0
  and, after primitivity, four at M>0. With THM-2426 these seven
  regime-typed shapes are the complete current deep-c_3 valuation
  residual. A strict primitive 91-unit example with all three blockers
  absent exactly partitions the full compositional 13-by-7 stalk on a
  positive base chamber; it lifts self-similarly through every
  septimal depth. Thus even frozen joint root geometry cannot remove
  the t=5 types. It removes no thirteen-adic row and does not prove
  LRC(14).
source: codex-2026-07-26-guard-top-thirteen-root-capacity
depends_on:
  - THM-2367-septimal-root-averaging-graft-and-cover-alignment
  - THM-2426-compositional-thirteen-root-final-septimal-lane-exclusion
related:
  - THM-2382-saturated-septimal-seven-bin-root-fibre-closure
  - THM-2391-blocker-caged-septimal-single-layer-address-reduction
  - THM-2430-guard-top-common-ninety-one-root-tiling-spectrum
  - THM-2431-repeated-step-rounding-exclusion-of-guard-top-zero-blocker-types
  - THM-2432-guard-top-pair-cage-and-low-blocker-residual-exclusion
  - THM-2435-top-blocker-essential-parent-and-punctured-stalk-carrier
script: 04-computation/lrc14_guard_top_thirteen_root_capacity_thm2427.py
output: 05-knowledge/results/lrc14_guard_top_thirteen_root_capacity_thm2427.out
script_sha256: 19a92396131daa3366debac4bddc43a25858a4c07e0fd0536e433b4a5960b7a6
output_sha256: efdd5896e0aad3b87c06fd92757bc1bd6ee76b731b2da611e2c392f270715039
hash_basis: working-tree bytes (LF)
cite_by_filename: true
---

# THM-2427 -- a deep blocker and top guard force all five top ordinary words

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2426 closes the whole subtop-guard/deep-`c_3` side of the septimal
valuation split. The remaining deep-`c_3` side has a top guard. The
same two-scale operation which closed the final lane now gives a short
capacity reduction:

```text
seven-bin cover
  -> every top bin has a top word
  -> pass to thirteen inverse roots with all blockers absent
  -> guard capacity 4 plus ordinary capacity 2t covers 13
  -> t=5.                                                        (1)
```

The lower ordinary words are not simply discarded from the scalar
cover. THM-2367 first proves that a top word already covers every
point; only then can the lower words be ignored.

## 1. Typed hypotheses

Retain THM-2367's primitive almost-everywhere scalar cover

```text
C_H subset
  union_(i=1)^5 D_(q_i)
  union D_(c_1) union D_(c_2) union D_(c_3),                      (2)

D_v={x:||vx||<1/14},                 C_H=E_H^c,

c_j=13C_j,
```

where `H,q_1,...,q_5` are units modulo thirteen. Put

```text
M=max(nu_7(H),nu_7(q_1),...,nu_7(q_5))
```

and assume

```text
nu_7(c_3)>M,                         nu_7(H)=M.                    (3)
```

Use THM-2367's counts

```text
k=#{j in {1,2}:nu_7(c_j)<=M},

t=#{i:nu_7(q_i)=M},

b=#{j in {1,2}:nu_7(c_j)=M}.                                    (4)
```

Because the guard is top, the total top arc weight is

```text
W=2+t+b.                                                        (5)
```

## 2. Top-bin coverage when `W>k`

Let

```text
N=7^(M+1).
```

On every generic `N`-orbit where all blockers of depth greater than
`M` are safe, THM-2367 equation (34) gives

```text
7m_r>=W-k                         for every r in Z/7Z,             (6)
```

where `m_r` counts top arcs occupying residue bin `r`.

Assume

```text
W>k.
```

Then (6) gives `m_r>=1` for every `r`. A depth-`M` ordinary danger
word is constant on its selected whole top bin, and the depth-`M`
guard is constant on its two selected bins. Consequently the union of
the top words already covers every point of the orbit. Subtop ordinary
words may also be dangerous, but no point relies on them for coverage.

## 3. A blocker-free thirteen-root fibre

Put

```text
A=D_(C_1)^c intersection D_(C_2)^c intersection D_(C_3)^c.
```

The union bound gives

```text
mu(A)>=1-3/7=4/7>0.                                               (7)
```

Choose `y in A` after removing:

- all quotient-blocker endpoints;
- for every `h in F_13` and `s in Z/NZ`, the affine pullback under
  `y -> (y+h)/13+s/N` of the scalar-cover and inherited THM-2367
  orbit-identity null exceptional sets; and
- the corresponding finite pullbacks of all top-word endpoints.

Only finitely many null sets are removed, so such a generic `y`
exists. For the thirteen inverse roots

```text
x_h=(y+h)/13,                         h in F_13,                    (8)
```

one has

```text
1_(D_(c_j))(x_h)=1_(D_(C_j))(y)=0

for j=1,2,3.                                                       (9)
```

Every blocker of depth greater than `M` is constant and safe on the
full `N`-orbit of each `x_h`. Section 2 therefore applies at every
root. Its mandatory top-word cover cannot use a top blocker because
all physical blockers are absent by (9). Each root lies in the guard
or in one of the `t` top ordinary danger words.

## 4. Thirteen-root capacity

On the roots (8), each speed unit modulo thirteen gives a permuted
thirteen-grid. The guard arc has length `2/7`, so it meets at most

```text
ceil(13*2/7)=4
```

roots. Each ordinary danger arc has length `1/7`, so it meets at most

```text
ceil(13/7)=2
```

roots. Covering all thirteen roots therefore requires

```text
13<=4+2t.                                                        (10)
```

Since `0<=t<=5`, equation (10) forces

```text
t=5.                                                            (11)
```

The bound is only a necessary capacity statement. The exact companion
retains an abstract unit-affine guard/five-ordinary root cover as a
positive control, so (10) does not falsely exclude `t=5`.

The boundary persists for actual integer speeds at one common phase.
At

```text
y=1/2,                H=1,                (q_1,...,q_5)=(2,3,5,11,19),
```

all six speeds are units modulo `91`, and their root masks are

```text
guard:       {0,1,11,12},

ordinary:    {6}, {4,8}, {2,10}, {3,9}, {5,7}.                 (11a)
```

These six sets partition `F_13`. Thus even common-base physical
thirteen-root geometry cannot eliminate the `t=5` rows.

There is a stronger compositional hostile. Put

```text
y_0=11/581,

H=1,

(q_1,...,q_5)=(547,1821,3095,4369,5643),

(c_1,c_2,c_3)=(91,1183,15379).                                 (11b)
```

All six top speeds are units modulo `91`; in fact every `q_i` is
congruent to `1` modulo `91`. The blockers have strict thirteen-adic
depths `1,2,3` and septimal depth one. On the product stalk

```text
x_(h,s)=(y_0+h)/13+s/7,

(h,s) in F_13 x F_7,                                             (11c)
```

the CRT coordinate

```text
r=7h+13s mod 91
```

turns the six top masks into the exact consecutive partition

```text
guard:       {0,...,12,78,...,90},

ordinary:    {13,...,25}, {26,...,38}, {39,...,51},
             {52,...,64}, {65,...,77}.                           (11d)
```

The thirteen-root section `s=0` is already the exact partition

```text
guard:       {0,1,12},

ordinary:    {2,3}, {4,5}, {6,7}, {8,9}, {10,11}.                 (11e)
```

Writing `c_j=13C_j`, the three quotient-blocker phases have norms

```text
||C_1 y_0||=11/83,       ||C_2 y_0||=23/83,
||C_3 y_0||=33/83.
```

They are all strictly outside their danger arcs. Since every `C_j` is
divisible by seven, every physical blocker is therefore absent on the
whole stalk (11c). The top words alone give exact one-fold coverage.
Moreover every mask and blocker truth value persists for

```text
|y-y_0|<1/635614.                                                (11f)
```

This is a positive base chamber, not an endpoint coincidence.

The same word has a toothpick-like all-depth lift. For `m>=0`, replace

```text
(H,q_i) by 7^m(H,q_i),

c_j by 13^j 7^(m+1),

y_0 by y_0/7^m,

7 by 7^(m+1) in the second stalk coordinate.
```

Then

```text
7^m u ((y_0/7^m+h)/13+s/7^(m+1))
 =u(y_0+7^m h)/13+us/7.                                        (11g)
```

Reduction of `h` by multiplication with `7^m` modulo thirteen and of
`s` modulo seven recovers (11c), repeated `7^m` times. The quotient
blocker phases remain the same three safe values. For `m>0` the
displayed row is imprimitive, so this lift is a self-similar hostile
to a proof mechanism, not an additional primitive residual type.

Consequently even the complete frozen `13 x 7` root/bin stalk, with
lawful blockers retained, cannot exclude the primitive
`(k,t,b,W)=(0,5,0,7)` shape. The missing information is cross-base
phase and endpoint evolution; in the blocker-bearing types one must
also retain owner/blocker state and the labelled valuation-zero
sidecar. Neither example is a global scalar-cover packet.

## 5. Complete residual type list

THM-2367 also proves

```text
W<=k             or             W>=7.                            (12)
```

If `W<=k`, equations (4)--(5) and `k<=2` force the single type

```text
(k,t,b,W)=(2,0,0,2).                                             (13)
```

If `W>k`, equation (11) gives `t=5`. Letting `0<=b<=k<=2` gives exactly

```text
(0,5,0,7),

(1,5,0,7), (1,5,1,8),

(2,5,0,7), (2,5,1,8), (2,5,2,9),                                (14)
```

together with (13): seven types.

If `M=0`, nonnegativity of valuations makes every blocker counted by
`k` have depth exactly `M`, so `b=k`. Only

```text
(0,5,0,7),                (1,5,1,8),                (2,5,2,9)     (15)
```

remain.

If `M>0`, any type with `t=5` and `b=k` makes the guard, all five
ordinary labels, and all three blockers divisible by seven: there is
no blocker below `M`, while every other blocker has depth at least
`M`. Primitivity excludes precisely

```text
(0,5,0,7),                (1,5,1,8),                (2,5,2,9).
```

The positive-depth residual is therefore exactly

```text
(1,5,0,7),

(2,0,0,2),

(2,5,0,7),

(2,5,1,8).                                                       (16)
```

Each positive-depth row carries an additional valuation-zero sidecar.
In the three `t=5` types, at least one of `c_1,c_2` has valuation zero;
in `(2,0,0,2)`, at least one lower ordinary or low-blocker label has
valuation zero. The tuple alone does not record which label.

## 6. Consequence and scope

THM-2426 proves

```text
nu_7(c_3)>M          implies          nu_7(H)=M.
```

This theorem therefore gives the complete current valuation residual
on the deep-`c_3` side: the three types (15) at `M=0`, or the four
types (16) when `M>0`.

It does not prove that any listed type is realizable. It does not
control owner-mask cancellation, provide a canonical target, remove a
thirteen-adic profile row, or prove LRC(14). The scalar ledger remains
`165`.

## 7. Exact companion

Run:

```text
python3 04-computation/lrc14_guard_top_thirteen_root_capacity_thm2427.py
python3 -O 04-computation/lrc14_guard_top_thirteen_root_capacity_thm2427.py
```

The dependency-free companion:

- reconstructs every strict/boundary thirteen-root mask and its
  unit-affine closure;
- verifies the sharp support caps `4` and `2`;
- retains abstract and common-phase physical `t=5` covers as positive
  controls;
- exhausts the `(k,t,b,W)` arithmetic; and
- checks the three-type `M=0` and primitive four-type `M>0` lists.

Normal and optimized modes must reproduce

```text
05-knowledge/results/lrc14_guard_top_thirteen_root_capacity_thm2427.out
```

byte-for-byte.

## 8. Independent hostile audit

An independent proof audit checked the quantifiers on the
high-blocker-safe `N`-orbits, the `13N` affine pullbacks of the
almost-everywhere exceptional set, and the identity

```text
1_(D_(c_j))((y+h)/13)=1_(D_(C_j))(y).
```

It independently recovered the root caps `4` and `2`, the implication
`W>k => t=5`, the three `M=0` types, the four primitive `M>0` types,
and the required valuation-zero sidecar. Normal and optimized
companions byte-match the stored transcript. The physical `t=5` root
partition was independently found during the frontier audit; it is
only a capacity hostile, not a scalar-cover packet.
