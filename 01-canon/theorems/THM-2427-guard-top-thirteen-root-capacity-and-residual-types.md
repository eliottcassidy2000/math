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
  residual. A common-phase 91-unit example exactly partitions the
  thirteen roots with one guard and five ordinary words, so root
  geometry alone cannot remove the t=5 types. It removes no
  thirteen-adic row and does not prove LRC(14).
source: codex-2026-07-26-guard-top-thirteen-root-capacity
depends_on:
  - THM-2367-septimal-root-averaging-graft-and-cover-alignment
  - THM-2426-compositional-thirteen-root-final-septimal-lane-exclusion
related:
  - THM-2382-saturated-septimal-seven-bin-root-fibre-closure
  - THM-2391-blocker-caged-septimal-single-layer-address-reduction
script: 04-computation/lrc14_guard_top_thirteen_root_capacity_thm2427.py
output: 05-knowledge/results/lrc14_guard_top_thirteen_root_capacity_thm2427.out
script_sha256: 36b407c2bb9935b141fff89f0e3e869b4b0db51da233912d1e6181bbe0a5e244
output_sha256: 2aad369e45e2cb84c9722cfed80d7ae7333a9a380d2b9cc4bb0dad740a5614d5
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
thirteen-root geometry cannot eliminate the `t=5` rows. Any such
elimination must use information discarded here: the simultaneous
septimal top-bin word, owner/blocker factors, or the labelled
valuation-zero sidecar. The example is a root-fibre hostile, not a
realizable scalar-cover packet.

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
