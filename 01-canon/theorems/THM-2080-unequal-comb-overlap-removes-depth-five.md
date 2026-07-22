---
id: THM-2080
title: "Unequal comb overlap removes the depth-five LRC14 tower"
status: >
  PROVED. If h is odd, the open 1/7 guard comb of h and the open 1/14
  danger comb of q overlap in measure at least 1/42; equality holds exactly
  when q=6h. Hunter's star bound then proves that the closed 1/14-safe set
  of any set of at most six distinct speeds cannot lie in one odd guard
  comb. Consequently a nontrivial THM-2073 terminal has at least seven
  speeds and the dyadic deletion tower has depth at most four. This removes
  the entire former r=5 lane, but does not close hereditary terminals of
  sizes 7 through 10 or LRC(14).
source: codex-2026-07-22-LRC14-unequal-comb
depends_on:
  - THM-638
  - THM-2073
related:
  - THM-2076
  - THM-2077
  - THM-2078
  - MISTAKE-239
script: 04-computation/lrc14_unequal_comb_overlap_referee_codex_20260722.py
output: 05-knowledge/results/lrc14_unequal_comb_overlap_referee_codex_20260722.out
script_sha256: 0ae2220c80f1eba0cf93819d4b442d88dc6ceeee06ae851a24eeb1f3ebd12696
output_sha256: 3cccf25892254dcc129448e0ab53060c26f7e10b04158b366de3b3b3095955d5
hash_basis: repository blobs with LF line endings
---

# THM-2080 -- unequal comb overlap removes depth five

Put

```text
E_h={t in R/Z:||ht||<1/7},
D_q={t in R/Z:||qt||<1/14},
G_Q=intersection_(q in Q) (R/Z minus D_q).              (1)
```

Thus `E_h` is the strict eligibility set of an odd dyadic-tower guard and
`G_Q` is the closed `1/14`-safe set of a quotient core.

The theorem has two parts.

1. For positive integers `h,q` with `h` odd,

   ```text
   measure(E_h intersect D_q)>=1/42,                    (2)
   ```

   with equality if and only if `q=6h`.
2. If `Q` is a set of at most six distinct positive integers and `h` is odd,
   then

   ```text
   G_Q is not contained in E_h.                         (3)
   ```

In fact the difference in (3) has positive measure. Applied to the last
guard in THM-2073, this proves that a nontrivial terminal quotient has at
least seven speeds. Since `|Q_r|=11-r`,

```text
r<=4.                                                    (4)
```

The proof is an unequal-radius companion to THM-638. Its faithful pairwise
observable is literal comb intersection mass; no tournament orientation is
introduced.

## 1. Exact unequal one-sided pair law

For `0<alpha,beta<1/2` and a positive integer `v`, write

```text
A_v(alpha)={t:{vt} in (0,alpha)}.                       (5)
```

Let `a,b` be coprime positive integers, and put

```text
x={a beta},                  y={b alpha}.               (6)
```

The offset enumeration in THM-638 works without requiring the two interval
lengths to agree and gives

```text
measure(A_a(alpha) intersect A_b(beta))
 =alpha beta+[min(x,y)-xy]/(ab),                        (7)

measure(A_a(alpha) intersect A_(-b)(beta))
 =alpha beta+[(x+y-1)_+-xy]/(ab).                       (8)
```

For completeness, decompose `A_a(alpha)` into its `a` half-open teeth and
`A_b(beta)` into its `b` teeth. Differences of their left endpoints run
bijectively through the `1/(ab)` grid by Bezout. At a point in the first
tooth, the number of translated second teeth covering it is

```text
floor(a beta)+1_{ {ab t}<x }.
```

The first tooth has scaled length `b alpha=floor(b alpha)+y`. Integrating
the indicator gives `floor(b alpha)x+min(x,y)`. Expanding the two floor
terms yields (7). Reversing the second tooth replaces `min(x,y)` by
`(x+y-1)_+`, proving (8). Endpoints are null and therefore do not affect
these identities.

For general `h,q`, put

```text
g=gcd(h,q),             a=h/g,             b=q/g.       (9)
```

Multiplication by `g` preserves Haar measure, so the intersection depends
only on `(a,b)`. Since a centered comb is the disjoint-a.e. union of its two
one-sided combs, (7)--(8) with

```text
alpha=1/7,              beta=1/14,
x=(a mod 14)/14,        y=(b mod 7)/7                   (10)
```

give the exact formula

```text
measure(E_h intersect D_q)
 =2/49+(2/(ab)) F(x,y),                                 (11)

F(x,y)=min(x,y)+(x+y-1)_+-2xy.                          (12)
```

This formula is useful beyond the tower: it records exactly when a fat
guard comb and an ordinary danger comb correlate or anticorrelate.

## 2. The sharp `1/42` floor

For every `x,y in [0,1]`,

```text
F(x,y)>=-1/8.                                           (13)
```

Indeed, if `x+y<=1` and `x<=y`, then

```text
F=x(1-2y),        0<=x<=1-y,
```

whose minimum is at least the minimum of
`(1-y)(1-2y)`, namely `-1/8`; the case `y<=x` is symmetric. If
`x+y>=1` and `x<=y`, then

```text
F=(2x-1)(1-y),
```

and the same quadratic bound applies after reflection. Again the other case
is symmetric.

If `ab>=15`, equations (11)--(13) give the strict estimate

```text
measure(E_h intersect D_q)
 >=2/49-1/(4ab)
 >=2/49-1/60
 =71/2940
 >1/42.                                                 (14)
```

It remains only to inspect coprime pairs with `ab<=14`. Pairs for which
`F>=0` already give at least `2/49`. The complete negative-correction table
is

```text
(a,b)       intersection             (a,b)       intersection
(1,4)          1/28                  (1,5)          1/35
(1,6)          1/42                  (1,11)         3/77
(1,12)         1/28                  (1,13)         3/91
(2,5)          1/35                  (3,4)          1/28
(8,1)          1/28                  (9,1)          2/63
(10,1)         1/35                  (11,1)         2/77
(12,1)         1/42                  (13,1)         3/91.       (15)
```

Every entry is at least `1/42`, and equality occurs only at `(1,6)` and
`(12,1)`. Because `h` is odd, both `g` and `a=h/g` are odd, excluding
`a=12`. The remaining equality `(a,b)=(1,6)` says exactly

```text
q=6h.                                                   (16)
```

This proves (2), including its equality statement. Notice that oddness is
used only to remove the reflected exceptional ratio in (15); it is exactly
the parity available for the tower guards.

## 3. Hunter's star rules out a six-speed terminal

Let `Q={q_1,...,q_s}` with `s<=6`, and use the star on the `s+1` events

```text
E_h,D_(q_1),...,D_(q_s),
```

whose edges join `E_h` to every `D_(q_i)`. Hunter's spanning-tree inequality,
proved pointwise in THM-638, and (2) give

```text
measure(E_h union union_i D_(q_i))
 <=2/7+s/7-s/42
 =(12+5s)/42.                                          (17)
```

For `s<=5`, the last quantity is at most `37/42<1`. For `s=6`, it is `1`,
but equality in every one of the six pair bounds would require all six
speeds to equal `6h`. Since `Q` is a set, at most one can do so. Hence the
sum of the six intersections is strictly greater than `6/42`, and (17) is
again strict. The complement has positive measure and is exactly

```text
{t:||ht||>=1/7 and ||qt||>=1/14 for every q in Q}
 =G_Q minus E_h.                                       (18)
```

This proves (3).

Now take a nontrivial THM-2073 tower. Its terminal safe set satisfies

```text
G_(Q_r) subset E_(h_(r-1))                             (19)
```

for the preceding odd guard, and `|Q_r|=11-r`. Equation (3) makes
`|Q_r|<=6` impossible. Therefore `|Q_r|>=7` and (4) follows. In particular,
THM-2077's former `r=5` six-comb saturation ledger is now a vacuous hostile
control rather than a live residual.

## 4. Scope and assumption challenge

This theorem removes all nontrivial tower depths `5,6,7,8`; THM-2076 had
already removed `6,7,8` by a scalar Haar-capacity tax. It does not treat the
depth-zero hereditarily primitive core, nor the remaining terminal sizes
`7,8,9,10` at depths `4,3,2,1`. Those lanes still require the original two
tail owners and their address/H-drift sidecar.

The challenged assumption is that the guard should be studied only through
its total capacity `2/7`. Pairing it separately with each ordinary danger
comb retains just enough arithmetic: a universal overlap floor plus one
rigid equality ratio. The useful graph is Hunter's undirected star, whose
edges carry exact intersection masses. Orienting those edges would add no
predicate and would discard the equality ledger; Tournament Analysis is
therefore not naturally available here.

## Correction lineage

The first pushed consumer reversed the containment implication and retained a
spurious `ab<=36` rank-six lane. MISTAKE-239 records the failure. From
`G_Q subset E_h`, it is `E_h^c` that must be covered by danger combs; the
overlap floor then rules out rank six completely, as Section 3 proves.

## Exact referee

The companion reconstructs (7)--(12) from exact interval atoms, checks the
formula on thousands of pairs, exhausts the small-product table and equality
classification, audits the Hunter arithmetic, and tests the containment
consequence on a finite hostile bank. It uses explicit runtime checks and
passes identically under `python` and `python -O`. QED.
