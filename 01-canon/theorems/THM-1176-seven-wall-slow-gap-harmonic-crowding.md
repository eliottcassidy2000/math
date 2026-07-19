---
id: THM-1176
title: Seven-wall slow-gap harmonic crowding and the finite toothpick ladder beyond the zero-density wall
status: PROVED analytic implication for every actual LRC(14) counterexample and every six-plus-seven partition.  If a is the slowest of the seven covering combs and m is the retained-core maximum, then either a<13m or a*sum(1/b_i)>1; in the latter branch distinctness sharpens b_1<6a to b_1<=6a-3.  More generally, r faster combs covering a slow gap force c*sum(1/d_i)>7-r, so r<=3 is impossible and recursive gap nesting stops after two steps with exact r=5 and r=4 integer cutoffs.  The normalized phase clock and equality-parity obstruction are exact.  This does not prove the universal six-comb slow-gap noncoverage conjecture or LRC(14)
source: codex-2026-07-18-S75 slow-gap carrier attack
depends_on: [THM-1094, THM-1153]
related: [THM-1015, THM-1123, THM-1166, HYP-7678]
script: 04-computation/lrc14_seven_wall_slow_gap_harmonic_crowding_codex_20260718.py
output: 05-knowledge/results/lrc14_seven_wall_slow_gap_harmonic_crowding_codex_20260718.out
---

# THM-1176 -- seven-wall slow-gap harmonic crowding

At the actual LRC(14) radius, write

```text
D_s={t in R: ||st||<1/14}.                              (1)
```

Let the seven deleted speeds be

```text
a<b_1<...<b_6.                                         (2)
```

The speed `a` is distinguished only as the slowest comb.  Its consecutive
danger teeth leave the closed safe gaps

```text
G_k=[(14k+1)/(14a),(14k+13)/(14a)],   k in Z,          (3)
|G_k|=6/(7a).                                          (4)
```

The central observation is that the other six combs have exactly enough bulk
duty to fill `6/7` of such a gap.  Their endpoint discrepancies can fill the
last `1/7` only if they carry one full unit of normalized reciprocal pressure.

## 1. Exact slow-gap survivor lemma, at every cardinality

> **Lemma 1.** For every positive integer packet (2), every integer `k`, and
> every `1<=r<=6`,
>
> ```text
> |G_k minus union_(i=1)^r D_(b_i)|
>   >=(6/49)((7-r)/a-sum_(i=1)^r 1/b_i).               (5)
> ```

**Proof.**  For an interval `J` of length `L`, the sharp one-comb periodic
discrepancy from THM-1094 is

```text
|J intersect D_b|<=L/7+6/(49b).                        (6)
```

For completeness, scale `J` by `b` and write `bL=N+s`, with integer `N` and
`0<=s<1`.  The `N` full periods contribute `N/7`, while the residual interval
meets the danger arc of length `1/7` in at most `min(s,1/7)`.  Hence its excess
over mean is at most

```text
max_(0<=s<=1)(min(s,1/7)-s/7)=6/49,                    (7)
```

which proves (6).  Apply (6) to the first `r` speeds on `G_k` and use the
union bound.  Since `|G_k|=6/(7a)`,

```text
|G_k minus union_(i=1)^r D_(b_i)|
 >= |G_k|-r|G_k|/7-(6/49)sum_(i=1)^r 1/b_i
  = (6/49)((7-r)/a-sum_(i=1)^r 1/b_i).                 (8)
```

This proves the lemma. ∎

In particular,

```text
sum_i 1/b_i<1/a                                        (9)
```

forbids the other six combs from covering **any** slow gap, uniformly in its
phase `k`.

There is also a component form.  Put

```text
U=D_a union D_(b_1) union ... union D_(b_6).            (10)
```

Under (9), no connected component of the open set `U` can join two
consecutive teeth of `D_a`, because doing so would cover the whole intervening
gap `G_k`.  A component therefore meets at most one `a`-tooth.  It can use at
most that tooth, of length `1/(7a)`, and portions of its two adjacent safe
gaps, each of length `6/(7a)`.  Consequently every component obeys

```text
component_length<13/(7a).                              (11)
```

The constant `13` is the exact carrier span `gap + period`:

```text
6/(7a)+1/a=13/(7a).                                    (12)
```

## 2. Strict LRC(14) dichotomy

Consider an actual LRC(14) counterexample `V`, so `V` consists of thirteen
distinct positive integer speeds and

```text
M(V)<1/14.                                              (13)
```

Partition it as `V=W disjoint_union S`, where `|W|=6`, `|S|=7`, let
`m=max W`, and order `S` as in (2).  The known six-speed lonely-runner bound
gives `M(W)>=1/7`.  At a maximizing time for `W`, the `m`-Lipschitz property of
`min_(w in W)||wt||` therefore supplies a closed interval `I` of length

```text
|I|=1/(7m)                                              (14)
```

on which every speed in `W` is safe at radius `1/14`.  By (13), the seven
combs in `S` cover all of `I`.

> **Theorem 2 (strict slow-gap crown).** Every such partition satisfies
>
> ```text
> a<13m
>       or
> a sum_(i=1)^6 1/b_i>1.                                (15)
> ```
>
> In the second branch,
>
> ```text
> b_1<=6a-3<6a.                                        (16)
> ```

**Proof.**  Suppose `a>=13m`.  Then (14) gives

```text
|I|>=13/(7a)=1/a+6/(7a).                               (17)
```

An interval of length at least one `a`-period plus one slow-gap length
contains some complete `G_k`: take the first slow-gap left endpoint at or
after the interval's left endpoint.  Since `D_a` is absent from the interior
of `G_k` and its strict endpoint convention does not cover the endpoints,
the other six combs cover this `G_k`.  Lemma 1 first gives

```text
sum_i 1/b_i>=1/a.                                      (18)
```

Section 3 proves that equality in (18) is impossible for a covered slow gap.
Thus the inequality is strict, proving (15).

Finally, distinct order first gives

```text
sum_i 1/b_i<6/b_1.                                     (19)
```

Combining this with the second branch of (15) gives `b_1<6a`.  Integrality
and distinctness give the sharper claim.  If instead `b_1>=6a-2`, then

```text
sum_i 1/b_i<=sum_(r=-2)^3 1/(6a+r)<1/a.               (19a)
```

For the strict inequality put `x=6a>=6` and clear the positive denominator
`x(x-2)(x-1)(x+1)(x+2)(x+3)`.  The numerator of
`6/x-sum_(r=-2)^3 1/(x+r)` is

```text
3(x-6)^4+62(x-6)^3+423(x-6)^2+988(x-6)+264>0.         (19b)
```

This contradicts the crowded branch and proves (16). ∎

For the top-seven deletion of ordered speeds `v_1<...<v_13`, Theorem 2 reads

```text
v_7<13v_6
  or
v_7 sum_(i=8)^13 1/v_i>1,                              (20)
```

and the crowded branch has `v_8<=6v_7-3`.  Thus the zero coefficient at `r=7`
in THM-1153 is not an absence of structure: beyond the carrier threshold
`v_7>=13v_6`, its former lower bound `0` becomes the strict bound `1`.
The finite nesting argument in Section 5 also gives, in that separated
branch,

```text
v_8/v_7<13/6  or  v_9/v_8<13/6  or  v_10/v_9<4/3.     (20a)
```

## 3. Equality rigidity at every cardinality

Strictness already has a short topological proof.  At equality every
one-comb discrepancy bound is saturated and the relatively open sets
`D_(d_i) intersect int(G_k)` are pairwise disjoint almost everywhere.  Two
overlapping open intervals have positive overlap, so these sets are actually
disjoint; a finite nontrivial disjoint family cannot cover the connected
interval `int(G_k)`.  The arithmetic proof below is retained because it gives
the stronger phase rigidity needed by the clock/owner viewpoint, not because
connectedness alone is insufficient.

More generally, let `1<=r<=6`, let `c<d_1<...<d_r`, and suppose the `r`
faster combs cover one complete `c`-slow gap.  Lemma 1 gives

```text
c sum_(i=1)^r 1/d_i>=7-r.                              (21)
```

Suppose equality holds.  The upper bounds in the proof of Lemma 1 then sum
to exactly the slow-gap length.  Since their union also has that measure,
all `r` one-comb discrepancy bounds must be equalities and the intersections
must be pairwise disjoint almost everywhere.

For a fixed `d=d_i`, the scaled gap has length

```text
d|G_k|=6d/(7c).                                        (22)
```

Equality in (7) occurs only when the fractional part in (22) is `1/7`.
Therefore there is an integer `h_i>=1` such that

```text
6d_i=c(7h_i+1).                                        (23)
```

Length saturation alone is not enough.  The residual interval of length
`1/7` must coincide with the danger arc.  Since the left endpoint of the
`c`-slow gap is `(14k+1)/(14c)`, this phase alignment is exactly

```text
d_i(14k+1)+c == 0  (mod 14c).                          (24)
```

Substituting (23) into (24) gives

```text
(7h_i+1)(14k+1)+6 == 0  (mod 84),                      (25)
```

or, after division by `7` and reduction modulo `12`,

```text
12 divides (h_i+1)(2k+1).                              (26)
```

The second factor is odd, so `4` divides `h_i+1`; hence

```text
h_i==3 (mod 4).                                        (27)
```

Write

```text
e_i=(7h_i+1)/2.                                        (28)
```

Every `e_i` is odd.  Equations (21) and (23) imply

```text
3 sum_(i=1)^r 1/e_i=7-r.                               (29)
```

Let `P=product_i e_i`, which is odd.  Clearing denominators in (29) gives

```text
3 sum_i P/e_i=(7-r)P.                                  (30)
```

Each `P/e_i` is odd.  The parity of the left side is therefore `r mod 2`,
whereas the parity of the right side is `(7-r) mod 2`.  Those parities are
opposite.  This contradiction proves the strict all-cardinality law

```text
c sum_(i=1)^r 1/d_i>7-r.                               (30a)
```

For `r=6`, (30a) is exactly the strictness needed in (18).

This also explains why a scalar mass calculation alone stops at a non-strict
inequality.  Strictness appears only after the endpoint phase is retained.

## 4. Exact normalized phase clock

Put

```text
x=at-k,              x in K=[1/14,13/14],              (31)
b_i=q_i a+r_i,        0<=r_i<a.                        (32)
```

Modulo one,

```text
b_i t == (q_i+r_i/a)x+k r_i/a.                         (33)
```

Thus `G_k` is covered by the six `b_i` exactly when

```text
K subset union_i {x: ||(q_i+r_i/a)x+k r_i/a||<1/14}.   (34)
```

Let

```text
g=gcd(a,b_1,...,b_6)=gcd(a,r_1,...,r_6),
T=a/g.                                                  (35)
```

The joint phase tuple in (33) has fundamental clock period `T`: `T r_i/a` is
an integer for every `i`, and the prime-adic gcd condition shows no smaller
positive integer advances every coordinate integrally.  Consequently the
binary cover-status word in `k` has period `T` (its least period may be a
proper divisor through accidental symmetry).

This is an exact finite reduction of the phase variable for a fixed speed
packet.  It is not a finite reduction of the speed ratios `q_i`.

For exact scalar replay, define

```text
F(y)=floor(y)/7+min(frac(y),1/7).                       (36)
```

If `A,B` are the endpoints of `G_k`, then, up to null endpoint conventions,

```text
|G_k intersect D_b|
  =[F(bB+1/14)-F(bA+1/14)]/b.                          (37)
```

Equations (32)--(37) are the quotient/residue/clock coordinates needed by an
exact slow-gap search.  A ratio-only or gcd-only classifier discards the
phase term `k r_i/a` and therefore cannot decide coverage.

## 5. Harmonic drift and the finite toothpick recursion

The all-cardinality law already has strong order consequences.  Since
`d_i>=d_1+i-1`, (30a) implies

```text
d_1/c<r/(7-r).                                         (37a)
```

For `r<=3` the left side of (30a) is strictly less than `r<=7-r`, so three
or fewer faster combs can never cover a complete slow gap.

Distinctness gives exact integer cutoffs at the two recursive rungs used
below.  Five faster combs covering a `c`-slow gap force

```text
d_1<=floor((5c-4)/2).                                  (37b)
```

Indeed, at the first excluded value
`d=floor((5c-4)/2)+1`, the consecutive envelope
`c sum_(j=0)^4 1/(d+j)` is already below `2`, and it decreases with `d`.
For `c=2u`, `u>=1`, the exact deficit is

```text
2-c sum_(j=0)^4 1/(5u-1+j)
 =2(625u^3+250u^2-75u-24)
  /[5(5u-1)(5u+1)(5u+2)(5u+3)]>0.                    (37c)
```

For `c=2u+1`, `u>=1`, the numerator over a positive denominator is

```text
625u^4+1000u^3+375u^2-58u-34>0.                       (37d)
```

The omitted case `c=1` has maximal envelope
`sum_(d=2)^6 1/d=29/20<2`.

Likewise, four faster combs covering a `c`-slow gap force

```text
d_1<=floor((8c-9)/6).                                  (37e)
```

At the first excluded value, write `c=6u+s`.  For `s=0,...,5` the first
excluded speeds are `8u+(-1,0,2,3,4,6)`.  The corresponding positive
cleared numerators of `3-c sum_(j=0)^3 1/(d+j)` are

```text
3(64u^2-8u-3),
128u^3-48u^2-34u-3,
320u^3+360u^2+125u+13,
3(64u^2+56u+9),
64u^3+72u^2+7u-8,
640u^3+1680u^2+1450u+411.                              (37f)
```

They are positive whenever that excluded speed is greater than `c`.  The
remaining `c=1,2,3,4` consecutive envelopes are respectively
`77/60,19/10,319/140,533/210`, all below `3`.

Now return to a six-comb cover of `G_k`.  If `b_1>=13a/6`, the length
`6/(7a)` of `G_k` is at least the carrier span `13/(7b_1)`, so it contains a
complete `b_1`-slow gap.  Both `D_a` and `D_(b_1)` are absent there, hence
`b_2,...,b_6` cover it and (37b) gives

```text
b_2<=floor((5b_1-4)/2).                                (37g)
```

If also `b_2>=13b_1/6`, that gap contains a complete `b_2`-slow gap;
`b_3,...,b_6` cover it and (37e) gives

```text
b_3<=floor((8b_2-9)/6)<(4/3)b_2<(13/6)b_2.             (37h)
```

Thus the self-similar nesting stops after two steps.  Every covered slow gap
obeys the concise ladder alternative

```text
b_1/a<13/6  or  b_2/b_1<13/6  or  b_3/b_2<4/3.        (37i)
```

The natural state along the ordered six-comb ladder is cumulative normalized
pressure

```text
H_j=a sum_(i=1)^j 1/b_i,       H_0=0.                  (38)
```

The crowded branch is precisely `H_6>1`.  If

```text
R_j=1-H_(j-1)>0,                                       (39)
```

then the remaining `7-j` teeth must cross that deficit.  Order gives

```text
R_j<a sum_(i=j)^6 1/b_i<=(7-j)a/b_j,                  (40)
```

with the conclusion

```text
b_j/a<(7-j)/R_j.                                       (41)
```

For `j=1`, this is `b_1/a<6`; (19a)--(19b) record the extra integral
three-unit saving.  Each subsequent tooth is sized by the entire deficit
left by the preceding stalk.  This is the slow-gap version of
THM-1153's H-drift/toothpick recursion, but it begins only after the geometric
carrier split `a>=13m`.

There is an all-`n` self-similar form.  At radius `1/(2n)`, let `a` be the
slow carrier and let `n-1` other danger combs have speeds `b_i`.  The slow safe
gap has length `(n-1)/(na)`, while the sharp one-comb endpoint excess is
`(n-1)/(n^2 b)`.  The same proof gives

```text
|slow gap minus union_i D_(b_i)|
 >=(n-1)/n^2 (1/a-sum_i 1/b_i).                        (42)
```

The LRC(14) coefficient `6/49` is the `n=7` tooth in this all-`n` law.

## 6. Carrier and tournament audit

The usual runner-order tournament is not the proof object here.  After
choosing `a`, the six remaining combs have the pairwise observable

```text
O(i,j)=a/b_i-a/b_j.                                    (43)
```

Use the positive sign as the switch/gauge and break a hypothetical tie by
the speed label.  Since the speeds are distinct, the resulting tournament is
transitive: score histogram `0,1,2,3,4,5`, no directed cycles, six singleton
SCCs, and one Hamiltonian path, namely speed order.  It preserves the order
needed in (19) and (40), but it destroys the magnitudes whose sum is `H_6`,
all residues `r_i`, every tooth overlap, and the clock `k`.  Tournament
fingerprints therefore cannot decide (34).

The challenged assumption is that vertices must be runners or arcs.  A more
faithful vertex set is the cyclic family of slow gaps

```text
Z/TZ,                                                   (44)
```

with six labelled tooth-interval families drawn over each vertex.  Equivalently,
one may use the ordered endpoint-crossing events inside `K`.  This carrier
preserves the exact local predicate "the six noncarrier combs cover `G_k`"
and its chronology.  Quotienting it to an owner tournament destroys interval
length and multiple ownership, so the correct object is a cyclic interval
hypergraph or owner word, not a bare tournament.

This locates the Fano/`chi_7` probe more sharply.  It is not needed in the
strict subcritical harmonic cone (9).  It can matter only in the crowded cone
`H_6>1`, where pair, triple, and endpoint incidence must spend the excess
pressure.  The Fano lines should therefore be laid over gap/endpoint events,
not imposed as a fixed runner tournament.

## 7. Verification and remaining frontier

The companion `Fraction` replay verifies the constants in (5), (12), the
all-cardinality parity split, and both exact integer cutoffs in (37b) and
(37e) through `c=1000`.  As non-proof telemetry, it also exhausts every
packet

```text
3<=a<=7,       a<b_1<...<b_6<=3a,       0<=k<a,        (45)
```

using exact rational closed teeth.  Closing teeth only makes coverage easier.
Across `4,166` packets and `27,730` slow gaps it finds zero closed-cover
candidates.  It also checks several larger consecutive and common-dilate
packets.  These rows support, but do not prove, the stronger conjecture that
six distinct integer combs never cover a full slow gap.

What is proved is the strict disjunction (15), its integral ratio consequence
(16), the all-cardinality law (30a), the finite nesting ladder (37i), the
exact phase clock, and the equality obstruction.  What remains open is
the harmonic-crowded region

```text
a sum_i 1/b_i>1.                                       (46)
```

Closing that region would yield the universal slow-gap theorem and force
`a<13m` in every seven-wall packet.  THM-1166's Fano/gcd inequalities and
adaptive forests are complementary constraints inside (46); neither theorem
currently excludes all mixed-residue phase clocks.
