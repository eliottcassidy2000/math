# AP8 with arbitrary high outliers: the optimal limiting margin

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
This is a routine analytic extension of the proved
[phase-packet supplier](lrc14_phase_packet_supplier_overnight_hexagon_sep05.md),
retaining the previously unpaid outlier branch divisible by nine. No theorem
ID or broad novelty claim.

## 1. Exact statement and the inherited branch

Let `v,w>8` be distinct positive integers, let `H=min(v,w)>=82`, and set

```text
C={1,2,...,8,v,w},
S=3C union T,
M(S)=max_(x mod1) min_(s in S) ||sx||,
L_lambda(S)=measure{x mod1: ||sx||>=lambda for every s in S}.
```

Here `T` consists of any three distinct positive integers not divisible by
three. Consequently the thirteen speeds are distinct; body/tail overlap
is impossible modulo three. There is no parity assumption and no bound
on any tail height.

Then

```text
1/9-9/H <= M(S) <=1/9.                                   (1)
```

More precisely,

```text
M(S)=1/9  iff  9 divides neither v nor w.                 (2)
```

The exact equivalence (2) holds for all distinct `v,w>8`; the assumption
`H>=82` is needed only for the positive quantitative margin below.

When at least one outlier is divisible by nine, put
`lambda_H=1/9-9/H>0`. There is also the tail-height-independent measure
bound

```text
L_(lambda_H)(S) >= 53/(4536H)>0.                          (3)
```

If both outliers are divisible by three, the right side of (3) improves
to `4507/(4536H)`. These are sufficient explicit measure floors, not sharp
measure constants. The limiting margin `1/9` in (1) is optimal and holds
uniformly over all three tail heights as both outlier heights tend to
infinity. Even in the strict branch of (2), the margins converge to `1/9`.

The point statement on the non-nine-divisible branch is already inherited
from the clock-27 degeneration mechanism in
[THM-492 — cprime14-two-head-reduction-not-exact-band-ladder](../../01-canon/theorems/THM-492-cprime14-two-head-reduction-not-exact-band-ladder.md).
Section 2 spells out the elementary unit-clock occupancy needed at the
endpoint `1/9`. That old clock does not handle the other branch, because
an outlier divisible by nine makes its physical speed divisible by 27.

The retained board is: AP8 endpoint rigidity; third-orbit body packets;
outlier valuation; one-comb discrepancy; physical safe measure. The map
from a body phase to its nine physical points keeps body safety and
replaces tail addresses by a uniform occupancy bound. Its sidecar is the
complete indexed packet and its factor-three integration normalization.
Outlier residue data are kept until the exact clock branch is separated.
The least-used operation is allowing the margin to approach the AP8
endpoint while retaining a strictly positive outlier density budget.

## 2. The AP8 ceiling and its complete equality set

For a body phase `y`, consider the nine circle points
`0,y,2y,...,8y`. Two have circle distance at most `1/9` by ordering them
and considering the nine cyclic gaps. Their difference is `k y` for some
`1<=k<=8`. Thus `min_(1<=k<=8)||k y||<=1/9` for every `y`, proving the
upper bound in (1) after setting `y=3x`.

Suppose all eight distances are at least `1/9`. The nine points are then
distinct, and every pair has circle distance at least `1/9`, because its
difference is again `k y` with `1<=|k|<=8`. In particular each cyclic gap
is at least `1/9`. Since the nine gaps sum to one, each equals `1/9`.
The point zero belongs to the configuration, so the configuration is
exactly `{0,1/9,...,8/9}`. Hence `y=a/9`; distinctness of all nine
multiples requires `gcd(a,9)=1`. Conversely each such phase has minimum
distance exactly `1/9`. The complete physical equality set is therefore

```text
x=u/27,       gcd(u,27)=1.                                (4)
```

There are eighteen such points modulo one. If `9` divides neither
outlier, every body speed `3c` is safe at all eighteen: `c` is nonzero
modulo nine, and multiplication by a unit preserves that property.
Each three-unit tail permutes the eighteen unit residues modulo 27.
At margin `1/9`, precisely four unit residues are strictly dangerous,
namely `+/-1,+/-2 mod27`. Thus each tail removes four points, and three
tails remove at most twelve. At least six points remain safe for the
whole row, proving the forward equality assertion in (2).

If `9|v` or `9|w`, its physical speed is divisible by 27 and vanishes at
every point (4). Consequently no full-row phase can attain `1/9`.
The finite minimum of continuous distance functions is continuous on
the compact circle and attains its maximum. Therefore `M(S)<1/9` in
this branch, completing (2), not merely showing failure of one clock.

## 3. Supplier at a variable margin

For an outlier `z`, set `omega(z)=1` if `3|z` and `omega(z)=3` otherwise.
The inherited supplier proves that, for `0<lambda<=1/9`, the AP8 core
has a common-third-orbit-safe interval

```text
J_lambda=[(2+3lambda)/21,(1-lambda)/8],
ell=|J_lambda|=5(1-9lambda)/168.
```

Its twenty-four addressed tooth inequalities are affine in `lambda`
and valid at both endpoints `lambda=0,1/9`; the exact controls retain
these whole intervals, not just sampled body phases. The third-orbit
danger comb of an outlier has duty `2lambda omega` and one-period
boundary discrepancy at most

```text
f_omega(lambda)/z,
f_omega(lambda)=2lambda(1-2omega lambda).
```

Subtract the two danger measures from `J_lambda`. If

```text
b=ell[1-2lambda(omega(v)+omega(w))]
     -f_(omega(v))(lambda)/v-f_(omega(w))(lambda)/w >0,    (5)
```

the residual has measure at least `b`. Its three translates by thirds
are disjoint and belong to the complete body-safe third-orbit
intersection, because `ell<1/3`. They supply at least `3b` body-phase
measure. The inherited nine-point packet then gives
`L_lambda(S)>=b`: each three-unit tail removes at most two of the nine
points, and the physical integration normalization is exactly the
factor `1/3` in Section 3 of the supplier. This deduction has no tail
height restriction.

## 4. The high-outlier budget

Set `delta=1/9-lambda`, so `ell=15delta/56`. If either outlier is
divisible by nine, at least one is divisible by three, leaving only
the two cases below. On the full interval `0<=lambda<=1/9`,

```text
f_1(lambda)<=14/81,
f_3(lambda)<=1/12.                                       (6)
```

The first function is increasing on this interval; the second reaches
its maximum at `lambda=1/12`.

For one three-divisible outlier and one three-unit outlier, the sum of
orbit multiplicities is four. Its density factor is at least `1/9`, so

```text
b >=5delta/168-83/(324H).
```

For two three-divisible outliers, the sum of multiplicities is two.
Its density factor is at least `5/9`, so

```text
b >=25delta/168-28/(81H).
```

At `delta=9/H`, both are strictly positive:

```text
mixed:  b >= [45/168-83/324]/H =53/(4536H),
both:   b >= [225/168-28/81]/H=4507/(4536H).              (7)
```

The hypothesis `H>=82` guarantees `0<lambda_H<1/9`, so all supplier
and discrepancy hypotheses apply. Equations (5)--(7) prove (3) and
the lower bound in (1) on the nine-divisible branch. The other branch
has already been closed exactly by Section 2. The upper bound and
these uniform lower bounds prove the stated optimal limiting margin.

## 5. Audit and bounded exact controls

[Source](../../04-computation/lrc14_ap8_asymptotic_empty_core_morning_sep06.py)
and [output](lrc14_ap8_asymptotic_empty_core_morning_sep06.out).

```bash
python3 -B 04-computation/lrc14_ap8_asymptotic_empty_core_morning_sep06.py
python3 -B -O 04-computation/lrc14_ap8_asymptotic_empty_core_morning_sep06.py
```

All **273 explicit gates** pass, with identical normal and optimized
output. The declared universe comprises 24 AP8 addressed tooth cells,
both endpoint inequalities at each of the two margin endpoints; all
eighteen unit tail residues modulo 27; 64 rational one-comb discrepancy
controls; and six named high-outlier budgets, including `H=82`.
The discrepancy controls independently compare the reduced primitive
comb with the literal union of the three observers `v(y+j/3)`, retaining
duplicates when `3|v` until the union operation. Both discrepancy maxima
and the final positive rational constants are checked as polynomial or
rational identities. No tail-height census is used or inferred.

The `three_ray_geometry` agent independently audited the root's supplier
substitution, all constants, quantifiers, and clock-27 overlap before
writing the corollary. The root then reread the complete Sections 1--4
and independently checked the pigeonhole equality rigidity, all eighteen
physical endpoint phases, the strict branch via compactness, the two
outlier budgets, and the physical measure factor: **PASS**. The source
uses only the Python standard library and exact `Fraction` arithmetic;
its finite checks support the analytic proof rather than replacing it.

```text
source SHA256 905d0879cb8be6e887b1617beb6a670087be30d20f78865edb14464358c2ccad
output SHA256 b7176f31b41a8c844e0c6989a8b738a391b49f8fc38be43706ef8a09491e5b6b
trace  SHA256 b77950345b4a64464a77c89f257362366f94913e33384dcace506977d6b8e151
```
