# Effective clock repair and the ten-subset gcd boundary

**Status: PROVED RELATIVE TO CITED LRC THROUGH THIRTEEN RUNNERS;
INDEPENDENTLY AUDITED; FINITE-EXACT controls.** This recovers
the gcd-sensitive clock mechanism, repairs a false specialization in
THM-737, and combines deletion inheritance with exact orbit orders.

For a primitive set `V` of thirteen distinct positive integer speeds,
write `M(V)=max_x min_(v in V)||vx||`. If `M(V)<1/14`, then:

* every twelve-subset has gcd one;
* every eleven-subset has gcd at most two;
* every ten-subset has gcd at most four.

For a ten-subset `P` with gcd `c>1`, let `w_1,w_2,w_3` be its complement
and sort `q_i=c/gcd(c,w_i)`. The only remaining scalar/deletion profiles are

| `c` | sorted effective orders |
|---|---|
| 2 | `(1,2,2)` or `(2,2,2)` |
| 3 | `(3,3,3)` |
| 4 | `(2,4,4)` |

These are necessary conditions on a hypothetical **strict counterexample**,
not realizations of unsafe rows. In particular, every primitive thirteen-speed
row with ten speeds of gcd greater than four is `1/14`-safe, at arbitrary
heights. Without primitivity the equivalent sufficient condition is
`gcd(P)>4 gcd(V)`, by common dilation. This does not prove LRC(14).

## 1. Inheritance and exact recovery

The closest mechanisms are
[THM-765 — hereditary primitivity, Part B](../../01-canon/theorems/THM-765-safe-component-tooth-deck-and-hereditary-primitivity.md),
[THM-761 — multi-exception sheet bound](../../01-canon/theorems/THM-761-multi-exception-sheet-covering-bound.md),
and [THM-4004 — three-detuned divisor combs](../../01-canon/theorems/THM-4004-lrc14-three-detuned-divisor-comb-profile.md).
The last already gives the exact open-arc count and prime-incidence cuts in
its typed `11+2` chart. [THM-860](../../01-canon/theorems/THM-860-primitive-hamming-six-finite-ramification-reduction.md)
already uses hereditary pair/complement lcms on a different, six-colour
AP-centred carrier at threshold `1/13`. Neither those old operations nor
their general use of gcd orders is a novelty claim here.

The corrected near miss is the coprime specialization in
[THM-737](../../01-canon/theorems/THM-737-pack-clock-sampling-measure-dispatch.md).
The canonical hostile below is a dilation of its own sharp tower. The
least-used sidecar is the actual orbit length `c/gcd(c,w)`, including the
literal branch multiplicity. The board is: common clocks; effective orders;
leave-one-out gcds; exact safe components; scalar-versus-joint coverage.

Targeted searches recovered these results, the old bounded
`lrc14_three_detuned_exceptional_kps_S127.out` order bank, and the coded
deletion interpretation of THM-2069. No prior exact global ten-subset-gcd-four
statement was recovered in that search. The elementary assembly below is
recorded as repository progress without an external-priority claim.

The only external input is LRC for up to twelve nonzero speeds. The
[current primary preprint, Theorem 1.3](https://arxiv.org/html/2604.23906v2)
states this result; v2 was posted September 1, 2026. Its Definition 2.1 also
retains leave-one-out gcds in the properness sieve. The theorem is cited,
not reproved here. `CORE-PAPERS.md` was consulted before checking the
primary abstract and theorem text.

## 2. Repair the scalar count before using it

At a body-safe phase `y` for `C`, all `c` lifts `(y+j)/c` are safe for
`cC`. For one exceptional speed `w`, let `g=gcd(c,w)` and `q=c/g`.
Its phases form an equally spaced order-`q` orbit, each point appearing
exactly `g` times. An open arc of length `1/7` contains at most
`ceil(q/7)` such points. This remains true when `7|q`: the open endpoints
give exactly the integer bound, without a further `+1`. Put

```text
beta(q)=ceil(q/7)/q.
```

The exact uniform branch budget is therefore

```text
# bad branches for w <= c beta(q),
mu(G_(cC union D)) >= mu(G_C) max(0,1-sum_(w in D) beta(c/gcd(c,w))). (1)
```

If the sum is strictly below one, every safe body phase has a free lift.
The coprime formula substitutes `q=c`; it cannot be applied to other gcd
strata. In particular the scalar cutoff `c>7d/(7-d)` is justified by the
coarser estimate only when every exceptional speed is coprime to `c`.
The original gcd-aware inequality in THM-737 stays valid.

For its exact hostile, take `C={1,...,12}`, `c=4`, `D={26}`. The exceptional
speed is not divisible by four, but its orbit order is two. Exact measures
are

```text
mu(G_C)=6617/194040,
mu(G_(4C union {26}))=mu(G_(2C union {13}))=6617/388080.
```

The falsely unqualified coprime floor is `3 mu(G_C)/4=6617/258720`, larger
than the actual measure. Formula (1) instead gives the attained floor
`mu(G_C)/2`. The first failed implication was removing `gcd(c,w)` after
counting the orbit. Nonmultiple detuning is weaker than coprimality.
This refutes the claimed general measure floor, not the actual row's
loneliness or the sound gcd-sensitive theorem.

## 3. Hereditary primitivity and the eleven-subset bound

Assume now that `V` is primitive with `M(V)<1/14`. For any twelve-subset
`A`, the cited lower-runner theorem gives `M(A)>=1/13`. Write its missing
speed as `w`. THM-765(B) says a single speed can destroy this margin only
if `gcd(A)|w`, since the target is below `1/4`. That gcd divides all of
`V`, so it is one.

For an eleven-subset `A` with gcd `c>1`, its two missing speeds are each
coprime to `c`: adjoining either gives a twelve-subset of gcd one. Their
two effective orders are both `c`. The divided body has a nonempty safe
set, so (1) forces `2 beta(c)>=1`. But `beta(c)<=1/2` for `c>=2`, with
equality only at `c=2`. Hence `c<=2`.

The elementary bounds used here and below are

```text
beta(q)<=1/2 for q>=2, equality only at q=2;
beta(q)<=1/3 for q>=3, equality only at q=3;
beta(q)<=1/4 for q>=4, equality only at q=4 or8.       (2)
```

Check `q=2,...,8` directly. For `q>=9`,
`beta(q)<= (q+6)/(7q)<1/4`, proving all tails of (2).

## 4. Ten-subsets and the complete four-profile residual

Let `P` be a ten-subset with gcd `c>1`. For each complementary speed,
`gcd(c,w_i)` is the gcd of an eleven-subset and is thus at most two.
Consequently

```text
q_i=c/gcd(c,w_i) >= c/2.                             (3)
```

Every twelve-subset has gcd one, giving the inherited pair-lcm condition

```text
lcm(q_i,q_j)=c for every pair i!=j.                  (4)
```

Indeed `lcm(c/g_i,c/g_j)=c/gcd(g_i,g_j)` and the last gcd equals
`gcd(P,w_i,w_j)=1`. Full failure also forces

```text
beta(q_1)+beta(q_2)+beta(q_3)>=1.                    (5)
```

If `c>=5`, (3) makes every integer `q_i>=3`. By (2), (5) forces all
three orders to equal three. Equation (4) would then give `c=3`, a
contradiction. This proves the unbounded gcd cap without any finite-height
extrapolation.

The table follows directly at the remaining three clocks. At `c=2`, (4)
permits at most one order-one tail. At `c=3`, (3) forces all orders three.
At `c=4`, the orders are two or four, and (4) permits at most one two.
All fours have budget `3/4<1`, leaving exactly `(2,4,4)`. The surviving
budgets are respectively `2`, `3/2`, `1`, and `1`; a budget of one alone
does not provide a free sheet.

At the order-one profile, the extra even tail is absorbed into an eleven-body
at clock two. The `(2,4,4)` profile also becomes a clock-two eleven-body
with two odd tails after absorbing its uniquely even exceptional speed.
Thus the phase problems are the existing clock-two ten/eleven-body and
clock-three ten-body branches. The table bounds clocks, not body heights,
tail heights, or the number of component-address states.

Incoming [THM-4442](../../01-canon/theorems/THM-4442-lrc14-bounded-ten-body-parity-free-scale-three-completion.md)
now closes the clock-three branch when the divided body is a ten-subset
of `[13]`. It does not close the arbitrary ten-bodies left here. The
original empty-hexagon paper motivates finding a smaller legal certificate;
the actual cross-lane map here is gcd inheritance, not a polygon analogy.

### Scale-three entry corollary and the concurrent reserved target

If `V=3C union T` is a primitive set of thirteen distinct positive speeds,
`|C|=10`, `|T|=3`, and `gcd(C)>1`, then `M(V)>=1/14`. Indeed the ten-pack
has gcd `3gcd(C)>=6>4`, so the proved cap applies directly. Thus any strict
counterexample already in this scale-three chart must have `gcd(C)=1`;
no bound on the entries of C or T is used. More generally a primitive row
`3dC union T` is safe for every integer `d>=2`.

Incoming [THM-4446, primitive ten-pack descent and dilation rays](../../01-canon/theorems/THM-4446-lrc14-primitive-ten-pack-descent-and-dilation-rays.md)
is still **RESERVED / UNPROVED EMPTY STUB** at the checked revision
`bf5a1cf355`. Its intended nonprimitive-ten-pack target is a direct consumer
of the theorem here. This paragraph proves the stated corollary independently;
it does not use or promote the reservation, nor assert arbitrary chart entry.
The independent clock-certificate referee checked both corollaries and the
reservation boundary after integration: **PASS**.

## 5. Exact controls and scope

The [standalone verifier](../../04-computation/lrc14_effective_clock_empty_core_sep06.py)
imports no repository geometry code. Danger-union measure and endpoint-cell
classification independently reconstruct the three named measures. Every
open-arc event and intervening cell is tested for `q=1,...,70`, including
the equality boundaries `7|q`. Every sorted divisor-order triple for
`2<=c<=420` is checked against the eleven-subset bound, all three pair
lcms and (5); exactly the four displayed profiles remain. This finite
control corroborates the proof in Sections 3--4; it is not its cutoff.

A separate phase-only hostile takes `C={1+8j:0<=j<=9}`, `y=3/8`, and
odd tails `(3,5,7)` at clock two. The body is safe, but the first two tails
kill opposite sheets, with masks `(2,1,0)`. This blocks declaring a
phase-uniform selector in the residual. It is not claimed to be an unsafe
full speed set.

```sh
python3 -B 04-computation/lrc14_effective_clock_empty_core_sep06.py
python3 -B -O 04-computation/lrc14_effective_clock_empty_core_sep06.py
```

The [frozen output](lrc14_effective_clock_empty_core_sep06.out) has 147
exception-raising gates; normal and optimized outputs agree. The source
and output SHA-256 hashes are
`12fde59aab96da6a312b358d058b7261834a4515a5497a64924a9ffa71aabecd`
and `f734022c70d8632417c1a60f8479243a7ad579c8a1b16dbaa501911d8b67e9b5`. The next decisive object
is the remaining phase-labelled component union, with the true effective
orders and the actual entry restrictions retained.

Independent final proof and source audit: **PASS**. The [audit sidecar](lrc14_effective_clock_audit_empty_core_sep06.md) checks the current primary lower-runner theorem, every inherited hypothesis, strict budget boundary, residual profile, common dilation, and the 147-gate replay.
