---
id: THM-1212
title: Nested first-safe windows close five and six carriers and complete the r=6 clustered stratum
status: PROVED analytically.  THM-1169's density residuals feed a nested sequence of closed first-safe windows.  Five carriers close in one step.  Six carriers force the exact ratio ladder y/x<35/12, z/x<1960/363<27/5, w/x<567/76, v/x<15876/2105<11, after which the final carrier has an open tooth strictly shorter than the protected closed interval.  Together with THM-1136, THM-1154, and THM-1169, every zero-through-six-carrier r=6 row is lonely
source: codex-2026-07-18-S77
depends_on: [THM-1136, THM-1154, THM-1155, THM-1169]
related: [THM-1135, THM-1151, THM-1137]
script: 04-computation/lrc14_r6_five_six_carrier_nested_window_referee_codex_S77.py
output: 05-knowledge/results/lrc14_r6_five_six_carrier_nested_window_referee_codex_S77.out
lean: 04-computation/lean/TournamentH7/TournamentH7/LRCNestedCarrierWindow.lean
---

# THM-1212 -- nested owner windows complete r=6

Let

```text
P={p1<...<p7} subset {1,...,12},       M=p7,
13M<k1<...<k6.
```

Suppose `rho` of the killers are divisible by `13`, and write their
normalized values as

```text
O={o1<...<o_rho},       x=o1,       carrier killers={13o:o in O}.
```

As in THM-1169, use the normalized carrier clock `u`.  A carrier-safe
boundary

```text
u=c/(14s),       s in O,                              (1)
```

produces the twelve owner charts

```text
t_a=a/13+u/13=a/13+c/(182s),       a=1,...,12.        (2)
```

If `u` lies in the least carrier's first safe window

```text
I_x=[1/(14x),13/(14x)],                               (3)
```

then the boundary readout of THM-1169 gives `cM<13s`.  For `rho>=4`, the
owner count is automatically `d<=7<2rho`, so THM-1154 absorbs all
`6-rho` noncarriers.  The only remaining task is therefore to prove that
the carrier combs do not cover (3).

No Covering hypothesis and no upper bound on a killer are used.

## 1. The nested first-safe window

For `x<=q<13x`, put

```text
J(x,q)=[1/(14x),13/(14q)],
L(x,q)=|J(x,q)|=(13x-q)/(14xq).                        (4)
```

> **Prefix-window lemma.**  Every carrier `o` with `x<=o<=q` is safe on
> the whole closed interval `J(x,q)`.

Indeed, throughout (4),

```text
1/14 <= o/(14x) <= ou <= 13o/(14q) <= 13/14.          (5)
```

There is no wrap and both equality endpoints are safe.  This elementary
window is the proof-bearing object below: each density contradiction adds
one carrier to the protected prefix and shortens only the right endpoint.

At radius `lambda=1/14`, THM-1155 says that if `m` danger combs of speeds
`W` cover a closed interval of length `L`, then

```text
(7-m)L <= sum_{w in W} 1/w.                            (6)
```

Because the carrier speeds are strictly ordered, if their least value is
`r`, the right side of (6) is strictly less than `m/r` whenever `m>=2`.

## 2. Five carriers close in one nested step

Let the five normalized carriers be

```text
x<y<z<w<v.                                             (7)
```

THM-1169 already closes the row when

```text
y>=14x/9.                                              (8)
```

Assume the complementary strict residual `y<14x/9`.  Since this is much
smaller than `13x`, the prefix-window lemma makes

```text
J(x,y)=[1/(14x),13/(14y)]                              (9)
```

safe for `x` and `y`.  If the later three combs covered it, (6) would give

```text
4L(x,y) <= 1/z+1/w+1/v < 3/y.                         (10)
```

But direct expansion gives

```text
4L(x,y)>3/y
  iff 4(13x-y)>42x
  iff 5x>2y.                                           (11)
```

The last inequality follows already from `y<14x/9<5x/2`, contradicting
(10).  Hence the five carriers have a common safe point in `I_x`; the
boundary readout (1) and THM-1154 finish the row.

> **Five-carrier theorem.**  Every r=6 row with exactly five
> `13`-divisible killers is `1/14`-lonely.

The proof uses the full metric interval (9).  A switch search on a bounded
integer range is unnecessary.

## 3. Six carriers: the exact nested ratio ladder

Write the six normalized carriers as (using `r` for the last carrier so the
carrier clock `u` from (1) remains unambiguous)

```text
x<y<z<w<v<r.                                           (12)
```

THM-1169 closes the row when `y>=35x/12`, so assume

```text
y<35x/12.                                              (13)
```

Suppose for contradiction that the five later combs cover `I_x`.  The
argument repeatedly restricts this cover to a prefix window.

### First handoff: protect x,y and bound z

The four later combs `z,w,v,r` cover `J(x,y)`.  Equations (4) and (6) give

```text
3L(x,y) <= 1/z+1/w+1/v+1/r < 4/z,
z < 4/[3L(x,y)] = 56xy/[3(13x-y)].                     (14)
```

The function `r -> 56r/[3(13-r)]` is increasing for `0<r<13`.  At the
strict endpoint in (13),

```text
z/x < 56(35/12)/[3(13-35/12)]
    =1960/363
    <27/5.                                             (15)
```

The last comparison is the exact one-unit inequality
`1960*5=9800<9801=27*363`.  In particular `z<13x`, so `J(x,z)` is a
nonempty closed interval safe for `x,y,z`.

### Second handoff: protect z and bound w

The remaining three combs cover `J(x,z)`.  Thus

```text
4L(x,z) < 3/w,
w < 3/[4L(x,z)] = 21xz/[2(13x-z)].                    (16)
```

Using the deliberately rounded strict bound `z/x<27/5` from (15),

```text
w/x < (21/2)(27/5)/(13-27/5)
    =567/76
    <13.                                               (17)
```

Therefore `J(x,w)` is safe for the first four carriers.

### Third handoff: protect w and bound v

The last two combs cover `J(x,w)`, so

```text
5L(x,w) < 2/v,
v < 2/[5L(x,w)] = 28xw/[5(13x-w)].                    (18)
```

Substitution of (17) yields

```text
v/x < 28(567/76)/[5(13-567/76)]
    =15876/2105
    <11.                                               (19)
```

Thus the first five carriers are all safe on `J(x,v)`.

### The last open tooth is too short

By (19), the length of this final protected interval satisfies

```text
L(x,v)=(13x-v)/(14xv)
      >2x/(14xv)
       =1/(7v)
       >1/(7r).                                        (20)
```

Every connected component of the last danger comb `D_r` is an **open**
tooth of length exactly `1/(7r)`.  The connected closed interval `J(x,v)`
is longer, so it cannot be contained in `D_r`.  This contradicts the
assumed carrier cover.

> **Six-carrier theorem.**  Every r=6 row whose six killers are all
> divisible by `13` is `1/14`-lonely.

The constants in (15)--(19) are not numerical reconnaissance.  They are
the exact monotone handoffs forced by the formal density inequality (6).
Endpoint strictness always points in the favorable direction: danger teeth
are open, prefix windows are closed, and distinct carrier order makes each
reciprocal estimate strict.

## 4. The entire r=6 stratum is closed

The carrier cardinalities now form a complete dispatch.

```text
rho=0:  THM-1136's five-noncarrier shadow;
rho=1:  THM-1136's unique-carrier corollary;
rho=2:  THM-1154's owner chart;
rho=3,4: THM-1169's owner windows;
rho=5,6: this theorem.
```

> **Uniform r=6 closure.**  For every seven-element
> `P subset {1,...,12}` and every six distinct killers
> `13 max(P)<k1<...<k6`, the thirteen-speed family
> `P union {k1,...,k6}` has a `1/14`-lonely time.

This removes the whole r=6 clustered branch from the LRC(14) proof map.
THM-1214 subsequently closes the adjacent clustered r=5 stratum by the
five-killer owner ledger; neither statement is a global LRC(14) proof.

This also supersedes the earlier all-shapes tail/enumeration residual recorded
in the THM-1121/1132/1134 horn programme.  There is no contradiction: that
programme treated all six danger combs inside the raw core-safe complex and
therefore faced an enormous shape bank.  The THM-1136/1154 quotient used here
first absorbs every noncarrier and keeps only normalized `13`-carrier phases;
the labelled prefix windows then solve that smaller compatibility problem
without enumerating raw killer shapes.

## 5. Tournament Analysis and assumption challenge

The speed-order tournament on five or six carriers is transitive and cannot
see any step of the proof.  The exact proof carrier is instead the chain

```text
J(x,y) -> J(x,z) -> J(x,w) -> J(x,v),                  (21)
```

with each vertex labelled by its closed endpoints, protected carrier
prefix, remaining-comb count, and reciprocal density budget.

For telemetry, orient two window obligations toward the one with the larger
protected prefix (equivalently, later in the containment proof order).
The four-vertex tournament is transitive, with score sequence `(3,2,1,0)`,
no directed cycle, four singleton SCCs, one Hamiltonian path, and six flips
under reversal of the containment gauge.  The pairwise observable is prefix
containment; the switch is the protected-prefix direction; the tie
Hamiltonian path is `J2 -> J3 -> J4 -> J5` (there are no metric ties).

This tournament is only an index.  It destroys the window lengths, endpoint
owner, strict open-cover predicate, and density slack that prove (20).
Candidate vertices explicitly challenged here are runners, carriers,
ratios, danger teeth, wall events, nested windows, and cover obligations.
The smallest faithful object is the labelled nested-window/cover ledger.

## 6. Exact replay and formal boundary

The dependency-free referee checks all rational identities and monotone
handoffs with `Fraction`, replays every integer first-pair residual through
large boxes, and independently constructs the exact danger-wall cells for
compact and adversarial 14-band carrier families.  Ordinary and optimized
runs are byte-identical, and every invariant uses a `require` that remains
live under `python -O`; the script contains no Python `assert` statements.

```text
04-computation/lrc14_r6_five_six_carrier_nested_window_referee_codex_S77.py
SHA-256 47aaa4e886571af2ae3e77c20712ba4596928ce3e068ad67b560bcfd4661b4ee

05-knowledge/results/lrc14_r6_five_six_carrier_nested_window_referee_codex_S77.out
SHA-256 93308505434db3aad56afecc0905a16dbd3c6704c8ecbb73e272f63d7d79b9ac
```

The only analytic input beyond elementary order arithmetic is THM-1155's
`multi_speed_density_bound`, already kernel-pure in
`LRCEssentialRegion.lean`.  Formalization now instantiates it on the actual
rho-five window and on all three actual rho-six tail Finsets.  It also proves
that a closed interval longer than one open tooth is not contained in that
comb.  The remaining geometric composition is narrower: derive each
restricted tail cover from a cover of the original window by excluding the
already protected carriers, extract the least common-safe point as a boundary
`c/(14s)`, and feed its broad-cap data to a future formal THM-1154 owner-chart
interface.  No finite certificate or unformalized asymptotic claim lies
inside the mathematical proof.

`LRCNestedCarrierWindow.lean` kernel-checks the dimensionless five-carrier
contradiction, all three six-carrier handoffs, the exact rounding constants,
the final tooth comparison, their composed ratio ladder, the concrete
rho-five non-cover theorem, all three concrete rho-six cover-to-ratio
handoffs, and final closed-interval/open-tooth topology.  Every public audit
uses only `[propext, Classical.choice, Quot.sound]`.
