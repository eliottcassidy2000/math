# Nine-subset gcd ceiling and the surviving clock-nine boundary

**Status: PROVED RELATIVE TO CITED LRC THROUGH THIRTEEN RUNNERS +
FINITE-EXACT profile classification + INDEPENDENTLY AUDITED.** Let `V` be thirteen distinct positive
integer speeds, `gcd(V)=1`, and suppose `M(V)<1/14`. For every nine-subset
`P` of `V`,

```text
gcd(P) belongs to {1,2,3,4,6,8,9}.                    (1)
```

In particular, any primitive thirteen-speed row containing nine speeds
whose gcd exceeds nine is safe. Without primitivity, the equivalent
sufficient condition is `gcd(P)>9 gcd(V)`. This is a reduction on a
hypothetical strict counterexample, not an unsafe example or LRC(14).

## 1. Recovered mechanisms and status

The inherited [effective-clock result](lrc14_effective_clock_empty_core_sep06.md)
proved twelve-, eleven-, and ten-subset gcd bounds `1,2,4` respectively.
Incoming [THM-4447, composite-clock capacity](../../01-canon/theorems/THM-4447-lrc14-composite-clock-capacity-and-small-clock-reduction.md)
now independently canonizes that ten-subset bound and its exact profiles.
It also supplies the divisor-absorption operation used below. Its literal
pointwise count retains orbit multiplicity and the strict endpoint.

[THM-4446, primitive ten-pack descent](../../01-canon/theorems/THM-4446-lrc14-primitive-ten-pack-descent-and-dilation-rays.md)
is now proved: its scale-three nonprimitive-body conclusion is covered by
the global ten-subset cap, while its phase-uniform fibre claim, dilation-ray
consumer and typed entry restrictions carry additional information.
[THM-4448, general shore attachment](../../01-canon/theorems/THM-4448-lrc14-general-shore-attachment-and-decoder-pair-cones.md)
retains one protected inverse sheet and proves attachment cones; it does
not give an unconditional nine-subset gcd bound. All three current canon
files were read at clean inherited `HEAD=6dd59c9c41`.

Targeted prior searches recovered
[THM-4079's nine-speed antipodal obstruction](../../01-canon/theorems/THM-4079-lrc14-antipodal-outlier-absorption-and-adaptive-clock.md),
which concerns a different phase predicate, not a nine-pack gcd ceiling.
No exact statement `(1)` was recovered in that bounded search. This note
records an assembly of known deletion and capacity mechanisms and makes
no external priority claim. Root independently recovered the same upper
bound and twenty initial profiles during this heartbeat.

The corrected near miss is to discard the effective order
`c/gcd(c,w)` or to turn a budget of one into a free-sheet certificate.
The least-used sidecar is the gcd after adjoining two or three tails.
The live concepts are deletion inheritance, labelled orbit capacity,
divisor absorption, actual phase partition, and safe realization.

## 2. All-height ceiling

Write `c=gcd(P)>1`, and let `w_1,...,w_4` be the four complementary speeds.
Set

```text
g_i=gcd(c,w_i),       q_i=c/g_i,
beta(q)=ceil(q/7)/q.
```

The inherited bounds on larger subsets imply

```text
g_i<=4,
gcd(g_i,g_j)<=2,
gcd(g_i,g_j,g_k)=1                  for distinct indices. (2)
```

Equivalently, `q_i>=c/4`, `c/lcm(q_i,q_j)<=2`, and every triple of
orders has lcm exactly `c`. At a safe phase for the divided nine-body,
its `c` physical lifts are all body-safe. Tail `w_i` has `q_i` distinct
orbit points, each with multiplicity `g_i`, so at most
`c beta(q_i)` labels are bad. The cited lower-runner theorem guarantees
such a body phase, with margin at least `1/10`. Therefore failure forces

```text
sum_(i=1)^4 beta(q_i)>=1.                              (3)
```

For every integer `q>=4`, `beta(q)<=1/4`, with equality only at
`q=4` or `8`. Check `4<=q<=8`; for `q>=9`, use
`beta(q)<=(q+6)/(7q)<1/4`.
If `c>=13`, all `q_i>=4`, so `(3)` forces all four orders to belong
to `{4,8}`. Every triple then has lcm at most eight, contradicting `(2)`.
Thus `c<=12`.

The three remaining clocks above nine are excluded analytically:

- `c=10`: `g_i` is one or two, so `q_i` is ten or five and
  every `beta(q_i)=1/5`; total budget is `4/5`.
- `c=11`: every `g_i=1`; total budget is `8/11`.
- `c=12`: possible orders are `3,4,6,12`. The pair-gcd bound allows at
  most one order three and at most one order four. Thus total budget is
  at most `1/3+1/4+1/6+1/6=11/12`.

All are strictly below one. At `c=5` or `7`, the bound `g_i<=4` forces
all four tails coprime to `c`, and the respective budgets `4/5` and
`4/7` again fail. This proves `(1)` without a speed-height census.

## 3. Complete finite order profiles and two absorption exits

The complete sorted order quadruples satisfying only `(2)-(3)` are:

| Exact body gcd `c` | Sorted effective orders |
|---|---|
| 2 | `(1,1,2,2)`, `(1,2,2,2)`, `(2,2,2,2)` |
| 3 | `(1,3,3,3)`, `(3,3,3,3)` |
| 4 | `(1,2,4,4)`, `(1,4,4,4)`, `(2,2,4,4)`, `(2,4,4,4)`, `(4,4,4,4)` |
| 6 | `(2,3,3,6)`, `(2,3,6,6)`, `(2,6,6,6)`, `(3,3,6,6)` |
| 8 | `(2,4,8,8)`, `(2,8,8,8)`, `(4,4,8,8)`, `(4,8,8,8)`, `(8,8,8,8)` |
| 9 | `(3,9,9,9)` |

These are **twenty abstract scalar/deletion profiles**. They are not all
necessary after the stronger exact ten-pack profile is inherited.
THM-4447's divisor absorption removes two of them:

```text
c=4, q=(1,4,4,4): absorb the order-one tail at clock four;
c=8, q=(2,8,8,8): absorb the order-two tail using divisor four.
```

In either case the absorbed body has ten speeds and the three residual
tails are odd. At clock four they kill at most three labels in total,
leaving a safe label. Equivalently, the resulting ten-subset has gcd
four but complementary orders `(4,4,4)`, which the exact ten-pack
classification excludes.

Exactly **eighteen** profiles survive after this absorption. For each of
the twenty profiles, the sidecar checks every divisor `d>1` of `c` and
computes

```text
Cap_d=sum_(i:d does not divide g_i)
          gcd(d,g_i) ceil(d/[7 gcd(d,g_i)]).              (4)
```

Because `d|c`, one has `gcd(d,w_i)=gcd(d,g_i)`, so the gcd word determines
this capacity exactly. Primitivity prevents all four tails being absorbed,
and the enlarged pack has at most twelve speeds. The cited lower-runner
input therefore applies. Only the two displayed profiles have any
`Cap_d<d`. The independent inherited-ten-profile check gives the same two
exits. Neither test proves that the eighteen surviving profiles contain
an unsafe row.

## 4. Safe realizations and a genuine fully-spoiled phase at clock nine

Every one of the twenty initial gcd profiles can occur in a **provably safe**
primitive thirteen-speed row. Given its `c` and `q_i`, put `g_i=c/q_i`,
index `i=0,...,3`, and set

```text
P=c*{1,...,9},
k_i=1+14(i+1),       w_i=g_i(1+c k_i),
x=1/(14c).                                               (5)
```

All tails exceed `9c`; their gcd with `c` is exactly `g_i`. If two tails
were equal, reducing modulo `c` first gives equality of their `g_i`
(the integers lie in `[1,c]`), and then equality of their distinct `k_i`,
a contradiction. The full row is primitive by the triple-gcd condition.
At `(5)` the body phases are `1/14,...,9/14`, while tail fractional phases
are

```text
(g_i+1/q_i)/14 in [1/14,5/14].                            (6)
```

Thus the full row is safe. This construction realizes the exact gcd
coordinates, not any supposed counterexample property.

The maximal residual clock also has an exact hostile to a phase-uniform
selector. Take

```text
c=9,       C={1,3,5,...,17},       T=(1,5,6,7),       y=1/2.
```

The divided nine-body is strictly safe at `y`. On the nine lifts
`x_j=(y+j)/9`, the tails' strict bad-label sets, in the displayed order, are

```text
{0,8},       {3,5},       {1,4,7},       {2,6}.             (7)
```

These sets partition all nine labels and attain the capacity of
`q=(3,9,9,9)`. All chosen lifts are spoiled. Nevertheless the full row
`9C union T` is safe at `x=1/14`: every body speed is odd, and the four
listed tails also have clearance at least `1/14` there. The failed
implication is therefore “every body-safe phase has a free lift,” not
lonely-time existence.

At divisor three, the six-tail is absorbed into a ten-body and leaves
three ternary-unit tails: the known clock-three component-address problem.
This is the smallest concrete next consumer. It must retain where the
body-safe components lie; capacities and hereditary gcd data alone do not
select a successful component.

## 5. Reproduction and stopping boundary

The [standalone source](../../04-computation/gcd_nine_audit_empty_core_next_sep06.py)
and [output](gcd_nine_audit_empty_core_next_sep06.out) use standard-library
exact fractions. Two separate complete classifiers enumerate gcd words and
effective-order words for `2<=c<=12`, with all pair/triple filters. They
agree on every row and the literal twenty-profile bank. The script also
checks every divisor, the exact inherited ten-pack signatures, all twenty
safe realizations, and `(7)` with its separate full-row witness.

```sh
python3 -B 04-computation/gcd_nine_audit_empty_core_next_sep06.py
python3 -B -O 04-computation/gcd_nine_audit_empty_core_next_sep06.py
```

There are 192 explicit gates. Normal and optimized outputs are byte-identical.
Semantic digest:
`1be9a4748b3e48af0bf32daf8d8f3d7f0765ddda3839e2ba24bbe73d9d253a32`.
Raw SHA-256 source/output:

```text
901ce3ec69ec771b56b5f8b5f79251aac171d3f7459215ab0758b9fdf43c0ebe
91ccda1b2fc6e4cbdc88faf8ca4172e22f7c51eb8d672d6a5b33d0ee176d6ed6
```

The connection sends a nine-body plus four tails to exact divisor orders,
preserving every inherited subset gcd and the uniform scalar branch budget.
It discards tail phase placement. The safe realization and fully-spoiled
phase above exhibit that loss explicitly. The global ceiling is proved;
emptiness of the remaining phase-labelled profiles remains **OPEN**.

Root final audit: **PASS** on the all-height proof, all twenty/eighteen profile distinctions, exact safe realizations and the clock-nine phase partition. The independent recursive hierarchy also reproduces all eighteen retained nontrivial profiles.
