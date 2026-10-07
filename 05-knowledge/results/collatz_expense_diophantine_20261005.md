> **CORRECTED 2026-10-05.** The exact type/rate formula and finite observations survive. The full rate tail jumps at every realized rate, not only record rates: `(l,A)=(4,7)` contributes `3/128` at the nonrecord rate 7. Stretched-exponential tail decay, logarithmic bank returns, and the assertion that no other adaptive deadline is available are **UNPROVED / RETRACTED AS CONCLUSIONS**. The proof below restores the missing density-tail argument and keeps coefficient stopping distinct from literal descent. Rate ceilings are now rationally certified; the DP's omitted mass is retained. This is a targeted correction, not a full audit of every inherited result.

# Diophantine expense: exact type rates, density tails, and scoped bank bounds

Original note: opus session `opus-2026-10-05-S9`; corrected by the concurrent
Codex audit on 2026-10-05. **PROVED:** Q's normalized rational-record identity,
T's natural-density series with the uniform-tail proof below, and B's stated
family-relative telescoping bound, conditional on its finite core being
certified. **FINITE-EXACT:** integer/Fraction DP and rate tests. The original
finite census tables are retained; the repaired run checks every observed
rate exactly and rejects int64 overflow. Decimal displays and floating
checks of the deadline bound remain numerical. **OPEN:** universal entry,
the full tail asymptotic, and asymptotic bank efficiency.

[Script](../../04-computation/experiments/collatz_expense_diophantine_20261005.py),
[output](collatz_expense_diophantine_20261005.out), and
[JSON](collatz_expense_diophantine_20261005.json).

## 1. Inheritance and definitions

The closest proved mechanism is [weak-reset families](collatz_weak_reset_family_20261005.md),
Proposition R and Lemma D1/Theorem D: the word carry supplies the actual
source threshold, and a uniform admissible rate pays the stated logarithmic
deadline on its domain. The source-identity and carry boundary is also
explicit in [THM-4512, coefficient-descent classes](../../01-canon/theorems/THM-4512-coefficient-descent-classes-one-member.md).
The former measurement-independence note's identification of coefficient
survival with literal non-descent is not used.

Write `alpha=log_2 3` and `U(n)=oddpart(3n+1)`. A valuation word of length
`l` and total valuation `A` has affine map `(3^l n+B_w)/2^A`, with
`B_w>0`. It is coefficient-descending when `2^A>3^l`; actual descent needs
`n>B_w/(2^A-3^l)`. Proper coefficient-rising prefixes certainly rise on
positive sources, but the converse is not asserted pointwise.

A first **literal** descent segment ends at the first odd iterate below
its positive start `n>1`. Its rate is the least positive integer `q` with
`2^(qA-l)>=3^(ql)`. The running-minimum chain concatenates these segments.
ROOT is handled separately. A chain rate or a finite census presupposes
that the relevant segments exist; a finite census proves only its range.

The live board is **type / normalized upper approximation / carry threshold /
cylinder mass / uniform remainder / bank-relative domain**. The canonical
hostile to a record-only staircase is type `(4,7)`. The missing sidecar in
the former density proof was a uniform tail bound, not merely finite
exceptions in each separate cylinder.

## 2. Theorem Q: rates and normalized approximation records

For every descending type,

```text
q(l,A)=ceil(l/(A-l alpha)),       A_min(l)=ceil(l alpha).
```

This follows by taking logarithms of the defining integer inequality.
Increasing `A` decreases the rate. The minimal type is realized by the
word `1^(l-1),A_min(l)-l+1`, so it really supplies the maximum at length `l`.

The strict record lengths of `q(l,A_min(l))` are exactly the strict records
of the rational upper bounds `A_min(l)/l`. To prove the ceiling step, put
`e_l=A_min(l)-l alpha` in `(0,1)` and `t_l=l/e_l`. A nonrecord rational
bound cannot give a new ceiling record. At consecutive strict rational
records `j<l`,

```text
t_l-t_j=(A_min(j)*l-A_min(l)*j)/(e_j e_l)>1.
```

The numerator is a positive integer. Thus no new rational record is lost
under ceiling. This elementary proof uses normalized error; the former
proof substituted the unnormalized error `e_l` without justifying the
equivalence. That equivalence is also true, with the following additional
argument. An `e_l` record is immediately a record of `delta_l=e_l/l`.
Conversely, suppose `delta_l` is a strict record but some `r<l` has
`e_r<=e_l`. Irrationality makes this inequality strict. Since
`delta_r>delta_l`, the positive upper approximation with denominator
`s=l-r` and numerator `A_min(l)-A_min(r)` has normalized error
`(e_l-e_r)/(l-r)<delta_l`. Replacing its numerator by `A_min(s)` can only
improve it, contradicting the record at `l`. Thus both record notions
coincide. This subtraction argument was recovered independently during the
correction audit; it is not a continued-fraction hypothesis. The repaired
comparison function still uses exact `Fraction(A_min(l),l)` values. The
usual continued-fraction interpretation is background, not a dependency here.

| `l` | `A_min` | `q_max` | `2^A - 3^l` | fraction |
|---:|---:|---:|---:|---:|
| 1 | 2 | 3 | 1 | 2/1 |
| 3 | 5 | 13 | 5 | 5/3 |
| 5 | 8 | 67 | 13 | 8/5 |
| 17 | 27 | 306 | 5,077,565 | 27/17 |
| 29 | 46 | 804 | `1.7e12` | 46/29 |
| 41 | 65 | 2,480 | `4.2e17` | 65/41 |
| 94 | 149 | 6,951 | `6.7e42` | 149/94 |
| 147 | 233 | 13,984 | `1.0e68` | 233/147 |
| 200 | 317 | 26,668 | `1.4e93` | 317/200 |
| 253 | 401 | 56,382 | `1.6e118` | 401/253 |
| 306 | 485 | 207,489 | `1.0e143` | 485/306 |

The former 60-digit ceiling gate has been replaced by exact rational log
intervals. For `0<t<1`, truncate
`2 sum_(j>=0) t^(2j+1)/(2j+1)` after `M` terms; its positive remainder is at
most `2t^(2M+1)/((2M+1)(1-t^2))`. Use `t=1/3,1/2` for `log2,log3`, divide
positive enclosing intervals, and increase `M` until both possible rate
endpoints have the same ceiling. Irrationality guarantees termination.
Small cases also receive an independent integer-power test. The displayed
large gaps are rounded scientific notation, not exact integers.

## 3. Theorem T: the full natural-density series

Let `N(l,A)` count words with positive valuations, coefficient-rising
proper prefixes `2^(A_j)<3^j` for `j<l`, and `2^A>3^l` at the last step.
Then the natural density among odd positive sources with first literal
descent type `(l,A)` is `N(l,A)/2^A`. Moreover the density with first-segment
rate greater than an integer `Q` is

```text
sum_(l,A with q(l,A)>Q) N(l,A)/2^A.                 (T)
```

Sources with no first literal descent form a density-zero set by the proof
below; assigning them infinite rate does not change (T). This is not a
pointwise assertion that literal and coefficient stopping agree.

**Fixed-type proof.** Each exact word of total valuation `A` has one odd
cylinder modulo `2^(A+1)`, of relative odd density `2^-A`. For a word whose
proper coefficient prefixes rise, deleting the finitely many positive
sources below its final carry threshold leaves exactly literal first
descent there. A fixed type has finitely many words. A word of that type
with an earlier coefficient exit can be a first literal descent only when
the source lies below that earlier exit's finite carry threshold. It adds
only finitely many sources. Thus the displayed fixed-type density follows.

**Uniform tail, needed for the countable sum.** Under normalized odd Haar
measure the valuations have product masses `P(a_i=a)=2^-a`. This follows
from the exact cylinders, not an assumption on a supplied integer. Put
`S_r=a_1+...+a_r`. Coefficient survival through `r` implies
`S_r<r alpha<8r/5`, since `3^5<2^8`. With `t=3/4`,

```text
E[t^a]=t/(2-t)=3/5,
P(coefficient survival through r)
 <= [(3/5)(4/3)^(8/5)]^r = theta^r,
theta^5=65536/84375<1.                              (1)
```

Indeed on that event `t^S_r>t^(8r/5)`, so the nonnegative moment bound gives
(1). The finite survival event is a finite union of valuation cylinders:
its partial sums are bounded by `r alpha`. Its natural density therefore
equals its Haar mass.

For fixed `r`, only finitely many sources can have a coefficient exit by
step `r` but fail literal descent at that exit. There are finitely many
proper coefficient-rising prefixes. For each such prefix the terminal
carry `B_w` is independent of the terminal exponent; when that exponent
is sufficiently large, `B_w/(2^A-3^l)<1`, leaving no positive exception.
Only finitely many remaining terminal exponents and bounded sources remain.

Consequently the set where first coefficient and first literal descent
differ has upper natural density at most `theta^r` for every `r`, hence
zero. The same applies to sources with no coefficient exit, and therefore
to sources with no literal descent. Finally the disjoint coefficient-exit
cylinders have total mass one: the remaining survival mass tends to zero.
For any finite selection of these cylinders, the complement has its exact
complementary density; enlarging finite selections until their total mass
tends to one proves countable additivity for every selected subseries.
Applying this to `q(l,A)>Q` proves (T). No unsupported exchange of natural
density with an arbitrary countable union is used.

**Finite DP and its bill.** The DP retains `l<=200` and
`A<=A_min(l)+40` (8,200 types). If its exact total mass is `M`, then every
reported partial tail lies below the full tail by an amount in `[0,1-M]`.
Here `1-M` is about `5.5e-8`. Threshold membership is now checked with the
exact rate routine. The following rounded entries are the retained
original observations, not an exact equality between the two rows:

| `Q` | 2 | 3 | 7 | 8 | 13 | 31 | 67 | 104 | 310 | 800 | 2000 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| truncated tail lower sum | 0.63006 | 0.25903 | 0.18743 | 0.17631 | 0.08182 | 0.05549 | 0.00696 | 0.00646 | 0.00093 | 0.00078 | 0.00018 |
| census `2^20` | 0.63009 | 0.25900 | 0.18711 | 0.17598 | 0.08138 | 0.05526 | 0.00686 | 0.00643 | 0.00088 | 0.00072 | 0.00015 |

The `(5,8)` type has seven words and mass `7/256`. Other record-type masses
include `N(17,27)=312455`, and the displayed approximate densities
`5.4e-4,1.7e-4,2.6e-6` at lengths29,41,94. These do not exhaust the jumps.
For example the three words `1114,1123,1213` have type `(4,7)` and exact
rate7; they cause a jump of at least `3/128` between thresholds6 and7,
although the neighboring worst-type record rates are3 and13.

No full asymptotic `exp(-c sqrt(Q))` is proved here. The heuristic needs
count asymptotics and comparability of successive approximation
denominators; the finite record list establishes neither. In particular
`q` being comparable to a product of successive denominators would not
make either denominator uniformly comparable to `sqrt(q)`.

## 4. Inheritance and the scoped bank bound

Every running minimum of a completed source starts one of its segments,
so its chain rate is inherited. The original finite census below `2^20`
found289666 sources with chain rate above30, against29072 with an own
segment rate above30. The following are finite counts:

| start | type | `q` | inheritors | share |
|---:|---:|---:|---:|---:|
| 47 | (34,55) | 31 | 49,927 | 0.095 |
| 31 | (35,56) | 67 | 48,898 | 0.093 |
| 91 | (28,45) | 46 | 45,951 | 0.088 |
| 95 | (5,8) | 67 | 27,694 | 0.053 |
| 379 | (5,8) | 67 | 11,384 | 0.022 |
| 103 | (26,42) | 33 | 11,201 | 0.021 |
| 847 | (5,8) | 67 | 6,664 | 0.013 |
| 455 | (5,8) | 67 | 5,349 | 0.010 |
| 71 | (32,51) | 114 | 4,479 | 0.009 |
| 335, 1243, 1711 | (5,8) | 67 | 4,025 / 3,766 / 3,321 | 0.008 / 0.007 / 0.006 |

**Proposition B.** Fix a threshold `Y` and rate `q_hi`, and let
`Y'=max(Y,10q_hi)`. Assume the finite odd core below `Y'` is certified.
For sources whose running-minimum segments starting at least `Y` are all
`q_hi`-admissible, the inherited segment inequality gives

```text
tau(n)<=floor(1.051 q_hi log_2 n)+max_{odd y<Y'} tau(y).    (2)
```

Proof: while segment starts remain at least `Y'`, Lemma D1 bounds total
length by the telescoping full logarithmic drop from `log_2 n`. The last
segment may land far below `Y'`, so this argument does not replace that
quantity by `log_2(n/Y')`. On entering the certified core, append its
receipt. If a finite certified bank supplies an earlier stop, append that
receipt instead; retain the small-core alternative and its bound. Finite
core certification is a premise for arbitrary parameters, not a conclusion
from the below-`2^20` census.

Original coverage among sources between `Y` and `2^20`:

| `Y` | `q_hi = 3` | 8 | 16 | 32 |
|---:|---:|---:|---:|---:|
| `2^6` | 0.112 | 0.216 | 0.474 | 0.571 |
| `2^10` | 0.204 | 0.346 | 0.605 | 0.700 |
| `2^14` | 0.373 | 0.515 | 0.731 | 0.804 |
| `2^16` | 0.494 | 0.623 | 0.801 | 0.857 |

Original bank-relative finite proportions, all sources below `2^20`:

| bank size | rate 3 | rate 8 |
|---:|---:|---:|
| 0 | 0.061 | 0.149 |
| 10 | 0.096 | 0.227 |
| 100 | 0.123 | 0.272 |
| 1,000 | 0.166 | 0.331 |
| 10,000 | 0.233 | 0.406 |

These tables do not prove logarithmic returns in bank size or impossibility
of shorter deadlines from other information. The claim that the bank
threshold is the only adaptive parameter is withdrawn. A source-aware
analytic bound or a different family grammar is not ruled out by (2).

A narrower positive-density obstruction does follow for the **named fixed
rate/finite-bank grammar**. Fix `q` and choose any `l>q`. The word
`1^(l-1),A_min(l)-l+1` has rate greater than `q`, because its excess is less
than one. Its proper prefixes rise. Outside its finite carry exception and
any fixed finite bank/core, its cylinder of density `2^-A_min(l)` has an
unpaid hard first segment. Thus that grammar cannot have density-one
coverage. No independence or renewal law is needed, and this does not
apply to unbounded adaptive rates or other methods.

## 5. Reproduction and correction boundary

```text
python 04-computation/experiments/collatz_expense_diophantine_20261005.py --census-bits 20 --json 05-knowledge/results/collatz_expense_diophantine_20261005.json
python -O 04-computation/experiments/collatz_expense_diophantine_20261005.py --census-bits 20
```

The corrected program uses exact normalized rational record comparisons,
rationally certified ceilings, and exact DP masses with an explicit missing
mass. The numba census now rejects before `3v+1` can overflow signed int64,
and every distinct observed type/rate is independently rechecked exactly.
Its finite successful run certifies only that range. Decimal proportions,
scientific-notation gaps, timing, and floating deadline comparisons are not
integer identities. Universal natural-density claims rest on the proof in
section3, not on the census. Universal convergence remains open.
