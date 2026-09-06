# Independent literal audit of the scale-three Haar consumers

**Status: INDEPENDENT AUDIT PASS.** The adjacent-family conclusion is
PROVED for all its stated heights. The finite-head statements are
FINITE-EXACT. The parity-free ceiling and proposed arbitrary ten-body floor
are REFUTED by literal controls; these are not counterexamples to LRC(14).

Reviewed the producer
[source](../../04-computation/overnight3_20260906_lrc_consumers.py) and
[report](overnight3_20260906_lrc_consumers.md). The audit
[source](../../04-computation/overnight3_20260906_lrc_consumers_audit.py)
imports no producer function or other repository module. Its
[preserved output](overnight3_20260906_lrc_consumers_audit.out) gives every
claimed finite count and the exact polynomial inequality certificates.

## 1. Independent physical predicates

For each tail `T`, the audit partitions the circle at all exact rational
points where

```text
||w(x+j/3)|| = 1/14,     w in T, j=0,1,2.
```

On each resulting open cell it tests the actual predicate

```text
for every j in {0,1,2}, some w in T has ||w(x+j/3)||<1/14.
```

It does not enumerate assignments of tails to sheets, construct their six
interval intersections, or use a previously derived tail-mass formula.
The common grid denominator is `42*lcm(T)`; midpoint comparisons are exact
integer comparisons on twice that denominator. Period `1/3` is used only
because translation by `1/3` permutes the three literal tests. The resulting
cell length is multiplied by three to obtain physical circle measure.

For bodies, the audit separately partitions at all points
`||cx||=1/14` and tests the literal conjunction `||cx||>=1/14`.
Midpoint tests give the positive-length pieces. Every breakpoint is also
tested independently, so isolated weakly safe points remain in the answer.
Closed safe intervals must have safe endpoints; this is an explicit gate.

## 2. Complete finite tail and body replays

The literal tail audit reproduces:

| Tail | Physical failure-set mass |
|---|---:|
| `(1,5,7)` | `8/245` |
| `(1,5,11)` | `6/77` |
| `(1,2,4)` | `0` |
| `(1,2,5)` | `0` |
| `(1,4,5)` | `1/28` |
| `(2,5,7)` | `22/245` |
| `(1,10,11)` | `6/55` |

I replayed **all 8,664 primitive triples** of distinct positive integers
at most sixty, none divisible by three. There are exactly 136 violations
of `6/77`; the smallest by maximum speed is `(2,5,7)`, and the unique
head maximum is `6/55` at `(1,10,11)`. No parity or additive-relation filter
was imposed. Of the violations, exactly one is nonadditive:

```text
(2,11,20),       mass=11/140.
```

Thus deleting only triples satisfying `a+b=c` does not repair the bound.
The census does not prove an arbitrary-height mixed-parity upper bound.

For `C=(1,2,3,5,7,8,9,11,12,13)`, the independently recovered body is:

```text
measure = 14249/252252;
12 positive-length closed intervals;
isolated points = {1,3,5,9,11,13}/14.
```

The six lower-half intervals match the producer's listed rational endpoints
exactly, and the remaining intervals are their reflections. The measure and
extremizer are recovered from
[THM-530 / lrc-gp-intersection-global-witness-floor, Section A](../../01-canon/theorems/THM-530-lrc-gp-intersection-global-witness-floor.md).
That section's universe is a subset of `{1,...,13}`. It supplies no floor
for arbitrary ten-speed bodies, and the producer properly preserves this
scope and provenance.

All 66 ten-subsets of `{1,...,12}` were checked: their minimum is
`16277/194040>6/77`. All 286 ten-subsets of `{1,...,13}` were checked:
exactly twelve have measure below `6/77`. No primitivity filter is needed.
This proves the claimed first maximum-speed boundary of thirteen.

For the full thirteen-speed set

```text
3C union {1,5,11} = (1,3,5,6,9,11,15,21,24,27,33,36,39),
```

the literal safe measure is `25331/756756`. There are eighteen
positive-length closed intervals and twelve isolated points:

```text
1/14, 5/42, 3/14, 11/42, 13/42, 5/14,
9/14, 29/42, 31/42, 11/14, 37/42, 13/14.
```

In particular `x=1/14` passes every weak-safety predicate. Since the body
measure is below the tail failure-set mass `6/77`, this independently
confirms that the scalar Haar gate is sufficient but not necessary. The
fixed-denominator witness is not a new LRC family closure.

## 3. Independent derivation of the all-height adjacent formula

Let `T=(1,b,b+1)`, where `b>=4` and `b=1 mod3`, and put `delta=1/14`.
The failure set is invariant under `x -> x+1/3` and reflection. Speed one
can fail only in the three disjoint `delta`-neighborhoods of those translates.
It therefore suffices to count `0<x<delta` and multiply by six.

On that interval speed one occupies sheet zero. Since each ternary-unit
tail can occupy at most one sheet, the other two tails must occupy sheets
one and two. The residues `b=1` and `b+1=2 mod3` show that `bx` and
`(b+1)x` must both lie within `delta` of a common number `t=k/3`, where
`k>=1` and `3` does not divide `k`. Their integer translates must agree:
their difference is `x<delta`, whereas different translates would require
a separation at least `1-2delta`.

For such a center the exact interval is

```text
((t-delta)/b, min(delta,(t+delta)/(b+1))).
```

The lower endpoint is positive. Distinct center intervals are disjoint
because center spacing is at least `1/3>2delta`. The interval is nonempty
exactly for `t<(b+1)/14`. Its length is

```text
((2b+1)/14-t)/(b(b+1))           if t<=b/14,
((b+1)/14-t)/b                  otherwise.
```

This proves the producer's physical formula, with no appeal to the
odd-tail theorem and with strict endpoints affecting only measure-zero
sets.

For a standalone exact summation, define

```text
N(K)=K-floor(K/3),
S(K)=K(K+1)/2-3 floor(K/3)(floor(K/3)+1)/2,
K0=floor(3b/14),        K=floor((3b+2)/14).
```

Here `N` counts positive integers not divisible by three and `S` sums
them. The complete numerator for the mass, over denominator `7b(b+1)`, is

```text
P(b)=3(2b+1)N(K0)-14S(K0)
     +(b+1)*(3(b+1)(N(K)-N(K0))-14(S(K)-S(K0))).      (1)
```

Writing `b=42q+r`, `r=1,4,...,40`, both cutoffs are `9q+constant`.
Their counts and sums are therefore polynomials of degree at most two,
and (1) is exactly

```text
P(b)=1134q^2+B_r q+C_r.
```

Independent expansion gives:

| r | B_r | C_r |
|---:|---:|---:|
| 1 | 54 | 0 |
| 4 | 204 | 5 |
| 7 | 396 | 31 |
| 10 | 630 | 84 |
| 13 | 738 | 120 |
| 16 | 846 | 156 |
| 19 | 1080 | 253 |
| 22 | 1188 | 307 |
| 25 | 1422 | 444 |
| 28 | 1530 | 516 |
| 31 | 1638 | 588 |
| 34 | 1872 | 769 |
| 37 | 2064 | 935 |
| 40 | 2214 | 1080 |

To prove the lower bound, check `11P(b)-6b(b+1)>0`. After the required
shift `q=u+1` for residues one and four, all three coefficients are strictly
positive; the smallest constant is five, in residue seven. For the upper
bound, check `42b(b+1)-55P(b)>=0`, shifting residue one to `q=u+1` for
`b>=4`. All coefficients are nonnegative; the leading and linear ones are
positive. The only zero constant occurs at residue ten, where the polynomial
is `11718q^2+2394q`. Hence

```text
M(4)=1/28;
6/77 < M(b) <= 6/55 for every admissible b>=7;
equality with 6/55 occurs only at b=10;
lim M(b)=1134/(7*42^2)=9/98.
```

The audit output preserves every shifted polynomial's coefficients. It
also compares the residue formula with the independent literal predicate
in every admissible residue at `q=0,1,2`. These numerical comparisons are
controls; the count/sum identities and coefficient signs prove the
unbounded statements. The broader additive-triple error bound is handled
by the separate additive audit and is not a dependency of this proof.

## 4. Reproduction and scope

```text
python3 -u -c 'import runpy; runpy.run_path("04-computation/overnight3_20260906_lrc_consumers_audit.py", run_name="__main__")'
```

Normal and optimized executions give identical output with 4,155 explicit
gates. The complete literal replay takes only a few seconds because exact
rational predicates are represented by integer grid numerators. There are
no floating-point constants or optimization-disabled assertions.

SHA-256 values over LF bytes:

```text
source:
87f1716b869cfdbe6d1ade4cf609b3d5dd8288112effe03402d86ea9a942f39c
output:
7023bbf288208ba46baaef85d894201e08dd2e0c69ffe4ec7f64266c5fe12ab7
semantic digest:
fe3e202c4b66aa7ae6c578adbbb4b40ef5a58d11792d36a6d516bd0ca38434f1
```

No producer files were changed. The accepted conclusions retain parity,
complete component addresses, and weak endpoints. They establish neither
LRC(14) nor a universal mixed-parity ceiling of `6/55`.
