# Independent audit: the nonadditive 11/140 cutoff consumer

**PASS ANALYTICALLY + FINITE-EXACT; inherited coefficient and low-sector proofs.**
This report audits a sufficient all-height reduction and its finite-universe
count. It does not repeat or independently certify the native physical census.
While this audit was running, incoming `058a8ded9` promoted THM-4437 (the
generic all-parity 6/77 bound outside the three low circuits) and THM-4441
(the sharp norm-five selected bound 46/665). Together with the fourth-round
norm-four theorem, these give the nonadditive 11/140 result more directly.
The reduction below remains an independently checked reusable consumer,
not a new priority claim or a necessary premise of that shorter proof.

## Scope and inherited mechanism

Let `1<=a<b<c`, `gcd(a,b,c)=1`, all three speeds be units modulo three,
and `c!=a+b`. Physical mass and every selected projection satisfy
`mu<=min E_i<=qN`, where `q=3/(7c)` and N counts complete owner-live raw
carriers. Thus `N<11c/60` suffices for strict `min E_i<11/140`.

The frozen fourth-round proof supplies a primitive relation v with
`S=||v||_1<4 sqrt(c/3)`, with no parity assumption. Its slice bound for a
nonempty admissible defect list is

```
N < F_v+B_v,   B_v<=2S/7+4/3,
F_v/c<=6/49+4/(7M),   M=max |v_i|.
```

The audited finite coefficient table covers exactly 747 patterns with
M<=18: 698 with three nonzero entries and 49 with an actual zero.
Its full ordered semantic digest is checked before using any record.
The source of those records is the previously independent error-polygon
referee, not a fresh claim that this audit re-proves its geometry.

The canonical hostile is the empty-defect pattern `(0,1,2)`, where N=0
directly. It must not pass through the strict inequality `0<0`. Among
the excluded primitive patterns, `(0,1,1)` forces repeated speeds;
`(1,1,1)` forces `c=a+b` and is excluded by the present nonadditive
hypothesis. Patterns `(1,1,2)` and `(1,2,2)` retain the norm-four and
norm-five roof sidecars. Coordinate permutations and global sign change
do not lose a positive-speed relation; a primitive relation has at most
one zero residue modulo three. These are exactly the inherited coverage
conditions, not a newly filtered coefficient universe.

## Exact cutoff audit

For every nonempty small-coefficient record `(p,alpha,B)`, the necessary
sufficient-consumer threshold is `B/(11/60-alpha)`. The maximum is exactly

```
35280/199 <178,
unique pattern (7,13,18), alpha=17/147, B=12.
```

All 746 nonempty records have positive denominator. Therefore c>=178
already gives strict N<11c/60 throughout this branch. The sole empty
record was separately discharged above.

For M>=19, put

```
alpha=142/931,
gap=11/60-alpha=1721/55860,
A=(2/7)/gap=15960/1721,
B=(4/3)/gap=74480/1721.
```

It suffices that c>=AS+B. If S<=53 and c>=535, then
`AS+B<=920360/1721<535`. If S>=54, the short-relation inequality gives
`c>3S^2/16`; furthermore

```
g(S)=3S^2/16-AS-B,
g(54)=18547/6884>0,
g'(54)=75561/6884>0,
g''(S)=3/8>0.
```

Hence c>AS+B also in this branch. The generic sector is strict for every
c>=535. The value `920360/1721` exceeds 534, so silently rounding this
sufficient cutoff down would be invalid. A head through c=534 suffices;
the H535 native census deliberately includes one redundant boundary layer.

The inherited norm-five bound uses continuum selected mass at most
`363/6272` and one-sided residue error `4/(7c)`. Its exact threshold is
`(4/7)/(11/140-363/6272)=17920/649<28`, leaving the already checked H27
head. The newer sharp 46/665 result is strictly below 11/140. The frozen
norm-four proof has `mu=min E_i<=11/140`, with unique equality `(2,11,20)`;
its tail c>=33 is strict, with gap 1/32340 at 33. No new assertion about
the additive sector, arbitrary entry, or LRC(14) follows from this audit.

## Independent universe count

Set `U(H)=H-floor(H/3)` and let mu(d) denote the arithmetic Mobius function.
Gcd inclusion-exclusion gives the number of primitive strict ternary-unit
triples through H as

```
T(H)=sum_(1<=d<=H,3 does not divide d) mu(d) binom(U(floor(H/d)),3).
```

Divisors divisible by three cannot occur in a unit-speed triple; scaling
by every remaining divisor preserves the unit test and strict order.

For raw additive pairs, the unit condition forces a and b into the same
nonzero residue r modulo three. Write `a=r+3i`, `b=r+3j`, `0<=i<j`, and
`J=floor((H-2r)/3)`. The count of `i+j<=J` is `floor((J+1)^2/4)` for
J>=0 and zero otherwise (sum the available j for each i). Summing r=1,2
defines A0(H). The primitive additive count is

```
A(H)=sum_(1<=d<=H,3 does not divide d) mu(d) A0(floor(H/d)).
```

The independently reconstructed counts are

| head | primitive unit triples | additive | nonadditive |
|---|---:|---:|---:|
| 534 | 6,453,474 | 10,838 | 6,442,636 |
| 535 | 6,514,176 | 10,908 | 6,503,268 |

The redundant c=535 layer contains 60,632 nonadditive triples. Controls
compare the Mobius formulas with direct count-only enumeration at every
height 1..80, compare A0 with an independent residue count at every height
1..535, and count primitive additive pairs directly by gcd at H535.
No physical roofs, owner contacts, or projection masses are recomputed.
The resulting universe matches the native output exactly. Its claimed
mass maxima still require the separately assigned raw/native comparison.

## Reproduction and freeze

From the repository root:

```
python -B 04-computation/overnight5_20260906_lrc_nonadditive_reduction_audit.py
python -B -O 04-computation/overnight5_20260906_lrc_nonadditive_reduction_audit.py
```

Normal and optimized outputs are byte-identical LF: **2,661 explicit
optimization-live gates PASS**, apart from the inherited polygon routine's
own gates. The source intentionally imports only the frozen polygon
referee routine as coefficient input; the new consumer and all universe
calculations are independent of both native census and prior head loops.

```
source 695f2081bbb92fdf386af29efb54cc862853ceb411b4bc26d9cda55821afa6ed
output fb1c1887e68c242f94958103cac754cec0fe69de061e9226810dd43db857eccd
inherited referee b5f52380e28bd3883f95090ff0c06a89e667d7aafd0e4c241b8fca68c76c7c00
coefficient records cf808062354debbefc1d8ead8ad0d10e9da5427cb42b8f083b6af24d0059c87c
native output e8142deeaea3278b6a07b07676570205fd7852c366e34088ec0973e597dc00a1
```

Source-target-loss contract: the source is the complete owner-live carrier
dictionary. Relation slicing and its frozen convex-width table map it to a
count upper bound, preserving an upper bound on every projection but losing
the individual roof values. The norm-four/norm-five sidecars recover those
values precisely where the generic count consumer is too weak. The new
operation changes the target constant and optimizes each inherited slope
with its own intercept, instead of discarding their correlation.

**Filing note.** Root filed this completed audit after checkpoint `07b2d91b2`.
Local paths were made repository-relative and output line endings normalized;
normal and optimized verification was rerun with matching output. The proof
and finite universes are unchanged.
