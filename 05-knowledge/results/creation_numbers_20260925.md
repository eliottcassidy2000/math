# The three fruit integers: literal audit, ratio bounds, and signed certificates

Status: **PROVED** for elementary algebra and carry identities; **FINITE-EXACT** for the recorded integer trajectories and elliptic-point controls; **CITED** for the primary-source provenance. Date: 2026-09-25. No universal Collatz convergence claim or global leastness claim is made.

## 1. Inheritance and literal preservation

Closest proved mechanism: [catalan_elliptic_20260921_elliptic.md](catalan_elliptic_20260921_elliptic.md) gives the symmetric fruit cubic, its exact elliptic correspondence, and the corrected primitive positive triple below. The prior transcription repair there concerns a different paste; it is not silently applied to this turn. [collatz_mod6_20260921_fruit_rank_positive_multiples.md](collatz_mod6_20260921_fruit_rank_positive_multiples.md) supplies the bounded positive-multiple census and exact group operations. [thirtysix_digits_20260925.md](thirtysix_digits_20260925.md) supplies the signed inverse-germ clock and its warning that an inverse-fibre index is not forward time.

Canonical hostile: a real ratio or an elliptic coordinate does not determine a Collatz valuation word. Corrected near miss: numerically similar decimal integers can change the cubic equation. Least-used sidecars: raw decimal input, initial powers of two, the entire halving word, and the selected positive elliptic component. Live board: projective fruit equality; normalized dominance ratios; exact binary carries; forward common futures; elliptic positive-point controls.

The literal input strings, without any correction, are:

```text
15,447,680,210,874,616,644,195,131,501,991,983,748,566,432,566,956,543,170,002,663,489,825,320,203,527,7999
4,373,612,677,928,697,257,861,252,602,371,390,152,816,537,558,161,613,618,621,437,993,378,423,467,772,036
368,751,317,941,299,998,271,978,115,652,254,748,254,929,799,689,719,709,962,831,374,716,372,246,340,555,799
```

Removing commas only gives `(a,c,d)`, of digit lengths `(81,79,81)`. In particular the first string's final group `7999` is retained. The relation to the previously verified triple is

\[
d=10b+9,\qquad b=(d-9)/10,
\]

where

```text
a = 154476802108746166441951315019919837485664325669565431700026634898253202035277999
b = 36875131794129999827197811565225474825492979968971970996283137471637224634055579
c = 4373612677928697257861252602371390152816537558161613618621437993378423467772036
```

Thus the three literal slots and three repaired slots contain **four** distinct integers. This note computes both sheets for every one of those four integers; it does not replace the literal third input without retaining its certificate.

For `s=a+b+c`, put

\[
F_N(a,b,c)=s^3-(N+2)s(ab+ac+bc)+(N+3)abc.
\]

For nonzero pair sums,

\[
\frac{a}{b+c}+\frac{b}{a+c}+\frac{c}{a+b}-N
=\frac{F_N(a,b,c)}{(a+b)(a+c)(b+c)}.
\tag{1}
\]

Exact integer arithmetic gives `gcd(a,b,c)=1` and `F_4(a,b,c)=0`. In contrast, `F_4(a,c,d)<0`, and the literal fruit sum is approximately

`2.743741758121343036518043593454996850803`.

The output retains the full integer residual and exact rational sum. The repaired fruit sum is exactly 4.

The [2014 Bremner–Macleod paper](https://ami.uni-eszterhazy.hu/uploads/papers/finalpdf/AMI_43_from29to41.pdf), section 2 and Remark 2.2, gives the elliptic family and the ninth-multiple example. Its [2025 corrigendum](https://publikacio.uni-eszterhazy.hu/8862/1/AMI_62_from26to27.pdf) corrects numerical table entries and retains the `N=4`, `m=9`, maximum-81-digit row; it explicitly identifies the three coordinate lengths as 79, 80, and 81. These primary documents were rechecked for this note. Our exact equality and trajectory certificates do not rely on a table-wide minimality assertion.

## 2. What the ratios do follow from

The corrected ratios are

\[
\frac bc=8.4312751287326437696965788642056\ldots,
\qquad
\frac ab=4.1891864406390187213487583959752\ldots.
\]

For the literal third input, `d/c=84.312751287326437696965788642056...`. Its change in scale is due to the explicit decimal relation `d=10b+9`.

There is a stronger analytic constraint than comparing those ratios to rounded constants. For any positive real solution of fruit sum 4, choose `a` as a largest coordinate and define

\[
r=\frac{a}{b+c},\qquad t=\frac{bc}{(b+c)^2},\qquad 0<t\le\frac14.
\]

Then the fruit equation is exactly

\[
4=r+\frac{r+1-2t}{r^2+r+t},\qquad
t=\frac{r^3-3r^2-3r+1}{6-r}.
\tag{2}
\]

Since the last fraction in the first equation decreases as `t` increases,

\[
r+\frac{2}{2r+1}\le4<r+\frac1r.
\]

Also `r>=1/2` because `a` is largest. Solving these inequalities therefore gives

\[
2+\sqrt3<\frac{a}{b+c}\le\frac{7+\sqrt{65}}4.
\tag{3}
\]

The upper equality occurs at `b=c` for positive real solutions. For rational or integer triples it cannot occur because the upper endpoint is irrational. The lower endpoint is approached when one of the smaller coordinates tends to zero relative to the other and is never attained by a positive triple. Numerically the interval is approximately `(3.73205,3.76556)`. The displayed solution has

\[
r=3.745006159239259220497607804489896\ldots.
\]

This is the exact reason its largest coordinate is close to four times the sum of the other two. It does not force `b/c` to be 8.431...: for real solutions the split parameter `b/c` varies, with `t=(b/c)/(1+b/c)^2`, and (2) gives the corresponding `r`. Retaining the split is necessary; collapsing to one approximate ratio loses the point on the cubic.

The constant `3/2` has another precise role here. For any positive triple, Cauchy–Schwarz gives

\[
\sum\frac{a}{s-a}=s\sum\frac1{s-a}-3
\ge s\frac9{2s}-3=\frac32,
\]

with equality exactly at `a=b=c`. This is a lower bound for a symmetric rational expression, not a transfer of the odd Collatz multiplier `3/2` into an elliptic evolution rule.

## 3. Exact trajectories on both sheets

Use the ordinary clock

\[
C_\sigma(n)=\begin{cases}n/2&n\text{ even},\\3n+\sigma&n\text{ odd},\end{cases}
\qquad \sigma\in\{+1,-1\}.
\]

For odd inputs use the accelerated odd clock

\[
U_\sigma(n)=\frac{3n+\sigma}{2^{v_2(3n+\sigma)}}.
\]

The first column below records the initial ordinary halvings needed before the odd clock starts. Odd preperiod counts stop at the first entry into the displayed odd cycle. Ordinary times to 4 refer to the original integer input.

| input | initial halvings | plus odd steps to 1 | plus ordinary steps to 4 | minus odd preperiod | minus eventual cycle minimum |
|---|---:|---:|---:|---:|---:|
| a | 0 | 681 | 2025 | 613 | 17 |
| b | 0 | 717 | 2116 | 548 | 5 |
| c | 2 | 608 | 1831 | 598 | 1 |
| literal d=10b+9 | 0 | 585 | 1778 | 633 | 5 |

The minus odd cycles, written from their minima, are

\[
(1),\qquad(5,7),\qquad(17,25,37,55,41,61,91).
\]

The certificates close each actual trajectory on its cycle and check the closing edge; they do not assume that these are all possible minus cycles. The plus trajectories supply genuine finite routes to the ordinary root 4 for all four numbers.

The first common odd states, with odd-clock indices measured from the odd parts of the inputs, include:

| sheet | pair | first common state | indices |
|---|---|---:|---|
| plus | a, b | 317 | 670, 706 |
| plus | c, d | 53 | 606, 583 |
| plus | a, c | 5 | 680, 607 |
| plus | a, d | 5 | 680, 584 |
| minus | b, d | 13 | 546, 631 |

All pairwise joins are retained in the output. Corrected `a,b,c` have no common future with one another on the minus sheet, because the complete certificates close on different cycles. This is a statement about these explicit inputs, not about the fruit curve as a whole.

## 4. Cheap hostile to the three-basin pattern

The corrected `9G` triple happens to realize all three displayed minus basins. The next inherited positive examples already disprove the suggestion that every positive fruit solution does so. On

\[
E:y^2=x^3+109x^2+224x,
\quad G=(-4,28),\quad T=(56,728),
\]

use the exact projective inverse from the inherited note, clear denominators, and divide by the common gcd. The resulting controls are:

| elliptic point | coordinate digit lengths | minus basin minima | odd preperiods |
|---|---|---|---|
| 13G+T | 167,168,167 | 5,5,17 | 1320,1303,1128 |
| 17G | 286,286,286 | 5,5,1 | 2164,2395,2051 |

All six coordinates and all six valuation certificates are stored in the output. Independently of the elliptic construction, the verifier checks their positivity, primitivity, and exact fruit sum 4. Thus a direct proposed implication “positive fruit solution implies one coordinate per minus basin” is **REFUTED** by `13G+T`. The strongest survivor is the specific finite basin assignment for the displayed `9G` triple. The first missing ingredient is an invariant connecting the elliptic point to forward Collatz itineraries; the fruit equality itself supplies none of that itinerary data.

## 5. The useful binary carry coordinate

Approximate real ratios do not encode exact halving counts. The actual initial binary coordinates are:

| input | odd part modulo 256 | v2(n_odd+1) | v2(n_odd-1) | initial plus one-halving run | initial minus one-halving run |
|---|---:|---:|---:|---:|---:|
| a | 175 | 4 | 1 | 3 | 0 |
| b | 155 | 2 | 1 | 1 | 0 |
| c | 33 | 1 | 5 | 0 | 4 |
| d | 23 | 3 | 1 | 2 | 0 |

Here `n_odd` means the input after its initial powers of two are removed. The row for `c` therefore starts at `c/4`.

For any positive odd `n>1` and fixed `sigma`, put `h=v_2(n+sigma)` and write `n=2^h u-sigma`, with `u` odd. Then for `0<=j<=h-1`,

\[
U_\sigma^j(n)=3^j2^{h-j}u-\sigma.
\tag{4}
\]

The first `h-1` steps have exact halving exponent one; the next exponent is at least two. This follows because, while `h-j>=2`,

\[
3(3^j2^{h-j}u-\sigma)+\sigma
=2(3^{j+1}2^{h-j-1}u-\sigma),
\]

and the parenthesized expression is odd. At `h-j=1` the corresponding parenthesis is even. These one-halving steps grow when the input is greater than one. Equation (4), checked against each recorded prefix, is an exact finite carry decoder. It does not control the itinerary after that run or turn the rounded ratios 4 and 8 into branch labels.

In the inverse-germ coordinates of [thirtysix_digits_20260925.md](thirtysix_digits_20260925.md), the initial plus ranks `v_2(3n_odd+1)` are `(1,1,2,1)` for `(a,b,c,d)`, and the minus ranks are `(2,4,1,2)`. These ranks describe a single edge pointing to its odd target. They differ from the prefix lengths in the table, from ordinary time, and from the elliptic multiplier 9.

## 6. Certificate format and reproduction

```text
python 04-computation/experiments/creation_numbers_20260925.py
python -O 04-computation/experiments/creation_numbers_20260925.py
python 04-computation/experiments/creation_numbers_20260925.py --verify
python -O 04-computation/experiments/creation_numbers_20260925.py --verify
```

The `.out` is machine-readable JSON. It stores raw strings, literal and repaired integers, exact fruit residuals, ratios, joins, the two hostile positive triples, and fourteen complete trajectory certificates. Each certificate contains its start, sign, initial halving count, exact valuation word, preperiod, and terminal cycle. The generator uses binary valuations; the separate verifier anchors the starts to the raw strings and stated repair, replays every ordinary numerator step and every individual halving, checks parity at each one, rejects premature repeated odd states, and checks the cycle closure. The `--verify` route reads only the saved data and does not run the trajectory generator or elliptic operations. Altering the first valuation is an explicit rejection control. All checks use exceptions and survive `-O`.

The bounded universe is four distinct inputs on both sheets plus six hostile coordinates on the minus sheet, each with a limit of 100000 odd steps. Every trajectory closed well inside that bound. No inference from a timeout is used. The inherited group-law module is imported only for the two control triples; its `main` and PARI computations are not called. Further curve-wide basin sampling is unnecessary to refute the universal three-basin proposal. A future structural claim would need an explicit map preserving the forward certificate predicate, including parity and valuation data, rather than just projective fruit equality.
