# The supplied sequence is a machine orbit; its two different trees matter

2026-10-26. **CITED:** the supplied sequence is the Marxen–Buntrock
five-state Busy Beaver machine's arithmetic macro-orbit. **Elementary
proofs below:** a normalization, exact integer predecessor forest,
survivor-cylinder count, periodic-address bound, and divisor valuation
identity. **FINITE-EXACT:** all bounded controls in `check.py` pass normally
and under `python -O`, with identical output. The
[independent review](INDEPENDENT_REVIEW.md) passes in these stated scopes.
No Collatz or universal machine-halting claim is made.

## 1. Exact origin, including its missing clock

The list is

```
0,6,16,34,64,114,196,334,564,946,1584,2646,4416,7366,12284,halt.
```

It occurs explicitly in [Aaronson, The Busy Beaver Frontier, p.10](https://eccc.weizmann.ac.il/report/2020/115/download/).
The primary configuration analysis is [Michel, Behavior of Busy Beavers,
Marxen–Buntrock winner](https://bbchallenge.org/~pascal.michel/beh),
which cites Buro (1991) and Michel (1993). It defines
`C(x)=...0(A0)1^x0...` and gives

| Input | Next configuration | Exact machine steps |
|---|---|---|
| x=3k | C(5k+6) | 5k²+19k+15 |
| x=3k+1 | C(5k+9) | 5k²+25k+27 |
| x=3k+2 | halt | 6k+12 |

Thus the nonhalting branches are `(5x+18)/3` and `(5x+22)/3`.
Summing these cited clock formulas on the displayed orbit gives exactly
47,176,870 machine steps; the final pattern has 4098 ones. This packet
checks that sum, not a new full machine simulation. Dropping the clock
would preserve reachability but destroy the Busy Beaver objective.

An exact-sequence repository search did not recover an earlier copy of
this orbit in the searched canon/result text. Its mechanisms do have
substantial ancestry below. An attempted exact OEIS lookup was not
accessible through the search interface; no OEIS identifier is invented
and no absence-from-OEIS claim is made.

## 2. A small normalization with exact arithmetic

Parity is preserved on each nonhalting raw branch. In the even sector,
put `z=x/2+5`. The original nonnegative sector becomes `z>=5`. Extend the
same formula to positive integers for a simpler structural statement:

```
F(z)=(5z+1)/3  if z=1 mod3;
F(z)=(5z-1)/3  if z=2 mod3;
F(z) undefined if z=0 mod3.
```

Equivalently, on its domain `F(z)=floor((5z+1)/3)` and
`3F(z)-5z=epsilon`, where epsilon is +1 or -1. The normalized source orbit
is

```
5,8,13,22,37,62,103,172,287,478,797,1328,2213,3688,6147,halt.
```

The floor formula must retain the domain guard: using it also at multiples
of three would create a different dynamical system.

### Proposition A: the actual integer graph is a forest of simple paths

The exact image of F on positive integers not divisible by three is

```
{y>0 : y=2 or3 mod5}.
```

Every such y has its unique predecessor

```
(3y-1)/5 if y=2 mod5;
(3y+1)/5 if y=3 mod5.
```

**Proof.** Reduce `3y=5z+epsilon` modulo five. The two signs give disjoint
image residues 2 and 3. Conversely, the displayed inverse is a positive
integer, and modulo three it is respectively 1 or 2, so its branch guard
holds. Every F step strictly increases z; every inverse step strictly
decreases y. Repeated inversion therefore reaches a unique positive root
whose residue modulo five is 0,1 or4. Every component is a simple finite
or infinite path, with no cycle or merger. In particular, 5 is a root.

This is not the binary tree described next. That tree classifies possible
residue itineraries across different sources, while the actual integer
inverse graph has at most one predecessor at each vertex. Ordinary
Collatz has actual merging; a correspondence based only on the word
“tree” cannot transport its orbit problem.

### Proposition B: the symbolic survivor set is an exact 3-adic Cantor set

For a sign word `epsilon_0,...,epsilon_(k-1)`, define

```
A_k = sum_(j=0)^(k-1) 5^(k-1-j) 3^j epsilon_j.
3^k z_k = 5^k z_0 + A_k.                         (1)
```

That word is realized by exactly one residue class of z_0 modulo `3^k`.
There are exactly `2^k` such classes, and their natural/Haar fraction is
`(2/3)^k`.

**Proof.** The inverse branches `h_e(y)=(3y-e)/5` are well defined on
Z_3. They contract by 1/3 and map onto the disjoint residue balls 1 and 2
modulo three. Repeated inverse composition gives one ball per word;
different words first separate at their first different sign. Alternatively,
(1) gives its representative `-A_k*(5^k)^(-1) mod3^k`; successive
integrality conditions enforce exactly the stated signs. Each surviving
ball has three refinements. Since `z_k(z_0+t3^k)=z_k(z_0)+t5^k`, precisely
one refinement halts on its next test and precisely two survive.

The infinite survivor set S is consequently conjugate to the full binary
one-sided shift by its itinerary. It has Haar measure zero and Hausdorff
dimension `log(2)/log(3)`: the `2^k` depth-k balls give the upper cover,
and assigning mass `2^-k` to each cylinder gives the matching lower bound.
The set of positive integers that never halt has natural density zero,
because it lies in every finite survivor union. This does **not** prove
that that set is empty, or that a specified integer is outside it.

Every nonempty finite survivor cylinder contains positive integers that
halt on the very next test, and positive integers that survive any specified
additional finite number of tests. Thus no such cylinder alone certifies
nontermination. This is not an impossibility theorem for adaptive,
source-labelled, or height-dependent proofs.

### Exact source hostile: one more digit changes the next event

The sources below have identical first fourteen signs:

| Source z | Its actual integer root | F^14(z) | Next test halts? |
|---:|---:|---:|---|
| 5 | 5 | 6147 | yes |
| 4,782,974 = 5+3^14 | 4,782,974 | 6,103,521,772 | no |
| 14,348,912 = 5+3^15 | 3,099,365 | 18,310,553,022 | yes |

The cylinder `5 mod3^15` certifies this fourteen-step halting itinerary;
`5 mod3^14` alone does not. Even matching finite signs and final halting
time does not preserve the actual root. This is the precise missing
coordinate in a proposed transfer to one fixed Collatz source such as 27.
As a separate numerical hostile, literal z=27 halts immediately here and
belongs to root16; sharing the label 27 is no dynamical correspondence.

### Proposition C: periodic symbolic survivors are not positive integers

A repeated sign word of length k has a unique periodic 3-adic point

```
z = -A_k/(5^k-3^k),       |z|_infinity <= 1/2.     (2)
```

**Proof.** Set z_k=z_0 in (1). The denominator is prime to three, so this
is a valid 3-adic point with the prescribed itinerary. Moreover,
`|A_k| <= sum_j 5^(k-1-j)3^j = (5^k-3^k)/2`.
Uniqueness of the infinite 3-adic coding implies that any point with that
periodic itinerary is this fixed point. A finite preperiod applies the
inverse maps h_e, each of which preserves the real interval `[-1/2,1/2]`.
Hence no positive integer survivor can have an eventually periodic sign
itinerary. A divergent positive-integer survivor, if one exists, must be
aperiodic. Excluding cycles alone does not exclude it.

There is a useful completion warning: the same two inverse affine maps
contract in the real metric as well, by 3/5. Their real images of
`[-1/2,1/2]` are `[-1/2,1/10]` and `[-1/10,1/2]`, which overlap and cover
the whole interval. In Z_3 the images are disjoint. The real and 3-adic
limits of a nonperiodic rational approximation sequence need not be the
same rational number. Treating them as equal would falsely finish the
positive-integer problem. Formula (1), with its denominator and actual
terminal integer, is the sidecar that must remain.

## 3. Sloane's ternary problem: a real bridge with a different endpoint

[Saye, JIS 25 (2022), Article 22.3.4](https://cs.uwaterloo.ca/journals/JIS/VOL25/Saye/saye3.pdf)
distinguishes Erdős's “no digit 2” question from Sloane's “no digit 0”
question. For Sloane the exceptional exponents are conjecturally
`0,1,2,3,4,15`; Saye verifies both conjectures through exponent
`2*3^45`. We use his exact lift mechanism, not a probabilistic conclusion.

The order of 2 modulo `3^k` is `2*3^(k-1)`. Indeed,
`v_3(4^m-1)=1+v_3(m)`, by the elementary lifting-the-exponent identity.
Consequently exponents modulo that order biject with all units modulo
`3^k`. Exactly `2^k` of those units have k trailing digits in `{1,2}`.
Each zeroless exponent class has precisely two good lifts among its
three next-level lifts. The relative fraction is `(2/3)^(k-1)`; the
initial parity/two-root factor explains the difference from Proposition B.

The abstract binary refinement pattern is shared with F. The operation
is not: shifting the exponent multiplies the value by 2, whereas F
consumes one ternary residue while multiplying by 5/3 and adding a carry.
A tree isomorphism preserves refinement and measure, not these operations.

Most importantly, Sloane's finite representation stops at its actual
leading digit. Padding `2^15` with a leading zero does not refute its
zerolessness. A candidate must be tested with the depth bound
`k<=floor(n log_3(2))+1`. The word length is a genuine archimedean coordinate.
This recovers the earlier repository warning in
[THM-4580, decimal zeroless powers](../../../../01-canon/theorems/THM-4580-zeroless-powers-of-two-the-5-adic-tree-and-verification-to-1-1e11.md)
and the [September sibling-dimension atlas](../../collatz_procgen_20260922_sibling_dimension_ladder.md):
fixed residue features and a thin exceptional set do not by themselves
settle the selected arithmetic trajectory.

## 4. Multiply-perfect numbers give an exact conservation law

For `n=product p^a`, define the prime valuation output of each prime-power
block by `v_q(1+p+...+p^a)`. Then, for a positive integer k,

```
sigma(n)=kn
iff, for every prime q,
sum_(p^a||n) v_q(1+p+...+p^a) = v_q(n)+v_q(k).     (3)
```

This follows from the multiplicative formula for sigma and unique prime
factorization. It is a finite weighted flow diagram, with loops and
multiple contributions permitted; making it a tournament would discard
its actual relation. The source is the list of prime-power blocks, the
map is factorization of each block's divisor sum, and the preserved
predicate is the entire equality sigma(n)=kn. Dropping exponents or
prime labels destroys it.

Controls are `n=6,28,120,30240`, with multipliers `2,2,3,4`.
Triperfect120 has aliquot image240, so multiply-perfect does not mean
aliquot fixed except at k=2. The older perfect/aliquot carrier is
[HYP-2221](../../../hypotheses/HYP-2221-perfect-aliquot-fixed-point-carrier.md);
its near-perfect odd control32445 has `sigma(n)-2n=6`. Small scalar defect
does not replace an exact valuation balance.

This suggests a construction, not a Collatz proof: assign to an excursion
certificate a labelled ledger of consumed and regenerated prime-adic
precision, and require conservation at every prime together with the
original source embedding. Formula (3) shows what a successful scalar-to-
vector refinement looks like in another problem. A positivity or descent
law for that Collatz ledger remains to be proved.

## 5. Reproduction, inheritance, and next targets

Run from the repository root:

```
python 05-knowledge/results/root27_machine_bridges_20261026/atlas/check.py
python -O 05-knowledge/results/root27_machine_bridges_20261026/atlas/check.py
```

[Normal](normal.json) and [optimized](optimized.json) outputs are identical.
The explicit universes are all residues through modulus3^10, all periodic
sign words of lengths1–9, positive integer forest controls through20000,
Sloane residue periods through10 digits and actual exponents0–512, plus
the named divisor and interpolation hostiles. Algebraic cylinder generation
is compared with direct iteration. These bounds validate the implementation;
the universal statements use the proofs above.

The closest mechanism is exact guarded affine coding and the inherited
zeroless lift tree. The hostile is the three same-prefix sources above.
The corrected near miss is mistaking an itinerary tree for the integer
predecessor graph, or confusing a padded zero with a represented digit.
The least-used sidecar is Michel's quadratic machine clock. The
[25-node historical atlas](HISTORY_ATLAS.md) reaches March and June
threads, not only the recent dominant route.

Anchor: preserve one arithmetic source through a certificate. Niche:
source-root and cylinder-tree disagreement. Wildcard: exact prime-flow
ledgers and modular-level defects. The live board has six objects:
literal source; affine carry; refinement depth; operation/clock; integer
root; prime valuation ledger. New analogies are retained only with an
explicit map and a cheap hostile test.

The strongest next question is concrete: can an orbit-dependent rank use
both the current affine carry and the original source-root information
to exclude indefinite certificate renewal? Proposition B rules out one
fixed-depth shortcut; it does not rule out that richer construction.
