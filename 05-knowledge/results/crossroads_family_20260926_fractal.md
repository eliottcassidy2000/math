# Occurrence rates, recursive graph carriers, and the missing ambient metric

**Status: PROVED elementary consequences of the independently audited
THM-4495; independently audited SOUND by the automata lane. FINITE-EXACT controls are
separate.** The dimension exponent is inherited, not a priority claim.
No positive-integer Collatz convergence claim follows.

## 1. Inheritance and the object

The closest mechanism is the parity-cylinder bijection and the exact
no-descent word order in [THM-4495, no-descent count](../../01-canon/theorems/THM-4495-no-descent-count-exact-order-spitzer-ballot.md).
The dimension exponent already appears in [the choice ladder](collatz_procgen_20260922_choice_ladder.md)
and [the strategy cube](procgen_cube_20260925_strategy_cube.md).
The corrected near miss is treating a dyadic exceptional set as a set of
positive integers. The hostile is a single negative periodic point, or the
two-completion example below. The least-used sidecar is the absolute
positive-integer realization of an infinite address.

Use the shortcut map T(n)=n/2 for even n, (3n+1)/2 for odd n, extended to
Z_2. A length-k binary parity word w has affine action

    T_w(n)=(3^a n+C_w)/2^k,

where a is its number of odd letters and C_w is a nonnegative integer.
It describes exactly one residue modulo 2^k. Two starting values have the
same first k parity bits if and only if they are congruent modulo 2^k:
induct on the first differing bit, since each branch divides the difference
by 2 and multiplies it by an odd unit. Consequently infinite parity coding
is an isometry between Z_2 and binary addresses with their dyadic metric.

Let W_k count words for which every nonempty prefix of length j has
3^(number of odd letters)>2^j. Let E be their infinite-address set in Z_2.
All allowed finite words extend by infinitely many odd letters, so the
number of occupied level-k cylinders is **exactly W_k**.

Put rho=log_3(2) and h=H_2(rho)=0.9499555... . THM-4495 proves

    W_k = Theta(2^(h k) k^(-3/2)).                     (1)

For fixed k, the positive integers in these cylinders have natural density
W_k/2^k. The same fixed-k density holds for actual no-descent through k
steps: a word with a contracting prefix has T^j(n)=a_j n+c_j with a_j<1,
so only finitely many positive sources on that cylinder avoid a drop;
there are only finitely many words at fixed k. The finite exceptions cannot
be discarded uniformly in k. In particular n=1 never falls below itself,
although one of its prefix multipliers is less than 1.

## 2. From frequency to a precise fractal statement

**Proposition.** In the dyadic metric, E has Hausdorff dimension and both
box dimensions equal to h, Haar measure zero, and h-dimensional Hausdorff
measure **zero**. Its exact dyadic covering numbers obey (1).

*Upper bound and critical measure.* The occupied level-k cylinders cover
E and have diameter 2^(-k). Their h-cost is
W_k 2^(-h k)=O(k^(-3/2)), tending to zero. This proves H^h(E)=0 and
dimension at most h. Their Haar measure is W_k/2^k, also tending to zero.
Equation (1) directly gives both box dimensions.

*Lower bound by a small graph.* Fix b. Form a graph with **one vertex and
W_b labelled loops**, one for each allowed b-bit word. A loop w acts by

    phi_w(y)=(2^b y-C_w)/3^a.                           (2)

These maps contract dyadic distances by 2^(-b), have disjoint images in
their different parity cylinders, and insert their indicated word before
the tail address. Concatenating allowed blocks stays allowed: the surplus
a log 3-j log 2 is positive within a block and accumulates positively
between blocks. Thus this graph's infinite paths encode a subset E_b of E.

The uniform independent block measure assigns mass W_b^(-m) to every
m-block cylinder. If s_b=log_2(W_b)/b, a dyadic cylinder of length mb+r,
0<=r<b, has mass at most W_b^(-m)<=2^(b s_b)2^(-s_b(mb+r)). The elementary
mass bound gives dim_H(E_b)>=s_b; block covers give the reverse inequality.
By (1), s_b tends to h. Taking b arbitrarily large proves dim_H(E)>=h.
QED. The b=1 case has dimension zero and causes no exception.

The same argument shows that every nonempty allowed prefix cylinder in E
has dimension h: append arbitrary allowed blocks to that prefix. This is
local recursion at every occupied scale, not merely a picture resembling
a fractal. The polynomial factor in (1) explains why the critical
h-dimensional measure vanishes despite dimension h. No claim of a positive
measure for a logarithmically corrected gauge is made.

The graph retains exact parity words and their affine carries through its
edge labels. Counting its loops alone loses those labels. The one-vertex
description also uses many loops as b grows; it is not a bounded finite
automaton recognizing the entire prefix language.

## 3. A deliberately hostile two-loop graph

The word maps (2) contract dyadically while their real behavior depends on
3^a/2^b. Even when they contract in both completions, real sign information
does not by itself exclude positive dyadic integers. Here is an explicit
control, **not the Collatz map**:

    phi_1(y)=(2y-1)/3, phi_2(y)=(2y-2)/3.              (3)

In the real metric their unique compact attractor is [-2,-1]: its images
are [-5/3,-1] and [-2,-4/3], whose union is the interval. In Z_2 their images
are respectively the odd and even residue balls, disjoint and covering
Z_2. The same two labelled loops therefore have dyadic attractor all Z_2.

For a concrete common address start x_0=1. Set d_j=1 when x_j is odd and
2 when it is even, and x_(j+1)=(3x_j+d_j)/2. Every x_j is a positive
integer. The rational prefix approximant A_m=phi_(d_0)...phi_(d_(m-1))(0)
satisfies exactly

    A_m=1-(2/3)^m x_m.                                (4)

In Z_2 its difference from 1 has valuation at least m, so A_m tends to the
positive integer 1. In the real metric the contractions give convergence
to a point in [-2,-1]. One rational sequence, with one graph and one
address, has different limits in the two completions. Any argument that
excludes positive integers just because a real attractor is negative must
identify and prove an extra compatibility condition.

For an allowed periodic Collatz word the fixed point is
C_w/(2^b-3^a)<0. That sign observation is correct, but applying it to every
aperiodic dyadic address is exactly the unsupported inference exposed by
(3). Large dimension or an everywhere-negative real limit is insufficient
to settle intersection with ordinary positive integers.

## 4. Comparison with the explicit family and square paths

[The orbit-family note](crossroads_family_20260926_orbits.md) constructs
disjoint positive progressions F_k with k+11 no-descent steps and a retained
three-prime triangle inside ten odd nodes. Their tail-union density is
2^(1-K)/19124224. The phase clock gives exact affine recursion, but the
nested closed tail families intersect only at -1. These families realize
arbitrarily long finite growth without one positive survivor at all scales.

By contrast E has dimension h and uncountably many dyadic survivors.
Neither fact asserts a positive survivor. The tiny constructed family has
density decreasing like 2^(-k); all positive-prefix cylinders together
have density Theta(2^(-(1-h)k)k^(-3/2)). Thus the explicit ten-node carrier
samples a much sparser part of the long-growth language. Its prime motif
does not account for the aggregate entropy.

[The square-sum note](crossroads_family_20260926_squares.md) supplies the
parallel lesson in a different arithmetic setting: a selected path lifts
by alternating offsets, but restoring the unused odd cycle in Q_15 forces
all offsets to zero. The exact shared mechanism is the signless-incidence
kernel, not an asserted trajectory map from square sums to Collatz.

| Source and map | Preserved predicate | Lost data / cheapest hostile |
|---|---|---|
| Words to residue cylinders via affine inverse | Exact finite parity and carry | Positive infinite realization; n=1 fixed-horizon exception |
| Allowed b-blocks to one-vertex loop graph | Infinite concatenation stays slope-positive | Completeness of language; finite b yields smaller dimension |
| Residue families to dyadic closures | Every finite dyadic prefix realized | Odd CRT constraints disappear; factor 4669 in orbit note |
| Rational IFS to real attractor | Real contractive limit | Dyadic limit and integer sign; graph (3) |
| Square path to lifted path | Consecutive square sums | All unused edges; the Q_15 odd cycle |

## 5. Reproduction and remaining proof target

Run `python 04-computation/experiments/crossroads_family_20260926_fractal.py`
and the same command with `python -O`. The [raw output](crossroads_family_20260926_fractal.out)
retains exact W_k at selected depths, with decimal dimensions explicitly
labelled approximations. Independent brute words to 16, Spitzer recurrence
to 128, all allowed blocks to 8, and exact two-completion controls accompany
the DP through 512. No floating comparison decides a parity barrier.

Independent audit re-derived all proof steps, checked all 2046 residue
words at lengths 1 through 10 above their exact carry thresholds, and
checked 100 direct rational compositions for (4). The separate paths
agreed. The four-lane reproduction command is
`python -B 04-computation/experiments/crossroads_family_20260926_audit.py`.

The strongest new target is an invariant retaining **one fixed positive
source** while increasing precision, together with carries or height.
Neither occurrence frequency, dimension, prime-support rank, nor a
finite graph description supplies this missing quantifier by itself.
