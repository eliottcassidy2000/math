# A 27-containing family with an exact recursive descent certificate

**PROVED elementary family, exponent recursion and rank / FINITE-EXACT
fixed suffix and controls / OPEN universal Collatz coverage.**
Date: 2026-09-27. No novelty claim for inverse-word completion, prime-power
lifting, or the classical Mersenne phase formula.

There is an explicit infinite family containing 27 whose members undergo
an arbitrarily long initial rise and then reach the fixed integer 47.
Its source parameters can be recovered directly from the integer. Along
the actual orbit, one parameter decreases by one every two odd steps.
The construction supplies a genuine recursive certificate on its stated
domain; it does not establish entry into that domain for arbitrary sources.

## 1. Inheritance and the live board

The closest proved mechanism is [inverse completion](arithmetic_braids2_20260917_inverse_completion.md):
a finite prescribed word can be completed to an admissible target.
The [guarded affine-port note](creative_descent_20260925.md) retains the
source cylinder and numerical descent threshold. The current construction
specializes those mechanisms to the repeated negative-cycle word `(1,2)`
and gives an explicit exponent recursion, source decoder, and actual
parameter rank. Its completion theorem is not independent evidence of
global coverage.

The canonical hostile is the fixed-coefficient family `4*8^k-5`, which
contains 27 and has an unbounded initial expanding phase. The corrected
near miss is confusing the eventual success of a deliberately completed
word with success for every source having its first few valuations.
The least-used sidecar is the exponent-index residue clock from
[the ternary-digit/sibling note](ternary_digits_20260925.md), section 3.

The live board is exact source, repeated valuation word, ternary exponent
index, terminal division, recursive rank, and family frequency. The anchor
is a source-decodable descent certificate; the niche is the prime-power
index clock; the wildcard is whether a small recursive certificate can
handle arbitrarily large actual time. The answer is affirmative for this
family, with exact boundaries.

Throughout, `U(n)=oddpart(3n+1)` on positive odd integers, and the first
descent time of n>1 means `min{j>=1: U^j(n)<n}`.

## 2. The expanding two-step block and its retained coordinates

If n=11 mod16, its next two division exponents are exactly 1 and 2, and

```text
U^2(n)=(9n+5)/8,       U^2(n)+5=(9/8)(n+5).          (F1)
```

The inverse block is `(8m-5)/9`, legal for an odd m precisely when
`m=4 mod9`; this inverse is then 11 mod16. These are guards on actual
integers, not freely chosen division words.

For a positive even integer a and k>=1, put `n=a*8^k-5`. Direct induction
gives

```text
U^(2j)(n)   = a*9^j*8^(k-j)-5,      0<=j<=k,
U^(2j+1)(n) = (3a/2)*9^j*8^(k-j)-7, 0<=j<k.         (F2)
```

Before each block, its shifted value is divisible by at least16, so (F1)
applies. Every one of these first 2k iterates is greater than n: the even
boundaries multiply n+5 by 9/8, and each intermediate is greater than its
own preceding boundary. The latter difference is `(n_current+1)/2>0`.

On these blocks, `v2(n+5)` decreases by3 and `v3(n+5)` increases by2.
Consequently `2*v2(n+5)+3*v3(n+5)` and the part of n+5 coprime to6
are conserved. This is an exact phase ledger. It is not a global rank:
the `(1,2)` guard eventually fails, and subsequent operations may create
new precision. For the original family `a=4`, the endpoint is `4*9^k-5`;
its large `v3(endpoint+5)=2k` records the elapsed phase, not a proved
future collapse.

## 3. Choose a terminal target, then solve its exact index guard

For the concrete target 47, define

```text
F(h)=47*4^h+7,                  h>=0,
K(h)=floor((v3(F(h))-1)/2).
```

F(h) is always divisible by3, and it is a positive integer, so K(h) is
a well-defined finite nonnegative integer. For `0<=k<=K(h)` set

```text
a(k,h)=2F(h)/3^(2k+1),
N(k,h)=a(k,h)*8^k-5.
```

The coefficient is a positive even integer. The admissibility condition
can also be written

```text
47*4^h = -7 mod3^(2k+1).                            (F3)
```

**PROVED recursive closure.** For every admissible k>=1,

```text
U^2(N(k,h))=N(k-1,h),
U(N(0,h))=47.                                      (F4)
```

The first assertion is (F1), with h unchanged and k decreased. For the
second,

```text
N(0,h)=(94*4^h-1)/3,
3N(0,h)+1=47*2^(2h+1).
```

Thus the terminal division exponent is exactly `2h+1`, and the endpoint
47 belongs to the actual orbit. No approximation of an exponent or of
the source was used. The whole k-block phase has total division exponent
3k; including its terminal step gives exactly `3k+2h+1` divisions by2
over `2k+1` actual U steps.

The main source family below has k>=1. The k=0 objects are included only
to close the recursion. In particular `N(0,0)=31`, and `31 -> 47` is
not a numerical descent; its role is a base transition to the fixed
verified suffix. No k=0 descent assertion is inferred from the main theorem.

The invariant h supplies the exact arithmetic data, while k is a
well-founded rank at the phase boundaries. A fixed controller can execute
this family with an unbounded k counter. This differs from claiming that
every state represented by that controller has the same rank.

## 4. A ternary digit recursion generates every admissible exponent

The elementary valuation identity

```text
v3(4^d-1)=1+v3(d), d>=1,                            (F5)
```

follows by factoring an exponent coprime to3 as a geometric sum, and noting
that replacing d by3d contributes exactly one more factor of3.
In particular, 4 has order `3^m` modulo `3^(m+1)` and runs through
all residues equal to1 modulo3.

Since `-7/47=1 mod3`, there is a unique residue `h_m mod3^m` satisfying
`F(h_m)=0 mod3^(m+1)`. Use its representative in `[0,3^m)`. It has the
following entirely integer recursion:

```text
h_0=0,
E_m=F(h_m)/3^(m+1),
d_m=(-2E_m) mod3,       d_m in {0,1,2},
h_(m+1)=h_m+d_m*3^m.                               (F6)
```

Indeed `4^(3^m)=1+3^(m+1) mod3^(m+2)`. After the update, the divided
error changes by `47*4^h_m*d_m=2d_m mod3`; the stated choice cancels it.
Large integer powers need not be materialized to find the next digit:
compute F modulo `3^(m+2)` first.

For a given phase length k, **all** the exponents are

```text
h=h_(2k)+t*9^k, t>=0.                              (F7)
```

The first minimal values of `h_(2k)` are

```text
k:       1   2    3     4      5       6       7       8
h_(2k): 0  45  369  6201  52128  465471  465471  465471.
```

The repeated last value is genuine extra divisibility, not evidence that
the exponent sequence stabilizes. Eventual stabilization at an ordinary H
would force the positive integer F(H) to be divisible by every power of3,
which is impossible. The digit recursion has unbounded information depth.

This is the same precise kind of ternary residue clock as the earlier
sibling mechanism: for h>=0 and d>=1,

```text
v3(F(h+d)-F(h))=1+v3(d).                            (F8)
```

Here that clock supplies an integrality guard for the actual orbit family,
and (F4) supplies the separate descending parameter. The clock alone would
not imply descent.

## 5. Exact descent and convergence, including the source 27

The pair k=1,h=0 gives a=4 and `N(1,0)=27`.
Every other member with k>=1 has n>47. For k>=2, the bound a>=2 gives
n>=123. For k=1 and h>=1 the defining formula already gives n>47.

Combining this comparison with (F2) and (F4) proves:

> Every member with k>=1 other than27 has exact first descent time `2k+1`, and
> its first smaller iterate is47. Its first2k actual odd iterates all
> exceed its starting value.

This is an infinite unbounded-depth result at the fixed insufficient-
precision interface `(q,R)=(1,1)`: every source in the family is3 mod8
and has `v2(3(n-2)+1)=2`.

The exceptional source27 has `27 -> 41 -> 31 -> 47`; 47 is larger than
27, so that three-step block is correctly rejected as a descent proof
for27. A directly checked fixed suffix is

```text
47,71,107,161,121,91,137,103,155,233,175,263,395,
593,445,167,251,377,283,425,319,479,719,1079,1619,
2429,911,1367,2051,3077,577,433,325,61,23,35,53,5,1.
```

Successive entries are U-images. Its first entry below27 is23 at index34;
it first reaches1 at index38. Thus27 has first descent at its37th U step,
and every member of the constructed family first reaches1 after exactly
`2k+39` U steps. This universal family claim depends only on the exact
formulas and this finite, replayable suffix, not on general Collatz.
There is no earlier1: every growing-phase state exceeds its source>=27,
every base `N(0,h)` is at least31, and the displayed suffix contains1
only at its final index38.

This also gives an explicit rank on the typed controller states. Assign
`2k+39` to phase boundary `(k,h)`, `2k+38` to its first-step intermediate,
and `38-i` to entry i of the fixed47 suffix. Each actual U transition
decreases this nonnegative integer by one, including the base transition
`N(0,h) -> 47`. The formula is specified by the finite controller and
its k counter, rather than by an unknown stopping-time oracle.

The first two values at k=1 are27 and7,301,195, with respective h=0,9.
The minimal k=2 value is30,647,931,178,357,398,081,453,215,355, at h=45;
its first descent occurs at step5 and ends at47.

In the recursive certificate language of
[the recursive-entry note](entry_20260927_recursive.md), the ordinary members
need a repeated `(1,2)` node and one terminal U step: two macro nodes,
although their actual time `2k+1` is unbounded. The finite27 exception
uses the fixed suffix. A bound on atomic orbit lookahead therefore does
not automatically exclude a recursive certificate grammar.

## 6. Recognize the family from one fixed integer

For h>=1, F(h) is odd, so `v2(a(k,h))=1` and

```text
v2(N(k,h)+5)=3k+1.                                 (F9)
```

For h=0 only k=0,1 are legal; the k>=1 source is27. This yields a
direct membership decoder for the k>=1 family:

1. Accept n=27 with `(k,h)=(1,0)`.
2. Otherwise require n positive odd, `e=v2(n+5)>=4`, and `e=1 mod3`.
   Put `k=(e-1)/3` and `a=(n+5)/8^k`.
3. Require `(3^(2k+1)*a-14)/94` to be an integer power `4^h`.
   The resulting h is the exact exponent; accept this pair.

Every accepted pair satisfies the defining equation and every family
member passes. The power-of-four test uses integer division and binary
valuation, not a floating-point logarithm or an orbit search. This is
a source-computable certificate interface with explicit arithmetic data.

The original one-parameter family `4*8^k-5` has shifted valuation `3k+2`.
Consequently its intersection with this new k>=1 family is **exactly27**.
We have constructed an infinite certified family through27, not proved
the conjecture for all later members of that original fixed-coefficient
family. No finite test is being extrapolated across that distinction.

## 7. Exact frequency of depth, and a sparse source family

Equation (F7) gives a precise recursive frequency statement:

```text
natural density of {h>=0: K(h)>=k} = 9^(-k),
natural density of {h>=0: K(h)=k}  = 8/9^(k+1).       (F10)
```

It concerns the ordinary exponent parameter h. Each additional available
two-step block fixes two further ternary digits. For example, over one
period of length `9^4=6561`, the counts at depths `0,1,2,3,>=4` are exactly
`5832,648,72,8,1`. This is an exact residue count, requiring no independence
or random-orbit assumption.

The associated sources form a zero-density subset of the natural numbers.
Indeed n<=X implies `2*8^k-5<=X`, so k=O(log X). The defining formula
then bounds h=O(log X), uniformly in those k. There are at most
O(log X) possible parameter pairs: with common bounds k<=L and h<=H,
the single residue class in (F7) gives at most
`sum_(k=1)^L(H/9^k+1)<=H/8+L`, and both H,L are O(log X).
This provides no entry theorem
for a fixed source outside the family, and the exponent densities cannot
be relabeled as orbit frequencies.

The [Bernoulli boundary selector](forest_20260926_bernoulli.md), section3,
has a literal use here: for `M=3^(2k+1)` and `d=F(h)`, its endpoint-aware
sawtooth identity detects exactly `M|d`. This retains the divisibility
guard, but adds no descent by itself. The new decreasing coordinate comes
from the actual recursive operation (F4). Primality of any displayed
power-plus-constant number is irrelevant to both mechanisms.
The selector works for this nondyadic M as well: its floor jump at d/M
is one precisely when d is a multiple of M.

## 8. Positive and hostile controls from adjacent power families

A simpler terminal control chooses target1 instead of47. Let

```text
4^h=-14 mod3^(2k+1),
a=(4^h+14)/3^(2k+1),     n=a*8^k-5.
```

The same ternary lifting supplies every h. The first minimal exponents
are `4,76,562,4207,37012,450355`; the first source is75. Its first2k
U iterates rise above the source, and its next iterate is1. Thus the
first descent time and hitting time of1 are both exactly `2k+1`.
This separates the long growing phase from any difficulty in the fixed
terminal suffix.

The Mersenne comparison is an exact hostile to paying exponential growth
with only an exponent-index valuation. For `M_H=2^H-1`,

```text
U^H(M_H)=(3^H-1)/2^e(H),
e(H)=1 for odd H; e(H)=2+v2(H) for even H.            (F11)
```

The first H-1 steps have division exponent1; the final step has exponent
`1+e(H)`. Factoring `3^H-1` proves the displayed valuation: odd factors
of H contribute no extra2, and each subsequent doubling adds one after
the initial `3^2-1=8`.

For every H>=9, `2^e(H)<=4H` and
`3^H-1>4H(2^H-1)`, so this phase endpoint is still above its source.
For the second inequality, `(3/2)^9>36` and the ratio
`(3/2)^H/(4H)` increases for H>2. Direct checks of H=1..8 show that
the endpoint descends exactly at H=2,4,8. No assertion is made about
the later orbit of the other Mersenne inputs.

## 9. Reproduction and next question

Run from the repository root:

```text
python 04-computation/experiments/entry_20260927_families.py
python -O 04-computation/experiments/entry_20260927_families.py
```

The [script](../../04-computation/experiments/entry_20260927_families.py)
and [frozen output](entry_20260927_families.out) check the exponent recursion
through60 ternary digits, independently enumerate exponent periods through
eight digits, count a complete depth period, decode fixed sources, replay
large actual integers for k=1..6 and three exponent-period offsets, and
verify all clocks and source comparisons: **7,263 explicit checks**, with
identical normal/optimized output. The target47 controls reach3,056,710-bit
inputs, and the target1 control reaches900,708-bit inputs.
Mersenne valuations are checked through H=1024 and
their actual U clocks through H=128. Explicit checks remain active under
optimized Python. Universal assertions above have direct proofs; these
ranges do not supply their quantifiers.

The constructive next question is whether an unresolved source can be
rewritten into a certified target family through an actual guarded orbit
segment with a decreasing recursive parameter. Simply completing an
arbitrarily chosen word to47 changes the source and cannot answer that
entry question. This is the boundary between an exact successful family
and a global Collatz proof.
