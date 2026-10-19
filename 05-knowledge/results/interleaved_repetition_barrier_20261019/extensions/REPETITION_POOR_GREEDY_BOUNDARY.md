# A fixed-rule first-crossing family with no macroscopic expansion powers

19 October 2026. **PROVED symbolic and guard-realizability hostile.**
Membership in any positive-displacement source set S_d remains **OPEN** for
unbounded members of this family. This note deliberately preserves that
arithmetic distinction.

## 1. What this refutes

A fixed bank alphabet, a fixed accepting grammar, and a short algorithmic
description do not force the consecutive *pure-expansion* repetition used
in the interleaved repetition theorem. This remains false for expansion-first
first-crossing words with unboundedly many episodes and with fully lawful
positive integer lifts.

The construction below does not refute the endpoint theorem. Its purpose
is to identify the arithmetic work that a grammar-based argument still owes.
It also does not refute longer matched-block methods that allow resets
inside a repeated factor. Those are a separate, stronger test.

## 2. An exact greedy construction

Use the chronological Terras bank from the primary note, and put

    v=s_V=5474401089420219382077155933569751763
          /5316911983139663491615228241121378304,
    mu=s_F=3/4.

Exact integer arithmetic gives

    1<v,  v^9<4/3<v^10.                           (1)

For an integer m>=1 construct W_m as follows. Start with cumulative
multiplier S=1. Repeat m times:

1. Append V, and replace S by vS.
2. If mu S>=1, append F and replace S by mu S.

Finally append one additional F. Only this last F is the first crossing.

Indeed after each iteration, 1<=S<1/mu=4/3. This holds inductively: before
an iteration S<4/3, and v<4/3, so multiplication by v can require at most
one reset to return to [1,4/3). Every conditional reset has output at least
one by its guard. The final reset has output below one. Expansion prefixes
inside V also have multiplier greater than one relative to V's input, so
there is no hidden crossing at a V prefix relevant to this bank stopping
rule. The word is expansion-first and all its earlier macroblock boundaries
have cumulative multiplier at least one.

Let theta=log(v)/log(4/3). The number of conditional resets after m iterations
is exactly floor(m theta). In fact their combined multiplier is

    S_m=v^m mu^floor(m theta) in [1,4/3).

Theta is irrational: a relation v^a=(4/3)^b with positive integers a,b
would give 3^(77a+b)=2^(122a+2b), contrary to unique prime factorization.
There are therefore floor(m theta)+1 explicit resets, tending to infinity.
This logarithm formula is explanatory; construction and all its decisions
use rational comparisons only.

Between successive conditional resets, V occurs either nine or ten times.
To prove this, immediately after such a reset the multiplier lies in [1,v):
the preceding state before its last V was below 4/3. From that interval,
eight further V's still cannot reach 4/3, while ten always do, by (1).
The initial run has length ten, and the final run may be truncated. A final
extra F can also occur directly after a conditional F; that creates no
expansion repetition.

## 3. Repetition and description bounds

In the displayed V,F parse every consecutive pure-expansion factor is a
power of V contained in one run of length at most ten. Thus every factor
U^k with U over {K,V} and k>=2 has duplicated Terras length

    r=(k-1)|U|<=9*122=1098,

while total expansion length is L=122m. For every m>900,

    100r<L

for every such factor. The primary one-percent repetition exclusion does
not apply to W_m in this parse.

A parse-robust version holds for m>915. Splitting V into K^4 F introduces
an internal reset, so each new expansion run has form V^q K^4, V^q, or K^4.
A power of a mixed K,V word cannot fit inside V^q K^4: it would require a K
before a later V. Thus duplicated length stays at most 1098, while total
expansion length in any such parse is at least 120m. For m>915 the strict
inequality 100r<L therefore persists. Conversely, K^4F=V cannot absorb an
explicit greedy F, because that F is immediately preceded by the final F
of a V or by another explicit F. These bounds cover all parses obtained
by this exact substitution; the m>900 threshold above refers to the
displayed V,F parse.

The entire family is specified by one fixed rational algorithm and the
integer m. A conventional binary description of a member uses O(log m)
bits in addition to the fixed program. No uncomputable Kolmogorov-complexity
claim is needed. All words also belong to the fixed regular bank language
V(V|F)*F. Neither fact implies a large pure-expansion power.

This statement concerns a fixed grammar that *accepts* the words, not a
uniformly bounded straight-line grammar generating each particular word.
Pumping an accepting automaton changes the word; it need not repeat a block
already present at a useful scale, and it does not preserve the source or
endpoint equation.

## 4. Lawful cylinders are not diagonal endpoint witnesses

Write the full branch of W_m as

    T^N(x)=(a x+b)/M,  M=2^N, a=3^h<M.

For every m, the residue x=-b a^(-1) modulo M has infinitely many positive
integer lifts, each following the entire displayed chronological word.
Thus this is not a fake parity language or a guard-inconsistent hostile.

For a fixed displacement d>=0, however, endpoint realization requires

    Delta=M-a>0,
    x=(b-Md)/Delta in the positive integers.      (2)

The exact missing test is equivalently

    b>Md and Delta divides b-Md.                 (3)

If (3) holds, (2) automatically gives a lawful source: a x+b=M(x+d), and
the standard odd-coefficient cylinder induction authenticates every prefix.
There is no additional parity assumption to hide in this equivalence.

An exact carry-residue formulation is useful. Since gcd(M,Delta)=1, set

    delta = b M^(-1) mod Delta,  0<=delta<Delta.

Then the possible nonnegative endpoint displacements are precisely

    d=delta+j Delta,  j>=0,  d<b/M.               (4)

In particular a non-descending endpoint exists on the lawful cylinder if
and only if delta<b/M. This is an exact scalar arithmetic filter, but its
modulus Delta is unbounded; a finite grammar does not carry it for free.

The step-aware carry lemma gives b/M<N/6. Consequently any admissible
residue in (4) must lie in a short initial interval [0,N/6) of the modulus.
The primary short-range lemma already excludes all such hits when N<=2^52.
Beyond that range, this family isolates a precise next target:

> Prove delta_m>=b_m/M_m eventually, or prove that for each fixed d the
> congruence and positivity in (3) cannot persist along an unbounded
> subsequence relevant to the proposed 2-adic boundary.

This is a statement about actual endpoint realizability. Merely building
W_m, verifying its cylinder, or knowing its compact description proves
none of it. A proof that this particular family is excluded would still
need an argument covering other repetition-poor families before resolving
the full S_d boundary.

## 5. Small exact controls

The companion checker verifies the two exact power inequalities, constructs
selected words m=1,10,901,1024 by integer comparisons, checks every crossing
boundary, reconstructs literal affine coefficients independently, and tests
the carry-residue endpoint filter. These are controls of the construction
and equivalence, not an enlarged episode-cap theorem. The all-m statements
above are proved by the interval invariant and the nine/ten-run argument.

## 6. The stronger matched-factor test still detects this family

There is a provable literal factor-complexity bound here, not just a short
program. Let f_n be the indicator of the conditional reset after the nth V.
The exact multiplier invariant gives

    f_n=floor(n theta)-floor((n-1)theta).

A length-k factor of this binary reset sequence is determined by the phase
z={(n-1)theta}. Its coordinates can change only at the k points
{-j theta}, 1<=j<=k. These points divide [0,1) into k+1 half-open intervals;
the floor convention assigns each endpoint to its right interval. Therefore
the infinite reset sequence has at most k+1 distinct length-k factors. This
is a direct rotation-partition proof; no classification theorem about
Sturmian words is needed.

Encode its symbols by

    0 -> V,  1 -> VF,

whose literal lengths are 122 and 124. The resulting infinite parity word
G is exactly the greedy stream before the extra final F. To describe an
r-bit factor of G, specify a length-k reset factor and an offset within its
first image, where

    k=ceil(r/122)+2.

There are at most k+1 reset factors and at most 124 offsets. The k images
have length at least 122k>=r+244, more than enough to cover the offset and
r desired bits. Consequently

    p_G(r)<=124(ceil(r/122)+3)
           <=(62/61)r+496.

Appending the final two bits F to a prefix of G introduces at most two new
length-r factors. (If that prefix is shorter than r, the whole word has at
most two such factors.) Hence for W_m,

    p_Wm(r)<=(62/61)r+498.                        (5)

The matched-block extension applies at r=ceil(122m/100). For every m>=5,

    (123/61)r+498<N+1,

because r<=61m/50+1, N>=122m, and
123m/50+123/61+498<122m+1 already at m=5 and increasingly thereafter.

Thus every non-descending lawful source of W_m, m>=5, would enter an actual
periodic Collatz orbit by time at most p_Wm(r)-1, with a period at most
p_Wm(r). In particular this entry time and period are each bounded above
by (62/61)r+498. The primary short-range theorem excludes the finite small
length range outright; the new content is an all-length, unbounded-episode
cycle obstruction derived from literal factor complexity.

This does not silently resolve the remaining endpoint question: an exact
internal cycle has a contracting intervening multiplier and is not ruled
out by the height contradiction. If a long matched factor is followed by
two different observed continuations within W_m, the return alternative
is eliminated and strict endpoint descent follows. Proving such branching
uniformly, or otherwise excluding the resulting actual cycles, is a sharper
next target than asking whether this short program has a long expansion
power. It has no such macroscopic power, yet it has plentiful long literal
matches.
