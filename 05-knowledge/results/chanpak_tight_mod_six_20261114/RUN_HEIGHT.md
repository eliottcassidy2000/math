# Ternary run invariant and an exact near-ceiling decoder

Let P be the binary address of a transitive-spine one-apex word,
p its final bit, and τ the number of adjacent 10 transitions.
The ordinary path count is g(P) from the [balanced reader](TIGHT_COMPLEMENT.md).
Then

    g(P) ≡ 1+τ+[P≡2 mod 3]  (mod 3),
    or equivalently H ≡ 1+τ+[R≡0 mod 3]  (mod 3).

The second form uses the actual root response R=P+1.
This is a run invariant on the same graph as the ordinary count.

To prove it, start with the weighted inversion formula for g−1:
a 10 pair at distance d contributes 2^(d−1) and one extra
when adjacent. Before removing the lower ones, a 1 at exponent
i≥1 contributes 2^i to all its lower positions, and the last
bit contributes zero. Put B=#ones, J=#adjacent 11 pairs and
U=Σ bit_i(−1)^i. Modulo 3 the subtracted weight of a 11 pair
before adjacency is −(−1)^(i+j); summing gives

    g(P) ≡ 1+U−p+(U²−B)/2−J.

If c is the number of one-runs, J=B−c gives B+2J≡c mod 3,
and c−p=τ. Also U≡P and U+U²/2=[U=2] in the three-element
field. This proves the invariant.

## The sharp fixed-run upper height

At normalized width n (first bit 1, last bit 0) and t≥1 drops,
the maximum is attained by

    1^a (01)^(t−1) 0^b,     a+b=n−2t+2,     |a−b|≤1.

An adjacent swap 01→10 strictly increases g: it adds an
adjacent inversion and weakly increases every other affected
inversion distance. Move each surplus one in a noninitial
one-run left across its preceding zero-run, leaving at least
one behind. Similarly move surplus zeros in nonterminal
zero-runs to the right. Each completed transfer retains t,
and repeated transfers give the asserted shape with its two
end lengths variable.

For m=t−1≥1, the alternating middle has addresses P=Q=D
with D=(4^m−1)/3, Q the reverse-complement address, and g
equal to (4^m−1+6m)/9 by the binary recurrence. Splitting its
two end seams gives the total

    (4^m−1+6m)/9 + D(2^a+2^b−2) + 4^m(2^a−1)(2^b−1)+2.

At fixed sum a+b the coefficient of 2^a+2^b is negative,
so balanced lengths maximize it. For t=1 the corresponding
formula is (2^a−1)(2^b−1)+2, also balanced. The maxima strictly
decrease with t: when t≥2 an extremal shape contains 1010,
and switching it to 1100 raises g and reduces t by one.
Thus (H mod 3,R mod 3) selects a drop-count class, and its
smallest allowed t gives a sharp upper height for that class
at fixed width.

## An exact constructive decoder near the top of a width

Given N=2^n−D with 0<D<2^n−1, g(P)≤P+1 and its converse
bound g(P)≤Q+1 imply P,Q≥2^n−D−1. Let
s=ceil(log₂(D+1)) and L=n−s. If L>0 the first L bits are
forced to be ones and the last L to be zeroes. Reject if they
overlap. Otherwise the word has shape 1^L Z 0^L with
ℓ=2s−n free bits (L=0 gives the unconstrained middle).
At deficit D=O(2^(n/2)) only O((D+1)²/2^n) middle words
must be tried, a constant at a fixed square-root height scale.

For L>0 and a nonempty middle, put u=2^L−1. The exact output is

    H=g(Z)+u(P_Z+Q_Z)+u²2^ℓ+1−Z_first+Z_last.

After expansion this is 2^n−α_Z·2^L+β_Z, with
α_Z=2^(ℓ+1)−P_Z−Q_Z and
β_Z=g(Z)+2^ℓ−P_Z−Q_Z+1−Z_first+Z_last.
For the empty middle H=u²+2, giving α=2, β=3.
When L=0 just evaluate g(Z) directly.
At fixed ℓ the residual β is bounded independently of L.
The decoder tests both the exact remainder and the coarse
quotient against the same literal Z. In particular for fixed ℓ
and growing L the residual is already smaller than the modulus,
so equality follows once the correct quotient and residue are
both selected. The run invariant optionally filters by a
rooted mod-3 demand. None of these facts guarantees a successful
middle for every odd N or in every width.

## Symmetric exact decoder near the mandatory lower height

For n≥4 the sharp lower height is

    1 + 2^(n−2) + 2^floor((n−2)/2) + 2^ceil((n−2)/2).

To see it, charge each interior vertex for its mandatory pair
with either the first 1 or the last 0. To minimize its cost,
choose zeros near the first end and ones near the last end,
because the pair weights are nondecreasing in the distance.
The balanced word 1 0^a 1^b 0 attains each minimum simultaneously
and creates no inversions within the interior. The two endpoint
sums are 2^a and 2^b, giving the formula.

For normalized width n≥3 the endpoint inversion contributes
2^(n−2), so put D=N−1−2^(n−2). If D≤0 there is no solution:
each interior vertex contributes at least one more inversion.
For D>0 let L=n−3−floor(log₂ D). If L≥1, every witness has
shape

    1 0^L Z 1^L 0,

or is impossible if the collars overlap. In fact an interior
one at right exponent j costs at least 2^(j−1) against the
last zero; an interior zero costs at least 2^(n−j−2) against
the first one. Values greater than D are forbidden. These
two inequalities give exactly the stated forced end bits.

The free middle width is ℓ=n−2−2L, and the actual count is

    H = 2^(n−2) + (2+P_Z+Q_Z)·2^L + g(Z).

This follows by treating 10^L and 1^L0 as the two outer
blocks, both of count 1+2^L, each with exposed end address
2^L toward the middle. Their direct cross term is 2^(n−2);
they have no adjacent correction into the middle or into each
other. The formula also holds with empty middle and g(empty)=1.

At D=O(2^(n/2)) only a bounded number of literal middles
remain. In any fixed such chart with ℓ<L, g(Z)≤2^ℓ<2^L.
Thus matching the exact quotient and residue at the prescribed
target forces equality, not just a candidate in the right
local class. The upper and lower charts both leave the
unbounded middle region and other tournament families open.

## The two height charts are exact reflections

Let W be any middle word of width ℓ and Z its complement.
The complement identity from the [bilateral reader](../chanpak_absolute_reflection_20261114/COMPLEMENTARY_READER.md)
simplifies the upper residual to β_Z=g(W)+2. Its coarse
coefficient α_Z is exactly q=2+P_W+Q_W, the same coefficient
as on the lower chart. Thus put B=2^(2L+ℓ), M=2^L. The two
physical words give

    1^L complement(W) 0^L     →  B−qM+g(W)+2,
    1 0^L W 1^L 0            →  B+qM+g(W).

They are reflected around B+g(W)+1, with the same bounded
middle count at both ends. If L>ℓ+1 both remainders are
strictly smaller than M, so choosing the coarse quotient and
the residue forces the entire absolute output. This is an
exact positive two-sided family. The constraint that q and
g(W) come from the *same* middle is essential; the reflection
does not choose a middle for every target.

For example L=5 and W=101 give q=9, M=32 and middle count
3. The two graphs have exactly 7909 and 8483 paths, with
remainders 5 and 3, both strictly below 32.

## The prescribed coarse coefficient identifies the asymmetric pairs

At fixed width ℓ, q fixes P_W−rev(P_W). In fact it determines
all the differences between bits in mirror positions, though
it leaves equal pairs 00 or 11 free. Suppose two choices of
mirror differences gave the same q. Subtract the corresponding
anti-palindromic polynomials in the address base, obtaining a
nonzero integer polynomial F with coefficients of magnitude at
most two and F(2)=F(1/2)=0. Then F=(2X²−5X+2)G with G an
integer polynomial by Gauss's lemma. Take the first coefficient
g_j of G of largest absolute value m>0. In F the coefficient
of X^(j+1) has magnitude at least
5m−2(m−1)−2m=m+2≥3, a contradiction.

There is also a direct greedy inverse. Put m=floor(ℓ/2),
k=ℓ−m, and write d=P−rev(P)=2^k A−B, where for outer-to-inner
mirror differences ε_i,

    A=Σ ε_i 2^(m−1−i),      B=Σ ε_i 2^i.

The congruence B≡−d modulo 2^k leaves at most two candidates
in [−(2^m−1),2^m−1]. For each, compute A. The outer difference
ε₀ is zero when B is even. When B is odd it must be the sign of A:
that highest digit dominates all lower digits combined. If A=0
in this odd case reject. Replace A by A−ε₀·2^(m−1), replace
B by (B−ε₀)/2, and shorten the list. A valid end must have
A=B=0. This uses O(ℓ) integer steps for each of the two
initial possibilities; no search among the signed profiles
is needed.

Consequently once q is chosen, at most 2^ceil(ℓ/2) words
remain to test for the small remainder. This does not make
that remainder free: all those choices still use the exact
path-count evaluation.

## What the middle parity does to the carry

With a middle of declared width ℓ, reversal sends P to
(−1)^(ℓ−1)P modulo 3, while Q=2^ℓ−1−rev(P). Thus P+Q≡1
mod 3 for odd ℓ, and P+Q≡2P for even ℓ. At odd ℓ, **both**
coarse coefficients α_Z above and 2+P_Z+Q_Z below vanish
modulo 3. At even ℓ those coefficients determine P_Z modulo
3, and the low count g(Z) then specifies the drop count
modulo 3 (in the upper chart first account for its displayed
boundary correction). This is a physical relation between
binary width, ternary phase, and the exact quotient.

Run **python3 replay_collatz_runs.py** for the checked invariant
and maximizers; the formulas give the all-width conclusions.
