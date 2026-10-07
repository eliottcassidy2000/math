# Whole-block routing can pass a guard that every individual swap fails

**Status:** PROVED elementary carry/guard rules and one infinite guarded
family; FINITE-EXACT census and ROOT controls. No priority claim. This is
a new tested representation of a composite routing move, not additional
coverage beyond the complete inherited adaptive selector: that selector
already covers the displayed seed and every tested family parameter.
Universal Collatz remains OPEN.

## 1. Inheritance and the changed primitive

The [four-slot compression note, section 5](collatz_four_slot_compression_20261004.md)
already gives the adjacent carry defect and warns that reordering valuations
changes their source cylinders. The
[partitioned-completion note](collatz_partitioned_completion_20261004.md)
already improves on locally decreasing moves by comparing the final child
with the original source. The [root-collision theorem, THM-4555](../../01-canon/theorems/THM-4555-uniform-switches-are-collisions-at-minus-one-trailing-ones-deletion.md)
is a separate uniform trailing-one mechanism.

The source inspiration is the trusted
[CrocSwap integer-multiplication routing work](https://github.com/CrocSwap/integer-mult-bounds):
its old adjacent-routing cost is bypassed using already allowed nonadjacent
field swaps. That paper does not permit swapping Collatz operations.
Here the corresponding question is whether an entire alternative word
can obey the integer guard even when no individual transposition can.

The board is **word order / exact carry / ternary divisibility / source
cylinder / original-source inequality / baseline coverage**. The canonical
hostile is accepting equal slope and total cost as permission to reorder.
The least-used coordinate is the total carry defect of a whole replacement.

## 2. Exact nonadjacent defect and its divisibility price

Apply a word w=(a1,...,aL) from left to right, with each ai>=1. Write

\[
F_w(n)=\frac{3^L n+B_w}{2^{A_w}},\qquad A_w=\sum_i a_i.
\]

For a block w=u,a,v,b,t, put P=cost(u), C=cost(v), and let B_v be the
middle word's carry (zero for the empty word). If w' swaps a and b, then

\[
B_w-B_{w'}=3^{|t|}2^P(2^a-2^b)(3B_v+2^C).                 \tag{1}
\]

To prove this, the carry of (a,v,b) is
3^(|v|+1)+3 B_v 2^a+2^(a+C); subtraction factors it as claimed.
Prefixing multiplies the defect by 2^P, and suffixing by 3^|t|.

Using zero-based positions i<j for the swapped letters, the source shift
that would give the same endpoint is exactly

\[
m-n=\frac{2^P(2^a-2^b)(3B_v+2^C)}{3^{j+1}}.             \tag{2}
\]

Since 3B_v+2^C is coprime to 3, this is an integer if and only if

\[
a-b\equiv0\pmod{2\cdot3^j}.                            \tag{3}
\]

Indeed, for unequal exponents the 3-adic valuation of their power-of-two
difference is zero when their difference is odd, and
1+v3((a-b)/2) when it is even, by the elementary lifting identity for 4.
Equal exponents give the identity move and also satisfy (3).
Thus a direct nonadjacent swap is not automatically cheaper arithmetically:
the later position determines its required ternary precision.

In particular, among words on {1,2,3}, **no nontrivial single transposition
has an integral source shift**. This is an all-length assertion, not an
extrapolation from the census.

## 3. A legal whole-block replacement

For two words with the same multiset, length L, and cost A, set

\[
d=(B_v-B_w)/3^L.
\]

If d is a positive integer, every positive source n of w with n>d has
the smaller positive source m=n-d obeying v and sharing its endpoint.
The shift is even: every nonempty word carry is odd, so the carry
difference is even and division by an odd power of 3 retains evenness.

For completeness the exact source cylinder of w is

\[
n\equiv(2^A-B_w)(3^L)^{-1}\pmod{2^{A+1}}.               \tag{4}
\]

This is an iff: final integrality and oddness in (4) successively imply
integrality and the exact valuation at every preceding step, by reducing
at each initial power of 2. Each positive actual source stays positive.
For the other word, F_v(m)=F_w(n) is the same positive odd endpoint and
m is a positive odd integer; the same argument proves its complete guard.
The program also replays both words independently.

A short carry collision over words on {1,2,3} at length four is
(2,3,1,3) versus (3,3,2,1): their carries are 223 and 547, with
547-223=4*3^4. It already has an immediate descent and is only a control.

The more useful all-growing source word is

\[
\begin{aligned}
w&=(1,1,1,1,1,2,2,2,3,1,1,3),\\
v&=(1,2,1,3,1,3,1,2,1,1,2,1).
\end{aligned}
\]

They have L=12, A=19, B_w=923953, B_v=3049717 and
B_v-B_w=4*3^12. Equations (2)-(4) therefore give

\[
\begin{aligned}
n_t&=257727+1048576t,\\
m_t&=257723+1048576t=n_t-4,\\
F_w(n_t)&=F_v(m_t)=261245+1062882t,\qquad t\ge0.
\end{aligned}                                             \tag{5}
\]

Every nonempty prefix u of w satisfies 2^cost(u)<=3^|u|. Since its
carry is positive, every one of the twelve actual source iterates is
strictly larger than n_t, for every parameter. Thus this dependency works
where direct descent within twelve odd steps does not. For the seed,
ordinary first descent occurs at step thirteen, while m0 first descends
at step four.

No single transposition connecting words over this multiset can pass the
integer guard. The entire replacement does: several nonintegral defects
cancel before the guard is checked. This is analogous to choosing a larger
permitted primitive, not to assuming the constituent operations commute.

## 4. What can actually be reused

A supplied first-hit ROOT word for m_t begins with v: its endpoint in
(5) exceeds one, so it cannot have reached the odd-map fixed point earlier.
Remove that prefix and prepend w. This produces n_t's literal first-hit
ROOT certificate, with the source and both exact guards checked. The API
does not search for a child proof; the independent finite controls do so
only to test the splice.

Appending a common suffix to both words preserves d and the common future.
Prepending a common prefix u multiplies the source shift by
2^cost(u)/3^|u|. Thus the integer translation is preserved only if the
extra ternary divisibility is supplied. For the shift four in (5), no
nonempty common prefix has this property. This precisely limits recursive
reuse: right continuation is free, left insertion costs a guard.

The inherited [adaptive selector](adaptive_boundary_selector_20261004.md)
already reduces the seed257727 via a sibling child206415. Its budget-eight
run also reduces every one of the independently tested t=0..255.
Consequently (5) adds to a **direct-twelve-step baseline**, not to the full
existing adaptive bank. No untested parameter is asserted covered or
uncovered by that stronger bank. The companion
[uncovered join routes](collatz_uncovered_join_routes_20261007.md)
owns the stronger-baseline improvement search.

## 5. Exact controls and finite minimality

The [script](../../04-computation/experiments/collatz_block_permutation_guards_20261007.py)
and [output](collatz_block_permutation_guards_20261007.out) declare the
complete universe. Run:

    python -B 04-computation/experiments/collatz_block_permutation_guards_20261007.py
    python -B -O 04-computation/experiments/collatz_block_permutation_guards_20261007.py

The census includes all797160 words of lengths1..12 over {1,2,3}, grouping
by the multiset and B modulo3^L. Among words satisfying the all-prefix
slope condition, none at lengths1..11 has a higher-carry partner in its
group; at length12 exactly one source word does, namely w above. The
algorithm retains the maximal-carry partner, which is sufficient for
existence of a smaller translated source. This is minimality only in the
declared finite alphabet and prefix-slope criterion.

Further controls cover35400 transpositions with lengths1..5 and letters1..5;
large differences that do satisfy(3); all4096 odd sources below8192 for
an independent cylinder/replay iff; 128 actual family parameters and
independent ROOT splices; and12 malformed/source/mismatch rejections.
The normal and optimized runs agree on79311 checks. The elementary proofs,
not these counts, establish the infinite family and swap criterion.
