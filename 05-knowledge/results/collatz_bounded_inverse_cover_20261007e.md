# A complete bounded inverse-word cover at the original source

**PROVED:** exact native inverse-word guards, positive smaller-child safety,
Mersenne phase pullback, and finite-union/CRT density formulas.
**FINITE-EXACT:** all 953 contracting positive valuation words of lengths
1 through 8, their 252 retained Mersenne-phase labels and the counting results
below. **CONDITIONAL:** a paid smaller dependency still needs its own ROOT
certificate. **OPEN:** covering every parameter and grounding all children.

## 1. Inheritance and the precise extension

The inverse-word algebra and source-owned common-future receipt are inherited
from [uncovered join routes](collatz_uncovered_join_routes_20261007.md), including
its already published F91 map `(64n-73)/81`. The three-generator normal-form
language in [child normal forms](collatz_child_normal_forms_20261007.md) is a
different bounded universe. The immediately preceding
[guard fusion bank](collatz_complement_guard_fusion_20261007e.md) contains the
first-crossing words `1^k,a` for even a, together with the inherited G17 word.
This packet completes a different finite language: **every ordinary positive
valuation word of length at most 8 whose inverse coefficient is below 1**.
It does not claim the inverse construction itself is new.

The hostile is `31 <--(1,2)-- 27`: a smaller inverse child is a paid dependency,
not an independently grounded point. The corrected near miss is counting a
phase without preserving its word, source and terminal obligation. The least
used sidecar is the full positive carry B. The live objects are the original
source, exact word, native guard, counting cell and still-open child.

The actual source chart is fixed throughout:

\[
 E(t)=924745897+2^{32}t,\qquad n(t)=2^{E(t)}-1,\qquad t\ge0. \tag{1}
\]

All density statements below concern the ordinary parameter t, not natural
density of this sparse set of Mersenne integers. No large source is expanded.

## 2. Every contracting inverse word is safe on its positive native class

For a nonempty word `w=(a_1,...,a_L)` of positive valuation letters, write

\[
 F_w(h)=\frac{Ph+B}{Q},\quad P=3^L,\quad Q=2^A,\quad
 A=\sum a_i,\quad B>0.
\]

Assume **Q<P**. For a supplied positive odd endpoint n, its exact inverse
guard and child are

\[
 n\equiv BQ^{-1}\pmod P,\qquad h=\frac{Qn-B}{P}.       \tag{2}
\]

**There is no lost finite positive native head.** Starting from the odd integer
n, apply the formal inverse letters `(2^a x-1)/3` in reverse order. If any
intermediate acquires a nontrivial power of 3 in its reduced denominator,
every later inverse letter retains that denominator and introduces one more
factor of 3: its numerator remains a 3-unit. Final integrality in (2) therefore
forces every intermediate to be integral. Such intermediate integers are odd.
An integral inverse of a positive odd integer is positive, since `2^a x-1>0`.
Thus all stages, including h, are positive odd integers and reverse to the
actual forward valuation word w. Finally Q<P and B>0 give **0<h<n**.

In particular n=1 cannot satisfy (2), since it would require a positive
integer h<1. A first-hit ROOT cannot occur inside the word: an actual ROOT
state cannot subsequently reach n>1. This proves strict native receipt safety
without relying on the first-crossing form of the preceding bank.

The least source is just the least positive odd representative of (2), with
ordinary period 2P. The implementation explicitly retains that least source
and checks the same result through the older general cut compiler for all
953 words. The 252 retained labels have maximum least source 13117; the
exponents in (1) exceed every such bound, so every parameter cut is zero.
There is consequently **no too-small positive native Mersenne hostile** for
this contracting alphabet. Dropping native integrality, or dropping Q<P, is
a genuine change of hypotheses: source 3 fails the `(1,2)` native guard, and
the one-letter word `(2)` includes ROOT padding at source 1 and is excluded.

The receipt is `n --()--> n <--w-- h`. Its original source never changes.
A supplied strict ROOT word for h must begin with w; cutting that prefix
exports the remaining ROOT word for n. The production API does no ROOT search.

## 3. Which words have Mersenne phases?

Substitute n=2^E-1 in (2):

\[
 2^{E+A}\equiv Q+B\pmod {3^L}.                         \tag{3}
\]

The carry recurrence gives `B=2^(A-a_L) (mod 3)`. Thus Q+B is a 3-unit exactly
when the final valuation a_L is even. If it is odd, (3) is impossible for
every E. If a_L is even, 2 generates the units modulo 3^L, so (3) has exactly
one exponent class modulo `2*3^(L-1)`. Reducing (3) modulo 3 shows E is odd.
The inherited exact digit-lifting logarithm computes that class.

Writing this exponent class as e modulo `2*3^(L-1)`, (1) gives the single
parameter cell

\[
 t\equiv\frac{e-924745897}{2}\,(2^{31})^{-1}
       \pmod {3^{L-1}}.                                \tag{4}
\]

This is an iff native guard at the same original source. `symbolic_child_mod`
retains precision `P*modulus` before dividing by P; it never divides a residue
computed at insufficient precision. The map preserves the ordered word and
the child obligation. Reducing to a counting cell loses those choices, so the
full labelled list is retained separately.

## 4. Exact bounded union, inherited overlap and residual

The complete universe is

\[
 1\le L\le8,\quad a_i\ge1,\quad 2^{\sum a_i}<3^L.
\]

The maximum cost is 12. There are **953 words**, of which **701** have odd
terminal valuation and no Mersenne phase. All **252** even-terminal records
remain in `entries()`. Removing contained cells only for counting leaves:

| t residue | Period | One shortest retained word |
|---:|---:|---|
| 1 | 3 | 12 |
| 3 | 27 | 1122 |
| 93 | 243 | 112122 |
| 225 | 243 | 111222 |
| 144 | 729 | 1211222 |
| 171 | 729 | 1121222 |
| 255 | 729 | 1212122 |
| 278 | 729 | 1111124 |
| 525 | 729 | 1122122 |
| 675 | 729 | 1112222 |
| 696 | 729 | 1112114 |
| 708 | 729 | 1111214 |

The twelve cells are disjoint and have total density

\[
 \tau_8=\frac{284}{729}=0.389574759945\ldots.             \tag{5}
\]

Length-eight words add labelled alternatives but no new counting cell beyond
the shorter representatives. A separate full-modulus census gives 852 covered
residues out of 2187 and checks every label against the counting view.

The earlier first-crossing/G17 bank through terminal reset 12 has mass
`131325004/387420489`. Combining its cells with all the present cells gives

\[
 \tau_{\rm combined}=\frac{150929272}{387420489},\qquad
 \tau_{\rm combined}-\tau_{\rm old}=\frac{332}{6561}.
                                                               \tag{6}
\]

This is an exact **5.060204... percentage-point increment in ternary parameter
guard mass**, not new ROOT coverage. The earlier F91 row is explicitly
recognized rather than counted as a new mechanism. The least parameter
outside this combined ternary bank is t=0. **t=23 is also outside every row**;
in fact the whole class `t=23 mod27` is disjoint from the combined ternary
cells. This is only a statement about these finite banks, not an all-depth
obstruction at those parameters.

There is also an ordinary-source boundary to a universal pure-inverse scheme:
source **7 has no smaller positive odd actual predecessor at any word depth**.
Its only possible smaller odd children are 1, 3 and 5; their forward orbits
are contained in `{1,3,5}`. Thus no contracting inverse-word certificate can
pay 7. Nevertheless the forward word `(1,1,2,3,4)` certifies 7 to ROOT. This
demonstrates why complementary constructions are necessary. Source 7 is not
a point of the fixed large-exponent chart (1), and this example supplies no
all-depth obstruction at t=23.

For any independently authenticated finite dyadic guard union of mass beta,
CRT gives the fused guard mass `beta+tau-beta*tau`. Therefore (6) adds only
`(1-beta)*332/6561` to that already fused parameter cover. The same formula
for an infinite binary grammar requires its separate uniform shrinking-tail
proof; finite CRT alone does not justify arbitrary countable limits. Different
labels at a covered point remain different child obligations.

## 5. Exact controls and reproduction

The script independently generates the 953-word universe by cuts between
cost-many ones, verifies the terminal-parity criterion, and compares all
positive native minima with the older source-cut compiler. It checks 756
strict ordinary receipts by forward replay and independent rational inverse
replay, exact parameter lifts and denominator-preserving modular readers,
the full 2187-residue counting universe, inherited overlap and the missing
t=23. A supplied child ROOT control cuts `11 --(1,2,3,4)-->1` to the source
13 word `(3,4)`; a wrong-source ROOT word is rejected. Boolean aliases,
nonintegral parameters, wrong guards and forged records are rejected.

Run from the repository root:

```text
python -B -X utf8 04-computation/experiments/collatz_bounded_inverse_cover_20261007e.py
python -B -O -X utf8 04-computation/experiments/collatz_bounded_inverse_cover_20261007e.py
```

The companion `.out` records the exact check count and rational ledger. No
production procedure discovers a ROOT path, and no individual representative
is promoted to a ROOT theorem for its infinite parameter class.
