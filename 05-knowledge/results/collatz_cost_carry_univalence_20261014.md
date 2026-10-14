# Cost–carry univalence through length 20

**Status, 2026-10-14.** A new exact range certificate proves cost–carry univalence for every slope-admissible valuation word of length at most 20. An independent exhaustive frontier confirms the result, including all 13,472,296 words at length 20. The unrestricted all-length hypothesis remains open here. External novelty is not claimed.

This sharpens the finite observation in [Polynomial carry fibres](collatz_polynomial_carry_fibres_20261014.md) and the range argument in the [exact-interface research roadmap](exact_interface_research_roadmap_20261014.md). The polynomial-fibre theorem and its entropy conclusion do not depend on this strengthening.

## 1. Definitions and the question

For an accelerated positive Collatz valuation word

\[
w=(a_1,\ldots,a_L),\qquad A_i=\sum_{j=1}^{i}a_j,\qquad A_0=0,
\]

write

\[
A=A_L,\qquad
B_w=\sum_{i=0}^{L-1}3^{L-1-i}2^{A_i},\qquad
F_w(n)=\frac{3^L n+B_w}{2^A},\qquad
r_w=B_w2^{-A}\pmod{3^L}.
\]

A word is **slope-admissible** when every prefix satisfies

\[
2^{A_i}<3^i,
\quad\text{equivalently}\quad
A_i\le b_i:=\lfloor\log_2(3^i)\rfloor.
\]

All logarithmic floors in the checker are computed by integer bit length. These inequalities are word conditions. They are not being substituted for a theorem about the entire orbit of one fixed positive integer.

The question is whether two slope-admissible words with equal `(L,A,r)` must coincide. For equal total cost, equal carry is equivalent to

\[
B_w-B_v=k3^L\quad(k\in\mathbb Z).
\tag{1}
\]

Exact equality of `(L,A,B)` already determines the word: if `L>1`,

\[
B=3^{L-1}+2^{a_1}B_{\rm tail},\qquad B_{\rm tail}\text{ odd},
\]

so `a_1=v_2(B-3^{L-1})`; divide and recurse, using total cost to recover the final letter. Thus the entire issue is whether a nonzero shift in (1) is possible.

## 2. Exact extremal range

For fixed admissible `L,A`, positivity of the letters and the remaining cost give

\[
i\le A_i\le \min(b_i,A-L+i).
\]

Both bounds can be attained simultaneously. Hence

\[
B_{\min}(L,A)=3^L-2^L,
\qquad
B_{\max}(L,A)=
\sum_{i=0}^{L-1}3^{L-1-i}2^{\min(b_i,A-L+i)}.
\tag{2}
\]

For the lower extremizer, put every surplus valuation in the final letter. For the upper extremizer, use cumulative costs `min(b_i,A-L+i)`. Their successive increments are positive because both candidate sequences increase by at least one, and their last value is `A`.

The largest upper endpoint across all allowed costs occurs at `A=b_L`, giving

\[
U_L=\sum_{i=0}^{L-1}3^{L-1-i}2^{b_i}.
\]

Let `D_L=U_L-(3^L-2^L)`. Then

\[
D_{L+1}=3D_L+2^{b_L}-2^L,
\qquad
\frac{D_{L+1}}{3^{L+1}}-\frac{D_L}{3^L}
=\frac{2^{b_L}-2^L}{3^{L+1}}\ge0.
\tag{3}
\]

The normalized range is nondecreasing. One endpoint check therefore certifies all preceding lengths.

## 3. The extra factor of two: every word begins with 1

Every slope-admissible word has `a_1=1`. For `L>=2`, its affine numerator therefore satisfies

\[
B_w=3^{L-1}+2B_{\rm tail}
\equiv3^{L-1}+2\pmod4.
\tag{4}
\]

Consequently any difference in (1) is divisible by 4. Since `3^L` is odd, `4 | k`, and a nonzero collision requires

\[
|B_w-B_v|\ge4\cdot3^L.
\tag{5}
\]

This improves the parity-only threshold `2*3^L`.

At length 20 the exact values are

\[
\begin{aligned}
B_{\min}&=3,485,735,825,\\
U_{20}&=17,290,915,285,\\
D_{20}&=13,805,179,460,\\
4\cdot3^{20}&=13,947,137,604.
\end{aligned}
\]

The margin is `141,958,144`. Equations (3)–(5) exclude every nonzero shift for all `2<=L<=20`; length 1 has only one admissible word. Exact `(L,A,B)` decoding then proves:

> **Finite univalence theorem.** For every `L<=20`, the map from slope-admissible words to `(A,r_w)` is injective.

At length 21, `D_21=43,561,973,452` exceeds `4*3^21=41,841,412,812`. This says only that this range certificate stops working. It does not produce a collision.

## 4. What any later collision must satisfy

Suppose distinct words share a prefix of length `h` and total cost `S`, and their first unequal letters are `a` and `b`. For a same-total-cost pair, that first difference cannot be the final letter. Splitting the affine numerator at the common prefix and then at the first differing letter gives

\[
B_w-B_v=2^S(2^a B'-2^b B''),
\]

with `B',B''` odd. Therefore

\[
v_2(B_w-B_v)=S+\min(a,b).
\]

Under (1), the exact constraint is

\[
\boxed{v_2(k)=S+\min(a,b).}
\tag{6}
\]

The slope height bound also gives `|k|<=L/3`; the sharper interval (2) is available at each cost. Large common prefixes therefore require collision shifts with large powers of 2, while the real range permits only a short list of small shifts.

For a hypothetical length-21 collision, the exact range permits only `k=+4` or `k=-4`. Equation (6), together with `a_1=1` and `a_1+a_2<=3`, forces the first divergence to be the second letter: the words begin `(1,1,...)` and `(1,2,...)`. Their final letters must have the same parity, since equal carry modulo 3 gives `r_w=(-1)^{a_L}`. These are concrete arithmetic restrictions for the next search or proof, not a claim that such a pair exists.

There is also a useful exact native reformulation. The native source cylinder is

\[
n\equiv(2^A-B_w)3^{-L}\pmod{2^{A+1}}.
\]

If `B_w-B_v=k3^L`, then shifting a source for `w` by `k` produces a source for `v`, and

\[
F_w(n)=F_v(n+k).
\]

Thus a collision would give two nearby native source cylinders that merge after the same number of odd steps and the same total dyadic valuation. Large positive representatives exist in both cylinders. This reformulation retains the native-word requirement; it does not establish an infinite orbit or a cycle.

## 5. Terminal slope alone is insufficient

The prefix hypothesis matters. The words

\[
w=(1,3,1,1,1,1,3),\qquad
v=(2,3,1,1,1,2,1)
\]

have the same length `7`, total cost `11`, and terminal slope `3^7/2^11>1`. Their affine numerators are `4591` and `8965`, whose difference is `2*3^7`, and both carries equal `911 mod 2187`. They are distinct, but they violate prefix admissibility: the first word fails at prefix 2 and the second at prefix 1.

This is an explicit counterexample to the weaker terminal-slope version of univalence.

## 6. Reproducibility and bounded census

The companion C++ checker independently generates every slope-admissible word through length 20 by the exact carry recurrence

\[
r'=(3r+1)2^{-a}\pmod{3^{L+1}},\qquad A'=A+a\le b_{L+1}.
\]

It compares the transition with the direct affine formula for every record, checks separate cost-only composition counts, tests the extremizers and normalized-width recurrence, sorts `(A,r)` keys for collisions, and checks bounded inverse-decoder samples. It retains no giant frontier in the published files. A fixed guard caps each frontier at 14,000,000 records and the census at length 20; two packed 16-byte frontiers require about 310.5 MB at the largest transition.

| Length | Budget | Words | Distinct `(A,r)` pairs |
|---:|---:|---:|---:|
| 17 | 26 | 663,535 | 663,535 |
| 18 | 28 | 1,900,470 | 1,900,470 |
| 19 | 30 | 5,936,673 | 5,936,673 |
| 20 | 31 | 13,472,296 | 13,472,296 |

The JSON includes all earlier levels and the length-21 range diagnostic, with no length-21 census. GNU C++17 builds at `-O0` and `-O2` produce byte-identical JSON. No external library is needed.

```sh
g++ -O2 -std=c++17 -Wall -Wextra \
  04-computation/experiments/collatz_cost_carry_univalence_20261014.cpp \
  -o /tmp/collatz-cost-carry
/tmp/collatz-cost-carry
```

**Remaining question:** does cost–carry univalence hold at every length for this prefix-slope language? The current result proves it through 20 and isolates a bounded dyadic-shift obstruction beyond that point. It gives no convergence conclusion for a fixed positive integer orbit.
