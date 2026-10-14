# Length 21: an eleven-candidate cost–carry certificate

**Status, 2026-10-14.** Cost–carry univalence holds at length 21. Together with the [range theorem through 20](collatz_cost_carry_univalence_20261014.md), this proves it through 21. The new step uses an exact shifted-numerator certificate with only 11 terminal candidates, rather than a full frontier census. The all-length hypothesis remains open; external novelty is not claimed.

## 1. Setting and a general maximal-cost reduction

For a positive integer valuation word `w=(a_1,...,a_L)`, put

\[
A_i=\sum_{j=1}^{i}a_j,\quad A_0=0,\quad
b_i=\lfloor\log_2(3^i)\rfloor,\quad
B_w=\sum_{i=0}^{L-1}3^{L-1-i}2^{A_i}.
\]

The word is **slope-admissible** if `A_i<=b_i` at every prefix. Its affine map and ternary carry are

\[
F_w(n)=\frac{3^Ln+B_w}{2^{A_L}},\qquad
r_w=B_w2^{-A_L}\pmod{3^L}.
\]

The question is whether equal `(L,A_L,r_w)` forces equal words. Exact `(L,A_L,B_w)` already determines the word by the odd-tail decoder described below.

**Maximal-cost normalization, valid at every length.** The numerator `B_w` does not involve the final valuation. Every admissible prefix `p` of length `L-1` and cost `A_p<=b_(L-1)` extends to an admissible word of total cost `b_L` by appending `b_L-A_p>=1`. Consequently the following are equivalent:

1. A same-cost, same-carry collision exists at length `L` at some allowed total cost.
2. Such a collision exists at the single canonical total cost `b_L`.
3. The map
   \[
   \Phi_L(p)=3B_p+2^{A_p}\pmod{3^L}
   \]
   is noninjective on the slope-admissible prefixes of length `L-1`.

For the forward implication, keep the two prefixes and change only their final valuations to reach `b_L`; the congruence between their numerators is unchanged. Conversely, normalize two distinct prefixes to this common cost and multiply their equal numerator residues by the same unit `2^(-b_L)`. This removes cost stratification from the all-length existence question. It does not prove injectivity.

## 2. Exact reduction at length 21

Here

\[
L=21,\quad q=3^{21}=10,460,353,203,\quad b_{21}=33.
\]

The attained numerator range across all admissible words is

\[
m=q-2^{21}=10,458,256,051,
\qquad
M=\sum_{i=0}^{20}3^{20-i}2^{b_i}=54,020,229,503.
\]

Every word starts with `a_1=1`, so all numerators have the same residue modulo 4. Equal carry at equal cost thus requires a numerator difference divisible by `4q`. Because

\[
4q\le M-m<8q,
\]

a distinct pair must have, after ordering its numerators,

\[
B_H-B_\ell=4q=41,841,412,812.
\tag{1}
\]

Therefore

\[
B_H\ge m+4q=52,299,668,863,
\qquad
B_\ell\le M-4q=12,178,816,691.
\tag{2}
\]

The exact fixed-cost ranges permit (1) only at costs 32 or 33. They are equivalent here, since every length-20 prefix has cost at most 31 and the last valuation can be changed by one. We use canonical cost 33.

If two words first differ after a common prefix of cost `S`, and their first unequal valuations are `a,b`, then oddness of the two remaining affine tails gives

\[
v_2(B_H-B_\ell)=S+\min(a,b).
\]

Equation (1) makes this valuation 2. Since both first valuations are 1, the first divergence is exactly at the second letter, between 1 and 2.

The orientation is forced as well. A word beginning `(1,2)` has numerator at least

\[
3^{19}\cdot5+8(3^{19}-2^{19})
=15,105,204,767,
\]

obtained by using ones thereafter before the final valuation. This exceeds the low ceiling in (2). Hence the high word begins `(1,2,...)` and the low word begins `(1,1,...)`.

## 3. An exact branch bound, not a census

Suppose a high-word prefix has length `i<=20`, cumulative cost `A`, and affine numerator `B`. Every completion has numerator at most

\[
U(i,A,B)=3^{21-i}B+3^{20-i}2^A
+\sum_{j=i+1}^{20}3^{20-j}2^{b_j}.
\tag{3}
\]

This bound is attained: jump to cumulative cost `b_(i+1)` at the next letter, follow the subsequent maximal cumulative costs, and finish at total cost 33. All increments remain positive.

Starting with `(1,2)`, the exact search tries each possible next valuation

\[
1\le a\le b_{i+1}-A,
\]

updates `(A,B)` to `(A+a,3B+2^A)`, and discards a branch when (3) is below the high threshold in (2). At depth 20, the numerator is exactly `3B+2^A`; the final valuation merely makes total cost 33.

The resulting certificate has:

- 168 visited prefix states;
- 55 discarded subtrees, each with its explicit attained upper bound;
- 11 retained length-20 prefixes.

The 55 discarded prefixes and 11 terminal prefixes are pairwise prefix-free. A separate cost-only recurrence counts the descendants of every discarded prefix. Their counts, plus 11, equal **2,427,209**, the exact number of all slope-admissible length-20 prefixes beginning `(1,2)`. Thus the small certificate exhausts the possible high branch; none of its discarded subtrees is being silently sampled.

## 4. All eleven shifted partners fail

For each retained high numerator, subtract `4q` and run the unique necessary affine decoder. For a residual word length `s>1`,

\[
B=3^{s-1}+2^aB',\qquad B'\text{ odd and positive}.
\]

Thus `R=B-3^(s-1)` must be positive, the next valuation is forced to be `a=v_2(R)`, and the next residual is `R/2^a`. Every decoded prefix must also satisfy its slope budget.

The table lists every candidate. The JSON supplies its full high word, exact carry at cost 33, decoded low prefix, and exact rejection data.

| High numerator | Shifted low numerator | First obstruction |
|---:|---:|---|
| 52,300,302,719 | 10,458,889,907 | Decoder remainder `-5278` at letter 12 |
| 52,387,642,751 | 10,546,229,939 | Decoder remainder `-166` at letter 16 |
| 52,409,616,767 | 10,568,203,955 | Decoder remainder `-1928` at letter 14 |
| 52,491,405,695 | 10,649,992,883 | Decoder remainder `-1304` at letter 14 |
| 52,569,041,279 | 10,727,628,467 | Decoder remainder `-166` at letter 16 |
| 52,642,734,431 | 10,801,321,619 | Prefix `(1,1,1,1,6)` has cost `10>b_5=7` |
| 52,661,275,007 | 10,819,862,195 | Decoder remainder `-8` at letter 14 |
| 52,730,284,415 | 10,888,871,603 | Decoder remainder `-326` at letter 15 |
| 52,795,789,439 | 10,954,376,627 | Decoder remainder `-166` at letter 16 |
| 52,812,269,951 | 10,970,857,139 | Decoder remainder `-586` at letter 15 |
| 54,020,229,503 | 12,178,816,691 | Decoder remainder `-76` at letter 17 |

Ten candidates cannot be the numerator of any positive-valuation word of the required length: their forced decoder eventually leaves a nonpositive remainder. The eleventh fails slope admissibility at its fifth letter. None is a collision.

> **Certified theorem.** No two slope-admissible words of length 21 have the same total cost and ternary carry. Combined with the preceding range theorem, cost–carry univalence holds through length 21.

## 5. Mechanism and remaining boundary

The mechanism couples an ordinary real bound to a dyadic decoding constraint. The only allowed ternary shift places one numerator very near the upper extremal range and the other near the lower range. Exact upper-completion bounds leave just 11 possible high prefixes; subtracting the forced shift makes each low candidate incompatible with its uniquely determined positive affine digits or prefix budget.

The maximal-cost formulation in Section 1 is a genuine all-length simplification. Further lengths can be investigated by bounded shifts and exact completion envelopes rather than by multiplying a separate census across costs. It supplies no general proof that every such shifted pair fails.

The question concerns a finite language of slope-admissible native words. It is not a convergence proof for one positive integer source, and it is not an identification of all finite-source no-descent words with the slope-admissible language.

## 6. Reproduction

The standard-library checker uses only exact integers, explicit checks that survive Python optimization, and a short memoized cost counter. It outputs the complete 55-subtree pruning certificate and all 11 rejection witnesses. Normal and optimized runs give byte-identical JSON and pass 7,381 checks.

```sh
python3 04-computation/experiments/collatz_length21_shift_certificate_20261014.py
python3 -O 04-computation/experiments/collatz_length21_shift_certificate_20261014.py
```

The published files contain the small certificate, not a materialized multi-million-word frontier.
