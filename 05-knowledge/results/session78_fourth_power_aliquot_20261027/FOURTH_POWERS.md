# Four fourth powers through 50: an exact budget theorem

**Status:** PROVED (elementary arguments below); FINITE-EXACT (independent integer checks in `check_fourth_powers.py`).

**Scope:** Nonnegative integer targets, fourth powers of integer bases, and at most four nonzero summands. Equivalently, use exactly four summands and permit `0^4`. Changing the zero convention changes the answer. This note proves a finite-range structure; it claims no new global Waring-type theorem and no significance from the primality of 13 or 37.

## 1. The observed pattern, completely determined

For a sum at most 50, every fourth-power summand is 0, 1, or 16: the next is `3^4 = 81`, and even powers cannot cancel. If `a` summands equal 1 and `b` equal 16, then

\[
 n=a+16b,\qquad a,b\ge0,\qquad a+b\le4.
\]

Thus the attainable values through 50 are exactly

\[
 \{0,1,2,3,4\}\ \cup\ \{16,17,18,19\}\ \cup\
 \{32,33,34\}\ \cup\ \{48,49\}.
\]

There are **14 nonnegative values**, of which **13 are positive**. The complement is

\[
 [5,15]\ \cup\ [20,31]\ \cup\ [35,47]\ \cup\ \{50\},
\]

with interval lengths **11, 12, 13, 1**, hence **37 missing values**. The complement count is 37 in either `[0,50]` or `[1,50]`, because zero is attainable and its removal changes neither complement. The identity is simply `50 = 13 + 37` for the positive interval.

The final singleton is a cutoff effect. If the endpoint is extended to 64, the final gap is the entire interval `[50,63]`, of length 14, followed by `64 = 4*16`. Through 64 there are 15 attainable nonnegative values and 50 missing values. Through 49 the positive partition is 13 versus 36, and through 51 it is 13 versus 38; the particular 13/37 counts are not cutoff-independent invariants.

### 37 versus 50: two different obstructions

Every integer fourth power is 0 or 1 modulo 16. Four such powers can therefore have only residues 0, 1, 2, 3, or 4. Consequently **37 is impossible globally with four fourth-power summands**, since `37 mod 16 = 5`.

In contrast, **50 has the allowed residue 2**, but its least number of fourth-power summands is five:

\[
 50=16+16+16+1+1.
\]

For targets below 81 the exact minimum is

\[
 L_4(n)=\lfloor n/16\rfloor+(n\bmod16).
\]

Thus `L_4(37)=7` and `L_4(50)=5`. The local residue condition loses the quotient, which records how much of the four-term budget is already spent. The complete finite-range condition is

\[
 n=16q+r,\quad 0\le r<16,\qquad q+r\le4.
\]

It is not valid to use this formula beyond the range without restoring the larger fourth powers: `81=3^4` needs only one term, whereas the restricted formula gives `5+1=6`.

## 2. General theorem: two positive denominations below the third power

Let `p >= 2` be even, let `B=2^p`, let `k >= 0`, and let `0 <= n < 3^p`. Write the Euclidean division

\[
 n=Bq+r,\qquad 0\le r<B.
\]

**Theorem A (budget criterion).** The minimum number of nonzero `p`-th powers of integer bases summing to `n` is `q+r`. Hence `n` is a sum of exactly `k` such powers with zeros permitted if and only if `q+r <= k`.

**Proof.** Every nonzero summand is either 1 or B, since all terms are nonnegative and the next positive power is `3^p > n`. Any count pair `(a,b)` satisfies `a+Bb=n`. Reduction modulo B gives

\[
 a=r+Bt,\qquad b=q-t,\qquad 0\le t\le q,
\]

where `t` is an integer. The number of nonzero terms is

\[
 a+b=q+r+(B-1)t.
\]

Its minimum occurs at `t=0` and equals `q+r`; that pair is always valid. Pad with zeros to reach exactly k terms. This proves both directions. No hypothesis `k<B` is needed for this minimum or existence criterion. □

This argument also applies to arbitrary positive integer p if bases are required to be nonnegative. The evenness hypothesis is what permits arbitrary integer bases without cancellation.

The finite-range reduction must not be confused with a general global congruence: it is not true for every even p that all p-th powers are 0 or 1 modulo `2^p`. For example, `3^6 mod 64 = 25`. The global modulo-16 argument above is specifically about fourth powers.

**Theorem B (all count pairs and collisions).** If `q+r <= k`, all possible count pairs in exactly k slots are indexed by

\[
 0\le t\le T:=\min\left(q,\left\lfloor\frac{k-q-r}{B-1}\right\rfloor\right),
\]

with

\[
 (a,b,c)=(r+Bt,\ q-t,\ k-q-r-(B-1)t),
\]

where c is the number of zero terms. If `q+r>k`, there are no such pairs. In particular, `k<B` forces `t=0`, so the count pair is unique. Conversely, once `k>=B`, the value `n=B` has the two count pairs `(a,b)=(0,1)` and `(B,0)`, both within the budget; thus `k<B` is the exact uniform injectivity boundary for this palette.

The unordered count of representations by nonnegative bases in k slots is `T+1` when attainable. The ordered count is

\[
 [x^n](1+x+x^B)^k
 =\sum_{t=0}^{T}\frac{k!}{a!b!c!}.
\]

Coefficients count representations, not distinct attained integers. If negative integer bases are counted separately, each summand in the last formula is additionally multiplied by `2^(a+b)`; this changes counts, not attainable values.

**Strict positivity.** If exactly k nonzero terms are required, then `a+b=k`, not merely `a+b<=k`. In the same finite range, attainability is equivalent to

\[
 n=k+(B-1)b\quad\text{for some integer }0\le b\le k.
\]

For four positive fourth powers through 50, the values are only **4, 19, 34, 49**. Thus an unqualified phrase “exactly four fourth powers” must state whether zeros are allowed.

## 3. Bands, arithmetic gaps, and counting formulas

First consider the restricted palette `{0,1,B}` without a cutoff. The b-th band of attained values is

\[
 I_b=[Bb,\ Bb+k-b]\cap\mathbb Z,\qquad 0\le b\le k.
\]

If `k<B`, these bands are pairwise disjoint. They may be adjacent when `k=B-1`. The gap after `I_b`, before `I_(b+1)`, has length

\[
 g_b=B-k+b-1,\qquad 0\le b<k.
\]

Indeed, subtract the prior upper endpoint and one from the next lower endpoint. Thus the lengths increase by one because each additional B consumes a slot and shortens the next band by one. For `B=16,k=4`, the lengths are 11, 12, 13, 14.

The map `(a,b) -> a+Bb` sends the integer triangle `a,b>=0, a+b<=k` injectively into the integers when `k<B`. Consequently, in the full palette interval `[0,Bk]`,

\[
 D_{B,k}=\frac{(k+1)(k+2)}2,
\]

and the number of missing values is

\[
 Bk+1-D_{B,k}
 =k(B-k-1)+\frac{k(k-1)}2
 =\sum_{b=0}^{k-1}g_b.
\]

For arbitrary k, use the canonical pair `(r,q)` instead of counting all pairs. The full palette count is

\[
 D_{B,k}=\sum_{j=0}^{k}\min(B,j+1)
 =\begin{cases}
 (k+1)(k+2)/2,&0\le k\le B-1,\\
 B(k+1)-B(B-1)/2,&k\ge B-1.
 \end{cases}
\]

For any integer cutoff `N>=0`, the restricted-palette count through N is

\[
 D_{B,k}(N)=\sum_{q=0}^{\min(k,\lfloor N/B\rfloor)}
 \left(1+\min\{B-1,k-q,N-Bq\}\right).
\]

Every summand counts allowable remainders in its own Euclidean quotient block, so this formula does not double-count even when `k>=B`.

**Range guard.** These palette formulas describe actual p-th powers through N only when `N<3^p`. In particular, full-interval formulas for `[0,Bk]` apply to actual powers if `Bk<3^p`; `k<B` alone does not imply that. At the first new denomination, `n=3^p`, the true one-term representation already invalidates the restricted minimum formula. For fourth powers with k=4, the whole palette interval ends at 64, so the guard is satisfied.

## 4. Exact representation counts through 50

For nonnegative bases in four ordered slots, the counts by attained n are:

| n | 0 | 1 | 2 | 3 | 4 | 16 | 17 | 18 | 19 | 32 | 33 | 34 | 48 | 49 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| ordered count | 1 | 4 | 6 | 4 | 1 | 4 | 12 | 12 | 4 | 6 | 12 | 6 | 4 | 4 |

The total is 80, not 14 or 13. There are `3^4=81` ordered choices from bases `{0,1,2}` and only `(2,2,2,2)` exceeds 50. Excluding the zero target leaves 79 ordered tuples and 13 distinct attained positive values. Counting signs gives 609 ordered integer-base tuples through 50, but the same 14 attained nonnegative values.

## 5. Reproduction and hostile controls

Run from this directory:

```sh
python3 check_fourth_powers.py
python3 -O check_fourth_powers.py
```

The checker uses no `assert` statements, no floating-point arithmetic, and no external dependencies. It compares independent ordered-tuple enumeration, dynamic programming, exact quotient/remainder formulas, factorial coefficients, and direct interval complements. Tests include:

- the exact 50 cutoff, the restored 64 endpoint, and all minima through 80;
- zero-allowed versus strictly-positive conventions;
- ordered/unordered/signed counting distinctions;
- the global modular obstruction for 37 and budget obstruction for 50;
- the first larger fourth power at 81 as a deliberate failure of an unguarded formula;
- general even exponents p=2,4,6, including budgets on both sides of k=B;
- the exact first collision at k=B and all full-palette counts for several B;
- deliberate wrong formulas that must be rejected.

The proof establishes the stated universal formulas. Finite tests independently audit representative parameters and boundaries; they are not substitutes for the proofs.
