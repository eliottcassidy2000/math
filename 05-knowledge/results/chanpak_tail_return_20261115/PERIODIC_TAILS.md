# Period, drift and bit reversal in exact tournament tails

The animal images motivate three levels of the same exact tail function: a literal period, the accumulating change over one period, and a possibly nonperiodic history seen through a moving readout. These are statements about the tournament model, not assertions about animal locomotion. The input recurrence is the [ordinal painted-cut transfer](../chanpak_selective_target_20261114/CUT_STATE.md); the [earlier binary suffix decoder](../chanpak_collatz_exact_addresses_20261114/BINARY_SUFFIX_DECODER.md) uses the same binary correction with a reversed/complemented address. Binary histories as coordinates also occur in [Bernstein–Lagarias's 2-adic shift conjugacy](https://websites.umich.edu/~lagarias/doc/bernstein.pdf); no claim of Collatz convergence is used here.

## A repeated physical block has an exact drift, and its extrema are rigid

For a painted transitive singleton with mark b, write the coordinates as (h,N,S,d). The exact updates are

```
N' = N + (1−b)(S+d),     S' = 2S+bh,     d' = bh.
```

Fix a word w of length l, with bits b_i (0 ≤ i < l), and set

```
A = 2^l,     B = sum b_i 2^(l−1−i),     C = sum_(b_i=0) 2^i.
x* = −B/(A−1).
```

Normalize the start by h; the block acts by x ↦ Ax+B. The fixed-point x* is a formal baseline, not in general a physical state. At step i in that formal orbit the normalized start is −θ_i, where θ_i is the binary fraction 0.b_i b_(i+1)… with the bits repeated cyclically. Use a cyclic incoming bias b_(i−1) at this baseline. Define

```
µ(w) = sum_(b_i=0) [b_(i−1) − θ_i].
```

Then k copies applied to a physical initial state give exactly

```
N_k = N_0 + C(S_0−hx*)(A^k−1)/(A−1) + h k µ(w)
      + [k>0 and b_0=0](d_0−h b_(l−1)).
```

Indeed the normalized start displacement from the formal orbit grows by 2 each step, and the bias matches the cyclic state immediately after the first bit. Sum the resulting geometric response and the k identical baseline contributions. This separates exact phase correction, exponential growth and linear drift. Runs of just 0's or just 1's have zero drift. The positivity of the mixed periodic drift is not claimed new; the following sharp **length-controlled extremal classification** is the refinement.

**Theorem.** At each length l≥2 over cyclic words using both letters:

- Minimum drift occurs precisely at one run of each letter, with lengths ⌊l/2⌋ and ⌈l/2⌉ (allow rotations and exchanges of 0 and 1). The minimum is (2^⌊l/2⌋ + 2^⌈l/2⌉ − 2)/(2^l−1).
- Maximum occurs precisely at the alternating words for even l, and the cyclic words with one same-bit edge for odd l. In the even case it is l/3. For every word, if E is the number of cyclic same-bit edges, then µ ≤ l/3−2E/9.

**Minimum proof.** For a zero run of length n followed by ones of length m, let θ denote the rotating binary fraction at the start of the following ones. The run gives exactly 1−θ+θ 2^(−n) to µ. Since θ starts with a 1 and its first zero is at position m+1, it gives at least (2^(−m)+2^(−n))/2. If there are r≥2 zero runs, pair each with its following one run; by convexity the sum is at least r·2^(−l/(2r)) ≥ 2·2^(−l/4). This is strictly greater than 2^(1−⌊l/2⌋) for l≥4. A single pair of runs of lengths m,n contributes exactly (2^m+2^n−2)/(2^l−1), strictly below 2^(−m)+2^(−n). Within that single-run class the balanced split uniquely minimizes the numerator, up to interchanging m and n.

**Maximum proof.** Write C_j for the number of cyclic pairs starting in 0 and ending in 1 at displacement j. It is half the Hamming distance to the cyclic j-shift, and C_(−j)=C_j, since 01 pairs and 10 pairs have equal totals. Expanding θ_i gives

```
µ = C_1 − sum_(j≥1) 2^(−j−1) C_j.
```

For odd j the final bit differs from the initial bit unless the intervening edges have an odd number of same-bit edges. A union bound gives 2C_j ≥ l−jE; for even j we may use C_j≥0. The odd weights sum to 1/3, and their j-weighted sum to 5/9. As 2C_1=l−E this yields the stated upper bound. With E=0 all inequalities are equalities; for even length every nonalternating word has E>0. At odd length the words with E=1 are all rotations or complements of each other. For those, at even j Hamming distance is at most j by the same union bound; at odd j use C_j≤l/2, sharpening j=1 by C_1=(l−1)/2. Since the even j-weighted sum is 4/9, this gives µ≥l/3−1/2−2/9+1/8=l/3−43/72. Every other odd-length word has E≥3, with µ≤l/3−2/3, strictly lower.

**Exact-return consequence.** On any common clock (take a common multiple of the two period lengths), two fixed periodic tails have constant gap at infinitely many repetitions only if both their exponential coefficients and their linear drifts after multiplying by h agree. Just matching the period length, or the leading slope after normalization, cannot remove a nonzero drift mismatch. The extrema quantify this mismatch independently of finite calibration. Changing the digital word itself or changing the attachment interface is not ruled out.

## Nonperiodic history: an exact distinction at the level of the limit set

Take a fixed evaluation from the [binary suffix decoder](../chanpak_collatz_exact_addresses_20261114/BINARY_SUFFIX_DECODER.md),

```
F(a)=N_0+S a+h p(a),  p(2u+1)=p(u),  p(2u)=p(u)+u+(u mod 2), p(0)=0.
```

Let a_n be the ordinary binary integer of the first n bits of an infinite stream beginning in 1. This is the **address stream**, not the chronological order of the literal marked tail (that address is reversed and complemented to obtain the marks). The decoder's floor estimate gives F(a_n)/a_n = S+h(1−ρ_n)+O(n/2^n), where ρ_n is the real fraction 0.b_n b_(n−1)…b_1. In particular ρ_(n+1)=(b_(n+1)+ρ_n)/2: old bits remain visible through an exponentially decaying memory. The slope accumulation set is finite **if and only if** the address stream is eventually periodic. It is the full interval [S,S+h] **if and only if** every finite binary block occurs infinitely often in the stream.

For a direct proof, take K to be the set of all subsequential limits in the compact binary sequence space of these reversed histories (append zero padding). Its image under real binary evaluation is the accumulation set of ρ_n. Each real number has at most two binary expansions, so a finite real set implies K finite. Deleting the first digit maps K onto K: predecessor histories give one inclusion and successor subsequences after n_j→∞ give surjectivity. On finite K the deletion is a permutation, hence K consists of finite cycles. Take a length r that distinguishes the elements. At all sufficiently large times the first r+1 digits must agree with an element of K, since every cylinder not meeting K occurs only finitely often. Consecutive histories then track the same cycle because deleting the first digit agrees for r digits. Their newest bits are exactly periodic from that point. Conversely an eventually periodic stream has finitely many reversed limit histories. For the full interval, recurrent occurrence of every finite block makes K the full binary space; if a block is absent recurrently, omit its cylinder, whose interval interior contains points with no alternative expansion outside that cylinder. This proves both claims without treating binary endpoint ambiguity as an extra cycle.

Thus a stream containing every block can approach all normalized slopes, while the same finite source continues to obey its exact absolute count and its capacity bounds. Full normalized slope coverage does not imply all odd targets, just as an exact constant gap in a compatible source window does not promise arbitrary repositioning of that window.

Replay: `python replay_periodic_tails.py` gives direct integer transfers and exact rational periodic coefficients for all words through length twelve; the infinite statements follow from the arguments above.
