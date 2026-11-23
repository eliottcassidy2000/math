# A 3-to-5 arc flip with two Collatz seeds exactly 21 − 7 apart

**FINITE-EXACT.** Label four tournament vertices by valuations \(1,7,3,2\). Use the arcs
\[
 1\to7,\quad 1\to3,\quad 7\to3,\quad 3\to2,\quad 1\to2,\quad 2\to7.
\]
Reversing \(1\to7\) raises the Hamilton path count from **3 to 5**. In the positive odd \(3x-1\) map, two of the affected spanning words have an exact common tail.

| Orientation | Hamilton path of valuation labels | Odd integer trajectory |
|---|---|---|
| Before | \((1,7,3,2)\) | \(541\to811\to19\to7\to5\) |
| After | \((7,1,3,2)\) | \(555\to13\to19\to7\to5\) |

The seed displacement \(555-541=14=21-7\) is forced by the two exchanged labels. Two other Hamilton words survive in both orientations, namely \((1,3,2,7)\) and \((1,2,7,3)\). The gain of two comes from the additional words \((3,2,7,1)\) and \((2,7,1,3)\). Each word represents one odd cylinder modulo \(2^{14}\). The [bucket formula](BUCKETS.md#hamilton-paths-are-equal-mass-residue-classes) also makes this an exact calculation of three versus five arithmetic residue classes.

## How deep can two adjacent valuations be exchanged?

**PROVED.** In a positive valuation word for either odd Collatz sign, transpose distinct adjacent letters \(u,v\). Let \(i\) be the number of preceding letters and \(s\) their sum. The two words admit seeds reaching the **same odd terminal integer** after the whole word if and only if
\[
                          2\cdot3^{i+1}\mid(v-u).
\]
When this holds, every allowable terminal integer for one is allowable for the other at the level of divisibility, and arbitrarily large choices give positive realizations. Their seed displacement in absolute value is
\[
                       2^s\,\frac{|2^v-2^u|}{3^{i+2}}.
\]
It is independent of the common terminal and every letter of the common suffix. One can in addition require the two orbits eventually to hit the same prescribed odd integer prime to three.

Indeed the inverse formula for a word of length \(m\) is \(3^m x=2^A y-\sigma C\). The transposition preserves \(A\). In \(C\) only the term immediately after the changed first letter varies, giving difference up to sign \(3^{m-i-2}2^s(2^u-2^v)\). The two ternary terminal congruences coincide exactly when \(2^{|u-v|}=1\pmod{3^{i+2}}\). The order of 2 is \(2\cdot3^{i+1}\) there, which proves the criterion. Division by \(3^m\) gives the seed formula and the [common basin construction](BUCKETS.md#lift-any-finite-word-into-a-specified-basin) supplies the final anchor. The two different words always have distinct dyadic seed residues, even when their ternary terminal classes coincide.

For \((u,v)=(1,7)\) the terminal classes can coincide with no preceding letter. The displacement is \((2^7-2)/9=14\). This explains the two initial integers in the example. Multiplying both labels by 3 gives \((3,21)\), whose difference is 18. It can be exchanged after **one** arbitrary prefix letter and fails the common-terminal test after two. For any pair of even difference, tripling the labels raises its available prefix depth by one because it multiplies that difference by 3. Multiplying by either 2 or 11 preserves its 3-adic valuation. This distinguishes one of the primes \(\{2,3,11\}\) within the exchange rule; their membership in the odd-square theorem follows from a separate size inequality. The path exchange requires both external incidences of the proposed paths to be legal in their respective tournaments.

Run `python3 -B collatz_buckets.py` for the exact words, both seed trajectories, and checks of the interchange rule against direct terminal congruences.
