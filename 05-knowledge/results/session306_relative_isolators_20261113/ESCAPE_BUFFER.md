# A logarithmic buffer for the parity words of a genuine escape

For the ordinary positive shortcut map x→x/2 (even), x→(3x+1)/2 (odd), suppose the orbit never repeats. Then all its values are distinct and it tends to infinity. Choose a time m at which the value y=x_m is the minimum of its entire forward tail; such future-minimum times occur arbitrarily late, with y tending to infinity.

Let S(k) be the number of odd letters among the first k parities after m, and p=log(2)/log(3). The exact accumulated identity is

    log(x_(m+k)/y) = (S(k)−pk) log 3 + ∑_{m≤i<m+k, x_i odd} log(1+1/(3x_i)).

The odd values used in the sum are distinct and at least y. Sorted, their j-th value (indexing at zero) is at least y+2j, whether y is even or odd. Thus, for any k,

    S(k)−pk ≥ −[1/(3y)+(1/6)log(1+2k/y)]/log 3.

Indeed log(1+t)≤t; the sum of the first k terms of 1/(y+2j) is at most 1/y plus the integral from 0 to k of 1/(y+2u). The nonnegative log-ratio at a future minimum permits dropping it. No Diophantine approximation assumption is needed.

Two consequences:

* Escape, if it occurs, gives actual **positive-integer** witnesses to the logarithmically buffered survivor language, for arbitrarily large starting heights y. This is stronger than only an asymptotic slope bound. The error budget is uniform over prefix length, not over an arbitrary collection of local subwords.
* The inherited [THM-4476](../../../01-canon/theorems/THM-4476-thin-divergent-orbits-reciprocal-sums-finite.md) proves that ∑ 1/x_i converges for every nonrepeating integer orbit here. Hence the parity discrepancy from time 0 tends to +∞ by the exact accumulated identity. Start at a global minimum of that discrepancy, or at a minimum over any sufficiently late tail. Every prefix from there is in the *exact* noncontracting survivor language. Consequently, excluding positive integers from the exact language would exclude every divergent positive orbit. Conversely, no point on a positive periodic orbit, trivial or not, can supply a witness in this exact set: on any positive cycle the exact period identity gives a strictly negative affine logarithmic drift, because the positive logarithmic correction must cancel it. A finite initial surplus cannot survive infinitely many repetitions of a negative-drift cycle. Thus, using THM-4476, the **no-divergence half** of ordinary Collatz is equivalent to the exact survivor set containing no positive integer. The cycle half remains entirely separate. In particular, the coefficient stopping-time problem cannot by itself exclude positive nontrivial cycles; compare the distinction in the [September directions note](../collatz_directions_20260926.md).

The distinction mirrors the input paper’s coverage condition: one can first isolate every possible obstruction with a buffered necessary constraint, but has not eliminated any positive counterexample just by calculating the size of that constraint. The first statement uses distinct integers rather than an independence model for the parity bits.
