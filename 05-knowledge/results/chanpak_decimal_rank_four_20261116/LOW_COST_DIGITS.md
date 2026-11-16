# The tower crossing, and exact least-wait lifting

Two historical ingredients worth keeping separate:

- The [March 20 computation](https://github.com/eliottcassidy2000/math/blob/ad21b2d13b10a215d0b1c8ff6a39dd76a048fe14/05-knowledge/results/deep_numerology.out) observes 2^10=10³+24=10³+4! with 10=binom(5,2), the number of arc choices in a five-vertex labelled tournament. The literal vertex-number equation `2^binom(n,2)=binom(n,2)^3+(n−1)!` indeed has only n=2,5 among n≥2. Check n=2,…,6; from n≥6 the left side grows by 2^n at the next n, while the two terms on the right grow by ((n+1)/(n−1))³ and n, both strictly smaller. The equality at n=5 is not an inflation rule.
- The [Gray selector and natural recursion note](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/gray_selector_recursion_20261018e.md) separates the full Fibonacci ladder state from its projected resistance coordinate. Here the exact rank-four image, the dyadic line and the integer scalar are likewise three different levels. The cheap path depends on which level one wants to return. [Old decimal windows](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/powers_two_decimal_windows_20261007.md) already keeps the exponent phase rather than forgetting it after a residue match.

## Exact tower-carry sieve

For v≥1, if `2^u=10^v+3·2^w`, then necessarily w=v. Otherwise dividing by the lesser power of 2 leaves an odd factor strictly larger than one. Thus `2^(u−v)=5^v+3`. The order of 2 modulo 5^j is exactly `4·5^(j−1)`: it is four at level one, and `2^4−1=15` has valuation one, so each power-of-five lift raises that valuation by one.

Let a_j be the least nonnegative discrete logarithm of 3 in this cyclic unit group. It starts with a₁=3. To lift, test only the five choices
```
a_(j+1) = a_j + c·4·5^(j−1),      0≤c<5.
```
Exactly one works. The first values are 3, 7, 7, 107, 607, 8107, 45607. This gives a fast *exact exclusion* sieve: if v≥j the required binary exponent u−v is at least a_j, whereas `2^(u−v)=5^v+3<2^(3v)` for v≥2. Thus for j≤v≤floor(a_j/3) a solution is impossible. At j=32 this excludes all v from 32 through **6,133,770,176,775,323,556,869**. Direct checks at smaller v give v=1 or 3, namely 16=10+6 and 1024=1000+24. This does not rule out larger v than the displayed finite cutoff. It illustrates how a simple order lift prunes enormous windows but does not by itself prove an all-scale exponential equation theorem.

## The genuinely least-wait prescription for one literal golden clock

In [the fixed golden return](GOLDEN_DECIMAL.md), the scalar orbit reads H_n=f M^n z, has exact period L_k=60·5^(k−1) modulo 5^k, and the mod-32 output is (25,29,17,1) in phases n mod 4. On the 60 base positions at a given such phase, exactly three give any chosen scalar mod 5. Each of those positions has unit five-adic derivative and therefore exactly one lift through each higher digit.

This proves an exact shortest-wait algorithm *within the one-word clock*. At the desired phase and residue, keep the three base time indices j∈[0,60). Lift each by at most five trials per digit using modular matrix powers for H_j. If the current time is n, the least positive delay is the minimum of the three positive distances from n to those lifted indices around the cycle of length L_k. No chronological simulation of the whole period is needed. This is an actual optimization result for one permitted repeated return, not a claim that no shorter route exists when all 168 words are available. The faster multiword [mod-800 controller](DECIMAL_CONTROLLER.md) illustrates the difference.

Run `python low_cost_digits.py` to verify the tower cutoff and the fixed-clock selector. The test computes a least wait to the 1001 decimal suffix at a compatible mod-32 phase, and checks the selector against direct enumeration at small precisions.

## Why “least resistance” still has an absolute bottleneck

The joint controller uses true physical words and can hold a dyadic choice while correcting the five-adic side. But it returns the **image**, not the absolute integer vector, and staying in the existing return bank forces an absolute word-versus-growth imbalance. The [golden and sparsity result](GOLDEN_DECIMAL.md) quantifies that obstruction even with occasional excursions. For Chan–Pak, the harder form of least resistance is to preserve enough choices at *small height*, or to pay a signed remainder from a growing carrier, not only to minimize time to the right suffix.
