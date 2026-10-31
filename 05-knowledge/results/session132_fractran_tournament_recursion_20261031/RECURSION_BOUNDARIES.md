# Recursive structure, output cost and an exact computation clock

**PROVED elementary consequences of the audited programs and counting models.** These clarify which kinds of recursion the session supplies. No arithmetic-circuit lower bound or fraction-minimality result is asserted.

## 1. Why a fixed FRACTRAN program does not imply a short linear recurrence

Let A_n be the number of unlabelled tournaments. There are 2^binom(n,2) labelled tournaments and every isomorphism orbit has at most n! members. Thus

`A_n >= 2^binom(n,2)/n!`.

Since n!<=n^n, its base-two logarithm grows at least as
`n(n-1)/2 - n log2 n`. This exceeds every linear function of n eventually.
A sequence satisfying a fixed finite-order linear recurrence with constant coefficients has at most exponential growth: its state vector evolves under a fixed finite matrix, whose powers are bounded by C R^n for some constants. Consequently A_n is **not C-finite**.

The same applies to merged counts M_n>=A_n/2. By the fold theorem, self-converse counts satisfy

`S_(2m) >= 4^binom(m,2)/m!`,

because an ordered pair of labelled m-tournaments has at most m! simultaneous relabelings. Thus S_n is not C-finite either. Literal-line counts also grow superexponentially.

In contrast the pure-blue sequence has the proved constant-coefficient recurrence

`p_(n+6)=p_(n+4)+p_n`.

The six-register window is therefore a complete linear counting state for that sequence. The ordinary/self-converse counters use bounded but input-dependent cycle-type states or paired objects instead. Finite program length and fixed recurrence order are distinct properties; unrestricted counters provide the additional memory.

## 2. A valuation-output lower bound

For any fixed FRACTRAN list and chosen output prime p, let

`c=max_i max(0,nu_p(numerator_i)-nu_p(denominator_i))`.

If c>0 and a run reaches an output state after t steps, then

`t >= (output valuation - initial valuation)/c`.

This is just telescoping the per-step increment bound. If c=0, that valuation cannot increase at all. It says nothing about a different output encoding.

The actual A000568 program has c=1 at prime3 and initial valuation0, so it needs at least A_n fraction steps. The16-fraction literal-line counter has c=1 at prime2 and initial valuation1, so it needs at least ell_n−1 steps. Its exact 6ell_n+O(n²) runtime has optimal asymptotic order within this fixed-list, bounded-initial-output, exponent-readout model. Direct binary arithmetic is a different model and computes the n<=25 table much more efficiently.

## 3. The pure-blue stream has an exact integer microstep clock

For a positive integerN, define G(N)=sum_p w_p nu_p(N), with weights

```
p:  2  3  5  7 11 13 17 19 23 29 31 37 41 43 47 53 59 61 67 71 73 79 83
w: 12 10  8  6 18 16 14  8  0 -7 -8 -6 -7 -5 -6 -4 -5 -3 -4 -2 -3 -1 -2
```

All other prime weights are zero. Direct dot products with the 22 fraction columns are exactly1, so every legal step satisfies

`G(N_(t+1))=G(N_t)+1`.

The seed has G=108. At checkpointn the exact elapsed time is

`12p_n+10p_(n+1)+8p_(n+2)+6p_(n+3)`
`+18p_(n+4)+16p_(n+5)+8n−108`.

This equals 10686 at n=24 and is checked against every literal checkpoint. The full valuation matrix has column rank22 with a unimodular22-by22 minor, independently confirming that no nonempty multiset of fractions has product1. In particular **no periodic positive-integer orbit exists from any seed**, not just the intended stream input. Arbitrary-seed termination and output interpretation are not classified.

There is also an explicit magnitude envelope. At p=2, w_p=12 log2 p. At every p>=3, any positive w_p is at most 18, while12 log2 p>=12 log2 3>18; zero and negative weights satisfy the inequality as well. Since all integer valuations are nonnegative,

`G(N)<=12 log2 N`.

Hence every trajectory, for as long as it exists, satisfies

`N_t^12 >= 2^(G(N_0)+t)`.

An infinite trajectory therefore grows at least exponentially in this lower-envelope sense. The integers themselves need not increase at every step. This is an exact consequence of the emitted fractions, not a statistical trend, and it explains why a strongly connected finite control graph is not a periodic integer orbit.
