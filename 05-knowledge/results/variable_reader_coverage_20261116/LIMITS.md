# The fixed window alphabets still cannot cover all odd counts

The two explicit binary-programmed banks give superlogarithmic
absolute coverage, but their combined capacity has a quantitative
ceiling. This concerns these seeds and ordinary/near-sink words
with at least one final near-sink, even allowing arbitrary word
length and any final near-sink length.

With seed coefficient Z after the common synchronization and
word length n, put D=2^n. For t≥1 final near-sinks the exact
program produces

    H = I + Z(2^t D−rev_n(b)) + g(2^t b−J(b)).

Since rev_n(b)≤D−1 and J(b)≤b, one has
H ≥ I + Z((2^t−1)D+1), uniformly in the word.
For fixed t, contributing powers D below a bound Y have
total word count at most 2Y. Consequently each intercept I
contributes at most

    2(X+|I|)/Z · Σ_(t≥1) 1/(2^t−1)

presentations up to X. The series is at most 5/3 because
the first term is 1 and the later terms are at most
(4/3)·2^(−t). Counting distinct values can only reduce this.

Thus the strong three-reader bank (Z=82) and the five-reader
bank (Z=334) together hit at most
(10/82 + 50/(3·334))X + O(1) < 0.173 X counts up to X.
In particular most odd counts cannot be covered by simply running
these two seed banks longer. The positive lower bound Ω(X^(2/3))
and the every-target square-root bracketing theorem do not change
that conclusion. This does not bound the spectrum of other seed
banks, different terminal instructions, or a target-dependent
rewire spanning the resulting gaps.
