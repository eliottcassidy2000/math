# Exact exchanges between different sources and different readers

The [singleton-cycle construction](README.md) with source A, a
transitive reader of order k, and s extra ordinal singleton
vertices inside the source has the absolute output

    N_A(s) = (2^k−1) f_A · 2^s + (3−2^k) h_A.

Here f_A is the two-slot cover count and h_A the ordinary
Hamiltonian count of the unpadded source. Both coordinates
multiply across the ordinal components of A; each padded
vertex doubles f_A without changing h_A.

## Exact synchronization criterion with independent widths

Write N_A(s)=a·2^s+b and N_B(t)=c·2^t+e, with a,c>0 and
s,t nonnegative. For a prescribed signed difference Δ from
A to B solve

    c·2^t − a·2^s = Δ+b−e.

If the right side is zero, the odd parts of a,c must agree.
If they do, equality of valuations fixes t−s and gives an
unbounded family, after imposing nonnegative widths. If the
right side is nonzero, the smaller of s,t is at most its
2-adic valuation: otherwise their common power of two would
not divide the right side. Try each smaller width and solve
the other by one exact division and power-of-two test, in
both orientations. This is a complete test for a given pair
of physical sources. It optimizes over both widths rather
than asking whether a congruence can approximate the target.

## Preserve the whole ordinary count while changing both inputs

The following rows are literal outgoing bit masks.

| Mark | h | f | Outgoing masks |
|---|---:|---:|---|
| A₁ | 91 | 2176 | (252,41,10,16,38,12,62,126) |
| A₂ | 59 | 940 | (12,17,2,102,13,23,55) |
| B₀ | 413 | 6392 | (204,81,34,38,13,19,188,62) |

Use A₁→A₂ with a reader of order 2, or B₀ followed by six
transitive singleton vertices with a reader of order 4.
The normalized profiles are (5369,2045440) and (413,409088).
Their two slopes agree because 3·2045440=15·409088, and their
intercepts agree because 5369=13·413. With s further padding
vertices in either source, the strong tournaments have orders
18+s and 19+s, both at the exact count

    6,136,320·2^s − 5369.

At s=0 this is the prime 6,130,951. The two literal full graphs
at that width are checked by subset path DP. The coefficients
are coprime, so this is not a rail forced to be composite by
a common numerical divisor. No infinitude of primes is claimed.

## Make an exact unit repair at the same slope

Use the following four cores.

| Mark | h | f | Outgoing masks |
|---|---:|---:|---|
| C₁ | 409 | 6352 | (204,49,50,70,73,25,38,126) |
| C₂ | 251 | 4520 | (252,105,234,144,70,152,168,18) |
| D₁ | 149 | 3176 | (188,73,10,16,38,14,189,62) |
| D₂ | 53 | 904 | (92,41,10,16,6,93,30) |

Use C₁→C₂ with reader 2 and s padding vertices, or D₁→D₂
with reader 4 and s+1 padding vertices. The slope common to
both is 86,133,120·2^s; the intercepts are −102659 and −102661.
Thus changing source and reader at the matched slope makes a
literal two-path correction at every width. At s=6 the outputs
are the two primes **5,512,417,021** and **5,512,417,019**.
The replay validates the literal cores, exact transfer and
trial-division primality (not full large-order DP for this
second pair).

This closes absolute equalities and unit repairs on specified
unbounded rails. In contrast to a modular carry, the offset and
the slope both match at actual integers. There is no proof yet
that every admissible odd target lies on such a rail. The two
physical examples came from the source-selection investigation;
the independent-width decoder and direct replay are included
here. Run **node replay_aligned_readers.js**.
