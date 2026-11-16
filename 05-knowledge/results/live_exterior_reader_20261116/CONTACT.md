# There is only one potentially small single-bit contact in a reader

The gap construction should not be interpreted as offering a second
independent small contact just by marking another bit in the same word.

For a synchronized binary reader let g=9 and suppose its normalized
slope is 82 <= Z/2^m <= 91. This includes the original reader (m=0)
and its programmable m-bit prefixes followed by synchronization,
when the exponent of the synchronizing copy of 2 is included in m.
Fix the n-bit word width and flip bit i from zero to one. Write
W=2^i, V=2^(n-1-i), and C=2^m V. The exact distance-kernel calculation
as in [gap selection](SWITCH.md) gives

    Delta = 18W - ZV - 9(e(X)+e(U)),
    0 <= e(X) <= W,    0 <= e(U) <= V.

Consequently

    9W - (Z+9)V <= Delta <= 18W - ZV.

If W/C <= 4, then Delta <= -10C. If W/C >= 16, then
Delta >= 9W-100C >= (44/16) W. Since D=2^n and 2^m D=2WC,
both bounds imply |Delta| >= 3 sqrt(2^m D).
A change among the synchronized three seeds may add at most four
to the change; this does not alter the asymptotic lower bound.

The ratio W/C is a power of two. Thus **a single-bit switch with
sub-square-root change must be at W/C=8**, equivalently
2i=n+m+2. There is at most one such site, and at widths with the
wrong parity there are none. The output Hamilton count is itself
of order 2^m D for all these normalized slopes and the inherited
prefixes: the I value from the prefix has magnitude O(2^m),
while the output contribution is between ZD and 2ZD+18D.

This is a constraint on one-bit physical derivatives only. Two large
derivatives can cancel, including with the mixed interaction term;
a different terminal reader need not have the same critical site.
It directs the next design toward multiple-bit commutators, multiple
seam types, or genuinely different readings of the same boundary,
rather than two small single-bit derivatives in one word.

## The tangent walk is physical even between samples

For fixed seed and width, changing an address bit reverses exactly
one arc, between that vertex and the retained root. Along the tangent
walk of [the crossing proof](EXTERIORS.md), the two exterior high
addresses have positive integer increments bounded independently of
N; the low phases stay fixed. Bit reversal just reorders the flips.
In L increments of size a, bit j changes at most La/2^j+1 times,
since each change crosses a multiple of 2^j. Summing over j bounds
the total number of arc reversals by O(L+n), including binary carry
steps. For the mod-9 walk sampled with uniformly bounded gaps, this
remains O(L+n) in the number of sampled steps L. Individual carries
can still require Theta(n) flips. This is amortized locality of the
walk, not a bound on the count jump of each individual arc.
