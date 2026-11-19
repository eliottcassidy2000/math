# Recurring four- and five-count windows from matched marked profiles

**Proved all heights from explicit finite source data.** No universal
coverage, optimal window size, or continuation of an apex-edge exchange
is claimed.

The previous two-source result 6*2^m±1 used a marked transitive blow-up.
Trying to widen just one copy by post-inflation recontact is not the only
way: change *which source is inflated*, requiring its nonconstant
coefficients to agree exactly. This implements the archive's separation
between complete local signatures and their projection onto one reading.

Use standard lexicographic upper-triangle edge codes, with bit 1 meaning
the smaller label beats the larger. In all sources the mark is vertex 0.
Here are two literal source banks:

| order | source codes | base H in corresponding order | common higher visit profile |
|---|---|---|---|
| 6 | 318,315,335,343 | 11,13,15,17 | c2=6, c3=0 |
| 7 | 1403,1399,11454,5374,1679 | 29,31,33,35,37 | c2=18, c3=c4=0 |

The proof uses the earlier finite marked reader (inherited substitution
kernel THM-3121), H(m)=sum_r c_r r! S(m,r), where c_r counts legal
outside-vertex orders with r separated visits to the marked block.
Here c1 is the displayed base H. The standalone replay literally
enumerates the outside orders and placements to verify every coefficient,
and independently counts paths on all nine inflated tournaments at m=1..6.
The all-height conclusion comes from the coefficients, not extrapolation.

Therefore for **every m>=1** there are actual tournaments for all odds in:

    [6*2^m - 1, 6*2^m + 5]        (order m+5),
    [18*2^m - 7, 18*2^m + 1]      (order m+6).

These are consecutive windows in count space. They need not be connected
by one-apex recontacts; changing to the matching source is the operation.
A target in either scope is decoded without searching for arcs: choose the
row and one of its listed offsets for which (N-offset)/k is a positive
power of two of exponent at least 1; uninflate to the marked source as
the base. For m>1 decrementing m gives genuine strict descent in the same
typed family, via deleting one twin. This does not provide a predecessor
for a target outside the lists.

The gain is repeatable *width* beyond two; the width is still fixed,
whereas contiguous coverage between the dyadic scales requires a width
growing proportionally to the center. The bounded-quotient theorem
in BOUNDED_QUOTIENT_BARRIER.md rules out solving that by merely combining
finitely many fixed uniformly inflated sources.

Run python3 reader_windows_replay.py.
