# A selective-contact reader breaks the factorial lock at unit scale

This experiment tests a mechanism, not a new size bound. Use the
order-nine and order-ten source pairs from [FACTORIAL_LOCK](
FACTORIAL_LOCK.md). Within each pair the sources have equal H and
P(2), while higher cover moments disagree. No common sequence of
whole-module singleton wraps and ordinal operations can separate
the source counts by two.

Now add one vertex v, with arcs v→U and (B−U)→v for a labelled
subset U of the source. This is *not* a whole-module context.
Using the same odd-cycle packing argument as for the common
readers gives the exact cut formula
[
 H(B+v)=H(B)+2\sum_{s\in U,t\notin U}K_B(s,t),
]
where K_B(s,t) sums H(B−P) over all even-order directed paths P
(in vertices) from s to t. This follows by removing the unique
odd cycle through v from each packing: once v is removed, its
remaining path has even positive order and the named entry/exit
contacts. The other packings live in B and contribute H(B).
The matrix is a genuinely selective observation of the two-port
path counts rather than the unrestricted P(2).

Here are two actual exact unit repairs. "Mask" refers to the
out-neighbors of v in the printed source row labels of FACTORIAL_LOCK.

| source order | mask B | mask B' | contact count per arm | cut weights | resulting H |
|---|---:|---:|---:|---|---|
|9|288|66|2|3331,3332|8931,8933|
|10|320|5|2|14165,14166|37689,37691|

In each case one may relabel B' to send its two-element mask to
the two-element mask of B. Then the same reader has exactly the
same labelled contact arcs to both modules. The difference is a
response to changing the source, not to changing its total or the
number of contacts from the reader. Both equalities can be checked
by literal Hamilton path subset DP, with no appeal to the cut formula.

This establishes a finite physical bypass of the lock. It does
not put a two-contact cut at every absolute target, nor claim an
unbounded common-slope continuation of these four particular
cuts. The next question is an inverse one: choose a physical
two-port matrix and marked two-element set for which its absolute
cut sum is the prescribed (N−H(B))/2.
