# Every single-pivot return on a transitive spine

Here a return means equality of ordinary Hamiltonian counts after changing
the orientation of a single apex arc. A gap is the signed version of the
same equality. This completely decodes the absolute target with a
prescribed gap for **all** one-apex single-pivot changes, not just the
earlier central motif 101→111. The classification supplies an infinite
nonconverse neutral family, but excludes this entire construction class
from full or eventual coverage by a neutral/constant-gap return.

## Exact coordinates

Write a contact word as `A b B`, with 1 meaning the apex beats the
spine vertex. Let x be the binary value of A and y the binary value of
the reverse-complement of B. Let p(t)=t mod 2 and let h(t) be the
one-apex ordinary count for the binary word t (h(0)=1). Leading zeros
or trailing ones are neutral, and x=0 or y=0 is allowed, including an
empty side. The count is
\[
 C_b(x,y)=h(x)+h(y)-1+2xy+(1-b)(x+p(x))+b(y+p(y)).
\]
Indeed, the A-to-B inversion terms from the standard path formula
have a factor 2 because the pivot sits in the separating position.
The inversions that use the pivot give the appropriate one-sided
address and the last/first-bit adjacency correction. No other
path term depends on the pivot. For either x or y zero we can remove
its whole neutral side, including its adjacency, so the same
formula applies. Thus a signed gap 2d occurs exactly when
\[
  \lceil y/2\rceil-\lceil x/2\rceil=d.
\]
Taking signed d (or using converse) handles both changes.

## Exact absolute-target decoder and global limitation

Fix any d and either endpoint b. Enumerate four parity pairs e,f.
In each take y=x+2d+e-f, with x≡e, y≡f modulo 2 and x,y≥0.
When both variables are positive, C_b is strictly increasing in x
in increments of 2. In fact 1≤h(t)≤t+1 by the binary recurrence,
so the two h terms can drop by at most x+y; the leading 2xy
increases by 4(x+y)+8, and the final linear term also increases.
The same proof applies starting at x=0 or y=0. Choose the minimum
permitted x of the required parity, binary search for the first
count ≥N, then accept only equality. An upper bracket x=2|d|+N+4
adjusted to parity suffices. Recover A and B from the binary
representations (empty allowed), and orient the actual graph.
This takes O(log²(N+|d|+2)) arithmetic operations, counting the
O(log(N+|d|+2)) binary recurrence steps per evaluation.

For fixed d each sequence has x,y=O_d(sqrt(N)); it has Θ_d(sqrt(X))
strictly distinct endpoints ≤X as X→∞. Consequently the entire class,
including all pivot positions and all orders, has Θ_d(sqrt(X))
distinct absolute baselines. One cannot repair this by choosing a
more complicated arrangement of the old transitive sources or
moving the one pivot: one must genuinely change the interface or
leave this fixed-gap model.

## Literal neutral switches that are not converse

For each odd t≥3 take x=2t, y=2t−1, and both central bits.
Their gap is zero. Both halves have equal binary width (since
2t is not a power of two), and each has at least two marked
bits in its own reading direction: t is odd and ≥3; 2t−1
cannot be a power of two. The resulting tournaments are strong
(first bit 1, last bit 0) and the apex is intrinsically unique:
there are two triangles through the apex using disjoint choices
of A-ones and B-zeros, so deleting any other vertex leaves a
triangle. Any isomorphism preserves the apex and the total
ordering of the deleted spine. The two words are different,
and a converse would exchange the two equal-width half
addresses, which would force x=y. Thus these are two genuine
nonconverse source classes at the same absolute baseline.
For t=3, the words are 1100010 and 1101010 and both give 73.

The same balanced nodes have different response to an added
reader: append r≥1 terminal zeros. The new reverse-complement
address is y+(2^r−1)2^len(B); its parity is unchanged. Consequently
the new count gap is exactly (2^r−1)2^len(B), not zero. This is
an explicit neutral change with a finite, physically measured
hidden coordinate. It is not closed under an arbitrary reader.

## Provenance and independent check

The contact-word count is inherited; the four-branch inverse and all-pivot
quantifier extend [the motif-specific decoder](https://github.com/eliottcassidy2000/math/blob/2d70eb9cd31437bbfa8cf418c6285f79fc8afaaa/05-knowledge/results/chanpak_11_13_26_ready_barrier_20261114/CUT_BALANCE.md).
Run `python3 replay_pivot.py` for literal tournament path counts.
