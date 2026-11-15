# A nonzero matched slide cannot carry a full binary address

This is a trade-off between two desiderata in a fixed-core slide:
many ways to place the tape, and nonzero ability to repair the count.
It concerns the exact [finite slide graph](SLIDE_COCYCLE.md), not all
one-apex constructions.

**Theorem.** In any strongly connected component of the slide graph
having two internal outgoing edges at every vertex, \(Q/h\) is a
power of two and \(h\) divides both \(s\) and \(r\) at every vertex.
Every physical return through such a component has zero displacement.
Consequently, at any fixed physical marked core, the number of
length-\(d\) closed addresses through components that admit nonzero
displacement is \(O(\lambda^d)\) for some \(\lambda<2\).

## Binary input and output

For any tournament \(h\) is odd by Rédei, and \(Q\) is even:
interchanging the two ordered cover slots is a fixed-point-free
involution on a nonempty vertex set. In the slide graph's update,
integrality of \(r'\) forces \(w_i\equiv r\pmod2\), independent of
\(z_i\). So there is at most one outgoing edge per input bit \(z_i\),
and its output \(w_i\) is determined by the preceding state. In reverse,
integrality of \(s=(s'-hz_i+Qw_i)/2\) forces \(z_i\equiv s'\pmod2\);
there is at most one predecessor per reverse input bit \(w_i\).

Assume the component has both forward edges internally everywhere.
Every one-sided infinite binary input \(z\) therefore has a literal
output \(w\) staying in the bounded graph. Regard the words as real
binary fractions \(x,y\). The bounded start recurrence gives
\[
 y=(s+h x)/Q
\]
from the initial state \(s\). Moreover, the first \(k\) output bits
are determined by the first \(k-1\) input bits. For every interior
threshold of this affine map at \(y=j/2^k\), its preimage
\(x=(Qj-2^ks)/(2^k h)\) must be a dyadic cut of the input: if it
were strictly inside a relevant input cylinder, the two open sides
would give conflicting output prefixes. This is independent of
the choice of binary expansion at the boundary itself.
For large \(k\) choose two consecutive such \(j\). Subtracting the
dyadic preimages and recalling that \(h\) is odd shows \(h\mid Q\).
Substituting one preimage then shows \(h\mid s\). By the degree bound
in reverse and the edge count \(2|V|\), the component also has both
reverse edges internally everywhere. The same argument in reverse,
interchanging the channels, gives \(h\mid r\).

Let \(m=Q/h\) now be an even integer. Take two infinite input words
differing only in the first bit and with an irrational common tail.
Their output real fractions differ by exactly \(1/(2m)\).
On the other hand, iteration of the end recurrence, for each \(d\),
shows that the two output **least-significant first** integers
(truncated at \(d\) bits) differ by \(m\) modulo \(2^d\).
Adding the ordinary integer \(m\) to that bit stream changes only
finitely many bits: above the digits of \(m\) any continuing carry
stops at the next zero, and the output stream is not eventually all
ones since its real fraction is irrational. Thus the difference
of the real output fractions is dyadic. This forces \(m\) to be a
power of two.

## Power-of-two delay means no count translation

At a physical return with \(m=2^b\) and \(h\mid s,r\), put
\(x=s/h\), \(y=r/h\). These are integers from 0 to \(2^b-1\).
Write \(u\) for the \(b\)-bit word coding \(x\). The unsigned start
matching equation says \(u z=w u\) as literal words including
leading zeros. Reversing this equality and comparing the unsigned
end matching equation gives \(y=\operatorname{rev}_b(u)\).

The block \(C\) now has exactly the ordinary and unsigned signature
of the transitive word \(u\) scaled by \(h\). Substituting \(u\) for
\(C\) in a transitive exterior can change the scaled apex count
only by three local corrections: the internal marked count, its
signed end difference times the last exterior bit on the left,
and its signed start difference times the complement of the first
exterior bit on the right. This follows by expanding the [cut
transfer](../chanpak_selective_target_20261114/CUT_STATE.md) twice.
For \(pCzq\) and \(pwCq\) those two exterior bits are the *same*,
since \(p=\operatorname{last}(w)\), \(q=\operatorname{first}(z)\).
The substituted words are identical, as are all three corrections.
The physical displacement is zero.

In each remaining strongly connected component with a nonzero
physical closed walk, there is a missing internal input branch.
From every state within that component one can reach such a state
in at most \(|V|-1\) edges. Thus at least one binary word of length
\(|V|\) cannot be read internally from any given state. Blocking
languages into those lengths gives the claimed exponential bound
with a base strictly below 2. There are only finitely many states
for the fixed ordinary core, so a common \(\lambda<2\) works across
the relevant components.

## Absolute capacity consequence

Restrict to the guarded, **strong** slide family (left guard 1, right
guard 0) using a fixed marked core and any active closed tape. One
may split that tape at any return position and move the core there.
For a tape of total length \(d\), this produces at most \(d+1\)
variants for each graph address. Their counts are at least
\(1+2^d\) for \(d\ge1\): from the apex go to the left guard, take any
even subset of tape vertices in ordinal order, then the right guard
and the apex. Each is a distinct directed odd cycle whose term in
the [odd-cycle packing identity](https://arxiv.org/abs/2307.05569)
contributes at least two. Ignoring other cycles gives the bound.
Consequently all these tape choices and slide sites together realize
only \(O(X^{\log_2\lambda}\log(2+X))\) counts through \(X\).
This capacity bound is for a fixed (or finite bank of) cores and the
closed-tape address grammar without unrelated exterior placements.

In particular even the unit return does not make the tape a universal
number-line address. The precise obstruction is a **binary delay
versus active carry** dichotomy: a seam that simply delays every input
bit is transparent to the Hamilton count, and a seam that produces
a nonzero displacement must lose some binary addresses. To obtain
both placement and repair, change the marked reservoir with the
target or assign those jobs to different interfaces.

The [exact branching replay](replay_slide_entropy.py) sharpens the
generic estimate for the sampled active components. Among 4096 possible
12-bit inputs, even the most permissive live state allows at most 256
for the earlier 13 profile and at most 32 for the new 45 profile.
For the two specified components the resulting absolute capacity
bounds have exponents (2/3) and (5/12), respectively, up to the
logarithmic factor for slide positions. These exponents are certified
upper bounds, not claims of optimal growth exponents.
