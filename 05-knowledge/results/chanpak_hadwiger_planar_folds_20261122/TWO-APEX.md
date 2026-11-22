# Paired cuts for an apex-to-apex flip

Let \(D\) be a tournament on the common deletion of two vertices
\(a,b\). Let \(A_i,B_i\in\{0,1\}\) say whether \(a,b\), respectively,
beat \(i\in D\). A single reversal changes only \(a\to b\) to
\(b\to a\). The following is an exact two-boundary pairing, not a
planar Jacobian criterion or a claim about Hadwiger–Nelson.
The odd-cycle comparison and the nested-contact classification below
are inherited from [Nested contacts](../chanpak_cech_perfect_seams_20261122/NESTED-CONTACTS.md),
recast as a matrix pairing. The adjacent-crossed decoder and the two
explicit nested growth rails are the new constructions.

## An exact signed bilinear contraction

For \(i,j\in D\) define the nonnegative kernel
\[
 W_{ij}=\sum_{\substack{S\subseteq V(D)\\ |S|\text{ odd}}}
              H(D-S)\,P_{ij}(D[S]),
\]
where \(P_{ij}\) counts spanning paths from \(i\) to \(j\) and
\(H(\varnothing)=1\). The signed change on the reversal, divided
by two, is
\[
 \frac{\Delta H}{2}
 =A^{\mathsf T}W({\bf1}-B)-B^{\mathsf T}W({\bf1}-A)
 =(A-B)^{\mathsf T}W{\bf1}
       +A^{\mathsf T}(W^{\mathsf T}-W)B .
\]
Indeed, the [odd-cycle packing formula](https://arxiv.org/abs/2307.05569)
counts every nontrivial directed odd cycle with weight two.
Only packings containing the flipped arc change. The unique
cycle in such a packing containing that arc consists of the
two apices and an odd-size ordered path on \(S\); the untouched
cycles off this support contribute \(H(D-S)\). A singleton path
is allowed, and in particular \(W_{ii}=H(D-i)\).
Thus the gap is affine in the contacts of either apex when the
other apex is pinned; the absolute height is a separate count.

### Inherited boundary: crossing contacts are necessary for any larger unit edge

If \(A\ge B\) coordinatewise, every summand with coefficient
\(W_{ij}\) is nonnegative. Each strict discrepancy at \(j\)
contributes at least one by \(W_{jj}\ge1\). Therefore an increase
of two forces exactly one difference, \(a\to j\to b\); the
dual case \(B\ge A\) is equivalent with orientation exchanged.
The diagonal term must be exactly one, so \(D-j\) is transitive:
any cyclic triangle there gives a positive weight-two term.
All longer changing terms must vanish.

Place \(D-j\) in its transitive order, and write \(t_i\) for the
common tags on the spine. Let \(\ell\) be the first spine index
beaten by \(j\) (the end sentinel if there is none), and \(r\)
the last one beating \(j\) (the start sentinel if there is none).
The absence of a changing path of three spine/deletion vertices
forces \(t_i=0\) for \(i<r\) and \(t_i=1\) for \(i>\ell\):
otherwise respectively \(i\to r\to j\) or \(j\to\ell\to i\)
contributes another positive term. Hence \(r\le\ell+1\).
If \(\ell>r\), the two are adjacent because a position strictly
between them could neither beat nor lose to \(j\). In either case,
apart from the at most two indices between \(r,\ell\) inclusive,
the spine consists of common exterior sources and sinks.
Removing them leaves at most five vertices.

For strong pairs the only surviving options have two active
spine vertices. To see this without a size-five search, in the
ordered case \(r=\ell-1\) the two neighbors of \(j\) must have
tags 1 and 0, respectively, to avoid a source and a sink.
In the crossed case \(r=\ell+1\) they have tags 0 and 1, by the
forced-tag condition. These two small graphs have lower rows
\((14,8,18,20,3)\), respectively \((22,16,10,19,4)\) in labels
\(a,b,j,\text{left},\text{right}\); both reverse \((0,1)\).
Literal counts are 13/15 and 9/11. With fewer than two active
spine vertices there are no strong pairs: direct enumeration of
the at-most-four-vertex cases gives none, and is replayed.
In particular *nested out-neighborhoods cannot support any
larger unit seam*, independently of the ordinary deletion.

## Adjacent opposite contacts on a spine have an exact address inverse

Now let the whole deletion \(D\) be a transitive spine.
Its \(W_{ij}\) is zero for \(i>j\), one for \(i=j\), zero
for \(j=i+1\), and \(2^{j-i-2}\) for \(j\ge i+2\):
choose an interior subset of odd size between the endpoints.
Suppose the tags are equal except at two adjacent positions
\(x,x+1\), where \((A,B)=(1,0),(0,1)\), respectively.
Define \(L\) from the common 1 bits before \(x\), and \(R\)
from the common 0 bits after \(x+1\), using weights
\(1,1,2,4,\ldots\) ordered nearest to the discrepancy first.
The diagonal terms cancel and the mixed discrepancy term
is zero. Subtracting neighboring columns/rows of the kernel
leaves exactly
\[
                       \Delta H=2(R-L).
\]
For example an earlier contact at distance \(d\ge2\) from
\(x\) contributes \(2^{d-2}-2^{d-1}=-2^{d-2}\) to the
half-gap; at distance one the contribution is \(-1\).
On the right these signs reverse, giving precisely the stated
weights. The orientation \(a\to b\) is the lower one
when \(R=L+1\). Swapping the apex labels and taking the
converse reduce every sign convention to this one.

Consequently this entire unit-gap architecture has just
one nonnegative integer \(L\) and two terminal-bit choices:
choose the nearest left tag \(\epsilon\) and nearest right
complemented tag \(\eta\), then the far left tags are the
ordinary binary expansion of \(L-\epsilon\), and the far
right tags are the reverse complement of the expansion of
\(L+1-\eta\). Outside zeroes on the left and ones on the
right are ordinal padding. A missing nearest bit is
equivalent to inserting the appropriate padding peripheral
vertex.

This also gives a complete bounded absolute target test.
Independently for each apex, a left common winning contact
and a right common losing contact are endpoints for odd
cycles using just that apex. If their distances from the
nearby discordances are \(d,e\ge1\), their combined weight
is \(2^{d+e}\), since the interior parity is chosen freely.
This is at least four times the product of the corresponding
weights in \(L,R\). Counting both apices' one-cycle packings
therefore gives \(H\ge 1+8LR=1+8L(L+1)\).
It suffices to test this bounded interval and at most four
contact pairs per integer \(L\).

For each choice count the absolute height by distributing
spine vertices among the three ordered path segments around
the two apices, in each of their two possible orders.
Each segment is automatically increasing along the spine.
A prefix must end at a vertex beating the first apex if
nonempty; the middle must start at a first-apex win and end
at a second-apex loss, and can be empty exactly when the
connecting arc points forward. The suffix must begin at a
second-apex win if nonempty. A streaming 12-state counter
records whether the current prefix end is admissible,
whether the middle is empty or has valid/invalid current
end, and whether the suffix has started; invalid initial
middle/suffix choices are never taken. This gives exact
heights in linearly many arithmetic operations in spine
length. The C++ replay checks this architecture for
\(p=3,\ldots,20\) and finds no alignment with either
perfect endpoint pair. This is not a universal exclusion.

The middle width-preserving comparison is a quadratic
transfer whenever the far tags are reverse-complements:
record for each of the three path slots whether it is empty,
or the first and last common tag bits of its assigned
subsequence. The state on the left determines the right
state sum by swapping and complementing boundary bits.
With the central contacts fixed, concatenate each assigned
slot across the three regions. Acceptance is a fixed
bilinear pairing in the two state vectors. Hence the
height is a fixed quadratic form in the transfer state of
one binary word; the signed edge difference is still the
much smaller paired-cut imbalance above.

## An explicit nested strong family

Take left common tags \(1^k0\), the two adjacent opposite
contacts, and right common tags \(0^{k+1}\) with \(k\ge1\).
Then \(L=2^k-1,\ R=2^k\), so each pair has gap two.
Inflate vertices 0 and 5 in the lower row-mask source
\[
                 (62,252,184,112,224,192,133,9)
\]
by transitive blocks of width \(k\), and flip only \((6,7)\)
for its upper partner. Both are strong: both apices beat
the first spine vertex and lose to the last, and the spine
itself points forwards. The family is nested as marked
induced pairs under enlarging the two blocks.

A path revisits the left block only after stepping in from
an apex, and departs the right block on each nonfinal visit
to an apex. Thus each has at most three runs. In a block of
width \(k\), the numbers of ways to fill one, two or three
ordered nonempty runs are \(1,2^k-2,3^k-3\cdot2^k+3\).
The lower skeleton signature matrix by the two run counts
is
\[
 \begin{pmatrix}259&236&66\\514&592&222\\258&362&162\end{pmatrix}.
\]
It is computed by appending an allowed next source symbol,
using each ordinary vertex once and the two marked symbols
at most thrice. Multiplication by the two filling vectors
and expansion gives the exact lower rail
\[
 162\cdot9^k-388\cdot6^k+298\cdot4^k+128\cdot3^k-250\cdot2^k+53 .
\]
Its upper companion adds two. The source itself gives
259/261 on eight vertices.

Despite the two independent physical block growths, the
balanced rail never straddles a Euclid center. For \(k\ge3\)
all non-leading mixed terms in the midpoint are divisible
by 16, while \(162\cdot9^k+54=54(3\cdot9^k+1)\) has
2-adic valuation exactly 3. Such a perfect midpoint would
have exponent 4, but the already smaller \(k=1\) midpoint
is 260, exceeding \(E_4=120\). The first two midpoints
260 and 4128 also fail, directly: their valuations
respectively prescribe exponents 3 and 6.

## A nested rail with all compatible dyadic phases, but no exact hit

The paired cuts do not impose bounded midpoint valuation.
For an explicit system in the same adjacent-opposite architecture,
take the far left word \(1^k01\), the far right its reverse
complement, and both nearest common tags zero. The two weighted
addresses differ by one because the far contributions agree,
while the nearest right zero contributes an extra one. The
lower source for width \(k=1\) has rows
\[
 (1022,4092,1016,4080,3040,1984,3968,3840,512,3072,2325,293)
\]
with ports 0 and 9 of width \(k\), and reversed pair \((10,11)\).
The paths have at most three runs in each port for the same
two-apex reason. The source's lower signature, with row and
column indexed by run counts \(1,2,3\), is
\[
 \begin{pmatrix}
 13919&20804&9250\\
 28900&45776&21626\\
 16090&26666&13122
 \end{pmatrix}.
\]
Put \(X=2^k,\ Z=3^k\). The exact midpoint (the lower height plus
one) is
\[
 M_k=18998X^2-30440XZ-10988X+P(Z),\quad
 P(Z)=1982+7488Z+13122Z^2 .
\]

**Every dyadic phase divisible by 32 occurs.** For even \(k\)
we have \(Z\equiv1\pmod8\), and
\[
 \frac{P(1+8u)}{32}=706+8433u+26244u^2.
\]
This polynomial is the identity modulo 2 and its derivative
is odd at every integer. In fact the difference between its
values at \(u\) and \(v\) is \((u-v)\) times the odd integer
\(8433+26244(u+v)\). Thus it is a permutation modulo every
power of two. Also \(9\) generates the units congruent to
1 modulo 8, at every dyadic precision: repeated squaring
gives \(v_2(9^{2^j}-1)=j+3\). Hence for every \(s\ge5\) and
every prescribed \(R\equiv0\pmod{32}\), choose \(u\) by the
polynomial permutation, and choose an even \(k\) with
\(3^k\equiv1+8u\pmod{2^s}\). Add a period to \(k\) so that
\(k\ge s\). The \(X\)-terms then vanish modulo \(2^s\) and
\(M_k\equiv R\pmod{2^s}\). This proves in particular that
every exact valuation \(v_2(M_k)\ge5\) occurs. The seam is
literal and strong for all \(k\ge1\); both apices beat the
first spine vertex and lose to the last.

**No member is an exact Euclid seam.** All even \(k\) give
\(M_k\equiv2\pmod3\), since \(X\equiv1,\ Z\equiv0\pmod3\).
Every Euclid center is 0 or 1 modulo 3. For odd \(k\ge5\)
the \(X\)-terms vanish modulo 32 and \(Z\equiv3\pmod8\).
Then \(P(Z)\equiv P(3)\equiv16\pmod{32}\), which forces
the Euclid exponent to be 5. All heights here exceed
\(E_5=496\) by induced monotonicity from \(k=1\).
The two short odd widths \(1,3\) have midpoints \(13920\)
and \(4323024\), neither an Euclid center.
So arbitrary dyadic precision is compatible with a fixed
obstruction at another prime, even on a genuine nested
family of single-arc flips.

Replay the [small and symbolic controls](check_two_apex.py) with Python. Compile [the finite scanner](hn-two-apex-scan.cpp) with C++17 and give it exponents 3 through 20 as arguments; its [saved output](hn-two-apex-scan.out) gives the address counts and exact results.
