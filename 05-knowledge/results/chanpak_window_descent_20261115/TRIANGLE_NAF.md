# Sparse signed binary decoding of triangle spines

This extends the count-three obstruction in [CRT_WINDOWS](CRT_WINDOWS.md)
to every ordinal spine whose strong components are **only** singletons
and directed triangles. It does not concern spines with other strong
components even when their Hamilton count is a power of three.

In the marked-spine interface, write $h$ for its ordinary Hamilton count,
$Q$ for its two-slot count, and $S,T$ for total weighted partial paths
starting in the marked set and ending outside that set, respectively.
A path has weight equal to the Hamilton count of its complement; the
paths counted by $S,T$ are nonempty. A triangle component has ordinary
count 3 and two-slot count 12, and its three individual start and end
weights are each 3. A singleton has counts 1 and 2. Thus $h=3^r$ for
$r$ triangle components, and $Q/h=2^m$ for **virtual width** $m$ counting
a singleton as one bit and a triangle as two.

Assign a singleton its physical mark bit and a triangle with $j$ marked
vertices the two-bit code $\operatorname{bin}_2(j)$. Let $w$ be their
concatenation from source to sink, including leading zeros. Exact
factorization of $S$ at its start block gives

$$
S/h=\operatorname{val}(w);
$$

indeed, each block's normalized local start weight is its own code
value, multiplied by the two-slot to Hamilton ratios of all blocks to
its right. Similarly $T/h$ is the reverse block word with local code
$1-b$ at singleton $b$, and $\operatorname{bin}_2(3-j)$ at triangle $j$.
Here "reverse block word" reverses the order of components, not the
digits *inside* each code. Therefore

$$
\frac T h-\operatorname{rc}_m(S/h)
=\sum_{\rm mixed\ triangles} 2^{a_i}
       \bigl(\operatorname{rev}_2(j_i)-j_i\bigr),
$$

where $a_i$ is the virtual width of the prefix before triangle $i$ and
$\operatorname{rc}$ is fixed-width reverse-complement.

Each nonzero coefficient is $+1$ or $-1$, and its $a_i$'s have pairwise
distance at least two since triangles occupy disjoint two-bit slots.
This is exactly the signed binary **nonadjacent form**, whose uniqueness
here needs no external theorem. An even discrepancy has zero lowest
digit; for an odd discrepancy its lowest digit is the unique $+1$ or
$-1$ congruent to it modulo 4, since its neighboring digit is zero.
Subtract that digit and divide by 2. The procedure recovers every
position and sign, also for negative discrepancies. Discrepancy zero
has no nonzero digits (its least nonzero digit would survive modulo
the next power of two).

The decoder thus reconstructs every mixed triangle, including its
mark count $j$; rotations identify its marked orientations. Every
remaining bit can be regarded as a singleton. In fact replacing a
monochrome triangle by two monochrome singletons divides its whole
cut-transfer matrix by 3, as shown in [the count-three proof](CRT_WINDOWS.md).
So the compressed marked spine, reconstructed from $(h,Q,S,T)$, has
physical apex count $H_{\rm compressed}$, and the original count is
$3^{r-r_{\rm mixed}}H_{\rm compressed}$. Here $r$ is read from $h$ and
$r_{\rm mixed}$ from the signed binary decoder. This proves:

**In this whole SCC grammar, $(h,Q,S,T)$ determines the local apex
count exactly. No nontrivial matched additive digit exists there.**

The parallel with unique sparse digit expansions is operational: it is
why even *combining arbitrarily many* individually rigid little cycles
does not manufacture a hidden switch. To get a genuine carry source one
must introduce a new strong component or give up these independent
exterior responses. The count-five pair does the former.
