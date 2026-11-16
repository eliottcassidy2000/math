# Birth fibres and a binary program for continued growth

**Main results.** The [birth atlas](births.tsv) lists **all** smallest
isomorphism classes for every count first attained through order eight:
320 distinct counts and 5,827 tied birth classes. It distinguishes
the order where a count is born from the order where a specified
growth rate is born. Independently, there is an exact binary
controller for the future rate of a marked tournament. It lands
on several prescribed earlier centers, including 30090, without
sacrificing the continuing +2 repair. It does not give universal
odd-count coverage or arbitrary simultaneous center/rate selection.

## 1. Every attainable count has a finite exact birth fibre

Write \(m(h)\) for the least order of a tournament with \(H=h\);
use \(+\infty\) for an empty fibre. Use nonempty tournaments for
the birth question, so \(m(1)=1\). The known
[exclusions 7 and 21](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md)
give empty fibres. We are not assuming that there are no other
empty fibres.

The strongly connected components of a tournament are linearly
ordered, and every full path traverses them in that order.
Their counts multiply; singleton components contribute 1.
Thus a birth graph of count \(h>1\) has no singleton component.
If \(m_{\rm str}\) is minimal strong order, then
\[
m(h)=\min_{h=\prod_i h_i,\ h_i>1}\sum_i m_{\rm str}(h_i).
\]
The complete fibre consists of all ordered component words whose
sum ties this minimum, with all choices of strong birth class in
each slot. Indeed replacing any component by a smaller graph
of its count would contradict birth, so each component is also
globally minimal for its own count. The factor order is part of
the graph: scalar prime factorization alone does not determine a
unique tournament.

A strong component of order \(n\ge4\) has at least
\(5^{(n-1)/3}\ge3^{n/3}\) paths by
[Busch's floor](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3);
at order three the count is 3. Consequently every realized count \(h>1\) has a singleton-free
witness of order at
most \(\lfloor3\log_3 h\rfloor\). Enumerating at these bounded
orders decides membership and gives all birth ties for any
*prescribed* count. This is a finite recipe, not a practical
large-count sampler and not a proof that every fibre is nonempty.

### Which induced supersets preserve the count?

The [odd-cycle expansion](https://arxiv.org/abs/2307.05569) gives
\(H(T)-H(T-v)=2\sum_{C\ni v}H(T-V(C))\), over directed odd
cycles. In a nontrivial strong component each vertex lies on a
triangle: otherwise all of its in-neighbors beat all of its
out-neighbors, contradicting strong reachability. Therefore
\[
H(T-v)=H(T)\quad\Longleftrightarrow\quad v
\text{ is a singleton strong component}.
\]
Applying the result one vertex at a time shows that two
same-count graphs are in an induced-subgraph relation only by
adding singleton components. Every fixed-count induced poset is
therefore a disjoint union of coordinatewise lattices
\(\mathbb N^{j+1}\), one for each singleton-free core word of
length \(j\): coordinates count singleton pads between and
outside its components. For \(h=1\) omit the empty graph under
the nonempty convention. This classifies the whole plateau
around a count, not just its smallest graph.

## 2. Exact census and the Mersenne test at 127

The [TSV](births.tsv) contains a representative (outgoing masks)
for every tied minimum, not one arbitrary choice. It is based on
the archived [order-seven and order-eight tables](../actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md).
Lower orders are recovered by deleting vertices from the complete
order-seven representatives and canonically deduplicating by
score cells. Counts are recomputed by subset path DP.

| Order | Newly attained counts | Classes in their full birth fibres |
|---:|---:|---:|
| 1 | 1 | 1 |
| 2 | 0 | 0 |
| 3 | 1 | 1 |
| 4 | 1 | 1 |
| 5 | 4 | 6 |
| 6 | 12 | 33 |
| 7 | 58 | 328 |
| 8 | 243 | 5,457 |

All odd counts through 609 except 7,21 occur in this table;
there are also higher attained counts. Across orders 4–8, every
birth class except the two for 63 has a vertex deletion whose
count was itself born at the previous order. For 63 the
one-deletion list is \(1,25,25,33,33,37,37,11\), and two
deletions can already reach the order-six birth of 17.
This is a bounded census of **birth-preserving ancestry**, not
an all-order rule. Its first exceptional step is exactly a
warning about enforcing single-vertex birth parents.

The 125 and 127 fibres at order seven each consist of a unique
strong self-converse class. Their canonical outgoing rows here
are
\[
(6,12,56,97,11,83,23),\qquad (10,20,33,70,77,27,39).
\]
The respective multisets of one-vertex deletion counts are
\(\{19,19,37,37,41,43,43\}\) and
\(\{29,29,29,31,31,37,37\}\); all are order-six births.
They share both of their actual count-37 parent classes,
not just the integer 37. The shared parents are two of the five
order-six H37 classes. The inherited aligned certificates cost
four old arc differences when preserving either common parent;
without that requirement some relabeling uses three.
This is an example of a nearer scalar neighbor that need not
be a nearer induced edit.

The Mersenne equality has different physical mechanisms at
different orders. The first 127 graph has counts 11, 23, 11 of
odd cycles of lengths 3, 5, 7, and nine disjoint triangle pairs:
\(127=1+2\cdot45+4\cdot9\).
For the first 125 graph these are 12, 28, 10 and six pairs:
\(125=1+2\cdot50+4\cdot6\).
Losing five individual odd cycles and gaining three pack pairs
changes the number by exactly two. In contrast the archived
[Fibonacci fan](../../../01-canon/theorems/THM-4748-two-reversed-arcs-exact-path-counts-the-seven-mod-eight-lines-and-the-mersenne-side-three-arcs-from-transitivity.md)
at exponent seven has nine vertices and one apex, so all its
odd cycles intersect at that apex: \(127=1+2\cdot63\).
Neither case says that all nonzero binary subsets index paths.
Also 255 first appears at order eight and 511 already appears
at order eight, whereas 63 first appears at order eight despite
its smaller exponent; the numerical exponent is not the birth
order function.

## 3. Birth of a growth rate is a separate invariant

For a marked graph write \((h,g,r,s)\) for the count, deletion
count, and weighted even/odd root-start path sums, as in the
[center-126 packet](../center_126_new_interface_20261116/README.md).
Let \(e\) be the full end count at the root. The near-sink
append gives both \(h+2r\) and \(h+(r+s)-e\), so \(e=s-r\).
Call a root stable when \(e=g\), and active when also \(r>0\).
Under the repeated near-sink append the growth slope is \(2r\).

An active root lies in the last SCC: otherwise no full path
ends there and \(e=0\). A preceding ordinal block with count
\(q\) multiplies \(h,g,r,s\) by \(q\); any nonempty suffix
destroys current stability. (It need not destroy the
possibility of a later different preparation.) It follows that
a fixed-count induced extension cannot create a new positive
stable slope. Earlier pads preserve a last-component port,
whereas a trailing pad can destroy stability. To gain a new
positive slope one must change a nontrivial core, not merely
increase the order.

Define \(m_{\rm act}(h,r)\) as the first order of an active
root with those parameters. Then the exact recurrence is
\[
m_{\rm act}(h,r)
=\min_{q\mid\gcd(h,r)}
\bigl(m_\emptyset(q)+m_{\rm str,act}(h/q,r/q)\bigr).
\]
Here \(m_\emptyset(1)=0\) allows no prefix, and other factors
use the ordinary birth order. All the tied pointed witnesses
are joins of minimum prefixes with tied minimum strong rooted
suffixes. Consecutive odd centers' endpoints \(h,h+2\) are
coprime, so their stable witnesses cannot share a preceding
ordinal multiplier \(q>1\). A multiplicative scalar trick by
itself will not build the paired interface.

Neither order-seven birth at 125 or 127 has an active root.
There are 12 classes at H125 and 10 at H127 at order eight.
Their first common positive slopes are exactly \(120,140,158\)
at that order. In particular the order-eight redesign for
slope 120 could not have been found by padding a birth graph.
The nine-vertex F7 fan above has \((h,g,r,s)=(127,1,107,108)\),
so slope 214 is first born at order nine (absent from the
complete order-eight 127 table).

## 4. An exact binary program for a matched center and a slope

Keep a marked old core. Append a transitive sequence of new
vertices which every old vertex except the root beats. A bit
\(w_i=1\) means the root also beats the new vertex; \(w_i=0\)
means that vertex beats the root. Require the word \(w\) of
length \(L\) to end in 0. The whole string is one literal
extension of the old graph, and the final root is stable
since every path of its deletion ends at the last new vertex.

Write \(x=X(w)\) for the forward binary value and \(z=Z(w)\)
for the reverse-complement binary value, so \(x\) is even,
\(2^{L-1}\le z<2^L\). Let \(f(z)\) be the one-apex path
count for the binary word of address \(z\) (the
[known binary recurrence](../prescribed_prime_30000_20261116/README.md)).
For a stable old root, the exact new coordinates are
\[
h'=h+2rz+g(f(z)-1),\qquad g'=g,\qquad
r'=2^L r+\frac{g x}{2},\quad s'=g+r'.
\]
More generally set \(b=s-g\). Replace \(2rz\) by
\((r+b)z+(r-b)[w_1=0]\) in the first formula and
\(2^Lr\) by \(2^{L-1}(r+b)\) in the third.

For the count, any odd cycle meeting the suffix must contain
the root. A cycle using only new vertices and the root gives
the one-apex contribution \(g(f(z)-1)\). A mixed cycle first
uses an old nontrivial root path and then an increasing new
segment ending at a losing suffix vertex \(j\). For \(j\ge2\)
its \(2^{j-1}\) choices split evenly by parity, giving weight
\(r+b\) per binary unit after counting the cycle factor 2.
The first position has only the singleton suffix path, so
uses \(2r\) rather than \(r+b\); this is the boundary term.
For the new \(r\), read the total split response \(R'=r'+s'\).
Each nontrivial old root path can choose any suffix subset,
contributing \(2^L(r+b)\). The singleton root can stop, or
enter at a 1-bit \(j\) and choose arbitrary later suffix
vertices, contributing \(g(1+x)\). Subtract \(e'=g\) from
\(R'\) to get \(2r'\). These arguments preserve all old arcs.

For **two seeds** with identical \(g,r,b\) and counts
\(h,h+2\), every common word preserves the +2 gap and
all three continuation coordinates. At fixed length \(L\),
the possible new \(r'\) from a stable seed form the exact
progression
\[
2^L r,\quad 2^L r+g,\quad \ldots,\quad
2^L r+g(2^{L-1}-1).
\]
Matching a slope selects the unique word via
\(x=2(r'-2^Lr)/g\), if that is an even integer in range.
There is **no additional binary choice left** to fix the
center at that same slope and budget. Conversely to select a
center \(C\) from stable \(h,h+2\), solve
\[
2rz+g f(z)=C-(h+1)+g
\]
by the exact binary inverse and reverse-complement its
word. This is a physical construction from an absolute
target whenever the inverse returns an address. For unequal
old \(g\), a common word with \(x>0\) generally breaks
slope equality; an equal-\(g\) source pair is a real resource.

### Explicit target hits

The supplied finite seed bank through order eight gives the
following literal examples. At 378 the old pair 69/71 differs
by a directed triangle reversal. At 30090 the 49/51 seed
pair differs by only two old arcs, and the same two edits
remain after all programmed and future near-sink growth.

| Center | Seed counts | \(g,r\) | Address \(z\) | Word \(w\) | New \(r'\) | Old edits |
|---:|---:|---:|---:|---|---:|---:|
| 378 | 69/71 | 13,32 | 4 | 110 | 295 | 3 |
| 10008 | 31/33 | 5,14 | 347 | 001001010 | 7353 | 3 |
| 30090 | 49/51 | 5,35 | 410 | 101001100 | 18750 | 2 |

The last row is a new **persistent** realization of the
prescribed 30089/30091 pair, different from the earlier
eleven-vertex witness. The slope at 378 is 590, not the
360 needed to continue a separate tripled-midpoint recipe.
The center finder scanned 2,808 stable keys and 134 matched
key pairs from the order-eight atlas; 30012 did not hit
in this restricted bank. None of these negative tests is
a nonattainment result for general tournaments.

## 5. When approximation selects a unique rate, and when it cannot

The programmable suffix has two independent coordinates
in a precise **normalized** sense. Put \(A=2r'/2^L\),
\(Y=C'/2^L\), \(U=x/2^L\), and \(V=z/2^L\). Then
\[
A=2r+gU,\qquad |Y-AV|\le (|h+1|+gL)/2^L .
\]
Indeed \(H(w)-1\) is the weighted forward \(1\)-to-\(0\)
pair sum plus the adjacent 10 correction. The same
product \(xz/2^L\) includes also the reverse-order pairs,
of total at most \(L\); the adjacent correction is also
at most \(L\), giving the stated error. Leading bits of
\(w\) prescribe \(U\), trailing bits prescribe \(V\),
independently once the budget separates them. The joint
limit region is \(2r\le A\le2r+g,\ A/2\le Y\le A\).
This does not imply rounding to an absolute integer:
the finite word is still unique at a given slope, and the
unscaled error bound grows with \(L\).

There is a different exact rounding window for an already
stable pair followed by ordinary near-sinks. For an active
stable root one always has \(g\le2r-1\): choose any out-neighbor
\(u\) of the root. Its contribution to \(r\) equals the
odd-start weighted sum at \(u\) in the deletion, namely
\((R(u)+e(u))/2\ge g/2\) by splitting full paths at \(u\).
As \(R=r+s=2r+g\ge h\), this gives \(h\le4r-1\).
Thus, for a stable pair with common \(r\) and initial counts
\(h,h+2\), after \(k\ge1\) ordinary appends the center
\(C\) satisfies, with \(L=2^k\),
\[
2r(L-1)+2\le C\le 2r(L+1)-2 .
\]
In particular if \(C\le L(L-1)\), the nearest integer to
\(C/(2L)\) is the **only** possible half-slope \(r\) and
equals it whenever such a pair exists. The arithmetic then
determines the required seed count \(h=C-1-2r(L-1)\).
The remaining exact test is the two active birth fibres
\((h,r)\) and \((h+2,r)\).

This selection does not solve arbitrary placement. Summing
over possible \(r\) gives only \(O(B^2/4^k)\) possible even
centers up to \(B\) for fixed depth \(k\) (each \(r\) gives
only \(O(r)\) offsets). Summing over
\(k\ge(1/2+\epsilon)\log_2 B\) yields
\(O(B^{1-2\epsilon})\). Long ordinary tails can identify
and maintain a rate, but cannot alone cover all centers.
The unresolved step is obtaining enough small active seed
fibres, or allowing different physical continuations when
the exact center and slope test fails.

## Replay

    python3 atlas.py ../actual_hamilton_multiplicities_20261021/session24-multiplicities
    python3 program.py
    python3 ancestry.py
    c++ -O3 -std=c++17 find_center.cpp -o /tmp/find-center
    /tmp/find-center ../actual_hamilton_multiplicities_20261021/session24-multiplicities/n8_classes.csv

The graph and profile checks in program.py use direct subset
DP on the actual masks, not the transfer recurrence. The
complete finite atlas replays recorded representatives and
induced subsets; all-order statements above are separately
proved from SCC order, insertion, and the displayed transfer.
