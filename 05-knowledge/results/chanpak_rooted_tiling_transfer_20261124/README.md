# A local classifier for fixed-path tilings

**PROVED** for the local classifier, strip recursion, and structural
consequences below. **FINITE-EXACT** for the marked H=11 and cycle examples.
This extends [the exact complement-line wiring](../session25_spectrum_wiring_20261021/SYNTHESIS.md)
and [eternal colours, THM-4930](../../../01-canon/theorems/THM-4930-eternal-colours-which-odd-path-counts-are-ever-self-converse.md).
It does not prove that every target multiplicity appears.

## Class membership by one-vertex moves

Normalize a directed Hamilton path to \(0\to1\to\cdots\to n-1\).
A tiling is the bit string of the remaining orientations. Let \(\tau\)
be positional reflection followed by edge reversal (the converse
normalized on the reversed path); let \(\phi\) flip all nonspine tiles.
The ordinary isoclass fibre has \(H(T)/A(T)\) tilings, where
\(A=|\operatorname{Aut}T|\). Indeed normalizing the various Hamilton
paths of one \(T\) gives the same string precisely on an automorphism
orbit, and the automorphism action on a linear order is free.
By [Rédei](https://acta.bibl.u-szeged.hu/13432/), \(H\) is odd, so both
\(A\) and this fibre size are odd. The merged fibre is single if \(T\)
is self-converse, and otherwise double with mass \(2H/A\equiv2\bmod4\).
These masses count tilings, not paths before dividing out symmetries.

**Local generation theorem.** A move deletes one vertex from the
spine, provided the shortened word is a directed path, and reinserts
it in any slot that gives a directed path. Renormalize using the new
word. Two tiling strings belong to the same ordinary isoclass iff
one can reach the other by such moves. In particular arbitrarily
permuting labels is not needed as a primitive operation.

*Proof.* For a labelled tournament consider the corresponding graph
of all Hamilton paths. Fix a vertex \(v\); let \(S_v\) be the set of
paths from which \(v\) is removable. Over each path of \(T-v\),
the insertions of \(v\) form a clique, and at least one insertion
exists: either \(v\) can go at one end or the string has a consecutive
pair \(a,b\) with \(a\to v\to b\). Inductively the graph on \(T-v\)
is connected. To lift an edge in that graph moving \(w\), delete \(w\)
first, leaving a common word of \(T-\{v,w\}\). Choose any admissible
gap for \(v\) in that common word. Each of the two original slots
for \(w\) can still be used, because if \(v,w\) use the same gap,
both fit there in the order of their mutual arc. This gives a lifted
edge moving \(w\) with both endpoints in \(S_v\), so \(S_v\) is
connected. The same independent-gap argument gives
\(S_v\cap S_w\ne\varnothing\) for every \(v,w\) (use the empty word
in the two-vertex base case). Every path has a removable end and
lies in some \(S_v\). Thus the labelled graph is connected.
A normalized move never changes the tournament isoclass, and
every tiling in the isoclass comes from a labelled path, proving
the claim.

**Exact classifier and transport.** On the rooted graph, take the
least bit mask of a component as its canonical root. For the merged
class use the lesser of the two keys for the component and its
\(\tau\)-image. Keep on every move the permutation \(q\) sending new
positions to old positions. If \(\sigma_x\) sends the positions of
mask \(x\) to the starting labels along a chosen tree, a directed edge
\(x\xrightarrow{q}y\) gives the loop permutation
\(\sigma_x q\sigma_y^{-1}\). These generate exactly the starting
automorphism group. Indeed every such permutation is an automorphism
by comparing the two transports of the same tiling; conversely every
automorphism sends the initial word to a labelled path which local
generation reaches, giving a loop with that permutation. Products
of the edge loops telescope to the transport of any loop. Thus the
physical path count is recovered from the rooted component **and**
its permutation labels. The replay implements this classifier.

**Sharp warning: \(\phi\) is not a map of components.** At order five,
use outgoing masks \((14,20,8,18,5)\) and move vertex 1 in
\(01234\) to make \(02314\). This keeps both endpoints and changes
the rooted tile mask \(19\to35\). Their complements have path counts
11 and 13, so these complements are not isomorphic. At order at
most four, a genuine interior move with fixed endpoints would swap
the only two interior vertices, requiring both orientations of the
same arc. This counterexample is therefore smallest for fixed endpoints.
The endpoints alone do not repair the loss. In fact on the original
labels, if \(P\) is the root, its tile-complement tournament has
outdegree
\[
 n-1-d^+_T(v)+{\bf1}_{v=P_{\rm start}}-{\bf1}_{v=P_{\rm end}},
\]
since it is the converse with the spine arcs restored. Hence the two
complements in the counterexample even have the same labelled score
vector. Under an interior move the complement changes by reversing
two triangles through the moved vertex, at the old and new gaps.
Each preserves the scores.

There is a further local color restriction: in a rigid self-converse
tournament the blue roots are an **independent set** in the local
re-rooting graph. Both blue roots would induce the same anti-involution.
The change of root must then commute with positional reversal. A
nontrivial single-vertex relocation is a cyclic shift on one interval;
commuting with reflection requires that interval to be reflection-fixed,
and conjugation by reflection inverts the shift. So it would have to
be a shift of order two, an adjacent swap, which cannot be a directed
path in both orders. In particular any local route between distinct
blue viewpoints of a rigid class makes at least one black visit.

## Blue/black outer-shell recursion on physical lines

A raw line is \(\{t,\phi t\}\). Blue means \(t=\tau t\) and black means
otherwise. For \(n\ge3\) a line is fixed under \(\tau\) iff it is
blue: \(\tau\) cannot swap \(t,\phi t\) because the end-to-end tile
is reflection-fixed. Keep separate the operation of merging converse
**nodes** and of identifying mirror **lines**. If identifying mirror
lines, each black line orbit has weight two and each blue orbit
weight one to recover incidence in the original tilings; loops have
two incidences in either convention.

Delete the outer spine pair from order \(n+2\). Restriction commutes
with \(\phi\) and \(\tau\); the new strip has \(2n-1\) tiles, one
fixed (the new closing tile) and \(n-1\) reflected pairs. Define
a defect as a reflected pair with different bits. Over a given raw
parent line, the number of raw children with exactly \(k\) added
defects is
\[
 2^n\binom{n-1}{k}.
\]
Choose the unique child endpoint restricting to a designated parent
endpoint. Pick the defect positions, one bit per reflected pair,
and the new closing bit. If the parent is black, every child is
black and the displayed counts also give the number over the parent
**mirror orbit** (choose one raw parent sheet). Over a blue parent
there are \(2^n\) blue children and, for each \(k>0\),
\(2^{n-1}\binom{n-1}{k}\) black child mirror orbits. The closing bit
is free on the strip, so over a raw parent the new closed endpoint
restricts half the time to each parent endpoint. In \(r\) raw strips
the defect-increment generating polynomial is
\[
 2^{r(n+r-1)}(1+z)^{r(n+r-2)}.
\]
This is a physical placement law, not a recursion on the counts in
a class alone. For a class-level obstruction, the order-five pair
of roots 4 and 9 in the H=11 fibre below restricts to order-three
masks 0 and 1, respectively, giving cyclic and transitive tournaments.
So the extension rule must keep the root port even inside a fibre.

## A path closes at one specific tile

A rooted tiling closes iff the end-to-end edge points from last to
first. This bit is \(\tau\)-invariant and flipped by \(\phi\).
Consequently every blue **and** every black physical line has exactly
one closing and one open endpoint, giving a directed line from the
first to the second. For a class \(T\), with \(C(T)\) directed Hamilton
cycles counted modulo rotation, closing incidences number \(nC/A\),
and open ones \((H-nC)/A\). These double for a chiral merged class.
Closing roots occur in cyclic necklaces of length
\(n/|\operatorname{Stab}_{\operatorname{Aut}T}(c)|\) per cycle orbit:
the cycle stabilizer acts by rotations and has odd order. In
particular every necklace length is divisible by \(2^{v_2(n)}\).

A finer blue count follows by fixing an anti-involution \(\alpha\)
of a self-converse tournament. Such an involution exists since the
extended automorphism group has an odd-order subgroup of index two.
It fixes at most one vertex (two would have an arc that must reverse),
so it fixes one vertex at odd order and none at even order. Its action
on a fixed directed Hamilton cycle, with order reversed as well,
reflects the cyclic word. At odd order it fixes one cut edge, at
even order it fixes two opposite cut edges. Thus an invariant cycle
gives respectively one or two \(\alpha\)-palindromic closing paths.
All anti-involutions are conjugate by an automorphism: their product
has odd order, and conjugation by a suitable power of that product
takes one involution to the other. A palindromic path determines
its anti-involution uniquely. Dividing the palindromic path count
for one chosen \(\alpha\) by its automorphism centralizer therefore
gives the number of blue rooted tilings. Its closing subset has size
\[
 \frac{(1\ {\rm if}\ n\ {\rm odd},\ 2\ {\rm if}\ n\ {\rm even})
          |\operatorname{Fix}(\alpha\text{ on directed cycles with reversal})|}
      {|C_{\operatorname{Aut}T}(\alpha)|}.
\]
In particular blue closing mass at even order is even, and the
odd blue fibre leaves an odd number of open blue roots. Summed
over self-converse classes, the blue closing count is exactly half
of all reflection-fixed strings (the closing tile is a free
reflection coordinate). This is a different aggregate from half of
the raw tilings, and makes the cycle restriction visible even after
merging.

For even-order self-converse tournaments a further invariant is
\(H\equiv1\bmod4\). By [Grinberg–Stanley's odd-cycle expansion](https://arxiv.org/abs/2307.05569),
\(H=\sum_{\mathcal S}2^{|\mathcal S|}\) over collections of vertex-disjoint
directed odd cycles, including the empty collection. A single invariant
odd cycle would require a fixed vertex of \(\alpha\), impossible here.
Thus odd cycles pair, and modulo eight the two-cycle collections
cancel the single-cycle pair contribution exactly when the two
reflected cycles are disjoint. More explicitly,
\[
 H \equiv 1+4\,\#\{\{c,\alpha c\}: c\text{ odd},\ c\cap\alpha c\ne\varnothing\}
 \pmod8,
\]
with reverse order understood in \(\alpha c\).
This rules out even-order blue-to-blue \(+2\) and \(p\to7p\)
moves, even with different roots.

**Distinct spectra.** Add two vertices \(b,c\) with
\(T\to b\to c\to T\). Every directed Hamilton cycle traverses
\(b,c\) consecutively and contains exactly one Hamilton path of \(T\),
so its cycle count is exactly \(H(T)\), and this correspondence
preserves any marked flips of old arcs. Conversely cycles already
attain 7 and 21, with outgoing masks respectively
\((18,36,9,19,38,77,31)\) and \((114,68,9,19,38,78,156,63)\),
checked by the replay. Thus those path holes do not persist under
changing the rooted observable. Conditional on Chan–Pak, all
positive odd cycle counts occur by the wrapper and those two seeds.

## The visible 2 and 4 at 11

Use lexicographically ordered nonspine tiles with 1 = forward.
The unique order-five \(H=11\) class has rooted masks
\[
 4,6,9,13,20,22,33,40,43,44,57 .
\]
Blue roots are \(4,22,33\); closing roots are \(9,33,40,43,57\),
so this class has one Hamilton cycle. Complementing with XOR 63
gives exactly **two internal black lines**, \(\{6,57\},\{20,43\}\),
which use **four** distinct path viewpoints. The lines are mirrors
of each other, giving one loop of weight two if line orbits are
identified. The unweighted fundamental black incident degree is
four: this loop contributes two, and two further black line orbits
leave the class. One leads to a doubled H=5 node (two converse
classes), and the other to a single self-converse H=13 node. So
one actual neighbourhood exhibits both odd and doubled-odd masses.
This is one concrete geometric reading of “11
becoming 2 and 4,” with both counts tied to physical placements.
It is distinct from replacing the 11 test in the prime sieve.

For any prime \(p\ge5\), every tournament with \(H=p\) is rigid,
at every order. If its automorphism group were nontrivial, its
order would have to be \(p\) (it acts freely on paths). On a
vertex orbit of its order-\(p\) generator, the induced tournament
is translation-invariant on \(\mathbb Z/p\mathbb Z\). Each of
its \((p-1)/2\) forward differences gives a distinct directed
Hamilton cycle, by repeatedly adding the difference. Their cuts
give at least \(p(p-1)/2>p\) paths. Such paths inject into the
full tournament by inserting the omitted vertices one at a time
in a deterministic valid slot; deleting the added vertices
recovers the original word. This is a contradiction.

In particular the exact prime-order distribution reduces to a core
polynomial. If \(B_p(z)\) counts strong isoclasses of path count \(p\)
by order and \(B_p^{\rm sc}(z)\) their self-converse part, the ordinary
class, undoubled merged-node, and doubled merged-node generating
functions are respectively
\[
 \frac{B_p(z)}{(1-z)^2},\qquad
 \frac{B_p^{\rm sc}(z)}{1-z^2},\qquad
 \frac12\left(\frac{B_p(z)}{(1-z)^2}
                  -\frac{B_p^{\rm sc}(z)}{1-z^2}\right).
\]
Exactly one strong component contributes \(p\), since path counts
multiply over the ordinally ordered SCCs, each nontrivial one having
count greater than one. The other SCCs are singleton padding on
the two sides. Converse swaps these sides and takes converse of
the core, so symmetry requires identical padding sizes and a
self-converse core. By rigidity, the two merged masses are \(p,2p\),
not \(p/A,2p/A\). At \(p=11\), the checker gives one core at order
five, self-converse, and none through order four. The
[Busch strong bound](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3)
excludes order at least six. Thus every odd order \(n\ge5\)
has one undoubled H=11 node, every even order has none, and
order \(n\ge5\) has \(\lfloor(n-4)/2\rfloor\) doubled H=11 nodes.
All padded copies at \(n>5\) are open at every root.

### Separate arithmetic interpretation

Removing 11 from the preceding [fixed-bank sieve](../chanpak_lrc_two_edit_ports_20261123/FILLER-KEITH-ULAM.md)
and testing divisibility by 2 or 4 instead adds no restriction (all
starts are odd). The support test then has exactly two maximal-size
supports, the two older 13-member supports each enlarged by \(8a,8b\).
The replay checks that at full support just one term survives at
each remaining prime and all three terms vanish modulo 11.
Thus the sharp dyadic exponent rises from 12 to **14**, but the
leading new families are all composite. Even separating odd counts
into phases 1 and 3 modulo 4 splits both maximal families equally
by the existing binary lifts (the seed coordinates \(0a,1a,2a\)
are retained). A parity split cannot substitute for the missing
11 test. Also, on the prior Kaprekar mod-11 sign quotient the
nonzero clock has \(1\to2\to4\) under doubling, while 11 itself is
zero and does not move.

**Next target.** The classifier finds exactly which placements are
one isoclass and detects the doubled odd masses. Neither its
connectivity nor the strip recurrence forces a component of any
prescribed size. A Chan–Pak induction would need to control the
sizes of the rooted components under extension, keeping the
physical port and its transport, not just the blue/black totals.

Replay: `python3 rooted-transfer-check.py`.
