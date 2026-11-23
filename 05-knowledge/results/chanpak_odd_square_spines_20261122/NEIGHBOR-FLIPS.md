# The two near brackets are isolated pairs in the one-arc metagraph

**FINITE-EXACT.** Put \(N=2^{60}(2^{61}-1)-1\). Consider the graph on
ordinary isomorphism classes of 118-vertex tournaments, adjacent when one
physical arc reversal suffices. For each of the [two explicit pairs of
near witnesses](NEAR-WITNESSES.md), take the induced subgraph on
\(\lvert H-N\rvert<B\), with the corresponding \(B\) in the table below.
The component of that pair consists of **exactly its two original classes**.
Enlarging the band to \(\lvert H-N\rvert\le B\) admits a third class adjacent
to one of them. In particular every sequence of flips that starts at the
given bracket and reaches \(E_{61}-1\) or \(E_{61}+1\) must first reach
a tournament at least \(B\) away from \(N\).

| Signed offsets of the original two counts from \(N\) | Sharp \(B\) |
|---|---:|
| \(-844{,}710,-844{,}708\) | 9,860,704,255,316,668,080 |
| \(184{,}814{,}588,184{,}814{,}590\) | 12,363,620,044,828,489,432 |

Each of the four labelled tournaments has 6,903 possible flips. Exactly 51
are neutral in each case. Every one of those 51 swaps tournament twins, two
vertices with identical incidence to every other vertex. Reversing the arc
between them therefore just transposes their names. Exactly one further
flip remains in the stated open band, namely the marked arc to the other
class. Distinct path counts distinguish the two classes. These facts
give the component theorem by transporting incident flips along any
isomorphism. The rows below count distinct Hamilton path counts among all
one-arc neighbors, including neutral ones.

| Source offset | Distinct neighbor heights | Smallest \(\lvert H-N\rvert\) outside the two original heights |
|---|---:|---:|
| \(-844{,}710\) | 6,831 | 9,860,704,255,316,668,080 |
| \(-844{,}708\) | 6,831 | 13,955,597,669,012,648,378 |
| \(184{,}814{,}588\) | 6,820 | 19,411,975,011,021,988,434 |
| \(184{,}814{,}590\) | 6,820 | 12,363,620,044,828,489,432 |

Another way to transport a unit flip is to choose a second arc and reverse
it on *both* original hosts. The resulting two opposite corners are still
separated by their marked physical flip. In the first pair 204 choices of
second arc preserve an absolute difference of two. In the second pair 203
choices do. In each census the 51 twins are exactly the choices nearest
the perfect-number midpoint. A two-flip search that first leaves the
induced band remains a different question.

## Counting every physical correction outside the original chart

**PROVED algorithm.** Suppose that deleting one distinguished vertex \(a\)
from a tournament leaves an ordinal sum of blocks of size at most \(s\).
For a candidate physical arc \(u\to v\), compute the number of spanning
paths using the adjacent ordered word \([u,v]\). Also compute the count
with \([v,u]\), formally enabling the reversed arc. The second count minus
the first is the exact change under reversal, since all paths avoiding both
words cancel. To contract an ordered pair \([i,j]\), incoming incidences
are tested at \(i\) and outgoing ones at \(j\). The contracted vertex can
have both or neither orientation with another vertex. [Inherited arc
contraction](../chanpak_framed_inverse_20261115/FRAMED_SIMPLE_GRAPH.md) ·
[THM-5267, arc-flip identity](../../../01-canon/theorems/THM-5267-seams-at-a-prescribed-centre-by-one-port-contact-the-arc-flip-lemma-deletion-monotonicity-the-order-10-extension-graph-and-the-heaviest-parent-tree-of-first-creators.md)

Delete \(u,v,a\) from the original blocks and try both orders of the
contracted pair and \(a\). For either order call the two separators \(P,Q\).
A path has three ordered slots before, between and after them. In each slot
the occupied blocks occur in increasing ordinal order. Within each small
block enumerate its vertex permutations with two cuts. Require each of the
three pieces to be an internal directed path. Across different blocks the
adjacent arcs are automatic. As blocks are read from left to right, the
state needs twelve possibilities:

* The prefix is empty or ends with a usable arc into \(P\), versus an
  unusable current end.
* The middle is empty, or has a valid start and one of two current end
  states for the arc into \(Q\).
* The suffix is empty or has already acquired a valid start after \(Q\).

Any new nonempty middle or suffix tests its start precisely when that slot
was empty. At the end accept an admissible prefix, and a middle whose end
is admissible or which is empty with \(P\to Q\). A suffix that stays empty
is allowed. If the tested arc touches \(a\), there is only one genuine
separator. A sentinel placed first and beating every vertex gives exactly
this two-slot calculation using the same transfer. Each path induces its
unique separator order, block allocation and local permutation. The transfer
therefore counts every permitted path once. Arithmetic uses GMP integers
throughout. The replay constructs the literal row matrices, covers all
\(4\binom{118}{2}=27{,}612\) physical reversals, and checks that the only
in-band arcs are the twins and the marked arc. The small control compares
378 individual reversals on 27 smaller hosts with direct subset enumeration.

With GMP on the compiler search path, run from this directory:

~~~sh
c++ -O3 -std=c++17 neighbor_flips.cpp -lgmpxx -lgmp -o /tmp/e61-neighbor-flips
python3 -B neighbor_flips.py /tmp/e61-neighbor-flips --check
python3 -B neighbor_flips.py /tmp/e61-neighbor-flips
~~~

For Homebrew GMP add \(\texttt{-I/opt/homebrew/opt/gmp/include}\) and
\(\texttt{-L/opt/homebrew/opt/gmp/lib}\) to the compile command.
