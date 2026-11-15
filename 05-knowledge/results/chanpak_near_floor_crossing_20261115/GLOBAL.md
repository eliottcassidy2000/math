# Cycle surplus in Busch's merge: an unrestricted near-floor theorem

**PROVED, for every strong tournament \(T\) of order \(n\ge3\).**
With \(H\) the number of directed Hamiltonian paths, \(C\) the number of
directed Hamiltonian cycles up to cyclic rotation, and \(f(n)\) Busch's
strong minimum,
\[
 H(T)\ \ge\ f(n)+3(C(T)-1).
\]
In particular this proves [HYP-9457(a)](../../hypotheses/HYP-9457-single-cycle-floor-band-and-lorentzian-join-crevasses.md)
at all orders, including its equality clause: equality in the conjectured
coefficient-two inequality holds exactly when \(H(T)=f(n)\).
This is a lower bound, not coverage of prescribed counts or a closure
of Chan–Pak.

## 1. Paths in a one-vertex gluing

Let \(A,B\) be induced strong subtournaments of \(T\), each on at least
three vertices, with union \(T\) and intersection \(\{v\}\).
Call a path of \(T\) broken if its linear restriction to either
factor is not a path. Other paths fall into fibres indexed by their
ordered pair of restricted factor paths.

Every fibre is nonempty. In fact, any two disjoint ordered directed
paths can be merged by comparing the first remaining vertex of each
list and taking the one that beats the other; after one list has
emptied, use the rest of the other. Each chosen vertex beats the
next remaining front of its own list and the current front of the
other, so the merged sequence is a path. Apply this separately to
the two prefixes before \(v\) and the two suffixes after \(v\) in the
chosen factor paths. The result links to \(v\) at each junction by
the corresponding factor arc. This is [Busch's merge argument, Lemma 2
and the proof of Theorem 1](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3),
also expressed as the first part of THM-4745. Consequently,
\[
 H(T)-H(A)H(B)=\#\{\text{broken paths}\}
                 +\sum_{\text{good fibres }F}(|F|-1).
\]

## 2. At least three paths of surplus per extra cycle

Call a cycle of \(T\) bad when its **cyclic** restriction to at least
one factor is not a directed cycle. If a factor cyclic restriction
has an invalid consecutive pair, then a rotation into a path can
avoid displaying that pair in its linear restriction only when its
cut is on the directed segment of the full cycle between the two
vertices in question. There are at least \(|A|-1\) cuts outside that
segment if the offending side is \(A\), one per other consecutive
arc segment in the restricted circular order. For \(|A|\ge4\)
this gives at least three broken paths. For \(|A|=3\), strength
makes \(A\) a directed triangle. The only invalid circular order is
its reverse, in which all three arcs are invalid; cutting at one
position removes only one, so every rotation is broken. The same
argument applies if the offending side is \(B\). Distinct cycles
produce disjoint sets of rotations, as the endpoints of any
closable path determine its unique closing arc. Thus each bad cycle
charges three distinct broken paths.

Group the good cycles by their pair of restricted factor cycles.
Within such a group, write the factor cycles starting at \(v\)
as \((v,a_1,\ldots,a_p)\) and \((v,b_1,\ldots,b_q)\), with
\(p,q\ge2\). A good cycle in the group is an admissible shuffle
of the two lists after \(v\). If the group has \(k\) members, then
rotating them to start at \(v\) gives \(k\) distinct paths in one
fibre. Rotating them to end at \(v\) gives \(k\) distinct paths in
a second fibre. This contributes \(2(k-1)\) surplus.

Now cut immediately after the first letter of each shuffle (after
\(a_1\) or \(b_1\)). There are \(t\in\{1,2\}\) different fibres so
obtained, determined by which factor supplied the first letter,
giving another \(k-t\) surplus. All these fibres have \(v\)
interior and are distinct from the preceding two.

If \(t=2\), suppose \(a_1\to b_1\), interchanging the factors
if needed. Choose a cycle beginning its shuffle at \(b_1\).
Immediately after the first occurrence of \(a_1\), its prefix is
\((b_1,\ldots,b_j,a_1)\) for some \(j\ge1\). This prefix can be
replaced by \((a_1,b_1,\ldots,b_j)\), since \(a_1\to b_1\)
and the \(b\)-arcs are factor-cycle arcs. Cut the full cycle just
after the original prefix. Its rotation is the complementary
shuffle suffix, then \(v\), then the original prefix. Replacing
the prefix produces another valid path in the same fibre:
\(v\) beats either possible first letter, and the restriction to
each factor is unchanged. The suffix is nonempty because \(p\ge2\).
The fibre is distinct from the first-step fibres, since both
factor prefixes are nonempty. This adds one more surplus, giving
at least \(3(k-1)\) whether \(t=1\) or \(t=2\).

No two groups charge the same good fibre. Every charged fibre's
two paths, when wrapped, recover the two original factor cycles.
And no good fibre is ever charged against the broken-path term.
Write \(b\) for the number of bad cycles and \(g\) for the number
of nonempty good groups. Adding all charges proves
\[
 H(T)-H(A)H(B)\ge 3b+3(C(T)-b-g)=3(C(T)-g).
\]
Since \(g\le C(A)C(B)\), we obtain the strong gluing inequality
\[
 \boxed{H(T)-3C(T)\ge H(A)H(B)-3C(A)C(B).}
\]

## 3. Apply Busch's split and induct

[Busch, Lemma 1 and Corollary 1](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3)
give the split above for every strong tournament of order at least
five unless it is nearly transitive, and give
\(f(3)=3,f(4)=5,f(5)=9,f(n+3)=5f(n)\).
A nearly transitive tournament is a transitive list with its
last-to-first arc reversed; all spanning cycles use that arc and
the remainder must be the unique increasing Hamiltonian path
from first to last. Thus \(C=1\) in this case. The orders 3 and 4
are also immediate: the only strong isomorphism class at each
order is nearly transitive. For order 4 one may pick a directed
triangle; the fourth vertex has both a win and a loss on the
triangle, and up to relabeling and converse gives the unique
strong class, which is self-converse.

For the split case put \(a=C(A)-1\), \(b=C(B)-1\), and
\(f_A=f(|A|)\), \(f_B=f(|B|)\). Inductively, \(H(A)\ge f_A+3a\),
\(H(B)\ge f_B+3b\), and \(a,b\ge0\) by Camion's theorem. Hence
\[
 H(T)-3C(T)\ \ge\
 f_Af_B-3+3a(f_B-1)+3b(f_A-1)+6ab
 \ \ge\ f(n)-3.
\]
The last step uses \(f_A f_B\ge f(|A|+|B|-1)\); by the displayed
period-three formula it is the six residue-pair comparisons,
with ratios 1, \(27/25\), or \(81/75=27/25\). Equivalently
each \(f(m)\) is the minimum of \(3^i5^j\) with
\(i,j\ge0,\ 2i+3j=m-1\), and multiplying preserves this
constraint. This completes the induction.

For the original coefficient two, \(C>1\) gives strictness.
If \(H(T)=f(n)\), the proved inequality and existence of a cycle
force \(C=1\). Conversely, equality with coefficient two
therefore has \(C=1\) and \(H(T)=f(n)\).
By parity the new lower bound also gives, for example, at least
\(f(n)+4\) paths when \(C=2\).

**Sharpness.** The order-five tournament with outgoing bitmask rows
\((10,4,25,18,3)\) on vertices \(0,\ldots,4\) has \(H=15\)
and \(C=3\) whereas \(f(5)=9\). The cycles with vertex 0 first
are \((0,1,2,3,4)\), \((0,3,1,2,4)\), and \((0,3,4,1,2)\).
The split into \(\{0,1,2\}\) and \(\{0,3,4\}\) gives two strong
triangles, so the local gluing coefficient and the unrestricted
coefficient are both sharp. [Independent DFS](verify_global.py) and its [replay](global_replay.txt)
check the witness and 7,936 literal strong-overlap controls.

The ribbon transfer route is still informative: the
[span-2 result](NEAR_FLOOR.md)
computes the *exact* ribbon floor and its full eventual
affine profile frontier. The unrestricted argument is instead
a positive merge-surplus certificate. Neither counts nor
classifies every integer in the global Hamiltonian spectrum.
