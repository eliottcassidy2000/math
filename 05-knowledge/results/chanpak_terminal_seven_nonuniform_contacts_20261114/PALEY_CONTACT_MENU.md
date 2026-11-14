# A 43-coordinate inverse palette for two nonuniform Paley contacts

The [one-contact obstruction](MINIMAL_PALEY_BREAK.md) is
structural: no one-apex choice escapes divisibility by 3.
For two contacts one can recover the *entire* variable-cut
response from exactly 43 canonical subset probes. The outside
reservoir can have arbitrary size, tournament structure, and
orientations to the two contact vertices.

Keep a Paley 7-orbit \(U=\mathbb F_7\) homogeneous to every
vertex except \(v_1,\ldots,v_k\). All outside arcs and the
directions among the \(v_i\) are fixed. The input bit \(x_{i,u}\)
records \(v_i\to u\); its complement records the opposite arc.
For every permutation of the vertices, write the indicator
that it is a directed path as a product of the directional
indicators on its consecutive arcs. A path uses at most
two arcs at \(v_i\). Therefore
\[
 H(\mathbf x)=
 \sum_{\substack{A_1,\ldots,A_k\subseteq U\\ |A_i|\le2}}
 c_{A_1,\ldots,A_k}\prod_{i=1}^k\prod_{u\in A_i}x_{i,u}
\]
as the unique multilinear Boolean polynomial with these
values. All coefficients are integral, and all nonconstant
coefficients are even (take Boolean differences of odd
Hamiltonian counts). No positivity of the signed coefficients
is claimed.

The affine group of \(U\) acts simultaneously on the sets \(A_i\).
It has 21 elements: the identity, six nontrivial 7-cycle
translations, and fourteen elements with one fixed point and
two 3-cycles. Acting on subsets of size at most two, these
classes fix respectively 29, 1, and 2 subsets. Hence the
number of coefficient orbits is
\[
                        (29^k+6+14\cdot2^k)/21.
\]
For \(k=2\), the orbit table by bidegree is
\[
 \begin{array}{c|rrr}
 |A_1|\backslash |A_2|&0&1&2\\\hline
 0&1&1&1\\
 1&1&3&7\\
 2&1&7&21
 \end{array}
\]
and the sum is **43**, independent of the outside size.
These are coordinates for a linear description, not a claim
that arbitrary integer assignments to the coordinates are
physically realizable.

To recover the coefficients, evaluate the ordinary count at
one representative of each input orbit with \(|A_i|\le2\).
The coefficient at \((A_i)\) is the Boolean Möbius difference
over its subtuples; all their count evaluations are in the
same bank after affine normalization. This exactly evaluates
every larger mask pair, even when neither mask has size at
most two. It also allows a direct inverse over mask pairs:
for a prescribed absolute target, evaluate one instance per
full mask orbit and reject that outside reservoir if none
matches. Among matches, lexicographically choose the least
pair or keep the full set of contact states. This is a
finite test per outside reservoir, not a proof that a
suitable outside reservoir always exists.

## Only 35 free coefficient directions modulo 3 are possible

A constant contact row (all 0 or all 1) is uniform to \(U\).
With at most two contact rows the one-contact theorem then
forces \(H\equiv0\) modulo 3 as a function of the other row.
Setting either row to 0 first removes the 5 coefficients with
at least one empty index set. Setting the first row to 1
and looking at the remaining row's singleton coefficient
forces \(c_{\{u\},\{u\}}\equiv0\): the stabilizer of \(u\)
has order 3, and all other nonempty subsets of size at most
two occur in triples under that stabilizer. Setting the second
row to 1 gives the same condition. For a fixed unordered
pair in the remaining row the stabilizer is trivial, so
its condition is a sum of the coefficients of total
bidegree \((1,2),(2,2)\) or \((2,1),(2,2)\), each with
multiplicity one. Those last two equations are independent:
one has a \((1,2)\) coordinate absent from the other,
and conversely. The resulting formal invariant space has
dimension 35 over \(\mathbb F_3\). Actual outside reservoirs
may realize a smaller set.

The three levels here separate: "one nonuniform contact" owns
the forced factor 3, "two contacts" owns a finite but genuinely
richer palette, and "choose the outside reservoir" still owns
the absolute target. The decoder makes that last obligation
testable without conflating a symmetric profile with coverage.

Run `python paley_menu.py` for the orbit count, Möbius recovery,
and literal/symbolic comparison including an outside vertex
with unequal directions to the two contacts.
