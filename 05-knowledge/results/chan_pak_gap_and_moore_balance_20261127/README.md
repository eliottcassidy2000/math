# Two kinds of Chan–Pak gap and the balance of rooted Moore encodings

2026-11-27. The [module-reversal theorem](../pants_diamond_module_reconstruction_20261127/README.md)
gives a useful way to separate one-vertex balance, two-vertex balance,
and the packing information that a Hamiltonian count can still see.
This note computes both local gap crossings and a sharp obstruction
in a dense tournament encoding of the degree-57 candidate.

## 1. Two distinct balance coordinates

For a tournament on \(n\) vertices let \(A\) be its signed adjacency,
\(G=A^2\), \(d_i^+\) the outdegrees, and \(c_3\) the directed triangle
count. Counting transitive triples by their unique source gives
\[
 c_3=\binom n3-\sum_i\binom{d_i^+}{2}
 =\frac{n(n^2-1)}{24}
 -\frac12\sum_i\left(d_i^+-\frac{n-1}{2}\right)^2.
\]
The transitive end has \(c_3=0\). The score-balanced end is regular
for odd order and has the two nearest degrees for even order.
This is a score statistic, not a complete ranking by Hamiltonian
path count.

There is a second axis that survives every cut switch. Let \(D\)
be the number of four-vertex diamonds (a directed triangle with a
source or sink). A four-set has signed Pfaffian of magnitude
three exactly for a diamond, and otherwise magnitude one; this
also follows from the earlier four-vertex pants calculation.
Writing \(\mathcal E=\sum_{i\ne j}G_{ij}^2\), the second and fourth
characteristic coefficients of \(A\) are
\(e_2=\binom n2,\ e_4=\binom n4+8D\). Newton's identity gives
\[
 \mathcal E=2\binom n2^2-4\binom n4-n(n-1)^2-32D .
\]
Thus fewer diamonds means greater square imbalance in this norm.
Zero diamonds includes transitive tournaments but is not restricted
to them. On even orders where a skew conference tournament exists,
\(\mathcal E=0\) is the perfectly balanced endpoint. On odd orders,
\(|G_{ij}|\ge1\) for \(i\ne j\), so
\(\mathcal E\ge n(n-1)\); equality characterizes a doubly regular
tournament up to cut switching (proof in Section 3).
The [four-cuff module theorem](../pants_diamond_module_reconstruction_20261127/README.md)
makes this same diamond local: in a regular module, its parity is
the sum of four rooted six-path detector bits.

## 2. What actually crosses the two missing counts

The all-order exclusions of 7 and 21 are already proved in the
[repository writeup](../../../06-writeups/research-document-for-mathematicians.md).
A Hamiltonian path visits strongly connected components in their
unique order, giving a product of the component counts.
[Busch's strong-component lower bound](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3)
is \(H\ge5^{(s-1)/3}\) on a strong component of size \(s\).
Thus a count of 7 would require a strong component of size at most
four with count 7. For 21, with 7 unavailable, it would require
a strong component of size at most six with count 21.
The replay enumerates the relevant labelled strongly connected
tournaments:

| Order | Strong-component spectrum |
|---|---|
| 3 | 3 |
| 4 | 5 |
| 5 | 9, 11, 13, 15 |
| 6 | 15, 17, 19, 23, 25, 27, 29, 31, 33, 37, 41, 43, 45 |

So the lower gap separates small strong cores. The upper gap
survives inside one strong-core size. Here are exact explanations
of what the two neighboring realizations can do differently.

On five vertices start with the increasing orientation and reverse
\(0\to3\). Its strong four-vertex core plus a singleton has count 5.
Now reverse \(2\to4\) to make the entire tournament strong, with
count 9. The triangle count goes from 2 to 3 and the directed
five-cycle count from 0 to 1, while the diamond count stays at 2.
This moves on the score axis without moving on the square-imbalance
axis.

On six vertices, start with the increasing orientation and reverse
\(0\to3\) and \(0\to5\). Reverse the directed triangle on
\(\{0,2,5\}\). The outdegree at **each labelled vertex** stays
fixed, and both tournaments are strongly connected. Both have five
directed triangles and four directed five-cycles. The first has no
disjoint pair of directed triangles; the second has exactly one.
By the [Grinberg–Stanley odd-cycle packing formula](https://arxiv.org/abs/2307.05569),
on at most six vertices
\[
 H=1+2c_3+2c_5+4d_{33}
\]
where \(d_{33}\) counts unordered disjoint pairs of directed
triangles. This move takes \(19\) to \(23\) entirely by creating a
pair, not by adding a triangle or a pentagon. The diamond count
increases from 3 to 5, so square imbalance falls from 144 to 80.
This illustrates the missing two units on this branch. The
all-order exclusion still uses the strong-core reduction above.

These gap crossings do not come from the invisible regular-module
reversal. The exact-square census at order six has only one
nonconstant Hamiltonian fiber type, \(\{13,17\}\), so keeping
all three-, four- and five-path *cut responses* fixed is stronger
than keeping just the score and ordinary triangle count fixed.
A proper nontrivial regular module in a small strong core must
be a triangle, since a quotient on only two vertices cannot
be strong and a regular module of size five would then need
order at least seven. At order five its quotient is a triangle
and the inflated count is 15; at order six the replay checks the
four-vertex strong quotients with a distinguished position,
getting 27 or 33. Thus regular-module signs cannot account for
the core of a 7 or 21 witness.

As a reusable internal construction, substituting a tournament \(M\)
for one vertex of a directed triangle gives
\(H=3H(M)+C_2(M)\), with \(C_2\) the number of ordered pairs of
nonempty directed paths partitioning \(M\). The three contiguous
placements supply the first term; otherwise the two outside
vertices must occur consecutively between the two parts of \(M\).
For a transitive chain of size \(k\), \(C_2=2^k-2\), giving
\(2^k+1\). For a four-vertex \(M\), \(C_2=14+4c_3(M)\) by
splitting into sizes \(1+3,3+1,2+2\); the result is
\(17+10c_3(M)\in\{17,27,37\}\).
This is an exact local grammar for how that branch skips 7 and 21,
not a claim that the grammar covers the full Chan–Pak spectrum.

## 3. Moore transport versus tournament balance

For a [diameter-two Moore graph](https://people.orie.cornell.edu/dpw/orie6334/HoffmanS60.pdf)
of degree \(d\), root at a vertex.
Its second neighborhood consists of \(d\) branches of size
\(m=d-1\), indexed by the first neighbors. Between each pair of
branches the edges form a perfect matching with matrix \(P_{ij}\),
because each vertex in one branch reaches the other branch's
first neighbor by a unique two-step path. For distinct \(i,j\),
\[
 P_{ij}+\sum_{k\ne i,j}P_{ik}P_{kj}=J_m .
\]
The summands are disjoint permutation supports. This is an exact
two-step balance: every slot is reached once, by the direct
matching or by one two-step route. The table and its identities
also suffice to reconstruct the graph: add the \(d\) parents and
the root, checking same-branch pairs through the parent and
the other cases by the matchings. The replay verifies these
identities on the Petersen graph (disjoint two-subsets of five)
and the [Hoffman–Singleton graph](https://aeb.win.tue.nl/drg/graphs/Hoffman-Singleton.html)
(five pentagons and five pentagrams with \(y=x+ij\bmod5\)
as the crossing rule). The degree-57 case remains open;
see the [recent research discussion](https://arxiv.org/abs/2606.29183).

Now make a tournament on the \(dm\) second-neighborhood vertices:
within branches use any tournament, choose any \(d\)-vertex
quotient tournament \(S\), and between \(i,j\) orient all arrows
according to \(S_{ij}\) **except** the direct matching, which is
reversed. The dense signed block is \(S_{ij}(J_m-2P_{ij})\).
This particular encoding includes local imperfections even
if the underlying matching system is perfect. Let \(G\) be
the square of its skew matrix. The exact off-branch block mean is
\[
 \frac{1}{m^2}\mathbf1^\top G_{ij}\mathbf1
   =\frac{(m-2)^2}{m}(S^2)_{ij}.
\]
Indeed, off-diagonal blocks of the skew matrix preserve the
block-constant and block-centered subspaces, and take a normalized
constant branch vector to \(S_{ij}(m-2)\) times the next one.
Only within-branch arrows can mix the two subspaces, and
the resulting two-step corrections return to the **same**
branch. So no internal orientation can repair an off-branch mean.
The formula also holds for partial search frameworks once all
cross blocks used in the computation are permutations.

For odd \(d\), all off-diagonal entries of \(S^2\) are odd.
They can all have magnitude 1 only for a switching of a doubly
regular tournament: take a nonzero vector \(x\in\ker S\).
Then
\[
 0=x^\top S^2 x
 \le -d\|x\|^2+\Big(\sum_i|x_i|\Big)^2\le0.
\]
Equality forces all \(|x_i|\) equal and
\((S^2)_{ij}=\mathrm{sgn}(x_i)\mathrm{sgn}(x_j)\).
Switching by those signs gives \(S^2=J_d-dI\) and a
regular tournament with \((d-3)/4\) common predecessors for
each pair. In particular \(d\equiv3\bmod4\).

This gives the exact optimal **maximum absolute off-branch mean**
over quotient and internal orientations, for any fixed matching
system at the indicated degree:

| Degree \(d\) | Branch size \(m\) | Best maximum block mean |
|---|---|---|
| 3, Petersen | 2 | 0 (the constant mode cancels) |
| 7, Hoffman–Singleton | 6 | \(8/3\) (Paley 7 attains it) |
| 57, conditional | 56 | \(2187/14\) (about 156.21) |

At 57 the quotient must have \(|(S^2)_{ij}|\ge3\) somewhere.
The lower bound is sharp by deleting two vertices of the Paley
tournament on 59 vertices. For its remaining signed matrix \(S\),
\[
 S^2=J_{57}-59I+uu^\top+vv^\top
\]
where \(u,v\) are the sign columns to the deleted contacts.
Thus all off-diagonal entries are \(-1,1,3\); the replay gives
406, 812, and 378 unordered pairs, respectively. Since the whole
dense encoding has even order, its off-diagonal squared entries
are even. A block of mean magnitude at least \(2187/14\) then
contains an entry of magnitude at least **158**.

This is a restriction on an encoding, not on existence of the
Moore graph. It explains why a perfect sparse Moore transport
does not imply a perfectly balanced dense tournament: the first
two examples can use doubly regular quotient signs, but a
57-vertex quotient cannot. In this encoding the mean defect is
a quotient obstruction; changing the internal arrows or the
matching labels cannot remove it. The Paley-59 deletion is a
concrete optimum for the mean calculation, not a proposed
matching completion.

## Replay

Run `node replay.js` to reproduce the strong spectra,
two exact gap crossings, square-fiber control, regular-module
inflations, Petersen and Hoffman–Singleton transport, and the
57-vertex Paley deletion. See [replay.out](replay.out).
