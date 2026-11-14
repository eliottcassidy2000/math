# To leave the factor 3, one nonuniform Paley contact is not enough

The [7-orbit module calculation](SEVEN_ORBIT.md) shows that a
homogeneous Paley 7-module forces the full Hamiltonian path count
to be divisible by 63. Here a **sharper robustness theorem**:
allow one external vertex arbitrary nonuniform arcs to the Paley
orbit, with all other vertices still uniform to it (and otherwise
completely arbitrary). The resulting ordinary count is *still*
divisible by 3. Two external vertices suffice to break this
obstruction for both ordinary prime and \(7p\) outputs.

## The apex profile depends only on the number of its exceptional arcs

Let \(G=T-v\), where the Paley set \(U\) is a module in \(G\).
Choose directions from an apex \(v\) to vertices of \(G\) using
\(x_i=1\) for \(i\to v\). Outside \(U\) all apex edges are
fixed and need not have any uniformity. For distinct \(i,j\)
let \(E_{ij}\) count ordered two-path covers of \(G\) having
first path end \(i\), second path start \(j\) (both nonempty).

Deleting \(v\) from a full Hamiltonian path leaves one of three
possibilities: \(v\) was first, last, or between two paths.
Thus the exact count has linear terms from the end positions
and the quadratic insertion term
\(\sum_{i\ne j}E_{ij}x_i(1-x_j)\).
This insertion identity is inherited; what follows is its
restriction to the Paley orbit.

Translations and square multipliers of \(\mathbb F_7\) are
automorphisms of \(U\) and extend by fixing \(G\setminus U\).
They act transitively on vertices and on unordered pairs,
since \(-1\) is not a square. On ordered pairs they have two
classes (the internal arc orientation). Hence the coefficients
of all linear terms on \(U\), including cross terms with outside
vertices, are the same. Also the symmetrized coefficient
\(\kappa=E_{ij}+E_{ji}\) is constant for \(i,j\in U,\ i\ne j\).
For a marked subset of size \(t=\sum_{i\in U}x_i\), the regular
score 3 of every Paley vertex makes the two directed cut counts
equal: each is \(t(7-t)/2\). Therefore
\[
                  H_t=H_0+\eta t+\kappa\,t(7-t)/2
\]
with integral \(\eta\), independent of the actual set of marked
vertices. This equality holds even though the distinct marked
sets need not be isomorphic. Taking \(t=0\) and 7 preserves the
whole homogeneous Paley module in \(T\). Both endpoint counts
are divisible by 63 by the module calculation, so \(9\mid\eta\).

It remains to count the total of \(E_{ij}\) over ordered distinct
pairs in \(U\). For such a two-path cover, the two vertices
beside its distinguished break both belong to \(U\), so they
belong to distinct occupied \(U\)-slots, even if the slots are
adjacent after concatenating the two paths. Specify all outside
vertices in order, the break, and the positions of the occupied
\(U\)-slots. Since the module is uniform, valid outside-slot
specifications depend on no internal choices. Their internal
weight with \(r\ge2\) slots is \(B_r(U)\).
Every \(B_r(U)\) for \(r\ge2\) is divisible by 126.
Hence \(126\mid\sum_{i\ne j\in U} E_{ij}=21\kappa\), so
\(6\mid\kappa\). Substitution in the displayed formula shows
\(3\mid H_t\) for every possible contact cut.

## Two nonuniform vertices really do suffice

Label the Paley orbit by \(0,\ldots,6\), use differences \(1,2,4\),
and add \(v=7,w=8\) with \(v\to w\). The bit masks below
record which orbit sites each of \(v,w\) beats.
\[
\begin{array}{c|r|r|r}
v\text{ mask}&w\text{ mask}&H&F_2\\\hline
96&115&2621\ \text{(prime)}&41140\\
36&63&2471=7\cdot353&38668
\end{array}
\]
Both quotient graphs are strong: each added vertex has both
incoming and outgoing arcs to the strong orbit. Counts and
primality are independently replayed from the literal arcs.

There is also an unbounded **two-contact** \(7p\) rail, neither
a homogeneous Paley module nor an order-7 symmetry argument.
Use the second nine-vertex core as a whole module in the cyclic
triangle with \(TT_r\) and a singleton. Outsiders beyond the two
special contacts are uniform to the original Paley orbit.
From the two-module count formula, the exact cofactor is
\[
                            5171\cdot2^r-4465,\quad r\ge1.
\]
Its only possible padding exponent is decoded by checking
\((p+4465)/5171\) for a positive power of two at least 2.
At \(r=6\) the cofactor is the prime 326479. This does not
establish arbitrary prime cofactors. It proves that two contacts
are the exact finite threshold for escaping the forced factor
3, and provides an explicit family with a real target inverse.

This cofactor rail is not merely another way to package an
order-7 automorphism. In the nine-vertex core no score occurs
seven times. In the rail, its nine scores shift up by \(r\),
and are all above the \(1,\ldots,r\) scores of the transitive
arm. The extra singleton can increase the multiplicity of
at most one class by one; it still stays below 7. Thus when
the cofactor is prime at least 11, the symmetry sieve of
[SEVEN_ORBIT.md](SEVEN_ORBIT.md) forces the whole graph to be
rigid. Both routes to the factor 7 (a genuine free orbit and an
arithmetic cancellation with no orbit) are physically available.

The two-contact profile also shows the cost of retaining the
next channel: for this core \((F_1,F_2,F_3)\equiv(0,0,1)\pmod7\).
One plain outer wrap therefore preserves 7-divisibility of \(H\);
the next wrap loses it, since the intermediate \(F_2\) is
\(4F_2+2F_3-2F_1\equiv2\pmod7\). In contrast every \(F_m\)
of a genuine 7-action is divisible by 7 (the group freely
permutes the fully ordered lists of colored paths), so such
a failure cannot occur on the natural symmetric extension.
