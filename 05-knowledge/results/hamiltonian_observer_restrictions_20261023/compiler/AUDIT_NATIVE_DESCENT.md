# Independent audit of native unmarked incidence descent

23 October 2026. **PASS: PROVED within the declared unmarked incidence
carrier.** Independently reviewed
[`NATIVE_CONNECTEDNESS.md`](../coupling/NATIVE_CONNECTEDNESS.md),
[`check_native_descent.py`](../coupling/check_native_descent.py), and the
definitions and normalization in
[`RESTRICTIONS.md`](../coupling/RESTRICTIONS.md). This audit addresses the
unbounded induction, not only the finite controls.

## Carrier and the isolated component

The compression is genuinely bijective: compressed pair `{i,j}` stores
the direction of native chord `(i,j+1)`, while all native consecutive arcs
are fixed forward. Any nonzero word contains a backward chord completing
a directed cycle with the intervening forward path. Thus its native
tournament is nontransitive. Conversely a transitive native tournament
has a unique Hamiltonian vertex order and yields only the zero word.

Every finite tournament occurs: a Hamiltonian path can be built by
inserting a new vertex immediately before the first path vertex it beats,
or at the end if none. The preceding path vertex, if present, beats the
new vertex by the choice of first position. Hence the target side of the
incidence graph includes every tournament class, not just a tested subset.

Ordinary compressed graph relabeling is a zero-cost change of representative
of a source class. It can change the native tournament. Reordering a
genuine Hamiltonian path inside that native tournament then witnesses one
source-to-source incidence adjacency. This distinction is essential and is
explicitly retained in the theorem.

## Strong-four reduction and its exact boundary

With compressed positions `(a,u,b,w)`, the forced compressed edges `01,12`
produce native backward arcs `02,13`. Let epsilon be compressed edge `02`,
which is native chord `03`.

- For epsilon zero, native prefix `(2,0,3,1)` is Hamiltonian. Its only
  backward internal pair is new positions `0,3`, namely old arc `1->2`.
- For epsilon one, prefix `(2,3,0,1)` is Hamiltonian and has the same one
  backward internal pair.

Both prefixes end at old native vertex one. The next vertex is four, and
the junction `1->4` is forward exactly when compressed edge `13={u,w}`
is absent. All prefix vertices remain before all suffix vertices, so the
*total* number of backward crossing arcs is unchanged. This is enough for
the edge-count drop only after the actual exit has separately been proved.
The internal counts drop from two or three to one, respectively.

A minimal boundary hostile is compressed `E={01,12,13}` on four vertices.
The unguarded path-type prefix is internally legal, but the proposed native
order `(2,0,3,1,4)` fails because `4->1`, not `1->4`. This is a spanning
star and is handled by the distinct star repair. It cannot invalidate the
guarded reduction because its center has no nonneighbor.

## Exhaustion of cases and the induction

For `m>=4`, if no degree lies between two and `m-2`, degrees lie in
`{0,1,m-1}`. With no universal vertex this is a matching. With one it is a
spanning star, since every other vertex already has degree at least one
and cannot have a second neighbor. With two, every remaining vertex has
degree at least two and must also be universal, so the graph is complete.
For `m=2,3`, the same three shapes exhaust the graphs directly.

The matching preparation is exact, including its boundary. Place two
matching edges at compressed `01,23`. Under native order
`(1,2,0,3,4,...)`, these become compressed `01,13`, and all remaining
matching edges keep their positions. The exit `0->3` is forward because
compressed `02` is absent. The weight stays constant and compressed vertex
one has degree exactly two. Since a two-edge matching needs `m>=4`, this
is an intermediate degree, so the very next step admits the guarded
strict reduction. The preparation cannot repeat forever.

For the spanning star, center-first compression yields a transitive native
suffix `1,...,n-1`, all vertices at least two beating zero, and `0->1`.
Rotating zero to the end is Hamiltonian with exactly that last arc
backward. For a complete compressed graph, consecutive increasing pairs
in decreasing block order are Hamiltonian. The only backward arcs are
the original consecutive pairs split between blocks, numbering
`ceil(n/2)-1=floor((n-1)/2)`. This is positive and strictly below
`binom(n-1,2)` for `n>=4`; order three is already the one-edge terminal.

Therefore a strict positive edge-count reduction costs at most two
source-adjacency moves. Induction gives distance at most `2(e-1)` from
every nonempty source to the one-edge class. Matching preparation changes
the source shape, and subsequent steps lower edge count; the chosen paths
can be taken simple. Combined with the isolated zero/transitive pair and
the existence of Hamiltonian presentations, this proves exactly two
incidence components for every `n>=3`.

## Conductance and spectral normalization

On the nonzero word space of size `Z=2^M-1`, let `W` have integer incidence
counts and let its row and column sums be `r_i,c_t`. The alternating chain
is `P=D_r^(-1) W D_c^(-1) W^t`, with stationary mass `pi_i=r_i/Z`.
For distinct adjacent source states,

    pi_i P_ij = (1/Z) sum_t W_it W_jt / c_t >= 1/Z^2.

At least one product is a positive integer and every `c_t<=Z`. The
Dirichlet form is one half the ordered-pair sum, equivalently the sum over
unordered pairs with conductance `pi_i P_ij`; hence there is no extra
factor of two in the manuscript's comparison.

For simple descent paths of length at most `R=2(M-1)`, first center at the
base state, apply Cauchy--Schwarz on each path, and sum using `sum_i pi_i=1`.
Each used undirected edge receives total path weight at most one. This
gives `Var_pi(f)<=R Z^2 E_P(f)` and precisely the stated gap lower bound

    1 / [2(M-1)(2^M-1)^2],  n>=4.

The chain is positive semidefinite because it is an alternating
conditional-expectation coupling. The Poincare gap convention therefore
also introduces no negative-eigenvalue qualification. This is a weak
order-dependent bound, not the auxiliary whole-cube gap `8/9`.

Complement-pairing source classes identifies the empty and complete
classes in the two components, so the paired incidence is connected.
Further target converse-merging cannot disconnect it.

## Computation and limits

Independent replay of both commands passed **30,882 explicit checks**:

    python -B 05-knowledge/results/hamiltonian_observer_restrictions_20261023/coupling/check_native_descent.py
    python -O -B 05-knowledge/results/hamiltonian_observer_restrictions_20261023/coupling/check_native_descent.py

The script covers every nonzero labelled compressed graph for native
orders three through six: `1,7,63,1023` words. It checks every actual
Hamiltonian junction, recompression, graph-class preservation under
compressed relabeling, and tournament-class preservation under the native
path change, using full permutation canonicalization. Maximum observed
lengths are `0,1,3,7`. These are meaningful controls of the local rules and
exception coverage; the preceding induction supplies the all-size claim.

No correction was required by this audit. The conclusion does **not**
assert connectedness while retaining one fixed native tournament across
the whole descent, indexed owner pins, attachment ports, or native
electrical mod-eight guards. It is also distinct from the compiler's
Hamming connectivity inside a single automorphism orbit, whose entire
ordered arc word and native isomorphism owner remain fixed.
