# One-sided periodic births at continuously fixed ports

## Exact physical scope

Fix any nonempty finite tournament T and any fixed ordered set of k distinct vertices as ports. A rightward birth with mask S adds a vertex which beats the selected old ports and loses to every other old vertex. Keep all the original ports forever; a newborn is never exposed, and no port is replaced or forgotten. Choose any nonempty periodic word of L such masks, and repeat that word. Let h_t be the actual Hamiltonian-path count after t whole periods.

**Spectral theorem.** The sequence h_t is annihilated by a polynomial all of whose roots are among `0^L,1^L,...,(k+1)^L`.

Consequently, an eventual nonconstant recurrence `h_(t+1)=b h_t+c`, with real b>1, requires `b=d^L` for an integer d in `{2,...,k+1}`.

**Ternary corollary.** If b=3, then L=1. No such system can satisfy `h_(t+1)=3h_t+4`, even eventually. A returning `3H+2` system is possible, as the independently audited three-port clock demonstrates, but its period must have length one.

All statements also hold for purely leftward births, by reversing every edge. They do not cover a period mixing rightward and leftward births, changing which physical vertices are exposed, graph rebuilding, or a state-dependent sequence of masks. The initial graph is arbitrary; it need not itself have been built by this language.

## Finite exact interface

Take unordered directed path covers of the entire old tournament with at most one component whose starting vertex is not a port. Encode each component by its start and end labels, using its port label when that endpoint is a port and a single symbol `*` otherwise. Distinct components have distinct start labels; a port cannot appear in two components, and simultaneous start/end use of one port is its singleton component. Thus a cover has at most k+1 components, and there are finitely many endpoint-connectivity shapes.

For each shape count the actual covers realizing it. The Hamiltonian-path count is the sum of the one-component coordinates. Port labels absent from the endpoints are internal vertices of their component; their identities remain part of the actual cover even though its shape does not expose them.

A rightward insertion can attach a predecessor component ending at a non-port or at an unselected port, and a successor component starting at a selected port. Either attachment may be absent. They cannot be the same component, which would create a directed cycle. Join these components through the newborn and discard covers with two non-port starts. This gives an exact nonnegative integral linear transfer M_S.

For the reverse bijection, delete the newborn from an output cover. At most two neighboring fragments are restored. A newly created start, if any, was an old out-neighbor of the newborn and hence a port. Therefore deleting the newborn cannot increase the number of non-port starts. Every output counted by the interface comes from exactly one old cover and one pair of attachments. This proves closure and exactness for all graph orders.

## Simultaneous triangularity, including the diagonal

For a shape sigma define

`Phi(sigma) = total used in/out incidences at the k ports + indicator(a non-port start exists)`.

Equivalently the first term is `2k` minus the number of port endpoint occurrences, counting a singleton port twice. The potential is an integer between zero and 2k+1.

Insertion never removes an old path-cover edge. Used incidences at a port can only increase. An old non-port start cannot be consumed as a successor, since the newborn beats only old ports. Thus the indicator of a non-port start cannot decrease either. Every transition weakly increases Phi.

If the potential does not increase, the insertion cannot attach to a port endpoint, and cannot create a new non-port start. It must therefore append the newborn after an existing non-port end, with no successor attachment. Conversely, each such choice is permitted for every mask S and leaves exactly the same endpoint shape. If `d(sigma)` is the number of non-port ends in the shape, the diagonal entry is exactly d(sigma), independent of S. There are no transitions to a different shape with equal potential.

Ordering the shapes by Phi makes all M_S simultaneously triangular, with the same diagonal values d(sigma) in `{0,1,...,k+1}`. A product for L births has diagonal d(sigma)^L. The Cayley-Hamilton theorem gives the claimed annihilating polynomial for every initial actual state and every linear observation, in particular h_t. Nilpotent behavior and possible Jordan blocks are included; no diagonalizability is assumed.

## Eventual affine recurrences and the +4 obstruction

If `h_(t+1)=b h_t+c` eventually, then on that tail

`h_t=A b^t-c/(b-1)`.

Nonconstancy means A is nonzero. Since the transfer's annihilating polynomial still annihilates every tail, b must be one of its roots. For b=3, the equation `d^L=3` forces L=1 and d=3.

For L=1 the mask never changes. Regard the first newborn as a distinguished vertex r in a fixed seed tournament. Every later newborn has the same incidences to the original graph as r, and all earlier newborns point to all later ones. Hence the growing newborn set is exactly a transitive module replacing r.

For any such fixed rooted tournament, let a_j count directed skeletons using every ordinary vertex once and r exactly j times without consecutive r's. Its clone counts are

`H_m = sum_j a_j j! S(m,j)`.

Permuting the j-1 nonempty internal fragments between occurrences of r gives a free action, so `(j-1)!` divides a_j. If H_m eventually satisfies the +4 recurrence, exponential independence forces the largest run index to be three and

`H_m = a_3 3^m - 2`.

But a_3 is even, while every tournament has an odd number of Hamiltonian paths. This is impossible. The theorem is therefore an all-order obstruction for arbitrary finite fixed port sets and arbitrary one-sided periodic words, rather than just a particular attempted gadget.

An immediate additional consequence is that for any odd prime multiplier p, a nonconstant affine recurrence in this exact language must have L=1, and its additive constant satisfies `c = p-1 (mod 2(p-1))`. Indeed the constant exponential coefficient is an integer, the leading coefficient a_p is divisible by `(p-1)!` and is even, and oddness forces the constant coefficient to be odd. This is only a necessary condition, not a realization claim.

## Verification

The argument above is symbolic for every finite k. The accompanying verifier checks every formal transition for k=0,...,4: potential monotonicity, the absence of equal-potential off-diagonal entries, and the exact mask-independent diagonal. It also checks the physical graph transfer against direct path-cover enumeration for all labelled tournaments through order four and all choices of up to three ordered ports. These are controls of the exact interface and indexing, not the source of the all-k theorem.

Run `python3 verify.py` and `python3 -O verify.py`; the receipts agree. No removable assertions are used. `interface.py` is the unchanged standalone interface implementation from the previously audited Session197 two-port packet, also supporting general k; this proof supplies the present all-k reasoning and fixed-port restriction explicitly.
