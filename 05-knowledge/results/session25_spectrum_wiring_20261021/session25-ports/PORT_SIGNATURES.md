# Exact port signatures, and why S3 is not an anonymous wiring quotient

## Scope and inherited mechanisms

PROVED elementary separator results, with FINITE-EXACT controls. These are
applications of standard boundary-colouring / connection-matrix ideas, not
claims of a new general theory. The new session connection is explicit:
state compression must preserve every allowed continuation, just as the
attached Busy Beaver paper preserves phase and the tournament wiring packet
retains exit/entry ports. No new plane chromatic-number bound, universal
Collatz claim, or full tournament classification follows.

Primary-source context: L. Lovász, *The rank of connection matrices and
the dimension of graph algebras*, arXiv:math/0408232v1, Theorem2.2 and
Lemma2.4 (pages5–6 in the PDF), gives the general automorphism-orbit rank
principle. Source checked directly: https://arxiv.org/pdf/math/0408232 .
The target complete graph K_k is twin-free; its colour-permutation orbits
on r-tuples are precisely P_r(k). Below is an independent elementary
proof of this special case, including boundary-only distinguishing tests.
No novelty is claimed for the general rank principle.

## 1. A complete counting signature for a named boundary

Fix integers k>=1 and r>=0. Let G be a finite loopless graph with r distinct named
boundary vertices 1,...,r. All colours are named and interchangeable; no
lists, precolouring, or colour-dependent weights are imposed. Let P_r(k)
be the set of partitions of the boundary into at most k nonempty blocks.
For pi in P_r(k), fix one assignment of distinct colours to its blocks,
and let f_G(pi) count its proper-colouring extensions to all of G.
Colour permutation proves independence from the chosen block colours.

If G and K share exactly this boundary, and have disjoint interiors, then

    Col_k(G union_boundary K)
      = sum_pi (k)_{|pi|} f_G(pi) f_K(pi).

Here (k)_b=k(k-1)...(k-b+1), and the empty-boundary convention is (k)_0=1.
The formula is exact even when boundary edges occur in both fragments:
restriction gives a pair of extensions and conversely each compatible pair
glues uniquely. A boundary violation makes the appropriate factor zero.

Thus the full partition-indexed vector is sufficient for every further
boundary gluing, including feasibility, which is positivity of the sum.
It discards internal geometry, vertex identities, and detailed colourings,
but none of these can change this particular abstract counting operation.
Geometry and embedding compatibility remain separate sidecars.

## 2. Necessity, not merely sufficiency

For fixed k,r, two fragments give the same glued colouring count against
EVERY graph on the named boundary if and only if their signature vectors
are equal. Consequently allowing arbitrary larger attachments cannot make
a smaller counting equivalence possible.

Proof. For a boundary edge set E let u_E(pi) be 1 when the partition pi
separates both ends of every edge in E, otherwise 0. These are exactly the
signatures of boundary-only graphs. Their linear span is all functions
on P_r(k). Indeed write d_ij(pi)=1 if i,j are in different blocks, else0.
To single out a partition sigma, multiply (1-d_ij) for all pairs within
its blocks and d_ij for all pairs in different blocks. This product is
the indicator of sigma. Expanding the factors1-d_ij expresses it as a
signed sum of square-free monomials in d_ij, each a u_E. Therefore equality
of all glued counts implies (k)_{|sigma|}(f_G(sigma)-f_H(sigma))=0 for every
sigma. Its prefactor is positive, proving necessity. Sufficiency is §1.
This is a spanning argument; it does not pretend that a signed linear
combination itself is a graph or a physical gadget.

In the corresponding universal boundary-state problem, the smallest exact
state partition has |P_r(k)| states: two colour assignments are equivalent
exactly when they have the same equality partition. Distinct partitions
are separated by a test adding one named edge. For k>=r this is Bell(r),
not r!, and for smaller k it is the truncated partition count. This is a
statement about distinct states, not a lower bound on bits or storage for
arbitrarily encoded vectors. Restricted physical attachment families may
have coarser equivalences; no geometric minimality is claimed.

## 3. The first S3 loss, versus harmless S2 port exchange

For two ports the two equality partitions are fixed by the port-exchange
S2 action. For three ports and k>=3 there are five equality partitions:
one all-equal, three specifying WHICH pair is equal, and one all-distinct.
S3 permutes the three pair-equality states and fixes the other two. Thus
colour-name quotienting and port-name quotienting are DIFFERENT operations.
For k=2 the all-distinct state is absent; the three pair states still differ.

Minimal hostile for k>=2: G has three boundary vertices and only edge12;
H has the same named boundary and only edge13. They are port-permutation
isomorphic and have identical anonymous type data. Attach K, the single
boundary edge12. Then G union K still has one edge, while H union K has
edges12 and13. The counts are respectively

    k^2(k-1)    and    k(k-1)^2,

whose difference is k(k-1)>0. Port names cannot be discarded while the
attachment port names remain fixed. This is not a contradiction to ordinary
graph isomorphism: relabelling G must relabel its interface as well.

The exact survivor is simultaneous equivariance: permuting the boundary of
EVERY operand in a gluing leaves the count unchanged. An independent port
permutation on one operand requires a permutation sidecar. This is the same
logical distinction as ordinary versus signed C3 symmetry and synchronized
signs between tournament modules; it is not an identification of those
objects. An S3 action on ports is not automatically a dynamical symmetry.

## 4. Six-colour forcing and metagraph consequences

The Session25 two-port support theorem reduces separating-fragment
assemblies exactly to their hard-inequality skeleton. Going to three or
more shared ports is a genuine way to escape that interface, but demands
the full named partition relation, rather than merely equal/different
counts or anonymous orbit sizes. At six colours a three-port counting
certificate has five entries; its feasibility version has five Boolean
entries. Counts retain more than support.

For tournament substitution the corresponding boundary is ordered entry
and exit vertices of path runs, not a colouring partition. The transferable
statement is the continuation test: a quotient is adequate only when all
allowed continuations agree. The agreement-graph/port packet supplies the
actual tournament-specific realization. No map from graph colouring counts
to Hamiltonian counts is claimed here.

For the paper's A.9 recursion, phase-forgetting fails the one-step
continuation test, while encoding phase as integer parity passes. This
provides an explicit common proof obligation, not a proof of Collatz
termination or a six-colour obstruction in the Euclidean plane.

## Verification

Run python3 verify_ports.py and python3 -O verify_ports.py. The checker
constructs all colour-equality partitions, reports their number, computes
exact ranks of the boundary-edge test matrix, checks the gluing formula by
full named-colouring enumeration, and checks the S3 hostile for palettes2..8.
