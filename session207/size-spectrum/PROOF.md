# Session 207: exact order spectra and the state needed to compute them

Status: FINITE-EXACT census through order 8; PROVED general recurrences below.
No claim of a closed-form spectrum at every order or of global odd-spectrum
completeness is made. No novelty claim is made for the previously recorded
small-order spectra. This packet independently regenerates them without
importing earlier representative files.

## Inheritance and scope

Closest inherited mechanism: strong-component multiplicativity and source/sink
padding (THM-1370, `01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md`).
Hostile example: 7 and 21 are permanently absent, while 63 is absent at order 7
but present at order 8. Corrected near miss: a gap below the maximum need not
be permanent. Underused sidecar: the full family of induced-subset endpoint
counts, rather than just the scalar H. Concept board: order support, strong
components, padding, vertex-extension profiles, metagraph quotient loss.
Anchor: exact order spectra. Niche: SCC product recurrence. Wildcard: whether
one-vertex continuation support descends to the count H (it does not).

Let H(T) count directed Hamiltonian paths as ordered lists of all vertices.
Only nonempty tournaments are considered, except for algebraic empty-state
conventions. Write S_n for their count support, and C_n for the support among
strong tournaments. A singleton is strong, so C_1={1}; C_2 is empty.

## Exact small-order answer

All displayed intervals below are odd integers only.

- S_1=S_2={1}.
- S_3={1,3}.
- S_4={1,3,5}.
- S_5=[1,15] minus {7}.
- S_6=[1,45] minus {7,21,35,39}.
- S_7=[1,189] minus {7,21,63,107,119,149,161,163,165,167,169,173,177,179,181,183,185,187}.
- S_8=[1,661] minus {7,21,611,615,617,619,623,625,635,647,655}.

Support cardinalities are 1,1,2,3,7,19,77,320. `census.json` also includes
all strong supports and one adjacency-code witness for every count at every
order (430 witnesses). Strong minima at orders 3..8 are 3,5,9,15,25,45.
These are finite exact results, not extrapolations.

### Why the census is exhaustive

Every tournament of order n has an induced tournament of order n-1. By induction
choose an isomorphism from that deletion to one representative already stored;
the deleted vertex has one of the 2^(n-1) possible incident orientation patterns.
Thus augmenting all stored representatives covers every isomorphism class.
For n<=7 representatives are deduplicated by an actual canonical form: sort the
vertices by outdegree, and minimize the adjacency code over every permutation
inside each equal-outdegree block. Outdegrees are invariant under isomorphism,
so this is constant on, and distinguishes, isomorphism classes. At n=8 no
deduplication is necessary: all 456*128=58368 extensions are evaluated. Duplicate
classes affect candidate multiplicities, not support. Candidate multiplicities
are NOT labelled or unlabelled tournament frequencies.

The DP builds paths by appending an unused out-neighbor. For n<=5 each evaluated
candidate is independently checked by all-permutation enumeration. The separate
Python checker uses the reverse endpoint recurrence and checks every published
witness. Class counts through n7 are 1,1,2,4,12,56,456, checked as secondary
controls but not needed for the induction proof of completeness. All generated
counts are checked odd. The singleton, transitive orientation and cyclic triangle
occur naturally as low-order controls.

Reproduction:

    g++ -O3 -std=c++17 verify.cpp -o /tmp/session207-spectrum
    /tmp/session207-spectrum > census.json
    python3 check_witnesses.py

The code supports precisely n<=8; its fixed arrays are not an all-order program.

## General exact recurrences

### Strong-component convolution

Put S_0={1} as an empty-product convention. For every n>=1,

    S_n = union over k=1..n of { a*b : a in C_k, b in S_(n-k) }.

Proof: the strongly connected components of a tournament have a unique linear
order with every arc pointing forward. A Hamiltonian path traverses each
component consecutively, in that order; restriction and concatenation give a
bijection with a tuple of component Hamiltonian paths. Its count is their
product. Taking the first component gives one inclusion in the recurrence.
Conversely an ordered join of a strong k-vertex tournament and any tournament
of order n-k realizes the stated product. This is an exact reduction to the
strong spectra, not a solution or closed formula for those strong spectra.

### Vertex extension with the necessary sidecar

For each induced subset A, let e_T(A,v) count Hamiltonian paths in T[A] ending
at v and b_T(A,v) count those starting at v. They obey the usual subset
recurrences, initialized to one on singletons. Given a new vertex x, put
I={v:v->x} and O={v:x->v}, a partition of V(T). Define

    L_I(empty)=R_O(empty)=1,
    L_I(A)=sum over v in A intersect I of e_T(A,v) for nonempty A,
    R_O(B)=sum over v in B intersect O of b_T(B,v) for nonempty B.

Then exactly

    H(T+x) = sum over A subset V(T) L_I(A) R_O(V(T) minus A).

Proof: cut a Hamiltonian path at x. A is the uniquely determined set before x;
the prefix ends in I, the suffix starts in O, and there is no further relation
between these two subpaths. Concatenation is the inverse map. Empty prefix and
suffix use the stated conventions. Augmenting every tournament with every
partition I,O gives an exact all-order procedure; retaining just H(T) does not.

### Scalar continuation is genuinely insufficient

The two five-vertex tournaments with adjacency codes 2 and 5 both have H=3.
Their sets of counts obtainable by adding one vertex are respectively

    {3,5,9,13,15,17,23,25,27,33},
    {3,5,9,13,15,27,29}.

Here bit k refers to the kth lexicographic pair i<j and equals one iff i->j.
`check_witnesses.py` exhausts all 32 extensions of both examples, using its
independent endpoint recurrence, and writes the exact counterexample JSON.
Structurally these are ordered joins `1 -> 1 -> C3` (code 2) and
`1 -> C3 -> 1` (code 5). They have the same multiset of strong components,
but place the triangle at different positions. Thus even the unordered strong-
component inventory loses continuation information: the current product is
commutative, but the extension geometry remembers component order.

Thus the map T -> (|T|,H(T)) preserves current count but destroys extension
support. The induced-subset endpoint profile is a sufficient restoration
sidecar. This does not assert it is a minimal sidecar.

## Relation to metagraphs

On the graded vertex-extension graph of tournament isomorphism classes, decorate
a vertex T by (|T|,H(T)). The row support is S_n. Source/sink padding provides a
vertical edge preserving H, hence S_n is contained in S_(n+1). The scalar
continuation counterexample shows that contracting each equal-H level is not a
lossless predictor of outgoing decorated edges: different vertices in the same
fiber have different one-step supports. The support quotient records that some
extension exists, not that every tournament in a count fiber admits it.

This graded graph is distinct from the repository's native fixed-order
metagraph, whose move keeps a directed Hamiltonian path and reverses every
other arc. Native-metagraph topology and fibers require that edge definition
and its sidecars; no claim that count support determines that graph is made
here.
