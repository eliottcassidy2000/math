# A folded path shadow retains the complete SCC skeleton, but not its context

**Session265, 11 November2026. PROVED statements below; finite controls separately labelled.**
**Chan–Pak remains OPEN.** No global reconstruction or literature-priority claim.

## Objects and inheritance

For a tournament T on n>=2 labelled vertices let H count directed Hamiltonian
paths, s(v) count their starts and e(v) their ends. The user's shadow is the
symmetric matrix W with zero diagonal and W(u,v) equal to the number of paths
using the pair consecutively. A zero-weight edge is retained as zero; omitting
it gives the same support graph. Path reversal proves W(T)=W(T*).

The nearest old mechanism is the consecutive-position convolution in
`01-canon/theorems/THM-050-consecutive-position-formula.md`, especially its total
edge-usage identity. A less-used sidecar is the split-strong-component Hamilton
cycle argument in repaired
`01-canon/theorems/THM-1900-the-insertion-response-calculus.md`, not its superseded
delete-a-path argument. The classical input is Camion's theorem: a strong
non-singleton tournament has a directed Hamilton cycle. Its source and imported
scope are recorded in THM-343-H7-impossible.md, references, item A10.
We use that input, rather than re-auditing it. Positive cycle packing below uses
THM-002-ocf.md, the inherited Redei–Berge odd-cycle formula.

The board is: folded edge counts / ordered SCCs / rank-one endpoint interfaces /
rooted even-prefix response / scalar insertion. The hostile is a labelled
six-vertex pair whose identical W hides a different scalar insertion count.

## Theorem1: tight-cut tomography

Write the uniquely ordered SCCs as C1 -> ... -> Ck. For every nonempty proper
vertex subset A define c(A)=sum_{u in A,v outside A} W(u,v). Then

    c(A)>=H,
    c(A)=H iff A is an SCC prefix or an SCC suffix.

In particular, W determines the SCC vertex partition and its order up to
reversal. It also determines whether T is strong.

Proof. A Hamiltonian path visits every SCC contiguously, in the unique SCC
order: leaving a component for a later one makes return impossible. Every
Hamiltonian path crosses any nontrivial cut at least once, giving c>=H.
A prefix or suffix is crossed exactly once. A union of whole components that
is neither prefix nor suffix has at least two transitions in its component
membership word, so every path crosses at least twice.

It remains to exclude a cut splitting a component C. Its size r is at least3.
Take a Hamilton cycle of C. Its binary A-membership word has an even number
q>=2 of transitions. Deleting a cycle edge produces a Hamilton path with at
least2 transitions: if q>=4 every deletion works; if q=2, since r>=3 there is
an edge that is not a transition, and deleting it retains both transitions.
Extend that path by arbitrary Hamilton paths in every other SCC. This is a
Hamiltonian path of T crossing the cut at least twice. All the other paths
cross at least once, so c(A)>H. This proves both directions.

Recovery: H=sum_{u<v} W(u,v)/(n-1). Declare u equivalent to v when no H-valued
cut separates them. The equivalence classes are exactly C1,...,Ck. Every two
vertices in different SCCs are separated by a prefix cut. Contract the classes
in the support graph of W. The result is the undirected path C1--...--Ck:
nonadjacent components cannot be consecutive in a Hamiltonian path, whereas
every consecutive-component pair has positive weights on ALL cross pairs.
The latter follows from Hamilton-cycle rotations, which give positive start
and end counts at every vertex of a strong component (and from the singleton
case). This path determines the order up to reversal. For k=1 it is one node.

This is algorithmic, not an exponential-cut requirement: pairwise minimum-cut
values in the nonnegative weighted graph give the equivalence relation via
lambda(u,v)>H. A standard undirected cut-tree may compress those queries, but
no cut-tree implementation is needed for the theorem or the replay.
For a strong T, the global minimum cut is strictly larger than H; for a
non-strong T it equals H and there are exactly k-1 unoriented H-valued cuts.

## Corollary: a Hamilton-cycle margin, not just strictness

Let T be strong of order n>=3 and let c_HC(T) count directed Hamilton cycles,
with cyclic rotations identified. Then EVERY nontrivial cut obeys

    c(A)>=H+(n-2)c_HC(T)>=H+n-2.

Indeed, a Hamilton cycle C has q_C>=2 cut transitions. Its n distinct deleted-
edge paths have total transitions (n-1)q_C and hence total excess above one
crossing at least n-2. Different Hamilton cycles give disjoint sets of such
paths, since a directed Hamilton path has at most one closing arc. All remaining
paths contribute nonnegative excess. Summing proves the inequality. Therefore
the folded shadow even gives the upper bound

    c_HC(T) <= floor((mincut(W)-H)/(n-2)).

This is not an exact Hamilton-cycle reconstruction. For a cut splitting an SCC
C of order r, the same proof and external SCC choices yield the local margin
c_T(A)-H(T)>=(r-2)c_HC(C) H(T)/H(C). This can certify that a nearly tight cut
cannot split a particular high-cycle core. The three-cycle attains equality.

## Theorem2: normalized boundary reconstruction

Put h_i=H(C_i), local starts s_i and ends e_i. Exact concatenation gives

    H=product_i h_i,
    W_T(u,v)=(H/h_i) W_Ci(u,v)                     [u,v in C_i],
    W_T(u,v)=H e_i(u)s_(i+1)(v)/(h_i h_(i+1))     [adjacent SCCs],
    W_T(u,v)=0                                   [nonadjacent SCCs].

Thus every boundary has total weight H and rank one, but its individual edge
weights are not generally uniform and are not generally H. For example
C3 -> C3 has H=9, each internal cyclic edge has weight6, and each of the nine
boundary edges has weight1. The sum, not each edge, has weight9.

Choose either recovered SCC direction and divide all weights by H. The row
and column sums of a boundary recover respectively e_i/h_i and s_(i+1)/h_(i+1).
Moreover local weighted degree gives

    s_i(v)/h_i + e_i(v)/h_i
      =2-sum_{u in C_i} W_T(u,v)/H.

Consequently if k>=2 the shadow recovers BOTH normalized endpoint vectors of
EVERY component: interiors have both from boundaries; at each extreme the
missing one is obtained by subtraction from the displayed degree identity.
It also recovers each local normalized shadow W_Ci/h_i. All of this is
simultaneously reversed when the global direction is reversed.

The denominator caveat matters. These rational vectors do not automatically
identify individual h_i; H gives their product, and normalized count data can
cancel common factors. We do not assert scalar factor recovery or internal
arc reconstruction from these data.

## Theorem3: one ordinary observer opens the converse fold

Append a new source z, so z -> every vertex of T. Every new Hamiltonian path is
z followed by an old one. Therefore

    H(z->T)=H(T), W_new|_T=W_T, W_new(z,v)=s_T(v).

A sink similarly reveals e_T. If T is non-strong, the marked source's positive
weights identify its first SCC and orient the entire recovered SCC line. The
source does not alter H. This is an exact endpoint measurement, not a claim
that an individual endpoint count is itself a scalar tournament H.

The measurement explains what an arithmetic constructor must preserve: an
undirected shadow can suffice for SCC concatenation after normalization, but
not for arbitrary rooted insertions.

## Theorem4 / FINITE-EXACT witness: a scalar context separates equal shadows

Use vertex set0,...,5, lexicographic pairs i<j, and bit1 to mean i->j. Let
A have mask1332 and B mask2744. Their full adjacency lists are in RESULTS.json.
Both are strong and have H=45. Their common W, in lexicographic pair order, is

    (13,13,13,13,23,13,23,13,13,13,23,13,13,13,13).

Their start vectors are respectively

    (6,9,6,6,9,9),  (6,9,9,6,6,9),

and their end vectors are 15 minus these coordinatewise. The SAME marked
vertex x=2 therefore has different start counts6 and9. These are distinct
rooted tournaments even if an unrooted relabelling happens to identify them.
The witness claims exact labelled shadow equality, not distinct unlabelled
converse classes.

Now attach z with sole outneighbor x, and with every other old vertex pointing
to z. The resulting scalar counts are

    H(A plus z)=93,    H(B plus z)=111.

Hence the shadow map, even retaining vertex labels and the marked vertex, is
NOT a congruence for tournament insertion. Its hidden arithmetic defect is18.

Here is a general positive explanation, not merely a search witness. Let
b_x(S) be the number of Hamilton paths of T[S] starting at x, taken as zero
when x is absent, and let H(empty)=1. Define

    K_x(T)=sum_{S containing x, |S| even} b_x(S) H(T-S).

Then the above apex operation obeys H_new=H+2K_x. By odd-cycle packing, a new
odd cycle must contain z, hence starts z->x and then follows an even-cardinality
path on S. Its last vertex is not x and therefore points to z. Every such
path closes uniquely; the remaining compatible odd-cycle packing contributes
H(T-S). At most one cycle uses z, and its cycle weight is2. This proves the
formula. It is the existing Session264 fusion carrier's coordinate

    K_x=(U-d)/2,  H_new=H+U-d,

where U is the complete rooted-start subset response and d its odd-minus-even
response. Thus no redundant new state is required: the old rooted carrier
restores exactly what W discarded in this context. For the two witnesses
K_x=24 and33. RESULTS.json splits these into |S|=2,4,6 contributions.

## Constructive consequences and boundaries

1. A proposed scalar constructor can be compiled first at the SCC skeleton
level. The shadow supplies a genuine certificate of which interfaces are
one-pass and which must handle recirculation.
2. Rank-one boundary margins give normalized start/end data for composing
existing SCCs. This is a concrete sufficient interface, not a visual analogy.
3. A search that quotients by W may discard useful arithmetic choices: the
same state has two legal rooted outputs93 and111. Keep (U,d), or a stronger
existing rooted response carrier, before searching insertion programs.
4. SCC detection supplies a clean strong-core reduction for inverse spectra:
prime H has at most one nontrivial SCC, and general H factors over the recovered
cores. It does not produce cores for7,21,63 or prove all large odd targets.
5. The user's path-reversal intuition is exactly right. The retained information
is stronger than a generic folded invariant: SCC topology survives completely.
The lost information remains arithmetically important inside strong components.

## Controls, source scope and reproduction

Run `python3 05-knowledge/results/session265_converse_dyadic_20261111/shadows/replay.py`.
It checks all labelled tournaments at orders2..5;12 seeded tournaments at each
order6..9; and both explicit six-vertex witnesses. Each test checks converse,
degrees, every cut, exact SCC reconstruction, local/boundary factorization, and
all single-mark apex formulas. Witness scalar counts are also independently
counted by enumerating all permutations, including all5040 permutations for
each seven-vertex apex extension. The general proofs above, not these finite
controls, supply the all-order statements.

No source claim from the supplied pasted texts is used as a proof dependency.
Their block/folding themes inspired the actual tight-cut theorem and the
contextual loss test. The old consecutive-position formula and repaired
insertion theorem are credited as nearest mechanisms; the (U,d) response is
explicitly inherited from Session264, `fusion/PROOF.md`.
