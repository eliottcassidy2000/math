# Contextual response interfaces: the difference between losing data and safely forgetting it

**Status: PROVED elementary framework and exact Hamilton-path realization; FINITE-EXACT controls.** The constituent mechanisms are elementary invariant-subspace observability and nilpotent matrix elimination. No literature-priority claim is made. The new research proposal is to use their common interface contract systematically: an observation is safe under composition only after closing its blind subspace under every permitted continuation. The topology obtained is explicitly relative to that response, not automatically the original topology.

## 1. Inheritance and concept board

S170 separates three facts that should remain separate: a profile can be constant on a physical flip family; nontrivial transport can remain on that family; and a choice of filled cells supplies relations not determined by the graph. Its holonomy packet further shows that promised subgroups admit much smaller separating probe banks than the whole transport group.

Relevant inherited sources are the S170 holonomy, monodromy, and meta packets as routed by [the S170 index](https://github.com/eliottcassidy2000/math/blob/2d137c00eaf36a6fbe9f06044ac56d2c3bf7c7d6/00-navigation/S170.md). Only those explicitly read packets are used here. The all-sheet obstruction, rather than a one-sheet return, is the closest proved mechanism. The hostile example is a presently invisible defect exposed by a later continuation. The least-used sidecar is the **allowed continuation algebra**: the same observer may be sufficient under one operation set and insufficient under another.

Five concepts remain live:

1. Current observation versus future contextual response.
2. Safe quotient versus arbitrary compression.
3. Observable holonomy versus actual holonomy.
4. Square-zero path memory versus ordinary walk memory.
5. Minimal observer completion versus minimal state dimension.

## 2. Definitions: response radical and observation horizon

Let k be a field, V a d-dimensional vector space, A_1,...,A_m linear operators on V, and C:V→k^p a linear observation. The operator alphabet is part of the specification. An allowed continuation here is every finite word in that alphabet, including the empty word. If applications restrict allowed words using a finite-state controller, incorporate that controller state into V first. The resulting finite-dimensional system has its own dimension and corresponding horizon bound. Arbitrary language restrictions can require infinite memory and are not covered by this finite-dimensional reduction; do not replace the admissible language silently by all words.

Write A_w for the corresponding composition. Define the **contextual response map**

    R(v)=(C A_w v)_w,

and its **response radical**

    N=ker R=intersection_w ker(C A_w).

Thus v and v' have the same response iff every future context gives the same observed output. Define

    N_0=ker C,
    N_(j+1)=N_0 ∩ intersection_i A_i^(-1)(N_j).

Equivalently N_j is the set invisible to words of length at most j. The **observation horizon** is the smallest j with N_j=N. The response interface is Q=V/N. These names describe elementary constructions, not asserted new mathematical objects without precedents.

## 3. Universal safe-forgetting theorem

**Theorem 1.** Put r=rank C. Then:

1. N is the largest common A_i-invariant subspace contained in ker C.
2. N=N_(d-r), and stabilization at any earlier step is permanent.
3. Q carries unique induced operators and observation making V→Q intertwine every operation and output.
4. If a surjective linear compression q:V→W supports operators B_i and observation D with q A_i=B_i q and C=D q, there is a unique surjection W→Q making V→W→Q equal V→Q. Consequently dim W≥dim Q, and Q is the unique minimal such interface up to the evident isomorphism.
5. Two states have unequal contextual responses iff some word of length at most d-r distinguishes them.

**Proof.** For v∈N and each i, every C A_w A_i v is zero; hence A_i v∈N. Conversely an invariant subspace inside ker C is killed after every word, so is contained in N. This proves (1).

The recursion gives the finite-word characterization by induction. Its subspaces descend from the (d-r)-dimensional N_0. Equality N_j=N_(j+1) makes N_j invariant and hence equal N by (1). Before equality each step reduces dimension by at least one. Thus by step d-r the chain has stabilized or reached zero, proving (2). Quotient invariance gives (3). For (4), ker q is invariant and contained in ker C, hence contained in N. Therefore q(v)↦v+N is well defined, unique and onto. Finally (5) follows from N=N_(d-r). ∎

The horizon bound is sharp: take C(v)=v_0 and A e_j=e_(j-1) for j≥1, A e_0=0. The state e_(d-1) first appears after d-1 steps. More generally selecting the first r coordinates gives horizon d-r for r≥1.

The theorem concerns a **linear compositional interface**, not the minimum number of bits needed by an arbitrary nonlinear representation. It preserves the specified observation in every specified continuation, not every property of the original object.

### The smallest warning

Take V=k², C(x,y)=x, and A(x,y)=(y,x). The states 0 and e_2 currently look identical, but their outputs after A differ. Quotienting by ker C does not define an A-action. Keeping only the current scalar output is therefore not an interface under the declared operation.

Conversely if every A_i preserves ker C, the present observation already is safe: N=ker C and no closure is needed. Lossy does not mean invalid; **unstable loss** is the obstruction.

## 4. Holonomy, relations, and response-flatness

Suppose now the A_i are invertible and include their inverses. Let F be the free group on the named generators, with representation ρ:F→GL(V). The induced representation on Q is denoted ρ_Q. A one-vertex cell complex is specified by attaching two-cells along words r_α. Its fundamental group is F modulo the normal closure of those words.

Call a relation **response-flat** when all contextual responses are unchanged by its holonomy.

**Theorem 2 (exact response-flatness criterion).** For a word r with h=ρ(r), the following are equivalent:

(a) R(hv)=R(v) for every v.
(b) (h-I)V⊆N.
(c) ρ_Q(r)=I.
(d) C A_w(h-I)=0 for every word w of length at most d-rank C.

Consequently the response local system extends across all the chosen two-cells iff these conditions hold for every attaching word. The original local system extends iff every h itself is identity.

**Proof.** The definition of R gives (a)⇔(b); the quotient gives (b)⇔(c); Theorem 1 gives (b)⇔(d). The induced representation factors through the presented group exactly when it kills every relator, since a kernel is normal. This is precisely the linear local-system extension criterion. ∎

This is an actual source-to-target map: the original transport representation maps to its induced response representation. Its destroyed information is the action on, and displacement into, N. The retained sidecars are the operation alphabet, observer, and attaching words. There is no conclusion that a response-flat loop is nullhomotopic in the original complex.

**Hostile curvature example.** Over k of characteristic different from two, let C=(1,0), h=diag(1,-1), and A swap the two coordinates. Although C h=C, applying A after the holonomy detects its effect on e_2: C A h e_2=-1 but C A e_2=1. For alphabet {A,h} and inverses, N=0. The relation h=1 is presently invisible but not response-flat. If A is forbidden and the alphabet is only {h}, then N=span(e_2), and the same relation really is response-flat. **Visibility is relative to future experiments, not merely the present statistic.**

For S170, a profile that remains constant throughout an allowed family can be a valid but extremely coarse response. It need not recover signed-sheet monodromy or determine which higher cells were selected. Conversely a claimed profile interface for a larger family of contexts needs a proof of closure; equality of seed profiles alone does not supply it.

## 5. An exact tournament realization: square-zero boundary response

Let G be a finite directed graph on labelled vertices 1,...,n, n≥1, with adjacency weights a_uv in k. Tournaments are the case a_uv∈{0,1}, a_uu=0 and a_uv+a_vu=1 for u≠v. Work over

    R=k[z_1,...,z_n]/(z_1²,...,z_n²).

This ring has basis z_S=product_(v∈S) z_v indexed by subsets. It retains exactly the visited-vertex sidecar: multiplying terms with overlapping vertex sets gives zero. Use characteristic zero if integer path counts, rather than their reductions, are required.

Adjoin a source s and sink t. Define M indexed by {s,1,...,n,t} by

    M_(s,v)=z_v,
    M_(u,v)=a_uv z_v,
    M_(u,t)=1,

for original vertices u,v, with all other entries zero. In particular M_(s,t)=0. There is no incoming edge at s and no outgoing edge at t. The powers satisfy M^(n+2)=0. Define M*=(I-M)^(-1)=sum_(j=0)^(n+1) M^j.

**Theorem 3 (Hamiltonian response and associative elimination).**

1. The entry (M*)_(s,t) equals the generating polynomial of all nonempty directed simple paths, with each path weighted by the product of its edge weights and by z_S for its visited vertex set S. In particular

       H(G)=[z_1...z_n](M*)_(s,t).

2. For any interior I⊆{1,...,n}, let B={s,t}∪(V\I) and define the **boundary response matrix**

       E_I=M_BB+M_BI (I-M_II)^(-1) M_IB.

   Then

       ((I-M)^(-1))_BB=(I-E_I)^(-1).

   Every inverse is a finite nilpotent sum. Thus E_I preserves the full simple-path polynomial and H(G).

3. Eliminating disjoint interiors successively gives exactly the same matrix as eliminating their union in one step. All elimination orders agree, retaining the original vertex variables.

**Proof.** A term of a matrix power is a directed walk with a factor z_v on every arrival at an original vertex. For an s-to-t walk, every original visit is an arrival, so a repeated original vertex gives z_v²=0. The surviving s-to-t terms are precisely the simple paths, proving (1). For an arbitrary starting index, a walk of length n+2 has at least n+1 arrivals at original vertices, because s cannot be entered and t can be entered at most once. Thus some arrival repeats, proving nilpotence.

For (2), uniquely divide any walk between boundary vertices into its successive excursions whose intermediate vertices lie in I. A single such excursion has generating matrix E_I; concatenating excursions gives its star. Any repeated arrival at an original vertex still vanishes, including repetitions between different excursions. Thus the two series are equal term by term, including the empty boundary walk. This also proves that E_I is nilpotent: any nonzero product of its entries lifts to a walk in M, whose length is at least the number of factors.

For (3), a walk whose intermediate vertices lie in the union admits the same unique decomposition whether the first interior is grouped before the second or both are grouped together. Equivalently block Gaussian elimination of I-M has invertible pivots here and gives the same Schur complement. ∎

An individual boundary entry can include a return to its initial boundary vertex: that starting vertex is not charged until an arrival occurs. This is intentional. In a complete s-to-t path, the first arrival at that boundary vertex already supplies its variable, so a later return vanishes. Do not interpret every isolated boundary entry as a simple-path polynomial with both endpoints already marked.

This is a genuine local-composition theorem, not just reconstruction from a complete graph encoding. An interior can be replaced by its matrix of boundary excursions, after which it interacts with the exterior only through that matrix. The variables and boundary labels are load-bearing sidecars. The theorem does **not** promise polynomial complexity: R has dimension 2^n, and a faithful response may be large.

### Interface contracts and a dense-graph caveat

If two pieces have disjoint interiors and no edges between their interiors, their excursion corrections add at the common boundary, after counting common boundary edges only once. Arbitrary tournaments do not have such a separator decomposition: two nonempty disjoint interiors necessarily have cross arcs. For arbitrary tournaments the always-valid operation is staged elimination, not unsupported parallel gluing. Module substitutions require their own operation law.

Replacing z_v by 1 is not a ring homomorphism out of R, since it would send z_v²=0 to 1. Ordinary adjacency resolvents count walks and do not preserve Hamilton paths. Replacing the individual variables by a single square-zero variable also fails: it annihilates every path of length at least two. This is an explicit reason the visited-set memory cannot simply be dropped.

## 6. Putting the three theorems in one interface

Since R is finite dimensional over k, every matrix or finite linear operation above can be expanded in the subset basis to obtain a finite-dimensional k-linear state space. Alternatively use the direct path automaton with states e_empty and e_(S,v), v∈S. A letter j appends j if unused and joined by a directed edge, and sends invalid extensions to zero. Initial appends send e_empty to e_({j},j). Let C select and sum states with S=V. Then

    H(G)=C (sum_j A_j)^n e_empty.

This gives an explicit response system to which Theorem 1 applies. Current Hamilton-count output is extremely sparse on partial states, while its contextual closure retains exactly what future extensions can observe. A particular elimination response need not be minimal; quotienting by its response radical is the canonical subsequent linear minimization.

Theorem 2 applies when the relevant operations are invertible transports. Vertex-appending maps in the path automaton are not invertible, so they do not themselves define a covering or a group representation. One must not call their cycles holonomy. Relabelling transports between labelled graph fibers, or separately specified signed-sheet transports as in S170, can supply the invertible part. The common theorem is safe compositional forgetting, not a false assertion that every operation is a symmetry.

## 7. New questions with precise tests

The following are **OPEN research directions**, not proved consequences.

- **Observer completion number.** For an existing pair (A,C), define κ(A,C) as the minimum number of additional scalar rows D needed to make the response radical of (C,D) zero. It measures experimental completion, while dim(V/N) measures memory. They need not coincide: with A=I and C=0, κ=d; with the nilpotent shift and C=0, one first-coordinate observer suffices. Determine κ for S170 signed transport under explicit admissible experiments, distinguishing snapshot probes from repeatable contextual probes.
- **Response holonomy filtration.** Compare the decreasing finite-window blind spaces N_j with a lower central filtration of a specified transport group. Determine the least continuation depth exposing each promised commutator layer. S170's Reed–Muller bank measures many initial sheets at a fixed word; it does not automatically answer this different depth-versus-probe question.
- **Tournament response width.** For a declared boundary and declared exterior operations, determine the minimal linear response dimension, rather than the size of the full square-zero matrix. The cheapest test is to compute ranks of context-versus-state evaluation matrices for small tournaments and compare equal-profile seeds.
- **Compatible higher cells.** Given a physical flip graph and response representation, select only cells whose boundary is response-flat. Compare this response-dependent complex with the actual all-sheet-flat complex. Theorems 1–2 provide a finite exact relation test; they do not supply a canonical choice of cells.

## 8. Reproduction and scope

Run `python3 check.py > observed.json`, then `cmp observed.json results.json`; repeat with `python3 -O check.py` to check that no assertion statements carry validation.

The pure-Python checker exhausts every labelled tournament of orders 1 through 5 and every interior subset, for 33,866 boundary cases. It independently compares the full polynomial against permutation enumeration and verifies one two-stage split for every interior: its sorted first half is eliminated before its remaining vertices. It does not enumerate every elimination order. It also checks the smallest unstable-loss example and sharp horizon controls through dimension eight. The proofs are all-order; these finite checks are hostile and transcription controls, not extrapolation.

Nothing here proves LRC(14), establishes a new tournament Hamilton-count extremum, identifies a scalar profile with a local system, or recovers higher cells from a graph. The deliverable is an exact criterion for when forgetting is safe, an all-order Hamilton-path interface that satisfies it, and a relation test explaining precisely which topology can survive a chosen observer.
