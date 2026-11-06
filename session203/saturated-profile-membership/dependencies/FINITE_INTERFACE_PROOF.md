# A finite weighted presentation of bounded-start profiles

**PROVED, independently audited.** The endpoint theorem and weak-tail
profile preservation are independently audited published dependencies,
frozen under dependencies/ with exact provenance. The independent review,
original candidate, and reproducible controls are preserved under audit/.

The conclusion is an effective finite weighted presentation of all
actual pairs (B,s) with1<=s<=K, for each fixed positive integer K. This is not a finite set of
integer profile states, an automatic recognition of the binary digits of B,
a scalar inverse rule, or a proof that prescribed B-values occur.

## 1. Interfaces and entry arity

For strong T, the endpoint theorem gives

    sum_{x in S} start_x(T) >= floor(|S|^2/4)+1,

for each nonempty S. Therefore at most M=4K-1 vertices have start<=K.
For the actual entry set A in the SCC decomposition of T-root,

    floor(|A|^2/4)+1 <= K/q,
    |A|^2 < 4K/q.

At K=26 every junction uses at most10 entry vertices, from an eligible
reservoir of at most103. No equality between these two different bounds is
being asserted.

A register construction has finitely many individually labeled active
vertices. A vertex introduction selects a free label u and a subset A of
currently active labels: the new vertex beats exactly those old vertices
and loses to every other old vertex. A forget operation removes an active
label, leaving the underlying vertex and all its edges intact. A forgotten
vertex can never be selected again. Thus it beats all subsequently
introduced vertices. Labels may be reused for newly introduced vertices.
The initial empty control has no root marker. A marker may be placed only
on an active vertex. The currently marked vertex cannot be forgotten until
the marker has moved to another active vertex. Final acceptance requires a
marked active root.

Put m(0)=0 and m(t)=floor(sqrt(4t-1)) for positive integer t. The endpoint
theorem says that any marked set whose TOTAL start count is at most t has
cardinality at most m(t). Define

    W(K)=1+max_{0<=a<=K} (m(a)+m(K-a)).

In particular W(26)=15, with a=13 attaining m(13)+m(13)+1=15.
The bound follows also from W(K)<=1+floor(2*sqrt(2K)).

**Proposition 1.** If T is strong and S is a prescribed nonempty marked set
with sum_{x in S} start_x(T)<=K, then T has a register construction using
at most W(K) active vertices simultaneously and retaining exactly S at its
completion.

Proof is by induction on order. A singleton uses one register. Otherwise
choose a root r in S and decompose T-r into C,D1,...,Dt as in the bounded
entry theorem. Let A be r's nonempty out-neighbor set in C and q the product
of H-values of the tail components. Define the required marks in C by

    S_C=(S intersect C) union A.

Endpoint counts of old vertices cannot decrease under induced strong
extensions, by rooted path insertion. Therefore

    sum_{x in S_C} start_x(C)
      <= sum_{x in S intersect C} start_x(T)+start_r(T)/q
      <= sum_{x in S} start_x(T) <= K.

The first inequality is an upper bound even when the two sets overlap;
there is no assumption that their contributions are disjoint. The inductive
hypothesis constructs C retaining exactly S_C. Introduce r, beating
precisely A. This needs at most m(K)+1<=W(K) registers. Immediately forget
every old C-vertex not belonging to S. Their arcs to every future tail
vertex point forward, so this is legal.

Process each tail block in order, retaining all its vertices temporarily
until that block is complete. New backward arcs lead only to r or to
previous vertices of the same block, all still active. Then forget the
block's vertices not in S.

Write a=start_r(T), so1<=a<=K. The retained marks other than r have total
start count in T at most K-a, hence there are at most m(K-a) of them. Every
current tail block D has H(D)<=a, including a singleton of H=1. Applying
the endpoint theorem to all of D gives |D|<=m(H(D))<=m(a). Thus at any tail
stage the number of active vertices is at most

    1+m(K-a)+m(a) <= W(K).

At completion precisely S remains. This is a target-dependent construction
using the actual graph and marked set; no online way to guess S is assumed.
Crucially, replacing the deleted marked root by its entry set preserves the
TOTAL marked start budget. It is unnecessary to retain every vertex whose
individual start count is small.

**Corollary.** Every actual positive rooted profile (B,s),1<=s<=K, has a
representative with a W(K)-register construction. The root lies in the
first SCC of its tournament. Remove singleton later SCCs, which preserves
both B and s by the published weak-tail product bijection. The root SCC has
start s0 with s=q*s0, where q is the product of later H-values. Construct it
using Proposition1 with S consisting only of the root. Append each
remaining nontrivial SCC as an ordered block, keeping that block's vertices
temporarily active. It has H<=q<=K and hence order<=m(K), so the root and
block need at most1+m(K)<=W(K) registers. All final B,s values are unchanged.

## 2. A finite nonnegative path-forest transfer

Fix W and a positive integer p. At a construction prefix, regard the active
labels as ports and all forgotten vertices as bulk. Count spanning directed
path covers by their unordered endpoint pattern: each nonempty component
is represented by (start label,end label), with a special symbol * for a
bulk endpoint. Active vertices that are not endpoints are internal. A
pattern also knows the current active-label set from the finite control.

Keep precisely those patterns with at most p components starting in bulk.
Every other component starts at a distinct active port, so there are at
most W+p components. There are only finitely many patterns for fixed W,p.
Repeated components with endpoint pattern (*,*) are allowed as a multiset.
Counts refer to actual unordered collections of vertex-disjoint directed
paths; two paths with the same endpoint type are still different paths.

For introduction of a new active vertex u with old out-neighbors A, form
all possibilities in which u is:

- a singleton component;
- appended to one old component whose end is * or a port outside A;
- prepended to one old component whose start is a port in A;
- placed between two distinct old components satisfying those end/start
  conditions, merging them into one component.

The two selected components must be distinct; using the same one would
form a cycle, not a path. Multiplicities from indistinguishable endpoint
patterns count the distinct choices of actual old components. They depend
only on the finite pattern, not on hidden vertex identities. Every child
cover has a unique inverse after deleting u and its incident path edges.
When u was internal, always split its former path into the two old
components, even if its old predecessor happens to beat its old successor.
No old predecessor-successor edge is added in this inverse. Thus this is an exact nonnegative integer linear
transfer on pattern counts.

On forgetting port u, replace its endpoint occurrences by *. If the number
of bulk-start components now exceeds p, discard the pattern. This does not
lose a cover that could contribute to a final cover with at most p paths:
future vertices cannot precede a bulk-start vertex, so the number of
bulk-start components can never decrease in a future legal construction.
All other patterns are merely regrouped. Forgetting is another exact
nonnegative integer linear transfer.

More strongly, the transfer maintains exactly the counts of all prefix
covers with at most p bulk-start components. Introduction cannot decrease
the number of bulk starts: a component merged on the right must start at
an active port. Hence discarded old covers cannot reenter the retained
class. The empty graph starts with one empty-cover state.

## 3. Physical B,s as nonnegative final readouts

Take p=2. At the end retain the root port r and forget any other ports.
Let w(pattern) be the exact pattern count. Then

    h = sum w over one-component patterns;
    s = sum w over one-component patterns starting at r;
    L = sum w over one- or two-component patterns
                          having a component ending at r.

In the two-component case the root-ending component is uniquely the
nonempty designated first path, and the other is the ordinary second path.
In the one-component case that second path is empty. Therefore these are
exactly the usual h,s,L, and B=h+L-s. Equivalently B has a nonnegative
readout: add every one-component pattern not starting at r, then add every
one- or two-component pattern ending at r. A one-component path counted
in both classes properly has coefficient2.

The same device supplies all fixed ordinary cover evaluations. With p>=k,
forget all ports and weight a q-component unordered cover by the falling
factorial k*(k-1)*...*(k-q+1), assigning its nonempty paths to distinct
labeled bins. Rooted designated-cover readouts use the analogous
falling-factorial weight for the ordinary components. Thus a fixed finite
cover requirement does not require an unbounded hierarchy.

## 4. What is finite and what remains unproved

The introduction/forget/root-mark alphabet is finite for fixed W. The
active-label set and root marker are finite control, and all forest
transitions are fixed nonnegative integer matrices. Their dimensions can
be enormous; no useful numerical implementation at W=15 is claimed.

Because the matrices and the s-readout are nonnegative, capping every
coordinate at K+1 commutes with all transitions and with the capped
readout. The finitely many capped vectors, together with control, recognize
exactly the legal words whose final root start satisfies1<=s<=K. For each
fixed s in that interval the accepted word language is therefore regular.
Keep the uncapped vector in parallel to read the exact B.

Every accepted word constructs an actual rooted tournament. Conversely,
Proposition1 and its weak-tail corollary represent every actual pair(B,s)
in the stated budget. Hence this is an effective finite weighted
presentation of exactly the bounded-start profile set, not merely a
sufficient subclass.

This does not show that the set of numerical B-readouts is semilinear or
regular in a positional base. Even a one-dimensional weighted system can
have the sparse image {2^n}. No arbitrary inverse branch is supplied here.
The remaining F(8,137) question is still whether every prescribed odd X
has one accepted word with (B,s)=(X-2j,r+3j) in each residue r=1,2.
