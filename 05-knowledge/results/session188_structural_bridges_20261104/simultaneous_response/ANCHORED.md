# Exact port averaging excludes every nested pair and exposed laminar family

**PROVED:** simultaneous reversals of a proper nested pair of distinct nontrivial modules cannot cancel in the complete named two-port response. More generally, every laminar family with an exposed anchor, defined below, has an actual query contrast equal to 2 modulo 4. This includes the disjoint theorem. **OPEN here:** fully covered laminar families and unrestricted crossing families. No exact collision is claimed.

## 1. Exact singleton-port averaging

Let A be a module in an old tournament R, with m=|A|>=1. Add the two fixed named ports x->y. Fix all attachments outside A, and assume y sees all of A uniformly. For each u in A let x's attachment bits on A be exactly e_u, meaning u->x and x->w for every w in A-u. Let R' converse only A.

**Averaging lemma.** Exactly over the integers,

    sum_(u in A) F_R(a_out+e_u,c)
      =sum_(u in A) F_R'(a_out+e_u,c).                       (1)

No condition is imposed on x's fixed incidences outside A, and y may uniformly dominate A or uniformly lose to A.

**Proof.** First ignore the attachment test at each x-A contact in a Hamilton word. Reverse each maximal contiguous A-run, retaining the order of runs and every outside name. This is an involutive bijection between the resulting partially tested words for R and R'. Every other boundary arc remains valid because all outside vertices, including y, see A uniformly.

Now sum over the m possible singleton attachments at x. A Hamilton word has at most one incoming and at most one outgoing x-A contact. The number of legal singleton choices is respectively

    no contact: m;
    incoming A->x only: 1;
    outgoing x->A only: m-1;
    both contacts: 1.

In the last case the two A endpoints are distinct, since Hamilton words never repeat a name. Reversing A-runs may change endpoint identities but not their contact types or distinctness, so the same weight is assigned to the paired word. Summing proves (1). The proof retains singleton A-runs and all possible positions of x and y. QED.

This is an equality of a **sum of actual queries**. It does not assert that an individual nonuniform query is invariant under a module converse.

## 2. An all-order form of the incoming anchor battery

For a tournament B on n>=3 names, put

    delta_B=(F_B-F_(B^op))/2 modulo 2,
    b_v=(S_v(B)+E_v(B)) modulo 2,
    d=outdegree_B(v).

The division is over integers. The incoming [E1 parity identity](https://github.com/eliottcassidy2000/math/blob/76b3b1a77d124100f03026420f6b821533cdec15/05-knowledge/results/two_port_parity_separation_20261104/PROOF.md) holds at every n>=3. Its finite-difference calculation gives the following explicit extension of incoming E2 to odd n as well.

Let c=e_v. Include the n-1 main queries

    a=e_u                 if u->v,
    a=e_u+e_v             if v->u,

one for every u!=v. Also include a=e_v if d is odd, and include a=0 if n-1+d is odd. These masks are all distinct. Call the resulting battery Q_v. Then

    sum_(a,c in Q_v) delta_B(a,c)=1+b_v modulo 2.           (2)

At even n there are exactly n queries. At odd n there are n-1 queries when d is even and n+1 when d is odd. If b_v=1, the single query (e_v,0) already has delta_B=1.

**Derivation.** Put f(a)=delta_B(a,e_v). Incoming singleton tomography gives delta_B(a,0)=sum_u b_u a_u and sum_u b_u=0. The actual mixed quadratic and shared cubic coefficients are

    D_uv=f(e_u)+f(0)+b_u,
    h_uv=f(e_u+e_v)+f(e_u)+f(e_v)+f(0).

Substitute these into E1: sum_(u!=v)D_uv+sum_(v->u)h_uv=1. The surviving coefficients of f(e_v),f(0) are d and n-1+d, and the remaining b sum is b_v. This is exactly (2). No even-order assumption is used. The n=2 boundary is excluded from E1; it is handled directly by its nonzero endpoint parity b_v=1. QED.

## 3. General nested-pair theorem

Let A properly contained in B be two nontrivial modules of T. Reverse both interiors by simultaneous edge XOR; equivalently, converse B and then converse A, since nesting preserves modularity. Call the result U.

**Theorem.** Some actual named two-port query satisfies

    F_T-F_U=2 modulo 4.

**Proof.** Choose any v in B-A. Apply the Section 2 anchor battery to the induced tournament T[B]. Lift every query to T by giving both ports B's exterior incidences. The incoming modulo-four substitution theorem multiplies every whole-B-converse contrast by the odd integer H(T/B). Thus the lifted battery sum is 2 modulo 4 when b_v=0. If b_v=1, the lifted endpoint query already has contrast 2 modulo 4.

It remains to show that conversing A after B cannot alter this certificate. In the battery, the queries indexed by u in A have c=0 on A and have the same x attachment outside A: a=e_u if A->v, or a=e_u+e_v if v->A. The module property makes the choice uniform over all u in A. Their sum is invariant under A converse by (1). Every other battery query, including the optional final queries, has both ports uniform on A, so each count is individually invariant by the ordered-path-cover substitution formula. The endpoint query (e_v,0) is likewise individually invariant.

Hence restoring A changes neither the battery sum nor the endpoint certificate, exactly as integer counts. Each actual query contrast is even because every extended tournament has odd Hamilton count. A battery sum equal to 2 modulo 4 therefore contains at least one individual query contrast equal to 2 modulo 4. The modular separator survives. QED.

This proof makes no assumption that B-A is a module or that its incidences to A all point in one direction as a whole. It uses only the actual module condition: each outside vertex separately sees A uniformly.

## 4. Exposed-anchor laminar families

A family of modules is **laminar** if two members are disjoint or one contains the other. A selected maximal module B has an **exposed anchor** if some v in B belongs to no proper selected module contained in B. The names here are definitions of the comparison class, not inferred structure.

**Theorem.** Any nonempty laminar family of nontrivial modules with an exposed anchor in at least one maximal selected module is separated from its simultaneous converse by an actual query contrast 2 modulo 4.

**Proof.** Choose such B,v and lift its anchor battery as above. First converse B. Every proper selected module A inside B avoids v, so the exact argument from Section 3 makes its reversal invisible to the battery sum, and individually invisible to the endpoint certificate. Apply this successively. Laminarity ensures that every selected set remains a module: reversing a containing set only reverses uniform boundary directions; reversing a contained or disjoint set does not alter its uniform boundary.

Every selected module outside B lies inside a disjoint maximal selected module. The lifted ports see all those names uniformly, so each such reversal is individually invisible to every actual query. The remaining contrast is therefore the whole-B certificate, which is nonzero modulo four. QED.

Every pair of distinct nested nontrivial modules has an exposed anchor in B-A. Every pairwise-disjoint family also qualifies. The theorem stops when every maximal selected module is covered by proper selected modules. The six-vertex whole-plus-three-pairs example in `PROOF.md` is such a fully covered laminar family; it is separated explicitly, but not by an asserted application of this theorem.

## 5. Exact checks and retained full responses

`anchored.py` consumes the complete integer response banks generated by `search.py`, including all cubic and quartic coordinates. Through five old names it checks every exposed laminar family and every available exposed anchor. It verifies both the modular certificate and the **exact equality** between the whole-B-converse battery sum and the all-reversals battery sum. No response is truncated at degree three.

It separately checks (1) through four old names for every original module, all exterior attachment masks, and both uniform orientations of y on A. The n=2 failure boundary of E1 is retained and handled only by endpoint parity. Fully covered laminar families are counted without being claimed covered by this proof.

    python anchored.py > anchored.json
    python -O anchored.py > anchored.optimized.json

The full response banks and parent census are in `search.py` and its receipts. Finite tests support the implementation and proof controls; the unbounded result is the weighted run-reversal bijection and the incoming parity theorem.
