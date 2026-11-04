# Uniform module-port compression preserves an exposed-anchor separator

**PROVED, alternate strengthening of the exposed-anchor method; independently reviewed.** The incoming E1 whole-converse parity identity survives identifying all port bits inside any disjoint module blocks, provided the anchor is an uncovered singleton. Consequently a separator can be chosen individually uniform on those modules. This independently settles arbitrary nested pairs and exposed-anchor laminar families, with a smaller compressed query battery. It does not cover a partition whose every block is nontrivial.

The companion independent [averaged proof](../simultaneous_response/ANCHORED.md) obtains the same family-separation scope by sums of singleton-port queries. This proof instead shows that uniform-port queries themselves suffice; it is not needed to justify that earlier averaging result.

## U1. The initial-run switching involution

Work in the incoming squarefree path algebra over F2. Let M be a module of a tournament B and let v lie outside M. Write s_u for all simple paths starting at u, and G_wv for all simple paths from w to v, retaining complete vertex support. Then, coefficientwise,

    sum_(u<w in M) (s_u G_wv+s_w G_uv)=0.             (1)

Products already discard intersecting path supports. Thus the left side counts ordered disjoint pairs (P,Q), both starting in M, with Q ending at v. Their starts are distinct automatically.

Decompose each path into its initial nonempty consecutive M-run and its remainder. Q has a nonempty remainder because v is outside M. P's remainder may be empty. Swap the two initial M-runs, keeping P's and Q's respective remainders in their old slots. If a remainder is nonempty, its first vertex is outside M. Since M is a module, every vertex of M has the same incidence to that first remainder vertex, so the new join is directed whenever the old join was. Internal run arcs and remainder arcs are unchanged.

The resulting paths remain disjoint and have the same union of vertex supports. Their initial M-runs are precisely the exchanged runs, even if a remainder later reenters M. Applying the operation twice restores the pair. The runs are disjoint and nonempty, so the pair cannot be fixed. This is a support-preserving fixed-point-free involution, proving (1).

For the whole-converse contrast delta_B=(F_B-F_(B^op))/2 modulo2, incoming T2 gives

    [a_u a_w c_v] delta_B
       = L_odd((U+1)(s_u G_wv+s_w G_uv)),  u!=w.

Consequently

    sum_(u<w in M) [a_u a_w c_v] delta_B=0.           (2)

This includes every such cubic coefficient, not only shared-index cubics. The involution applies before the support functional and therefore preserves its weights exactly.

## U2. Aliased coefficients and the surviving identity

Assume |B|>=3. Partition the names of B into pairwise disjoint module blocks and the remaining singletons. Choose an uncovered singleton anchor v. A block may also have size1. Let q be the total number of blocks. Identify a_u=a_C and c_u=c_C whenever u lies in block C, obtaining the Boolean polynomial bar_delta.

The original whole-converse contrast has degree at most3. In reducing a_u a_w to a_C when u,w lie in C, cubic terms can in principle contaminate the mixed quadratic coefficient [a_C c_v]. Equation (2) proves that their total contribution is zero. Therefore, for every block C other than {v},

    Dbar_Cv := [a_C c_v] bar_delta = sum_(u in C) D_uv,
    hbar_Cv := [a_C a_v c_v] bar_delta = sum_(u in C) h_uv. (3)

There is no other aliasing route to these monomials: c_v belongs to the singleton anchor, and the whole-converse polynomial has no quartic part. Higher terms of an arbitrary ambient internal-reversal contrast are not being discarded.

Each block C has one uniform incidence to v, since it is a module. Grouping the incoming E1 identity at v gives the exact compressed relation

    sum_(C!={v}) Dbar_Cv + sum_(v->C) hbar_Cv = 1.      (4)

The compressed block count q may be2: (4) is inherited from the original tournament on at least3 names, not asserted as E1 for the quotient tournament itself. Thus bar_delta cannot vanish identically. Some actual query constant on each selected block separates B from B^op by2 modulo4.

The hypotheses are necessary to this proof: v must remain a singleton, and blocks must be modules. When every block has size at least2, the coarse C3[TT2,TT2,TT2] example in PROOF.md has all uniform-query contrasts0 modulo4, so no such unconditional extension is asserted.

## U3. A compressed explicit query battery

Let b_C=sum_(u in C)(S_u+E_u) modulo2. Incoming singleton tomography gives

    bar_delta(a,0)=sum_C b_C a_C,   sum_C b_C=0.

At anchor v, let d be the number of quotient blocks C with v->C. This is a number of blocks, not the original outdegree in vertices. Set c=e_v. Include the q-1 main compressed queries

    a=e_C      when C->v,
    a=e_C+e_v  when v->C.

Also include a=e_v if d is odd and a=0 if q-1+d is odd. All listed masks are distinct. Their sum satisfies

    sum bar_delta(a,e_v)=1+b_v.                       (5)

Indeed put f(a)=bar_delta(a,e_v), substitute

    Dbar_Cv=f(e_C)+f(0)+b_C,
    hbar_Cv=f(e_C+e_v)+f(e_C)+f(e_v)+f(0)

into (4), and collect terms. The coefficients of f(e_v),f(0) are d and q-1+d. The remaining b sum is b_v.

If b_v=0, one query in this battery separates. If b_v=1, the single additional query (e_v,0) separates. The battery has q queries when q is even; q-1 or q+1 when q is odd. Including the endpoint query gives at most q+2 candidates. This is detection for this comparison, not reconstruction of the complete response.

## U4. General nested pairs and exposed laminar families

For A properly contained in B, both nontrivial modules of an ambient tournament, choose v in B\A and use A as one compressed block, with every other name of B left a singleton. U2-U3 provide a whole-B-converse separator uniform on A. In that actual port extension, A remains a module, so conversing A has no effect on the integer Hamilton count: all its ordered path-cover counts are converse-invariant. Thus the same individual query separates the result of simultaneously reversing A and B from the source. Give the ports B's exterior incidences to lift through the ambient quotient; T1 multiplies the core residue by its odd H, preserving2 modulo4.

More generally, in a laminar selected family choose a maximal selected module B and an exposed anchor v belonging to no smaller selected module inside B. Compress the maximal proper selected submodules inside B, which are pairwise disjoint, and leave other names singleton. The resulting query is uniform on every proper selected submodule, including all descendants. Each of their reversals is individually invisible to that query over the integers. Reversals in other maximal selected modules are likewise invisible because the lifted ports are uniform there. Only the whole-B contrast remains, and is2 modulo4.

This is stronger at the query level than averaging: the separator can be chosen uniform on all proper selected modules. The same family boundary remains. Fully covered laminar families and arbitrary larger crossing families are not settled here.
