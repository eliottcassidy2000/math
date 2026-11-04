# A marked pair exposes an ordered cut after quotient deletion

**PROVED:** quotient converse with unchanged child interiors is detected whenever one child is a homogeneous pair and deleting its quotient vertex leaves a non-strong quotient, provided at least three old exterior names remain. The other child tournaments are arbitrary. More generally, a grouped quartic recovers every exterior ordered-cover skew between unions of child blocks. This gives a new criterion inside the strong-prime quotient frontier, not a solution for all prime quotients.

## 1. Physical setting and exact coefficient

Let Q be a tournament on fixed child indices, and let T=Q[A_i] be its substitution by nonempty, disjoint, physically named child tournaments. Let U=Q^op[A_i], retaining every child's internal orientation and every physical label. Fix one child A_z={p,q}, with p->q. Write R=T-{p,q}.

Choose disjoint nonempty unions A,B of the other child blocks. They need not be modules of T as unions. Attach x->y and use variables a_v=[v->x], c_v=[v->y]. Identify all a variables on A to a_A and all c variables on B to c_B. All other exterior first-port variables are zero and all other exterior second-port variables are zero. Keep the pair variables separately, setting a_p=c_q=0 and varying a_q,c_p. Thus every exterior child is individually seen uniformly by both ports.

Let C_AB(R) count ordered spanning covers (P1,P2) of R by two nonempty directed paths, with P1 ending anywhere in A and P2 starting anywhere in B. The connecting arc is not tested. The grouped Boolean coefficient satisfies

    [a_q a_A c_p c_B](F_T-F_U)=C_BA(R)-C_AB(R).       (1)

This coefficient is an alternating sum of sixteen actual physical attachment queries. It retains degree four and uses four distinct Boolean variables. Since A and B are disjoint and exterior, every underlying old-name quartic has four distinct physical names.

### Proof of (1)

Let W be T with only the pair p->q reversed. Whole converse W^op has quotient Q^op and pair p->q; all other child interiors are conversed. Since the selected ports are uniform on each of those other children, restoring each child's orientation leaves its actual Hamilton count unchanged, by reversing every maximal child run. Therefore

    F_U = F_(W^op)

on every query in this grouped four-variable family, exactly over the integers.

The complete quartic part of any tournament response is unchanged by whole converse. Substituting grouped variables and constants cannot increase degree, so the degree-four part of F_(W^op) equals that of F_W on this family. Hence the left side of (1) is the grouped quartic coefficient of F_T-F_W.

For distinct exterior names i,j, the incoming exact twin derivative gives

    [a_q a_i c_p c_j](F_T-F_W)=C_ji(R)-C_ij(R).      (2)

For completeness, every surviving word uses the changed pair edge. The marked port neighborhoods force the block j,y,p,q,x,i in T and i,x,q,p,y,j in W. Removing this block's pair and ports leaves the indicated two-path cover of R. Both quartic signs are positive. No pair-to-exterior incidence is tested in these blocks, so the formula is independent of Q's incidences at z. Summing (2) over i in A,j in B proves (1). No quartic truncation of the actual quotient-converse comparison is used.

## 2. Exact positive skew across an ordered cut

Suppose R=A->B, with A and B nonempty. They may have arbitrary internal tournaments. Then

    C_AB(R)=H(A)H(B),
    C_BA(R)=(H(A)+q_2(A))(H(B)+q_2(B)),             (3)

where q_2 denotes ordered covers by two nonempty paths. Consequently

    C_BA-C_AB=q_2(A)H(B)+H(A)q_2(B)+q_2(A)q_2(B).  (4)

To prove the first equality, a path ending in A contains no B vertex, and a path starting in B contains no A vertex. They must therefore be Hamilton paths of A and B. For the second, the first path must contain B and the second must contain A. Each path either stays on its required side or crosses once from A to B. The four possibilities give respectively H(A)H(B), q_2(A)H(B), H(A)q_2(B), and q_2(A)q_2(B). These cover slots remain ordered, so there is no extra factor two.

Every nonempty tournament has a Hamilton path. Moreover q_2(D)>0 exactly when |D|>=2: split a Hamilton path for sufficiency. Thus (4) is strictly positive whenever |R|>=3. For |A|=|B|=1 it is zero; that boundary is not claimed to be detected by this quartic.

Combining (1) and (4), at least one of the sixteen actual queries separates T and U whenever the exterior admits such an ordered cut into unions of child blocks. The statement gives an exact integer separation, not necessarily a nonzero residue modulo four.

## 3. Consequence for a strong prime quotient

If Q-z is not strongly connected and has at least two vertices, its strongly connected components are linearly ordered. Split that order into a nonempty prefix and suffix, and let A,B be the corresponding unions of full child blocks. Then R=A->B. If R has at least three physical names, the criterion applies.

In particular, this holds for any strong prime quotient of order at least five with a vertex z whose deletion is not strong, whenever the child at z is a pair. All other children may be arbitrary nonempty tournaments; if they are also nontrivial, this is a fully covered laminar operation and is not an exposed-anchor case.

A concrete prime quotient is

    Q=(14,28,8,16,5)

in outgoing-row bitmask notation on names 0,...,4. Its first four names form the transitive order 0->1->2->3; name 4 beats 0,2 and loses to 1,3. It is strong and prime, and Q-4 is transitive. Inflate name 4 by a pair and the other names by arbitrary tournaments. The theorem separates Q[A_i] from Q^op[A_i] with all interiors unchanged. For all five children equal to a transitive pair, take A as the first two blocks and B as the next two. Then A and B are transitive of order four, and the quartic coefficient is

    14+14+14*14=224.

The checker confirms this through its entire sixteen-query signed battery on ten physical old names. The quotient's primeness and strong connectivity are checked directly; neither property is assumed from its conventional name.

## 4. Scope and dependencies

The general coefficient criterion (1) remains useful when Q-z is strong, but (4) no longer supplies automatic positivity. The cyclic five-vertex regular quotient is a concrete residual example: deleting any quotient vertex remains strong, so this corollary does not apply. No collision is claimed there.

The special child must be exactly a pair in this proof. An arbitrary larger child is not silently replaced by a pair. The result concerns the actual substitution T and U; arbitrary strong outer ambient lifting is not asserted. Ordered SCC padding may be handled by the existing exact localization theorem.

The twin derivative and whole-converse quartic invariance are inherited from the independently verified `../response/PROOF.md`, Sections 2 and 5. The new step is applying them to quotient converse with unchanged interiors, using individual child uniformity, then evaluating the grouped skew across an ordered cut. All physical labels and exact degree-four information are retained.

## 5. Reproduction

The standard-library verifier recomputes ordered-cover matrices by subset endpoint convolution and recomputes every signed battery term as a full augmented Hamilton count. It checks general grouped skews on a declared seeded sample, the ordered-cut formula on arbitrary small tournaments, and strong-prime W5-type quotient controls with both transitive and cyclic children. Normal and optimized receipts are identical. Finite controls supplement the all-order arguments above.
