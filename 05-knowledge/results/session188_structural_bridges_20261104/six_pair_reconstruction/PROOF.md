# Named reconstruction on a known partition into six pairs

**PROVED:** Fix twelve physical names partitioned into six named pairs. Among tournaments in which every pair is a module, the complete named two-port response distinguishes every quotient orientation and every internal pair orientation. There are 2^(15+6)=2,097,152 such tournaments. The same statement holds when this twelve-name core is substituted into any fixed rooted ambient tournament, with its exterior unchanged.

This is reconstruction on a specified, known pair partition. It is not reconstruction for arbitrary twelve-name tournaments, an advance in unrestricted named reconstruction, or an exhaustive response-bank claim for the 2,097,152 states.

## 1. Physical observer and inherited identities

Let B_i={2i,2i+1}, i=0,...,5. A tournament in the class is Q[B_0,...,B_5], where Q is an arbitrary named six-vertex tournament and each B_i has either internal direction. The external direction between two pairs is uniform. Add physical ports x->y. A bit a_v=1 means v->x; c_v=1 means v->y.

A coarse query has a and c constant on each B_i. Write f_Q(a,c) for its actual Hamilton-path count. This is independent of all six internal pair directions: reversing every maximal run inside one module gives an exact path bijection when the ports see that module uniformly.

The companion `../even_pair_refinement/PROOF.md`, equation (6), proves that the coarse quotient-converse contrast satisfies

    [a_i a_j c_j](f_Q-f_(Q-op)) = 4 h_ij (mod 8),
    h_ij = sum_(S subset V(Q)-{j}, |S| even)
                      (S_i(Q[S])+E_i(Q[S])) (mod 2),      (1)

for i!=j. Here S_i and E_i count Hamilton paths starting and ending at i, zero when i is absent. This is a proved coefficient of the full physical response, including the correct cubic and quartic contributions; no degree truncation is used.

Whole path reversal followed by exchanging x,y gives

    f_(Q-op)(a,c)=f_Q(1-c,1-a).                            (2)

Indeed whole converse reverses every pair as well as the quotient, but pair directions are invisible to these coarse queries. Thus the data h_ij in (1) can be computed from the response of the single unknown tournament, by exact finite differences of the difference in (2).

## 2. A five-name endpoint identity

**Lemma.** For every tournament R on five names,

    sum_(S subset V(R), |S| even)(S_i(R[S])+E_i(R[S]))
       = outdegree_R(i) (mod 2).                          (3)

The empty set contributes zero. Every two-name support containing i contributes one, and there are four, so only the four-name supports matter.

To evaluate one four-name support {i,a,b,c}, pair each candidate path starting at i with its reversed candidate ending at i. Write A_uv=1[u->v]. If the starting path is i,a,b,c, its two validity indicators are

    A_ia A_ab A_bc + (1-A_ia)(1-A_ab)(1-A_bc).

Modulo two this equals one plus the three linear and three quadratic terms in those factors. Sum over the six permutations of a,b,c. The constant and linear terms cancel. The quadratic terms give

    E_i+S_i
      = 1 + sum_(u in {a,b,c}) A_iu
          + sum_(u!=v in {a,b,c}) A_iu A_uv (mod 2).      (4)

For clarity, the terms A_iu A_vw with u,v,w distinct sum to sum_u A_iu, since A_vw+A_wv=1. The terms A_uv A_vw sum to the Hamilton-path count of the other three names, which is odd, giving the constant one in (4). The remaining terms are the displayed ordered quadratic sum.

Now sum (4) over the four choices of a three-subset of V(R)-{i}. The constant occurs four times. Each A_iu occurs three times. Each A_iu A_uv occurs twice. Exactly the degree parity remains, proving (3).

Applying (3) to R=Q-j gives the explicit six-name formula

    h_ij = d_i + A_ij (mod 2),
    d_i=outdegree_Q(i) (mod 2),       i!=j.                (5)

## 3. Recovering the quotient

First, h cannot vanish identically. If it did, (5) would make every row of the tournament adjacency matrix constant off the diagonal. Every vertex would then be a source or a sink. A tournament with at least three names cannot have this property: at most one vertex can be a source and at most one a sink.

Second, h determines Q up to whole converse. Suppose Q and Q' have the same h. Equation (5) gives

    A'_ij+A_ij=d'_i+d_i.

The left sides for (i,j) and (j,i) are equal, since both adjacency matrices are tournaments. Hence all six bits d'_i+d_i are one common bit. If it is zero the quotients coincide; if it is one every quotient edge is reversed.

An explicit recovery is available. Choose a trial value d_0 in {0,1}, then set

    d_j = 1+h_0j+h_j0+d_0,
    A_ij = h_ij+d_i,                     all modulo 2.

The two trials produce the two opposite quotient candidates. Since h has a nonzero entry, (1) guarantees that these two candidates have different actual coarse responses. Their ambiguity is therefore resolved by the same query data used to compute h.

## 4. A fixed query battery and pair decoding

Let L be the following 94 coarse mask pairs:

- c=0 and a=0, any singleton, or any two-element mask: 22 queries;
- for each j, c=e_j and a=0, any singleton, or any two-element mask containing j: 6*12=72 queries.

Every coefficient in (1) is a finite difference on eight queries in L. Add the reflected queries (a,c)->(1-c,1-a). They are disjoint from L: their a-weight is five or six, whereas queries in L have a-weight at most two. The resulting fixed battery has188 actual coarse queries and determines every h_ij via (1)-(2). If the two opposite quotient candidates agreed on this battery, their nonzero shared-cubic contrast could not be recovered from it. Thus the battery recovers the quotient orientation, not merely its converse class.

After Q is known, six further physical queries read the pair directions. For B_i={p,q}, use internal a_p=c_p=1 and a_q=c_q=0, with x->y. Give both ports B_i's exterior incidences: on every other pair their two attachment masks are all zero or all one according to its quotient edge to i. The augmented four-name pair core has Hamilton count1 for p->q and5 for q->p. It is a common module, and every other pair remains uniformly seen by the ports. Its orientation bit is therefore independent of all five other pair directions and remains visible modulo eight by the transfer below.

Consequently194 adaptive actual queries suffice:188 fixed coarse queries followed by six pair decoders. A fully fixed battery can instead include every possible exterior incidence pattern for each pair: six pairs times2^5 patterns, hence192 decoder queries. Together with the coarse battery this gives380 fixed queries. These are explicit nonoptimal upper bounds, not separation lower bounds.

## 5. Arbitrary fixed ambient context

Let K be any fixed rooted tournament and replace its distinguished root by the twelve-name core. The exterior labels and directions are identical for all states. For each coarse query above, give x,y the root's exterior incidences, so the actual fourteen-name augmented core is a common module.

The separate even-core theorem in `../mod8_transfer/PROOF.md` states that for same-size even-order A,B,

    H(K[A])-H(K[B])
        = H(K[TT2]) (H(A)-H(B)) (mod 8),                 (6)

and H(K[TT2]) is odd. The multiplier depends only on the fixed context, not on the attachment query or Q. Only the six core attachment bits are reflected; all exterior port incidences stay fixed. Subtracting a query from its reflection, and then taking the coefficient finite differences, multiplies4h_ij by this same odd integer and therefore leaves it unchanged modulo eight. The quotient-recovery argument is unchanged. A nonzero coarse contrast also remains nonzero, resolving its whole-converse ambiguity in the ambient context.

For each pair decoder, the common augmented pair core has even order four. Its exterior rooted quotient is fixed once Q is known, apart from other pair interiors, which are exactly invisible to this query by module run reversal. Equation (6) sends its core gap4 to the same gap4 modulo eight. Thus each pair bit is independently decoded in the full ambient response.

This proves that all2^21 members of the known-partition class have distinct named responses in every fixed ambient context. No assumption that Q or the ambient quotient is strongly connected is needed. The exterior must be the same for the compared states; arbitrary changing ambient tournaments are outside the statement.

## 6. Reproduction and evidence boundary

The standard-library verifier checks (4) by literal four-name Hamilton paths, checks (3) by integer induced-subset endpoint DP on every five-name tournament, and checks (5) on all32,768 six-name quotients. It retains the complete30-bit h vector for every quotient and groups the exact vectors, obtaining16,384 classes, each exactly one whole-converse pair. Hashes bind the retained data but are not used as equality substitutes.

A separate C++17 literal physical subset DP retains full exact query values for the188-query battery on selected twelve-name states, and all64 pair-direction states for one six-name quotient under the six decoding queries. It also checks a genuine strong ambient substitution. These controls use the full augmented tournament and do not implement the endpoint formula or a truncated response polynomial.

Run `python verify.py > normal.json` and `python -O verify.py > optimized.json`, then compare the receipts. Checks raise explicit exceptions under both modes; no software is installed. The all-order-in-the-ambient theorem follows from the proofs above, not from extrapolating the finite controls. The entire2^21-state response bank is not enumerated.
