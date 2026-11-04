# Equal-profile cyclic children through order seven are visible in every ambient quotient

**PROVED:** let three nonempty physically named child tournaments have the same complete ordered-path-cover profile and common order n<=7. Reversing their cyclic quotient, with arbitrary profile-preserving changes inside the children, is separated by one explicitly lifted physical two-port query in every fixed ambient tournament. The contrast has a strict integer sign, not merely a nonzero modular residue.

## 1. Exact hypotheses and query

Let A1,A2,A3 be disjoint, nonempty, fixed named sets, each of size n, where1<=n<=7. Put M=A1 union A2 union A3. Let the source core be

    T=C3[A1,A2,A3], with A1->A2->A3->A1.

Each A_i denotes its specified internal tournament. Assume all three have the same complete ordered-cover profile q, where q_k counts ordered spanning covers by k nonempty directed paths, q0=0, and q_k=0 for k>n.

Let the target core U have the opposite cyclic quotient on these same sets. Its internal tournaments may be unchanged, individually conversed, or otherwise replaced, provided each target child still has the same complete profile q. Physical names are not relabelled.

Let Q be any fixed tournament with a distinguished root r. Replace r by T or U; all other vertices and arcs of Q remain fixed. There is no restriction on the order, strong connectivity, or modular decomposition of Q.

Attach actual ports x->y. On the core, use

    a=1_(A1),     c=1_(A2),

where a_v=[v->x], c_v=[v->y]. For each exterior vertex z of Q-r, set both port incidences to match the root incidence:

    if z->r in Q, use z->x and z->y;
    if r->z in Q, use x->z and y->z.

These assignments define one actual named query on the entire source and target. They are identical for both. The claimed result is

    F_(Q[T])(a,c) - F_(Q[U])(a,c) < 0.                  (1)

The core masks in the notation above are extended by the exterior bits just specified; they must not be read as zero masks on every exterior vertex.

## 2. Ordered-cover domination of the augmented cores

Let T+ and U+ denote the port-augmented cores only, each on3n+2 vertices. Because the ports are uniform on each child, every ordered-cover count of these cores depends on each internal child solely through its complete q profile. Profile-preserving target changes therefore have no further effect.

The exact all-coordinate identity proved in `../covered_ambient/PROOF.md` is

    q_h(T+) - q_h(U+)
      = D_h(q)
      = (h/3) L_q[K E0^(h-1)/(1-XYZ)^(h+1)]            (2)

for every h>=1, with its stated K,E0 and monomial functional L_q. This identity retains paths in which the two ports occupy the same or different ordered slots. It is not a cubic truncation of a named response.

The classical odd-cycle expansion of Grinberg and Stanley, recorded with its primary citation in `../profile_cone_obstruction/PROOF.md`, writes every actual order-n tournament profile as

    q=t^(n)+sum_(s<n, s=n mod2) w_s t^(s),
    t^(s)_k=k!S(s,k),     w_s>=0.

The leading coefficient is one. Moreover t^(s)_1=1, so H(A_i)=q1=1+sum_s w_s. Counting vertex permutations gives H(A_i)<=n!. Thus q lies in exactly the continuous domain certified in `../profile_cone_bounded/PROOF.md`.

For1<=n<=7 that exact polynomial certificate proves

    D_h(q)<=0 for every h>=2.                           (3)

The independently proved scalar positive-square identity in `../covered_response/PROOF.md` and `../covered_positivity/PROOF.md` gives

    D_1(q)=-P(q)<0.                                    (4)

It applies because q is finite, nonnegative and nonzero. No log-concavity assumption is needed. Combining(2)--(4), every ordered-cover coordinate of T+ is at most the corresponding coordinate of U+, with strict inequality for their Hamilton counts.

## 3. Arbitrary ambient lift and quantitative sign

The union M union{x,y} is a module in the full augmented source and target, precisely because the ports were given the root's exterior incidences. Collapsing this module gives the same rooted ambient Q in both cases.

The exact module-run substitution formula supplies coefficients c_h(Q,r)>=0, independent of the internal augmented core, such that

    H(Q[V])=sum_(h>=1) c_h(Q,r) q_h(V).

Its first coefficient is c1(Q,r)=H(Q)>0: use one run of the root module in a Hamilton word of the quotient. Consequently the contrast of the actual physical query is

    sum_(h>=1) c_h(Q,r) D_h(q)
      <= H(Q) D_1(q)
      = -H(Q)P(q) < 0.                                (5)

This proves(1) for every fixed ambient tournament, including arbitrary strong ambient contexts. In particular the literal gap has magnitude at least H(Q)P(q).

## 4. Scope and relation to earlier boundaries

The theorem requires a common complete child profile and common child order at most seven. It does not assert coordinatewise domination for arbitrary tournament profiles of larger order or for unrelated child profiles. The order bound comes from a finite exact continuous-domain polynomial certificate, not a census of all order-seven tournaments.

Unlike an even-core modulo-eight argument, this proof also covers odd n, and it can detect a nonzero integer gap when the chosen query is zero modulo eight. For example equal transitive order-four children have a nonzero core gap-407568 divisible by eight; inequality(5) nevertheless forces a strict gap in every fixed ambient.

All changes here are inside the distinguished module, and all target children retain the source common profile. Arbitrary overlapping operations involving exterior names, arbitrary changes to the ambient quotient, unrelated quotient replacements, and unrestricted named-response reconstruction are not claimed. An exact response collision is not produced.

## 5. Reproduction

The verifier recomputes every displayed D_h coordinate on explicit transitive profiles through order seven, plus cyclic and nontransitive controls, and checks the complete positive-square scalar formula. Its separate C++ subset-DP checker retains actual physical source/target counts, old adjacency rows, extended masks and the predicted quantitative bound in several fixed outer contexts. Those controls supplement the proof; arbitrary ambient size follows from(5), not from sampling.
