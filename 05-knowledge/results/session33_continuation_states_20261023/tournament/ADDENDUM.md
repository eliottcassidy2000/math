# Triangle lattice and Hall obstructions for module observers

Status: PROVED elementary consequences and explicit hostile. No priority claim. Uses the definitions and main theorem in BOUNDARY_OBSERVABILITY.md. The main files were left unchanged while this addendum was audited. The statements concern directed path counts, not a general recovery of tournament isomorphism.

## A. The high path-cover coordinate is a triangle census

Let S have order n>=3, c3(S) directed triangles, and F_S(d)=d! pc_S(d). Then

    pc_S(n)=1,
    pc_S(n-1)=binom(n,2),
    pc_S(n-2)=binom(n,3)+2c3(S)+3binom(n,4).                 (A1)

Proof. An (n-2)-cover has either one three-vertex directed path or two two-vertex directed paths, all other components singleton. A transitive triple has one spanning directed path, whereas a directed triangle has three. This contributes binom(n,3)+2c3. On four chosen vertices each of the three pairings gives two directed edges, uniquely oriented in the tournament. This contributes 3binom(n,4). These cases are disjoint and exhaustive.

For the inherited centered polynomial

    E_S(u)=sum_(d=1)^n F_S(d)(u-1/2)^(d-1),

we obtain

    [u^(n-3)] E_S(u)
      =2(n-2)! c3(S) - n!(n-2)/12.                         (A2)

Indeed the coefficient is

    F(n-2) -(n-2)F(n-1)/2 + binom(n-1,2)F(n)/4;

substitute (A1) and simplify. The expression may be rational, as the centered polynomial itself may have rational coefficients. The coefficient is affine in the triangle census, with exact lattice spacing 2(n-2)!.

### Necessary form of a hypothetical second-depth twin at order seven

Suppose S,T have order seven and F_S(1)=F_T(1), F_S(2)=F_T(2). The leading u^6 coefficient cancels, their difference is even, and the endpoint-jet identity forces a double root at u=1/2. Evenness forces a double root at -1/2 too. Therefore

    E_S(u)-E_T(u)
      =240(c3(S)-c3(T))(u^2-1/4)^2.                        (A3)

This is a necessary identity, NOT existence or nonexistence of such a pair. If the profiles differ then their triangle counts differ. The sampled order-seven negative does not establish that (A3) cannot occur.

One can make its observable consequence explicit. Since E''(1/2)=2 F(3), equation (A3) gives

    F_S(3)-F_T(3)=240(c3(S)-c3(T)).                         (A4)

The four-exterior-vertex probe has R_2=11F(1)+17F(2)+4F(3), so a nonzero hypothetical difference would be

    R_2(S)-R_2(T)=960(c3(S)-c3(T)).                         (A5)

These normalize all quantities consistently: F is ordered-path-cover count, not pc. Equations (A3)-(A5) constrain an unresolved deeper hostile; they do not manufacture one.

## B. Hall coverage is exactly the top observer's support

For a context with r outgoing and r incoming neighbors of the module, let B be its zero-one outgoing-to-incoming adjacency matrix. The main theorem proves

    a_Q(r+1)=r! permanent(B).                              (B1)

Thus this coefficient is positive exactly when the bipartite graph has a perfect matching, equivalently when its coverage function obeys Hall's inequalities

    |N(X)|>=|X| for every X subset of the outgoing shore.   (B2)

This is a literal bridge to donor-mask/coverage data: the observer's first possible deep coefficient counts ordered systems of distinct representatives. Singleton degrees and total margins do not retain the subset-union information needed by (B2). No other clock predicate or number-theoretic theorem transfers automatically.

## C. Same labelled tournament scores, absent versus present top coefficient

Take r=3 and the matrices

    B_bad  = [[1,0,0], [1,0,0], [0,1,1]],
    B_good = [[1,0,0], [0,1,0], [1,0,1]].

Both have labelled row sums (1,1,2) and labelled column sums (2,1,1). But

    permanent(B_bad)=0,
    permanent(B_good)=1.

The first two rows of B_bad jointly cover only the first column, violating Hall at a two-element set although each singleton is covered. B_good has the unique diagonal perfect matching.

Construct quotient tournaments Q_bad and Q_good on

    *, o1,o2,o3, i1,i2,i3.

Orient *->oj and ik->*; orient each shore transitively in the displayed order; orient oj->ik exactly when the corresponding matrix entry is 1. Both quotient tournaments have the SAME labelled outdegree vector

    (3,3,2,2,4,4,3).

Indeed each outgoing vertex sees its fixed internal outdegree plus its row sum, and each incoming vertex sees its fixed internal outdegree, the module, and three minus its column sum. The module has outdegree three in both.

Nevertheless their top response coefficients differ:

    a_(Q_bad)(4)=0,
    a_(Q_good)(4)=6.                                      (B3)

The script independently counts their entire directed quotient-word coefficients, not merely the permanents, and confirms (B3). Their full labelled score data therefore do not determine even support of the top module-observer coefficient. They also have the same directed-triangle census, because for a tournament q the identity

    c3=binom(q,3)-sum_v binom(outdegree(v),2)

counts each transitive triple exactly once at its source. Thus the triangle census in (A2), though a real path-cover coordinate, does not restore this missing interface information.

This is a sharp small obstruction for the specific row/column-margin compression: with r=1 there is no ambiguity. For r=2, a zero-one matrix has positive permanent exactly when it has no zero row or column; in that case its permanent is two only for the all-one matrix, and otherwise one. The labelled margins distinguish these cases. Hence r=3 is the smallest shore order with same row/column margins but different permanent. We do NOT claim this is the smallest arbitrary quotient-score obstruction under other context conventions.

## D. Exact replay

    python3 tournament/check_addendum.py > tournament/addendum.json
    python3 -O tournament/check_addendum.py > tournament/addendum-optimized.json

The transcripts agree byte for byte. The checker verifies (A1)-(A2) on all 1,096 labelled tournaments of orders three through five, the two seven-vertex quotient word kernels and scores, and all 18 zero-one incidence matrices at shore orders one and two. The exact certificates contain script hashes. The elementary proofs above are all-size where quantified; these finite tests are independent controls, not replacements for the proofs.
