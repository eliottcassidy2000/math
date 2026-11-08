# Graft determinant, deletion-ratio obstruction, and scalar observation

**Status: PROVED, with independent FINITE-EXACT controls.** No claim of literature priority, arbitrary nonlinear state minimality, or Hamilton-count coverage is made.

## Inheritance and scope

The dependency is the all-order rooted graft closure in `05-knowledge/results/session225_graft_recovery_20261108/PROOF.md`, Sections 1–3. Its source state is `x=(D,E,U,d)`: total Hamilton paths, paths ending at the root, ordered two-path covers with the second nonempty path starting at the root, and root-deleted Hamilton paths. Empty unmarked path bins count once. Grafting nonempty rooted B onto nonempty rooted A directs B to the old root a, all old nonroots to B, and marks B's root r.

The closest prior mechanism is that general four-state matrix. The hostile is the strong directed triangle: strong connectivity will not imply invertibility. The corrected near miss is inferring scalar coverage or a common low-dimensional observer from command rank. The least-used sidecar is the ratio of two-bin to one-bin path counts after deleting the gadget root.

The concept board here comprises determinant factorization, root-deletion ratios, positive split counts, and scalar context observation. The main new structural conclusion is the factorization and rank/image criterion. Scalar observation is a language-specific refinement of existing work, not a first four-dimensional-minimum claim: `05-knowledge/results/adaptive_state_induction_20261103/BLOCK_SWITCH.md`, B2, already proves a rational four-dimensional minimum for its R/N/J1 language. Its conventions reverse orientations and use different coordinate ordering. We keep Session225's conventions throughout.

## 1. Exact determinant factorization

Use the inherited gadget statistics `(H,Q,s,e,alpha,delta,T,q,b)`. Thus Q and b count ordered two-bin covers of B and B-r; H and q count Hamilton paths there; s,e are root-start and root-end Hamilton counts; alpha,delta count rooted two-bin covers; T counts rooted three-bin covers. For |B|>=2 put

    Delta = q Q - H b.

Then the graft matrix satisfies

    det K_B = s e Delta.                                      (1)

**Proof.** Change source coordinates to `(p,E,t,d)=(D-E,E,U-d,d)`. The change has determinant 1. The output matrix, still in `(D',E',U',d')` coordinates, becomes

    [H      0   Q-H       H       ]
    [e      0   delta-e   0       ]
    [alpha  s   2T-alpha  2alpha-s]
    [q      0   b-q       q       ].

Expansion along the second column gives minus s times the determinant of

    [H  Q-H      H]
    [e  delta-e  0]
    [q  b-q      q].

The latter determinant is e(Hb-qQ), proving (1). In particular alpha, delta, and T cancel from the determinant even though they affect the operator.

The singleton gadget is exceptional: its inherited J has rank 3 and determinant 0. Formula (1) is not asserted there because its last row is different.

## 2. Strongness, exact rank, and smallest full-rank gadget order

For every nonempty tournament B, s>0 iff r reaches all vertices. To prove the nontrivial direction, take a longest directed path beginning at r. If an outside vertex v is beaten by any path vertex, insert v immediately after the last such path vertex; this extends the path because v beats every later path vertex. Thus every outside vertex beats the entire maximal path. No directed path from r can then leave its vertex set, contradicting reachability unless it is Hamiltonian. Reversing all arcs gives the corresponding statement for e.

Consequently s e>0 iff B is strong. Therefore for |B|>=2:

- K_B is invertible iff B is strong and Delta is nonzero.
- If B is strong and Delta=0, K_B has rank exactly 3.
- No nonstrong B gives an invertible K_B.

For the rank claim, the second column of the displayed matrix is nonzero only in row 3. The remaining three rows have rank 2 when Delta=0: row 4 is q/H times row 1, whereas row 2 is independent of row 1 since its final entry is 0 and row 1's final entry H is positive. Their rank is 3 when Delta is nonzero. Adding the independent row/column involving s proves the dichotomy.

Order 4 is the smallest possible full-rank gadget order. Order 1 is the exceptional singular J, order 2 is nonstrong, and at order 3 the only strong tournament is the directed triangle. Its `(H,Q,q,b)=(3,12,1,4)` has Delta=0. The Session225 gadget `[8,1,3,6]` at order 4 has s=e=1 and Delta=-8, so attains full rank.

Strongness is not a sign condition: the rooted order-four gadget `[2,8,3,5]` has `(H,Q,q,b,s,e)=(5,32,3,12,1,2)`, hence Delta=36 and determinant 72. The strong order-five gadget `[18,8,3,5,14]` has `(H,Q,q,b,s,e)=(13,104,3,24,3,2)`, hence Delta=0 and rank 3. In particular the triangle is not the only strong singular gadget order. No infinite strong singular family is claimed.

## 3. Deletion-ratio geometry and an exact inverse

Set w=p+d. The first and last rows say

    D' = H w + (Q-H)t,
    d' = q w + (b-q)t.

Eliminating w yields the all-order image identity

    q D' - H d' = Delta (U-d).                              (2)

The source split count U-d is nonnegative. It is zero exactly when the old root a is a sink: if a has an outneighbor v, use path (a,v) in its rooted bin and any Hamilton path on the remaining vertices in the other bin; conversely a sink supports no nontrivial rooted path. Existence of a Hamilton path in a tournament follows by the usual insertion argument, with the empty tournament allowed.

Thus Delta records the direction in which the output deletion ratio d'/D' moves from q/H. For a nonsink source root it is strictly below q/H when Delta>0 and strictly above q/H when Delta<0. For a sink source root there is equality. When Delta=0 all source states satisfy the same exact output relation

    H d' = q D'.                                           (3)

For a strong singular gadget, (3) is exactly its one-dimensional linear information loss: its rational image is the entire three-dimensional hyperplane specified by (3), not a claim that every point of that hyperplane is a physical tournament state.

Equivalently,

    Delta = H q (Q/H - b/q).

So invertibility requires that deleting the root change the ratio of ordered two-bin covers to Hamilton paths, as well as preserving both endpoint access conditions. This is a ratio comparison, not an asserted probabilistic covariance formula.

For an invertible gadget and any output state `(D',E',U',d')`, recover the input rationally by

    t = (q D' - H d') / Delta,
    w = ((Q-H)d' - (b-q)D') / Delta,
    p = (E' - (delta-e)t) / e,
    d = w-p,
    E = (U' - alpha p - (2T-alpha)t - (2alpha-s)d) / s,
    D = p+E,   U = t+d.                                    (4)

These equations invert K on all rational vectors, not merely on physical ones. They do not assert that an arbitrary integral output has a physical predecessor. Chains of invertible grafts can therefore be reversed at the state level by applying (4) in reverse gadget order. This does not recover the underlying tournament from its state.

## 4. Four scalar observations, with a necessary second layer

Let J denote the singleton graft, and L any directed-edge graft. The choice of new root for L does not affect its total Hamilton count. For a source A record

    h0 = H(A),
    hJ = H(J(A)),
    hJJ = H(J(J(A))),
    hL = H(L(A)).

The inherited first-row formula and J matrix give the observation matrix

    [1  0  0   0]
    [1 -1  1   0]
    [1  0  2   0]
    [1 -1  3  -2],       determinant 4.                    (5)

All four source coordinates can be recovered as

    D=h0,
    U=(hJJ-h0)/2,
    E=h0+U-hJ,
    d=(hJ+2U-hL)/2.                                        (6)

Hence two rooted tournaments have equal four-state vectors iff every finite rooted-graft context produces the same final Hamilton count. The forward direction follows from inherited closure under each graft; the reverse direction uses these four contexts only. No statement about arbitrary root changes or arbitrary graph contexts is intended.

This also supplies a minimal rational *linear* observer statement. Physical states span four dimensions: singleton, source-rooted edge, rooted directed triangle, and source-rooted transitive triangle have state-row determinant -4. Their responses against (5) therefore give a matrix of determinant -16 and rank 4. Any model assigning each source a vector in Q^k with every context response a linear functional factors this response matrix through Q^k; hence k>=4. The inherited state achieves k=4. Arbitrary nonlinear or integer encodings are outside this assertion.

There is a sharper depth distinction. Every one-graft observation row has form

    H [1,-1,0,1] + (Q-H)[0,0,1,-1].

These rows span exactly two dimensions, as singleton and edge gadgets already show. Adding the empty context gives dimension exactly three. Thus no collection, even infinite, of one-layer graft observations plus the original count gives a full-rank rational linear observer on the source-state space. A nested depth-two context is necessary and (5) shows it sufficient. This concerns linear observation rank; it does not assert a collision between physical tournaments or exclude nonlinear recovery on a restricted source family.

## 5. Finite controls, illustrative bridge, and boundaries

The portable verifier includes an unchanged copy of Session225's verifier as `inherited_verify.py`; the new script imports its independently specified literal subset DP, direct permutation enumeration, and graft constructor. It checks all 1,099 labeled tournaments of orders 1 through 5 with root fixed at 0, without filtering. For every gadget it independently compares all nine statistics with direct permutation/split counts. For each non-singleton gadget it checks the factorization, rank conditions, and literal graft outputs on the four physical spanning sources. Every invertible case is inverted on those literal outputs. It also checks all four scalar observations literally on the entire 1,099-source bank, including source singleton and sink-root boundaries.

A secondary finite bridge uses the two Session224 H=9 birth seeds

    A=[14,28,24,16,1], B=[10,28,25,16,1].

Applying the fixed Session225 four-vertex gadget at each of their five possible old roots gives respectively

    A: [423,605,357,185,261],
    B: [315,541,551,217,207].

These sets are disjoint, so this particular inherited scalar collision is separated regardless of root choice. This is a finite illustration, not a universal graph-separation result or a new birth-language coverage theorem.

Reproduce in this directory:

    python3 -I -B verify.py > replay.json
    cmp replay.json CHECK.json
    python3 -I -B -O verify.py > replay.optimized.json
    cmp replay.optimized.json CHECK.json

The all-order conclusions come from the proofs above, not from the census. The all-odd Chan–Pak completeness question, any improved unrestricted scalar coverage, and nonlinear observer minima remain unproved here.
