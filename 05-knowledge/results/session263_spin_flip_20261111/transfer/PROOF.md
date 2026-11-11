# A free sign-flip source compiler and three-query state tomography

**Status: PROVED from the inherited THM-5130 transfer law.** The new finite controls test only the compiler, decoder, and tomography below. No inherited census or transfer proof is re-audited. Chan–Pak and universal odd-coefficient supply remain OPEN.

## Sources, inheritance, and discovery scope

Source pin: THM-5130, `01-canon/theorems/THM-5130-sign-flip-towers-gadget-transfer-law-and-seeds-of-both-signs.md`, file commit `7c14d6741eab74550c9337f6143751399d36053e`. Its parts A–C are accepted at their audited PROVED status. The closest prior construction is `05-knowledge/results/session262_interface_supply_20261111/INCOMING_CONNECTIONS.md`, file commit `83f0f4ee1eb2eaed4ab09d3f4f1caa74be933727`, which already gives F3,F5, paired flips, four first singleton states, and indefinite fixed-offset reuse. We do not claim those results. Its warning that reuse does not imply target control is the corrected near miss. The endpoint coordinate E is the least-used necessary sidecar. The hostile example is scalar slope/count projection, which forgets h and E.

Narrow prior search: the two source packets and Session260–262 packets for free semigroup, ping-pong, projective interval, word decoder, and return-state terms. No prior version of the disjoint-cone decoder was located. This is a bounded novelty search, not a global priority claim.

Concept board: physical moving marks / projective separation / odd slope supply / scalar tomography / target coverage. Reflection-versus-rotation and chirality snippets motivated looking at noncommuting products; they are not mathematical evidence. These matrices are not involutions or rotations: det(F3)=-4, det(F5)=-8.

## 1. The literal grammar and its guard

For a pair (B,O), use the inherited state s=(r,M,h,E), with r=H(B)-W_B(O), M=U_B(O), h=H(B), E=End_B(O). Write q=h+r. A letter a in {3,5} appends the transitive tournament of order d+1, where d=1 for a=3 and d=2 for a=5, closes it by a new vertex as in THM-5130, and moves O onto its first and last vertices. It adds a total of 3 or 4 vertices. The exact state law is

    r'=-r,
    M'=a q+(a²+1)M-E,
    q'=(2a-2)M,
    E'=aM+2r-q.                                      (1)

The admissible cone K consists of LITERAL selected tournaments with |O|>=2, r in {-1,1}, and 0<=q<=M/2. Their counting interpretations give

    h>=1, M>=2h, 0<=E<=h<=M/2.                       (2)

We use literal objects, not arbitrary integer points satisfying inequalities. A decoder can reject a tuple without determining general tournament realizability.

An explicit initial object is B=TT2 with BOTH vertices marked:

    s0=(-1,3,1,1), q0=0, order 2.                    (3)

The first letter states are s3=(1,29,11,7) and s5=(1,77,23,13). Thus the grammar starts with actual odd slopes, not a hypothetical interface.

## 2. Disjoint images and an all-length free action

**Theorem 1.** Both letters preserve K. Their images have disjoint rational intervals:

    a=3: 8/23 <= q'/M' <= 8/19;
    a=5: 16/57 <= q'/M' <= 16/51.                   (4)

For either letter, 19M/2 <= M' <= 57M/2. For a=3 specifically the multiplier lies in [19/2,23/2]; for a=5 it lies in [51/2,57/2].

*Proof.* From 0<=q<=M/2 and 0<=E<=M/2, (1) bounds M'/M between a²+1/2 and a²+1+a/2. Dividing q'/M=(2a-2) by these bounds gives (4). Both intervals lie in [0,1/2], and 16/51<8/23 (368<408), so they are disjoint. Literal construction preserves |O|=2; (1) flips unit offset and (2) remains valid by the counting meanings. This proves invariance and all inequalities. No E=0 reset is used. ∎

**Theorem 2 (unique decoder).** Fix ANY initial literal state in K. Different words in the two letters give different full states. They even give different triples (r,M,q,E), which are just coordinates for the same state. There is an exact linear-time-in-word-length decoder once the initial state is supplied.

*Proof and algorithm.* Given terminal (r',M',q',E'), the ratio in (4) uniquely determines the last letter, unless the state is the initial state and decoding stops. With c=2a-2 recover

    r=-r', M=q'/c, q=aM+2r-E',
    E=a q+(a²+1)M-M'.                               (5)

These formulas invert (1) exactly. Reject a nonintegral M or a failed cone/counting guard. Iterate. M decreases by at least a factor 19/2 at each inverse step. Thus no nonempty word can fix the initial state. If two words had equal outputs, (4) makes their last letters equal; (5) cancels them. Cancellation ends at equal empty words, since a remaining nonempty word cannot return to the initial M. Conversely, a tuple decoding to the initial state with all exact forward equalities is realized by that decoded literal construction. ∎

The last converse is a complete membership test for this INITIALIZED STATE ORBIT, not all literal tournament states. Guards alone are not sufficient. Bit complexity with elementary integer arithmetic is polynomial in terminal bit length: word length is at most log_(19/2)(M/M0), and every intermediate integer has O(log M) bits on K.

**Corollary 3 (legal paired-flip compiler).** Group successive letters into 33,35,53,55, with letters listed in order of application. These four macros act freely on every initial state in K and preserve r. They give all-length, exactly decodable fixed-offset physical source families, with moved selected vertices and no consumed finite reserve of trade slots. This sharpens the repeatability result of Session262 by proving absence of word collisions and exact return-word membership. It does NOT return to the original marked vertices or restore the original state.

Parity is also preserved in a useful form. If M and E are odd initially, then they remain odd: q=h+r is even because every tournament has odd H, the first two terms of M' are even, and E'=aM+r-h is odd. The seed (3) has this property. A word of length k gives

    r=(-1)^(k+1), M odd, order 2+3k+#(letters 5),
    3(19/2)^k <= M <= 3(57/2)^k.                   (6)

For every t>=1 its literal closing has

    H=(-1)^(k+1)+2^t M, order 3+3k+#5+t.           (7)

Thus even lengths give anti-seeds and odd lengths seeds, without odd-part normalization or a missing finite head. As illustrative homogeneous branches, M for 3^k begins 3,29,319,3461,37599; for 5^k it begins 3,77,2109,57551,1570759. Their recurrence/generating-function descriptions follow from (1), not new scalar digit maps. There are 2^k distinct FULL STATES at length k; (6) alone does not imply 2^k distinct slopes.

## 3. Three physical scalar queries recover the complete source

**Theorem 4 (tomography with a known sign).** Suppose r is known (for example from word parity). Query the following three actual Hamiltonian path counts, all with a one-vertex transitive closing tail:

    A=H(T1(B,O))=r+2M,
    B=H(T1(F3(B,O)))=-r+2(3q+10M-E),
    C=H(T1(F5(B,O)))=-r+2(5q+26M-E).               (8)

Then recover

    M=(A-r)/2,
    q=(C-B-32M)/4,
    E=3q+10M-(B+r)/2,
    h=q-r.                                         (9)

Hence the three-query observable distinguishes every pair of words of a fixed length and determines its exact word by Theorem 2. Queries B,C are distinct constructions from the SAME source, not simultaneous contradictory choices. Their orders exceed the source by 5 and 6, respectively; A exceeds it by 2.

*Proof.* Subtract the final two equations of (8), and solve. The determinant in coordinates (M,q,E) is nonzero (absolute value 16). All quantities are finite-head counts, and no limit or asymptotic query is substituted for them. ∎

If r is not known, four queries suffice: the count r+2M on s,F3s,F3²s,F3³s. Their matrix on (r,M,h,E) is

    [1,    2,   0,    0]
    [5,   20,   6,   -2]
    [59, 218,  62,  -20]
    [637,2368,674, -218],

with determinant 16. Thus parity information can be replaced by one further physical query.

**Corollary 5 (honest scalar complexity bound).** Let D_k be the union of the THREE count sets in (8), for all 2^k words from (3) of length k. Then

    |D_k| >= ceil(2^(k/3)),
    max D_k <= 1+171(57/2)^k,

and every member is realized by a tournament of at most 4k+8 vertices.

*Proof.* The map from words to ordered triples in D_k³ is injective by (9) and Theorem 2, so 2^k<=|D_k|³. Equation (6), the largest one-letter slope multiplier 57/2, and the closing law give the maximum bound. The order bound follows from k additions of at most four vertices and the largest query overhead six. ∎

This is a guaranteed scalar diversity bound for THIS narrow compiler, not a claim of a new global spectrum-density record. In cutoff form it gives a polynomial lower bound along X_k=1+171(57/2)^k with exponent log(2)/(3 log(57/2)); it gives neither positive density nor arbitrary target reachability. The 2^k-state claim MUST NOT be promoted to 2^k distinct counts.

## 4. Failure boundaries and next obligation

The source-to-state map preserves the exact four-coordinate interface; passage to q/M loses scale, E, and sign, and is used only to classify a last letter, not as a complete state. Inversion restores the necessary sidecars from the full state or three counts. A scalar slope is not asserted to determine a word. Word freeness is arithmetic state freeness, not pairwise tournament nonisomorphism after forgetting marks (equal underlying tournaments could have different selected sets).

For a concrete warning about missing state coordinates, (r,M,h,E)=(1,10,3,1) and (1,10,3,3) both satisfy the elementary numeric inequalities of K and have equal first closing counts but different F3 slopes (111 versus109). We DO NOT assert these two tuples are physically realizable. The example shows that the proposed *inequality-only algebraic* scalar quotient does not close; literal scalar collisions are not needed for the proofs above. The executable controls also test that a forged tuple with the same ratio cannot bypass exact inverse/end-state checks.

The genuine next target is a physical grammar with an all-target scalar cover: a finite menu of actual target-count operations with inverse residues and smaller realized predecessors, including heads. The present exact reverse STATE decoder does not supply such a scalar cover. It establishes lossless information transport, odd-slope source supply, and quantitative diversity within one legal compiler. Chan–Pak remains open.

## Reproduction

Run `python3 -B 05-knowledge/results/session263_spin_flip_20261111/transfer/controls.py` from repository root. The receipt records exhaustive new-word checks, exact integer interval checks, reverse-word uniqueness, physical subset-DP counts for new composites, tomography, and hostile forged-state rejection. The all-length conclusions rest on proofs above, not finite enumeration.
