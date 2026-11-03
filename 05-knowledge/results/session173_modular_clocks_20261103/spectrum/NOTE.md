# Arbitrarily resilient full-profile-blind tournament storage

**Status: proved.** Transitive feedback modules remove the single-articulation obstruction from S172. A particularly explicit two-feedback family has 5M+2 vertices, survives every single-vertex deletion, retains all M hidden bits, and has a one-vertex Hamilton read head. More generally, any prescribed vertex resilience r is attainable provided the two old sides each have more than r vertices. These are constructive families, not a solution of the missing Hamilton-count spectrum {7,21}.

**Resilience scope.** Vertex-connectivity resilience means that strong connectivity survives deletions, not that the stored bits survive vertex erasure. The two seeds differ only on the arc between local vertices 0 and 2. Deleting either endpoint in an addressed seed makes its two bit versions identical, so even one vertex deletion can erase that bit. Canonical bit recovery is asserted for the intact storage family.

## 1. A substitution lemma for the complete path-cover profile

For a nonempty tournament X, write p_X(k) for its number of unordered covers by k nonempty directed paths, and

    F_X(u)=sum_k p_X(k)(u)_k,      q_X(k)=k!p_X(k).

Here (u)_k is the falling factorial; q_X(k) counts covers assigned bijectively to k ordered positions. The complete profile determines all p_X(k), including the order of X.

**Lemma (uniform-module substitution).** Fix a tournament D on module positions. Replace every vertex i of D by a nonempty tournament X_i, with all intermodule arcs directed according to D. The resulting complete profile depends only on D and the individual complete profiles F_(X_i).

**Proof.** For a positive integer u, a cover in u labelled slots consists of a directed path, possibly empty, in each slot. Collapse every maximal consecutive run from a single module to its module name. This leaves a uniquely determined run skeleton: a word in each labelled slot, with no equal consecutive module names and with allowed consecutive names directed by D. If module i occurs k_i times, its occurrences are distinguished by slot index and position within the word. They can be filled in exactly q_(X_i)(k_i) ways. The contributions from different modules multiply. Sum over all skeletons, restricting k_i≤|X_i|. This restriction makes the sum finite, and every original cover occurs exactly once. Equality for all positive integers u proves equality of the polynomials. ∎

This statement concerns uniform modules, not arbitrary vertex-by-vertex attachments. In particular, it does not contradict the existence of distinguishing read heads. It establishes a full-profile congruence for every fixed tournament substitution context, beyond the particular S172 singleton wrapper.

## 2. Resilient storage and canonical recovery

Use the verified five-vertex seeds

    B_0=(16,1,3,7,14),     B_1=(20,1,2,7,14),
    F(u)=u^5+6u^3+2u,     H=9.

The tuples are outgoing-neighbor bit masks. The seeds are strong, nonisomorphic, and differ on one labelled arc. The checker verifies strongness and nonisomorphism independently, exhausting all 120 vertex permutations; no new seed classification is assumed here.

For a binary word b of length M, split its ordered seeds after K positions:

    A=B_(b0)⊕...⊕B_(b(K−1)),
    L=B_(bK)⊕...⊕B_(b(M−1)).

Add a transitive tournament Z_r on r vertices, with uniform cyclic orientation

    A → L → Z_r → A.

Denote the resulting tournament V_b(r,K).

**Theorem.** Suppose r≥1, 1≤K<M, 5K>r and 5(M−K)>r. Then:

1. V_b has 5M+r vertices and remains strong after deletion of any fewer than r vertices.
2. Z_r is the unique r-element vertex set whose deletion destroys strong connectivity. In particular, vertex connectivity is exactly r.
3. The 2^M choices of b are pairwise nonisomorphic and have identical complete profiles.
4. On the declared labelled vertex set, arc distance is exactly Hamming distance between the bit words.

**Proof.** After removal of at most r vertices other than precisely Z_r, each of A,L,Z_r remains nonempty: the first two sides originally have more than r vertices. Any three nonempty tournaments with uniform cyclic intergroup orientation form a strong tournament, regardless of their interiors. Deleting Z_r instead leaves the nonstrong ordinal sum A⊕L. This proves (1)–(2). An isomorphism must preserve the unique minimum separating set Z_r. After deleting it, the ordered strong-component sequence is exactly the original sequence of five-vertex seeds. This recovers b and proves nonisomorphism. The substitution lemma proves profile equality. Only the one seed arc for each changed bit changes, proving (4). ∎

For any desired r, take K=floor(r/5)+1 and M≥2K. The feedback overhead is exactly r vertices. For r=2, K=1 works for every M≥2. Unlike S172, no vertex of this family is a strong articulation vertex. This does not claim maximal connectivity or an optimal vertex budget.

## 3. A positive Hamilton arithmetic formula for two feedback vertices

For arbitrary nonempty A,L and r=2, put x_i=q_A(i), y_i=q_L(i). Then

    H(V)=5x_1y_1+3x_2y_1+3x_1y_2+7x_2y_2
                  +2x_3y_2+2x_2y_3+2x_3y_3.             (1)

**Proof.** Collapse the Hamilton path into A-runs, L-runs and the two individually labelled feedback vertices z_1,z_2, with z_1→z_2. The allowed transitions are A→L, L→z_i, z_i→A, and z_1→z_2. Both feedback vertices must appear exactly once. There are at most three runs of each old side. The numbers of legal skeletons with (i,j) old-side runs are

    (1,1):5; (2,1):3; (1,2):3; (2,2):7;
    (3,2):2; (2,3):2; (3,3):2.

One may verify this finite list by beginning with the empty word and appending each permitted transition until both feedback vertices have appeared; ending is allowed precisely when both old sides have appeared. The accompanying checker exhausts this exact finite tree. Filling the ordered run positions contributes x_i y_j, proving (1). ∎

For a proof-ready closed equivalent, if f_X=(F_X(1),F_X(2),F_X(3))^T, then

    H(V)=f_A^T [[15,−11,2],[−11,13,−4],[2,−4,2]] f_L,    (2)

because (q_1,q_2,q_3)=(F(1),F(2)−2F(1),F(3)−3F(2)+3F(1)). Formula (1), rather than cancellation in (2), makes positivity and response monotonicity transparent.

For K=1 and the equal-profile seeds, F(1)=9,F(2)=84,F(3)=411. Ordinal sums multiply F, so for M≥2

    H(V_b(2,1))=33·9^(M−1)−651·84^(M−1)+504·411^(M−1). (3)

At M=2 this is 152757. Formula (3) is a three-exponential strong construction, distinct from S172's two-exponential singleton expression.

## 4. One-vertex reads survive every resilience level

Replace the seed at address y by Q_(b_y), obtained by adding a vertex which beats only the seed's local vertex 0. All external arcs follow that seed's module orientation, including its old side's orientation against Z_r. This is a fixed one-vertex induced extension, not an attachment chosen after learning the bit.

Literal seed enumeration gives

    p_Q0=(0,17,95,138,73,15,1),
    p_Q1=(0,19,109,150,75,15,1).

Thus

    F_Q1−F_Q0=2(u)_1+14(u)_2+12(u)_3+2(u)_4.            (4)

Products of polynomials with nonnegative falling-factorial coefficients again have nonnegative coefficients. Indeed,

    (u)_a(u)_b=sum_t binom(a,t)binom(b,t)t!(u)_(a+b−t),

where 0≤t≤min(a,b); count two injections by their t shared images. The remaining ordinal factors therefore preserve coefficientwise domination in (4). For the changed old side, all q_k differences are nonnegative and the q_1 difference is exactly 2·9^(side_length−1).

The Hamilton run-skeleton expansion for any r has nonnegative coefficients. In particular the skeleton A,L,z_1,...,z_r occurs and contributes q_A(1)q_L(1). Consequently the two bit responses have gap

    Δ_y ≥ 2·9^(M−1)>0,                                (5)

for every r,K satisfying the theorem. Exact response distinguishes the addressed bit; the other bits remain irrelevant by the substitution lemma. The extension retains resilience r and the unique separating r-set, by the same three-nonempty-groups argument. The count-error decoding guarantee is 2|error|<Δ_y; this does not translate into an arc-error bound.

For r=2,K=1 the exact response formulas are, with b=b_y:

    R_0(b)=(689+2b)9^(M−1)−(3071+254b)84^(M−1)
                                      +(2036+200b)411^(M−1),

    R_y(b)=(561+66b)9^(M−2)−(145824+20832b)84^(M−2)
                                      +(730296+81648b)411^(M−2),  y≥1.

At M=2 both addresses have responses 585033 and 645915, gap 60882. The positive-basis proof above establishes positivity uniformly, without extrapolating these examples.

## 5. Separate modular clocks and spectrum boundary

For (3), the exact 3-adic valuation is M. Divide by 3^M: the middle term contributes −217·28^(M−1), nonzero modulo 3; the first and third terms are divisible by 3 when M≥2. Thus all these counts are 0 modulo 9. Independently, modulo 8, they are 5 when M=2 and 1 for M≥3, since 9≡1,84≡4,411≡3 and 504≡0 modulo 8.

These two clocks must not be conflated: a residue modulo 8 does not undo the exact 3-divisibility obstruction. This family does not realize either missing candidate 7 or 21, does not produce a general prime-lane theorem, and does not justify a Collatz relation involving 527,44,5,8,13 or the proposed cycle {4,2,5,3,7,13}. The constructive advance is removal of the articulation obstruction with exact full-profile and readout guarantees, rather than a claim that a suggestive numerical cycle proves spectrum completeness.

## 6. Reproducibility and scope

Run `python3 check.py > certificate.json` and `python3 -O check.py > certificate.optimized.json`, then compare the outputs. `common.py` contains standalone literal subset-DP and least-vertex path-partition recurrences copied from the inspected S172 checker, not its wrapper formula. The checker verifies all 4 two-bit storage graphs and both addresses on every graph: 12 complete profiles on orders 12–13, exact Hamilton formulas, all one-vertex deletions, every separating pair, and labelled distances. It also checks three literal M=3 Hamilton counts (one storage graph and two read addresses), on orders 17–18. It separately checks r=1,3,4 resilience and unique minimum cuts on literal graphs; arithmetic clocks are tested for 2≤M≤30. General statements are proved above, not inferred from finite tests.

Sources inspected: `/tmp/s172_spectrum/NOTE.md` and its checker; `/tmp/s171_synthesis/NOTE.md`. The inherited nonisomorphic seeds can alternatively be checked exhaustively over their 120 vertex permutations; the checker includes this verification. No repository-wide novelty claim is made.

Final normal and optimized outputs agree. SHA-256:

    check.py: 0c6d335db04b87539c42caf3b688b8443a6fdc08bf5684bc68741b6bd4a5c78c
    common.py: 3f8f89bde8fa82c0d654c8189b9e74df3445f54d4178471e96bd721d97f287d7
    certificate.json: a2195267a87dcf7077d7821bb66036ef22706ab700ef2f12656b8d934c99d1b4
