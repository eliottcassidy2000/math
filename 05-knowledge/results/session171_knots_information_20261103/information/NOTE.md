# A hidden holonomy memory with a one-vertex Hamilton read head

**Status: proved.** The construction combines S170's Reed–Muller holonomy identification and constant-profile tournament seed compiler. Its new feature is a compact, distance-preserving tournament storage medium and an explicit contextual Hamilton-count readout. This is not an autonomous dynamics theorem or a classification of all Hamilton counts. Coding facts are inherited, not claimed new in coding theory.

## 1. The distinction that makes the bridge useful

An invariant of an object and the responses of that object to controlled extensions are different observers. A whole complete path-cover profile can be blind to exponentially many stored states while a one-vertex extension reads one specified state bit. Information loss is therefore relative to the allowed interface, not an absolute feature of the tournament.

Put d=n−1≥2, M=2^d, and 2≤r≤d. S170 identifies the promised depth-r signed transport returns with

    C_r = RM(d−r,d) = {h:F_2^d→F_2 of Boolean degree ≤d−r}.
    D_h(s,y,t)=(s+h(y),y,t).

Its dimension and minimum distance are

    K = sum_{j=0}^{d−r} binom(d,j),       delta=2^r.

S170 proves these returns are attainable by scheduled closed binary circuits in a sufficiently large fixed controller cube. No claim is made that every return is realizable by prescribed runner speeds. The present construction stores the resulting transformation D_h, not its word, timing, controller, or geometric origin.

## 2. A physical seed with two indistinguishable bulk profiles

Use vertices 0,...,4 and outgoing-neighbor bitmasks

    B_0=(16,1,3,7,14),       B_1=(20,1,2,7,14).

Thus B_0 and B_1 differ only by reversing arc {0,2}. Both are strong and nonisomorphic. They have the same unordered directed-path-cover polynomial

    P(z)=9z+33z²+31z³+10z⁴+z⁵.

For r labelled path slots, including empty slots, the corresponding counting polynomial is

    F(u)=sum_j P_j (u)_j = u⁵+6u³+2u,

where (u)_j=u(u−1)...(u−j+1). In particular F(1)=9 is the number of directed Hamilton paths. The symbol u here is a polynomial variable, not commutator depth r.

These finite seed facts are inherited from S170 and checked independently here. Strongness is checked by directed reachability; nonisomorphism by all 120 permutations; P by enumerating paths and cuts and removing ordering duplicates.

## 3. The constant-profile tournament code theorem

List y∈F_2^d in binary integer order. Define

    T_h = B_{h(0)} ⊕ B_{h(1)} ⊕ ... ⊕ B_{h(M−1)},

where ⊕ is ordinal sum: every arc between earlier and later blocks points forward.

**Theorem 1.** The family {T_h:h∈C_r} has the following properties.

1. Every member has 5M vertices, complete labelled-slot profile F(u)^M, and Hamilton count 9^M.
2. It has exactly 2^K pairwise nonisomorphic members. The codeword is recoverable from the ordered strong-component types without vertex labels.
3. On the fixed labelled vertex set, the number of arc reversals between T_h and T_g is exactly wt(h+g). Consequently its minimum labelled arc distance is 2^r.
4. From a received labelled tournament with at most e adversarial arc reversals and s erased unordered pairs, the original member is uniquely determined whenever

       2e+s < 2^r.

   This allows errors on any pairs, not merely the M designated data arcs. It is an information-theoretic decoding guarantee, not an efficient decoding claim.
5. For uniform h, the original complete profile has zero mutual information with h, and the hidden entropy is exactly K bits.

**Proof.** A path in an ordinal sum visits component blocks in increasing order. A cover in u labelled slots restricts to an independent u-slot cover in each block, and conversely these restrictions concatenate uniquely within each slot. Thus F multiplies. Hamilton paths are the u=1 case. Since every seed is strong, these are exactly the strong components. The strong components of a tournament have a unique linear order, so an isomorphism preserves the position and type of each block. Their two types differ, proving injectivity on isomorphism classes.

All interblock arcs and all but one intrablock arc per block are fixed. The remaining arc is precisely h(y), up to the fixed orientation convention; this proves the exact metric identity. Two candidate codewords consistent with e errors and s erasures can disagree on at most 2e+s pairs: outside erasures every disagreement must be charged to an error in at least one candidate. This contradicts the minimum distance. Finally the profile is constant on 2^K equally likely codewords. ∎

This is a finite binary code embedded isometrically in a tournament invariant fiber. The words “hidden memory” mean hidden from that specified profile observer, not hidden from the graph itself.

The table representation costs 5M vertices. By contrast, storing M full signed-sheet registers independently would use 5M(n+1) vertices in the inherited compiler. The saving is possible because the addresses y and the shared unchanged index action are external to the stored sign table. It is not a minimum-vertex theorem.

Examples:

    d=3,r=2: 40 vertices, 16 members, distance 4, H=9^8.
    d=4,r=2: 80 vertices, 2048 members, distance 4, H=9^16.
    d=4,r=3: 80 vertices, 32 members, distance 8, H=9^16.

Deeper holonomy thus means fewer stored messages and greater protection against arc corruption, while the entire bulk profile is unchanged.

## 4. A one-vertex context reads one bit

For a seed B_b, add vertex a which beats only old vertex 0. Write Q_b for the resulting six-vertex tournament. Its rows are

    Q_0=(16,33,35,39,46,1),
    Q_1=(20,33,34,39,46,1).

The Hamilton counts by starting vertex are respectively

    Q_0: (1,2,4,8,1,1), total 17;
    Q_1: (1,2,3,7,5,1), total 19.

These are exact finite identities. One way to verify each entry is the endpoint recurrence

    N(S,j)=sum_{i∈S\{j}, i→j} N(S\{j},i),
    N({j},j)=1,

with a fixed start retained; the checker independently uses the unrestricted endpoint recurrence, and the listed small rows also permit direct enumeration.

For a selected block y in T_h, place this new vertex within that block: all earlier blocks beat it; it beats all later blocks; within its own block it beats only the labelled vertex 0. Denote this one-vertex induced extension by Q_y(T_h).

**Theorem 2 (Hamilton read head).**

    H(Q_y(T_h)) = (17+2h(y))9^(M−1).

Consequently the bit is read exactly by

    h(y) = [H(Q_y(T_h))/9^(M−1) − 17]/2.

**Proof.** The extended tournament is the same ordinal sum with its y-th seed replaced by Q_{h(y)}. Apply Hamilton-count multiplication and the two seed counts. ∎

The context needs the chosen block and the local label 0; it uses the same attachment rule for both possible bit values. It does not first inspect which seed is present. This is an operational separation: profile-equivalence is not a congruence for these controlled one-vertex attachments.

**Corollary 3 (optimal local tomography).** Using this read head at the K addresses of Hamming weight at most d−r reconstructs the full holonomy. Exactly K local binary read-head queries are necessary in the worst case, even adaptively, for full reconstruction or deterministic identity certification.

For sufficiency, the coefficient of a monomial indexed by A is

    c_A = sum_{B⊆A} h(1_B) mod 2.

For necessity, each local read provides one binary linear evaluation of the K-dimensional code. For adaptive identity testing follow the all-zero transcript: fewer than K evaluations have a nonzero common kernel. Full reconstruction also has the elementary 2^K-message binary decision-tree lower bound. This is the same probe theorem as S170, now realized by physical Hamilton counts of one-vertex extensions.

The sibling S171 information-geometry packet further identifies the exact residual entropy K−rank(E) for a partially read probe bank E and the statistical directions it loses. The present packet supplies a concrete tournament implementation of that interface.

## 5. Why the observer bandwidth must be stated

Choose q different blocks S and add the local read vertex in each of them simultaneously. This is a q-vertex extension, not a one-vertex read. If w=sum_{y∈S}h(y) as an ordinary integer, then

    H(Q_S(T_h))=9^(M−q)17^(q−w)19^w.

For known q, different w give different counts. Thus one such scalar count can have q+1 outcomes and reveals subset Hamming weight, not just one bit. If S is the full cube, it reads the codeword weight; it still does not generally identify h because distinct codewords can have equal weight.

Therefore the K-query lower bound is **not** a lower bound for arbitrary extension contexts or arbitrary integer-valued Hamilton queries. Any claim of universal “information per probe” must specify the allowed probe and its response alphabet. A scalar integer can encode many bits.

Likewise the arc-error guarantee concerns decoding the received labelled adjacency structure. It does not assert that a corrupted graph's Hamilton read-head response differs in at most e output bits. A single erroneous arc can alter many path counts. Label uncertainty, lost block addresses, corrupted promises, and arbitrary graph isomorphism are different noise models.

## 6. Relation to resets and the boundary of the synthesis

S170's guard language shows that Hamilton specialization can collapse an endpoint transfer to rank one, even though the full polynomial profile remains injective. Here the loss is complementary: even a complete unpointed profile can collapse a whole code family, while a controlled pointed extension recovers its state. Together these say that neither “Hamilton count” nor “complete profile” alone defines an absolute information ordering without the allowed contexts.

The exact synthesis map is

    closed scheduled circuit → diagonal transport h → tournament T_h
                              → local contextual Hamilton response.

The first map forgets words, clocks and geometry. The second preserves codeword identity and Hamming distance as labelled arc distance, but forgets the transformation's origin. The third reads one specified coordinate. No map here makes the physical schedule intrinsic or proves a new runner existence theorem. The construction is generally not strongly connected because it deliberately uses ordered strong components as an address space.

A concrete next question is whether one can retain both the distance guarantees and a constant complete profile in a uniformly strongly connected storage family with comparably cheap local Hamilton read heads. This packet does not resolve that stronger realization problem.

## 7. Reproduction

    python3 check.py > certificate.json
    python3 -O check.py > certificate.optimized.json
    cmp certificate.json certificate.optimized.json

The standalone checker imports no repository files, uses explicit failure checks, and verifies seed profiles and nonisomorphism, literal two-block contexts, all small codewords where K≤11, optimal probe banks through d=5, literal labelled-arc isometry, and all 12,480 single-arc corruptions of the sixteen 40-vertex d=3,r=2 code tournaments. Every corruption has the correct unique nearest codeword. Normal and optimized outputs agree.

SHA-256:

    check.py: e99e459f54eec5ac04208ad941ec283dca82d0a2997cc4327c36411e846d1182
    certificate.json: 26bcabfc1c6a09f363a16c3179dc3becad6876ea6ca9904e11ad609704f6d91a

Inherited source: the S170 holonomy packet linked by the session index at https://github.com/eliottcassidy2000/math/blob/2d137c00eaf36a6fbe9f06044ac56d2c3bf7c7d6/00-navigation/S170.md.
