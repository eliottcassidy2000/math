# Strong full-profile-blind memory with one-vertex Hamilton reads

**Status: proved.** A single fixed feedback vertex upgrades the S171 ordinal storage to a uniformly strongly connected family. It preserves the complete-profile blindness, isomorphism separation, labelled arc-distance guarantees, and one-vertex coordinate readouts. The construction has a canonical strong articulation vertex, so it does not settle the stronger 2-vertex-strong storage question. The resulting additive-exponential Hamilton counts escape the storage family's pure powers of nine, but are not asserted to be previously unrealized odd integers.

## 1. Physical wrapper and its canonical recovery

Let A_1,...,A_M be nonempty tournaments, M≥2. Take their ordinal sum, meaning all interblock edges point from earlier to later blocks. Add a vertex a with

    a → A_1,        A_i → a for i≥2.

Call the result W(A_1,...,A_M). Set A=A_1 and L=A_2⊕...⊕A_M. The three groups have uniform cyclic orientation

    A → L → {a} → A.

Thus W is strong, regardless of the internal orientations of the blocks: every vertex can reach every other by going around this three-group cycle. If |A|,|L|≥2, deleting any old vertex still leaves all three groups nonempty and therefore leaves a strong tournament. Deleting a leaves A⊕L, which is not strong. Consequently **a is the unique vertex whose deletion destroys strong connectivity**.

If each A_i is strong, deleting this canonically identifiable a recovers the ordered strong-component sequence A_1,...,A_M. In particular, nonisomorphic ordered seed sequences give nonisomorphic wrapped tournaments. This needs no unproved claim about a modular decomposition.

## 2. An exact complete-profile formula

For a tournament X, let F_X(u) count covers of all vertices by u labelled directed path slots, allowing empty slots. It is a polynomial, expressed by the unordered path-cover coefficients in the falling-factorial basis. Put

    f_i=F_(A_i)(u),    l_i=F_(A_i)(u−1),    g_i=F_(A_i)(u+1).

Then the full profile of the wrapper is

    F_W(u)=u { (l_1−f_1) ∏_(i≥2) l_i
               +(3f_1−g_1−l_1) ∏_(i≥2) f_i
               +(g_1−f_1) ∏_(i≥2) g_i }.                 (1)

**Proof.** For positive integer u choose which of the u slots contains a. Split that slot into a prefix P before a and suffix Q after a. Other slots retain their labels. A nonempty P must end in L, and a nonempty Q must start in A_1. The count is obtained from the following three cases and one subtraction, per choice of the slot:

- P=Q=empty: l_1∏l_i.
- Q=empty, P nonempty in L: ∏f_i−f_1∏l_i.
- Q nonempty in A_1, with initially unrestricted P: (g_1−f_1)∏g_i.
- Subtract the invalid assignments in the preceding case where P is nonempty only in A_1: (g_1−2f_1+l_1)∏f_i.

In all products without a first-block subscript the range is i≥2. Each original slot restricts independently to blocks and concatenates in block order. Summing these terms gives (1). Polynomial identity follows from agreement at every positive integer u. ∎

This proves that equality of the full block profiles is preserved by this particular strong wrapper. It is a controlled operation, not a statement that full-profile equivalence is preserved by every context.

At u=1, every nonempty block has F(0)=0. If h_i=H(A_i) and t_i=F_(A_i)(2), formula (1) becomes

    H(W)=(t_1−h_1) ∏_(i≥2)t_i +(3h_1−t_1)∏_(i≥2)h_i. (2)

For two blocks A,B this is the symmetric bilinear arithmetic interface

    H(W(A,B))=t_A t_B−h_A t_B−t_A h_B+3h_A h_B.       (3)

It uses the two-slot count as well as H. Replacing the interface by H alone is not justified.

## 3. A strong constant-profile code family

Use S171's five-vertex strong seeds, with outgoing-neighbor masks

    B_0=(16,1,3,7,14),       B_1=(20,1,2,7,14).

They are nonisomorphic, differ on exactly one labelled arc, and share

    F(u)=u^5+6u^3+2u,       H=9,       F(2)=84.

For any binary word h of length M≥2, put

    U_h=W(B_(h(0)),...,B_(h(M−1))).

**Theorem.** The 2^M tournaments U_h have 5M+1 vertices; they are strong and pairwise nonisomorphic; they have one common complete profile, given by (1) with every block profile F; and

    H(U_h)=75·84^(M−1)−57·9^(M−1).                     (4)

On the fixed labelled vertex set,

    arc_distance(U_h,U_g)=Hamming_distance(h,g).       (5)

**Proof.** Strongness and canonical recovery are Section 1. Formula (1) gives full-profile constancy, and (2) gives (4). All wrapper and interblock arcs are fixed; only one data arc per differing block changes, giving (5). ∎

In particular, use S171's holonomy code C_r=RM(d−r,d), with M=2^d, d≥2, 2≤r≤d. Then K=Σ_(j=0)^(d−r) binom(d,j), and the subfamily indexed by this code has 2^K members, minimum labelled arc distance 2^r, and unique decoding under 2e+s<2^r for e adversarial arc errors and s erased pairs. This follows from the inherited code distance and (5); the present theorem does not claim efficient decoding. The common profile loses all K bits of uniform source information.

The upgrade costs exactly one vertex beyond the ordinal code. Its storage rate is K/(5M+1) bits per vertex, compared with K/(5M) before wrapping. A useful structural resource measure is deletion resilience: these U_h survive every one-vertex deletion except the unique wrapper a. They are therefore strong but not 2-vertex-strong.

## 4. The read head still costs one vertex

For a seed B_b, S171 adds a vertex z beating only local vertex 0, giving Q_b. Direct subset enumeration verifies

    (H(Q_0),F_(Q_0)(2))=(17,224),
    (H(Q_1),F_(Q_1)(2))=(19,256).                      (6)

To read block y, replace that block by Q_(h(y)) using the same local extension rule as S171. Keep the old wrapper a: it beats z if y=0 and is beaten by z otherwise. Interblock orientations remain ordinal. Thus this is a fixed one-vertex induced extension of U_h, without first learning the bit.

If y=0, its Hamilton response is

    R_0(b)=(207+30b)84^(M−1)−(173+26b)9^(M−1).        (7)

If y≥1, the response is

    R_y(b)=(16800+2400b)84^(M−2)
                         −(969+114b)9^(M−2).          (8)

These formulas follow immediately by inserting (6) into (2). The integer response gaps are

    Δ_0=30·84^(M−1)−26·9^(M−1)>0,
    Δ_y=2400·84^(M−2)−114·9^(M−2)>0.                 (9)

Hence b=(R_y−R_y(0))/Δ_y recovers the bit exactly. Every queried tournament remains strong.

For the code subfamily, the same K evaluation addresses of Hamming weight at most d−r recover the codeword. Fewer than K adaptive binary coordinate queries cannot universally reconstruct it or certify the zero word. This statement concerns only this declared two-outcome read head, not arbitrary integer-valued contexts.

An additional interface metric is the **count-read separation** Δ_y. If a measured Hamilton count has adversarial integer additive error ε with 2|ε|<Δ_y, nearest-response decoding still identifies the bit. This is a count-measurement noise guarantee, not an arc-error guarantee: one bad arc need not cause a small count perturbation. The distinction prevents transferring the adjacency code's noise model into its observer without proof.

## 5. Arithmetic gain and its exact boundary

The first two block lengths give

    M=2: order 11, H=5787=9·643;
    M=3: order 16, H=524583=27·19429.

These are strong tournaments, whereas the unwrapped memory had H=9^M. Formula (4) is an actual subtraction of two exponential products with a positive combinatorial interpretation, not just another ordinal multiplicative family. Its exact 3-adic valuation is

    v_3(H(U_h))=M,

because H=3^M[25·28^(M−1)−19·3^(M−1)] and the bracket is 1 modulo 3. All these counts are odd. None of this proves a missing-prime-lane theorem or the full odd spectrum, and the numbers are not claimed globally new.

The guard/reset framework and this wrapper use different interfaces. Guard resets collapse endpoint channels; the present construction carries (H,F(2)) through a strong feedback attachment. Bulk memory is still invisible, but a specified one-vertex context can access the original structure. A conceptual connection to reverse arithmetic systems such as Collatz is that a scalar projection can hide the state needed for legal continuations. There is no Collatz conjugacy, orbit theorem, or inherited return-certificate theorem here.

Concrete next target: remove the canonical articulation while retaining an audited full-profile identity and bounded-cost bit read. Two feedback vertices are a natural candidate, but their exact profile and isomorphism recovery must be proved rather than inferred from this one-vertex result.

## 6. Independent verification

Run:

    python3 check.py > certificate.json
    python3 -O check.py > certificate.optimized.json
    cmp certificate.json certificate.optimized.json

The standalone checker uses literal adjacency subset dynamic programming. It computes every unordered path-cover coefficient by a canonical least-vertex partition recurrence, not by the wrapper formula. It checks all four two-block storage words, and both possible one-head addresses on each, for 12 literal graphs of orders 11 and 12. For each graph it verifies every profile value u=0,...,order against (1), strongness, and uniqueness of the deletion articulation. It also checks the seed and query profiles, all labelled arc distances in the two-block family, and response formulas (7)–(8). Polynomial values at order+1 distinct arguments suffice to verify the degree-at-most-order identity on these literal controls. General formulas are proved above, not extrapolated from enumeration.

Normal and optimized output agree. SHA-256:

    check.py: 6e59b5c4e66a8bd7d1cb4842d407423e3ed9efa2f417626e921ddbd94a1b8a77
    certificate.json: e78d47d59e5a06a78e68e2d07b183822626e14b4c090fa61a535bf6eaed8ff83

Inherited sources: [S171 information NOTE](../../session171_knots_information_20261103/information/NOTE.md) sections 2–5 for the seeds, code embedding and old read head; [S171 synthesis](../../session171_knots_information_20261103/synthesis/NOTE.md) section 5 for the explicit unresolved strong-storage target; [S170 guards NOTE](../../session170_holonomy_homotopy_20261103/guards/NOTE.md) for the distinct reset interface and its limitations. No repository-wide novelty audit is claimed.
