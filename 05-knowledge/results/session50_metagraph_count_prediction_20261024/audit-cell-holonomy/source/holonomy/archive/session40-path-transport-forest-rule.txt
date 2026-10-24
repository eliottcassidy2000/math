# All-size actual-path transport and a sharp bounded-core wiring rule

23 October 2026. **PROVED.** These are sufficient rules with explicit retained data. They are not a bounded-number-of-marks theorem and do not classify arbitrary anonymous metagraph margins.

## 1. Exact directed refinement transport

Let A be a tournament and P=(p_0,...,p_(n−1)) an ACTUAL directed Hamiltonian path. Let B reverse exactly the arcs not used by P. Write P_uv=1 when u immediately precedes v along P. For distinct u,v,

    B_uv = 1−A_uv + P_uv−P_vu.

For every named vertex subset C,

    sum_(w in C) B_vw
      = |C|−1_(v in C)−sum_(w in C) A_vw
        +1_(successor_P(v) in C)−1_(predecessor_P(v) in C).

Absent predecessor/successor indicators are zero. Therefore target colour refinement can be computed iteratively from source outgoing counts into each current TARGET colour cell, plus the actual predecessor/successor owners along P. No target-class isomorphism oracle is needed. The formula works for any chosen marks, provided the same actual vertices and mark roles are retained.

In particular target degrees are n−d_A(start), n−2−d_A(end), and n−1−d_A(v) at interior vertices. Independent histograms of predecessor and successor colours are insufficient: the correction belongs to the SAME vertex whose source count is being updated. This is the exact shared-incidence requirement.

Converse with reversed P commutes with the construction. Transport every mark by its actual vertex and reverse its path-position role. All native parity/electrical mod8 guards are still checked separately. A decoder is not permission to discard those guards.

## 2. A sharp all-size bounded-core theorem

**Theorem.** In a tournament on n vertices, individualize r distinct actual vertices. If n−r≤5, the stable directed equitable quotient and canonical refinement record determine the marked tournament up to isomorphism. If n−r≤6, at most one additional shared-selector chirality bit is sufficient. Equivalently, use J mod8 from ENDPOINT_CHIRALITY.md when two mixed cyclic-three cells survive.

**Proof.** Every nontrivial stable cell has odd size, at least three. At most six unmarked vertices remain, so there is at most one regular five-cell, one triple, or two triples. The unique regular-five and C3 structures, uniform singleton attachments, and the one possible matching/co-matching block give exactly the order-eight proof. The number and arrangement of distinguished singleton vertices do not enter that argument. If only five unmarked vertices remain, two triples are impossible and no chirality bit is needed.

For Hamiltonian ports, the two endpoints already supply two actual marks. At order n≥8, marking any further n−8 actual path vertices guarantees the n−6 total-mark budget. At smaller orders the two endpoints suffice. These can be fixed path positions; retain their reflected positions when checking converse. This is a linear, high-information fallback, not a claim of uniformly small marking cost. An adaptive version stops earlier whenever the stable residual cells already satisfy the proved small-core condition.

For a complete literal source/path inventory, transport these observations by Section1, decode the target ordinary class, and aggregate into the required class-level metagraph edges. Preserve literal color, converse, loop and mod8 conditions. This is a constructive all-size sufficient class-wiring rule with an explicit marking budget. It does not separately distinguish literal target ports that are already isomorphic as marked tournaments, nor turn arbitrary feasible margin tables into complement-realizable wirings.

## 3. Sharp obstruction at seven unmarked vertices

The bound cannot be replaced by n−r≤7 with the SAME quotient-plus-triple-chirality observer. Let R and R' be the circulant tournaments on Z/7 with positive step sets {1,2,3} and {1,2,4}. Both are regular and self-converse; their nonisomorphism is certified by different Hamiltonian-path counts. Insert either as a uniform module between a specified source s and sink t, with s→R→t and s→t. Mark s,t.

Both order-nine tournaments have stable marked partition {s},R,{t}, quotient identical, and no triple block to which the chirality sidecar applies. The explicit path (s,0,1,2,3,4,5,6,t) is valid in both. Their marked tournaments are nonisomorphic because deleting the two marks leaves the nonisomorphic R,R'. They also remain distinct after converse merging, since each R is self-converse and the outer source/sink exchange identifies each whole tournament with its converse.

Add arbitrarily many individually marked vertices to a fixed transitive outer chain, uniformly before/after R, to obtain the same seven-unmarked-vertex obstruction at all larger orders. This is an observer counterexample, not a counterfeit pairing with one common fixed source class, and not a proof that an external global electrical mod8 guard also agrees.

## 4. Lower marking cost on a proved structural class

There is a different all-size extension when stable cells have size one or three. Their nonuniform regular blocks form a support graph on triple cells. Record each block's polarity and sign. If that support is a forest, these data and the stable quotient determine the whole marked tournament: cyclic chart rotations gauge away all translation phases along a tree. No cycle phases remain. This can use far fewer than n−6 marks, but its structural hypothesis must be checked.

For support containing cycles, retain the JOINT tuple of based affine-F3 cycle holonomies from Session39, modulo a single common root translation in each component. Signs alone do not suffice. This gives a complete decoder only for the specified equitable singleton/triple carrier. A regular seven-cell is outside its scope; a nontrivial cycle holonomy is a valid tournament, not a realizability failure. The three levels are distinct: a sufficient small-core budget, a potentially cheaper forest certificate, and the general shared-phase sidecar within the triple-cell carrier.

## Loss ledger and stopping boundary

Source: actual source tournament, actual Hamiltonian path, actual marked vertices. Target: transported target refinement record plus a shared rectangle sign or joint phase sidecar. Preserved: marked ordinary target isomorphism in the stated scopes, hence aggregated class wiring from a complete inventory. Lost by a quotient alone: internal regular-cell structure and phase/chirality alignment. Restored by the small-core bit: the only possible mixed C3-to-C3 orientation. Not restored: a regular seven-cell or general larger-cell structure. Cheapest obstruction: the explicit order-nine regular-module pair above.

## 5. Coprime-cell decoupling and a cheaper adaptive stop

**Lemma.** If an equitable tournament partition has cells X,Y of coprime sizes a,b, their adjacency block is uniform. Each row has q ones and each column p ones, so a q=b p. Coprimality forces q=0 or b (and p=0 or a). Equitability is essential; arbitrary partitions do not obey this conclusion.

Consequently every five-cell couples uniformly to every singleton and triple. The triple-forest decoder in Section4 extends to partitions with cell sizes1,3,5 provided all five-to-five blocks are uniform: each internal regular-five tournament is unique; all its external attachments are then quotient-determined. The same extension works with the joint triple-holonomy sidecar when the triple support has cycles. This gives an adaptive sufficient stopping certificate that can leave arbitrarily many vertices unmarked.

The five-to-five hypothesis cannot be silently removed. Even a degree-one block between two regular five-cells is a permutation in S5, not generally affine overF5. Their internal automorphism groups are cyclic C5, and the permutation block has eight C5-left/right double orbits: four of size5 and four of size25. A single permutation-parity bit leaves four orbits of either parity. Degree-two and degree-three blocks introduce still more data. The checker enumerates all120 permutation blocks and their25 chart rotations to certify this precise boundary; no general five-cell synchronization theorem is claimed.
