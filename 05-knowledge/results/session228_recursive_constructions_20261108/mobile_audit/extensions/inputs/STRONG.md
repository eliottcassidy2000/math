# Strong physical fibers with fixed deletion count three

Status: PROVED; independent audit requested. This is an addendum to PROOF.md and does not modify its frozen theorem or receipt.

For every r≥3 and p=1,...,2r−4, start with T((01)^r,p) from PROOF.md. Delete the initial singleton0 and final singleton1. Call the rooted result V_(r,p). Then V_(r,p) is a strongly connected tournament on2r vertices, and

    x(V_(r,p),a)=(B_r+(-2)^p−2, 3, 2^(2r−1)+1, 3),
    B_r=(4^r+6r−1)/3.

In particular one physical sidecar fiber contains2r−4 distinct scalar values, all witnessed by strong tournaments. Strongness therefore does not restore a bounded-fiber theorem.

Proof: the initial0 is unreplaced because p≥1, and is a global source. The final1 is unreplaced because p≤2r−4, and is a global sink. Removing either does not change D, since its unique position in every Hamilton path is forced. Neither operation changes d, since it removes a global source or sink in the root-deletion graph. Removing the initial source does not change U by the bijection in MEMBERSHIP.md Section2. Removing a global sink divides U by2: the sink can be appended either to P (including its empty case) or R, and these two possibilities are bijective with covers of the smaller graph. Thus U=(4^r+2)/2. After deletion, the last block is the original position2r−2, an unreplaced0 because p+1≤2r−3. Its vertex beats the root; every root-deletion Hamilton path ends there, so E=d=3.

For strongness, the root-deletion graph remains an ordered sum of strongly connected blocks (singletons and one cyclic triangle). The root beats a vertex in the first block: the first block is either the original singleton1 at position1, or, for p=1, a triangle replacing10 and having two outneighbors from the root. The root loses to the final singleton0. It can therefore reach every vertex: enter the first block, traverse its strong connectivity when needed, and move forward. Conversely every vertex can move forward to the last block and then to the root. The paths through the root show that every vertex reaches every other. This reasoning also covers a first triangle and makes no assumption that all root-to-triangle edges agree.

The scalar values are distinct because the powers(-2)^p are distinct. The boundary r=3 yields two strong six-vertex states (23,3,33,3) and(29,3,33,3). No prescribed interval or all-odd realization statement follows.
