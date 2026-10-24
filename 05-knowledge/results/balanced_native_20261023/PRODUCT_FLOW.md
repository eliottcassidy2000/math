# Sharp refresh congestion with retained kernel coordinates

**PROVED.** Let X=product_(i=1)^d A_i, with finite alphabet sizes q_i>=2 and the uniform law. At each integer event choose coordinate i with probability w_i>0, sum w_i=1, and resample that coordinate uniformly. Assume every such one-coordinate change is an authenticated primitive; no missing edges are silently filled. Put a_i=1-1/q_i and A=sum a_i.

For independent uniform endpoint refresh demand, the optimal length-weighted path congestion is exactly

    Lambda = max_i (1+A-a_i)/w_i.                       (1)

This is an all-route lower bound, including fractional/randomized routes, and is attained by changing differing coordinates in any fixed order. The native event length counts actual primitive transitions. Other work costs must remain separate ledgers.

## Proof

Any endpoint pair differing in coordinate i must use at least one i-edge, and every route has length at least D, the number of differing coordinates. Sum length-weighted loads over the directed edges of coordinate i. Their total load is at least

    E[D 1_(x_i != y_i)] = a_i (1+sum_(j!=i) a_j).

Their total directed conductance is exactly w_i a_i. Taking the ratio proves the i-th lower bound in (1). This uses a coordinate cut, strengthening a bound that sums all edges and forgets which coordinate carries the load.

For the fixed-order shortest routing, every directed i-edge has unweighted endpoint flow1/(|X| q_i), and conditional expected route length1+sum_(j!=i)a_j. Its native conductance is w_i/(|X| q_i). Thus its load ratio equals the i-th expression in (1), proving optimality.

The best choice of coordinate probabilities for this demand is

    w_i=(1+A-a_i)/[d+(d-1)A],
    Lambda_min=d+(d-1)A.                               (2)

Indeed max c_i/w_i>=sum c_i for sum w_i=1, with equality at proportional weights. This optimization targets this congestion, not every mixing objective. The exact product-chain spectral gap is min_i w_i (its coordinate resampling projections commute); the weights maximizing that gap are equal. These objectives need not choose identical weights.

## Actual native kernel application

The concurrent `native_packet_decomposition_20261023/workspace/ONE_BIT_WORKSPACE.md` authenticates S_r=T_r -> C3. Its full records are (b,k) in F2^r times Z/3. The three k-values are retained owner lifts over ONE normalized resource b. Cube moves change two owner clocks; a nontrivial kernel rotation changes three. Appending one genuine rigid block supplies a resource bit without turning the three lifts into three capacities.

For a single full-record chain on S_r, all three kernel transitions are native (C3 has both cyclic directions), so (1) applies with r binary coordinates and one ternary coordinate. For equal coordinate probabilities1/(r+1), r>=1,

    gap=1/(r+1),
    Lambda=(r+1)(3r+7)/6.

The bit-coordinate cuts determine the maximum; the kernel-coordinate load ratio is (r+1)(r+2)/2. The optimal probabilities for refresh congestion instead give

    Lambda_min=(3r^2+10r+6)/6,

an improvement of1/6 over equal coordinate selection. This small difference is exact. It is not evidence of a new resource: every path and its kernel coordinate remain in the record, while multi-slot capacity is still tested after mapping (b,k) to b.

This result concerns one full-record chain and its uniform marked endpoint law. It is not a product law for several simultaneously occupied slots, a balanced interface for arbitrary nonuniform weights, or full-atlas coverage. Actual nonuniform source laws need their own conductance and demand comparison. In particular a quotient that forgets k cannot charge kernel rotations correctly merely because it preserves resource counts.

`product_flow.py` checks (1) and (2) by exact directed path loads on small mixed alphabets, including the native family through r=5, with unequal q and rational probabilities. The general proof is the coordinate-cut argument above.
