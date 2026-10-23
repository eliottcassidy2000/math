# The order-eight Z/256 factor belongs to the bridge-free core

23 October 2026. **PROVED:** the weighted support-bridge decomposition at
all sizes. **FINITE-EXACT:** its full order-eight localization and the
transport of every saved dyadic certificate into the core. No all-size
claim that metagraph support bridges are leaves is made.

## 1. Exact bridge splitting over the integers

Let G be a connected finite multigraph with positive integer conductances.
Loops do not enter its electrical Laplacian. Suppose {a,b} is a bridge of
the underlying simple support, carrying total conductance w. Removing this
whole bundle gives connected sides G_A and G_B, with a in G_A and b in G_B.
Parallel edges in that bundle are not individual literal bridges when w>1.

Ground the potential at b. Use new integer coordinates

    z = p_a,
    q_i = p_i-p_a for i in G_A except a,
    r_j = p_j for j in G_B except b.

This is an invertible integer change of coordinates: p_a=z, p_i=z+q_i
on the first side, and p_j=r_j on the second. Both maps have integer
coefficients, so their matrices are unimodular. The energy is exactly

    p^T A_G p = w z^2 + q^T A_(G_A,a) q + r^T A_(G_B,b) r.

There are no other cross-side edges. Equality of these quadratic forms
over the integers gives the unimodular congruence

    A_G ~ [w] direct-sum A_(G_A,a) direct-sum A_(G_B,b).  (1)

Therefore, as finite abelian groups and tree counts,

    K(G) = Z/w direct-sum K(G_A) direct-sum K(G_B),
    kappa(G) = w * kappa(G_A) * kappa(G_B).              (2)

A one-vertex side has the zero-dimensional reduced matrix and trivial
critical group. This proves the leaf case and the genuinely nonleaf case.
It also works componentwise without an all-order connectedness assumption.
Repeated splitting over every support bridge decomposes the group into
its support-bridge-free blocks and the cyclic bridge-weight factors.

## 2. The actual order-eight bridge inventory

The checker reconstructs undirected support directly from the full literal
edge ledger, summing both closed-to-open directions before a fresh low-link
bridge search. It finds exactly eight support bridges, each with one leaf
endpoint:

| Leaf node | Other endpoint | Total weight | Color |
|---:|---:|---:|---|
| 0 | 3376 | 1 | blue |
| 1 | 2388 | 2 | black |
| 3 | 2540 | 2 | black |
| 10 | 1856 | 2 | black |
| 40 | 1140 | 2 | black |
| 329 | 2383 | 2 | black |
| 331 | 2432 | 1 | blue |
| 3357 | 2387 | 1 | blue |

Node IDs refer to the complete saved `nodes.tsv.gz`, not to a standalone
canonical naming convention. After removing precisely these eight leaves,
the induced core C has 3,520 vertices, is connected, has **no support bridge**,
and has minimum support degree **three**. Thus there is no additional leaf
cascade and no hidden bridge bundle left in the core. Loops are retained in
the ledger but ignored when measuring support degree.

Equation (2) yields the exact integral decomposition

    K(G_8) = (Z/2)^5 direct-sum K(C),
    kappa(G_8) = 32 * kappa(C).                         (3)

Combining this with the certified full-network tower gives

    K(C)_(2) = (Z/2)^3353 + (Z/4)^4 + Z/256,
    v_2(kappa(C)) = 3369.                              (4)

The large factor is therefore not an effect of pendant doubled black
lines. Those account for exactly five factors of order two, and nothing
larger. Equation (3) also implies that the full network and core have the
same odd-primary critical group, without computing that group.

## 3. Independent core projection and transported integer certificates

The checker does not rely only on subtracting five from a published count.
It rebuilds the core Laplacian and computes its first half-form afresh:

    blue components = 3348,
    blue bicycle dimension = 11,
    dim V_1(C) = 3358,
    rank q_1(C) = 3353,
    dim V_2(C) = 5.

It then restricts the five saved full-network mod-four lifts to the core.
Their binary reductions have rank five. Their integer core Laplacian images
are all divisible by four, and their complete q_2 matrix is unchanged,
with rank four. The surviving potential is likewise restricted at every
modulus 8,16,32,64,128,256. Direct integer core-edge sums give

    q_3 = q_4 = q_5 = q_6 = q_7 = [0],   q_8 = [1].

All those restrictions retain the same nonzero binary class and the same
grounded core root. This independently transports the exact Z/256
certificate into the bridge-free core.

Why restriction works is also elementary. For a leaf l attached to v with
weight w in {1,2}, full divisibility A_G u = 0 modulo 2^k implies
w(u_l-u_v)=0 modulo 2^k at l. Removing the leaf changes the v equation by
that same divisible term, so the core restriction still qualifies.
For k>=2, the omitted bridge-energy term satisfies

    w(u_l-u_v)^2 / 2^k = 0 modulo 2.

Indeed u_l-u_v is divisible by 2^k for w=1 and by 2^(k-1) for w=2.
Thus every higher half-energy pairing is preserved; the cross-pairings
follow by the same product-divisibility argument. In particular a full
binary lift supported only on those leaves must vanish at k>=2, so no
nonzero higher dyadic mode can disappear entirely into the deleted leaves.

## 4. Reproduction and scope

From this directory, with the complete ledger in `../wiring`:

    python check_bridge_localization.py --base ../wiring
    python -O check_bridge_localization.py --base ../wiring --out bridge_localization_receipt.optimized.json

The checker also directly verifies (1) for six weighted nonleaf controls
with three vertices on each side and bridge weights 1,2,3,4,5,16. All
checks use explicit exceptions and survive optimized Python. Normal and
optimized receipts agree exactly. Dependencies are Python's standard
library and the already frozen `check_electrical.py` helper.

This addendum localizes, rather than explains a uniform combinatorial
origin for, the order-eight high-torsion factor. It proves no all-order
bridge count, no all-order exponent bound, and no assertion that every
future metagraph bridge has a leaf endpoint.
