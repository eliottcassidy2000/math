# A new target reader genuinely using a third slot

Build a strong tournament F with a forward spine of blocks

    1 → 1 → 1 → B → 1 → 1 → 1,

where B is any nonempty tournament and the first/last groups of three are singleton vertices. An additional apex beats the first three and all of B, and loses to the last three. Mark the first and third spine vertices. Let g=H(F), W be the weighted start count, U the joint two-cover count, and E the full marked end count. Then

    g = 56P₂(B)−5H(B),
    W = 494P₃(B)−80P₂(B)−7H(B),
    U+E = g.

Consequently the normalized same-source theorem gives the ordinary count and its strong sevenfold lift at

    p = 3W+2g = 1482P₃(B)−128P₂(B)−31H(B).

All arcs of the growing B are retained. The graph F is strong because the apex reaches the early blocks and the late blocks return to it. The rooted seed M for the normalized construction is obtained by adjoining y beating just the two marked vertices of F. Then M is strong (there are both tagged and untagged vertices in F), C in the menu proof is strong by the return through x,z,y, and the recontact with tags y,s restores strength after the singleton s and the two return vertices are added.

**Proof of the coefficient identities.** A full path in F can cross B twice only by using the apex between the two visits: the spine quotient alone is ordinal. Likewise each of two paths in a cover can meet B once without the apex, but only one of the two can use the apex to meet it again. Hence at most three B-runs occur. Every read g,W,U,E is a fixed linear combination of ordered nonempty one-, two- and three-run cover counts of B. This is equivalent to a linear read of H,P₂,P₃ by binomial inversion. Transitive B on 1,2,3 vertices has profile (1,2^k,3^k); those three vectors are linearly independent (determinant 12). Literal subset DP at those anchors gives:

| k | g | W | U+E | p |
|---:|---:|---:|---:|---:|
| 1 | 107 | 1315 | 107 | 4159 |
| 2 | 219 | 4119 | 219 | 12795 |
| 3 | 443 | 12691 | 443 | 38959 |

These exact data prove both the normalization and the displayed laws for all B. A nontransitive triangle B has (H,P₂,P₃)=(3,12,33), giving the independent test p=47277. The replay verifies the profiles and certifies 4159 and 38959 prime by trial division. It does not show that the set of readers or possible profiles hits every prime. For readers such as this, higher cover costs must remain visible; fitting the lower two coordinates alone misses a genuine independent slot.

## What the third slot does and does not buy

The new read is at least 1323P₃(B), using H≤P₂≤P₃. Hence it also satisfies the previous branching obstruction whenever its inputs belong to that grammar. For one-apex B specifically, the obstruction strengthens: the same extreme inversion in an active spine word of length m contributes 3(4^(m−2)+2^(m−2)) to P₃. In the odd-cycle expansion, all cycles are through the apex; an inversion at distance d can include any even subset of the d−1 intervening spine vertices. With slot weight 6 for each odd cycle, each intervening vertex selected into the cycle costs one singleton weight of 3. The even binomial sum gives 3^(L−d)(4^(d−1)+2^(d−1)) for that inversion in a spine of length L. As there are only 2^(m−2) active words, one-apex input to this new reader yields O(X^(1/2)) possible p≤X after summing singleton padding geometrically. Arbitrary ordinal sums of these inputs yield O_epsilon(X^(1/2+epsilon)) via multiplicativity of P₃ and the same multiset Euler product. The new coordinate improves expressive power of the interface, but exact prime placement still requires a richer source pool.

For an explicit nontransitive prime example, B with outgoing rows (118,88,26,97,104,6,36) has profile (113,1444,9531), giving prime p=13936607 and its sevenfold companion 97556249. The same odd-cycle expansion gives P₂(B)≡(−1)^|B| H(B) (mod 3): each odd cycle changes the parity of singleton count once. Thus this third-slot reader is automatically 0 mod 3 on even-order inputs or inputs with count divisible by 3, since p≡P₂−H (mod 3). The remaining odd-order, 3-free inputs avoid this particular local obstruction, not every obstruction.
