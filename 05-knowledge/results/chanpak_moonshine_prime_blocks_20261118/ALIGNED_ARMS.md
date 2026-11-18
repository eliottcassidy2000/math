# A continuation-preserving positional inverse without the common sevenfold state

The [two-axis design](PAIRED_AXIOMS.md) need not transport the same state on both sides. Here an ordinary target chooses where a distinguished rooted arm goes in a word, with the continuation exactly unchanged by its position. An independently chosen upper graph supplies 7p in the examples.

## Shared root, different ordinal order

A rooted arm consists of an ordinary tournament plus a distinguished root. For its deletion D put g=H(D), b=f_D(2). For each nonempty S⊆D give each directed Hamiltonian path of D[S] weight H(D−S). Define U by summing the weights of paths whose first vertex is beaten by the root, and V by signing the same weights by (−1)^|S|. Define Ubar,Vbar by the corresponding last vertex beating the root.

Glue two arms by identifying the roots and directing every deletion vertex of the earlier arm X toward every deletion vertex of the later Y. The exact rules are
```
g' = g_X g_Y, b' = b_X b_Y,
U' = b_Y U_X + g_X U_Y, V' = g_X V_Y,
Ubar' = b_X Ubar_Y + g_Y Ubar_X, Vbar' = g_Y Vbar_X,
H' = g_Y H_X + g_X H_Y − g_X g_Y + U_X Ubar_Y + V_X Vbar_Y.
```
For U, a nonempty marked start can first enter X then visit an optional Y path, or first enter Y. The signed sum for an optional Y path is f_Y(0)=0 for a nonempty deletion, so the first term drops in V; ends reverse the roles. For H use odd-cycle packing: a cycle through the root can be confined to one arm, or run through a nonempty X path then a nonempty Y path. Its exterior has even size; the two parity reads select it, and the half in the selector cancels the weight two. The other cycles lie entirely in the deletions.

## A physical companion for every normalized arm

Assume V=0 and 0<rho=U/(b−g)<1. These are exact constraints, not statements about arbitrary arms. Write rho=a/d in lowest terms. The numerator is even because V=0 balances odd-size and even-size support, and the denominator is odd since b is even and g is odd. Choose k with 2^k≡1 modulo d. For a transitive deletion of size k, the binary number encoded by the root's outcontacts (first vertex highest bit) equals U. Its signed read equals minus the last bit. Therefore choose the contact integer a(2^k−1)/d; it is even and less than 2^k. The resulting actual arm has g=1, b=2^k, V=0 and exactly the same rho.

Any word of aligned arms (same rho, V=0) has U=rho(b−g), by the gluing law. Thus their order can change H while **leaving (g,b,U,V) invariant**, the reads needed for any future right attachment. In particular let A be the core, B a transitive-deletion companion, r=b_B, and delta=U_A·Ubar_B−U_B·Ubar_A. In the word B^j A B^(m−j), moving the core one position to the right subtracts exactly r^j delta from H. Locally swapping A,B changes the end read by delta/rho; the unchanged incoming prefix amplifies the local difference by g+U/rho=b=r^j. All suffix B arms have g=1 and the start reads agree. Therefore
```
H_j = H_0 − delta (r^j−1)/(r−1).
```
For chosen core and length, a prescribed absolute target now directly tests one power and chooses the physical position. This is a recursive positional decoder even though the two arms do not have proportional states.

## Two independent lifts

Take A=(26,12,9,16,6) rooted at 4 and B=(6,12,8,1) rooted at 3. Their (g,b,U,V,Ubar) are (3,24,12,0,15) and (1,8,4,0,6); rho=4/7 and delta=12. The two orders already give 97 and 85, with identical continuation. For m companions the j=0 height is (576·8^m+84m+61)/49, from iterating the gluing law. At m=4 the prescribed primes 48,157 and 48,049 select positions 0 and 2 respectively.

On the upper axis use a different ordinary family, not a putative shared-root replica: one apex over a transitive spine of length 20 with low-bit-first outmask 402849 gives 337,099=7·48,157; length 19 with outmask 15509 gives 336,343=7·48,049. Every directed cycle uses the apex. An inversion in its contact bits at i<j contributes 2^max(1,j−i−1) to H−1 by choosing an even number of middle vertices. This verifies the exact upper counts without a large path census.

This bypasses universal proportional states in two actual prime pairs. The general remaining step is to find enough cores and companion words (or further independent upper readers) to place each prescribed target. Exact position choice once the core and length are selected does not itself solve that supply problem.

Replay: `python3 aligned_arms.py`.
