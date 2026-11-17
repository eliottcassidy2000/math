# Prime-producing three-cover menu with a physical 8-to-1 gate

Let B,C be arbitrary nonempty tournaments, X=TT3⊕B, and form F=C3[X,C,a], with arcs X→C→a→X and singleton a. Mark the first and third vertices of the TT3 prefix of X. Write g=H(F); E counts full paths ending in the two tags, W counts ordered two-path covers with second path starting tagged and first optional, and U counts those covers whose first path also ends tagged. The exact seed criterion for the strong ×5/7/9 menu is D:=U+E−g=0. Its source count when D=0 is J=3W+2g: append y beating only the tags in F, so h=H(F+y)=W under normalization, then apply the two-ear and sink construction with p=3h+2g.

Put v(T)=(H(T),P₂(T),P₃(T)), where P_j counts j labeled path slots, with empty slots allowed. Exact formulas are
    D = v(B) [[8,−1,0],[0,0,0],[0,0,0]] v(C)^T
      = H(B) (8H(C)−P₂(C)),
    3W+2g = v(B) [[9,−5,0],[−46,91,−30],[0,−78,78]] v(C)^T.

**All-order proof.** In a full path, removing a if present leaves one or two monotone runs through the ordinal order TT3→B→C, so any one of B,C occupies at most two nonempty consecutive runs. In an ordered two-cover only one path can contain a; hence either slot occupies at most three nonempty runs in total across the paths. All designated endpoint tags are outside B,C. Contract maximal runs to slot symbols in paths or covers. A single such skeleton contributes the product of the nonempty j-slot counts of B and C for its run multiplicities. These counts are H, P₂−2H, or P₃−3P₂+3H. Thus both quantities really are bilinear forms in v(B),v(C), regardless of size and internal orientation. Determine their nine coefficients by substituting TT_i,TT_j, i,j=1,2,3, where v(TT_i)=(1,2^i,3^i), three independent rows. Direct subset DP gives evaluation tables
    D = [[6,4,0],[6,4,0],[6,4,0]],
    3W+2g = [[325,1255,4159],[885,3691,12795],[2473,10903,38959]].
Multiplication by the proposed coefficient matrices reproduces these nine values uniquely. The replay also uses C3 and nontransitive C below as separate controls.

**Physical gates.** The exact normalization condition is P₂(C)=8H(C). The injection cutting each Hamilton path of C at its |C|+1 positions implies |C|≤7. The attached independent catalog fixes a Hamilton path using Rédei's existence theorem and enumerates each unused arc through order seven. It finds only the following four profiles (frame presentations, not isomorphism classes): (1,8,27), (3,24,99), (13,104,471), (15,120,525), with counts 1,2,13,3. The all-order bound is what makes this exhaustive. At order two and above P₂ is a multiple of four, by pairing complementary subsets; H is odd. Thus 9 is not an alternative integer ratio on any tournament, whereas 8 is the exact ratio here. For singleton C the ratio is 2.

Substituting the four profiles shows that only the first and third can give prime p. The second and fourth have all read coefficients divisible by 3 and yield values larger than 3. This cuts the normalized reader exactly to two possible prime-producing profiles before asking for the variable B. It is the physical obstruction on the receiver, not a direct application of Catalan's theorem.

For C=TT3, the result is the genuine three-coordinate p=1482P₃(B)−128P₂(B)−31H(B). For B=TT_k this means p=1482·3^k−128·2^k−31; in particular TT2 gives the third-slot 9 next to the gate's two-slot 8. No number-theoretic theorem about powers is being imported. With B=TT1 and TT3 the prime p's are respectively 4159 and 38959. A nontransitive seven-vertex B with masks (118,88,26,97,104,6,36) has (113,1444,9531) and gives prime 13936607. Changing the gate core to C=(22,12,24,17,2), profile (13,104,471), gives p=28626P₃(B)−5264P₂(B)−403H(B) and on the same B gives the prime 265187651. Trial division is independent of both structural formulas.

The third coordinate cannot be dropped even as a nonlinear read of the first two. For C=TT3, B masks (58,356,81,198,490,332,129,295,77) and (178,428,97,309,390,400,315,76,133) have the same (H,P₂)=(2207,35824) but P₃=300621 vs 300861. They give distinct exact normalized source counts 440866433 and 441222113.

Every F in this construction is strong by the directed three-block cycle; y has outgoing tags and incoming untagged vertices and therefore preserves strength. The marked-ear menu then gives strong witnesses for 5p,7p,9p. This is a genuine larger physical profile bank, not a selector for arbitrary prescribed primes. Replays: `python replay_receiver_gate.py`, then `g++ -O3 -std=c++20 -Wall -Wextra -Werror ratio8_catalog.cpp -o ratio8 && ./ratio8`.
