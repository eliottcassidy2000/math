# A uniform counting boundary for finite-alphabet incidence signatures

Status: PROVED elementary counting obstruction. It concerns class separation, not uniqueness of an entire metagraph from an authenticated inventory.

Suppose each outgoing arc receives one of L fixed symbols, each vertex retains only the multiset of its outgoing symbols, and the tournament signature is the multiset of these vertex rows. A row is an L-component nonnegative integer count vector with each coordinate at most n. Hence there are at most(n+1)^L possible rows, and at most(n+1)^(Ln) signatures on n vertices, even before taking the outer multiset quotient.

For n>=3, orient one fixed Hamiltonian cycle and choose all other arcs freely. These2^(n(n−1)/2−n) tournaments are strong. An isomorphism class has at most n! labelled members. Therefore strong isotype counts are at least2^(n(n−1)/2−n)/n!.

If

(n−3)/2 > (L+1)log_2(n+1),

then this lower bound exceeds(n+1)^(Ln), using n!<=(n+1)^n. Two nonisomorphic strong tournaments must share the signature. This proof is independent of the particular symbol-assignment rule.

For the12-symbol E_(4,3) menu, the inequality holds at every n>=256. Indeed choose2^k<=n<2^(k+1), k>=8. The left side is at least2^(k−1)−3/2 and the right side at most13(k+1); the strict inequality holds at k8 and its gap increases thereafter. This is an explicit sufficient threshold, not the first collision order.

The separate15-vertex E collision is much sharper for the actual incidence rule and even preserves H-distinguishing information loss under explicit extensions. The counting result remains useful because it rules out every fixed L-symbol row-multiset scheme at all sufficiently large sizes, not only that one construction.

It does not prove different H in the counting collision, does not rule out an alphabet growing with n, and does not imply that capacity constraints cannot resolve signature-colliding ports. A source inventory contains more than a single node signature. Any stronger all-size claim requires its own proof.
