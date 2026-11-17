# What a homogeneous tournament seam can hear

A nonempty tournament B is substituted at one vertex of an exterior tournament Q. An exterior vertex has the same direction to every vertex of B. A *run* of B is a maximal consecutive nonempty segment in a Hamiltonian path of the resulting tournament. Let C_r(B) count ordered spanning covers by r nonempty directed paths, and let P_q(B) count ordered covers in q slots allowing empty slots. Thus

    P_q = sum_(r≥1) binom(q,r) C_r,        P_1 = H.

**The natural messages for Hamiltonian subscribers are J_r(B)=(r−1)! C_r(B), not C_r or P_r alone.** There are explicit physical subscribers which detect exactly this tuple over the integers or modulo any M.

## Factorial quotient and sharp readers

Compress every B-run of a Hamiltonian path to one slot symbol. The remaining symbols are distinct exterior vertices, and successive symbols must follow Q-arcs (in particular two slot symbols cannot be adjacent). Let a_r(Q) count these skeletons with r slot occurrences. Any skeleton with r occurrences has exactly C_r(B) lifts.

Between consecutive slot occurrences there are r−1 nonempty exterior excursions. Each excursion starts at an out-neighbor of the slot and ends at an in-neighbor. They occupy disjoint sets of labelled vertices. Permute just these r−1 excursion words. This never changes legality, and is a free symmetric-group action. Therefore

    H(Q[B]) = sum_r b_r(Q) J_r(B),        b_r(Q)=a_r(Q)/(r−1)! ∈ Z≥0.

The exact factorial cannot be improved for all Q. For each r choose 2r−2 exterior vertices in transitive order a₁,b₁,…,a_(r−1),b_(r−1), with the slot beating all a_i and losing to all b_i. There are at most r slot runs. At the top degree each excursion must be exactly a_i,b_j, covering all exterior vertices; the forward edges force i≤j, hence the matching is the identity. Permute the pairs arbitrarily. Thus a_r=(r−1)!, giving b_r=1. For r=1 the context is just the slot.

The resulting coefficient matrix is unit lower triangular over Z. Consequently two nonempty cores are indistinguishable to **every** Hamiltonian uniform subscriber with at most R runs if and only if their first R values J_i agree. This is true over Z or modulo **any** M; reversal/converse is just one way to make such equality. In particular independence of low-slot observations is not assumed: these are sufficient statistics *on physically realized cores*, not surjectivity of the statistics.

## Modular cutoff and automatic transmission

For M≥2 a convenient cutoff is the least r≥1 for which M divides r!(r−1)!. Indeed C_r itself is divisible by r! by freely permuting the ordered nonempty paths, so J_i is automatically zero modulo M at every i≥r. This is a safe cutoff, not a claim that every coordinate below it is nonzero on physical cores. The preceding calibrated readers give an exact finite test for indistinguishability even at unbounded depth and order. Uniform substitution is associative at one marked leaf: a chain of wrappers is itself a single uniform context for that leaf. Thus the equality passes through every depth and every choice of wrapper, without rebalancing at intermediate stages. This does not apply to a reader that subsequently makes selective vertex contacts inside the hidden core.

For an odd prime p one can shrink the packet further for two cores of the same order parity. The primary [Grinberg–Stanley odd-cycle identity](https://arxiv.org/html/2307.05569#S1.SS8.SSS3) at q variables equal to 1 writes P_q as an integer polynomial of the order's parity with P_0=0. Therefore P_(p−q) ≡ (−1)^|B| P_q and P_p≡0. For r<p, factorial division from J_r to C_r is invertible modulo p; binomial inversion then recovers all needed C_r from P_1,…,P_(p−1). At r≥p the ordered-cover factorial makes J_r vanish. Thus at a fixed parity only P_1,…,P_((p−1)/2) modulo p need to be published; this test is necessary too. Parity cannot simply be discarded.

## An obstruction to fine higher-channel switching

Suppose B,B' have the same positive order and exactly equal P_1,…,P_(r−1). Their slot-polynomial difference has all coefficients even: in each odd-cycle expansion the only coefficient not guaranteed even is the common all-singleton leading term. Its parity and zero constant term give the polynomial divisor

    2 q ∏_(i=1)^(r−1) (q²−i²)

in Z[q] (the remaining factor may have an additional q). For integer q≥0, the product without the 2 is (2r−1)! binom(q+r−1,2r−1). Binomial inversion gives 2(2r−1)! | ΔC_j for every j; ΔC_j=0 for j<r. Every exterior skeleton at j≥r contributes an additional factor (r−1)! by its independently permutable excursions. Hence every uniform subscriber obeys

    H(Q[B'])−H(Q[B]) ≡ 0 mod 2 (r−1)! (2r−1)!.

For r=2 this gives 12 when H is exactly equal; for r=3 it gives 480 when H and P₂ are exactly equal. At even order the polynomial divisor contains one additional q, and the resulting modulus has an additional factor r: q binom(q+r−1,2r−1) = r[binom(q+r−1,2r−1)+2 binom(q+r−1,2r)]. Thus at r=3 the modulus is 1440 at even order. It is sharp: the order-eight cores (136,89,75,128,45,15,57,118) and (42,32,163,150,39,200,31,83) both have H=345, P₂=5640, and their P₃ values are 44505 and 45225. The calibrated three-run reader has top coefficient 2 on C₃, producing the difference 2·720=1440. They are indistinguishable to all one- and two-run subscribers, but not to the three-run subscriber. These bounds concern differences of two outputs of the *same* exterior reader and same-order input blocks. They are not obstructions on the unrestricted path spectrum, and not bounds on readers which change at the same time. They explain why retaining more exact topics while trying to fill small gaps can make things worse.

## Physical producers that stay invisible modulo any chosen M

Let feedback block R_m be the transitive (m+2)-vertex tournament with its extreme arc reversed. Each nontrivial odd cycle consists of the extreme arc and an odd nonempty subset of the m intermediate vertices, in increasing order; two such cycles cannot be disjoint. Therefore

    P_q(R_m)=q^(m+2)+q[(q+1)^m−(q−1)^m].

Take e at least the largest prime-power exponent occurring in M (e=M.bit_length() is more than enough), choose m=2φ(M)e, t=φ(M)e. For p | M and p ∤ q, P_q(R_m) is a p-adic unit. Indeed reduce modulo p: if q±1 are both units, both powers are 1; if exactly one vanishes, q≡±1 and the result is 2; if both vanish (p=2), the result is q²=1. Euler's congruence now shows P_q(R_m)^t ≡ q^((m+2)t) modulo the relevant prime power. If p | q, each P_q contains a q factor, and raising to t≥e kills both sides modulo that prime power. Thus the ordinal sum of t copies of R_m and the same-size transitive tournament publish identical slot values modulo M for every q≥0. These are genuine different Hamiltonian counts: (2^m+1)^t versus 1. Any positive multiple of t works as well. If strong cores are needed, substitute the pair into any common strong outer quotient, e.g. the directed triangle. The new pair remains indistinguishable to every further homogeneous context modulo M.

This is a reusable way to keep a low residue without engineering each downstream seam. It does **not** meet the absolute placement scale: the construction is intentionally coarse. For fine correction the negative result above says which assumption to release: keep congruence of the subscribed topics rather than exact equality, or let the next graph contact selected vertices unequally.

Run `python replay.py`: factorial rows through 8, 80 arbitrary literal substitutions, finite congruence and feedback-block controls. The all-order statements follow from the arguments above, not the replay bounds.
