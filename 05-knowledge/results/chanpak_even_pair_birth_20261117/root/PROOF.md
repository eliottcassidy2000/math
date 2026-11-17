# What an even-order pair remembers, and what it loses

The user's even-only idea has a precise boundary. Even hosts can recover odd first orders if they retain marked structure. Unmarked even scalar counts are not enough. Pairs do furnish a simple exact growth at every size, but it acts on the entire slot polynomial, not just the Hamiltonian total.

## Even data and the marked interaction

Write P_q(T) for the number of ordered covers by q paths, with empty paths allowed, and P_1=H. The [Grinberg–Stanley odd-cycle expansion](https://arxiv.org/abs/2307.05569) views P_q as a polynomial in q: a spanning packing has a factor q for each singleton or nontrivial odd directed cycle, and an extra factor 2 for every nontrivial cycle.

In an even-order tournament E and for marked vertices a,b, set

    Γ_q(a,b)=P_q(E)−qP_q(E−a)−qP_q(E−b)+q²P_q(E−{a,b}).

The subtracted terms count packs where a or b is a singleton. Inclusion–exclusion shows that Γ_q has nonnegative coefficients and counts packs where *both* are on nontrivial cycles. It includes either one common odd cycle or two disjoint odd cycles, not merely adjacency of the pair. The same inclusion–exclusion for any marked subset counts packs in which every marked vertex is on a nontrivial cycle.

At q=1 the genuine even-to-odd midpoint is therefore

    H(E−a)+H(E−b)=H(E)+H(E−{a,b})−Γ_1(a,b).

This is not arithmetic interpolation without a term: Γ_1 is an actual marked interaction. If all these marked pair data are retained, three distinct vertices suffice to recover an individual odd deletion, by adding the two pair sums involving a and subtracting the third, then dividing by two. The full polynomial version holds in the same way after division by q. This gives a pair representation of the midpoint, not a proof that the needed interactions are simple.

A literal hostile is

    A=(6,12,48,5,11,27),   B=(2,12,25,49,35,7).

Both have H=33. For **every even subset size** their multisets of induced Hamiltonian counts coincide: at size four there are three 1s, three 3s and nine 5s. Nevertheless the five-vertex deletion multisets are (3,11,11,11,15,15) and (5,9,9,11,13,15), respectively. So even the entire even scalar deck, not just the two adjacent scalar levels, does not identify the intervening odd counts. This says nothing against a deck retaining which pairs were deleted, or their actual quotient graphs.

## A single physical operator adds a pair

For a nonempty T, attach two vertices x,y with x→T→y→x. Call this strong envelope W(T). It adds two vertices, so stays within either parity. A direct odd-pack computation gives the full transfer

    P_q(W(T)) = L P(q),       L f(q)=q[f(q+1)−f(q−1)+qf(q)].

Indeed a mixed odd cycle must contain x,y and one odd-length-in-vertices path in T, with the complement weighted by P_q. For the unrestricted optional chosen path the sum is P_(q+1)(T). For the sum signed by its order, assign each odd cycle or singleton in a packing to either one negative color or one of q positive colors, giving P_(q−1)(T). Their difference picks out the odd paths twice. If x,y are both singletons, the contribution is q²P_q(T). This proves the formula at integral q≥1, hence polynomially.

The first two reads are

    H(W(T))=P₂(T)+H(T),
    P₂(W(T))=2P₃(T)+4P₂(T)−2H(T).

So a pair step does not close on H or even (H,P₂). More exactly the shell sequence H(W^k(T)) determines P_(k+1)(T) from its first k+1 terms: the coefficient of P_(k+1) in the k-th read is k!, from the one all-upward path through the q shifts, and the remaining terms only use earlier positive slots (P₀=0). This gives a triangular exact decoding. At known order n, P_q has degree n, parity n, and zero constant term. Consequently ⌈n/2⌉ initial shell readings reconstruct the entire slot polynomial by interpolation. No bounded-length profile is claimed sufficient as n grows.

Padding by one neutral singleton is multiplication by q on P. It does **not** commute with taking a pair step:

    [L,q] f(q) = q[f(q+1)+f(q−1)].

Thus pad-then-pair and pair-then-pad differ exactly by P₂(T) in H despite having the same final order. This is a specific obstruction to a parity-blind blending rule.

## Finite modular recurrence and its prime obstruction

Modulo an odd prime p, write d=(p−1)/2 and ε=(−1)^|T|. Since P_q is a polynomial of parity n, P_(p−q)≡ε P_q and P₀≡0. The L update on its first d evaluations is tridiagonal: diagonal q², upper q, lower −q, with last diagonal d²+εd. Therefore **every** pair-shell sequence modulo p satisfies the same order-d characteristic recurrence for its parity, regardless of the starting tournament. Compute the characteristic polynomial by tridiagonal continuants, starting C₀=1, C₁=t−a₁ and C_j=(t−a_j)C_(j−1)+j(j−1)C_(j−2), where a_j=j² except a_d=d²+εd. In particular at p=3,

    H(W(T)) ≡ 2H(T) mod3 if |T| is even;
    H(W(T)) ≡ 0 mod3 if |T| is odd.

So repeating this very natural pair operation in the odd sheet cannot construct any prime count >3. The right use of pairs on hard targets is to introduce different contacts or change the pair interface, not only to extend the same envelope.

## The 3,4,5 example and the repeated pair

The [earlier F=U+S work](../../05-knowledge/results/divisor_layers_pair_geometry_20261004.md) separates labels, multiplicative letters and geometric pair data for p²qr. The analogous useful distinction here is between **vertex order**, **Hamilton count**, and **slot multiplier**. A triangle has (order,H)=(3,3), and the unique strong four-vertex block has (4,5). Their ordinal product has H=15 at order 7. But putting the triangle inside W adds just a pair: its P₂ is 12, and the new five-vertex tournament has H=12+3=15, genuinely at the first order for 15. Repeating this smaller block ordinally gives 15^k at order 5k, while the [largest singleton-free order](../chanpak_birth_cubes_smooth_20261117/PROOF.md) is 7k. The 5k construction is an upper bound on first order, not an exact formula for all k.

A neutral two-vertex pad contributes q² to P_q: it contributes 4 at q=2 but 1 to H. The product 5·3·2² is therefore a useful *mixed* reading of these roles, not a tournament Hamilton count. This is why the repeated pair from the divisor analogy is useful when its two channels are kept distinct.

Replay `python replay.py` for the literal equal even decks, conflicting odd midpoint, full pair transfer at q=1,…,4, commutator, and order-five compression.
