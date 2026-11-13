# An unbounded clock with changing, nontransitive modules

**Proved, restricted chart.** The [Π/Γ clock](../target_collatz_transitions_20261113/collatz/CLOCK_BOUNDARY.md) admits a refill by an arbitrary actual module, not just a singleton. One can transfer a whole nontransitive ordinal factor across the same separator, preserve its internal arcs, and then expose another arbitrary module for the next round. An explicit choice gives infinitely many exact valuation-one seven multiples with no fixed common prime divisor apart from 7 on the indicated subfamily. This does **not** prove infinitely many prime cofactors or arbitrary prime/7p reachability.

## 1. Whole-factor crossing and re-entry

In the notation of the inherited clock, put M in slot u and N=N₀→U in slot v, where U is any nonempty tournament and N₀ can be empty; all arcs in the displayed ordinal join point from left to right. Reverse exactly the arcs from the vertices of U to the separator c. Because U is an entire final ordinal factor of the v-slot, the moved vertices now form the first ordinal factor of the u-slot. The new graph is

    Π[U→M,N₀]  when N₀ is nonempty,
    Γ[U→M]    when N₀ is empty.

This can be checked incidence by incidence in the nine-vertex skeleton: the neighbors of v and u differ only at c (apart from the two slots themselves, where N₀→U→M already holds); no internal arc of U, M, or N₀ changes. The old names a and c remain fixed.

Write H0=H(M)H(N), and Q_W(t)=P_W(t) for the ordered t-color cover polynomial. Multiplicativity under ordinal joins gives the interior transfer exactly:

    H(Π[U→M,N₀])−H0 = [Q_U(3)/Q_U(2)]·(H(Π[M,N])−H0).

If N=U, the boundary becomes

    H(Γ[U→M]) = [Q_U(3)/Q_U(2)]·(H(Π[M,U])−H0)−17H0.

In particular the inherited −18H0 deviation from the formal continuation still applies. The boundary does not use a fictitious empty module in Π.

Now add a disjoint copy of an arbitrary nonempty tournament R, letting a beat all its vertices and letting all its vertices beat the remaining old vertices. This is exactly insertion into the now-vacant v-slot, because the only in-neighbor of the v-site outside that site is a. It creates Π[U→M,R], with H exactly

    60 Q_M(3) Q_U(3) Q_R(2) + H(M) H(U) H(R).

Thus the same marked recharge site is ready again, but the carried module has genuinely changed and may grow arbitrarily. The formulas make no assertion that another chart's rooted statistics, unit-root condition or primitivity at another port are retained. The return is physical and indefinite, not neutral and not a directed path along ordinary Collatz values. For singleton U=R this specializes to the inherited native recharge; mixed U,R provide distinct growth factors, with the operation defined on literal tournaments.

## 2. A primitive two-coordinate seven-adic example

Let U be the tournament with outgoing bit rows (40,1,3,6,15,30), in the order 0,…,5 and bits numbered from 0. Its counts are (H(U),Q_U(2),Q_U(3))=(19,256,1611), by the finite subset path-color calculation supplied with the note. Here Q_U(t)=sum over assignments of vertices to t labeled classes of the product of the Hamilton path counts in those classes, with empty classes contributing 1; this is a direct exact counting definition.

Start with one copy of U in each slot and use U again for every recharge. The resulting readouts for the carried module consisting of n≥1 ordinal U factors are

    H_n=15360·1611^n + 19^(n+1).

The two summands are coprime, since 19 does not divide 15360 or 1611. This is only a two-coordinate gcd statement, not a claim about marked primitivity. Modulo 7 the ratio 1611/19 is 3 of order 6, and the starting sum at n=0 is 15379≡0 (mod 7). Thus 7|H_n iff n≡0 mod6. For n≡0 mod42, the two clock bases are units modulo 49 and both raised to powers divisible by 42, giving H_n≡15379≡42 mod49. In particular v₇(H_n)=1 for an infinite subsequence of physical outputs, rather than just divisibility by 7 at an unspecified valuation.

On that subsequence the gcd of all H_n/7 is 1. To check this, first note gcd(H_n, 2·3·5·19)=1 for all n≥1. If a remaining prime ℓ divides every H_(42t)/7, then ℓ divides H_(42t) for all t≥1. If ℓ divides either 1611 or 19, it cannot divide any H_n because the two summands are coprime. Otherwise the recurrence in t has invertible constant coefficient modulo ℓ, and vanishing at all t≥1 forces H_0≡0 mod ℓ. Since H_0=7·13³, only 7 and 13 need checking; exact v₇=1 removes 7, and checking t=1 removes 13: 42≡6 mod12, 1611≡12 mod13 and 19≡6 mod13; 12^42=1 and 6^42=12, so H_42≠0 mod13. This rules out any fixed prime divisor of the cofactors in this subfamily. Nonexistence of a common divisor does not assert even a single prime cofactor beyond individually checked examples, much less infinitely many.

The seven-adic map is actually an isometry at *every* precision. Set
A=1611^6 and B=19^6. Modulo 49, A≡8 and B≡1. For t>s≥0, binomial lifting
gives v₇(A^(t−s)−1)=1+v₇(t−s) and v₇(B^(t−s)−1)≥2+v₇(t−s).
Subtracting the two affine sums at t,s leaves unequal valuations, so
the cofactor map satisfies

    v₇(H_(6t)/7 − H_(6s)/7) = v₇(t−s).

Consequently as t runs through any 7^k consecutive indices, these
cofactors cover all residue classes mod 7^k, each once. At precision
one the cofactor is 6+2t modulo 7, so v₇(H_(6t))=1 exactly when
t≠4 mod7. This all-depth local coverage coexists with the sparsity of
the exact integer family; it is not global target coverage.

## 3. Source connections and limits

The older [Lonely Runner–Collatz relation packet](../collatz_lrc_relations_20260930.md) emphasizes that cancellation depends on which integer relation is actually in the live fibre. Here the exact multiplier Q_U(3)/Q_U(2) and offset H0 both travel with the moved block; keeping just H would erase them. In the neighboring [local/global packet](../procgen_localglobal_20260930_lrc_collatz_primes.md), a divisor sieve certifies the easy region, not the hard dynamics. Similarly the explicit seven-adic rule certifies an infinite exact residue family, not coverage of prescribed cofactors. [The nearby medíant-tree session](../collatz_mediant_tree_K_20260930.md) is a contrasting arithmetic address tree: determinant-preserving address movement is not a substitute for proving the retained module can be moved through the physical separator.
