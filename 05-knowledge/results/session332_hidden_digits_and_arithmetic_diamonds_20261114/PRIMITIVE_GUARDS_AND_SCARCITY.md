# Exact gcd guard closure, but a thin old/first language

**Status:** proved all-order physical divisibility guard and absolute sparsity under infinitely many allowed Fermat letters. No unrestricted spectrum exclusion.

## The first/old sidecars cannot develop different gcd guards

For input (H,h,L,E), choose the old root O, first interior I₁, terminal T, or any other interior I_k. Use the corrected S331 all-port transfer. The gcd of the new four entries is exactly

    O: gcd(H,L,h)
    I₁: gcd(H,L,E)
    T or I_k with k≥2: gcd(H,L).

At O, h'=E'=h and L'=2·2^j L−h, so odd divisibility of the new tuple is precisely old divisibility of H,L,h. The terminal assertion follows from h'=E'=H and L'=H'. At I_k for j≥2, h'=H+2^(j−1)L and E'=2^(k−1)L, so H', h', E' expose the same odd divisors as H,L when paired with H'. The first interior adds −E in L'; all other interiors do not. At j=1, the first interior has E'=L and L'−h'=H+L, while H'=H+2L. These identities give the same gate gcd(H,L,E). All the gcds are odd since H,H' are odd.

More importantly, **after every move**

    gcd(H',L',h') = gcd(H',L',E') = gcd(H',L',h',E').

This is immediate at O,T from h'=E'. At j≥2 interior it uses precisely the H',h',E' combinations above, and at j=1 interior the difference L'−h'=H+L exposes the same primes and prime powers as E'=L. Thus after a primitive move, every further arbitrary sequence of O or I₁ remains primitive, with no exact residue test, regardless of genuine length choices. Terminal or other interior can only preserve primitivity if the pair (H,L) is coprime. The claim is about the physical four-state and not about whether a primitive count is prime. It is a native version of a shared guard closure; see [S323](../session323_computability_and_native_descent_20261113/GUARD_CLOSURE_AND_HALTING.md) for why guard compatibility alone does not prove a computation halts or reaches an intended target.

The pair gcd itself is visible by integer elimination:

    O: gcd(H',L') = gcd(H+2^j L, 2H+h)
    I₁: gcd(H',L') = gcd(H+2^j L, (2^(2j)−2·3^(j−1))H−2^j E).

The second formula uses multiplication by 2^j, a unit modulo any divisor of H'. At interior k≥2 put t=j−k and assume gcd(H,L)=1. Then gcd(H',H)=1, so the same elimination gives

    gcd(H',L') = gcd(H', G_t),  G_t = 2^(2t+1)−3^t.

For p≠2,3, a prime can divide G_t exactly when (4/3)^t=1/2 mod p, an exact discrete-log condition on the port distance. Three never divides G_t (check t=0 separately), seven never does ((4/3)≡−1 and 1/2≡4), and eleven never does (4/3≡5 generates {1,5,3,4,9}, not containing 6). For five precisely t≡1 mod4 is bad. Therefore such interior selectors, regardless of length, never contaminate a primitive pair at the three small axes 3,7,11. They can contaminate it at other primes.

## Nevertheless absolute outputs from this protected language are sparse

**Theorem.** Start from any finite bank. Allow every genuine exponent j=1,2,4,8,…, and at each step select only O or the first interior. Then the count of distinct Hamiltonian outputs ≤X is O(X^σ) for some σ<1. No primitivity assumption is needed for this upper bound.

At any state L≥H and E≤H,h≤H. After the first step in this subgrammar h/H≤2/3 and E/H≤1/2. For the old row this follows because h'=E'=h and H'≥3H. For a first interior j≥2 use h'=H+2^(j−1)L and E'=L; j=1 only subtracts E from the displayed deletion count. Thereafter the same bounds persist.

For the next three letters let d_i=2^j_i and c_i=2·3^(j_i−1)≥d_i. Normalize the input H=1; write the initial (h,L,E) as otherwise variable with the preceding bounds. We show H₃ is minimized termwise at the universal anchor (H,h,L,E)=(1,2/3,1,1/2). Put H₁=H+d₁L. For the first choice O, L₁=2d₁L−h, h₁=E₁=h; for first I, L₁=d₁H+c₁L−E, E₁=L, h₁=H+(d₁/2)L−[j₁=1]E. The final selector does not affect H₃. For second choice O or I the linear expressions are respectively

    O: H₃=H₁+d₂(1+2d₃)L₁−d₃ h₁;
    I: H₃=(1+d₂d₃)H₁+(d₂+c₂d₃)L₁−d₃ E₁.

Substitution shows nonpositive coefficients of h,E and nonnegative coefficient of L in every case. In the only subtractive L cases, first I then O subtracts d₃(d₁/2) against the larger coefficient d₂(1+2d₃)c₁; first I then I subtracts d₃ against the larger coefficient (d₂+c₂d₃)c₁. The sign for the I,O E term is −d₂(1+2d₃)+[j₁=1]d₃<0. This proves the uniform anchor bound.

For all 216 three-letter words restricted to exponents {1,2,4}, summing reciprocal expansion at that anchor gives 0.958050926545...<959/1000, checked using exact rational arithmetic in the replay. For words containing any longer exponent use the simpler product lower growth ∏_i(1+d_i). The full alphabet weight ∑_j2/(1+2^j) is below 6/5, and the tail from j≥8 is less than 79/10000, using 2/257+2/65537+2^(2−32). By the union bound the entire long-word weight is at most 3·(79/10000)·(6/5)². Thus the three-letter reciprocal sum at exponent one is strictly below one. Continuity in σ near one follows by uniform geometric domination of the explicit tail; choose σ<1 so the three-letter weight at σ remains below one.

Partition any word into its first move, three-letter blocks, and a tail of length at most two. The first and tail have finite inverse expansion sums because each step grows H by at least 1+2^j. Sum the geometric series over three-letter blocks. Any physical output ≤X costs at least a constant times X^(−σ) in this sum, hence there are O(X^σ) such words and no more distinct outputs. This is not a claim about the entire all-port or all-tournament grammar. It proves that removing every modular obstruction on the protected old/first lane is still not the required absolute-target bridge.

## A primitive source can nonetheless be dead before the first move

The balanced law is an **image** property. It is not satisfied by all strong primitive rooted tournaments. The literal order-seven rows

    (6,24,42,81,69,27,39), root 1,

have (H,h,L,E)=(135,33,165,25). The full gcd is one, but the old port forces a factor 3, the first interior forces a factor 5, and every other port forces both. No exponent changes any of these new full-state divisors. A complete scan at orders at most six finds no strong primitive dead root; the displayed witness therefore has minimum strong source order in this sense.

This is an unbounded physical obstruction, not an isolated small example. Replace vertex 0 by a transitive tournament TT_m, with uniform exterior directions and root still the old vertex 1. The blow-up is strong since its quotient is strong and every block is nonempty. The only exterior out-neighbors of the replacement block are the old vertices 1 and 2. In a path after deleting vertex 1, or a path ending there, an exit between two separated visits must use vertex 2. So at most two visits are possible, and the internal filling counts are 1 and 2^m−2. The direct m=1,2 controls determine

    h(m)=33+16(2^m−2)=16·2^m+1,
    E(m)=25+16(2^m−2)=16·2^m−7.

Their difference eight ensures gcd(h,E)=1 for all m, hence the overall state is primitive independent of all other factors.

For each of H,h,L,E decompose a witness into its itinerary through the fixed six other vertices and separate visits to the clone block. The number of visits is bounded by the fixed exterior size plus the number of path slots. All start/end guards depend only on which module is adjacent to the marked vertex. For k visits the number of internal fillings is the number of ordered nonempty cuts into k transitive-subset paths, exactly

    F_k(TT_m)=Σ_(i=0)^k (−1)^(k−i) binom(k,i)i^m.

Thus every measured coordinate is an integer linear combination of finitely many i^m. For m≡1 mod4, i^m≡i mod15 even in the nonunit cases modulo 3 or 5. So all four coordinates reduce to the m=1 ones modulo 15. Consequently at every order m+6≡7 mod4 with m≥1 there is a **strong primitive rooted source outside the image of every Fermat port**, and every allowed move from it locks in 3 or 5. No later count can be prime or 7p: inherited divisibility persists, while H is already at least 135 (path insertion on the induced base). This exhibits an exact separation between primitive *measurements* and physically reachable primitive *guards*.
