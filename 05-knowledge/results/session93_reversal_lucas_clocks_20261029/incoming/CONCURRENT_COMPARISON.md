# Concurrent reflection towers: provenance and a composite-clock control

**Status:** source comparison; independently checked finite composite/port controls below. The source packet arrived during Session93 publication at commit
`4368e3847072d9a5d88ddb7d80470f9861c97b93`, under
`05-knowledge/results/lucas_prime_reflection_towers_20261028/`.
Its source-reported general extensions are not silently counted as this session's discoveries.

## 1. Overlap and distinct contributions

The concurrent packet proves the strict prime-Lucas reflection classification and full recursive centralizer holonomy. These were also reconstructed and independently audited in this session. It supplies the general identity H=C∩<K>, a constructive loop-length bound, all-r port moments, prime-divisor valuation records, and a period-nine composite chart.

Our exact minimum k(k−1)/2 **individual triangle generators**, pointwise port-base minimum, complete LCA/orientation tuple invariant, blue-quotient/certificate distinction, exact trace-clock kernels, and exceptional12-state decoder minimum have their own proofs and scopes. The finite Jordan certificates A29/S199/A521 are supplied here explicitly. No priority competition or equivalence between different loop-cost measures is inferred.

The general H=C∩<K> proof is consistent with our examples: square-root edge labels belong to<K>; based triangle corrections are precisely Schreier stabilizer generators for its conjugation action on marks. Paley7 has<K>=C7 and C=C3, so H is trivial; in the recursive family our constructive basis generates all C. This is a proof comparison, not a claim that an arbitrary odd group has H=C.

## 2. Port moments and a specified probability law

The source `PORT_MOMENTS.md` views a uniform A_k automorphism through its number X_k of fixed vertices. Its root rotation is zero with probability1/3; then its three child automorphisms are independent and uniform. Therefore

    X_0=1,
    X_k = 0 with probability2/3,
          X_(k-1)^(1)+X_(k-1)^(2)+X_(k-1)^(3) otherwise.

This is a critical0-or-3 Galton–Watson population law. Burnside identifies its rth raw moment with the number of ordered r-port orbits. The centralizer product decomposition gives fixed-count variable

    Y_k = 1 + 2 sum_(j=0)^(k-1) X_j,

with independent X_j in this expression. In particular

    E Y_k=2k+1, E Y_k²=8k²+1,
    Var(Y_k)=4k(k−1).

This is exactly the excess of joint two-port orbit count over the product of marginal orbit counts. The sampling law is uniform group elements, not uniform arithmetic sources, random Collatz trajectories, or a quantum measurement law.

As an independent finite check, `check_concurrent.py` enumerates the literal centralizers through k=3 and compares the first three fixed-count moments with the recurrence. It obtains the third moment

    48k³−40k²+18k+1.

The source's full all-r polynomial theorem is credited rather than reconstructed here. The complete tuple invariant in our port note explains what distinguishes these joint classes beyond their counts.

## 3. The period-nine ring is not a cyclic76-clock

At n=9, the inherited odd-word reader has76 states. Since

    phi^9−1 = 20+34phi,

its quotient is

    Z[phi]/(phi^9−1) ≅ F4 × F19,
    (u,v) -> (u+v phi modulo2, u+5v modulo19).

The coordinate ideals are comaximal because2 is a unit modulo19, so CRT gives surjectivity. Both displayed coordinates kill the ideal generator and the quotient has index76, so the map is an isomorphism. Its additive group is C2×C38, not C76. Addition by1 has order38. A full-cycle prime-degree theorem therefore cannot be applied by treating this ring translation as a76-cycle.

Let R reflect cyclic word positions and let T be addition by1 through the reader. Reflection acts on the F4 coordinate by Frobenius x→x². For each nonzero F19 coordinate y, its four possible F4 tags give four distinct reflected F19 coordinates. This source claim is checked independently from all76 legal words, not from a chart receipt.

Consequently the F19 coordinate alone is not reflection-closed. More precisely, the observation profile

    (pi19(x), pi19(Rx), pi19(RTx))

distinguishes all76 states. For y≠0 the second coordinate already distinguishes the four tags. For y=0, T sends the four states bijectively to the y=1 fibre, where the third coordinate separates them. Thus the coarsest future-observation quotient for pi19 under these named operations is the full76-state carrier. The missing F4 tag is exactly two bits, not a free conjugate choice.

This is a concrete composite counterpart to the trace-clock theorem, but the observations and operations differ: reflection of decoded words is not the linear trace recurrence. We do not identify their kernels merely because both need a sidecar.

## 4. Prime records versus prime-valued terms

There is an exact conversion between its projective clock and our norm-one clock. Over an odd unramified prime, let h=ord(a) for a=phi². The eigenvalue ratio of the Fibonacci matrix is theta=phi/bar(phi)=−a. Since the norm-one group is cyclic,

    ord(theta)=2h if h is odd;
    ord(theta)=h/2 if h≡2 modulo4;
    ord(theta)=h if4 divides h.

For even h its unique involution is a^(h/2)=−1, so theta=a^(1+h/2) and the formula follows by a gcd. For odd h, −1 lies outside the odd cyclic subgroup and doubles the order. At11 and29 the orders5,7 become10,14; at47 order16 remains16, locating the Lucas midpoint at8. These are different distinguished clocks, not contradictory order computations. The checker verifies this conversion for every odd unramified prime through100 by direct ring powering.

The concurrent `PRIME_ATLAS.md` retains projective order d_r, initial valuation e_r, and splitting character chi_r. It proves the Lucas-divisor midpoint and valuation tower and cites Moree's density2/3 for primes dividing some Lucas number. The primary [Moree paper, Theorem2](https://www.numdam.org/item/JTNB_1996__8_2_449_0.pdf), was reopened here and explicitly gives the ordinary Lucas case P=1, crediting Lagarias. That is a different population from indices at which L_n itself is prime. We read the exact source scope; we do not reproduce its whole prime5000 census or claim a new density theorem.

The atlas also distinguishes three actions at eleven: word reflection fixes three words, J_6 fixes zero of twelve matching frames, and the Lucas projective half-turn fixes two points. These are not interchangeable involutions. Our direct finite-field order controls and the source's projective records are complementary rather than competing definitions of one clock.

## Reproduction

`python3 check_concurrent.py` and `python3 -O check_concurrent.py` use exact standard-library arithmetic, with no producer imports. They independently enumerate the period-nine words and CRT action and the small centralizer moments. Source-only all-r and prime-density extensions remain explicitly separated from these finite replays.
