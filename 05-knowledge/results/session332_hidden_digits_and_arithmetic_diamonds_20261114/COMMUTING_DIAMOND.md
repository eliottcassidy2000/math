# Two different strong constructions with exactly the same ready state

**Status:** proved infinite all-order multiplicity, exact uniqueness of the elementary commute, and a two-slot observer that opens the hidden collision. This does not supply an arbitrary prescribed count.

Use the corrected [all-port Fermat transfer](../session331_returns_sections_and_switch_surfaces_20261114/FERMAT_INTERIOR_PORTS.md). Write (H,h,L,E), where h deletes the marked root, E counts full paths ending there, and L is its relaxed split response. At length one denote by O the old-root choice, by I the first (also penultimate) interior choice, and by T the terminal. A ready state is (H,e,H,e), with e<H and a **literal ready terminal** has outdegree one. The j=1 rows are

    O(H,h,L,E)=(H+2L, h, 4L−h, h)
    I(H,h,L,E)=(H+2L, H+L−E, 2H+2L−E, L)
    T(H,h,L,E)=(H+2L, H, H+2L, H).

Read each word chronologically. From any ready state both OI and IO have the same height and split response (H₂,L₂)=(11H−2e,14H−3e), but the omitted coordinates (h₂,E₂) are (7H−2e,4H−e) and (2H−e,2H−e), respectively. T discards those omitted coordinates and returns the entire identical state

    (39H−8e, 11H−2e, 39H−8e, 11H−2e).

The IOT scalar branch (39H−8e) and its 109 seed already occur in S331; the parallel OIT branch, equality of the complete return, degree decoding, and observer opening are the new content here. Because the next state is ready again and the whole word adds six vertices, this is a genuine repeatable diamond. Every composition of such words may make its choice independently without changing future active values. It applies to physical marked tournaments and preserves strongness.

## No other elementary old/first commute from ready is available

For arbitrary lengths j,k≥1 put d=2^j and c=2^k. Let I_k always denote the *first* interior, with its −E correction. The difference in the Hamiltonian outputs of the two length-two paths O_j I_k and I_k O_j is exactly

    (d−c)(H+e) + d(c−2·3^(k−1))H.

This follows by projecting the first move to (H',L') and using H''=H'+d_second L'. If k≥2, then 2·3^(k−1)−c≥2. For d≤c the difference is strictly negative. For d>c the positive term is less than d(H+e)≤(4/3)dH for literal terminal-ready input e/H≤1/3, so the difference remains strictly negative. For k=1 the second term is zero and the first is zero exactly when d=c=2. So the all-length identity is isolated precisely at the triangle old/first pair. This is an assertion about the two named moves from readiness and first-coordinate equality before the terminal, not a global word-relations classification.

## Unrooted isomorphism genuinely distinguishes the histories

Start with the physically ready triangle C3: (H,e)=(3,1), with any marked vertex. At each new diamond call the three pairs (a,b),(c,d),(f,g), with b,d,g their terminals in time order. In OIT the new pair degrees are (5,5),(2,3),(1,1); in IOT they are (3,5),(3,3),(1,1). To verify, each pair vertex has outdegree one just after addition, gains two outgoing edges for each later pair, and loses one of those two exactly when it was the attachment root of that later pair. The initial marked root v was used twice in OIT and once in IOT. It had old degree one, so it now has degree five or six, respectively. Every other old vertex has gained six and, since the old tournament was nontrivial strong, has degree at least seven.

The presence of degree two detects the last letter OIT. In IOT the old vertex set is exactly the degree-at-least-six set, and its mark is the unique degree-six vertex. In OIT the old nonroot set is exactly the degree-at-least-six set. The remaining mark is the unique degree-five vertex beating any member of that set: the old mark had one old out-neighbor, whereas the new first pair's only possible return into the entire old tournament is to the old mark, not an old nonroot. Thus the induced previous *marked* tournament can be recovered from the unrooted output, and iterating the inverse recovers the entire history. At order 3+6m there are 2^m pairwise nonisomorphic strong graphs with one and the same exact ready pair. No claim of separation up to converse is needed.

Their common count obeys H₀=3, H₁=109, H_(m+2)=37 H_(m+1)−10 H_m; this follows by Cayley–Hamilton from ready matrix [[39,−8],[11,−2]]. The next count is 4003. The pairs remain primitive: from primitive odd (H,e) the only possible new common factor is five (determinant ten), and it occurs exactly when e≡8H mod5. The initial pair avoids this line. Under any allowed step the new ratio e/H modulo five is four (check that the second matrix row is four times the first modulo five), avoiding it permanently. Outputs increase strictly, but this remains one sparse exact ray, not a cofinite family. The geometric payoff is many graphs per exact H and per exact future observer state, not many distinct counts at that order.

## An exact two-slot observer opens the same hidden collision

For an arbitrary terminal-ready input A rooted at v, put C=Q_(A−v)(2), Q=Q_A(2), and R± the signed rooted responses of S328 at z=2. View the six new vertices as a fixed rooted right context B sharing v: all A−v dominate all new vertices. In the [general right attachment formula](../session328_fermat_gluing_and_phase_repair_20261114/FERMAT_CLOCK_MEMORY.md), the needed context data (Q_(B−v),Q_B,E_B^+,E_B^−), where E_B reads paths *ending* at v, are

| context | data at z=2 |
|---|---|
| OIT | (128,596,395,87) |
| IOT | (256,644,395,183) |

They are verified on literal seven-vertex buffers starting from a singleton root. Equivalently, compute E± by the starting observer on the converse. Substitute in the right attachment formula (4) to obtain

    Q2(A via OIT)−Q2(A via IOT) = −128Q+256R⁺+64R⁻−112C.

Since v has exactly one out-neighbor, the identity Q/2=R⁺+R⁻−C reduces this to 48(3C−4R⁻). Suppose A came from a terminal selection of length j on U rooted at u. Then C=2^j Q_U(2), since deleting the terminal leaves U followed by a transitive chain. By the unsigned interpretation of R⁻, an optional distinguished path in A−v may be empty or end anywhere except the unique out-neighbor u. The covers excluded for ending at u put the entire new chain into the other slot, so their count is precisely B_u(U), the number of two-slot covers of U whose distinguished slot is nonempty and ends at u. Consequently

    Q2(A via OIT)−Q2(A via IOT) = 48(4 B_u(U)−2^j Q_U(2)).

Slot swapping gives 2B_u(U)≤Q_U(2), and this is strict when u has an out-neighbor: take the path u→x in one slot, and any Hamilton path on the remaining vertices in the other; u ends neither slot. In particular the difference is strictly negative for every terminal constructed from a nontrivial strong source. For the triangle seed its predecessor is a singleton at j=1, so equality at depth one is expected. At depth two the four exact Q2 readings are 1671956,1698644,1631444,1671956, all with H=4003. This is pointwise visibility of the final diamond given its predecessor, not injectivity of Q2 on the full word history.

The replay directly computes all counts and the observer constants on literal low-order graphs, and exercises the inverse degree decoder through depth five. The proofs supply the unbounded conclusions.

## One lacunary exact ray still covers every scalar 7-adic address

Let D=[[39,−8],[11,−2]] and x=(3,1), so the common exact ready state is D^m x. Modulo seven the D orbit runs through all 48 nonzero pairs. The replay verifies that finite field statement and

    D^48 ≡ I + 7 [[3,5],[1,5]] (mod 49).

For any desired H modulo seven, choose a nonzero pair with that first coordinate and with nonzero derivative 3H+5e modulo seven; there is always a choice because the e coefficient is a unit. Let m be its first occurrence in the base orbit. For each t≥1, changing m by a multiple c·48·7^(t−1), c∈{0,…,6}, changes the next digit of H by c times that same unit, and fixes all earlier digits. This follows from D^(48·7^(t−1))≡I+7^t[[3,5],[1,5]] modulo 7^(t+1), by binomial lifting. Consequently this single positive exact recurrence realizes every scalar residue modulo each 7^k. It nevertheless has only O(log X) exact counts through X by uniform growth, even while its exact counts have exponentially many strong graph realizations.

There is also an exact zero-phase law. The first-coordinate zeros modulo seven occur precisely at m≡3 (mod 8): D^8 is scalar 3 modulo seven and only m=3 in the first eight has zero first coordinate. In the localized ring with denominators prime to seven let Y=D^8/3. The replay computes

    Y ≡ I+7B (mod 49), B=[[3,2],[6,1]];
    D^3 x ≡ (0,4) (mod 7).

Put f(ℓ)=(Y^ℓ D^3 x)_1/7. The denominator is legitimate since the first coordinate is always divisible by seven, and the factor 3^ℓ removed from H_(3+8ℓ) is a unit. For any r≥0,

    f(ℓ+c·7^r)−f(ℓ) ≡ c·7^r (mod 7^(r+1)),

because Y^(7^r)≡I+7^(r+1)B mod 7^(r+2), and the first component of B(0,4) is one modulo seven. Hence f is a bijective digit selector modulo every 7^a, and v7(H_m)≥a holds exactly on one phase modulo 8·7^(a−1). The first lift of the zero phase is m≡35 (mod 56); the literal recurrence gives v7(H_35)=4. This is an exact scalar modular success and an exact absolute sparsity example in the same grammar, not a surjection onto all large odd integers.
