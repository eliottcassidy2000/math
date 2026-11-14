# Binary cylinder decoder and word capacity for marked tournament sources

**Proved:** exact prescribed-target test for a given source; a sharp harmonic law for suffix presentations; an obstruction for finite source banks. **Open:** choosing a bank of physical sources that covers all admissible targets. The familiar one-apex/transitive density bound (Session26) is inherited, not rediscovered as a new obstruction.

## 1. Address and count

Freeze a marked ordinal spine with one apex. Let `h>0` be the Hamilton count of the old spine without the apex, `N0` its count with the apex, and `S>0` the total weighted nonempty paths starting at marked (apex-outneighbor) vertices, each weighted by the Hamilton count of its complement. Assume its last ordinal block is a singleton receiving a 0 mark, so the signed start bias is zero. Arbitrary earlier blocks are allowed. The exact cut-state law from [CUT_STATE](../chanpak_selective_target_20261114/CUT_STATE.md) then gives, for every appended transitive word `w`,

    N(w)=N0 + S t + h p(t).

Here the suffix can be empty or end in 0, and `t` is its integer address: bit i (least significant first) equals 1 exactly when the (i+1)st mark is 0. Nonempty words ending in 0 are in bijection with positive integers; t=0 is empty. The empty convention has p(0)=0. The suffix's unrestricted Hamilton count is 1 and the total end weight into 0 marks is t. The remaining coefficient p(t) is the count of the apex over that suffix minus 1.

For a pair i<j (zero-based), the marks 1 at i and 0 at j permit an odd apex cycle of weight two per even interior subset. This contributes 2 if j=i+1 and 2^(j−i−1) otherwise. All complements are transitive and have count one; all nonempty cycle packs have exactly one cycle. Prepending a 0 mark adds no such pair, and prepending a 1 mark adds u plus one extra if the next mark is 0. Thus for u≥0,

    p(2u+1) = p(u),
    p(2u)   = p(u)+u+(u mod 2).

The second equality at u=0 is just the boundary p(0)=0, not a new nonempty canonical suffix. The recurrence implies by induction that p(t) is even and 0≤p(t)≤t; when u odd in the second branch, p(u)=p((u−1)/2)≤u−1.

This is a digital inverse analogue of the [Bernstein–Lagarias Collatz parity-word coordinate](https://websites.umich.edu/~lagarias/doc/bernstein.pdf): both read successive low bits by the shift. Here the arithmetic retains the absolute address quotient and its count. No assertion about Collatz orbit termination or positive Collatz preimages is needed.

## 2. Exact decoder for a prescribed N

Let Y=N−N0. Reject Y<0, and take t=0 if Y=0. For Y>0, cut the address at a power q=2^r (r≥1). Write t=q u+s, 0≤s<q. Let rev_r(s) be the r-bit reversal (including padding), and delta_s=1−bit_(r−1)(s). Applying the preceding recurrence to the first r bits gives the exact split

    p(q u+s) = p(s) + (q−1−rev_r(s))u + p(u) + delta_s (u mod 2).

Indeed, the i-th zero bit of s contributes 2^(r−1−i)u on dividing the address by 2^(i+1). The only adjacency bit crossing the cut is the bit u mod 2 read after a zero top bit of s. Define

    A_s=q S+h(q−1−rev_r(s)),      M_s=Y−S s−h p(s).

Reject M_s<0. The possible remaining tail u is in the explicitly computable integer interval

    max(0,ceil((M_s−h delta_s)/(A_s+h))) ≤ u ≤ floor(M_s/A_s).

This follows from 0≤p(u)≤u. Test equality M_s=A_s u+h(p(u)+delta_s (u mod 2)) on that list. The interval and test are together necessary and sufficient and return literal marked suffixes, not compatible residues. Its real length is at most hY/(qS)^2 + h/(qS). Taking q comparable to sqrt(Y) gives O_h,S(1+sqrt(Y)) candidate tests; each evaluation of p takes O(log(Y+2)) arithmetic operations. One can choose a smaller q when the transported S is large. This is a deterministic decoder for each marked source, not a guarantee that any source works.

## 3. Harmonic capacity of a source

Let rho(t)=sum_i bit_i(t)2^(-i−1), the binary-reversal fraction in [0,1). For t≥1 the recurrence expresses p(t) as a sum over zero bits i in its finite expansion of

    floor(t/2^(i+1)) + bit_(i+1)(t).

Replace floor by its argument and include the infinite trailing zero bits. The floor and adjacency error has size at most 1+floor(log_2 t) and the added infinite tail is at most 1. Consequently

    p(t) = t(1−rho(t))+O(log(2+t)).

Fix h,S and write W(Y) for the number of literal suffixes with N(w)−N0≤Y (including empty). Then

    W(Y) = (Y/h) ln(1+h/S) + O_h,S(sqrt(Y+1)).

To prove this, partition t by s modulo q=2^r. In this class, rho(t) differs from rev_r(s)/q by at most 1/q. Every eligible address satisfies t≤Y/S. On that range N(w)−N0 differs from a_s t by O_h,S(Y/q+log(2+Y)), with a_s=S+h(1−rev_r(s)/q)≥S. Thus membership differs from the interval cutoff Y/a_s only in a boundary interval of the same error size. Both the cutoff and genuine candidates are in t≤Y/S. Counting addresses of this residue in the interval and accounting for the boundary costs O_h,S(Y/q²+log(2+Y)/q+1) per residue. The main term is Y/q times the sum of 1/a_s. Bit reversal permutes the grid points, so this Riemann sum is the integral of 1/(S+h(1−z)) over 0≤z≤1, with error O_h,S(1/q). Choose q comparable to sqrt(Y+1). This proves the claim and does not assume distinctness of counts.

In fact conditioned on membership, rho(t) converges to the distribution with density proportional to 1/(S+h(1−z)); restrict the same sum to any interval of reversed addresses, paying at most two boundary cells. The distribution is not uniform: conditioning on a bounded **absolute** count biases the address digits toward a low slope.

## 4. Consequence for closed banks

A finite bank of actual marked sources, indexed by j, realizes at most

    X sum_j [ln(1+h_j/S_j)/h_j] + O_bank(sqrt(X+1))

distinct positive counts at most X by these right continuations, irrespective of additional cross-source collisions. If the sum is <1/2, a positive fraction of odd targets are absent from this entire bank. This is strictly sharper than summing the individual coarse address bounds X/S_j, since ln(1+z)<z. For example, the strong source consisting of apex over the marked word 10 has h=1,S=2,N0=3; it reaches all odd-compatible congruences but misses at least 1−2ln(3/2) of positive odd targets in lower relative density. The all-order, all-source **one-apex/transitive** bound at least 1−ln 2 omitted is stronger and already in [Session26](../session26_random_bridges_20261022/session26-tournament/APEX_CLASSIFICATION_AND_SPECTRAL_LAW.md). This new bank estimate handles arbitrary fixed cycles and nontransitive cores.

In particular if all h_j/S_j are small, modular controllability by increasingly long words cannot substitute for total capacity from genuinely different sources. Conversely, exceeding the capacity threshold is only a necessary condition, not a coverage proof. Even knowing the exact number of candidate words does not count distinct witnessed targets.

## 5. Changing every cut on one core cannot repair everything

Fix an arbitrary nonempty core R of order m. Form each marked seed by appending R between a leading 1 singleton and a terminal 0 singleton, using any of the 2^m choices U of marked vertices of R. Continue arbitrarily to the right in the same manner. Let Q=P_R(2) be the number of ordered path two-covers allowing empty parts, and h=H(R). If x_U is the weighted total of nonempty paths in R starting in U, the seed has S_U=2(Q+x_U) and x_U+x_(R\U)=Q−h. The first equality follows by splitting paths in the marked seed at the cut after the first singleton, then doubling at the last singleton.

The convex function x↦1/(Q+x), paired at x and Q−h−x, has its maximal sum at the endpoints 0 and Q−h. Hence

    sum_U 1/S_U ≤ 2^m/(4Q) · (1+1/(2−h/Q)).

The OCF two-cover expansion includes the empty pack term 2^m. If R contains a directed triangle, one extra pack contributes 2^(m−1), giving Q≥(3/2)2^m. Also Q≥(m+1)h: each Hamilton path generates m+1 ordered pairs by splitting at a cut, and these are distinct across paths and cuts because concatenating each ordered pair recovers both. Divide the right side of the preceding bound by 1/2 to express the bound relative to odd targets. For every R the result is at most 5/6, since m≥1; for cyclic R it is at most (1+4/7)/3=11/21 because m≥3. In particular, after exhausting every possible choice of cut on one cyclic core, at least 10/21 of the odd targets are still omitted in lower relative density. The sharp harmonic bound from §4 is stronger for each individual core. This is a scope limit on a literal family, not on cut strategies with freely changing cores or more general suffixes.

## 6. No residue can choose even the first branch

For any chosen initial string of r suffix marks, the counts of *longer* suffixes extending this string have the same projected residue image modulo every M as all suffixes of the seed. This includes any prescribed address cylinder t≡s mod 2^r. To see it, append the chosen string and then a forced 0 singleton. By the singleton cut-state transfer, h stays fixed, S is iterated by S↦2S+bh for a mark b, and the signed bias after the forced 0 is zero. Thus gcd(h,S) stays fixed since h is odd. The start bias at each step is 0 or h, so each count increment is a multiple of gcd(h,S). Now apply the [fixed-cycle congruence theorem](../chanpak_cross_problem_theorems_20261114/FIXED_CYCLE_CONGRUENCES.md) to this new literal prefix. Its full residue set is exactly the original set: odd-compatible and N≡N0 modulo gcd(h,S,M). The new family is a subfamily of the requested cylinder.

Contrast the Collatz parity itinerary: one input residue modulo 2^r there determines the first r parities, whereas no single finite target modulus here can even eliminate one prescribed suffix prefix. This is a statement about necessary residue tests, not about high residues used jointly with an absolute interval or the explicit decoder in §2. Different cylinders solve the *same* congruence independently; equal residues are not equal integers.

## 7. Capacity is balanced across every fixed modulus

Let W_a(Y) count suffix presentations of size at most Y such that N(w)≡a modulo M. Write M=2^b m, m odd, d=gcd(h,S,m), and j=2^max(0,b−1) (the number of allowed odd binary classes). For an odd-compatible a with a≡N0 mod d, there is eta_M>0 such that

    W_a(Y) = (d/(m j)) (Y/h) ln(1+h/S) + O_(h,S,M)((Y+1)^(1−eta_M)).

All other classes have W_a(Y)=0. The proof uses the finite-state version of the low-bit recurrence. Generate the binary address from high bits down. Work modulo B=lcm(2,M) and keep the pair (x,z)=(u,p(u)) modulo B, where z is even. Under an independent fair bit beta the transitions are

    beta=0: (x,z) ↦ (2x,z+x+(x mod 2)),
    beta=1: (x,z) ↦ (2x+1,z).

Each state has exactly two incoming choices of the old x under its forced branch, and each specifies exactly one old z of the correct parity. Thus uniform probability on the B²/2 states is stationary. For irreducibility, first go to an even x by appending 0. The blocks 0 and 1^j0 on the even-x layer apply exactly the A and B_j recursions of the [one-apex joint congruence compiler](../chanpak_least_missing_jacobian_20261114/MODULAR_ESCAPE.md) to (x,1+z). That compiler sets any admissible joint state from any initial even x and even z: the optional low-bit toggles and forced anchor do not require a specific initial pair, and the odd-modulus moves apply once anchored. Thus all even-x states communicate. To reach an odd-x state, prescribe the required low b' bits of x by appending a fixed b'-bit word, with b'=v_2(B); modulo the odd m the initial x is freely chosen because the multiplier 2^b' is a unit. The effect on z is addition of a known even increment once the initial x is fixed, so the initial even z can be chosen to place the desired z. The state (0,0) has a self-loop, proving aperiodicity. Hence the finite chain mixes exponentially, say with total-variation error O_M(lambda_M^l) after l fair bits for some lambda_M<1.

For the absolute restriction, pad t to L bits with 2^L comparable to Y/S and with t≤Y/S included. Read the top r and bottom r bits, where r=floor(L/4), before sampling the middle bits. These ends determine t up to O(2^(L−r)) and rho(t) up to 2^(-r). Thus the ends alone determine the truth of N(w)−N0≤Y except for an uncertainty band of width O_(h,S)(Y2^(-r)+log(Y+2)) in the actual value. By §3 this band contains O_(h,S)(Y^(3/4)+sqrt(Y+1)) words: apply the asymptotic at the two actual interval bounds. In each block of fixed ends the middle bits mix for L−2r steps. Reading the fixed bottom bits, the residue x is now fixed on the binary part (when r≥b') and uniform modulo m. The uniform even z remains uniform and independent of this x: it is merely translated by an even increment determined by the old x. The total-variation error cannot increase during this fixed readout. Since the count is N0+Sx+hz, h odd and S even give uniform odd-compatible binary classes, and the linear projection on the odd part gives uniform classes a≡N0 mod d. Sum over the fixed-end blocks; the mixing error is O_M(2^L lambda_M^(L−2r)), not multiplied by an extra boundary count. Combine the two power savings with §3.

Consequently a congruence-complete finite source bank in which every gcd(h,S)=1, but whose sum of harmonic capacities is less than 1/2, omits a positive lower proportion of the odd targets in **every individual odd-compatible progression modulo every fixed M**. The word count above is still not a distinct-value count. For multiple banks with differing gcds, use the corresponding d/(m j) multiplier only in their allowed progressions. No claim is made for a varying modulus.
