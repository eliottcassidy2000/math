# A full physical offset alphabet at slope 99

**Theorem.** For every integer m≥1, each of the 99 consecutive odd integers

    99·2^m+1, 99·2^m+3, …, 99·2^m+197

is the Hamiltonian-path count of a strong tournament. In particular the previously discussed target 99·2^m+23 has three distinct strong realizers from the constructions below, not just a match in residue. This is a repeated band of fixed absolute width, not a cofinite spectrum theorem.

## 1. The phase bank and its physical decoder

Use the [one-return star](../session324_busch_inverse_and_gap_structure_20261113/FERMAT_SEVEN_RAIL.md): form I→U→W and add apex v, which loses to all of U and beats marked sets in I and W. In I let h=H(I), let α count ordered two-covers whose second part starts in the marked set, and let λ count those covers in which the first part also ends in the marked set. A path with an endpoint requirement is nonempty. In W let k=H(W), d=H(W+v), and β count the ordered two-covers whose first part is empty or ends at a vertex beating v. Other slots are allowed empty. The exact law is

    H(star)=αβ Q_U(2)+(h d−k λ)H(U).

Indeed split a full path at the apex. If both parts meet U, the I suffix supplies α, the legal W prefix supplies β, and U has two nonempty parts. If U lies wholly before the apex, the nonempty I suffix again gives αβ, and the empty I suffix gives hd. If U lies wholly after the apex, the I suffix must be nonempty; excluding a bad I prefix when no W prefix is present subtracts kλ from αβ. Adding the three cases and using Q_U(2)=F₂(U)+2H(U) proves the formula.

The [bank](corridor_witnesses.json) gives literal outgoing masks and cuts for each odd c from 1 through 197. Its finite identities are

    αβ=99·2^r_c,  h d−k λ=c,  0≤r_c≤4.

At U a singleton, each template has exactly one nonsingleton SCC, containing that singleton. All other SCCs are singleton padding. This is a finite structural condition on the actual mask, checked along with the profile identities in the replay. Inflating the marked singleton by any nonempty U preserves the strong core: the quotient of that core is strong, so every vertex of U can reach every other through the exterior, and the remaining SCCs are still singleton padding. Deleting those singleton SCCs does not change H. Now set U=TT_(m−r_c) when m>r_c. Its count is one and Q₂ is 2^(m−r_c), so the resulting strong count is exactly 99·2^m+c. The 115 remaining prescribed targets (m≤r_c) have explicit strong masks in the same bank. The search is not a proof of bank minimality, and no search cutoff is used for its infinite tail.

The first two bands adjoin, covering every odd number from 199 through 593. The bands for later m need not adjoin. These facts are about this supported band's entry, not the least known global cutoff.

## 2. The offset can be chosen before the physical reservoir

A stronger fixed-interface reading is available. Let a nonempty tournament U of profile (H,Q)=(H(U),Q_U(2)) have at least five singleton SCCs. For each c remove r_c of these singleton SCCs before inflating its chart. Then H is unchanged and Q becomes Q/2^r_c. Thus, for **every** allowed c, a strong tournament has the exact count

    N=99Q+cH.

This is one uniform physical profile alphabet: it does not require H=1 or U transitive. The offsets typically give spacing 2H, not two.

For a prescribed odd N and a proposed H coprime to 99, choose the unique odd c in {1,3,…,197} satisfying cH≡N mod 198. The count now forces an exact even slot budget

    Q=(N−cH)/99.

So for this alphabet one can remove the final residue/phase ambiguity completely. What remains is a lift to an actual U with that exact (H,Q), enough singleton padding, and Q>0. A residue solution for Q is not enough.

There is an explicit obstacle to trying finitely many fixed H. For each H, the chain of nonsingleton strong components has finitely many possibilities by Busch's bound and the SCC product (see [S324's fixed-count fiber](../session324_busch_inverse_and_gap_structure_20261113/FIXED_COUNT_FIBER.md)). After fixing those components, singleton entries only multiply Q by powers of two, independent of where they are placed. The obtainable Q for fixed H is therefore a finite union of dyadic rays. No finite fixed H bank plus this digit alphabet can cover a positive-density set of odd targets. The remaining inverse has to choose an unbounded profile, or route the omitted targets to a different class.

## 3. Three genuinely different cuts at 99·2^m+23

First use I=TT₂, with v beating both I vertices. Its (h,α,λ)=(1,3,2). The two W choices below give the same star state (αβ,h d−kλ)=(99,23):

| W rows | v's outcut in W | (k,β,d) |
|---|---:|---|
| (0,9,3,5,15) | 22 | (3,33,29) |
| (8,1,3,6,15,31) | 62 | (5,33,33) |

Both stars are strong without removing any padding: v enters I, I reaches U and W, U beats v, and all of W can reach an in-neighbor of v. They give exactly 99Q_U(2)+23H(U) for **every** nonempty U. Their sizes differ by one.

The third chart uses [the unit-end Busch construction](../session324_busch_inverse_and_gap_structure_20261113/UNIT_END_GRAMMAR.md) with terminal component C a singleton and prefix cut word 1010 followed by m copies of 1. Its bits record y→prefix vertex in forward order. The prefix has outgoing relaxation L_P=11·2^m, ordinary count H_P=9, and end E_P=0. The only descending pairs of the bit word are (1,2),(1,4),(3,4), of weights 2,4,2 in the one-apex formula; the suffix adds none. The full unit-defect source therefore has

    H_A=3L_P+H_P−E_P−2=33·2^m+7,  D_A=1.

The right triangle at the old root then gives 3H_A+2=99·2^m+23. This Busch graph has order m+9. The W6 star has the same order but just one vertex each of outdegree 1 and n−2; the Busch graph has two of each. In the Busch graph the new triangle's two off-root vertices are the degree-one vertices, and the first two prefix vertices are the degree-(n−2) vertices. In the W6 star the first I vertex is the unique vertex of degree n−2, and W vertex 1 is the unique one of degree one. Hence the two graphs are distinct even after taking converse. The W5 star is smaller.

**Orientation correction to S324.** In the single-apex forward-chain inversion formula the word x records that the apex beats the selected chain vertex. In the notation u_i=[y→i] of S324's unit-end transfer, the weighted inversion count uses x_i=u_i, not x_i=1−u_i. The exact transfer there retains the actual H_P and its replay counts H_P directly, so the transfer and classification are unchanged. The corrected orientation matters in choosing the 1010 prefix here.

## 4. A second, completely elementary three-lane crossing

Use C a singleton again and take the entire b-bit prefix to be ones. Now L_P=2^b, H_P=1, E_P=0, yielding H_A=3·2^b−1 and D=1. Attaching the right R_j realizes 3(2^j+1)·2^b−1. The fully transitive selective chart realizes (γ3^u+6·2^u−3)·2^b−1 with source size u. Within these parameter choices, the common absolute targets are exactly the rails from (γ,u,j)=(2,1,2),(6,1,3),(2,3,5): slopes 15,27,99 and intercept −1. Different b cannot agree by the valuation of N+1, since both coefficients are odd.

The coefficient test is 2^(j−1)=3^(u−1)+2^u−1 for γ=2, and 2^(j−1)=3^u+2^u−1 for γ=6. For odd e>0, v₂(3^e−1)=1; for even e, v₂(3^e−1)=2+v₂(e), by factoring each successive square (the new factor has valuation one after the initial square). For γ=2 and u≥2 this valuation is less than u except u=3; consequently the right side cannot be a pure power of two except that case, which gives 16. For γ=6 the valuation is less than u except u=2,4 among u≥2; their right sides are 12 and 96 rather than powers of two. Check u=1 directly. These meetings do not contradict [the disjointness of the two fully transitive selective charts themselves](EXPONENTIAL_OVERLAP.md), because the third route is Busch's unit-end chart.
