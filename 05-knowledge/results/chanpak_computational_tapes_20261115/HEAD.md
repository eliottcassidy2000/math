# A physical head with two exponentially separated registers

Let a nonempty tape be an ordinal sum U = D_1 → ... → D_m of nonempty tournaments (arbitrary internal structure). For a chosen block D = D_i write U=A→D→B and require both A and B nonempty. Select a nonempty subset V of D. Add a last sentinel t with U→t and a control vertex x with t→x; x beats U except for V→x. Denote the result by W(i,V). Since x beats A and every tape vertex beats t, all vertices are in a single strong component (each tape vertex reaches x via t, and x reaches the tape via A).

Write h(D) for the Hamiltonian path count and q(D)=Q_2(D) for the number of ordered pairs of possibly empty directed paths partitioning D. Define:
- L(D,V): ordered two-slot covers in which the second slot starts in V, with no constraint on the first slot.

Let E count covers whose first slot ends in V, and J count those among E whose second slot starts in V. Then the exact readout is

`H(W)=h(U)+q(U)+h(B)[q(A)E(D,V)-h(A)J(D,V)]-h(A)q(B)L(D,V).` (1)

Proof: classify Hamiltonian paths by the predecessor of x. With x first, the whole tape must be traversed before t and starts in A, giving h(U). With t preceding x, distribute U into two ordered slots, both potentially empty. Exclude the q(B)L(D,V)h(A) choices whose second slot starts in V: all A is forced into the first slot and B is arbitrary. If a vertex in V precedes x, the first slot receives no B or t and must end in V. Both halves of A are unrestricted; B and t enter the second slot in that order. If the second half of D is nonempty and starts in V, it is forbidden precisely when no A is present after x. All A is then in the first slot, hence the factor h(A) on J rather than q(A).

### Adjacent head windows do not overlap

Set e=E/h(D), j=J/h(D), l=L/h(D), and r_i=q(D_i)/h(D_i). Put p_i=product of r_a for a<i, s_i=product of r_a for a>i. Relative to the constant h(U)+q(U), and normalized by h(U), (2) is

`F_i(V) = p_i e - j - s_i l.`

The positive term minus j fits the band `[p_i, p_i(r_i−1)]`:
- E−J≥h(D), by cutting each Hamiltonian path immediately after its last vertex in V, an injection;
- E≤q(D)−h(D), since the first slot cannot be empty;
- J≥0.
In fact p_i e−j=(p_i−1)e+(e−j)≥p_i since e≥1. Thus it lies in `[p_i, p_i(r_i−1)]`. Each nonempty D has q(D)≥2h(D) by the two placements with an empty slot.

Likewise L≥h(D) by cutting a Hamiltonian path just before its first V, and L≤q(D)−h(D) because its second slot is nonempty. So s_i l lies in `[s_i,s_i(r_i−1)]`. When the head advances one block, p_{i+1}=r_i p_i and s_i=r_{i+1}s_{i+1}. The positive band moves strictly up; the band being subtracted moves strictly down. Therefore for every pair of nonempty masks V,V' on consecutive permitted blocks,
`H(W(i+1,V')) − H(W(i,V)) ≥ h(U)(p_i+s_{i+1}) > 0.`

This is a physical strong-tournament obstruction, not a restriction on abstract tape machines. Changing every unmarked module to one with the same (h,q) keeps the band boundaries unchanged. A target, if it belongs to this architecture, determines at most one head location by the band inequality; the signs of the two registers prevent the right sentinel from compensating the left register. It leaves open simultaneous marked blocks, a compensating change in the tape's (h,q), or a nonuniform contact that changes the sign geometry. Those are precise ways to revise the architecture rather than try more masks.

For a transitive tape of n vertices with both end guards and singleton head k (zero-based), E=L=1 and J=0. The explicit count is `1+2^n+2^k−2^{n−k−1}`. This is a physical strong two-sided readout at every n≥3, but its counts do not form a neutral head return.


## Multiple marks: an exact tape and a saturated quaternary address

Now let each block i have its own possibly empty V_i, and keep the first and last blocks unmarked. Use the same x,t orientation with V the union of all V_i. For each i define e_i=E_i/h_i, j_i=J_i/h_i, l_i=L_i/h_i as above, and t_i and a_i as the fractions of Hamiltonian paths in D_i ending and starting in V_i. Keep p_i,s_i,r_i as above.

The exact count of this strong tournament is

`H=h(U)+q(U)+h(U) Σ_i [e_i B_i − j_i − l_i s_i − t_i a_{i+1}]`,

where a_{m+1}=0, B_1=1 and `B_{i+1}=r_i B_i−l_i`.

Proof: with predecessor x absent or equal to the sentinel, the preceding proof applies, except that the forbidden starts in V are summed over the disjoint choices of their block. The provisional count of paths with predecessor in V_i puts no restriction on the first vertex after x and contributes h(U)p_i e_i. Then remove paths whose next vertex lies in V_j. For j=i the correction is h(U)j_i. For j=i+1 it is h(U)t_i a_{i+1}: all earlier tape vertices go in the first slot, all later ones in the second. For j>i+1 it is zero, since a block strictly between them can enter neither slot. For j<i it is h(U)e_i l_j times the product of r's strictly between j and i: the left exterior belongs to the first slot, the right exterior to the second, and the intervening blocks can occupy both. In particular these forbidden cases are disjoint, without a hidden triple term. The sum of the last kind is h(U)e_i(p_i−B_i), by expanding the displayed recurrence. The last unmarked block ensures a vertex after x whenever the predecessor is in V.

Notice the digit identity `B_{i+1}−1=r_i(B_i−1)+(r_i−1−l_i)`. All digits on the right are nonnegative because L_i≤q_i−h_i. They are not necessarily integers or independently free. But for an actual C3 block, h=3, q=12, and L_V=3|V|: fixing the first vertex of slot 2, its lengths 1, 2 and 3 each give exactly one pair. Thus the four normalized digits 0,1,2,3 are all physically possible. Also E_V=3|V|, J_V=|V|(|V|−1), and a=t=|V|/3. For J, each internal arc from u in V to v in V yields exactly two covers with first ending u and second starting v (one cover for each choice of which slot has two vertices); nonarcs yield none.

Take one unmarked singleton at each end and d consecutive C3 blocks between. After the left guard B=2. At the end of the d C3 blocks, any integer in `[4^d+1,2·4^d]` is obtained by reading the base-four digits of B−(4^d+1) and taking the corresponding complement masks. This is literal target-directed writing of one register using strong tournaments, not just residue coverage.

It is also a clean failure mode. Up to rotation of each C3 module there are only 4^d distinct mask programs. At least half of q(U)=4·12^d contributes paths having t immediately before x: simply place the left singleton in the second slot, where it gives a legal restart. Thus H≥2·12^d and there are O(X^(log_12 4)) distinct counts ≤X across all depths in this exact architecture. A universal address register does not by itself control the full scalar readout. What needs to change is the cost of a physical digit or the projection from registers to H, not the base-four address decoder.

## The local carry controls the first nontrivial 3-adic digit

For the C3 tape with d≥2, writing z_i=|V_i| on the C3 blocks, the same formula gives

`3^(2−d) H ≡ −Σ_{i=1}^{d−1} z_i z_{i+1} (mod 3).`

In fact 9H/h(U) is integral. The terms involving the running digit, L and q(U) are multiples of nine, the same-block correction is a multiple of three, and only the neighboring term with denominator nine can survive modulo three. In particular H is always divisible by 3^(d−2), and its valuation is exactly d−2 if and only if that nearest-neighbor sum is nonzero modulo three. This lower bound is attained at every depth d≥2 by choosing two consecutive singleton masks and all other masks empty. For d=2 that same rule realizes counts not divisible by three. This result is about the actual H, not just the register.

A version does not depend on using C3. At each block put g_i=gcd(h_i,q_i). For any prime p, let w_i=v_p(g_i). The tournament's count has p-valuation at least the sum of w_i minus the two largest w_i (with absent weights taken as zero). In the unnormalized proof, every summand multiplies one h_i or q_i from each unexposed block; only the independently chosen predecessor and successor at x can expose up to two blocks. This shows why increasing head bandwidth without choosing the bank's arithmetic is not enough. A mutually coprime bank eliminates this particular growing divisor, but still needs an absolute readout and exact carry correction.

## A balanced bank swap changes H without changing the address

The contact bound above recovers the lower floor for contiguous triangle chains from the [earlier all-k triangle filtration](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/berggren_target_return_20261111/arithmetic/TRIANGLE_ALL_K_PROOF.md). The nearest-neighbor valuation test is for this particular programmable tape, not a new general extension-content classification.

The bank can be changed without touching the digit output. In addition to C3, take TT2, whose normalized radix is again q/h=4. In TT2, the second-slot start weights at its source and sink are (2,1), and the first-slot end weights are (1,2). Thus both have exactly the same digit l in {0,1,2,3}. At digit 0 both have all marked readings zero. At digit 3 both have (e,j,l,t,a)=(3,2,3,1,1), even though their unnormalized h,q differ by a factor of three. At an isolated nonzero digit (both neighboring masks empty), convert C3 to TT2 at that site and simultaneously convert TT2 to C3 at any digit-0 parking site. Tape order, h(U), q(U), and **every** prefix B and suffix s stay fixed. Here “tape order” means the order of the block positions (the number of vertices is also fixed in aggregate). The effect of the swap is exactly

| kept digit l | change in H, divided by h(U) |
| --- | --- |
| 1 | B_i |
| 2 | 2/3−B_i |
| 3 | 0 |

For digit 1, E/h switches from 1 to 2 and J=0 in both. For digit 2 it switches from 2 to 1 and J/h drops from 2/3 to 0. Digit 3 has no normalized change. Disjoint pairs of active and parking sites, with empty neighbors at every active site, give independent signed shifts on the same address. The formula certifies exact control in one coordinate while keeping a reusable address readout in another. It does not give a unit interval: the B_i separate multiplicatively on base-four banks, and the factor h(U) can retain a large power of three. A concrete next problem is to build a second signed bank in the same radix with flat (ideally ±2) shifts after the height-dependent B term is canceled, not merely to add more digits.

## Flat unit repair: change pair interaction, not either register

Use the regular cyclic five-vertex tournament D, lexicographic forward-arc code 76 (rows 00011, 10001, 11000, 01100, 00110). Compare masks {0,2} and {0,1}. Both have (h,q,E,L,T,S)=(15,120,42,42,6,6), where T,S count the Hamiltonian paths ending/starting in the mask. The only differing coordinate is J, which is 8 for {0,2} and 10 for {0,1}. Thus for any nonempty transitive left and right guards of lengths a,b,

`H(W(V))=15+120·2^(a+b)+42·2^a−42·2^b−J(V).`

The two physical strong tournaments differ by exactly two at every height, have the same entire programmed digit register, and are joined by reversing the two incident arcs to 1 and 2. At a=b=1 the counts are 485 and 487. In fact the flat shift persists with arbitrarily many other live masks: in the multihead formula only j_i changes, so the unnormalized gap is twice the product of the h_k for all other blocks. This conditional extension states exactly when the unit survives neighboring live blocks: their H product must be 1.

Why the bank works: in any vertex-transitive D, the four coordinates E,L,T,S are proportional to mask size. At fixed mask size, all freedom lies in the symmetric pair energy `Σ_{u<v in V}(J_uv+J_vu)`, where J_uv counts two-slot covers whose first ends at u and second begins at v. The energy coefficients are constant on unordered-pair orbits. The five-cycle is a two-point spatial correlation changing the count while all one-point statistics stay fixed. On arbitrary ordinal tapes it remains embedded within the exact two-register formula; it is not a construction for choosing the baseline from an arbitrary target.

Sensitively test the mechanism rather than high-count coverage: for the standard cyclic tournaments of orders 5, 7, 9, 11 and 13, compute the attainable pair energies for each fixed mask size. The gcd of differences (excluding sizes yielding only one value) is respectively 2, 2, 2, 6 and 2 at each of those sizes; the longest run of consecutive step-two energies over all sizes has lengths 2, 2, 1, 1 and 1. Thus even when the residue obstruction is absent, just enlarging this symmetric core loses the exact two-step repair in three cases; it does not grow a window. This is a finite architecture test, not a conjecture about all cyclic tournaments. `check_energy.py` reproduces the table and verifies the order-five flat switch at four different left/right heights independently.
