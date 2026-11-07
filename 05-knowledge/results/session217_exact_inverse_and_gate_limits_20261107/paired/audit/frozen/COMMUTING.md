# Q/Z is an abelian two-counter construction, not a binary digit language

Status: PROVED elementary result; independent audit pending. The unrestricted odd Hamiltonian-path spectrum, existence of a full paired Q/Z source, and eventual coverage remain OPEN. This result changes the amplification interpretation of a hypothetical pair, not the validity of the inherited pair criterion.

## Inheritance and scope

Input: the independently audited seven-state transfers in `session216-structural-release/LEFT_RESET.md` and its audit/REVIEW.md. Source order is nontrivial, deletion is nonempty, and all graphs are physical tournaments. The vector is (h,s,e,L,R,d,b); h and d count whole and root-deleted Hamilton paths, s/e are root-start/end counts, L/R are the inherited rooted two-cover coordinates, and b is the deletion two-cover count. Q adds a vertex dominated by the root and dominating all old nonroots. Z adds a vertex dominating the root and dominated by all old nonroots. Both keep the same root.

The nearest positive mechanism is the unequal-deletion Q ray. The hostile near miss is treating two named moves as two independent digital choices. The least-used sidecar here is the isomorphism type of the entire rooted graph, not just its nine-state projection. Live concepts: chronological insertion, ordered root amalgamation, reset transients, dyadic counters, interval occupancy, and normalized menu defects. No new negative bounded search is used in the theorem.

## 1. Graph-level commutation, including transients

For EVERY initial rooted tournament G, and every word w in {Q,Z}, the resulting rooted graph is determined up to root-preserving isomorphism by q=#Q(w), z=#Z(w).

Proof. Delete the root. The added Q vertices form a transitive prefix before G-root, with the most recently added Q first. The added Z vertices form a transitive suffix after G-root, with the earliest added Z first. Every Q vertex dominates every Z vertex regardless of which was inserted first. The root dominates all Q vertices and loses to all Z vertices. The old rooted incidence is unchanged. These statements specify the resulting graph completely as a labeled block construction with canonical order inside each new block. Thus chronology is erased even before taking any counting interface. In particular QZ and ZQ are isomorphic for arbitrary starting states, including those not satisfying s=e=d.

The exact seven-state transfers are

Q(h,s,e,L,R,d,b)=(h+L-s,d,e,2L-d,R+b,d,2b),
Z(h,s,e,L,R,d,b)=(h+R-e,s,d,L+b,2R-d,d,2b).

Their common composite is

(h+L+R+b-s-e,d,d,2L+2b-d,2R+2b-d,d,4b).

Thus one Q and one Z, in either order, reset to s=e=d, and all remaining moves preserve this condition. At a fixed word length k there are at most k+1 rooted isomorphism types, not 2^k. Some types can coincide further.

## 2. Exact stable closed forms

For a physical stable source s=e=d, put a=L-d and rho=R-d. For all q,z>=0, write u=2^q and v=2^z. Then

b(q,z)=uv b,
a(q,z)=u[a+(v-1)b],
rho(q,z)=v[rho+(u-1)b],
h(q,z)=h+(u-1)a+(v-1)rho+(u-1)(v-1)b.

Proof: apply q Q moves first, then z Z moves, using the geometric-series sums. Section 1 permits this reordering for every word. These formulas include q=0 or z=0.

For the width-eight compiler, suppose d=r+3j with r in {1,2}, 0<=j<=8. The adjusted base X=h+a+2j becomes

X(q,z)=2^(q+z+1)b+2^(q+1)(a-b)+2^z(rho-b)+C,
C=h-a-rho+b+2j=3d+2 kappa+2j.

The second equality for C uses the inherited conserved defect. Coefficients a-b and rho-b can be negative; do NOT infer positivity from this expanded expression. The unexpanded h formula supplies the required positivity. Compiled counts are 3X(q,z)+2r.

## 3. All-length counting bound

For every nontrivial physical stable source, a>=d>=1 and rho>=d>=1. Indeed L contains the d covers with singleton root and the d full paths ending at the root with the ordinary bin empty; these families are disjoint. Similarly R>=2d using full root-start paths. Also h>=d, since the d root-start paths are among all Hamilton paths.

Consequently h(q,z)>=d 2^q and h(q,z)>=d 2^z. Since X(q,z)=h(q,z)+a(q,z)+2j>=h(q,z), every compiled count H<=N satisfies

q,z <= floor(log2(N/(3d))).

If N<3d there are no such outputs. Otherwise the number of distinct compiled counts at most N is bounded above by

(1+floor(log2(N/(3d))))^2.

This is an upper bound on output cardinality, not a statement that each pair of exponents gives a distinct count.

For an arbitrary fixed starting source, words containing both letters reduce to the single stable source QZ(G), followed by q-1 and z-1 additions. The two remaining axes are affine dyadic sequences. More explicitly, if q>=1,

h_q=h+L-s+(2^q-2)(L-d),
L_q=d+2^q(L-d), s_q=d.

Hence its adjusted base is an affine function of 2^q with nonnegative coefficient 2(L-d). The coefficient may vanish, yielding only one distinct value. Here L>=d because the singleton-root cover is present. For z>=1,

h_z=h+R-e+(2^z-2)(R-d),
L_z=L+(2^z-1)b, s_z=s.

The adjusted base is affine in 2^z with coefficient R-d+b>0. The old s (on the Z axis) or d (after Q) must satisfy the actual compiler residue/budget condition; inadmissible axes supply no compiled outputs. Each admissible axis contributes O_G(log N) distinct outputs, or one constant output. Therefore a fixed finite source bank, allowing arbitrary Q/Z words and their valid compiler outputs, supplies only O_bank((log N)^2) distinct counts <=N.

In particular its natural density is zero and it cannot contain an eventual interval of odd integers. A full paired five-state seed would still be a useful synchronized TWO-COUNTER family, but would not by itself be a binary interval-amplification certificate. An infinite bank, changing the root, other moves, or rebuilding sources is outside this theorem.

## 4. Executable controls

Run `python3 verify_commuting.py` and `python3 -O verify_commuting.py` from this directory. No external inputs or dependencies are needed. The checker independently implements endpoint subset DP. Its universe is every labeled tournament of orders 2..4 with root 0. It checks literal Q/Z transfers on arbitrary transient sources, every length-five word against its sorted word, closed forms after stable reset, and literal stable descendants with at most two further letters. The finite checks support but do not replace the all-length proof. Results: 3,552 algebra/word controls and 592 literal DP controls, normal and optimized outputs identical.
