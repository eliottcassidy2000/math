# Exact named endpoint transfer and the finite-suffix-bank obstruction

2026-11-04. PROVED algebraic transfer and counting obstruction; FINITE-EXACT independent controls. No claim of closed ternary grammar or global tournament spectrum completeness.

## 1. An augmented interface which genuinely retains joint ports

Let T be a tournament with physically named vertex set V and adjacency matrix A, with A_ij=1 when i→j. Work over the commutative integer algebra

R_V = Z[z_v : v in V] / (z_v² : v in V).

Set Z=diag(z_v) and G=(I−ZA)^−1 Z. The inverse is the finite geometric series: every positive-degree monomial of degree greater than |V| vanishes. Then G_ij is exactly the sum of weights of all directed simple paths from i to j, including the one-vertex path. A path has weight the product of its vertex variables. Repeated vertices vanish; every surviving walk is a simple path, counted once.

Write [V] for coefficient extraction of the full squarefree monomial, P=Σ_ij G_ij, e_r=Σ_i G_ir and s_r=Σ_j G_rj. The following are exact, with all names retained:

- H(T)=[V]P.
- The number of Hamilton paths ending/starting at r is [V]e_r or [V]s_r.
- E_r(T,2)=[V] e_r(1+P), for ordered two-slot covers whose first slot ends at r and whose second slot may be empty.
- S_r(T,2)=[V] s_r(1+P), using the analogous start condition.
- P_T(2)=[V](1+P)² for nonempty V.
- For example, the joint response with first slot ending r and second starting s is [V]e_r s_s.

More generally products of specified entries, rows, or columns count ordered disjoint slots with prescribed joint endpoints; the nilpotent variables enforce disjointness without losing correlations. Constraints on internal named vertices can be retained by coefficient restriction. This is an exact sufficient interface, not a bounded-size compression: its worst-case size grows exponentially with |V|.

### Birth transfer

Add x, let a_i=[i→x] and b_j=[x→j], so a_i+b_i=1. In R_{V∪{x}},

G'_{V,V}=G+z_x(Ga)(bᵀG),
G'_{V,x}=z_xGa,
G'_{x,V}=z_xbᵀG,
G'_{x,x}=z_x.

Proof: a simple path either avoids x, ends at x, starts at x, or visits x internally. In the last case its unique split at x gives a left path ending at an in-neighbor and a right path starting at an out-neighbor. The squarefree product kills precisely overlapping left/right paths. These cases are disjoint and exhaustive. In particular every new named joint-port field can be extracted after the update; agreement of separate marginal fields is unnecessary.

### One edge reversal

Suppose i→j is reversed, with every other edge unchanged. Then

G'_{ab}=G_{ab}+G_{aj}G_{ib}−G_{ai}G_{jb}.

A new path which uses j→i splits uniquely into a j-ending and i-starting old path. Those pieces cannot themselves use i→j, since they would repeat i or j. Conversely the removed paths split at i→j. Products with overlapping vertices vanish. All remaining paths are unchanged. Equivalently this is the resolvent identity; higher correction terms vanish because each would repeat i or j.

For multiple reversals use the formula successively with the UPDATED matrix. Applying all first-order corrections to the original matrix drops genuine joint terms. Vertex deletion is the substitution z_v=0 followed by removing its row/column; physical re-marking selects a different existing row/column and is not a state reset. This provides exact transfer for the elementary operations, but it does not prove the legality of any native Hamilton-path-complement guard; such guards require their separately supplied literal witnesses.

## 2. A finite fixed-suffix bank cannot close both ternary branches

For a fixed strong marked suffix (B,w), the established zero-ready grammar with transitive prefix TT_s has

H=h+2^sF,  s≥1,

where h=H(B)≥1 and F=E_w(B,2)≥1. Let a bank consist of k fixed such suffixes. Allow arbitrary clock resynchronizations, arbitrary switching among those k suffixes, arbitrary selection among duplicate representations, and arbitrary finite phase labels. None of these enlarge its set of numerical counts beyond the union of the k clock families.

For X≥2 the bank contains at most k floor(log_2 X) distinct counts at most X: each admissible 2^s≤H≤X, so there are at most floor(log_2 X) possible s per suffix. Duplicate representations only decrease the number.

Suppose this bank were closed under both H→3H+2 and H→3H+4 starting from H_0. At depth n it must contain 2^n distinct integers, all at most

X_n=3^n(H_0+2)−2.

Indeed a branch word d_1,...,d_n∈{2,4} produces
3^nH_0+Σ_{j=1}^n d_j3^{n−j}.
After subtracting the common all-2 sum, the residual has unique base-3 digits in {0,2}, proving distinctness. Therefore

2^n ≤ k floor(log_2 X_n)

would be necessary for every n. The right side is O(n) and the left side exponential: contradiction.

This is stronger than failure to synchronize two particular clocks. NO finite bank of fixed marked suffixes, however chosen and however richly phase-labelled, can represent a fully branching reusable ternary grammar solely by changing transitive-prefix lengths. It is not an obstruction to an infinite family of genuinely changing suffixes, a different prefix grammar, or a single selected branch.

Quantitative control: for H_0=56935 and k=2, depth6 requires64 distinct counts at most41507071, but any two fixed clock families offer at most50. This bound does not require the particular suffix values or phase coordinates.

## 3. Literal independent controls and phase scope

Source: pinned commit 4ad839b20dd280cdcfbadf86cb62b2be3c9615ad, directory `05-knowledge/results/triangle_return_170807_20261103/`, specifically RETURN.md, ZERO_READY.md, audit.json and zero_ready.json. Frozen input JSON is included here. The verification code is independently written standard-library Python, using backward endpoint recursion and direct path enumeration for the small algebra tests.

Verified literal H sequence:56935,68533,67599,170807.

Verified canonical marked suffix fields:

- Old suffix: (h,F)=(18977,18979); first slot ending at physical39 and second starting at physical6 gives joint response229.
- New suffix: (h,F)=(2195,42153); first slot ending at physical31 and second starting at physical6 gives joint response371.

These two newly calculated joint quantities concern the SUFFIX tournaments, not an unlabelled identification of the two full states.

The inherited81/27 hostile pair is recomputed from its actual two adjacency matrices: both have H=462915, while the same actual cut256 gives1061775 versus1395279. A separate exact backtracking computation now also enumerates ALL adjacency-preserving permutations, independently establishing full automorphism orders81 and27 (361 and121 search nodes, respectively). Every candidate image is constrained by outdegree and all already mapped adjacencies; the search branches over every eligible unused image, so no automorphism is omitted. See automorphism.py and automorphism_receipt.json. This control demonstrates that equality of H is insufficient for literal birth transfer; the symmetry orders are not substituted for joint incidence data.

The inherited F49/charge and F343 phase data remain frozen in input_audit.json. In particular the final F49 pair is fixed while the F343 codec has its three-step orbit. The algebra above takes the physical matrix and names as state; it does not identify those distinct observer actions or infer tournament automorphisms from them. Further candidate phases must be recomputed from the candidate's marked endpoint fields, not copied from an ancestor.

Exhaustive algebra controls: every labelled tournament of orders1 through4, every birth cut, and every edge reversal. Compared EVERY coefficient of EVERY endpoint matrix entry against fresh enumeration:1098 birth identities and410 flip identities. These are finite controls of the preceding proofs, not the basis for an all-order extrapolation.

Run: `python verify.py`. See receipt.json and verify.out. No third-party software is required.

## 4. Independent replay of the next concrete return

The separately constructed19-vertex candidate has H=512423=3·170807+2. Input names, full rows, flips and birth names are in the sibling returns/candidate.json. Independent C++ subset DP in joint.cpp verifies its count and the following ACTUAL physical endpoint fields, in the order [E1,S1,H(T−r),E2,S2]:

- root39: [45,5655,5655,512423,20042207].
- spare6: [21382,270,173211,33915934,1588204].
- port51: [90,126692,259039,782772,94607628].
- inport4: [45,477,477,1013059,2455969].
- old root50: [45,253384,259039,518078,144393953].

The full20 directed joint end/start two-slot responses are in candidate.out. The relation retained is J(r,s)=[V]e_r s_s, not a pair of separately tabulated marginals. In particular J(39,6)=351 and J(6,39)=1493231. Root39 has E2=H, verifying zero-readiness directly. Its golden vector is [512423,10021081,10021126]; root50 instead gives [512423,72196954,72196999]. Phase conversion is done by the separate phase packet using these recomputed fields.

This independent numerical check does not establish native guard compilation or an iterable family; use the separate return packet's exact scope and certificates for those claims.

The new suffix clock is5655+31673·2^t. Even this new bridge has just its observed+2 clock intersection with the170807 suffix clock, and no+4 intersection. Indeed the equation

5655+31673·2^t =3(2195+42153·2^s)+d

is equivalent to31673·2^t−126459·2^s=930+d. A positive right side forces t>s. Both coefficients are odd, so the left side has2-adic valuation s. For d=2 this gives s=2 and then2^t=16, hence exactly(s,t)=(2,4). For d=4 it gives s=1 and then2^t=253852/31673, which is not even an integer. Thus the new literal return is additional finite progress, not closure by clock resynchronization.

Reproduce the candidate control with the existing system compiler (no software installation):

`g++ -O2 -std=c++17 joint.cpp -o joint`

`./joint < candidate.in > candidate.out`

The arithmetic is uint64_t and the actual order19 control is safe: the total number of ordered two-slot path covers is at most20·19!=2432902008176640000<2^64. All nonnegative partial sums are bounded by that number. Do not extrapolate this implementation to larger orders without changing its integer type or proving a new bound. The Python matrix and automorphism checks use arbitrary-precision integers.
