# Physical cut intervals, boundary-balanced sumsets, and their cofactor obstruction

Session 192 continuation, 2026-11-05. **PROVED elementary statements; FINITE-EXACT certificates below; independently audited.** The global Hamiltonian-path spectrum conjecture remains OPEN. This packet constructs actual tournament families, rather than interpreting convex or toric support as a construction.

## 1. Inheritance and the question being changed

The inherited [weighted-cut proof](../weighted_cuts/PROOF.md) proves, for every tournament U and root attachment S,

    H(U+x_S) = H(U) + Delta_U(S),
    Delta_U(S) = sum_(i in S,j outside S,i->j) c_ij,
    c_ij = Q_ij+Q_ji > 0, with c_ij even.                 (1)

Here x beats exactly S. The H=13 core has strong-child support

    {23,25,29,31,33,37,43,45}.                            (2)

Its normalized increment support has gcd one, but holes. Its three-fold abstract sumset is an interval, but that addition is not a tournament operation. We retain the actual common parent and specify when simultaneous switches really add.

The source inspiration is THM-2552 / `01-canon/theorems/THM-2552-flat-common-base-gain-versus-unjoinable-wall-spectrum.md`: individually available gains cannot be added without one common physical carrier. The useful difference here is that a cut response is quadratic, so after a physical Boolean cube is supplied, pairwise interaction tests do suffice. No LRC conclusion or arbitrary higher-order independence is transported.

The second inspiration is THM-3541 / `01-canon/theorems/THM-3541-keller-grade-monoid-factorization-lengths-and-elasticity.md`: numerical factorization intervals must be distinguished from realizations on one object. Below, the ordered-join construction supplies the missing realization and also forces a previously invisible cofactor scale. These are methodological sources, not dependencies for the elementary proofs.

## 2. A constructive integer-interval criterion on a physical cut cube

Put w_ij=c_ij/2, and let W be the symmetric matrix with these off-diagonal entries and zero diagonal. For a cut indicator z in {0,1}^n,

    f(z)=Delta_U(z)/2
        =sum_i a_i z_i - sum_(i<j) w_ij z_i z_j,
    a_i=sum_(i->j) w_ij.                                 (3)

Choose disjoint coordinate sets and, on each, two binary patterns. Equivalently choose z^0 and disjointly supported vectors d_1,...,d_r with entries in {-1,0,1}, such that

    z(epsilon)=z^0+sum_j epsilon_j d_j                    (4)

is Boolean for every epsilon in {0,1}^r. Thus all choices in (4) are actual attachments to this same U. Suppose

    d_i^T W d_j=0 for every i != j.                      (5)

Let b_j=f(z^0+d_j)-f(z^0). Substitution in (3) gives the exact identity

    H(U+x_(z(epsilon)))=H(U)+2f(z^0)
                               +2 sum_j epsilon_j b_j.  (6)

There are no untested higher interaction terms: (3) has degree two.

Discard zero b_j and sort their absolute values as q_1<=...<=q_s. If

    q_1=1,  q_j<=1+sum_(i<j)q_i for j>=2,               (7)

then this one cut cube realizes every odd integer from

    H_min=H(U)+2f(z^0)+2sum_(b_j<0)b_j

through H_min+2sum_j|b_j|. For a strong U, these children are strong provided all cuts are nonempty and proper. For a weak U, use the exact source/sink SCC test in the inherited packet.

Proof of the last assertion: reversing a switch with negative b_j converts the values to H_min plus twice the subset sums of the q_j. If the previous coins fill [0,M], a coin q<=M+1 joins [0,M] to [q,q+M] without a gap. Conversely, at the first violation of (7), the number M+1 cannot be made. This is the standard complete-sequence subset-sum criterion, applied only after (4)-(5) certify its physical meaning.

For a pair exchange d=e_b-e_a, condition (5) with e_d-e_c is simply

    c_ac-c_ad-c_bc+c_bd=0.                              (8)

Ordinary independent root bits cannot produce an additive cube of dimension two: for d_i=e_i,d_j=e_j, the mixed coefficient is -w_ij<0 by (1). Exchange or more complicated balanced switches are genuinely required; dropping these interactions is not justified by capacity positivity.

### An explicit four-odd interval

Take the strong six-vertex parent with outgoing bit rows

    U=(50,8,3,5,14,30), H(U)=25.

Start at attachment mask 37={0,2,5}, and independently exchange 0 with 4 and 2 with 3. The relevant capacity rectangle is

    c_02-c_03-c_42+c_43 = 14-10-18+14 = 0.

The two changes in H are +4 and -2. The four masks and actual counts are

    37 -> 131,  52 -> 135,  41 -> 129,  56 -> 133.       (9)

Consequently this one strong parent realizes the entire odd interval [129,135]. This is an explicit certificate, not a claim that the interval itself was previously unknown. The verifier counts the parent and all four children independently by literal permutations, as well as subset DP and capacities.

The H13 hostile (2) does not meet a universal interval conclusion. In particular gcd one, even positive capacities, and a correct minimum/maximum are insufficient. Criterion (7) is sufficient for a supplied cube, not necessary for every interval that might occur elsewhere in the cut support.

## 3. Ordered-join capacities have an exact four-coordinate boundary

Let U=V_1 -> ... -> V_m be an ordered join of nonempty tournaments. Write h_a=H(V_a), P=product_a h_a. For a vertex i of a block V, let s_i,e_i count Hamiltonian paths starting and ending at i. Define

    B(V)=sum_(A subset V) H(V[A])H(V[V\A]),
    R_i(V)=sum_(A subset V,i in A) s_i(V[A])H(V[V\A]),
    L_i(V)=sum_(A subset V,i in A) e_i(V[A])H(V[V\A]).  (10)

Empty paths in the complementary factor are allowed, with H(empty)=1. Thus R_i is an ordered spanning two-path response whose second path starts at i and is nonempty; L_i is its ending analogue. These are physical two-path counts, not formal scalar states. Sum_i R_i=sum_i L_i=B(V)-H(V).

Internal capacities are

    c_ij(U)=(P/h_a)c_ij(V_a),  i,j in V_a.               (11)

For i in V_a and j in V_b with a<b, put

    O_ab=product_(t<a or t>b)h_t,
    M_ab=product_(a<t<b)B(V_t).

Then

    c_ij(U)=O_ab [ R_i(V_a)L_j(V_b)M_ab
                   +1_(b=a+1)e_i(V_a)s_j(V_b) ].       (12)

Proof. The valid exposed step i->j can occur in a Hamiltonian path only for adjacent blocks. It gives the second term. For the reverse exposed step j,i, all other steps are valid. Each block strictly between a and b splits into its portion before and after the exposed step, giving B(V_t). The initial block contributes a nonempty second path starting at i, giving R_i; the final block contributes a nonempty first path ending at j, giving L_j. Blocks outside [a,b] each occur as one path. The ordered component direction makes these decompositions unique and their converses valid. The same argument inside one block proves (11). This works even when a named block is not itself strong.

### Physical sumset theorem

For every block choose a nonempty family F_a of attachment subsets. Require each of the following four sums to be constant on F_a:

    R(S)=sum_(i in S)R_i,   L(S)=sum_(i in S)L_i,
    e(S)=sum_(i in S)e_i,   s(S)=sum_(i in S)s_i.         (13)

Let D_a={Delta_(V_a)(S)/2:S in F_a}. Formula (12) makes the total cross-block cut contribution K constant on the full Cartesian product of these families. Explicitly its a,b term is

    O_ab [ R(S_a)(L(V_b)-L(S_b))M_ab
             +1_(b=a+1)e(S_a)(s(V_b)-s(S_b)) ].         (14)

Every independently chosen local subset is a real subset of U. Thus the EXACT support, not just an inclusion, is

    {H(U+x_S): S intersect V_a in F_a for every a}
       = P+K+2 sum_a (P/h_a)D_a.                       (15)

For a vertex-transitive block on n_a vertices, s_i=e_i=h_a/n_a and R_i=L_i=(B_a-h_a)/n_a. Therefore fixed cardinality |S|=k_a automatically gives (13). If these blocks are strong, k_1>0 and k_m<n_m ensure every resulting child is strong.

This gives a broad physical setting in which sumset arguments are legitimate. It also identifies exactly what must be retained when vertex transitivity is unavailable: the four boundary sums (13), rather than H alone.

## 4. The cofactor lattice is an all-order obstruction to this slice method

Every difference between two values in (15) is divisible by

    2g,  g=gcd_a(P/h_a)=P/lcm(h_1,...,h_m).              (16)

The identity follows prime by prime: the minimum valuation among P/h_a is the total valuation in P minus the largest valuation in one h_a. Therefore a single slice of the form (15) cannot contain consecutive odd values when g>1, regardless of how broad its local supports or abstract sumsets are.

If the nontrivial blocks come from a fixed finite library, their h values have bounded lcm. As the number of blocks with h>1 grows, P grows, so g grows. Adding any number of h=1 singleton blocks does not remove this factor. Hence repeating a finite independent block library cannot amplify one boundary-balanced slice into longer and longer consecutive-odd intervals.

This is explicitly a SLICE obstruction. Different boundary slices can have different constants K, and their union may escape it. It does not exclude arbitrary root cuts, unions of parents, all tournament constructions, or the global spectrum conjecture. It says which retained scale prevents the tempting independent-sumset proof.

## 5. Concrete all-order families

Let C5 be the regular cyclic tournament i->i+1,i+2 mod 5, with rows

    (6,12,24,17,3).

Direct enumeration gives

    h=15, B=120, s_i=e_i=3, R_i=L_i=21;
    {Delta(S)/2:|S|=2}={13,14}.                         (17)

For example S={0,1} gives 13, while S={0,2} gives 14. These are small finite input constants verified by exact path enumeration; all-order formulas below follow from (11)-(15).

### 5.1 Repeated C5 blocks: actual sumsets, wrong interval scale

Let U_m be the ordered join of m copies of C5, and require |S|=2 in each block. Every child is strong. The exact support is

    { A_m+2*15^(m-1)t : t=0,...,m },                    (18)

where A_1=41 and, for m>=2,

    A_m=15^(m-2)(54*8^m+66m+117).                       (19)

Every t is constructed by using {0,2} in t blocks and {0,1} in the rest. No semigroup-normalization step is being used.

To derive the constant, cross capacities at block separation d are

    k_(m,1)=450*15^(m-2),
    k_(m,d)=441*120^(d-1)*15^(m-d-1), d>=2.

Each block pair contributes 2*3 times this capacity to Delta. Thus

    A_m=15^m+26m15^(m-1)+6sum_(d=1)^(m-1)(m-d)k_(m,d).

The sum simplifies to (19), using
sum_(d=1)^(m-1)(m-d)8^(d-1)=(8^m-7m-1)/49. Initial supports are

    m=1: {41,43};
    m=2: {3705,3735,3765};
    m=3: {419445,419895,420345,420795}.

Their difference gcd is exactly 2*15^(m-1). All m>=2 values are composite: (19) and (18) have a factor 3 (and a factor 3*15^(m-2) for m>=2). The construction is real, but repeated independent modules do not improve prime or odd-interval coverage.

### 5.2 A general sink-padding law and a dyadic four-odd window

Take ANY core V, let U=V->T_t where T_t is a transitive chain, t>=1, and let the new root lose to every chain vertex. For a fixed attachment S within V, formulas (11)-(12), with singleton B=2, give

    H(U+x_S)=H(V+x_S)+R(S)(2^t-1)+e(S).                (20)

Indeed the j-th chain vertex receives capacity R_i*2^(j-1) from old vertex i, with the extra valid-adjacency contribution e_i only at j=1. Summing this geometric series proves (20). For strong V and nonempty S, every resulting child is strong. Two local cuts with equal R(S) and e(S) therefore retain their exact spacing at every t. Full vertex transitivity and all four signatures are unnecessary for this one-sided transport.

The inherited S191 construction packets already use physical two-path responses and exponential transitive-wrapper laws; no novelty is claimed for the broad idea of an affine dyadic wrapper. The addition here is the explicit cut-capacity identity (20), its boundary-fibre interpretation, and the following concrete common-slope bank. Unlike the single-expanding-run one-apex architecture in [Session 191 one-apex proof](../../session191_tournament_growth_20261104/one_apex/PROOF.md) section 6.1, these cores are not required to be one-apex.

Use two strong five-vertex cores:

    V9 =(20,1,2,7,14), H(V9)=9;
    V13=(18,8,3,5,14), H(V13)=13.

The four exact boundary records are

    core  mask  H(core+x)  R(S)  e(S)  H(core+x)-R(S)+e(S)
    V9     10       29       32     4             1
    V9     17       33       32     2             3
    V13     6       31       32     6             5
    V13     3       33       32     6             7.       (21)

The last column is the intercept in (20). Consequently, for EVERY k>=6, all four odd values

    2^k+1, 2^k+3, 2^k+5, 2^k+7                          (22)

are realized by strong tournaments of the SAME order k+1. Construct them by appending t=k-5 singleton sink vertices to the named parent, then adding x with the exact mask in (21) on the core and no outneighbors in the padding. All old core arcs remain unchanged. The H13 core in (21) is isomorphic to the hostile core in (2); its holes do not prevent this exact transported two-point fibre.

For example k=6 gives 65,67,69,71 on seven vertices, while k=7 gives 129,131,133,135 on eight vertices. The separate cube (9) realizes the latter four values on seven vertices with one common parent. The constructions prove attainment, not minimum realizing order, new global membership, unbounded interval length, or recursive closure.

As another instance of the same law, a vertex-transitive V of order n, with parameters h,B and |S|=k, has common shift

    k[((B-h)/n)(2^t-1)+h/n].                             (23)

Thus C5 and k=2 give strong tournaments of order t+6 with

    H=42*2^t+5  and  H=42*2^t+7,  t>=1.                  (24)

Their first pairs are (89,91), (173,175), (341,343), (677,679); the unpadded pair is separately (41,43). These are exact consecutive odd pairs, not a claim of primality. Equations (20)-(24) explain why h=1 padding transports local spacing while nontrivial independent blocks dilute it. Neither construction yields intervals of unbounded length.

## 6. Verification, scope, and next test

`verify.py` is standalone and uses explicit exceptions rather than removable assertions. Run

    python3 verify.py > receipt.json
    python3 -O verify.py > receipt.optimized.json
    cmp receipt.json receipt.optimized.json

It checks the explicit interval cube by literal permutations and subset DP; the H13 hole control; (11)-(12) on every ordered pair of labelled blocks of orders 1..3; all eight triples from {singleton,C3}; the C5 interfaces; all 10^m local size-two cuts for m<=3; direct DP for each displayed progression value through parent order 15; the padded C5 pairs through t=8; all four dyadic-window core signatures by literal permutations; and all four padded children for k=6..12. Finite data do not assert a general parent-bank interval or solve the global conjecture.

The proof now suggests a precise next target: find large boundary-balanced local cut intervals with cofactor one, or arrange several different boundary slices whose constants K cover the cofactor residues. The first task is a genuine constrained support problem; the second must retain (14), so arbitrary free sums are again unavailable. A successful certificate must output the actual masks, boundary sums, cross constant, and scale, rather than only a normalized integer interval.
