# Signed controller frustration, exact star-phase optimization, and half-path extrema

Status: proved structural identities plus a finite exhaustive certificate through five pairs. No general maximum formula is claimed. In particular, the initial proposed 'all labels coherent' domination theorem is false. The unsigned relaxation remains valid.

## 1. Setup and the directional sign that must not be dropped

Let sigma be an anti-involution of a tournament with m two-element orbits and either zero or one fixed vertex. Choose names i+ and i- so the final half-path fold is permitted precisely for i+. In even order this means i+ -> i-; in odd order it means i+ -> c -> i-, where c is the fixed vertex.

Between distinct pairs i,j the four arcs have matrix

    [[x,y],[1-y,1-x]],   x,y in {0,1}.

Exactly one pair controls the allowed sign: the matrix is constant along rows or along columns. Define the **controller tournament** Q by v -> w when pair v controls their interaction. For each arc v -> w, define b_vw in {+1,-1} to be the sign required of v when w occurs immediately BEFORE v in the half-ordering. If v occurs before w, the required sign is -b_vw. This reversal of the pin is essential.

Every Q and every assignment b defines a valid anti-invariant pair-interaction system. Thus there is no additional realizability condition at this level. The missing condition in the initial relaxation was compatibility of pins within a half-ordering, not existence of the tournament.

## 2. Exact wedge-factor theorem

For a permutation p=(p_1,...,p_m), let D(p) consist of internal positions j for which p_j controls both neighbors:

    p_j -> p_(j-1),  p_j -> p_(j+1) in Q.

Let L(p) be 1 if p_m controls p_(m-1), and 0 otherwise. Put c(p)=|D(p)|+L(p). For m=1 use D empty and L=0.

The exact number of signed half-paths with pair order p is

    W_b(p) = product_(j in D(p)) (1-b_(p_j,p_(j-1))*b_(p_j,p_(j+1)))
             * (1+b_(p_m,p_(m-1)))^L(p).                         (1)

Consequently the fixed Hamiltonian path count is

    F(Q,b) = sum_p W_b(p).                                      (2)

Proof. Each of the m-1 adjacency constraints pins its controller, and the final fold pins the final pair. There are m pin occurrences. The only repeated pins are the two pins at an internal double-controller or the adjacency pin plus fold at the final controller. Each repeated pin subtracts one distinct pinned pair; hence exactly c(p) pairs are free. Internal repeated pins agree iff the two b labels are opposite. The final repeated pin agrees iff b=+1. If all agree the count is 2^c(p), otherwise zero. This is exactly (1).

This is a literal signed-interaction model: each outgoing star wants its adjacent selected edge labels opposite, while the fold supplies a positive boundary field. 'Antiferromagnetic wedge with fold boundary' describes the actual factors, not merely an analogy.

The polynomial (2) is multilinear of degree at most m-1. A path edge belongs to only one controller and appears in at most one factor. No constant term is generated from nonempty products, so its constant coefficient is m!.

## 3. Consequences: exact means, a Hamiltonian-count subfamily, and parity

For EVERY fixed Q, independent uniform b labels give

    E_b F(Q,b) = m!.                                            (3)

Each factor in (1) has mean1, and its edge variables are disjoint from those in every other factor of that permutation.

The terms c(p)=0 are exactly directed Hamiltonian paths of Q. Indeed, along p a double-controller is a backward-then-forward transition, and the fold detects a final backward edge. Avoiding both forces every edge forward. Therefore

    F(Q,b) >= H(Q),   F(Q,b) = H(Q) mod 2.                       (4)

Setting every b=-1 kills all terms with c(p)>0 and yields

    F(Q,-1) = H(Q).                                              (5)

Thus every ordinary tournament Hamiltonian-path count has an exact realization as an anti-involution-fixed Hamiltonian-path count on both 2m and 2m+1 vertices. The converse does NOT hold: fixed counts7 and21 occur in the census, so these are not subject to the ordinary Hamiltonian-count forbidden-set claim.

Double-controller positions and the final-controller position cannot be adjacent. They lie in positions2,...,m. Hence c(p)<=floor(m/2), giving the proved bound

    F(Q,b) <= U(Q) <= 2^floor(m/2) m!,                           (6)

where U(Q)=sum_p 2^c(p). More precisely,

    U(Q) <= 2^floor(m/2)m! - (2^floor(m/2)-1)H(Q).

This improves the earlier elementary bound2^(m-1)m! for m>=3, without assuming the Springer extremal conjecture.

## 4. Exact obstruction to saturating the unsigned relaxation

Always F(Q,b)<=U(Q). Equality holds if and only if:

- every vertex of Q has outdegree at most1; and
- every edge label b is+1.

The empty-edge cases are interpreted vacuously.

Proof. To saturate every relaxed permutation, put any controlled edge w,v at the end of a permutation. Its fold condition forces b_vw=+1. Thus all labels must be+. If v controls two distinct neighbors u,w, the ordering with consecutive u,v,w has incompatible pins, contradicting saturation. Conversely, outdegree<=1 eliminates all internal double-controllers, and labels+ satisfy every possible fold condition.

For m>=4, average outdegree=(m-1)/2>1, so equality is impossible for every Q. For m=3 the cyclic controller attains U=9; the transitive controller has U=11 but cannot attain it.

This disproves the initial realizable-coherence reduction. The failure is structural and occurs before any large enumeration.

## 5. Exact star-phase elimination

For each v with nonempty outgoing star, fix a representative bipartition of its outgoing edges, encoded by base signs a_vw, modulo simultaneous sign reversal. Write

    b_vw = epsilon_v a_vw,   epsilon_v in {+1,-1}.

All wedge factors are independent of epsilon_v, since both labels belong to the same star. A permutation has at most one fold factor. Therefore

    F(Q,b) = A(Q,a) + sum_v epsilon_v B_v(Q,a),                  (7)

with

    A = sum_p product_wedges (1-a_e a_f),
    B_v = sum_(p ends at v and v controls predecessor w)
                    a_vw product_wedges (1-a_e a_f).

Thus

    max_epsilon F = A + sum_v |B_v|,                            (8)
    min_epsilon F = A - sum_v |B_v|.

There are no phase-phase interaction terms. If E=binom(m,2) and r is the number of nonempty outgoing stars, this reduces the sign search from2^E assignments to2^(E-r) star bipartitions, with all2^r phases optimized algebraically. A tournament has at most one zero-outdegree vertex, so r is m or m-1 except the trivial case.

This is an exact computational reduction, not a proof of which Q or bipartition is globally extremal.

## 6. Even/odd order transfer, including multiplicities

The entire attainable fixed-count spectrum is the same on orders2m and2m+1. In even order the final fold signs are the internal mate arcs. In odd order the final fold signs are the arcs to the fixed vertex, while all m internal mate arcs are irrelevant to the fixed half-path count. Pair-interaction data are identical.

For canonical fixed involutions on the two labelled vertex sets, if N_even(r) and N_odd(r) count invariant orientations with fixed count r, then

    N_odd(r) = 2^m N_even(r).                                   (9)

Indeed each gauge-fixed controller/sign datum has2^m even realizations and2^(2m) odd realizations. This is a labelled fixed-involution census identity, NOT an assertion that isomorphism-class multiplicities differ by2^m. The normalized distributions and every normalized moment coincide.

Use F_m^max for the common maximum here. Do not identify it with the user's M_k: physical repository searches found several unrelated M_k usages and did not recover a unique relevant definition.

## 7. Finite exact extrema and a sharp transitivity counterexample

Complete enumeration of every controller Q and every edge-sign assignment for m<=5 gives:

| m | F_m^max | max_Q U(Q) | maximizing controller/sign gauges |
|---|---:|---:|---:|
|1|1|1|1|
|2|3|3|2|
|3|9|11|14|
|4|39|57|56|
|5|195|361|40|

The last column counts labelled pair-controller/sign gauges, not tournament isomorphism classes.

Independently, phase elimination (8) gives the same max-per-Q histograms:

- m3:9 on all8 labelled controllers.
- m4:37 on32 controllers;39 on32 controllers.
- m5:185 on24;189 on480;191 on240;193 on240;195 on40.

At m5 the transitive controller has maximum193. The global maximum195 is achieved exactly when Q is C3 -> T2. There are40 labelled such controllers; each has one maximizing sign assignment. Therefore **transitive controllers do not maximize the genuine signed problem**, already at five pairs.

A simple195 witness: label the cyclic triple0,1,2 with arcs0->1->2->0, append vertices3->4, and direct every triple vertex toward3 and4. Set all cyclic-triple labels+, all v->3 labels+, all v->4 labels- for v in{0,1,2}, and b_34=+. Its phase decomposition is

    A=156,   (B_0,B_1,B_2,B_3,B_4)=(10,10,10,9,0),
    F=156+10+10+10+9=195.

Explicit order10 and11 adjacency masks and anti-involutions are in independent_checks.json. They agree with the independently implemented ternary half-path DP, not just the new pin-count code.

For m3 a short algebraic maximum proof replaces enumeration. If Q is cyclic, F=6+sum_edges b<=9. If Q is transitive and a,b label the two edges controlled by its source and c labels the third edge,

    F=6+a+b+c-2ab<=9.

Opposite a,b and c=+ attain9. Orders8/9 and10/11 maxima above are finite exhaustive certificates, not all-order formulas.

## 8. The uniform-plus family and its limits

If every b=+, only permutations with no internal double-controller survive. Their Q-edge directions are a forward run followed by a backward run, with a unique sink. If I(Q) counts such permutations whose sink is internal, then

    F(Q,+1)=3H(Q)+2I(Q).                                        (10)

The sink-last and sink-first cases contribute H and2H respectively. Equivalently, this count is twice the convolution of two directed paths ending at a common sink, minus H(Q).

For transitive Q, F(Q,+1)=2^m-1. Each nonmaximum vertex chooses which side of the maximum to occupy, giving2^(m-1) sink-shaped permutations; one has sink last and weight1, the rest weight2.

Maxima of this uniform-plus subfamily over Q at m2,3,4,5 are3,9,27,105. It is not the optimizing family for m>=4. Its m5 maximum is attained by the regular five-vertex controller, whereas the actual195 maximum uses the nonuniform signs above.

## 9. What the Springer numbers really count here

For transitive Q, the relaxation U(Q) is the sum over permutations of2^(number of left peaks). A left peak includes the first position when it exceeds the second, but excludes the final position. This is the same free-pair statistic in the original unsigned description; reversing the total order if needed only exchanges peaks with valleys.

Let P_m(t) be this left-peak polynomial. Inserting the new largest entry into a permutation with k left peaks leaves k unchanged in2k+1 slots and increases it in m-2k slots. Hence

    P_(m+1)(t)=(1+mt)P_m(t)+2t(1-t)P'_m(t),  P_0=1.

The exponential generating function is

    sum_m P_m(t) z^m/m!
       = sqrt(1-t)/(sqrt(1-t)cosh(z sqrt(1-t))-sinh(z sqrt(1-t))).

Specializing t=2 gives

    sum_m U(T_m) z^m/m! = 1/(cos z-sin z),

so these are precisely the Springer numbers1,1,3,11,57,361,... . This identity is classical peak-enumeration mathematics, not a newly discovered generating function. A primary reference is Ma, Fang, Mansour, Yeh, 'Alternating Eulerian polynomials and left peak polynomials', arXiv:2104.09374v3, especially its introduction and left-peak sections:

https://arxiv.org/html/2104.09374v3

The claim that transitive Q maximizes U(Q) is verified here only for m<=5. It remains a separate relaxation extremal problem. Even a proof would give an upper bound, not the genuine maxima, because Section4 proves the relaxation cannot be saturated once m>=4.

## 10. Verification, lineage, and remaining obligations

Files:

- enumerate.cpp: exhaustive all-controller/all-sign count for m<=5, with complete fixed-count histograms.
- results.txt: exact output; includes every fixed-count frequency in the gauge-fixed labelled ensemble.
- check_model.py: constructs actual even/odd tournaments and compares the factor formula to the independent ternary DP.
- independent_checks.json:146 independent DP comparisons, explicit maximum witnesses.
- phase_opt.py: an independent algebraic optimization over star bipartitions.
- phase_results.json: complete max-per-Q histograms and maximizing controller codes.

At m5, direct enumeration has1,048,576 controller/sign gauges. Phase elimination examines43,008 star bipartitions and reproduces all max-per-Q values. No m>5 exhaustive enumeration was performed.

The prior repository reflection 'the-reversal-involution-and-the-half-principle-opus-20260629.md' already interprets fixed paths as half-paths and records Paley values1,9,185 at orders3,7,11. This report extends that existing half-path language; it does not rename or rediscover it. The195 order11 witness shows that the Paley185 value is not the unrestricted fixed-half-path maximum.

Open obligations:

1. Prove or refute transitive extremality for U(Q) for arbitrary m.
2. Determine F_m^max beyond the finite certificate, using star bipartitions and exact phase elimination rather than a false coherence assumption.
3. Develop endpoint/port recursions for structured controller extensions, especially the C3 -> T2 maximizing witness.
4. Relate fixed-count multiplicities to isomorphism-class and centralizer-weighted blue-tile counts; do not substitute the gauge histogram for a metagraph-node histogram.
5. Recover the repository's intended M_k definition before making any comparison of extrema or asymptotics.

The additions below close the global blue-spectrum existence question, but not the larger fixed-order fullness conjecture.


## 11. Binary universal-controller construction: every positive odd blue count

Let A(Q,b) be the wedge-only sum from Section5, omitting every fold factor. Add a new controller u with u->v for every old vertex, and label every outgoing edge of u by the same epsilon in{+1,-1}. Every ordering with u internal vanishes, since its wedge factor is1-epsilon^2=0. If u is first, the new count is the old F. If u is last, the old final fold is replaced by1+epsilon, and the contribution is(1+epsilon)A. Thus the exact extension laws are

    F_new=F+(1+epsilon)A,   A_new=2A.                            (11)

Starting from one pair with F=A=1 gives

    A_m=2^(m-1),
    F_m=1+sum_(j=2..m)(1+epsilon_j)2^(j-2).                     (12)

As the row signs epsilon_j vary, these are exactly the odd integers1,3,...,2^m-1, each once in this row-uniform subfamily. The controller Q is transitive.

To pass from fixed paths to blue tilings, note that a sigma-commuting automorphism must permute sigma-pairs, preserve the controller tournament, and preserve the positive fold sign. A transitive tournament has no nonidentity automorphism; the fold condition forbids an individual pair flip. Hence the centralizer is trivial. By the already proved blue-count formula B(T)=F_HP(T,sigma)/|C_Aut(T)(sigma)|, the blue count is exactly F_m.

Therefore every positive odd h is a blue count at both orders2m and2m+1 for

    m=ceil(log2(h+1)).                                          (13)

It is also attainable at every larger pair count by choosing later epsilon=-1. Literal adjacency/DP/centralizer checks verify blue7 at n6/7, blue21 at n10/11, blue31 at n10/11, blue63 at n12/13, and blue127 at n14/15. The latter checks use explicit individual constructions, not exhaustive enumeration past m5.

In particular the global positive blue-count spectrum is exactly the positive odd integers, using the established oddness of blue fibres. Non-self-converse classes additionally supply blue count0. The repunits2^m-1 are the all-plus row construction and are not rare in this observer.

This does NOT settle the ordinary Hamiltonian-count spectrum. It also does NOT prove THM-643's fixed-order C4 interval up to3^floor((n-2)/2), which is larger than the binary interval in relevant orders. THM-643 already records blue{1,3,5,7,9} through n7 and its earlier3-power upper cap was separately refuted by blue39 at n8/9. Those are lineage, not findings newly supplied here.

## 12. Odd-order mate-arc colouring and the exact blue maximum

Let G=Aut(Q,b). In even order the fold gauge fixes all mate arcs, and

    C_Aut(T)(sigma) is isomorphic to G.                          (14)

In odd order the m internal mate arcs form an arbitrary bit-vector d, while all other data are fixed, and

    C_Aut(T)(sigma) is isomorphic to Stab_G(d).                  (15)

Proof. The centralizer fixes the unique central vertex in odd order and permutes pairs. Preservation of the fold orientation forces each positive representative to map to a positive representative, so pair flips are impossible. The remaining preservation conditions are exactly the controller orientation, its edge signs, and, in odd order, the mate-arc vector. Conversely any such pair permutation lifts to a sigma-commuting automorphism.

The group G is a subgroup of a tournament automorphism group, so has odd order: an involution would transpose two vertices and reverse their mutual arc, and Cauchy's theorem excludes even group order. Gluck's Corollary1 states that an odd-order permutation group has a subset with trivial setwise stabilizer. Applying it to G's faithful action on the pairs, choose d to encode a regular subset. This kills the centralizer without changing a single half-path count.

Consequently EVERY signed-controller fixed count is realizable as a blue count at order2m+1, and in particular

    max_(|T|=2m+1) B(T) = F_m^max.                              (16)

This theorem uses the ordinary odd-order Gluck corollary, not the converse-extended regular-subset statement with exceptions C3,C5,P7. It therefore genuinely uses the repository's palette theme without identifying those exceptions with Hamiltonian-count gaps.

Primary source: David Gluck, 'Trivial Set-Stabilizers in Finite Permutation Groups', Canadian Journal of Mathematics35(1)(1983),59-67, Corollary1 on p66:

https://doi.org/10.4153/CJM-1983-005-2

The cited corollary relies on the odd-order theorem and Gluck's classification/induction. This report does not claim an elementary replacement proof. The repository's THM-4743 already discusses Gluck and its converse extension; the application here is specifically to unused mate-arc bits and the blue-count denominator.

At m5, the unique maximizing sign assignment on Q=C3->T2 has G=C3. Thus its even-order f195 gives blue65, not195. All40 maximum gauges have this symmetry. The transitive controller witness f193 has G=1. It follows from the finite complete fixed-count census that

    max_(|T|=10) B(T)=193,
    max_(|T|=11) B(T)=195.                                      (17)

The second value is attained by reversing one mate arc belonging to the cyclic triple, breaking C3 while preserving f195. check_blue.py verifies the actual adjacency automorphisms, not just the signed-controller description. At even orders2,4,6,8,10 the blue maxima certified here are1,3,9,39,193; at odd orders3,5,7,9,11 they are1,3,9,39,195.

Additional files: check_blue.py and blue_checks.json give all binary witnesses above and the193/195 centralizer controls. All claims concerning maxima above five pairs remain open.


## 13. Strong binary realizations and portable verification

The binary construction need not hide inside disconnected tournaments. When the new universal controller u has row sign epsilon=+, its cross arcs are exactly

    u- -> every old paired vertex -> u+.

In even order the fold arc u+->u- closes this into a cyclic wrap. Every old vertex is mutually reachable with u+ and u-, regardless of the old tournament's strong components. In odd order the fold passes through the central vertex,

    u- -> every old paired vertex -> u+ -> center -> u-,

so all vertices are mutually reachable whenever there is at least one old pair. These orientations follow directly from the convention that the controller has sign+ when it appears later.

For every odd h>=3, the minimal binary size m=ceil(log2(h+1)) has its highest binary digit positive, so its final row has epsilon=+. Consequently h has strong blue-count realizations at both orders2m and2m+1. For h=1, use C3 on the odd side: reverse the single mate arc relative to the positive central fold. Thus every positive odd h has a strong odd-order blue realization on at most2ceil(log2(h+1))+1 vertices. No assertion is made about the possibility or impossibility of a strong even-order blue1 realization.

The independent audit packet session29-signed-audit proves this extension as well as a faster pair-mask fixed-path algorithm using O(m^2 2^m) arithmetic operations. That compression is a separate audited algorithm; this packet's verification still calls the earlier ternary DP and does not relabel that implementation as the faster one.

The DP dependency in check_model.py now resolves the sibling directory session29-halfpath relative to this packet, or an explicit HALFPATH_DP_PATH environment override. It no longer depends on an absolute workspace path. A copied sibling layout was executed successfully before refreshing the hashes.

check_strong.py and strong_checks.json additionally verify literal strong connectivity and fixed-path counts for all63 odd h from3 through127 at both minimal binary orders, plus the h1 C3 case:127 constructed graphs in total. This is not an exhaustive tournament enumeration.
