# Capacity reorientation transport and score-fibre Ehrhart certificates

Session 23b, 2026-10-21. Status: **PROVED** in the stated scopes; the accompanying integer controls are separately **FINITE-EXACT**. This is an integral-capacity lift of the original fixed-Hamiltonian-path nonpath-complement operation, with an explicit change of demand. It is not an assertion that the partner is regular or that a circulation stays a circulation. Classical total unimodularity, lattice-polytope counting, and reciprocity are imports, not novelty claims.

## 1. Objects, conventions, and exact reorientation transport

Fix a labelled undirected graph with edge set E and an orientation T. Every coordinate remains attached to its undirected edge. Let B_T have column e_v-e_u on the arc u->v, so B_T f is inflow minus outflow. For integer M>=0 and integer demand vector b, put

    F_T(M,b) = {f in Z^E : 0<=f_e<=M, B_T f=b}.

The same notation with real coordinates denotes its polytope. An element of F_T(M,0) is a bounded nonnegative circulation. It is not an Euler tour and need not have connected support.

Let D be any edge subset, and reverse exactly those arcs, obtaining T'. Define

    (Phi_D^M f)_e = M-f_e  if e in D,
                     f_e  otherwise.

Write B_D 1 for the sum of the columns of B_T indexed by D. Then

    B_T' Phi_D^M(f) = B_T f - M B_D 1.                 (1)

Indeed, outside D the contribution is unchanged. On D its new contribution is -B_e(M-f_e)=B_e f_e-M B_e. Thus Phi_D^M is a bijection

    F_T(M,b) <-> F_T'(M,b-M B_D 1).                    (2)

It is an ambient integral affine isometry: its linear part is a diagonal matrix with entries +/-1 and its translation is integral. In the reverse orientation B_D'1=-B_D1. Therefore applying the same coordinate rule a second time returns both f and b. This proves involutivity with the changed incidence and demand, not with an incorrectly held-fixed demand.

For a tournament on n vertices let s(T) be its labelled **outdegree** vector. Since

    B_T 1=(n-1)1-2s(T),
    B_T'1=B_T1-2B_D1,

we have B_D1=s(T')-s(T), and (2) becomes

    b'=b+M(s(T)-s(T')).
    M s(T)+b = M s(T')+b'.                            (3)

Thus the exponent vector M s+b is the invariant of the decorated transport. In particular, tournaments with the same labelled score vector have affinely lattice-isomorphic bounded-circulation polytopes at every M. Relabeling gives the corresponding conclusion for the same score multiset. The ambient affine maps are not generally tournament isomorphisms and need not conjugate their automorphism-group actions.

The same argument permits unequal integer capacities c_e, replacing M-f_e by c_e-f_e and the demand correction by -B_D c. Uniform capacities are what make the score-vector specialization (3) especially simple.

## 2. The original fixed-path nonpath-complement operation

Fix an actual directed Hamiltonian path

    P=(v_0 -> v_1 -> ... -> v_(n-1))

in T, and let D=E(T)\E(P). Reverse D and preserve every path edge, obtaining T'. Since

    B_D1=B_T1-B_P1,
    B_P1=e_(v_(n-1))-e_(v_0),

(2) gives the exact operation map

    F_T(M,b) <-> F_T'(M,b-M B_T1+M(e_last-e_first)).    (4)

The selected Hamiltonian path and its labelled endpoint ports are literally unchanged. For a regular tournament, B_T1=0, and therefore

    F_T(M,0) <-> F_T'(M,M(e_last-e_first)).            (5)

This sends circulations to capacity-M flows with M units of net output at the first endpoint and M units of net input at the last endpoint. They may also contain cyclic flow; no assertion identifies a flow vector with one path or one traversal.

For n=2r+1, the partner's labelled outdegrees are

    s(T')=r1-B_P1:
    r+1 at first, r at internal vertices, r-1 at last. (6)

So T' is not regular for n>=3. For general T the score law is

    s(T)+s(T')=(n-1)1-B_P1.                           (7)

This is the existing repository score-complement skeleton, expressed in the present first-to-last path convention. The older THM-646 uses the opposite base-path direction; its endpoint constants consequently appear reversed. Equations (4)-(5) are a correctly typed integral lift of that operation.

A minimal hostile is a directed triangle with its two-edge Hamiltonian path fixed. Nonpath complement makes the triangle transitive. The source has M+1 circulations; the target has only one circulation but exactly M+1 flows of the endpoint demand in (5). At n=5 the exact controls give 24 source circulations and 24 endpoint-demand flows at M=1, but only 14 target circulations. At M=2 the corresponding numbers are 219,219,90.

## 2A. The zero-on-path face: the actual Hamiltonian-presentation count

Retain the distinguished directed Hamiltonian path P and now impose f_e=0 on every path edge. For a regular tournament T define

    R_(T,P)(M) = |{f in F_T(M,0): f|P=0}|.

If T and U are regular and both contain this same labelled P, their disagreement set avoids P. Hence (2) preserves the zero-on-P condition and proves that this count depends only on n and M, not on the regular base. Write it R_n(M). A regular base containing any prescribed labelled path exists: use the cyclic tournament in the path order.

At M=1, reversing the support of f gives a regular tournament U that contains the same P, and every such U occurs exactly once. Thus

    R_n(1) = number of regular labelled tournaments containing fixed P
           = sum_(regular isomorphism classes T) H(T)/|Aut(T)|.    (20)

For the last equality, label an unlabelled representative along each of its Hamiltonian words. Two words give the same labelled tournament precisely when an automorphism carries one to the other. The action on Hamiltonian words is free, since fixing the ordered list fixes every vertex. This is the actual path-normalized presentation count, not an analogy between flows and paths. The bijection depends on the selected regular base and does not canonically identify one circulation with one Hamiltonian path of that base.

For odd n>=5, the digraph T minus P is strongly connected. Indeed, if S is any nonempty proper vertex subset with |S|=a, regularity gives exactly

    a(n-a)/2

outgoing cut arcs. A directed Hamiltonian path uses at most min(a,n-a) of them: each tail in S and each head outside S can occur at most once. Since max(a,n-a)>=3,

    a(n-a)/2 - min(a,n-a)
       = min(a,n-a) [max(a,n-a)/2-1] > 0.

Every proper cut therefore has an outgoing arc after deleting P, which is equivalent to strong connectivity. Every remaining arc lies in a directed cycle; summing one such cycle flow for every arc and then scaling gives a circulation strictly between zero and one on every remaining edge. Consequently the zero-on-P polytope has full dimension in the remaining circulation kernel:

    delta=(m-(n-1))-n+1=m-2n+2=(n-1)(n-4)/2.                     (21)

It is an integral face of P_T (or follows from the same incidence integrality argument). Ehrhart theory proves that R_n(M) is a polynomial of exact degree delta and

    sum_(M>=0) R_n(M) z^M = a_n(z)/(1-z)^(delta+1).               (22)

The all-ones vector on the remaining arcs is generally not balanced. Therefore the full-polytope interior shift by two, and its specific palindromic degree law, do **not** automatically transfer to this face.

There is nevertheless an exact replacement: its **codegree is three**, so a_n has exact degree delta-2. At dilation one there is no interior integer coordinate assignment. At dilation two the only possible one is 1_D, but B_D1=-p, where p=e_last-e_first, so it is not a circulation. Strong connectivity supplies a directed first-to-last path Q using only D. The vector 1_D+1_Q has divergence -p+p=0 and all coordinates equal to 1 or 2, so it is interior at dilation three. Series reciprocity then makes the lowest interior exponent equal to delta+1-deg(a_n), proving the degree claim. No palindromicity follows. This refinement was independently derived during the packet audit.

For cyclic R5 with positive differences {1,2}, delete P=(0->1->2->3->4). The remaining graph has arcs 02,24,40,41,13,30. Its circulation semigroup is free on the cycles

    A=(0->2->4->0),
    B=(0->2->4->1->3->0).

To verify this explicitly, put a=f_40 and b=f_41. Balance forces f_13=f_30=b and f_02=f_24=a+b, with no other coordinates. Thus capacities are equivalent to a,b>=0 and a+b<=M, proving

    R_5(M)=binom(M+2,2),
    sum_(M>=0) R_5(M) z^M=1/(1-z)^3.                            (23)

For this representative the total arc-use weight is 3a+5b, which is extra information not captured by (23). Its exact joint capacity/use series is

    sum_(M>=0) sum_(a,b>=0, a+b<=M) z^M q^(3a+5b)
       = 1/[(1-z)(1-z q^3)(1-z q^5)].                           (27)

The other two regular path-normalized n=5 presentations already differ at capacity one: the three total-use polynomials are 1+q^3+q^5, 1+q^4+q^5, and 1+q^3+q^4. Thus the same underlying regular tournament with different selected Hamiltonian paths can have different weighted face data, although the unweighted Ehrhart count is universal. Direct edge-coordinate controls give R_5(M)=1,3,6,10,15,21,28,36,45 for M=0,...,8. At n=3, deleting P leaves no directed cycle, so R_3(M)=1; neither the strong-connectivity claim nor formula (21) is asserted there.

At n=7, delta=9 and R_7(1)=91. Independent enumeration of all 2^15 orientations of the nonpath edges finds 91 regular presentations, split by Hamiltonian-path count as

    H=171: 57 presentations;
    H=175: 25 presentations;
    H=189:  9 presentations.

The exact flow enumeration for each of three displayed regular bases containing P gives the same 91. There is also a decisive nonpalindromicity control: a_7 has constant coefficient 1, linear coefficient 91-10=81, degree 7, and leading coefficient 34. To obtain the last value, subtract 1_D from an interior lattice point of the capacity-three face. The result is a binary zero-on-P flow of demand p in T minus P. Independent star-coordinate and Hamiltonian-path-tree eliminations both count exactly 34 such flows. Series reciprocity identifies this number with the leading h* coefficient. Thus this face is not palindromic, in contrast to the full regular circulation polytope. The top coefficient also has a direct score-fibre interpretation: reversing the support of that binary flow g in the original regular base gives a tournament still containing P with score vector r1+p, namely outdegree r-1 at the first endpoint and r+1 at the last. A separate fixed-path orientation enumeration gives 1 such presentation at n=5 and 34 at n=7, agreeing with both interior counts. This verifies the value in (20); a census of all regular isomorphism classes is not needed to prove the affine invariance.

Finally, the original operation D=E\P preserves zero path coordinates literally. Restricting (5) therefore gives

    {f in F_T(M,0): f|P=0}
        <-> {g in F_T'(M,M(e_last-e_first)): g|P=0}.              (24)

This refines the integral lift to the actual path-presentation carrier while retaining the correct nonzero target demand and both labelled endpoint ports.

## 2B. Centered signed currents: an exact electrical-coordinate statement

Set j=f-(M/2)1. Under coordinate complementation on D,

    j'_e=-j_e on D and j'_e=j_e otherwise,
    B_T' j'=B_T j=b-(M/2)B_T1.                                  (25)

Thus the centered signed current and its divergence represent the same edge-based physical quantities after reorientation. The invariant divergence equals b+Ms-M(n-1)1/2, precisely the exponent invariant (3) minus a constant vector. If fixed undirected-edge resistances r_e>0 are retained, then

    sum_e r_e (f_e-M/2)^2                                       (26)

is preserved. Ordinary uncentered energy sum r_e f_e^2 need not be preserved. In the regular-to-nonpath-complement case, centered divergence is zero on both sides: the target's nonzero f-demand is offset by its changed baseline incidence. This is a genuine change-of-orientation statement, not a claim that every feasible current obeys Ohm's law, arises from vertex potentials, or minimizes energy.

## 3. The score-fibre coefficient formula and the M=1 bijection

For every tournament T, every M>=0, and every integer b, the bounded-flow count is

    |F_T(M,b)| = [product_i x_i^(M s_i+b_i)]
                product_(i<j) (sum_(k=0)^M x_i^k x_j^(M-k)).       (8)

To prove this, give an original arc u->v carrying f_e the monomial x_u^(M-f_e)x_v^f_e. The exponent at vertex i is M times its outdegree, minus its total outgoing flow, plus its incoming flow: M s_i+b_i. Changing which endpoint is called the tail merely reorders the M+1 terms of the factor, so the product in (8) is orientation-independent. Terms have coefficient one. A coefficient with any negative exponent is zero.

At M=1 a circulation is the indicator of an arc set D with B_D1=0. Reversing D gives a tournament with the same labelled score vector, and every tournament in this score fibre arises uniquely by comparing its arcs with T. Thus

    |F_T(1,0)| = number of labelled tournaments with score vector s(T). (9)

These are directed Eulerian subsets (balanced at every vertex), not merely undirected even-degree subsets. In particular the fibre is not thereby a binary linear code or a coset. Its size can be 24, which is not a power of two.

For M>1, (8) counts multiplicity allocations **without binomial weights**. It is not the number of orientations of M distinguishable parallel copies of every edge. That different count replaces each factor by (x_i+x_j)^M. At n=3,M=2, the unweighted circulation count is 3, whereas the distinguishable-parallel-edge count is 1+2^3+1=10.

## 4. Regular tournaments: one capacity polynomial for each odd order

Let n>=3 be odd, T regular, and m=n(n-1)/2. Write

    P_T={f in R^m : 0<=f<=1, B_T f=0},
    d=m-n+1=(n-1)(n-2)/2.

The incidence matrix has rank n-1 because the underlying graph is connected. The point (1/2)1 belongs to its kernel and lies strictly between all coordinate bounds. Hence P_T has affine hull ker B_T and dimension exactly d.

P_T is an integral polytope. One standard proof uses total unimodularity of oriented incidence matrices and integral bounds. A direct vertex proof is also short: if a feasible vertex had any nonintegral coordinate, the undirected graph of its nonintegral edges could have no vertex of degree one, since the integral demand equation there would express one noninteger as an integer. It therefore contains an undirected cycle. A sufficiently small signed circulation perturbation on that cycle preserves all balances and bounds in both directions, contradicting extremality. Thus every vertex is integral. This also proves integrality for integer demands and unequal integer capacities.

Ehrhart's theorem now gives a polynomial

    N_n(M)=|M P_T intersect Z^m|,       degree exactly d,             (10)

independent of the chosen regular T by (2). Its ordinary series has the form

    E_n(z)=sum_(M>=0) N_n(M) z^M = h_n*(z)/(1-z)^(d+1).             (11)

Here h_n* is the standard lattice-polytope numerator; its nonnegative coefficients are classical. This is a polynomial rather than a nontrivial quasipolynomial because P_T is integral in its circulation lattice.

Every coordinate is nonconstant on the affine hull: 0 and 1 are both feasible. Since (1/2)1 is strictly feasible, the relative interior of M P_T is exactly

    {f : B_T f=0, 0<f_e<M for every e}.

For integral f and M>=2, subtracting 1 from every coordinate gives the bijection

    relint(M P_T) intersect Z^m <-> (M-2)P_T intersect Z^m.          (12)

It preserves balance specifically because B_T1=0. The M=1 interior count is zero; at M=2 the all-ones vector is its unique lattice point. Ehrhart-Macdonald reciprocity therefore yields

    N_n(-M)=(-1)^d N_n(M-2)                       for M>=2,         (13)
    N_n(-1)=0.

Since both sides of (13) are polynomials and agree at infinitely many integers, the identity extends polynomially to all arguments. In the centered variable t=M+1, N_n(t-1) is even when d is even and odd when d is odd.

The series of interior counts is z^2 E_n(z). The series reciprocity law gives

    z^2 E_n(z)=(-1)^(d+1) E_n(1/z).

Substituting (11) and clearing denominators proves

    h_n*(z)=z^(d-1) h_n*(1/z).                                    (14)

As h_n*(0)=1, it has **exact degree d-1**, is palindromic, and has leading coefficient 1. No separate Gorenstein theorem is needed for this deduction. At n=3, d=1, N_3(M)=M+1 and h_3*=1, agreeing with (13)-(14). The n=1 edgeless exception has N_1(M)=1 and h_1*=1; it has no edge inequalities and must not be assigned degree d-1=-1 or codegree two.

These arguments also apply to a nonempty graph with an Eulerian orientation: d=m-n+c where c is the number of connected components. Tournament regularity is the complete-graph specialization.

## 5. Exact five-vertex polynomial, with two counting routes

For n=5 the dimension is six. The exact answer is

    N_5(M)=(M+1)^2 [3(M+1)^4+5(M+1)^2+4]/12                       (15)

or, in the ordinary monomial basis,

    M^6/4 + 3M^5/2 + 25M^4/6 + 20M^3/3
            +79M^2/12 + 23M/6 +1.

Consequently

    E_5(z)=(1+17z+72z^2+72z^3+17z^4+z^5)/(1-z)^7.                (16)

This is a proof by finite exact determination after the proven degree bound, not an extrapolation from an unproved pattern. Here is the coefficient recurrence used to calculate the determining values. For a complete graph on r remaining vertices, let C_M(a_1,...,a_r) denote the coefficient with those target exponents. Set C_M(empty)=1 and set it to zero whenever any a_i<0, any a_i>M(r-1), or sum a_i != M r(r-1)/2. Then

    C_M(a_1,...,a_r)
      = sum_(0<=k_j<=M, sum_(j=2)^r k_j=a_1)
          C_M(a_2-M+k_2,...,a_r-M+k_r).                         (17)

The recurrence follows by choosing all r-1 factors incident with the first vertex and then deleting that vertex. Targets can be sorted for memoization because the remaining complete-graph product is symmetric; this does not divide by any symmetry factor or identify coefficient terms.

Apply (17) to (2M,2M,2M,2M,2M). The resulting values for M=0,...,10 are

    1, 24, 219, 1136, 4175, 12216, 30429,
    67264, 135621, 254200, 449031.

The first seven determine the polynomial of degree at most six uniquely and give (15). The last four are independent capacities beyond the interpolation range. A second implementation assigns the six non-star edge coordinates directly and solves the four star coordinates from the balance equations; it agrees through M=8. It does not use recurrence (17) or the polynomial.

## 6. What the certificate forgets, and what restores information

All regular tournaments on a fixed odd order have the same N_n(M), despite different Hamiltonian-path counts, automorphism groups, and weighted flow distributions. For two explicit regular tournaments on Z/7,

- Cyclic tournament with positive differences {1,2,3}: H=175.
- Paley tournament with positive differences {1,2,4}: H=189.

Both have N_7(1)=2640. The coefficient of q^4 in their bounded-capacity-one total-use polynomials is 28 and 21 respectively. These controls therefore show genuine information loss, not just a logical possibility.

The full multivariate capacity polynomial does transport exactly. For

    K_T(M,b;y)=sum_(f in F_T(M,b)) product_e y_e^f_e,

(2) implies

    K_T'(M,b';y) = (product_(e in D) y_e^M)
                  K_T(M,b; y_e^-1 on D, y_e otherwise).          (18)

The total use changes by

    |Phi_D^M f|=|f|+M|D|-2 sum_(e in D) f_e.                     (19)

Therefore simply setting every y_e=q does not generally preserve total-use series. Retaining the reversed-edge load is a sufficient sidecar for this one-step weight transformation. Retaining T, the selected path P, D, and labelled demand is necessary to interpret the operation as the original metagraph transport rather than just an equality of scalar counts.

Equation (9) counts a **labelled score fibre**, not tournament isomorphism classes. It does not solve A000568, classify the original metagraph's nodes or ports, establish connectivity, or count Euler tours. Polytope equivalence does not restore the forgotten tournament isomorphism type. Its useful positive role is an exact score-fibre certificate and an integral transport across the original path-complement edge.

## 7. Inheritance and primary mathematical imports

The targeted repository inheritance search found:

- `01-canon/theorems/THM-646-line-metagraph-score-complement-skeleton.md`: the pointwise labelled-score complement law and score-fibre skeleton of the original line operation. Our formula (7) is its existing law in the explicitly stated reversed path convention.
- `05-knowledge/results/tournament_even_transfer_20261011b.md` and the current session's oriented-flow packet: binary cycle-space versus integral nonnegative-flow distinctions. Their loss-of-orientation guardrail remains active.
- `01-canon/MISTAKES.md`, 2026-10-11 Vitali-origin correction: score fibres are not binary code cosets; the labelled regular n=5 fibre has size 24.
- `01-canon/theorems/THM-087-worpitzky-tournament.md`: an earlier Ehrhart analogy, not the present bounded-incidence polytope.

The contribution here is the explicit affine capacity/demand map on the actual nonpath-complement operation and the resulting correctly scoped score-fibre polynomial certificate. This local search does not establish literature priority for general reorientation transformations or for the explicit polynomial.

Classical imports, checked at primary sources:

1. A. J. Hoffman and J. B. Kruskal, *Integral Boundary Points of Convex Polyhedra*, in *Linear Inequalities and Related Systems* (1956), pp. 223-246. The original article is reproduced with an author introduction in *50 Years of Integer Programming* (2010), pp. 49-76. The introduction explicitly states the integer-right-hand-side / total-unimodularity integrality theorem. Our elementary vertex argument above also proves the special case directly.
   https://www.cs.umd.edu/~gasarch/BLOGPAPERS/kruskalhoffman.pdf
2. R. P. Stanley, *Combinatorial Reciprocity Theorems*, Advances in Mathematics 14 (1974), 194-253, Section 6, especially Proposition 6.1 and the degree discussion on printed pp. 208-210; Section 5, Proposition 5.2 relates reciprocal evaluations and generating functions. These are the reciprocity and polynomial-degree imports used here.
   https://math.mit.edu/~rstan/pubs/pubfiles/22.pdf
3. M. Beck and F. Sottile, *Irrational proofs for three theorems of Stanley*, European Journal of Combinatorics 28 (2007), 403-409, arXiv v4. Introduction and Theorem 2 give the Ehrhart denominator and numerator nonnegativity, with lattice period p=1 here; Theorem 1 gives cone reciprocity.
   https://arxiv.org/html/math/0501359v4

## 8. Reproduction and finite scope

Run from the packet directory:

    python3 check.py > normal.json
    python3 -O check.py > optimized.json
    cmp normal.json optimized.json

The program has explicit failure checks, not optimization-removable assertions. It computes the exact n=5 counts through M=10 and independently through M=8; verifies all 291 labelled score vectors across all 1024 n=5 tournaments at capacity one; checks all 262,144 triples (K4 orientation, reversed edge set, capacity-one vector) including their arbitrary demands; compares complete source and target flow sets for every regular path-normalized n=3 and n=5 orientation at capacities 1,2,3, totaling 4,146 transported source flows; and independently counts the n=7 regular score fibre by coefficient recurrence and star-coordinate flows for the displayed representatives. It also checks the zero-on-path face through capacity 8 at n=5 and capacity 1 at n=7, and enumerates all 32,768 path-normalized n=7 tournaments to verify the 91 regular presentations and their H distribution. It additionally counts the 34 binary endpoint-demand flows certifying the n=7 path-face leading h* coefficient. It does not enumerate large n=7 capacities or claim a full n=7 isomorphism census.
