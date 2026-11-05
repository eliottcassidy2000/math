# Every single-module tournament context is positive in central-factorial coordinates

Session 197 follow-on to `PROOF.md`. **PROVED all-order operator positivity and explicit physical H13/H19 escapes; inverse targeting remains OPEN.** No literature-priority claim is made. The fixed bin-profile conventions, odd-cycle identity, q/r bases, and S193/S195 source links are in the companion proof. This file does not alter its restricted-grammar omission statements: it adds new physical generators outside those grammars.

## 1. The statement is about the full operator matrix

For a rooted tournament (Q,r), let C_Q[S]=Q[r:S] replace r by the entire tournament S, with uniform exterior incidences inherited from r. Empty S means delete r. There is a unique linear map A_Q on Q[x] such that

    A_Q B_S = B_(C_Q[S])

for every tournament S, including the empty tournament. To see full-profile linearity, cut each output cover into maximal S-runs. These form an unordered cover of S. All actual runs are distinct fragments, and their exterior incidence types are the same because S is an intact module. For a fixed number j of fragments the possible recombinations with the named exterior vertices have a fixed output bin polynomial K_j(x), independent of which actual S-cover supplied the fragments. Thus

    A_Q (x)_j = K_j(x)

defines the linear map on the falling-factorial basis. The case j=0 is the bin profile of Q-r. This is the one-module specialization of the inherited S195 cover kernel. Equivalently, A_Q(x^m)=B_(Q[r:Tr_m]) for all m>=0, which determines the map since monomials form a basis. No scalar-H-only compression is being assumed.

Consider the full basis

    1, q_0,r_0,q_1,r_1,...,
    q_k=x product_(j=1)^k(x^2-j^2),
    r_k=x^2 product_(j=1)^k(x^2-j^2).

**Theorem A (universal single-module positivity).** For every actual rooted tournament (Q,r), the matrix of A_Q in this basis is entrywise a nonnegative integer matrix. If Q has N exterior vertices, A_Q shifts polynomial parity by N modulo two. Consequently its restriction from the appropriate q- or r-state space to the output parity is nonnegative.

The theorem applies to every individual basis vector, most of which are not physical tournament profiles. It is therefore stronger than merely observing that physical inputs give positive physical outputs. It is not a realization theorem for arbitrary vectors or arbitrary positive matrices. Nor does it assert simultaneous multivariable positivity for several independently virtual input modules.

## 2. Three elementary positive operators

Let X multiply by x, let S shift x to x+1, and put

    delta = S-S^(-1),
    Delta = S+S^(-1)-2I.

The symbol S here is an operator, not an input tournament. In this section operator composition acts from right to left.

Besides X1=q_0 and the positive X actions in the companion proof, one has, for k>=1,

    Delta q_k = 2k(2k+1) q_(k-1),
    Delta r_k = 2(k+1)(2k+1) r_(k-1),
    delta q_k = 2(2k+1) r_(k-1),
    delta r_k = 4(k+1)q_k + 2k(k+1)(2k+1)q_(k-1).                 (1)

The boundary formulas are

    Delta 1=delta 1=0,
    Delta q_0=0,       delta q_0=2,
    Delta r_0=2,       delta r_0=4q_0.                            (2)

The constants in (2) are actual multiples of the separate basis element 1. One must not omit them, or extend the formulas to a fictitious q_(-1).

**Proof.** In the first three identities, evaluate the shifted polynomial at the consecutive integer roots shared by the relevant basis element, then compare degrees and leading coefficients. For example Delta q_k has degree 2k-1, is odd, and vanishes at +/-1,...,+/-(k-1), so it is a multiple of q_(k-1); the second finite difference fixes the coefficient 2k(2k+1). The even identities also vanish at zero when k>=1, with a double zero by parity. Finally r_k=xq_k gives

    delta r_k=x delta q_k + q_k(x+1)+q_k(x-1),

and (1) follows using Delta q_k and Xr_(k-1)=q_k+k^2q_(k-1). Direct calculation proves (2). QED.

Thus X, delta, and Delta all preserve the cone of nonnegative integer coordinates in the extended basis. Delta preserves polynomial parity and delta and X reverse it.

## 3. The exact odd-cycle skeleton formula

The following establishes Theorem A, including all signs and multiplicities. Begin with the actual graph Q[r:Tr_m]. Use the inherited identity

    B_T(x)=sum_Gamma 2^|Gamma| x^[|T|-sum_(C in Gamma)(|C|-1)]        (3)

over vertex-disjoint odd directed cycle packings.

Every cycle either stays entirely outside the clone module or contains clones. A cycle entirely inside the transitive clone module is impossible. In each crossing cycle replace every maximal nonempty clone-run by a marked gap. The remaining cyclic word contains distinct exterior vertices; its marked gaps record where the clone-runs were. Unmarked consecutive exterior vertices must have the corresponding arc in Q. A marked gap from u to v requires u->r->v in Q.

Call a collection of such exterior-only cycles, crossing cyclic skeletons, and unused exterior vertices a skeleton packing. Exterior vertex sets are disjoint. Exterior-only cycles must already be odd. A crossing skeleton need not have odd exterior length. If it has l exterior vertices and t marked gaps, its eventual cycle is odd precisely when the total number of clones assigned to its t gaps has parity

    epsilon = 1-l (mod 2).                                     (4)

There are finitely many skeleton packings, depending only on Q. Canonical rotation starts every cycle/skeleton at its least exterior vertex; ordering the objects by those anchors removes any ambiguity about their unordered collection. The marked gaps of a skeleton are distinct positions. Assigning a nonempty set of named clones to each gap determines that run uniquely by its increasing transitive order. Thus cutting runs and reconstructing them are inverse, with no rotation factor, path-order factorial, or unlabelled-clone division omitted.

For t>=1 define the Laurent shift operator

    L_(t,epsilon) = (S-I)^t + (-1)^epsilon (S^(-1)-I)^t.           (5)

For an integer d>=0, half of L_(t,epsilon)(x^m) evaluated at d counts assignments of m named clones to d unused singleton bins and t specified nonempty run bins, with the total assigned to the run bins having parity epsilon. Indeed inclusion-exclusion for the required run bins gives (S-I)^t. Giving each occupied run-bin choice sign -1 instead gives (S^(-1)-I)^t, whose contribution for an assignment is (-1) to the total number assigned to runs. Averaging with (-1)^epsilon is the exact parity projector. This argument is a signed polynomial evaluation, not a claim that d-j is a negative number of bins.

For several crossing cycles, the corresponding operators multiply. Their shifts commute: each clone chooses exactly one of the distinct marked-run bins or an unused singleton bin. Inclusion-exclusion and the parity sign are applied separately to each cycle's group of bins. Each crossing cycle's factor 2 from (3) cancels the factor 1/2 in its projector. If a skeleton packing has

- c_0 exterior-only cycles;
- c_1 crossing skeletons;
- u unused exterior vertices;

its operator contribution is therefore

    2^c_0 X^(u+c_0+c_1) product_(i=1)^c_1 L_(t_i,epsilon_i).       (6)

Crucially X acts **after** all shift operators. The x exponent records the unused exterior vertices and one bin-profile factor for each packed cycle; unused clone vertices are already counted inside x^m by their d available labels. Moving X through a shift would give a different and generally wrong operator.

Summing (6) over all skeleton packings proves the exact identity

    A_Q = sum_(skeleton packings) 2^c_0 X^(u+c_0+c_1)
                                      product_i L_(t_i,epsilon_i). (7)

First it proves equality on every monomial x^m, including m=0, by the physical clone construction and (3). The full-profile linearity of Section 1 then proves equality on all Q[x], in particular every nonphysical basis vector q_k or r_k. This is the step that upgrades a graph-count identity to a full operator-matrix theorem.

### Why all the signs in (5) become positive operators

Put y=S-I and z=S^(-1)-I. These commute and satisfy

    y+z=Delta,        yz=-Delta,        y-z=delta.

Let P_t=y^t+z^t. Then

    P_0=2, P_1=Delta,
    P_t=Delta(P_(t-1)+P_(t-2)),       t>=2.                      (8)

Let R_j=(y^(j+1)-z^(j+1))/(y-z), understood as the polynomial geometric sum. Then

    R_0=1, R_1=Delta,
    R_j=Delta(R_(j-1)+R_(j-2)),       j>=2.                     (9)

Both are polynomials in Delta with nonnegative integer coefficients. Hence

    L_(t,0)=P_t(Delta),
    L_(t,1)=delta R_(t-1)(Delta).                               (10)

For example L_(1,0)=Delta, L_(1,1)=delta, L_(2,0)=Delta^2+2Delta, and L_(2,1)=delta Delta. This checks both exterior parities at t=1 and t=2 rather than implicitly fixing one sign. Equations (1)-(2), (7), and (10) prove entrywise nonnegative integrality.

For parity, the shift contributed by a skeleton is

    u+c_0+c_1 + sum_i epsilon_i
      = u+c_0 + sum_i l_i (mod 2).

Every exterior-only cycle has odd length, so this equals the total number N of exterior vertices modulo two. Thus every summand has the required output parity. This completes the proof of Theorem A.

One immediate consequence for every physical input T is the coefficientwise central-factorial inequality

    state(C_Q[T]) >= H(T) state(C_Q[Tr_p]),

where p=1 for odd |T| and p=2 for even |T|. It follows by keeping only input coordinate u_0. This is a coordinate inequality, not deletion or a claimed physical decomposition. The usual rooted-template inverse problem remains.

## 4. Explicit new generators reach H13 and H19

For a nonempty core A, a nonempty port set P in A, and clock T, form U_(A,P)[T] with

    w -> T -> A,
    q -> w for q in P,       w -> q for q outside P.

All old arcs are retained. The S193 multi-port theorem gives

    B_U=x[(B_A-K)f - L f(x-1)+R f(x+1)],                      (11)

where f=B_T, L sums first-bin starts at ports in x bins, R sums first-bin ends at ports in x+1 bins, and K sums the two-distinct-bin port responses in x+1 bins. These polynomials are pinned to actual cores, not fitted to a few values.

### A triangle with two ports

Use the directed triangle rows [4,1,2], ports {0,1}. For each vertex the start response is x^2-x+1, and the end response at x+1 is x^2+x+1. For the two ordered port pairs the mixed responses are x+1 and x-1, giving K=2x. As B_A=x^3+2x, (11) gives

    Ff=x^4 f+2x[(x^2+x+1)f(x+1)-(x^2-x+1)f(x-1)].             (12)

A particularly transparent positive operator form is

    F=X^4+4X^2+2X^2 Delta+2(X^3+X)delta.                       (13)

Alternatively F=2XDX+2D-X^4-2X^2, consistent with the companion tridiagonal operator. At x=1,

    H(U)=h+6 B_T(2)
        =13h+36u_1   for odd input,
        =25h+72u_1   for even input.                          (14)

For T=K1 the actual strong five-vertex tournament is

    [20,17,2,7,12],       H=13.

Thus it physically escapes the H13 omission of the previous D/E grammar while preserving positive coordinates.

### A transitive four-vertex core with two nonconsecutive ports

Take A with rows [0,1,3,7], whose order is 3->2->1->0, and ports {0,2}. The responses are

    L=(x-1)^3+x^2(x-1)=2x^3-4x^2+3x-1,
    R=(x+1)^3+x^2(x+1)=2x^3+4x^2+3x+1,
    K=x(x+1)+x(x-1)=2x^2.

For L, a first bin starting at sink 0 excludes all three other vertices; a first bin starting at 2 allows vertices 1,0 in that bin or the other bins but excludes vertex 3. The end cases give R. For mixed port order (0,2), the remaining vertices have x+1 and x choices; for (2,0) they have x-1 and x choices. This proves K without fitting.

With B_A=x^4, (11) yields

    Gf=x[(x^4-2x^2)f
          +(1-3x+4x^2-2x^3)f(x-1)
          +(1+3x+4x^2+2x^3)f(x+1)].                         (15)

Again the positive form is short:

    G=X^5+6X^3+2X+(X+4X^3)Delta+(3X^2+2X^4)delta.            (16)

At x=1,

    H(U)=-h+10 B_T(2)
        =19h+60u_1   for odd input,
        =39h+120u_1  for even input.                          (17)

With T=K1 its actual strong six-vertex rows are

    [32,1,35,7,15,26],       H=19.

Every core vertex reaches port 0, then w, then T, while T reaches the core, proving strongness for any nonempty clock in this second family. The triangle core gives the same conclusion in the first family.

These are exact positive physical additions to the grammar, not a proof of arbitrary prime targeting. A theorem allowing an arbitrary rooted Q does not solve the inverse problem of finding Q with a specified count. In particular it does not justify adding arbitrary positive state transformations or treating the entire positive cone as realizable. The negative-Hankel witnesses in the companion proof still apply; Theorem A concerns matrix entries, not positivity of every continued moment matrix.

## 5. Reproduction and hostile controls

Run alongside the companion verifier:

    python3 verify_extension.py > extension_normal.stdout
    python3 -O verify_extension.py > extension_optimized.stdout
    cmp extension_normal.stdout extension_optimized.stdout

The standard-library checker uses explicit -O-safe exceptions. It independently enumerates the rooted-template cycle skeletons by exterior cyclic permutations and marked gap sets, then builds the Laurent-shift operator and compares against full bin profiles of actual clone graphs. That expected path uses neither the positive lowering theorem nor the port formulas.

Its finite controls include all 285 rooted labelled tournaments of orders 1 through 4; actual clone sizes m=0,...,4 give 1,425 full physical profile comparisons, including the deletion boundary. There are 1,140 matrix-entry sign checks on the first q/r basis vectors. Literal assignment enumeration separately checks 168 parity-projector cases, with m=0,...,6, t=1,...,3, both parities and d=0,...,3 (the exact count is emitted by the verifier). Elementary lowering identities and F/G positive actions are checked through index 30. Both new port operators are checked on every labelled clock of order at most four, 150 additional actual graphs, and the explicit H13/H19 rows are independently counted by endpoint DP.

The theorem is proved by the finite skeleton formula and positive recurrences, not inferred from these finite controls. The companion packet and this follow-on together separate physical constructions, exact state transport, arbitrary-formal-state positivity, and the still-unproved inverse-targeting step.
