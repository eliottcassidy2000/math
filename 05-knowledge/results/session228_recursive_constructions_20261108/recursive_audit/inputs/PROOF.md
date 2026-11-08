# Exact target-bounded recursive cyclic construction

**Status: PROVED local theorems; FINITE-EXACT compilation, not yet independently audited.**
Chan–Pak completeness remains OPEN. No unrestricted coverage improvement or literature-priority claim is made.

## 1. Inheritance and precise new question

The current source checkout was main `12fb69799f72a83fdbfa5e7900090fa7ee240913`.
The complete path-cover substitution kernel is inherited from
`01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md`.
The full cyclic profile / endpoint-jet calculus is inherited from
`01-canon/theorems/THM-3134-tournament-endpoint-jet-and-c3-newton-profile-transform.md`.
The multiplicative path-colour coordinate is inherited from
`01-canon/theorems/THM-3166-falling-factorial-order-join-path-colour-transform.md`.
The older THM-1975 discussion is not used in place of those proved results.

Session227's multi-root finite-menu theorem explicitly leaves recursive reuse outside its scope. This report treats a genuinely recursive grammar: either child of a new cyclic node can be an arbitrarily nested earlier output. The new scoped consequences are an exact **target-bounded decision procedure**, exact minimum cyclic-depth compilation, and a depth-dependent minimal rational linear summary theorem. The binary polynomial identity itself is a short specialization of the inherited kernel, independently obtained in the sibling arithmetic lane. It is included with a direct proof for self-containment, not presented as a new general substitution theorem.

Inheritance board: physical graph construction; nonlinear polynomial multiplication; recursively growing moments; target rather than graph-order bounds; exterior SCC products. The closest mechanism is the inherited path-cover transform. The canonical hostile is equal low moments with unequal deeper cyclic responses. The corrected near miss is applying a strong-node order bound to arbitrary exterior joins. The least-used sidecar is the next colour value Q(3), then the full ladder Q(4),Q(5),....

For a nonempty tournament T, Q_T(j) counts ordered j-bin spanning path covers, allowing an empty bin once. Equivalently Q_T(t)=sum_d pc_T(d)(t)_d, where pc counts unordered nonempty-path covers. Thus Q_T(0)=0, Q_T(1)=H(T). For the transitive tournament TT_n, Q(t)=t^n. All polynomials below use ordinary powers of t.

## 2. Exact physical recursive operations

Let J(A,B)=A=>B be the ordered join. Let C(A,B)=C3[A,B,K1], with arcs A->B->z->A. Both inputs must be nonempty. The latter graph is strong regardless of the internal connectivity of A and B. Write p=Q_A, q=Q_B and define

    Delta+ p(t)=p(t+1)-p(t), Delta- p(t)=p(t)-p(t-1).

The exact laws are

    Q_J(A,B)(t)=p(t)q(t),
    Q_C(A,B)(t)=t[p(t)q(t)+Delta+ p(t) Delta+ q(t)
                              +Delta- p(t) Delta- q(t)].       (1)

The first is inherited multiplicativity. To verify the second at a positive integer j, choose the colour of z in j ways. The other j-1 colours contain an optional A path followed by an optional B path. In z's colour, the only possible layout is

    A_pre, B_pre, z, A_post, B_post.

The pre/post A slots and the j-1 ordinary slots provide j+1 A bins; likewise for B. The illegal assignments are A_pre nonempty with B_pre empty, or B_post nonempty with A_post empty. Inclusion-exclusion gives

    p(j+1)q(j+1)
    -[p(j+1)-p(j)]q(j)-p(j)[q(j+1)-q(j)]
    +[p(j)-p(j-1)][q(j)-q(j-1)].

The subtraction terms are independent exactly because deleting one specified empty bin leaves j ordinary labelled bins. Simplification gives the bracket in (1). Every legal assignment yields a unique coloured cover and conversely. Equality at every positive integer proves the polynomial identity.

Put h_T=Q_T(1) and r_T=Q_T(2)-Q_T(1). Evaluating (1) gives the positive scalar rule

    H(C(A,B))=2h_A h_B+r_A r_B.                            (2)

This is a physical nonlinear arithmetic operation, not an arbitrary operation on formal pairs. Updating r_C needs Q_A(3),Q_B(3), so (h,r) is not asserted to be a closed state.

For the unary collar U(A)=C(A,K1), (1) specializes to

    Q_U(A)(t)=t Q_A(t+1)+t^2 Q_A(t)-t Q_A(t-1).            (3)

In particular H(U(A))=Q_A(1)+Q_A(2).

## 3. Sharp strong-node order bound and exact target compiler

Let |A|=a, |B|=b. Partitioning into two bins gives

    r_T=H(T)+sum_(empty != S != V(T)) H(T[S])H(T-S)
        >=2^|T|-1.                                      (4)

Every tournament, including an empty induced subgraph under H(empty)=1, has a Hamilton path. Therefore

    H(C(A,B)) >= (2^a-1)(2^b-1)+2
                 >=2^(a+b-1)+1.                         (5)

For the second inequality set x=2^(a-1), y=2^(b-1); the difference between the two displayed lower bounds is 2(x-1)(y-1)>=0. Both bounds hold for every nonempty A,B. The first is attained by transitive inputs. The second is attained, for every total output order n>=3, by one singleton and one transitive input of order n-2. Thus the uniform output-order bound is sharp:

    H(C(A,B))<=X  ==>  |C(A,B)|<=N(X),
    N(X)=2+floor(log_2(X-1)), X>=3.                       (6)

**Theorem (exact finite decision procedure).** Fix a finite bank S of actual nonempty seed tournaments containing K1, and close it under arbitrary finite applications of J and C. For each integer X>=3 the entire attainable scalar set at most X, and each attainable value's minimum cyclic nesting depth, can be computed by a finite procedure:

1. Compute each seed polynomial. Enumerate, by increasing graph order n<=N(X), all polynomial profiles of grammar expressions of order n and H<=X. Include seeds of that order; combine all smaller ordered-join pairs whose orders sum to n, and all cyclic pairs whose orders sum to n-1. Deduplicate polynomials while storing minimum cyclic depth. Seed depth is 0, join depth is max of child depths, cyclic depth is one plus their maximum.
2. Form the scalar generator bank from every seed's H<=X (including a seed larger than N(X), if present) and every cyclic output obtained from the enumerated input profiles with total output order at most N(X). Store the minimum possible depth for each generator scalar.
3. Take the multiplicative closure of these generators through X. A dynamic program over increasing positive integers records reachability, and minimum of max depths over factorizations. Ignore generator 1 when descending to a strict smaller integer; it cannot improve the minimum depth below 0.

**Proof.** All profiles at each fixed finite order are obtained by their final grammar operation; induction on n proves step 1 complete. H is nondecreasing from a child to its parent, because join multiplies positive integers and (2) is at least either child's H. Thus pruning H>X cannot remove a successful target expression.

Every cyclic node in an expression whose final H is <=X satisfies (6). Its children have smaller order and occur in step 1, irrespective of their nesting depth. An arbitrary expression whose outermost operation is a join can be flattened only across its outer joins, until the factors are seeds or cyclic nodes. Their H values multiply, and each factor appears in step 2. Conversely every listed profile and generator has a literal graph expression; joining chosen generator witnesses constructs any accepted product. This proves both directions of scalar membership. For the depth refinement, profile equality preserves every future J/C polynomial response by (1), so replacing a profile representative by a minimum-depth representative never loses a smaller-depth continuation. The max/min dynamic programs therefore give the exact minimum depth. ∎

This is a total exact algorithm for this specified grammar, not a runtime bound efficient in log X, a recognition algorithm for arbitrary tournaments, or an existence proof for all odd targets. It does not assume a fixed list of gadget letters or a bound on recursive depth.

**Essential hostile.** No order bound like (6) holds for the entire grammar: C3=>TT_m has H=3 at arbitrarily large order. The compiler handles these arbitrary outer joins by scalar multiplication, rather than falsely discarding every expression larger than N(X). The literal four-vertex case already violates the naive N(3)=3 bound.

## 4. Exact finite recursive-depth findings

The canonical verifier uses as seeds **all labelled tournaments of orders 1 through 5**, reduced only by equality of their complete Q polynomials. There are 1,099 labelled presentations and 15 distinct order/profile pairs. At X=4095, N(X)=13. Arbitrary recursion and arbitrary exterior joins yield exactly **577** distinct attainable odd scalars. Their exact minimum cyclic-depth distribution is

    depth 0: 79; depth 1: 390; depth 2: 107; depth 3: 1.

The first values requiring depth 1,2,3 are respectively 17,241,4001. These depth statements are **FINITE-EXACT consequences of the complete target algorithm**, not unrestricted minimum graph orders or new scalar attainability results.

Two particularly simple recipes are

    T0=TT_4, T1=U(T0), T2=U(T1), T3=U(T2),
    H(T0),H(T1),H(T2),H(T3)=1,17,241,4001.

T2 has 8 vertices and outgoing-neighbour masks

    [126,188,248,240,224,192,130,1].

T3 has 10 vertices and outgoing-neighbour masks

    [510,764,888,1008,992,960,896,772,514,1].

The script constructs and directly counts these graphs. Despite full recursive reuse, this grammar still misses 19, alongside 7,21 and many others. This is a scoped membership result, not a nonexistence claim for tournaments. There is no coverage conjecture asserted for this grammar.

## 5. Recursive depth requires successively deeper linear moments

For k>=0, let f_k(B)=H(U^k(B)), where U^0(B)=B. Formula (3) implies

    f_k(B)=sum_(j=1)^(k+1) a_(k,j) Q_B(j),
    a_(k,k+1)=k!.                                      (7)

The top coefficient arises only from k successive upward shifts, of weights 1,2,...,k. The evaluation at zero never causes a new unknown because Q_B(0)=0 and the operator (3) vanishes at t=0. For example

    f_0=Q1,
    f_1=Q1+Q2,
    f_2=-Q1+5Q2+2Q3.

**Theorem (exact linear summary dimension).** For arbitrary physical nonempty B, the collection f_0,...,f_k has minimal rational linear-summary dimension k+1. Precisely, assigning B a vector in Q^d with each of these responses a rational linear functional requires d>=k+1. The moment vector (Q_B(1),...,Q_B(k+1)) attains equality.

**Proof.** The coefficient matrix in (7) is lower triangular with nonzero diagonal 0!,1!,...,k!, hence invertible. Evaluate its moments on the physical tournaments B=TT_n, n=1,...,k+1. The evaluation matrix j^n is invertible: factor j from column j (or transpose and factor the corresponding row) to obtain a Vandermonde matrix on distinct nonzero j. Consequently the physical response matrix has rank k+1. A d-dimensional linear summary factors that matrix through Q^d, forcing d>=k+1. ∎

Thus no one fixed finite-dimensional rational **linear-response** model works for all these recursive depths. This sharpens the local bounded-depth statement for this particular chain of physical contexts; it does not contradict the inherited general full-profile calculus. It does not rule out nonlinear encodings into one integer, adaptive branch-dependent state selection, or other constructions. It is not a Chan–Pak obstruction.

### Concrete physical low-moment collision

The compiler supplies two order-nine physical realizers with

    p(t)=t^9+14t^7+56t^5+28t^3,
    q(t)=t^9+16t^7+46t^5+36t^3.

Their outgoing rows are respectively

    [510,252,504,496,480,384,288,320,2],
    [508,505,506,240,224,384,288,320,24].

Both have (H,Q2)=(99,4320), but Q3 is respectively 64665 and 66825. One collar gives H=4419 for both; two collars give respectively 150831 and 155151. These are literal graph counts, not merely different formal polynomials. Their recipe strings and seed rows are included in CHECK.json. This refutes treating (H,Q2) as sufficient for arbitrary recursive continuation, while preserving its sufficiency for a single next scalar evaluation.

## 6. Reproduction, scope, and next positive question

Run in this directory:

    python3 -I -B verify.py > replay.json
    cmp replay.json CHECK.json
    python3 -I -B -O verify.py > replay.optimized.json
    cmp replay.optimized.json CHECK.optimized.json

The verifier is standalone standard-library Python; it reads no repository files, imports no scratch module, and uses explicit exception gates active under -O. The receipt's scalar/depth digest is over sorted exact `(value, minimum depth)` pairs, not graph isomorphism types.

Independent counting paths inside the checker:

- All 1,099 labelled seed graphs through order 5 have profiles counted by permutation backward-adjacency cuts, then checked against subset Hamilton DP plus independent bin convolution at every integer 0,...,n.
- All 121 pairs of labelled nonempty A,B through order 3 test the entire literal binary profile, scalar rule, and order bound. All 75 labelled graphs through order 4 separately test the entire unary profile.
- Exact recursive compilation uses the proved operation laws, with no graph-isomorphism or unproved feasibility filter. Minimum depth and missing targets are claims solely within this grammar and seed bank.
- The two depth witnesses, low-moment collision and both closing collars are checked by literal Hamilton subset DP. Physical response matrices have exact rational rank k+1 through k=8; (7) and the dimension theorem are proved for every k rather than extrapolated.

The target compiler is the constructive survivor: it replaces incomplete fixed-height searches with complete target membership in a truly recursive language. The positive question now is which additional physical prime seeds or operations admit an inverse target rule that covers the remaining integers. Merely adding more finite seeds or quoting scalar factorizations does not prove that rule. Chan–Pak completeness remains OPEN.
