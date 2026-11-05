# An injective three-response signature and sharp bilinear source bounds

Session 197. This is a separate corollary to the complete start-one source grammar and the two-source C2 normal form. It replaces four numerical source responses by three different linear observables, with a nonlinear exact inverse on the physical class. It also sharpens product pruning by the sources' Hamilton-path counts. No scalar coverage, source-graph injectivity, or closed linear three-dimensional transition is asserted.

## 1. Physical sources, conventions and inherited inputs

A source is a nonempty rooted tournament (V,b) with exactly one directed Hamilton path starting at b. Write

    s=(D,E,U,d)=(H(V), e_b(V), R_b(V), H(V-b)).

Here R_b counts ordered covers (P,Q), where Q is nonempty and starts at b and P may be empty. Counts are literal directed vertex words; H(empty)=1. The singleton has response o=(1,1,1,1). A nontrivial source means |V|>=2, not necessarily strong.

The reviewed ROOTED and BACKWARD_PROFILE_INVERSE inputs give the complete strong start-one grammar. For its ordinary gadget M_L,

    (D,E,U,d) -> (g(L)D-(L mod2)E+LU,
                  D-[L=1]E, Q(L)D+B(L)U, D).         (1)

The constants have g(1)=Q(1)=B(1)=1 and, for l>=1 in the even line and l>=2 in the odd line,

    (g,Q,B)(2l)=(g(l),Q(l),2B(l)),
    (g,Q,B)(2l-1)=(g(l)+l-(l mod2),Q(l)+B(l),2B(l)).

Appending a sink while retaining the root gives O(s)=(D,0,2U,d). The reviewed weak-source theorem says every weak source is uniquely a strong start-one core (possibly the singleton) followed by a nonempty transitive tail. These operations therefore generate every source response and only physical ones.

All source D and d are odd. This can be checked within the grammar: g and Q are odd, B is even except at L=1, and the invariant E=U modulo2 is preserved by (1), by O, and at the singleton. Thus D remains odd; the new deletion count is an earlier D or an unchanged d. No parity assertion about arbitrary formal tuples is used.

For strong sources with D>=3 the inherited cone is

    D<=U<=2(D-E),  E>=1.                            (2)

They have one of these last-gadget descriptions:

- a nonsingleton last gadget: d=E;
- a singleton last gadget: d=U-D+E and 3E>=D.

The second inequality follows from the exact last-step theorem for D>3. The only D=3 source is the triangle, which has E=d=1 and U=3 and can be placed in the first description. The singleton also has d=E. A weak source has E=0.

The named copied dependencies are ROOTED.md, BACKWARD_PROFILE_INVERSE.md, WEAK_START_ONE_ORACLE.md and its unchanged exact oracle. The separate BILINEAR_C2_INVERSE.md supplies the physical meanings and completeness of the N/M grafts used in Section5. Their hashes are in RELEASE_FILES.sha256; these are inputs, not claims inferred from finite tests here.

## 2. Exact recovery from three linear observables

Define

    A=D+2d,  B=E+d,  C=U-d.                         (3)

**Theorem 1.** On the entire physical source-response set, including weak sources and the singleton, (3) is injective. Its inverse is the following exact piecewise-linear readout:

    if B is odd:
        d=B             when A>=3B;
        d=A-B-C         when A<3B;
    if B is even:
        d=B/2           when 2A<=3B+2C;
        d=A-B-C         when 2A>3B+2C.

Then recover

    D=A-2d, E=B-d, U=C+d.                           (4)

**Proof.** In a weak source E=0, hence B=d is odd and A-3B=D-d>=0. In a strong source with d=E, B=2d is even. In a strong singleton-last source, a possible odd B is distinguished from every weak source by

    A-3B=D-3E-d<=-d<0.

For even B put S=2A-3B-2C. On the sheet d=E,

    S=-2(U-D)<=0.

On the sheet d=U-D+E,

    S=U-D>=0.

At S=0 both sheet formulas give the same d. These inequalities and the preceding odd-B distinction exhaust the inherited physical source alternatives. The singleton (A,B,C)=(3,2,0) is recovered by the even branch, as is the triangle. This proves every branch and injectivity. QED.

Theorem1 concerns numerical response quadruples. Many nonisomorphic rooted graphs can have the same response; no graph-isomorphism conclusion follows. Also the formula is only a candidate inverse on unrestricted integer triples. The coordinate inequalities alone are not reachability conditions.

**Exact membership.** Compute the candidate in (4), reject impossible signs, oddness or elementary bounds, then run the reviewed weak-source oracle on (D,E,U). Accept precisely if d belongs to the returned set of possible deletion counts. The oracle returns one physical source recipe for that d. Theorem1 and the soundness/completeness of that oracle prove that this is exact membership in the image of (3), with witness reconstruction. It remains a numerical-target algorithm, not a polynomial-bit-time claim.

This changes the observables rather than claiming the old triple (D,E,U) always determines d. For example the two reviewed weak sources (9,0,18,1) and (9,0,18,3) collide in that old triple but have respective signatures (11,1,17) and (15,3,15). The known arbitrarily large weak deletion ambiguity is likewise retained by (3).

## 3. Rank, integral content and bounded inventory

In the original coordinates the symmetric N form has matrix

    [[1,-1,1,0],[-1,1,0,-1],[1,0,3,-1],[0,-1,-1,0]].

Its kernel is spanned by (-2,-1,1,1), which is exactly the kernel of the linear map (3). Consequently the form descends to

    N(s,t)=(A-B)(A'-B')+AC'+CA'+3CC'.                (5)

The descended symmetric matrix is

    [[1,-1,1],[-1,1,0],[1,0,3]],

with determinant -1. Thus N has rank3. Equivalently, in X=A-B,B,C the matrix is [[1,0,1],[0,0,1],[1,1,3]]. Its entries are nonnegative, but its entire integer cone is not asserted to be physical.

On physical sources the signature preserves exact integral content:

    gcd(A,B,C)=gcd(D,E,U,d).                         (6)

One divisibility direction follows from (3). For the other, gcd(A,B,C) is odd because A is odd. Each recovery branch makes d an integer combination of A,B,C or B/2; an odd common divisor also divides B/2 when that branch is used. Equation(4) then proves the other direction. This is a gcd identity, not a representability or prime-coverage theorem.

For U<=K the elementary bounds d<=D<=U and E<=D give

    1<=A<=3K, 1<=B<=2K, 0<=C<=K.

By injectivity, there are O(K^3) physical response quadruples under this cap. An exact forward inventory needs no graph-order cutoff: in (1), U_new=Q D+B U>U and B(L)>=L, so any accepted last gadget satisfies L<=K/U. Try those finitely many L, retain U_new<=K, and include weak dyadic tails. This enumerates the whole capped physical inventory. Predecessor recipes can be retained. The estimate improves the crude four-coordinate state bound, but does not turn (1) into a linear three-coordinate action: recovering d is essential and piecewise.

## 4. Elementary physical inequalities with no signed cancellation

For every nontrivial source,

    E<=d,  D-E>=d,  U>=D,
    f:=U-D+E-d>=0.                                 (7)

Here are direct injections, independent of numerical source enumeration.

Deleting b from a path ending at b injects the E paths into paths of V-b, proving E<=d. Since a source has a b-start Hamilton path and |V|>=2, b has an out-neighbor. In any path of V-b, insert b immediately before its first out-neighbor. Every earlier vertex points to b, so the preceding vertex, if present, points to b. The insertion is valid and not at the end; deleting b recovers the path. Hence D-E>=d.

Cutting each full path immediately before b injects all D paths into the U covers. Covers whose two pieces concatenate to a valid path are exactly those D images. Among the other covers, take Q to be the singleton b and P to be a path of V-b. Exactly E of those d choices have P's last vertex pointing to b, since those and only those extend to a path ending at b. Therefore at least d-E covers fail to concatenate, proving U-D>=d-E and the last inequality in(7).

In particular E<=D/2. Since D is odd, E<=(D-1)/2. This sharper integer bound, and the nonnegative slack f, are crucial below. Neither argument assumes a global upper bound on U at fixed D; weak tails can make U arbitrarily large.

## 5. Sharp bilinear Hamilton-count product bounds

For two **nontrivial** physical sources s_i=(D_i,E_i,U_i,d_i), the independently derived two-source C2 grafts have

    N=(D1-E1)(D2-E2)+D1 U2+U1 D2+3U1 U2
        -d1(E2+U2)-d2(E1+U1),

    M=D1 D2+2(D1 U2+U1 D2)+6U1 U2-E1 U2-U1 E2
        -d1(E2+U2)-d2(E1+U1).                       (8)

Their physical construction and complete C2 interpretation are inherited from the separate normal-form proof; the following inequalities also follow algebraically from (7).

**Theorem2.** Every such pair satisfies the exact integer bounds

    4N >= 15D1 D2+5D1+5D2-1,
    2M >= 17D1 D2+3D1+3D2-1.                       (9)

The leading constants 15/4 and 17/2 cannot be increased uniformly.

**Proof.** Put b_i=U_i-D_i>=0 and f_i=U_i-D_i+E_i-d_i>=0. Direct expansion gives

    N0=6D1 D2-2(D1 E2+E1 D2)-E1 E2,
    N=N0+b1(3D2-2E2)+b2(3D1-2E1)+b1 b2
            +f1(E2+U2)+f2(E1+U1);

    M0=11D1 D2-2(D1 E2+E1 D2)-2E1 E2,
    M=M0+b1(7D2-3E2)+b2(7D1-3E1)+4b1 b2
            +f1(E2+U2)+f2(E1+U1).

All omitted slack terms are nonnegative by(7). Both N0 and M0 decrease separately with E1 and E2. Substituting Ei<=(Di-1)/2 gives exactly(9). This proves the bounds without capping or cancelling signed coordinates.

For sharpness, let j>=1 and x=2^j+1. Apply gadget L=2^j to the singleton, then a singleton gadget. Their literal source has

    (D,E,U,d)=(2x-1,x-1,2x,x).

Indeed the first source is (x,1,x,1), as follows from g(2^j)=Q(2^j)=1 and B(2^j)=2^j; then (1) gives the displayed source. Two identical sources in(8) give

    N=15x^2-2x, M=34x^2-6x+1.

Divide by D^2=(2x-1)^2 and let j tend to infinity. The ratios approach15/4 and17/2. These are actual C2 constructions by the graft theorem, so the limiting constants are physically sharp, not merely sharp on a relaxed cone. QED.

The nontrivial-source restriction is necessary: the singleton fails D-E>=d and is deliberately not admitted in(9). Weak sources such as a transitive two-vertex source are allowed.

For target n, scalar symmetry permits D1<=D2. Inequalities(9) then bound the smaller source by exact integer comparisons

    15D1^2+10D1-1<=4n   for N,
    17D1^2+6D1-1<=2n    for M.

The two-source normal form also has independent U-product bounds; retain those when enumerating sources, since a bound on D alone does not bound weak tails. Together these inequalities are sound inverse pruning, not a construction for every numerical target.

## 6. Reproducibility and scope

Run `python3 verify.py` and `python3 -O verify.py` in a copied packet directory. The verifier writes VALIDATION.json. It uses explicit checks, not removable assert statements. It checks:

- all32549 physical response quadruples in the complete U<=2000 grammar bank: exact recovery, collision freedom, gcd identity and coordinate bounds;
- every184512 integer candidate in the stated U<=31 signature rectangle against exact physical membership;
- every1695 start-one root among all labelled tournaments through order5, using actual path/cover counts;
- all10609 ordered nontrivial source pairs in the U<=41 bank: both slack identities, both product bounds, and the signature form;
- 1275 instances in the unbounded weak-ambiguity construction through parameter50;
- three literal sharpness sources and both resulting graft graphs each, through order14, including their exact C2 counts.

Normal and optimized receipts agree byte-for-byte. These finite checks supplement the proofs. They do not certify that every cone point is physical, that gcd1 implies coverage, that the signature identifies a graph, or that a fixed linear signature transfer exists. The exact oracle remains responsible for source reachability, and the all-order low-cycle scalar coverage question remains open here.
