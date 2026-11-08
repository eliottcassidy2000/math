# Binary cyclic synthesis and a finite physical inverse

**Status: PROVED candidate, pending independent audit.** All-order statements below
are proved combinatorially. The companion receipt is FINITE-EXACT. Chan–Pak
all-odd completeness remains OPEN. No unrestricted coverage interval or
literature priority is asserted.

## Inheritance and purpose

The general positive walk-content substitution kernel is inherited from
`01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md`;
the complete cyclic profile calculus is inherited from
`THM-3134-tournament-endpoint-jet-and-c3-newton-profile-transform.md`.
The order-join product for Q is inherited from
`THM-3166-falling-factorial-order-join-path-colour-transform.md`.
Session227's multiroot packet gives positive covers and sharp gadget-size
bounds. The four-state convention comes from Sessions225–227. The
Session199 balanced-pair product recorded in THM-4916 is a different
bilinear observation and is not claimed anew. The formulas here are useful
explicit specializations of the inherited substitution kernel; the main
point is an exact **physical numerical inverse** with a finite terminating
search, rather than an additional linear feasibility test.

The closest proved mechanism is positive path-cover substitution. The
hostile is a numerical cone/lattice point without a physical representative.
The corrected near miss is assuming (H,Q(2)) is closed under recursive
cyclic composition: its next Q(2) needs Q(3). The underused sidecar is the
ordered two-nonempty-path cover count. Concept board: positive substitution,
full colour polynomial, factor arithmetic, physical inverse, graph modules.

## 1. Physical operation and complete colour law

For nonempty tournaments A and B form C=C3[A,B,{u}] by directing
A -> B -> u -> A, retaining both induced graphs. This is strong, since
paths through the other two nonempty blocks connect every ordered pair of
vertices. Let H(T) count directed Hamilton paths. Let Q_T(j) count ordered
j-bin spanning path covers, allowing empty bins, and put
F_T(2)=Q_T(2)-2H(T). Thus Q_T(0)=0 and Q_T(1)=H(T).

For every integer j>=1, and hence as a polynomial identity,

    Q_C(j) = j [ Q_A(j) Q_B(j)
               +(Q_A(j+1)-Q_A(j))(Q_B(j+1)-Q_B(j))
               +(Q_A(j)-Q_A(j-1))(Q_B(j)-Q_B(j-1)) ].       (1)

**Proof.** Choose the bin containing u in j ways. The other j-1 bins
contain A then B, with either part possibly empty. Split the distinguished
bin into the portions before and after u. In either portion, A precedes B.
Ignoring boundary constraints gives Q_A(j+1)Q_B(j+1) assignments to the
j-1 ordinary slots and the two special slots for each factor. The illegal
case “A-before nonempty but B-before empty” has count
(Q_A(j+1)-Q_A(j))Q_B(j). The illegal case “A-after empty but B-after
nonempty” has count Q_A(j)(Q_B(j+1)-Q_B(j)). Their intersection has count
(Q_A(j)-Q_A(j-1))(Q_B(j)-Q_B(j-1)): each factor has one forbidden special
slot and its other special slot is required nonempty. Inclusion–exclusion
and rearrangement give (1). Slot assignment and path concatenation are
inverse operations. Each Q is a polynomial (sum of nonempty cover counts
times binomial polynomials), so agreement at all positive integers proves
the polynomial identity.

All differences in (1) are nonnegative at positive integers, with a direct
meaning: covers in which a specified extra bin is used. This is a positive
composition, despite its expanded shifted-polynomial form. It needs one
additional colour evaluation for each recursive layer. Setting B=K1
recovers the independently developed collar law
j Q_A(j+1)+j^2 Q_A(j)-j Q_A(j-1).

## 2. Rooted state, not just the scalar

Use x(T,u)=(D,E,U,d), where D=H(T), E counts Hamilton paths ending at u,
d=H(T-u), and U counts ordered pairs (P,R) spanning T with P possibly
empty and R nonempty starting at u. Put

    a=H(A), b=H(B), r=F_A(2), s=F_B(2).

Then

    d=E=ab,
    D=3ab+as+br+rs,
    U=D+br.                                                (2)

Equivalently, k_A=Q_A(2)-H(A)=a+r and k_B=b+s give

    D=k_A k_B+2ab.                                         (3)

**Proof.** A spanning path uses one of the six block-run words
uAB, ABu, BuA, ABuA, BuAB, ABuAB. Their counts are respectively
ab, ab, ab, rb, as, rs. A path ending at u must have word ABu,
so E=ab; deletion leaves A->B, giving d=ab.
For U, initially distribute A and B between the P-bin and the R-bin
following u. The R-bin can contain A then B. If its A-part is nonempty,
there are (Q_A(2)-a)Q_B(2) choices. If its A-part is empty, its B-part
must be empty as u loses to B; this contributes ab. Thus
U=(a+r)(2b+s)+ab=D+br.

For example, A=TT2 and B=TT3 give state (23,1,25,1),
whereas swapping the factors gives (23,1,29,1).

The scalar symmetry in A,B does not extend to U: switching A and B
changes U-D from br to as. This is a physical, asymmetric sidecar.

## 3. Exact constructive inverse for a prescribed rooted state

Let P be the set of all pairs (H(T),F_T(2)) of nonempty finite tournaments.
For an integer vector x=(D,E,U,d), define

    t=U-D,       z=2D-U-3d.

**Theorem.** There exists a physical C3[A,B,{u}] with state x iff
E=d>0, t>=0, z>=0, and at least one positive divisor b of gcd(d,t)
satisfies

    a=d/b,        r=t/b,        s=bz/(d+t),                 (4)

with s an integer and both (a,r),(b,s) in P. Here gcd(d,0)=d.

**Necessity.** By (2), t=br, d=ab and
z=as+rs=(a+r)s. Hence b divides d and t, a+r=(d+t)/b,
and (4) follows. All coordinates in the factor pairs are actual counts.

**Sufficiency.** Choose physical representatives A,B for the two pairs.
Their construction (2) has deletion count ab=d and br=t. Also
(a+r)s=z. Thus its Hamilton count is 3d+t+z=D, its U is D+t=U,
and its E=d=E. This supplies a literal physical representative.

There are at most tau(gcd(d,t)) arithmetic candidates. Physical pair
membership is essential: cone, integrality and divisibility alone do not
replace it. The t=0 boundary is harmless. A nonempty T has F_T(2)=0 iff
T is a singleton: if its order is at least two, any nontrivial vertex
partition gives two spanning paths. Thus t=0 forces A=K1, a=1 and b=d.
Likewise z=0 forces B=K1.

**Physical hostile to arithmetic feasibility alone.** The literal tournament
with outgoing masks (20,1,2,7,14), rooted at vertex 1, has state
(9,3,9,3), verified independently by subset DP in the receipt. It is strong:
0->4->3->2->1->0 is a spanning cycle. Here t=z=0, so the inequalities,
parity and rational divisibility conditions admit the candidate factor pairs
((3,0),(1,0)) and their swap. Neither is physical because F_T(2)=0
forces H(T)=1. Thus this genuinely physical rooted state is outside the
binary cyclic state image, although its scalar count 9 belongs to that image.
Physical pair membership is a load-bearing step, not redundant arithmetic.

This theorem tests existence of some graph with the requested state,
not membership of a given graph in the construction. A graph-level test
would also require u's out-neighbours and in-neighbours to be nonempty
modules A,B with all A->B. State equality alone does not certify those
edges. No uniqueness of physical representatives is claimed.

## 4. Removing the apparent infinite membership oracle

If T has n vertices, every nonempty proper subset and its complement
have Hamilton paths. Consequently

    F_T(2)>=2^n-2,   H(T)>=1,
    Q_T(2)-H(T)>=2^n-1.                                   (5)

For factors of orders m,n, (3) gives

    D >= (2^m-1)(2^n-1)+2 >= 2^(m+n-1)+1.                 (6)

The second inequality is exactly
2(2^(m-1)-1)(2^(n-1)-1)>=0. In particular a representation with D<=X
has total order at most 2+floor(log_2(X-1)); each factor separately
has order at most floor(log_2(X-1)). These statements require X>=3.
The bound is sharp for one singleton and one transitive factor.

Thus (4) is an effective finite physical inverse: enumerate all labelled
tournaments up to the individual order cap, retaining one representative
of each pair in P; test each divisor candidate against that bank, and
return C3[A,B,{u}] when it succeeds. The graph search may be expensive,
but terminates and is exact. There is no unproved realizability oracle.
Do not apply the order cap to arbitrary tournaments with H<=X; it is a
cap for factors in this particular strong cyclic construction.

For a scalar target N the corresponding exact inverse enumerates
1<=d<=floor(N/3), factorizations d=ab and N-2d=k_A k_B,
and tests (a,k_A-a),(b,k_B-b) in the same finite physical bank.
All factors are positive; actual membership enforces all further parity
and size restrictions. It returns witnesses or proves absence within
this operation, not absence among unrestricted tournaments.

## 5. Useful boundaries and exact controls

The operation produces 11 from A=B=TT2 and 33 from A=C3,B=TT2.
It does not produce 13, even with arbitrarily large factors: (6) implies
m+n<=4 for D=13. The possibilities are two singletons (3), a singleton
and TT2 (5), a singleton and TT3 (9), a singleton and C3 (15), or two
TT2 factors (11). No case gives 13. This is not an unrestricted
nonexistence result; 13 has other tournament realizations.

The standalone verifier enumerates all labelled tournaments of orders
1 through 6, retaining physical (H,F2) pairs. The proof cap makes its
computed scalar image through 128 exhaustive for arbitrary factors, not
merely for a fixed seed menu. Literal tests use every ordered pair of
labelled factors of orders 1 through 3; subset path DP checks the full
polynomial law at j=1,2,3 and the entire rooted state. A separate inverse
implementation is compared with direct pair-image enumeration, including
integer vectors without a physical preimage. Exact universes and counts
appear in CHECK.json. The receipt does not prove all-order quantifiers;
those come from the proofs above.

Reproduce:

    python3 -I -B verify.py > replay.json
    cmp replay.json CHECK.json
    python3 -I -B -O verify.py > replay.optimized.json
    cmp replay.optimized.json CHECK.json

## 6. What this changes and what remains open

The new object is a finite **physical** inverse for a genuine bilinear
composition, with explicit witnesses and rejection certificates inside
its scope. It separates arithmetic factor choice from realization and
makes both finite at a target. The full colour law lets recursive outputs
become new gadgets without truncating away the coordinates needed by
later composition. It does not furnish an all-odd addition law: 13
already requires a seed or operation outside this cyclic class.

An adaptive all-odd proof would still need a family of operations and
physically realizable inverse branches covering every admissible target.
The divisor inverse gives one honest branch, not that covering theorem.
