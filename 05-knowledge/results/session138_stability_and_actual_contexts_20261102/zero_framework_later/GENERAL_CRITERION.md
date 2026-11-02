# Finite cyclic-response certificates for actual tournament modules

**Status: PROVED general sufficient framework; independent Session138
reconciliation proof and witness audits ACCEPT.**
Failure of a certificate is explicitly inconclusive. The criteria are not
claimed necessary for a sign on all realizable side modules.

## Typed carrier

Let M0,M1 be nonempty tournaments of the same order n. They may, but need not,
be two completions of one vertex-deletion card. Let

    d_r = f_(M1)(r)-f_(M0)(r),

where f_W(r) counts ordered covers of all vertices by r nonempty directed
paths. Extend each profile by zero beyond its order, and set d_0=0. Define
the finite symmetric tridiagonal matrix K on coordinates1,...,n+1 by

    K_(r,r)=d_(r-1)+3d_r+d_(r+1),
    K_(r,r+1)=K_(r+1,r)=d_r+d_(r+1),                   (1)

with all other entries zero. For every pair of nonempty side tournaments A,B,

    H(C3[M1,A,B])-H(C3[M0,A,B]) = f_A^T K f_B.         (2)

This is exact: contract each maximal module run. Outer C3 run-words are
cyclic sequences; for r occurrences of M there are three words at side
multiplicities(r,r), and one word each at(r-1,r-1),(r+1,r+1), and the four
adjacent pairs(r-1,r),(r,r-1),(r+1,r),(r,r+1). Regrouping gives(1).
Coordinates beyond n+1 in a large side module do not contribute. This loss
of high side coordinates is exact, not an approximation.

## Certificate P: universal positivity

If every entry of K is nonnegative, then the response in(2) is nonnegative
for every nonempty A,B. If additionally K_11>0, then

    response >= K_11 H(A)H(B)>0.                       (3)

This is an exact finite sufficient test on the actual module profiles.
Within the entrywise-nonnegative class, K_11>0 is also necessary for strict
positivity for **all** nonempty A,B: singleton sides isolate that entry.
Entrywise nonnegativity itself is not claimed necessary over actual sides.
A negative matrix entry therefore means `NOT_CERTIFIED`, not a realizable
negative context.

## Certificate N: a negative polynomial-multiplier chamber

Choose any nonempty actual base tournament V, with order N and profile
`v=(v_1,...,v_N)`. Suppose the following finite conditions hold:

1. v is log-concave: `v_r^2>=v_(r-1)v_(r+1)` for interior coordinates.
   All its coordinates are positive, as for any tournament profile.
2. There is a pivot c with1<=c<=N such that every entry in row i<c of K is
   nonnegative and every entry in row i>c is nonpositive. The pivot row
   may have either sign or mixed signs. Pad K by zero rows as necessary.
3. The exact quadratic form satisfies `v^T K v<0`.

Write the colored path-cover polynomial

    P_W(m)=sum_r f_W(r) binom(m,r).

Then the response in(2) is strictly negative for **every actual pair** A,B
satisfying

    P_A=P_V R_A,  P_B=P_V R_B,                          (4)

where R_A,R_B are nonzero polynomials with nonnegative monomial coefficients.
The quantifier is over actual side tournaments in this polynomial chamber;
the theorem does not assert that every formal product in(4) is realizable.

### Why the chamber is nonempty and structurally useful

The colored odd-cycle-family expansion gives P_W nonnegative monomial
coefficients for every tournament W. Indeed color each ordinary odd-cycle
family: a family C contributes `2^|C| m^(|V(W)|-sum_C(|cycle|-1))`.
Ordinal joins multiply P. Thus(4) holds for A=V=>W_A and B=V=>W_B, with
arbitrary tournaments W_A,W_B, including empty ones. It also holds when
these factors appear in any ordinal order.

For an eventual-power family, fix a tournament F and a proposed integer k.
Take `V=TT_k[F,...,F]`, whose color polynomial is P_F^k, and check the three
finite conditions above. An accepted certificate then proves negativity
for all side multiplicities a,b>=k in `TT_a[F]`, `TT_b[F]`, and for further
nonnegative polynomial factors. It is not asserted that some k must pass;
the procedure certifies a supplied base, rather than performing an unbounded
search or deciding every eventual sign problem.

### Proof of Certificate N

The following multiplier lemma is also used in the preceding
`zero_balance_later/NO_ZERO_THEOREM.md`; its proof is included here. If v has positive
contiguous support and is log-concave, multiplication by m in the binomial
basis sends

    v_r -> (Lv)_r=r(v_r+v_(r-1)).

For convolution by(1,1), its adjacent log-concavity difference is

    (v_r²-v_(r-1)v_(r+1))
    +(v_(r-1)²-v_(r-2)v_r)
    +(v_r v_(r-1)-v_(r-2)v_(r+1)),

which is nonnegative by log-concavity and monotone adjacent ratios. Zero
boundary coordinates are harmless. Multiplication by the positive
log-concave sequence r also preserves log-concavity. Moreover `(Lv)_r/v_r=r(1+v_(r-1)/v_r)` is
nondecreasing. Iteration and nonnegative linear combination prove that
multiplication of P_V by any nonzero monomial-nonnegative polynomial yields
positive ratios `f_A(r)/v_r` nondecreasing on1,...,N, with nonnegative tail.
No assumption that the final sum remains log-concave is used.

For any nonnegative vector w, the row-sign condition makes Kw nonnegative
below c and nonpositive above c; its pivot coordinate can fall on either
side. Thus Kw has at most one change from positive to negative, and no
positive coordinate beyond N.

If a vector h has this one-crossing shape and `v^T h<0`, the increasing-ratio
property implies `f_A^T h<0`: choose a positive scale at the last positive
coordinate, bound positive terms from above and negative terms from below
in magnitude; any tail contributes nonpositively. If h has no positive
coordinates the conclusion is immediate because the original support
remains positive.

First use h=Kv and condition3 to obtain `f_A^T K v<0`. Next use
h=K f_A; symmetry gives `v^T h<0`, and the same one-crossing argument gives
`f_B^T K f_A<0`. This proves Certificate N.

## Algorithm and truth boundaries

`certify.py` exposes the exact tests. `certify_actual(M0,M1,V)` validates the
literal tournament masks, computes complete ordered profiles by backward-H
and anchored-partition DP, and then runs both certificates. The lower-level
functions accept profile data but certify only their algebraic predicates;
realizability must not be inferred from a formal vector. The certificate
contains the pivot, base quadratic, and scope; failed tests return
`NOT_CERTIFIED` rather than `FALSE`.

Computing the profiles can be expensive; once they are supplied, checking K,
log-concavity, a proposed pivot and the quadratic form is finite exact
integer arithmetic. No numerical eigenvalue or sampled parameter test is
load-bearing. Broader one-crossing patterns, bases without log-concavity,
and arbitrary side color polynomials remain outside the certificate.

## Minimal actual hostile to omitting the row condition

Let the common six-vertex card have masks

    B=(32,1,3,7,15,30).

Add its seventh vertex at cuts5 and14, respectively:

    M0=(32,65,3,71,79,94,5),
    M1=(96,1,3,7,79,94,14).

Both modules are strong. Their exact ordered profiles are

    f_M0=(51,678,3798,11280,18240,15120,5040),
    f_M1=(39,678,4038,11760,18480,15120,5040),
    d   =(-12,0,240,480,240,0,0).                       (5)

Take actual base V=K1, with log-concave profile(1). Then

    v^T K v = -36 < 0,
    H(C3[M1,TT2,TT2])-H(C3[M0,TT2,TT2])=828>0.           (6)

The new side polynomial is `P_TT2=m^2=P_V(m)m`, with a nonnegative
monomial multiplier. Triangle sides give another actual positive control,
with response67932 and multiplier `m^2+2`. Thus all multiplier/base conditions hold, but the
negative conclusion fails when the row condition is removed. In fact
K_22=228>0 while the only permissible pivot for V is c=1, so the forbidden
positive high row is explicit. This is an actual same-card counterexample,
not a formal network or a conjectural profile.

**Minimal module order7, FINITE-EXACT certificate through order6.** The
independent literal program `minimality_literal.cpp` enumerates every
labelled tournament of orders1,...,6 and all its vertex permutations:
33,867 tournaments,23,717,429 words. At each order it recovers all distinct
ordered profiles and checks every pair's cyclic K. The counts are

    orders:             1,2,3,4,5,6
    distinct profiles:  1,1,2,3,8,24
    profile pairs:      0,0,1,3,28,276.

Every one of these308 matrices is entrywise one-signed. If its base
quadratic is negative, K is nonpositive and a negative entry already meets
the positive base support; multiplying by a nonzero nonnegative polynomial
retains positive support on those coordinates. Thus no negative-to-positive
propagation hostile exists at smaller module order, for any actual base or
admissible multipliers. This proves minimality, not merely the absence of
one selected small context. All39 small profiles also pass log-concavity.

## Log-concavity boundary

Dropping log-concavity from the multiplier lemma is unsafe even for a short
formal profile: v=(2,1,2) gives Lv=(2,6,9,8), whose ratios on the original
support are(1,6,9/2), not nondecreasing. This is explicitly a formal control,
not an asserted tournament profile. No actual non-log-concave tournament
profile is claimed; none occurs through order6 in the complete literal
control above. The actual hostile(5),(6) keeps log-concavity and isolates
precisely the missing row-crossing condition.

## Reproduction notes

`framework_audit.py` independently counts the seven-vertex hostile profiles
both by anchored partition DP and by all5,040 permutations per module. It
then counts all Hamiltonian paths of the six literal outer contexts by a
separate forward-extension subset DP, confirming absolute counts831,795
(singleton sides),27,279,28,107 (TT2 sides), and1,408,455,1,476,387
(triangle sides). It also exercises
the reusable certifier on the previously accepted actual Q4 controls and
checks that the hostile is correctly reported `NOT_CERTIFIED`. Normal and
optimized runs agree.

The failed order-six discovery precursor is retired under `exploratory/`
and is not part of the public manifest. The accepted hostile is fully
specified by `actual_row_hostile.json` and reproduced by `framework_audit.py`;
its minimality is separately reproduced by `minimality_literal.cpp`.
