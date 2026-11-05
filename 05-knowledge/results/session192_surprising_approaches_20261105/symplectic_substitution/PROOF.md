# Heisenberg patch substitution in the affine symplectic rank-four carrier

**PROVED construction and closure; FINITE-EXACT certificates. Independently reviewed.** Session 192 follow-up, 2026-11-05. This supplies geometric
local surgery outside every bounded-apex language and additional explicitly
constructed counts. It does not prove the full `{7,21}` spectrum, a scalar
ternary recursion, or a new literature-priority claim.

## 1. Inheritance and the unexpected source

The starting point is [the symplectic compiler](../symplectic/PROOF.md): every two-apex word has
an integral representation by points `(1,p,q,r)`, with

    omega(u,v) = u1*v2-u2*v1+u3*v4-u4*v3,
    u -> v iff omega(u,v)>0.

The previous closure was ordered join. The neglected source is
**THM-3243, `contact-deformation-blowup-equivariance-and-full-orbit-resolution`**,
specifically the distinction between moving a marked centre and preserving
the marked geometry under a Heisenberg action. Its finite-field blowup,
Singer orbit counts and LRC interpretations do not transfer here. The useful
idea is to retain the central compensation when placing a local configuration
at a new centre. The following real integer construction is proved directly.

The path-count sidecar comes from **THM-3121,
`path-cover-walk-content-substitution-kernel`**, and **THM-3134,
`tournament-endpoint-jet-and-c3-newton-profile-transform`**. Their path-count
formulas are inherited, not new. The new conclusion is that their arbitrary
recursive substitutions can remain in the same four-dimensional affine
symplectic carrier whenever all inputs are in that carrier.

The corrected near misses are `session191-grammar/END_COMMUTATION.md`,
`SPARSITY.md`, and `INDUCED_CORE.md`. Their obstructions remain valid.
Live board: central compensation; unequal geometric scales; block-run
profiles; feedback deletion; affine chart; literal versus gauge commutation.

## 2. Exact integer local surgery

Write `v=(p,q,r)` for the homogeneous point `(1,p,q,r)`. Define

    L_v(P,Q,R) = (p+P+rQ-qR, q+Q, r+R),
    D_M(p,q,r) = (M^2 p, Mq, Mr).

Direct expansion gives

    omega(L_v(u),L_v(w)) = omega(u,w),
    omega(D_M(u),D_M(w)) = M^2 omega(u,w).                 (1)

Thus these are Heisenberg translations and positive conformal symplectic
dilations on the affine chart. The central coordinate has weight two, not
weight one. The translation law is

    (p,q,r)*(P,Q,R) = (p+P+rQ-qR, q+Q, r+R).

Its commutator has central coordinate `2(rQ-qR)`. This noncommutativity is a
coordinate fact; a common translation still preserves every tournament arc.

**Theorem 1 (substitution closure).** Let a quotient tournament `Q` and
nonempty factor tournaments `S_a` each have integral representations in the
chart `(1,p,q,r)` with no zero off-diagonal pairing. Then `Q[S_a]` has such a
representation in the same four-dimensional space. Every internal pairing
can be retained exactly; every cross-block arc is retained from `Q`.

**Explicit compiler.** Let the outer points be `v_a=(1,p_a,q_a,r_a)` and
the points of block `a` be `w_ai=(1,P_ai,Q_ai,R_ai)`. For each `a<b` and
each `i,j` in the respective blocks, put

    A_ab = omega(v_a,v_b),
    B_aibj = (r_b-r_a)(Q_ai+Q_bj)+(q_a-q_b)(R_ai+R_bj),
    C_aibj = omega(w_ai,w_bj).

Let `D=max |B_aibj|`, `E=max |C_aibj|`, taking both maxima as zero if
there is only one outer vertex, and choose

    M = D+E+1.

The output points, with their labelled block membership, are

    z_ai = (1,
            M^2 p_a + M(r_a Q_ai-q_a R_ai) + P_ai,
            M q_a + Q_ai,
            M r_a + R_ai).                              (2)

**Proof.** Equation (1) gives, within one block,

    omega(z_ai,z_aj) = omega(w_ai,w_aj).                 (3)

Between distinct blocks, direct expansion of (2) gives

    omega(z_ai,z_bj) = M^2 A_ab + M B_aibj + C_aibj.     (4)

The integer `A_ab` has absolute value at least one. Moreover

    M^2-MD-E = M(E+1)-E >= 1.

Consequently (4) has the sign of `A_ab` and absolute value at least one.
Equations (3)--(4) are precisely all the labelled arcs of substitution. All
coordinates are integers, the first coordinate stays one, and the alternating
pairing matrix is `V^T J V`, so its rank is at most four. QED.

This is an actual local replacement operation, not merely a count identity.
Before clearing denominators, the equivalent small-patch map is

    w_ai -> L_(v_a)(epsilon^2 P_ai, epsilon Q_ai, epsilon R_ai).

It places each whole factor in an arbitrarily small neighbourhood of its
outer point, retains the factor pairing up to the positive scale
`epsilon^2`, and preserves all cross arcs for sufficiently small positive
epsilon. Formula (2) is its integral version after a common dilation.
Numerical coordinates of old outer points need not remain frozen; their
labelled tournament relation does. Selecting one point from each block gives
an induced copy of the quotient.

The common affine chart is a hypothesis. No claim is made for arbitrary
rank-four alternating sign realizations whose vectors cannot be put in one
positive affine chart. No finite-dimensional rank-two oriented-matroid
chirotope interpretation is being asserted.

## 3. Exact H response and the necessary sidecar

For a tournament `S`, let `pc_S(k)` count spanning covers by `k` unordered,
nonempty directed paths, and put `F_S(k)=k! pc_S(k)`. Let `W_Q(k_1,...,k_q)`
count directed quotient words in which letter `a` occurs exactly `k_a`
times; adjacent letters are distinct. Then

    H(Q[S_a]) = sum_(1<=k_a<=|S_a|)
                  W_Q(k_1,...,k_q) product_a F_(S_a)(k_a).       (5)

This is the `d=1` law of THM-3121. To see it directly, split a Hamilton path
at every change of block. The runs in each block form a spanning path cover.
The quotient word fixes the ordered occurrences of that block, and `F_S(k)`
counts the assignments of its actual path components to those occurrences.
Conversely those choices concatenate uniquely, proving (5).

For the directed triangle, a quotient walk follows the cyclic order, so the
three occurrence counts differ by at most one. For equal blocks this gives

    H(C3[S,S,S])
       = 3 sum_(r>=1) (F_S(r)^3
                     +F_S(r+1) F_S(r)^2
                     +F_S(r+1)^2 F_S(r)),                    (6)

where `F_S(r)=0` outside `1..|S|`. With two singleton factors,

    H(C3[S,K1,K1]) = 3H(S)+F_S(2).                          (7)

Equation (7) is not a scalar `3H+b` operation with a fixed b. The singleton
and the transitive two-vertex tournament both have H=1, but their F(2)
values are 0 and 2, producing 3 and 5 respectively. Thus even this elementary
geometric patch requires a sidecar that H alone does not determine.

## 4. A rank-four family outside every bounded-apex language

Set `T_0=K1`, `T_(d+1)=C3[T_d,T_d,T_d]`. Start with the triangle points

    (1,0,1,0), (1,0,0,1), (1,0,-1,-1).

Theorem 1 recursively gives an integral rank-at-most-four representation of
every `T_d`. All `T_d` for `d>=1` are strong.

**Theorem 2.** The largest transitive induced subtournament of `T_d` has
exactly `2^d` vertices. Its feedback vertex number is therefore

    tau(T_d)=3^d-2^d.                                     (8)

For `d>=2`, its minimum possible alternating sign-realization rank is
exactly four.

**Proof.** A transitive subset cannot meet all three outer blocks, because
one point from each is a directed triangle. It meets at most two blocks, with
at most `2^(d-1)` points in each by induction. Equality is attained by taking
maximum transitive subsets in two blocks, which form an ordered join. This
proves (8). Rank at most four follows from Theorem 1. For `d>=2`, a vertex
dominates a block containing a triangle. Hence its out-neighbourhood is not
transitive. Every planar determinant tournament has transitive out- and
in-neighbourhoods, since those vectors lie in an open half-plane and are
linearly ordered by angle. Thus rank two is impossible. Alternating ranks
are even, so the minimum is four. QED.

In particular `T_2` has nine vertices, H=3159, and tau=5, so it cannot be
represented by a two-apex word under any relabelling. The compiler finds
`D=4,E=1,M=6` and the following explicit nine points:

    (1, 0, 7, 0), (1,-6, 6, 1), (1, 6, 5,-1),
    (1, 6, 1, 6), (1, 0, 0, 7), (1,-6,-1, 5),
    (1,-6,-5,-6), (1, 6,-6,-5), (1, 0,-7,-7).

Their outgoing bit rows are

    (58,60,57,464,480,456,135,263,71).

The verifier checks every arc, H by subset DP and all 9! permutations,
tau=5 by exhaustive subsets, and alternating rank four by exact elimination.
The full ordered path-cover profile, independently counted from vertex
subsets and the component containing the least name, is

    F_(T2)(1..9) = (3159,42606,253206,874800,1934280,
                   2812320,2630880,1451520,362880).

For `T_3`, the next compiler has `D=28,E=49,M=78`. All 27 points and 351
physical pairings are checked in `receipt.json`. Formula (6), with the
independently counted preceding profile, gives the inherited target

    N_star = H(T_3) = 417382500592116023859
                    = 3^13 * 1309949 * 199849717.          (9)

The two latter factors are prime, verified by trial division. Thus (9) is
outside the original Session 192 theorem's 139-smooth domain. The new result
is its explicit rank-four geometry, not the previously known integer count.
The rank is exactly four and (8) gives tau=19. No 27-vertex subset-DP claim
is made: its H is certified using the proved substitution law.

Combining this construction with the earlier ordered-join theorem gives
another explicit infinite family:

    N_star^k * N,  k>=0,

where N is any of the previously constructed positive odd 139-smooth
integers other than 7 and 21. This is an extension of that sufficient family,
not permission to extract the two prime factors in (9) separately.

## 5. Hostiles, commutation, and the remaining spectrum boundary

**Ordinary translation is not enough.** Translate the triangle's `(p,q,r)`
coordinates by `(0,1,0)` without the central correction. The three upper
pairings become `(2,-2,-1)` instead of `(1,-1,1)`: a directed triangle turns
transitive, so H changes from 3 to 1. The Heisenberg correction restores all
three pairings exactly. The central correction cancels first-order centre
effects; the quadratic scale keeps central differences and planar areas at
the same order, preserving the desired internal tournament.

**Noncommuting gauge transformations are not a ternary grammar.** Common
symplectic transformations, Heisenberg translations and positive dilations
preserve all arcs, hence H. Their group commutator cannot create a nonzero
Hamilton-count increment. Likewise substitutions at two fixed, disjoint
module addresses commute as physical tournament operations. They cannot
satisfy incompatible universal affine laws `3H+2` and `3H+4`: Session 191's
commutation obstruction applies unchanged. Changing a later active address
or nesting module substitutions can be order-sensitive, but no suitable
scalar digit laws are proved here.

**The fixed-template sparsity theorem is not contradicted.** Recursive
substitution by newly constructed nontransitive blocks is expressly outside
`SPARSITY.md`'s finitely many fixed templates with boundedly many transitive
clocks. The tower (8) demonstrates that the geometric class really exits
bounded feedback-vertex languages. A single tower is still a sparse target
sequence; exiting a prior hypothesis does not prove density or completeness.

**The induced-core obstruction survives.** The operation retains an induced
copy of each factor and of the quotient, so it cannot yield a strong H=39
tournament containing an H=13 core. The existence of rank-four coordinates
does not relax that purely combinatorial prohibition from `INDUCED_CORE.md`.

The remaining concrete task is to design nested/address-dependent patches
whose complete path-cover sidecar supports inverse target selection. Geometry
is now closed under those patches; selecting arbitrary admissible odd H
values is still open.

## 6. Verification and connection contract

Run from the workspace root:

    python3 verify.py > normal.stdout
    python3 -O verify.py > optimized.stdout
    cmp normal.stdout optimized.stdout

The standard-library verifier has no optimization-removable assertions. It
checks the explicit nine- and 27-vertex constructions; 101 further labelled
geometric substitutions (all two-apex words of length at most three with
one triangle patch, plus all 16 ordered pairs of one-letter outer/inner
words); 21 small outputs by two independent path algorithms; 21 instances
of (7); the same-H hostile; and the ordinary-translation hostile. The
finite universe is explicitly a control bank, not all rank-four tournaments.

Source: affine symplectic point configurations and complete block-run
profiles. Target: the labelled tournament substitution and its exact H.
Map: centre-dependent Heisenberg translation after unequal coordinate
dilation. Preserved: every internal arc, every quotient cross arc, and the
common affine chart, with integer unit margins. Lost by a scalar quotient:
all higher path-cover counts. Sidecar: `F_S(k)`, retained in (5). Cheapest
decisive tests: the translated triangle with H=1 versus3, and same-H singleton
versus transitive edge with wrapper counts3 versus5. No LRC physical state,
finite-field orbit interpretation, or full-spectrum implication transfers.
