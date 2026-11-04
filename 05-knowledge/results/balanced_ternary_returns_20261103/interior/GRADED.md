# Positive rooted defects and a bridge from balanced sources

**PROVED scoped attachment identities and rigidity; FINITE-EXACT named
controls.** The H35 source that supports both ternary branches suggests
allowing the rooted defect to increase. This produces an exact preparation
chart and an actual lift of every balanced source, but neither is a closed
ternary induction. The free parameter needed for closure is constrained by
an additional source-realizability condition.

## The inherited prefix law retains three exterior quantities

Let W be an actual tournament with two ordered ports a->b and J={a,b}.
Write u=H(W), let E=E_J,W(2) count ordered spanning two-path covers whose
first nonempty path ends in J, and let C count covers whose first path ends
in J and second nonempty path starts in J. Empty second paths are allowed
in E. Prepend a transitive k-chain to W, k>=1, then add R which loses to
a,b and beats all other vertices. The inherited
[ATLAS equation A4](../../ready_return_pentagon_20261103/ATLAS.md) gives

    H = 2^k E + u - C,
    s_R=q_R=u,       F_R=H+C.                    (G1)

These identities hold for arbitrary actual W, not just the transitive
residual examples in that atlas. Removing R splits a Hamilton path into
two paths; each prefix vertex can be assigned to either path, giving 2^k E.
The forbidden second-path starts in J require all prefix vertices in the
first path, and subtract C. Paths starting at R contribute u. Deleting R
gives the ordinal sum of the prefix and W, also of count u. The two-in-port
weighted-cut formula then gives F_R-H=C. This proves the extension of G1
directly; it does not infer arbitrary W from a fitted coefficient pair.

If W= a -> P is a new global source a followed by a tournament P with a
named vertex b, write its rooted data as (h,s,q,F). Then

    u=h,       E=h+2F,       C=s+F.              (G2)

For E_a the first path must be the singleton a, leaving h choices. For
E_b, a may prefix either of the two paths of a cover counted by F, giving
2F. In C_ab the first path is (a), leaving s choices; in C_ba removal of
the initial a from the second path gives precisely F choices. These
bijections retain the ordered ports and all old vertices.

## The full two-parameter H35 preparation chart

Take P=b->TT_r->v->b, r>=1. Its rooted profile at b is
(u,1,1,u), u=2^r+1. Equations G1 and G2 give

    H = 3(2^r+1)2^k-1,
    q_R=s_R=2^r+1,
    F_R-H=2^r+2,                       k,r>=1.   (G3)

In the standard two-apex notation this actual tournament is
`3 0^k 3^r 1`, with R the first chain vertex and a,b the apices. In the
READY in-port notation its residual word is `11^k 00^r 01`. This makes
the dictionary to ATLAS explicit instead of treating similar words as
the same representation.

Changing r=1 to r=3 by two interior 3 births retains the physical ports
and root and gives

    H -> 3H+2,       s -> 3s,       (F-H) -> 3(F-H)-2.

At k=2 the named path is H35 -> H107 with profiles
(35,3,3,39) -> (107,9,9,117). The growing defect is real attachment data.
The numerical family and fixed-prefix mechanism are inherited; the
receipt retains the actual old-name embedding from the conditional lane.

This source is already an actual intermediate of the preceding balanced
construction. For m=2, its B_2 word2032 followed by only the first NR birth
has H35 and old physical x profile (35,3,3,39). Re-marking at (c_1,c_2)
gives word30031, with old x as its first chain root. The second NR birth
would instead seal it to the balanced H59 source. Thus keeping the
intermediate's positive defect is a literal alternative to that sealing
choice. This identification is specific to m=2: at m=3 the old-x profile
is (71,3,3,87), whereas the graded k=3,r=1 source has (71,3,3,75).
Equal H,s,q do not identify their attachment data.

This chart has an exact iteration barrier. Any ternary plus-two step
between two parameter pairs would require

    2^K(2^R+1)=3*2^k(2^r+1).

Oddness forces K=k, and then 2^R=3*2^r+2. The only r,R>=1 solution is
r=1,R=3. Every chart value is -1 modulo 3, so no plus-four branch can
land anywhere in the same chart. Allowing the defect to grow therefore
does not by itself close the source family.

The [conditional lane's H109 sibling](../conditional/NOTE.md) has profile
(109,5,5,121) and exterior data (u,E,C)=(5,29,12). Its prefix family is
29*2^k-7, which equals 3(9*2^k-1)+4 only at k=2. The successful sibling
changes the exterior attachment data, and its first numerical success
cannot be promoted to a uniform branch rule.

## Assuming a free exterior parameter reveals a rigidity obstruction

A tempting hypothesis is to replace u=2^r+1 by a freely selected attained
odd u while retaining the exterior requirements

    H(W)=u,       E=3u,       C=u+1.              (G4)

That is an actual joint-realization hypothesis, not an algebraic definition
of a tournament. In the factorization W=a->P with P strong it is rigid.
The singleton P gives u=1. For P of order at least three, G2 and G4 force
F=h and s=1. The inherited nonnegative decomposition

    F-h = (q-s)+delta,
    q>=s,       delta>=d(d-1),

where d is the indegree of b in P, forces q=s=1, delta=0 and d=1.
Thus P-b is transitive. Strongness forces the sole in-neighbour v of b
to be last in that transitive order, so P is exactly b->TT_r->v->b.
The proposed free u has collapsed back to 2^r+1. Releasing the global-source
factorization or changing E and C is necessary for this approach to gain
new values of u.

The singleton boundary supplies an actual first step before G3:
`3 0^k -> 3 0^k 3 1 -> 3 0^k 333 1`. These are two successive plus-two
steps with u=1,3,9 and positive defects 2,4,10. At k=2 they give
11->35->107. The next required u=27 is excluded by the rigidity theorem
for this factorization. The first step recovers the singleton preparation
mechanism of READY section 6; it is not a new numerical family.

Nor can G4 be assumed for all odd u: u=7 and 21 are excluded by
[THM-1370, the all-order exclusion](../../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md).
Even hypothetical coverage of every allowed u in G4 leaves targets such
as the attained primes 41 and 167 outside this chart: their unique
decompositions H+1=3*u*2^k require u=7. Actual word witnesses are retained
as hostile controls. A covering proof therefore needs several compatible
charts, in addition to a source-realization theorem within each chart.

## An actual bridge from every balanced source

Let P be any actual strong tournament with named root b and balanced
profile (B,t,t,B). Apply G2 and G1. The resulting strong tournament
retains every old vertex and arc, with fresh root R, and has

    H=3B*2^k-t,
    q_R=s_R=B,
    F_R-H=B+t,                         k>=1.     (G5)

This is a constructive specialization of the inherited prefix kernel.
The receipt directly checks the old C2 pair (59,17), whose k=1 lift has
profile (337,59,59,413), and two further named literal graphs. Thus the
positive-defect strategy can use the earlier balanced sources without
assuming that arbitrary exterior signatures exist.

The boundary k=0 is excluded. Then a is a global source of the whole
output, and the remaining graph is the ordinary native birth at b in P.
Its H is 2B-t, not 3B-t, and the new root has s_R=0 and q_R=B. For the
pair (59,17) the literal profile is (101,0,59,236), whereas illegitimate
substitution in G5 would give H=160. In particular k=0 cannot be used
below to manufacture a formal carry of four.

For a change of actual balanced source (B,t)->(B',t') at fixed k, G5
gives the exact ternary carry condition

    3*2^k(B'-3B) - (t'-3t) = c,       c in {2,4}. (G6)

The already-proved interior map (3B+2,3t+2) instead gives carry
6*2^k-2, at least 10 for k>=1. It therefore does not meet G6. This is the
specific missing compensation when transferring that construction into
the graded chart.

The bridge has a useful 7-adic control. When B is a unit modulo 7,
the lift is divisible by 7 exactly when t/B lies in
3<2>={3,5,6} modulo 7 and k selects the corresponding phase. A primitive
pair with 7|B cannot give a multiple of 7 through this lift. The original
pairs (15*2^m-1,3*2^m+5) never land in the 7p lane by G5. Their interior
plus-two children do: for m=0 mod3 and k=1 mod3 the lift is divisible by7.
At k=1,m=3j the actual counts and quotients are

    H_j=261*8^j-23=7q_j,
    q_(j+1)=8q_j+23,       j>=1.

The first quotient is the hostile composite 295=5*59. The next is the
verified prime 2383, giving a literal 17-vertex H16681 source. Furthermore
49 divides H_j exactly when j=4 mod7. These are explicit source and phase
controls, not an infinitude or spectrum-coverage claim.

Run `python3 -B graded.py` and `python3 -B -O graded.py` in this directory.
Both outputs must equal graded.json. The all-size formulas rest on the
displayed bijections; the finite receipt retains literal graphs, physical
maps, root profiles, positive and hostile controls, and source hashes.
