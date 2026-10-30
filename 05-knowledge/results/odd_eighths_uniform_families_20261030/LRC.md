# Uniform tail banks with directional owners retained

**PROVED scoped**, same-agent proof review. L1--L3 are local labels in this
file. LRC(14) and arbitrary-core classification remain OPEN.

Write ||x|| for distance to the nearest integer, S_V(r) for the closed
survivor set {t mod1: ||vt||>=r for all v in V}, and M(V)=max_t min_v||vt||.
All speed sets are finite, nonempty, positive and integral unless stated.

## Inherited mechanism and exact scope

[Session114 GENERAL_CORE.md](../session114_continuation_symmetry_20261030/lrc/GENERAL_CORE.md)
proves the complete top chamber for

    C_q={1,...,q-1,q+1}, q>=3,
    rho=1/q, r0=1/(q+2), s=s(q)=min{a>=2:gcd(a,q)=1}.

For r0<r<=rho, write epsilon=rho-r. The surviving components are

    [a/q-epsilon/p_a, a/q+epsilon/n_a],  gcd(a,q)=1,

where p_a=a^-1 mod q in {1,...,q-1}, n_a=q-p_a, and every occurrence
of owner1 is replaced by q+1. In particular the minimum directional owner
is s, occurring on two reflected sides. This is stronger than a list of
maximizers: its proof supplies the entire radius chamber. We inherit it.
The same source gives the single-tail maximum and its two strict-branch
maximizers. The extension here treats a whole simultaneous tail bank.

## L1. A first-arrival test for every finite resonant bank

For a finite multiplier bank T, define its first arrival at radius r by

    tau_T(r)=min{x>=0: ||mx||>=r for every m in T},

with infinity when the set is empty. The minimum may be taken in [0,1].
For every q>=3, every such T, and every r0<r<=1/q,

    S_(C_q union qT)(r) is nonempty
       iff  q*r+s(q)*tau_T(r)<=1.                         (1)

More precisely, on a side of a/q with owner e in {p_a,n_a}, its entire
surviving part is parametrized by

    t=a/q +/- x/q,
    0<=x<=(1-q*r)/e,   x in S_T(r) as a periodic subset of R. (2)

Proof. A tail qm has integer phase at a/q, so its distance at the displaced
point is ||mx||. The inherited core window has exactly the displacement
bound in (2). At least one side has owner s, and every other owner is at
least s. Existence on any side is therefore equivalent to the first arrival
fitting in the largest available interval. This is precisely (1).

The same argument works for ANY core whose complete top windows are
[c_j-(rho-r)/p_j,c_j+(rho-r)/n_j], with Qc_j integral, and a tail bank QT:

    min_j min(p_j,n_j) * tau_T(r) <= Q*(rho-r).             (3)

This conditional general theorem requires the full top-window hypothesis.
It does not infer windows from symmetry, nor classify radii below the chamber.
All joint tails enter the SAME tau_T. Separate tail arrival times cannot be
minimized independently and then assumed to coexist.

## L2. Complete classification of maximum-preserving finite attachments

Fix q>=3 and m>=1 with 2m>s(q). Put

    D=qm+s, a=s^-1 mod q in {1,...,q-1}, h=(as-1)/q,
    K=am+h, r=m/D, V0=C_q union {qm}.

The inherited strict-branch theorem gives

    M(V0)=r, argmax(V0)={K/D,1-K/D}.                       (4)

For ANY finite set Z of additional positive speeds, including nonresonant
ones, the following is an exact iff classification:

    M(V0 union Z)=r
       iff  m <= (vK mod D) <= D-m for EVERY v in Z.       (5)

When (5) holds, the two points in (4) remain all the maximizers. If any
speed fails (5), the new maximum is strictly smaller than r.

Proof. Adding speeds cannot increase M. Equality requires a maximizer from
the two-point list (4). Reflection gives the same distance for every added
speed at both points, and their exact safety condition is (5). Compactness
of the circle makes failure strict, not just nonattainment of the supremum.
The identity qK=aD-1 shows gcd(K,D)=1. Consequently (5) admits exactly

    D-2m+1=(q-2)m+s+1

residue classes modulo D. This is a classification of ALL finite attachment
banks preserving this maximum, not a finite search over a height cutoff.

For resonant additions v=qz, use qK=-1 mod D to simplify (5) to

    m <= (z mod D) <= D-m.                                (6)

For example every finite T with min T=m and max T<=(q-1)m+s has
M(C_q union qT)=m/(qm+s). The residue version (6) includes unbounded
and widely separated multiplier sets, beyond that interval family.

## L3. Exact terminal mass and a sharp loss boundary

Suppose Z satisfies (5). Define the actual fastest binding owners

    P=max({s} union {v in Z:v=s mod D}),
    N=max({qm} union {v in Z:v=qm mod D}).

Then for every sufficiently small epsilon>0,

    measure S_(V0 union Z)(r-epsilon)
        =2 epsilon (1/P+1/N).                            (7)

Proof. At K/D, sK=m mod D and qmK=-m mod D. Because K is a unit,
an added speed ties the positive wall exactly when v=s mod D, and the
negative wall exactly when v=qm mod D. Every other speed has strictly
positive slack. Among tied walls on one side the largest speed gives the
smallest physical expansion epsilon/v, so the local interval length is
epsilon(1/P+1/N). There are exactly two isolated maxima, and reflection
swaps the two owners. Compactness excludes additional components sufficiently
near r. This proves the exact germ (7).

Taking added speeds s+jD and qm+jD preserves the maximum and both maximizing
times while making the mass coefficient tend to zero. Thus neither the
maximum value nor its two-point symmetry orbit controls a positive uniform
terminal mass floor. The retained sidecar is the actual largest tied speed,
not just its residue or the original core owner. This is an infinite-family
obstruction to that proposed compression, not to LRC.

## The literal odd-eighth family

For q=8, C_8={1,2,3,4,5,6,7,9}. Its maximum times and owners are

| Time | Positive owner | Negative owner |
|---|---:|---:|
| 1/8 | 9 | 7 |
| 3/8 | 3 | 5 |
| 5/8 | 5 | 3 |
| 7/8 | 7 | 9 |

The strict branch begins at m=2. At that value,

    M(C_8 union {16})=2/19,
    argmax={7/19,12/19},
    admissible added speeds = all positive v NOT in {0,8,11} mod19.

Hence adjoining ANY four distinct new admissible speeds gives an infinite
thirteen-speed family with exact clearance2/19>1/14. Taking even one forbidden
speed loses2/19; it does not imply failure at1/14. For instance adding8
gives maximum1/10, independently checked by all opposite-slope clocks.
At m=3 and4 the exact maxima are1/9 and4/35; the checker records their
complete forbidden residue lists as further controls.

The number8 here is an exact denominator and the four primitive phases form
its unit orbit. The unequal directional owners show why that orbit cannot
be substituted for the full window geometry.
