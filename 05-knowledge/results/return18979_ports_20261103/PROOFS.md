# Grandparent defects, port completion and an exact ternary carry boundary

**PROVED scoped elementary mechanisms; FINITE-EXACT supplied-source certificates.**
Same-agent proof review, with independent counting paths described in AUDIT.
Local labels RC1--RC5 are scoped to this file. Full source closure and the
Hamilton spectrum remain OPEN. A return of the actual18979 source has not
been constructed in this packet.

All tournaments have actual Boolean orientations. H counts directed spanning
paths, H(empty)=1. F_r counts ordered spanning covers by a nonempty path ending
at r and a second possibly empty path. A root is b-ready when F_r-H=b.
The inherited two-birth reset produces H+2F_r=3H+2b.

## RC1. The complete ready locus of transitive root growth

Let T be a strong tournament, r a named vertex, and replace r by an ascending
transitive block r1,...,rm, with m>=2. Among all marks in this block:

* a1-ready mark exists iff r has a unique in-neighbour w and T-{r,w} is
  transitive; in that event it is precisely the first clone;
* no mark is2-ready.

**Proof.** Put I=N-(r), d=|I|. Deleting the first clone leaves a strong
tournament C, namely T with r replaced by m-1 clones. Deleting the initial
r1 from a spanning path shows q-s=sum_(i in I) start_i(C)>=d. The inherited
rooted defect identity gives

    F_r1-H = q-s + delta,       delta >= d(d-1).

See [READY equation7](../two_entry_funnels_20261103/READY.md) and its exact
weighted-cut dependency. Hence the first defect is at least d^2. For d>=2
it is at least4. For d=1, every C-path starting at w must immediately traverse
all m-1 remaining clones in order: their only external incoming vertex is w,
which cannot be reused. The remainder is an arbitrary path of T-{r,w}. Thus

    F_r1-H = H(T-{r,w}),                       m>=2.       (RC1a)

It is1 iff the residual tournament is transitive, and is otherwise at least3. The latter
equivalence follows from H=1 iff the tournament has no directed triangle:
a triangle already has three paths and insertion of each other vertex
preserves distinct restrictions. Conversely a transitive tournament has
one path.

For adjacent twins u->v, the exact difference F_v-F_u is the positive F
of their contraction, as proved in [O4](../return6325_observer_20261103/PROOFS.md).
The contraction here is strong, so this increment is at least H>=3.
Consequently later clones cannot have defect1 or2. QED.

For the single-clone case, C=T-r and w its unique in-port, direct deletion
gives s_r(T)=H(C)-s_w(C), H(T)=s_r(T)+F_w(C), and
F_r(T)=H(C)+F_w(C). These also prove F_r-H=s_w(C).
The concurrently published [CL1/CL2](../return6325_completion_20261103/CLONE_LAW.md)
contains this grandparent law for every parent, including weak boundaries.
Our derivation was made before reading that commit; the identity is credited
as overlapping work. The new extensions here are the arbitrary-indegree
clone classification and the concrete repair/word-carry analysis below.

## RC2. Which old ports must a repair touch?

Fix the vertex set, active root r and a set S of other vertices whose incident
arcs may be changed. Permit arbitrary changes of the root star too. A necessary
condition for a corrected strong tournament to have1-ready transitive root
growth is

    there exists w!=r such that T-( {r,w} union S ) is transitive.       (RC2a)

Indeed RC1 forces a unique final in-port w and an entirely transitive
grandparent. The displayed restriction is untouched by every allowed edit.
This is an endpoint-support obstruction; helpers or a different edit order
do not remove it if the final edits have the same support.

For n>=3, nonempty S, and a prospective w for which an active vertex remains
outside {r,w}, the condition is also sufficient. Keep the untouched chain,
insert active vertices to form a full transitive order on V-{r,w}, and place
one active vertex last. Set the root star to w->r->V-{r,w}, and set the last
grandparent vertex to beat w. Every vertex reaches w along this order, and
w reaches all vertices through r, proving strongness. If S={w}, the same
construction uses the editable w star to orient the last vertex into w.
These constructions preserve every forbidden-to-edit arc.

Consequently, if a nonempty old spare set S0 is already editable, the exact
number of additional editable old ports needed is

    min_(w!=r) tau( T-( {r,w} union S0 ) ),               (RC2b)

where tau is the minimum number of vertices meeting every directed triangle
(equivalently, making the tournament transitive upon deletion). Necessity
is RC2a; a minimum triangle transversal supplies the additional editable
vertices, and the preceding construction proves sufficiency. This has the
same set-cover form as finite observer completion, but the covered objects
are untouched cycles, not unresolved pairs of sources. No H-neutrality is
included in this minimum.

**Actual source.** Use the returned6325 source, old root31, spare6, in-port21.
Its grandparent H is27. If only root31 and spare6 can be edited, RC2a fails
for every prospective w. The receipt supplies an untouched directed triangle
for each of the twelve choices. Thus **arbitrary mixed root/spare corrections
cannot restore clone readiness**, even without demanding H-neutrality.

Exactly one additional old port is sufficient. Exhausting the eleven possible
extra labels gives precisely these two choices:

|Additional editable port|Only possible final in-port|
|---|---|
|2|21|
|21|2|

If the original root star is retained, only the first choice remains. This
is a repair-port completion number, not the query completion number of O1:
changing an arc and observing a response are different operation menus.

The second choice also has an explicit finite witness: reverse the five
unordered pairs {2,31},{5,21},{8,21},{14,21},{21,31}. The source is strong,
with profile (1837,778,779,1838), unique in-port2, and clone law
1058*2^m-279. Its first-clone defect is1 for every m>=2. Five is minimal
within this specified root/spare/extra-port support, as certified by all
orders inserting6 and21 into the fixed uneditable chain and all four
allowed incidence choices with2. The numerical normalization
A=1058 ->529 shifts m to m+1 under RC3; it retains the clone hierarchy.

In the fixed-star choice two reversals suffice, both through old vertex2:
reverse 2->32 and 4->2. The grandparent becomes the transitive order

    32,22,2,1,6,8,4,13,5,14,7.

A subset DP minimizes the number of backward arcs in a linear order and
returns2, proving arc-edit minimality for this grandparent. The corrected
base is strong with profile (1425,181,181,1425). Its full old-root clone family is

    H_old(m)=5407*2^m-4489,
    H_repaired(m)=1243*2^m-1061,                         m>=1.

The unique in-port permits at most two clone runs, so these counts lie in
the span of1 and2^m; exact counts at m=1,2 determine both coefficients.
BC1 supplies the matching certificate, and the word formula in RC3
independently recovers the repaired coefficients.

For m>=2 the repaired first-clone profile is

    (1243*2^m-1061, 1243*2^(m-1)-1062,
     1243*2^(m-1)-1061, 1243*2^m-1060).

The net count change is -4164*2^m+3428. This is a readiness repair with
an explicitly nonzero count debt. It is not a neutral return at6325/18979.

The [BC2 common-frame compiler](../boundary_carry_returns_20261103/PROOFS.md)
implements the two exterior reversals uniformly in m: four supplied-path
native complements, six temporary births/deletions, peak three helpers,
no old deletion. Root31, spare6, extra port2 and every clone address survive.
At source order n+m-1 with n=13, the expanded arc-reversal clock is
4*binom(n+m+1,2); constant native event count is not constant arc work.
Literal guards at m=1,2,4 and all five Q4 phases are recorded.

## RC3. A complete binary carrier for this ready-growth locus

By RC1, every qualifying rooted base consists of r, its unique in-port w,
and a transitive grandparent Q. Put Q in its unique order and encode
beta_j=1 iff w->Q_j. Strongness is exactly beta_last=0. The w cut on the
chain r1,...,rm,Q is 1^m beta. Conversely every such binary word gives an
actual strong tournament. Distinct words give distinct rooted base classes,
because the root determines w and Q's order is intrinsic.

At base order n>=3 there are2^(n-3) classes, with ordinary GF

    x^3/(1-2x).

Encode beta by the positive integer A whose little-endian bits are 1-beta_j;
the last bit is1. Write B(A)=H(beta)-A+(A mod2), where H(beta) is the one-apex
path count of the cut beta alone. Then, for every A>=1,m>=1,

    H_A(m)=A*2^m+B(A).                                      (RC3a)

**Proof of the count.** Use the direct weighted inversion formula of
[feedback-word Theorem1](../feedback_word_hamilton_spectrum_20261020/proofs.md):
H(b)=1+2 sum_(i<j) b_i(1-b_j)w(j-i), with w(1)=1 and
w(d)=2^(d-2) for d>=2. The inversions internal to beta contribute H(beta).
A zero at beta_1 contributes2^m; a zero at beta_j, j>=2, contributes
2^(j-1)(2^m-1). Summing proves RC3a. This recovers the earlier one-apex
mechanism; the target here is its complete readiness/carry classification.

Deleting a leading zero from a one-apex cut leaves H unchanged. A leading
one is an extra root clone. These facts give the exact recursion

    B(0)=1,
    B(2k)=B(k),
    B(2k+1)=B(k)-k-(k mod2).                               (RC3b)

Here A=0 is only a recurrence boundary, not a strong family. In particular
B(A) is odd, B(A)<=1, and B(A)=1 iff A is a power of two. Also

    B(A)>=1-A+(A mod2).

For even A>=2 the stronger lower bound B(A)>=3-A holds, because beta
begins in1, ends in0, and contains a directed triangle.

Structural halving is exact: H_(2A)(m)=H_A(m+1). The removed leading beta
letter is a real transitive twin of the root. This is a change of clone
presentation, with addresses retained; it is not division of H by two.

## RC4. Uniform ternary offsets have a sharp obstruction

Put Delta(A)=B(3A)-3B(A). A uniform identity between these families at the
same clone parameter, H_C(m)=3H_A(m)+c for all m>=2, forces C=3A and
c=Delta(A), by independence of2^m and1.

Let a be the odd part of A. Then

    a=1 mod4:              Delta(A)<=-4;
    a=3:                   Delta(A)=0;
    a=3 mod4, a>=7:        Delta(A)>=4.

Moreover Delta(A)=4 iff a=7. In particular **Delta(A) is never2**.

**Proof.** Define E(k)=B(3k+1)-3B(k), F(k)=B(3k+2)-3B(k), and
epsilon=k mod2. Direct substitution in RC3b gives

    Delta(2k)=Delta(k),       Delta(2k+1)=E(k)-2+4epsilon,
    E(2k)=Delta(k)-3k-epsilon,
    E(2k+1)=F(k)+3k+3epsilon,
    F(2k)=E(k),              F(2k+1)=F(k)-2+2epsilon.       (RC4a)

The lower bound on B(k) and upper bound B(3k)<=1 imply
E(2k)<=-2-4epsilon, so Delta(4k+1)<=-4.
E(1)=-2 gives Delta(3)=0. For k>=1, examine

    E(2k+1)=B(3k+2)+3k+3epsilon-3B(k).

For odd k>=3, B(k)<=-1 and the general lower bound imply E>=6;
for k=1 its exact value is2. For even k which is not a power of two,
the stronger even lower bound and B(k)<=-1 imply E>=4.
For k=2^j, the remaining exact values are4 at j=1,2, and
9*2^(j-2)-4>=14 for j>=3. Therefore E(2k+1)>=2 for k>=1, with
equality exactly at k=1. Substituting into Delta(4k+3)=E(2k+1)+2
proves every assertion. QED.

This excludes a **uniform family identity**, not an individual H value,
and not a return through larger-feedback sources. Even allowing a fixed
integer shift of the clone parameter does not repair +2: powers of two
transfer between A and2^m while B is invariant under that transfer.

The carry has a closed six-coordinate integer representation
(Delta,E,F,k,epsilon,1). On appending a binary digit it transforms by

    0: (Delta, Delta-3k-epsilon, E, 2k, 0, 1),
    1: (E-2+4epsilon, F+3k+3epsilon, F-2+2epsilon, 2k+1, 1, 1).

These are exact linear integer operators on the displayed coordinates.
Their binary-reading clock is distinct from physical attachment clocks;
no nonnegative orientation or native-move claim is made for the matrices.
On the repaired source A=1243, B=-1061, B(3A)=-1999, so Delta=1184.
The would-be +2 target misses by1182 in the constant mode at every m.

## RC5. The surviving +4 branch has a two-birth implementation

For A=7, beta=000. For3A=21, beta=01010. Therefore

    H_7(m)=7*2^m-5,
    H_21(m)=21*2^m-11=3H_7(m)+4.                         (RC5a)

Embed the old three Q vertices as positions0,2,4 of the new five-vertex
chain; birth positions1 and3 in the gaps, with both beaten by w and by
all root clones. Orient all chain relations by that order. Every old arc
is retained. Thus two literal births, zero native complements and no
old deletion implement RC5a uniformly in m. Old root and any old spare
survive. For m>=2 both endpoints have a1-ready first clone; the +4 identity
comes from this specified two-cut operation, not the1-ready reset.

The examples include23->73,51->157,107->325. The code retains complete
source/target rows, joint-port polynomials and Q4/F49 phases. This branch
does not iterate inside this carrier: oddpart(21)=1 mod4 puts its next
coefficient step on the negative side of RC4. A closed ternary source
family must enlarge the representation or operation/source class.
