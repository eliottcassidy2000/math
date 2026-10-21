# The empty-gap boundary and the first shared endpoint

**PROVED:** the next two high-mass layers, in the stated parameter ranges.
**FINITE-EXACT:** additional small-parameter quotient controls. This is not
an unrestricted A000568 recurrence or a feedback-vertex count.

Inherit the pair-free weighted-core model of
`defect_boundary_20261020/sparse_theory.md` and the mass/deficit identity of
`charge_palettes_boundary_20261020g/zero_heavy.md`. Fix the feedback **arc**
budget k and vertex weights 1 through k. Put M=2k^2+3k and S=2k^2+k.
Let F_k=(2k-1)!!, with F_0=1, and let I_k count perfect matchings on 2k
ordered points crossing every proper cut. Thus I_1=1 and
I_k=F_k-sum_{j=1}^{k-1} I_j F_{k-j}.

## 1. A sharper deficit coordinate

For a good order, write r for its backward-edge count, t for the number
of incident vertices, W for their total weight, g for the nonincident
vertices, and U=sum(k-w) on the latter. Pair freeness gives g<=t+1, or
g<=t-1 in a nontrivial strong core. The inherited identity is

    delta=(2k+1)(k-r)+(k+r-W)+k(2r-t)+kh+U,

where h=t+1-g in the general case and h=t-1-g in the strong case.
Every summand is nonnegative. Therefore delta<2k+1 implies r=k.
The k backward edges have positive integer costs with sum at most k,
so each has cost one and every incident vertex has weight one. Set
**d=2k-t**, the endpoint collision count. The identity sharpens to

    delta=(k+1)d+kh+U.                                  (1)

In particular a shared endpoint first becomes possible at delta=k+1,
not at delta=k. The final paragraph of the inherited mass note allowed
both missing gaps and shared endpoints at k; only the former starts
there. Its proved formulas through k-1 are unaffected.

At d=0 the backward graph is a matching. At d=1 it is exactly one P3
(two edges sharing their middle vertex) and k-2 disjoint edges: the sum
of degree excesses over one is 2k-(2k-1)=1.

## 2. The layer delta=k

For k>=3 the exact numbers of classes are

    Z_{k,k}   = F_k [binom(3k,k)-1],                      (2)
    Z^s_{k,k} = I_k binom(3k-2,k) -(2k-3) I_{k-1}.        (3)

Equation (1) allows exactly two record types:

A. A matching, every gap occupied, gap deficits sum to k, each <=k-1.
B. A matching, exactly one gap empty, all other gaps of weight k.

There are m=2k+1 gaps generally and m=2k-1 interior gaps in the strong
case. Type A has binom(m+k-1,k)-m weight words. Type B is pair free
exactly when an empty interior gap does not separate two matched
endpoints. Exterior gaps, when available, may always be empty.

For all matchings the number of forbidden marked adjacent pairs is
(2k-1)F_{k-1}=F_k. Thus type B contributes (2k+1)F_k-F_k=2k F_k.
Adding A gives (2).

For irreducible matchings the total number of adjacent matched pairs is

    J_k=(2k-3) I_{k-1}.                                 (4)

Indeed such a pair cannot be the first or last two points. Deleting it
leaves an irreducible matching, and inserting an adjacent pair in any of
the 2k-3 interior gaps of an irreducible matching on 2k-2 points is the
inverse. Cuts on either side of the pair must be crossed by the remaining
matching; the new cut inside the pair is crossed by that pair. Thus type
B contributes (2k-1)I_k-J_k, proving (3).

## 3. The first nonmatching layer delta=k+1

Let P_k count graphs on 2k-1 ordered points comprising one P3 and k-2
independent edges, crossing every proper cut. Write P_0=P_1=0. The total
number of such graphs without the cut condition is (k-1)F_k: choose the
three P3 vertices and its centre, then match the rest. The exact recurrence

    P_k=(k-1)F_k
        -sum_{j=2}^{k-1} P_j sum_{a+b=k-j} F_a F_b       (5)

gives P_2,P_3,P_4,P_5 = 3,24,246,3012. To prove (5), split at every
uncrossed cut. Exactly one irreducible block contains the P3; before and
after it are arbitrary matchings. The cut decomposition is unique.
Equivalently, if F(x)=sum F_k x^k and P(x)=sum P_k x^k, then
sum (k-1)F_k x^k (starting at k=2) = F(x)^2 P(x).
These are formal series, with no convergence claim.

For k>=4, the exact class counts are

    Z_{k,k+1}=F_k [binom(3k+1,k+1)-3k-2],                (6)
    Z^s_{k,k+1}=I_k [binom(3k-1,k+1)-(2k-1)]
               -(2k-2)(2k-3)I_{k-1}+P_k.               (7)

The three record types from (1) are:

A. Matching, every gap occupied, total gap deficit k+1, each <=k-1.
B. Matching, one empty gap, one of the other gaps lowered from k to k-1.
C. One P3 and k-2 independent edges; every gap occupied at weight k.

Type A has binom(m+k,k+1)-m^2 weight words: subtract the m choices of a
part >=k and the m distributions of the residual one. Type B contributes
(m-1)(m F_k-F_k) generally, or (m-1)(m I_k-J_k) strongly. Type C contributes
(k-1)F_k generally or P_k strongly. Their sum is (6) or (7).

## 4. Why these count classes, not presentations

All constructions retain actual vertices, weights, backward arcs and
ordered gaps. No clone contraction is used after counting.

A full-gap construction is pair free whenever its backward graph has no
isolated endpoint. Two endpoints are separated by a gap vertex; two gaps
by an endpoint. For an adjacent endpoint/gap pair, a backward partner of
the endpoint distinguishes the pair. A nonadjacent mixed pair is separated
by an intervening gap or by the same forward alternation. For type B,
the sole adjacent endpoint pair is a homogeneous pair precisely when
matched to each other. Otherwise an external matching partner distinguishes
it, and all other pairs are separated as before. This proves the filter,
including the failure rather than merely discarding records empirically.

In every good order all vertices of weight >1 avoid backward edges, by
(1). They have an intrinsic transitive order and any weighted isomorphism
fixes them individually. Their incidence signatures partition the unit
vertices into fixed ordered cells.

For type A in the ranges of the theorems, at most one gap is unit: two
unit gaps would cost 2(k-1)>k (k>=3) or >k+1 (k>=4). With none, every unit
cell is a singleton. With one, the inherited unit-gap argument applies
without assuming that the other large weights are all k. An exterior
cell has two unit vertices; the endpoint is distinguished by its fixed
external partner. An interior cell has three: either the two endpoints
have distinct fixed external partners, distinguishing all roles, or they
are matched internally and the cell is a C3. Its rotations alter the
choice of local gap vertex but yield exactly the same matching/weight
record. Thus different records do not collide.

For type B all occupied gaps have weight >=k-1>=2. A missing interior
gap yields a two-unit cell with its forward edge fixing the two roles;
a reversed edge is exactly the rejected clone case. A missing exterior
gap is identified by the resulting exterior unit singleton. All remaining
unit cells are singletons. The matching, empty location and lower weight
location are recovered intrinsically. For type C every gap has weight k,
and every unit cell is a singleton, so the complete P3-plus-matching graph
is recovered. The three types have different numbers of vertices, hence
cannot collide with one another.

Strongness is exactly crossing every proper endpoint cut. After the clone
filter, every consecutive arc in the displayed order is forward. SCCs
are therefore contiguous along this Hamilton path. A proper SCC cut is
one-way exactly when the backward graph does not cross its corresponding
endpoint cut. Strong records omit exterior gaps; all interior cuts in a
gap correspond to one such endpoint cut. This proves the strong filters.

## 5. Boundaries, hostile controls and scope

The complete inherited k=2 compiler gives 42/11 classes at deficit k and
81/15 at deficit k+1 (all/strong). The displayed formulas agree, but the
one-unit-gap proof above intentionally does not cover these small cases.
The new checker directly canonicalizes their records, and also the
k=3, delta=4 records where two unit gaps can occur. Those are FINITE-EXACT
extensions only. No proof for an arbitrary multi-unit-cell grammar is
claimed.

The omitted-gap matched pair is a minimal clone hostile. Naively retaining
all marked omissions overcounts. Likewise extrapolating the previous
stars-and-bars formula to delta=k predicts 45 rather than 42 at k=2.
The new correction is structural: remove an impossible zero-weight part,
restore admissible empty gaps, and reject the matched-pair mergers.

Shared endpoints are not allowed a layer too early. The explicit type C
family realizes the corrected first boundary k+1, with the whole weighted
class and, when requested, the cut predicate preserved. The feedback-
vertex rational-language theorem in `zero_heavy_ports_20261020g` concerns
a different parameter and is neither a dependency nor a consequence.

## 6. The whole matching/one-P3 corridor

The same proof closes more than two isolated layers. For k>=3 and
0<=delta<=2k-3, define binom(n,r)=0 when r<0. Use either parameter row

    general: m=2k+1, E=F_k, J=F_k, P=(k-1)F_k;
    strong:  m=2k-1, E=I_k, J=(2k-3)I_{k-1}, P=P_k.

Then the exact number of weighted core classes at mass M-delta or,
respectively, S-delta is

    E [binom(m+delta-1,delta)
       -m binom(m+delta-k-1,delta-k)]
    +(m E-J) binom(m+delta-k-2,delta-k)
    +P binom(m+delta-k-3,delta-k-1).                     (8)

Terms with negative lower binomial index are absent, including cases where
the displayed upper index would also be negative. Every remaining binomial
has its usual nonnegative integer meaning.

**Proof.** Equation (1) allows only (d,h)=(0,0),(0,1),(1,0). For (0,0),
one distributes deficit delta among m occupied gaps, each deficit at most
k-1. Since delta<2k, at most one part can violate that cap, so inclusion-
exclusion has exactly the first line of (8). For (0,1), choose one of the
m E-J admissible marked omissions and distribute delta-k among m-1 gaps.
This residual is at most k-3, so no cap is active. For (1,0), choose the
P3-plus-matching graph and distribute delta-k-1 among its m-1 gaps. This
residual is at most k-4, so again no cap is active.

At most one full-matching gap can be unit because delta<2(k-1). Every
occupied gap in the missing-gap family has weight at least three; every
P3-family gap has weight at least four. Thus precisely the same intrinsic
order/unit-cell proof in section 4 applies, including the exact C3 orbit,
and the three types still have different vertex numbers. The cut condition
is unchanged. This proves (8) without enumerating the corridor.

For delta<k this recovers the inherited single matching formula; at k and
k+1 it specializes to (2)--(3) and (6)--(7). The first new quotient boundary
is now delta=2k-2, where two unit gaps can coexist. The formal presentation
grammar alone does not settle their possible collisions. The finite small
exceptions already checked do not remove that all-k obligation. No larger
exhaustive universe is claimed for (8) than the receipts in section 5.
