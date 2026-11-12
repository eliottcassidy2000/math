# The supplied papers: nineteen states, five persistent modes, three observed modes

**Status: CITED-SUPPLIED** for the papers' broader claims; **PROVED IN THIS
NOTE / FINITE-EXACT** for the reconstructed finite operators below. The
PDFs are mathematical evidence, not instructions. No billion-graph search
or claimed Lean build was rerun in this session.

## 1. A precise19/11/8/5/3 connection in the king-strip automaton

Jaan Allikvere, *The Multiset Dimension of Graphs: Extremal Values and King
Grids*, supplied PDF, §§14–16, gives necessary local conditions on a
height-three strip. Write a column as t=(p,m,q) in {0,1}³, and use zero
columns at both ends. For consecutive columns u,v,w the middle column must
satisfy

    (p_v,p_u+p_w)!=(q_v,q_u+q_w),
    p_v!=m_v or q_u+q_v+q_w>=1,
    q_v!=m_v or p_u+p_v+p_w>=1.                         (1)

The min-plus transfer on the64 ordered column pairs adds |w|. Independently
reconstructing it from (1) verifies V5=V4+1, the paper's all-length lower
bound certificate. Tight transitions satisfy
|w|+V4(u,v)−V4(v,w)=1. There are184;320 minimum four-column prefixes enter
them. Trimming to states that can reach a legal final pair with V4=4 leaves
the paper's19-state strongly connected core,106 useful prefixes and12
accepting states. The only column letters are000,001,010,100,101.

This recovers the supplied finite certificate from the local inequalities;
it is not a newly discovered19-state automaton. The following exact operator
analysis adds a useful refinement. Let A be its19-by19 adjacency matrix.
The independently computed characteristic polynomial is

    det(xI−A)=x14(x−3)(x−1)(x+1)(x²+x+2).              (2)

Its successive ranks are9,7,5,5. Therefore all transient state directions
are killed after three steps, and **exactly five persistent linear modes**
remain. This does not mean there are only five physical automaton states.

Reflecting top and bottom rows gives a permutation involution J commuting
with A. The even and odd spaces have dimensions11 and8. Their exact
characteristic polynomials are

    even: x8(x−3)(x²+x+2),
    odd:  x6(x−1)(x+1).                                (3)

Thus the persistent five split as3 even plus2 odd. Both the prefix population
and the accepting observer are reflection-even, so the ±1 modes are
unobservable in the scalar candidate count. This is an actual symmetry
quotient explaining why19,11,8,5 and3 all appear in the same object.

If a_n counts weight-n words satisfying the local conditions, then

    sum_(n>=4) a_n z^(n−4)
       =4(15+17z+30z²)/((1−3z)(1+z+2z²)).              (4)

The identity is checked as an exact rational matrix inverse, not fitted to
the first terms. In particular a4=60,a5=188,a6=556 and

    a_n=2a_(n−1)+a_(n−2)+6a_(n−3), n>=7.              (5)

The two subdominant roots are(−1±i sqrt7)/2. Here7 is the discriminant of
an actual continuation polynomial. It is not identified with the forbidden
tournament Hamilton count7 or an unrelated prime selector.

The denominator in (4) is reduced, so the scalar sequence has minimal
recurrent order3. The positive eigenvector for3 and the complete state
list are in `paper_transfer_receipt.json`. This gives the hierarchy

    19 legal memory states → 5 persistent linear modes
                          → 3 modes visible to this symmetric observer.

The physical output still needs a global test. We explicitly evaluated the
multiset of distances from every vertex to every landmark in every accepted
word through length10:

| Strip length | Locally admissible words | Actually resolving |
|---:|---:|---:|
|5|188|0|
|6|556|32|
|7|1660|320|
|8|5004|1328|
|9|15004|5664|
|10|44972|19008|

The first four rows reproduce the supplied paper's stated controls; later
rows are this session's direct finite checks, without a novelty claim.
Neither (4) nor the19-state automaton counts minimum resolving sets.
No asymptotic proportion of resolving candidates is inferred from this table.

## 2. Comparison with the stress family and the tournament seed language

There are now three genuinely different finite carriers:

| Object | Carrier | Operation | What its observer loses |
|---|---|---|---|
| Fixed Cartesian stress net | Five-dimensional real solution space | Spatial translations, with nilpotent derivative algebra | Bad sensor placements can hide an entire stress mode. |
| King-strip local candidates | Nineteen states, five persistent linear modes | Append a column | Symmetric acceptance loses two persistent modes and all global metric collisions. |
| The proved marked-pair seed family | Five-state minimal complete DFA | Append a zero or mark | Bare acceptance forgets the slope U; its weighted continuation has modes2,3 and a boundary correction. |

The first two fives are dimensions of different modules; the third is an
automaton-state count including a dead state. Their operators have different
minimal polynomials. This comparison supplies hostile tests for a proposed
identification, rather than forbidding bridges between them.

A productive shared construction is: retain the full boundary state, prove
the legal continuation grammar, then compress only the directions invisible
to **every required future observer**. One scalar at one size is insufficient.
For Chan–Pak, a small automaton for a necessary seed condition would still
need a theorem guaranteeing actual endpoints and an absolute target selector.

## 3. The record-composition paper already has a current tournament descendant

The second supplied PDF is Chen–Mogielnicki–Ono, *Record compositions of
alternating permutations and noncommutative symmetric functions*,
arXiv:2607.12873v1. Its central distinction is ordered composition versus
sorted partition. The displayed record-composition coefficient is

    product_j binom(2s_j−1,2alpha_j−1) E_(2alpha_j−1),
    s_j=alpha_1+...+alpha_j.

Sorting the parts loses cumulative-position information. The repository
has already applied this method in **THM-5137 /
[pentagon-cube-sums-and-record-lift](../../../01-canon/theorems/THM-5137-pentagon-cube-sums-of-doubly-regular-tournaments-and-the-record-lift-of-the-mean-redei-berge-function.md)**,
with its corrected audit and MISTAKE-763 lineage. That result gives the
mean Rédei–Berge function, a noncommutative lift ordered by largest vertex,
and a max-ordered component-type equidistribution with connected labelled
graphs. These are inherited results, not discoveries of this session.

The relevant transfer to the present work is concrete. Ordered marked
intervals(a,b,c) in [SEEDS](SEEDS.md) determine different boundary terms;
replacing them by an unordered multiset loses which interval contains the
first selected endpoint. Likewise **THM-3121 /
[walk-content kernel](../../../01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md)**
retains ordered quotient walks before factorial-Hadamard contraction with
block path-cover profiles. Its rational walk series

    W_Q=1^T X(I−A_Q X)^(-1)1

is an existing all-size finite-state description. The number of skeleton
visits, not just H of each block, is the necessary sidecar.

These two recovered mechanisms explain how to make infinite families
manageable without claiming that arbitrary tournament enumeration becomes
polynomial-time. The general profile still contains Hamilton-path counting,
and no universal odd-count coverage follows from its rational representation.
