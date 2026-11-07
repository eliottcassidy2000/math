# Arbitrary-apex chain words: a positive block compiler

Status: PROVED elementary construction and identities; FINITE-EXACT controls.
Chan–Pak completeness remains OPEN. No literature-priority claim.

## Inheritance and research board

Read-only source: eliottcassidy2000/math at
0f8750677dd20b594a0e032a3b56aa23f58a950c. The closest proved mechanism is
Session191's two-apex run transfer, with its 36 states and rank-17 observable
quotient. Its homogeneous-letter polynomial t(t-1)(t-2)(t-3) is inherited,
not a new empirical discovery. THM-4878 already gives more general finite
spine transfers; the novelty claimed here is the elementary commuting-
idempotent explanation, explicit all-k bound, and positive block expansion.

The hostile is Session214's all-order 7477 omission for one restricted exit;
allowing just one more exit supplies physical witnesses but not a complete
language. The corrected near miss is confusing a single block's scalar
count with its continuation state. Least-used sidecars: Stirling surjections
and THM-3134's centered permutation polynomial. Our live board is: arbitrary
apex masks; run endpoint states; homogeneous block powers; path-cover jets;
physical cyclic substitution. Anchor: positive richer constructions. Niche:
structural explanation of integer clocks. Wildcard: descent-polynomial EGF.

## 1. Exact all-k construction

Fix any tournament A on k>=0 labeled apexes, not necessarily transitive.
A word w over {0,...,2^k-1} gives a tournament T_A(w): its chain vertices
are increasing, and bit a of letter c means apex a beats that chain vertex.
All other incidences are reversed. These are actual tournaments and retain
arbitrarily many chain exits at each apex. Their union over k includes all
tournaments, so no completeness follows just from allowing unbounded k.

For each permutation p of the k apexes, a Hamilton word uniquely has form
L_0,p_0,L_1,...,p_(k-1),L_k. Each L_r is an increasing chain subsequence.
A state s has s_r in {0,1,2} for r<k: empty, nonempty with last vertex
beating p_r, nonempty with last vertex beaten by p_r. The last slot has
s_k in {0,1}: empty/nonempty. Thus the full space has k! 2*3^k states
(an inessential redundant state for k=0).

Define E_(c,r) on state basis vectors by assigning the next chain vertex
to slot r. If r>0, s_r=0, and bit p_(r-1) of c is zero, the result is zero:
the required incoming apex edge fails. Otherwise update only s_r, to 2 if
r<k and bit p_r is one, to 1 otherwise. The first slot has no incoming
restriction. Set M_c=sum_(r=0)^k E_(c,r), independently on each apex order.
Initialize all slots empty once for every p. At the end accept precisely
states with no s_r=2 for r<k, and, for each empty internal slot 0<r<k,
the apex arc p_(r-1)->p_r. Empty initial/final slots impose no edge.

Then H(T_A(w))=f M_(w_last)...M_(w_first) v. Proof: assign every chain
vertex to its unique run. First-edge checks occur upon opening the run;
last-edge checks occur only at the end, allowing a later vertex to repair
a bad last boundary. Empty internal runs use the actual apex arc. All
within-run edges hold because the chain is transitive. This is a bijection.
For k=0 the empty tournament is assigned one empty Hamilton word.

## 2. New structural all-k block identity

For fixed c and p, each E_r is idempotent, including its zero action:
a forbidden empty slot stays forbidden; an allowed update makes the slot
nonempty and further updates write exactly the same state. Distinct E_r
commute because they affect separate coordinates and their tests depend
only on those coordinates. These assertions hold on the full state space,
not merely on reachable states. Consequently, for every integer m>=0,

 M_c^m = sum_(U subset {0,...,k}) onto(m,|U|) product_(r in U) E_(c,r),
 onto(m,j) = sum_(i=0)^j (-1)^(j-i) binom(j,i) i^m = j! S(m,j).

Use 0^0=1, so the empty subset contributes identity exactly at m=0.
Proof: expand the m-fold sum into words of slot choices. Idempotence and
commutation reduce each word to its used-slot subset U. Exactly onto(m,|U|)
words have that subset. Although inclusion-exclusion evaluates the number,
the coefficient itself is nonnegative. Every term is a physical assignment,
not virtual subtraction of tournaments.

Commuting idempotents are simultaneously diagonalizable over Q: split into
image and kernel successively, which are preserved by every other E_r.
Thus every eigenvalue of M_c is a sum of k+1 bits, and

 product_(j=0)^(k+1) (M_c-jI) = 0.

Every fixed-context family T_A(u c^m v) therefore has, for m>=1, an exact
expression sum_(j=1)^(k+1) alpha_j j^m with rational alpha_j, and a possible
separate m=0 value. No polynomial-in-m factors occur. At k=2 this explains
the inherited Session191 identity rather than claiming it anew. For any
number of variable homogeneous blocks, multiply the explicit block
operators in order; this gives multivariate finite sums of products of
integer clocks. Different letters need not commute. The full boundary
vector, not H alone, is the compositional state.

## 3. A descent-polynomial window into cyclic substitution

For k>=1 take w=(2^k-1)^r 0^s. This is exactly C3[A,TT_r,TT_s] when
r,s>0. Thus general substitution is INHERITED from THM-3121 and its
path-cover/centered-permutation equivalence from THM-3134. The following
EGF is a direct specialization/presentation, not a new substitution law.
It is also valid at r=0 or s=0 as a chain-word construction.

Let d(p) count backward adjacencies in permutation p of A. Put
Z=(exp(x)-1)(exp(y)-1),
F=1+exp(x+y)-exp(x), G=1+exp(x+y)-exp(y). Then

 H(T_A((2^k-1)^r 0^s))
 = r!s! [x^r y^s] F G sum_p Z^d(p) (1+Z)^(k-1-d(p)).

Proof: each internal run is empty or contains vertices of both blocks,
with EGF Z. An empty internal run is allowed exactly at an apex forward
adjacency. The initial run is either empty or has at least one second-block
vertex with any number of first-block vertices, giving F. The final run is
either empty or has at least one first-block vertex with any number of
second-block vertices, giving G. Exponential coefficients account for
assigning the labeled ordered vertices of each block to runs; within-run
orders are forced. Summing apex orders proves the formula. In THM-3134
notation the sum is exactly E_A(Z+1/2), the inherited centered polynomial.

Connection contract: map apex A to cyclic two-cap tournament; preserved
information is its complete descent polynomial (equivalently full ordered
path-cover profile at fixed order), not its adjacency matrix. The destroyed
information is individual apex labeling/edge detail beyond this profile.
Needed sidecar is block lengths r,s. Cheapest test is the two independent
slot and exponential-polynomial calculations in verify.py. The companion
fixed-cap observer theorem is the more general
profile-recovery framework; this is its transparent homogeneous-block view.

## 4. Explicit fourth-clock families and mixed masks

Take cyclic apex rows A=[2,4,1] (0->1->2->0), r=3. Its descent counts are
3 at d=0, 0 at d=1, 3 at d=2, so the polynomial in Z is 3+6Z+6Z^2.
Expanding the displayed EGF gives, for every m>=1,

 H(T_A(7^3 0^m)) = 87 - 189*2^m + 72*3^m + 36*4^m.

The m=0 value is 3, and must not be obtained by this positive-clock formula
(which would give 6). This boundary is an explicit hostile control.
For m=3 the nine-vertex outgoing bit rows are
[58,60,57,496,480,448,391,263,7], and H=2823.
The 4^m coefficient is nonzero, so this sequence cannot have a fixed-context
homogeneous-clock representation with <=2 apexes: such sequences use only
1,2,3, and distinct exponential sequences are linearly independent. This
is a sequence separation, NOT a claim that individual counts require three
apexes. The cyclic symmetry forces divisibility by 3 here, so this family
alone plainly cannot solve Chan–Pak.

More generally any k-apex A with initial block length exactly k and final
block length m has highest-clock coefficient k!^2 at (k+1)^m. To use all
k+1 slots during the 0-block, the initial k all-ones vertices must occupy
each slot 1,...,k once (k! assignments); there are k! apex orders. Conversely
these assignments allow the high eigenmode, with final nonempty conditions
contributing inclusion-exclusion lower clocks. No other initial assignment
can access k+1 active slots, because a zero-letter cannot open slots 1,...,k.
The argument establishes coefficient k!^2 for k>=1; k=0 is the constant
chain boundary. These are positive families with arbitrarily high clocks.

The compiler also allows nonuniform masks, beyond the homogeneous cap
specialization: for A=[2,4,1] and word [7,1,2,4,0], the receipt gives an
explicit eight-vertex tournament, checked by independent subset DP. No
coverage or isomorphism-minimality claim is attached to this example.

## 5. Verification, scope, and next experiment

Run python3 verify.py > RESULT.replay.json and compare with RESULT.json;
repeat under python3 -O. All tests raise explicit exceptions, not assert.
The receipt states exact finite universes: all labeled apex tournaments
k<=3 and all words to lengths 5,6,5,3 respectively for k=0,1,2,3;
all states for algebra checks k<=4; all apex tournaments k=1..4 and two
block lengths 0..3 for the EGF; fixed-context block checks with k=3; and highest-clock extraction for every
labeled apex tournament at k=1..4 (75 cases), by applying
(shift-1)...(shift-k) to positive-clock values m=1,...,k+1.
Independent vertex-subset DP does not use the run decomposition.

The all-k proofs are not inferred from these finite tests. The full transfer
method is inherited, the all-k idempotent/Stirling proof is the present
structural contribution, the cyclic EGF is a direct specialization of
inherited substitution, and displayed coefficients/witnesses are exact
corollaries with finite independent controls. This does not settle a
fixed-k spectrum or supply an all-odd inverse compiler. The next useful
experiment is arbitrary-mask three-apex word search with continuation-state
synthesis, retaining the multi-coordinate boundary rather than fitting H.

## Source routes

All paths below are pinned, with git blob and SHA256 in SOURCES.json:

- 05-knowledge/results/session191_tournament_growth_20261104/two_apex/TWO_APEX.md
- 01-canon/theorems/THM-4878-graded-spine-grammars-cut-state-automaton-hankel-and-seam-ranks-cutwidth-censuses.md
- 01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md
- 01-canon/theorems/THM-3134-tournament-endpoint-jet-and-c3-newton-profile-transform.md

Use the actual pinned paths in SOURCES.json if a theorem's maintained slug
changes; IDs alone are not source identity.
