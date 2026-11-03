# Adaptive insertion: exact root changes, an inverse detector, and a sharp state loss

**Status: PROVED for the identities and structural consequences below;
FINITE-EXACT for the stated minimality and restricted-language exclusions.
No all-target prime or seven-prime construction is proved.**

The operation is always on an actual named tournament. Given a root v in T,
adjoin a new vertex x with v→x and x→u for every old u≠v; retain every old
arc. The adaptive language starts from the named triangle
0→1→2→0 and allows any old root at every step. Names of new vertices are
successive integers. H counts directed Hamiltonian vertex words, with
H(empty)=1. All the restricted languages below use this same operation.

## 1. Inheritance and the live board

The closest proved mechanism is Session 150's
[rooted block formula](../session150_hamiltonian_spectrum_20261103/gadgets/RESULT.md)
and [fixed-root profile dynamics](../session150_hamiltonian_spectrum_20261103/prime_induction/RESULT.md).
Its [adaptive computation](../session150_hamiltonian_spectrum_20261103/adaptive/RESULT.md)
already reaches every odd count 3 through 199 except 7 and 21. That is a
complete bounded computation for its language, not an all-order induction.
The present packet does not extend that numerical coverage record.

The canonical hostile is the regular strong tournament on five vertices:
it has no indegree-one vertex, so cannot be an adaptive child. Its count 15
is nevertheless attained by the adaptive history (0,0,4). Thus failure to
generate a class is different from failure to realize its count.

The corrected near miss is the proposed insertion identity 2h−s. The true
identity is h+F−s: deletion can create a broken adjacency. Likewise,
strict growth follows from F>q≥s; the stronger assertion q>s is false.

The least-used relevant sidecar here is the mixed endpoint path-cover
response. Its general mechanism is inherited from
[THM-4099, squarefree gap transfer and mixed insertion](../../../01-canon/theorems/THM-4099-squarefree-gap-transfer-and-mixed-insertion-boundary.md),
and its pair positivity is inherited from
[Session 124, canonical weighted-cut law](../session124_tournament_cut_response_20261031/README.md).
Section 5 specializes that mechanism to the actual adaptive operation;
it does not claim the general path-cover hierarchy as new.

| Live concept | Representation and operation | Result of the hostile probe / next boundary |
|---|---|---|
| Anchor: prime and 7p induction | Named adaptive histories, any old root | Exact class inverse exists, but scalar target inverse remains open. |
| Closed arithmetic subsystem | Keep the root or mark the newborn; four coordinates | Exact matrices and a reset exist; this language has permanent count gaps. |
| Missing switching data | Same rooted profile, then arbitrary old-root switch | Sharp strong adaptive order-six pair has different attainable next counts. |
| Inverse sidecar | Rooted defect and deletion count | The profile detects legal last vertices; the parent profile still needs a two-deletion coordinate. |
| Wildcard: grammar ancestry | Newest vertex in the earlier five-letter unique-cycle grammar | The predecessor mechanism extends to unrestricted adaptive histories; the languages are not identified. |

The last connection recovers the narrower newest-maximum argument in
[Session 103, forbidden words and strong bridge cuts](../session103_words_necklaces_20261029/singularity/FORBIDDEN_WORDS_ARE_STRONG_BRIDGE_CUTS.md).
The general unique-inverse and rigidity theorem is in [ANCESTRY.md](ANCESTRY.md),
from the parent lane and independently checked here. No external priority
claim is made.

## 2. Exact root changes and a positive reset

For a rooted actual tournament (T,v), write

    h=H(T),  s=# Hamiltonian paths of T starting at v,
    q=H(T−v),
    F=# ordered spanning two-slot covers (P,Q) of T with P ending at v.

P is nonempty and Q may be empty. This F equals Session 150's definition:
remove v from the end of P; the remaining path is empty or ends at an
in-neighbor of v. Empty slots have one empty path, never an extra ordering.

After insertion at v, let R mean that v stays marked and N mean that the
new vertex x becomes marked. The exact maps on z=(h,s,q,F) are

    Rz = (h+F−s, q, q, 2F−q),                         (2.1)
    Nz = (h+F−s, h−s, h, h+F).                       (2.2)

Equation (2.1) is inherited. For (2.2), a path starting at x is exactly an
old Hamiltonian path not starting at v, prefixed by x, so s_x=h−s. Deleting
x gives q_x=h. In a two-slot cover with the first path ending at x, that
path is either the singleton x, contributing h, or is an old path ending
at v followed by x, contributing F. Hence F_x=h+F. This proves (2.2)
for every nonempty actual T, without a strongness assumption.

Applying N and then R gives the actual two-birth reset

    R(Nz) = (h+2F, h, h, h+2F).                      (2.3)

Call a profile (b,t,t,b) balanced. R preserves balanced profiles, sending
(b,t) to (2b−t,t). From a balanced source with count b, the same two-birth
reset gives

    R(N(b,t,t,b)) = b(3,1,1,3).                      (2.4)

These are constructions, not assertions that arbitrary numerical profiles
are realized. Every map is applied to its named source and specified cut.
Because R and N are linear, every further keep/new word from (2.4) has
count b times the count obtained by that same word from the rooted C3.

**Factor obstruction, with its guards.** If b>1 and b≠7, no such reset
followed only by keep/new moves can attain an odd prime p or 7p. Both
factors are integers at least 3; neither can equal 7. The C3 factor cannot
be 7: the first insertion has count 5, the next has count 9, and counts
strictly increase. For an adaptive source the condition b≠7 follows from
the same inherited all-root argument. The b=1 exception is real: the
singleton reset constructs C3. For an arbitrary source, retain the b≠7
guard rather than silently asserting source membership or a global gap.

In the C3 keep/new language, every R-state is balanced. Consequently a
single N between the start (or an R) and the next R creates this factor
obstruction. This is a necessary word-pattern exclusion for prime/7p
outputs in that restricted language, not a sufficiency criterion and not
an obstruction to switching to an older vertex.

### A positive three-coordinate form on the zero-defect language

The coordinator extracted this corollary after reading the concurrent golden
single-bridge theorem. On the invariant locus delta=0 from Section3, write

    (a,q,s)=(h-s,q,s),    h=a+s,    F=a+q.

These are nonnegative physical counts on actual sources. Equations(2.1)–(2.2)
become

    R(a,q,s)=(2a,q,q),
    N(a,q,s)=(a+q,a+s,a).                              (2.5)

Starting at the rooted C3 vector(2,1,1), these two nonnegative integer
matrices describe exactly the degree-cap-one adaptive language. After a
birth the only indegree-one roots are the protected root and the newborn,
so its two choices are precisely R and N. This is a three-dimensional
weighted representation, not three possible numerical states.

It has minimum rational linear dimension3. Use actual prefixes empty,R,N
and continuations empty,R,RR. The Hamilton-response matrix is

    [[3,5,9],[5,9,17],[5,9,15]],    determinant=-4.

Every linear representation factors this matrix through its state space,
so dimension at least3 is necessary. The three observations also recover
the coordinates: a=(H_RR-H_R)/2, q=2H_R-H_RR, s=H-a. Thus they give the
coarsest deterministic future-H observer within this operation language.
The explicit15-coordinate path-forest representation in ACTIVE_PORTS is
redundant; its coordinate count was never a minimum-rank assertion.

There is a literal connection to the incoming
[golden single-bridge recurrence](../golden_hamilton_transfer_20261103/RECURRENCE.md):
on the converse source T^op its coordinates(A,B,C) equal(h-s,s,F).
Joining one singleton by its exceptional backward bridge and then taking
converse is exactly the adaptive birth with the newborn marked. The golden
map(C,A,A+B+C) agrees with N under this identification. Keeping the old root
is a different command, and(2.5) uses delta=0. General old-root switching
does not inherit this three-coordinate closure; Section4 is its hostile.

## 3. The defect detects indegree and the legal inverse

Define

    δ_v = F_v−h−q_v+s_v,    S=N⁻_T(v),    d=|S|.

For distinct a,b in U=T−v, let W_ab(U) count ordered covers (A,B) by two
nonempty paths, A ending at a and B starting at b. Then

    δ_v = Σ_{a,b∈S, a≠b} W_ab(U)
        = 2 Σ_{unordered {a,b}⊂S} b_ab(U),             (3.1)

where b_ab=(W_ab+W_ba)/2 is Session 124's positive pair weight. In
particular its proved bound gives

    δ_v ≥ 2 Σ_{{a,b}⊂S} H(T−{v,a,b}) ≥ d(d−1),
    δ_v=0  iff  d≤1.                                 (3.2)

To prove (3.1), splitting a Hamiltonian path just after v injects the h
paths into the F covers. The covers left over have Q nonempty and starting
in S. Those with P=(v) number q−s. In every remaining cover, remove the
last vertex v from P. Both surviving paths are nonempty; the first ends
in S and the second starts in S. This is the stated bijection. Positivity
and integrality in (3.2) are inherited from the weighted-cut theorem,
including the empty residual tournament convention.

The maps have exact defect laws

    δ(Rz)=δ(z),     δ(Nz)=0.                          (3.3)

These follow by substitution and explain the reset: N gives an
indegree-one root, and the following R sets s=q as well. For a strong
tournament, δ=0 is exactly indegree one, since indegree zero is excluded.

**Profile criterion for the legal last vertex.** In any strong T of order
at least three, v has indegree one and T−v is strong if and only if

    δ_v=0  and  q_v>s_v.                             (3.4)

Indeed, let u be the sole in-neighbor of v. Then q−s counts the Hamiltonian
paths of U=T−v starting at u. Strongness of T implies that every old vertex
reaches u within U, so u belongs to U's terminal strong component. If U
has more than one strong component, no Hamiltonian path starts there.
If U is strong, a Hamiltonian path starting at u exists: choose a longest
path with fixed start u. Any vertex that cannot be inserted after the
start or at the end must dominate the whole path; if the path were not
spanning, all outside vertices would dominate it, contradicting strongness.
This proves both implications of (3.4), using (3.2).

[ANCESTRY.md](ANCESTRY.md) proves that at most one vertex satisfies this
criterion when |T|>3. Repeated deletion down to C3 recognizes the adaptive
language. The parent count is exactly q_v, so this supplies an intrinsic
decreasing count attached to the unique object inverse. It does not
construct a source for an arbitrary prescribed integer.

At a legal newest vertex x, (2.2) recovers from its child profile

    h_parent=q_x,   s_parent=q_x−s_x,   F_parent=F_x−q_x.

The missing fourth coordinate is
q_parent=H(T_child−{x,u}), where u is x's unique in-neighbor. Thus the
profile detects the legal inverse and its count, but its own inverse
profile still needs the named root and a two-deletion sidecar.

## 4. A sharp failure of arbitrary root switching from four coordinates

**FINITE-EXACT minimal order; exact named counterexample.** Let an adjacency
row be the bitset of outgoing neighbors, with vertices numbered 0,...,5.
Use the two strong adaptive tournaments

    A=(40,17,3,6,13,30), marked root 2,
    B=(34,24,3,5,13,30), marked root 3.

They are not ordinarily isomorphic, but both marked profiles are

    (n;h,s,q,F)=(6;37,5,13,71).

Insert vertex 6 at the marked root. Both children have H=103. A further
insertion at each named root 0,...,6 has the following count:

| Root for the second insertion | 0 | 1 | 2 | 3 | 4 | 5 | 6 (newborn) |
|---|---:|---:|---:|---:|---:|---:|---:|
| From A | 287 | 321 | 219 | 329 | 275 | 233 | 179 |
| From B | 289 | 307 | 339 | 219 | 271 | 227 | 179 |

Hence neither the set nor the multiset of the next attainable counts is a
function of (n;h,s,q,F). In particular, the named A source reaches 233 and
287=7·41 here, while B does not; B reaches 307, while A does not. These are
ordinary prime and 7p examples, not new spectrum attainments. The keep-root
value 219 and newborn-root value 179 agree, precisely as (2.1)–(2.2) require.

Actual construction histories make membership and port names reproducible:

| Source | Root history from C3 | Map from insertion labels to displayed labels | Marked insertion label |
|---|---|---|---:|
| A | (0,1,2) | (3,1,0,2,4,5) | 3 |
| B | (0,2,1) | (3,0,1,2,4,5) | 0 |

For example, complete histories (0,1,2,3,5), (0,1,2,3,2), and
(0,2,1,0,2) attain 233, 287, and 307, respectively. The checker verifies
the named transports, ordinary nonisomorphism, and the final counts by
all 8! vertex orders as well as subset DP.

Minimality is relative to this precise claim: among every ordinary strong
tournament of each fixed order 3,4,5 and every marked root, equal profiles
have equal multisets of counts after the first marked insertion followed
by a free root choice. The first failure is at order six, and both sources
already lie in the adaptive language. No minimality over arbitrary larger
statistics or other operation families is asserted.

The first failed implication is merging actual sources by their four
coordinates before allowing arbitrary root changes. The strongest survivor
is exact closure under keep/new moves. The repaired state is a named source
with mixed endpoint response data, as follows.

## 5. The exact mixed sidecar for an older root

All covers in this section are of the same old named T. For distinct u,v,
define:

* A_uv: ordered two-slot covers (P,Q), P starting at u and ending at v;
* B_uv: ordered two-slot covers (P,Q), P ending at v and Q starting at u;
* G_uv: ordered three-slot covers (P,Q,L), P ending at u and Q ending at v.

Every constrained slot is nonempty; Q in A and L in G may be empty.
Thus G_uv=G_vu by exchanging the first two slots.

Insert x at old root v and mark a different old root w. Its exact profile is

    h' = h+F_v−s_v,
    s'_w = A_wv,
    q'_w = H(T−w)+F_v(T−w)−s_v(T−w),
    F'_w = 2F_w−A_vw−B_vw+2G_vw.                    (5.1)

In the deletion row, v retains its name in T−w. In particular, the next
count from switching to w is h'+F'_w−s'_w. No alignment of ports or
deleted labels is discarded.

For s'_w, x cannot start the path; its predecessor is v. Removing x leaves
exactly a cover counted by A_wv. The q row is the inherited insertion
identity applied to T−w. To prove the F row, locate x in the two-slot cover
whose first path ends at w. If x starts its path, placing it in the first
slot gives F_w−A_vw and placing it in the second gives F_w−B_vw. Otherwise
x immediately follows v. Cutting around x gives G_vw in either slot,
with the first two old slots interchanged in the second case. These four
cases are disjoint and exhaustive, including an empty final old slot.

This is an exact one-step repair, not a claim that the displayed finite
list closes under indefinitely many arbitrary old-root changes. Updating
the mixed sidecar itself invokes higher marked cover data. THM-4099 gives
the inherited general gap formulation of that hierarchy.

## 6. Two small grammar probes and their stopping reasons

First restrict every move to R or N starting from (3,1,1,3). It misses
11 and 13 at every order: the first child has H=5, every second child has
H=9, and the third step has H at least 15. All sources are strong and H
strictly increases thereafter. This disproves the initial proposal that
marking only the current root or newborn could retain the full adaptive
target language.

The second probe gives that proposal every rooted version of the actual
four-vertex child as a seed, then permits R/N. Exhaustion under the cutoff
199 visits 352 distinct profiles and misses exactly

    7,21,39,55,67,83,103,113,121,127,139,157,163,173,181,187

among odd values 3,...,199 (with the triangle's count 3 included). Thus
39 is the first missing nongap and 67 the first missing odd prime other
than 7 in this enlarged restricted language. Every retained profile has
a replayed actual named construction. Merging profiles is valid here
because (2.1)–(2.2) are closed maps; Section 4 proves why the same merging
would be invalid for the unrestricted adaptive search.

These exclusions are permanent for the specified languages: Session 150's
strict-growth proof makes cutoff pruning complete. They say nothing about
unrestricted all-order absence of those counts; the inherited adaptive
bank already realizes all the listed nongaps through 199.

## 7. Connections, verification, and the remaining induction obligation

| Source → target | Map and preserved predicate | Lost data / required sidecar | Decisive control |
|---|---|---|---|
| Actual rooted T → four coordinates | Count h and exact keep/new response | Other old-root responses; mixed endpoints and named deletions | The strong adaptive A/B pair |
| Strong T → canonical adaptive parent | Criterion δ=0 and q>s, then delete its named vertex | A target integer alone does not supply T; parent profile needs two-deletion count | Both nonadaptive strong H15 classes at order five |
| Balanced actual state → scaled C3 dynamics | Two actual insertions, then linear R/N maps | Original internal structure; only this operation language is controlled | b=1 exception and b≠7 factor guard |
| Fixed cutoff profile graph → target support below cutoff | Exact R/N transitions and strict growth | No unrestricted root choices and no conclusions above cutoff | Replay all 352 stored recipes |

`adaptive_probe.py` has no repository imports. It constructs every ordinary
tournament through order six by exhaustive one-vertex augmentation and
canonicalization over all permutations within equal-outdegree cells. The
class counts are 1,1,1,2,4,12,56 at orders 0,...,6. Direct endpoint-subset DP
on actual graphs checks 247 strong rooted profiles and 1,188 ordered pairs
of distinct old roots, including (2.1)–(2.3), (3.2)–(3.4), and (5.1).
The legal-inverse recognizer accepts 1,1,4,20 classes at orders 3,...,6,
agreeing with the all-order ancestry theorem. Both earliest rejected
strong order-five classes have count 15: the regular class has no
indegree-one vertex, while C3[C3,1,1] has one with nonstrong deletion.

Reproduce from the repository root:

```sh
python3 -I -B 05-knowledge/results/adaptive_metagraph_reformulation_20261103/adaptive_probe.py
python3 -I -B -O 05-knowledge/results/adaptive_metagraph_reformulation_20261103/adaptive_probe.py
```

Both outputs agree with `adaptive_probe.json`; on the session machine each
run took about 0.5 seconds. Explicit exceptions enforce checks under -O.
The complete ordinary-class universe supports the stated minimality; the
three independent full-permutation counts support the named target ports.
The separate `adaptive_independent_audit.py` imports no producer code and
checks all 285 rooted labelled tournaments of orders one through four by
literal permutations and path covers. It also recomputes both six-vertex
sources and all fourteen eight-vertex children by literal Hamiltonian
words. Its receipt is `adaptive_independent_audit.json`; it audits the
identities and displayed witness, not the producer's minimality census.

The remaining target induction must construct an actual adaptive source
and a named root whose value h+F−s is the prescribed prime or 7p, or prove
a target reduction that supplies such a source. A formal numerical
solution of the profile equations does not suffice. The unique object
inverse, exact reset, and mixed sidecar make this obligation precise;
none eliminates it.
