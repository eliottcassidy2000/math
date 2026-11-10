# Anchored inverse selectors and the exact loss in the compact even-gap family

**Status: PROVED elementary image/selector bounds; FINITE-EXACT seed-bank census.**
Chan–Pak remains OPEN. No full-count nonexistence outside the explicitly stated
construction classes is asserted. The finite bank computation is useful for
classifying this proposed repair, not a record for finite Chan–Pak coverage.

## Inheritance, scope, and research board

Input commit: `c6b359b9c0ca`. The inherited physical state is `(D,E,U,d)`:
full Hamilton path count, paths ending at the root, distinguished rooted
suffix/complement path covers, and root-deletion count. The physical moves are

    O: (D,0,2U,d)
    Z: (D+U-E,d,2U-d,d).

Dependencies: `session199/physical-interpolation/PROOF.md`, especially §§1,4;
`05-knowledge/results/session221_computation_to_construction_20261107/physical/PROOF.md`;
and `05-knowledge/results/session243_prime_gap_bridges_20261110/constructions/GAP_UNIVERSALITY.md`.
The last source already gives the unique binary response selector and the
strong induced-deletion pairs `L(N)=g(N)+2N`, `R(N)=g(N)+4N`, including the
special pair 3,5 at N=1. The all-even gap theorem is inherited, not new.

Anchor: select absolute counts from this gap family. Niche: enumerate *all*
retained-root tails of a finite seed bank by target rather than order.
Wildcard: exploit 2-adic valuation to quantify how much the gap projection
forgets. Board: response uniqueness; endpoint-image redundancy; valuation
counting; unit-deletion rerooting; target-bounded seed banks. The canonical
hostile is a target such as 13, which is a legitimate full tournament count
but impossible in the one-apex subfamily. Corrected near miss: a complete
response selector is not a complete full-count selector. The least-used
sidecar here is the absolute source count retained beside its response.

The packet makes no literature-priority claim. It leaves existing files alone.
The message processor could not run because `.machine-id` is not present;
no identity file was created in this shared worktree.

## 1. Exact recurrence and pair-image redundancy

Put g(1)=g(2)=1, with g(1) an auxiliary convention for the exceptional pair.
For m>=2,

    g(2m)=g(m),
    g(2m-1)=g(m)+m-(m mod 2).                    (1)

The even rule also holds at m=1. These follow directly from O and Z and the
fact that a canonical source with response m has E=m mod 2. The state after
the unique word selecting n>=2 is `(g(n),n mod 2,n,1)`.
Induction gives 1<=g(n)<=n. For odd n>=3, (1) gives

    g(n) >= (n+1)/2.                             (2)

Indeed g(m)>=1 and m mod 2<=1. All g(n) are odd.

**Proposition 1 (exact image, PROVED).** The union of all lower and upper
counts in the compact Session243 gap pairs is exactly

    S = {g(N)+2N : N>=1}.                        (3)

For every N>=1, R(N)=L(2N), since g(2N)=g(N). Thus the upper endpoints add
no new absolute counts. For h>3, an exact inverse for S is

    ceil(h/3) <= N <= floor((h-1)/2),
    N>=2, and g(N)+2N=h.                        (4)

Both bounds follow from 1<=g(N)<=N, and testing all N in (4) is sufficient
and necessary. This is a complete target-bounded selector, not a claim that
the selector succeeds for every h.

A collision family, independently checked through k=30, is

    L(2^k+1)=L(3*2^(k-1))=3*2^k+3,   k>=2.     (5)

All-Z words give g(2^k+1)=2^k+1; stripping O gives
`g(3*2^(k-1))=g(3)=3`. The two inputs are distinct for k>=2.

## 2. Positive-density loss despite every even gap

**Proposition 2 (PROVED).** Define

    c = sum_{j>=0} 1/(2^(j+2)+1).

Then

    limsup_{X->infinity} |S intersect [1,X]|/X <= c,
    c = 0.431166447015110... .                   (6)

Consequently the lower limiting fraction of positive odd integers omitted
by S is at least `1-2c = 0.137667105969779...`.
This concerns only the compact compiler, not all tournaments.

Write N=2^j m with m odd. By (1), g(N)=g(m). Equation (2) also holds for
m=1 under the auxiliary convention. Hence

    L(N) >= (2+2^(-j-1))N + 1/2.                (7)

If L(N)<=X, necessarily `m <= (2X-1)/(2^(j+2)+1)`.
The number of positive odd m satisfying this is at most
`X/(2^(j+2)+1)+O(1)`. Only O(log X) valuations can contribute since
L(N)>=2N+1. Summing yields `|S intersect [1,X]| <= c X+O(log X)`;
counting inputs can only overestimate the image because of collisions.
For a rational numerical enclosure, sum j=0,...,J-1 exactly and bound the
remaining tail above by `2^(-J-1)`. The verifier uses J=16 and records the
exact fractions; this certifies c<0.432 without floating-point dependence.

Thus the even-gap projection is quantitatively inadequate here: it misses
at least a fixed positive proportion of absolute odd targets, even though
it selects every even difference and both counts lie below 5N.

## 3. All one-apex counts have a complete target-bounded inverse

A one-apex tournament is one in which deleting a distinguished vertex leaves
a transitive tournament. Write that chain in its transitive order. Any initial
chain vertices beating the apex are global sources; remove them until the
first chain vertex is beaten by the apex. Removing a global source preserves
the path count: every Hamilton path must start there. If nothing remains,
the original tournament was transitive and its count is 1.

The remaining tournament is exactly a word from the rooted edge: append
chain vertices in order, using O when the apex beats the new vertex and Z
otherwise. Every old nonroot vertex beats the new vertex, as required.
Trailing O vertices are global sinks and may also be stripped without
changing the count. A nonempty reduced word ends Z, so its response U is
odd. Conversely each such word is a physical one-apex tournament.

**Proposition 3 (PROVED).** For every odd h>=3, the following are equivalent:

1. h is the path count of some one-apex tournament, of any order.
2. There exists odd U with `h<=U<=2h-1` and `g(U)=h`.

The inequalities follow from g(U)<=U and (2). This is an exact all-order
classification by a finite target-dependent test. The verifier applies it
through h=2001. For each target, a negative answer excludes the entire
one-apex class, not merely graphs up to a tested order.

In particular, h=13 requires U in `{13,15,17,19,21,23,25}`. The corresponding
counts are `{11,9,17,15,19,15,23}`, so 13 is impossible in this class.
Nevertheless the five-vertex tournament with adjacency masks
`[18,8,3,5,14]` has exactly 13 Hamilton paths, independently counted by
subset DP. Root 0 has profile `(13,2,18,3)`. This explicitly distinguishes
the subfamily obstruction from the desired universal statement.

The condition d=1 is exactly transitive root deletion: a tournament has one
Hamilton path iff it is transitive. To see the nontrivial direction, a
nontransitive tournament contains a directed triangle, with three Hamilton
paths. Inserting each remaining vertex into each path is always possible
(the usual first-insertion argument), and deletion of the inserted vertex
recovers the old path, giving an injection and preserving the lower bound 3.
Consequently rerooting while maintaining d=1 cannot repair the h=13 hole.
One must leave the unit-deletion class somewhere, or use a different physical
operation whose result leaves that class. This does not forbid rerooting
with d>1, J moves, composition, or multi-apex constructors.

## 4. A complete inverse for finite seed banks, with unlimited tail length

Let B be a finite bank of actual rooted tournaments, each with its root
reaching every vertex, at least two vertices, and state `(D,E,U,d)`.
Allow arbitrary O/Z words, retaining the marked root throughout.
All states satisfy `0<=E<=D<=U` and `U>=2d`; d stays fixed.
In particular each letter strictly increases U:

    O: U'=2U;    Z: U'=2U-d >= U+d.

After k letters from a seed with response U0,

    U_k >= 2^k (U0-d)+d.                         (8)

Suppose a nonempty word has full count h. Strip any trailing O letters;
these are sinks and do not affect h. If any letters remain, its last move
is Z, from a state s satisfying

    h = D(s)+U(s)-E(s) >= U(s).                  (9)

**Proposition 4 (PROVED exact algorithm).** To find every count <=h in the
bank with arbitrary O/Z tails, include the seed counts, then enumerate all
reachable states with U<=h and test their final Z counts. Deduplicate equal
four-coordinate states. There are finitely many such states, bounded in
word length by (8). Every permitted count <=h is found by (9); no discarded
branch can later produce it because U is strictly increasing.

This is both an exact inverse and a proof obligation suitable for an anchored
coverage theorem. To use a bank to prove all-odd coverage, prove that for
every odd h outside {7,21}, either a seed has D=h or some state in this
finite target-bounded reachable region satisfies `D+U-E=h`. Merely producing
states with every possible difference does not establish this intersection.
The algorithm retains the needed absolute-count coordinate rather than
subtracting a generated offset algebraically.

### Complete finite-value experiments, not an unbounded coverage theorem

We enumerate every labelled tournament of orders 2 through b, and every
root reaching the entire tournament; deduplicate their exact profiles.
For b=3,4,5 there are respectively 3,10,42 distinct seed profiles. Subset
DP computes D,E,d and a separate reversed-graph DP computes U from its
rooted-path/complement definition, without using the O/Z recurrence.

For every odd target at most 2001, the all-tail inverse gives:

- b=3: 545 covered values; the same image as the one-apex inverse.
- b=4: 789 covered values; first extra holes beyond 7,21 are 55,67,83.
- b=5: 985 covered values; holes through 201 are exactly 7,21,105,173.

There are 1001 positive odd targets <=2001. Seed counts include 1.
These are exact statements about the specified bank family, not all
order-b tournaments after arbitrary rerootings. In particular **no length
of retained-root O/Z tail on any seed of order <=5 can realize 105 or 173**.
That negative assertion is all-order inside this family because the response
cutoff proves exhaustive exclusion. Larger banks, midword root changes,
and other operations are not excluded. The improvement from b=3 to b=5
shows that restoring a non-unit deletion coordinate is genuinely useful;
the remaining holes show why merely changing the initial root is insufficient.

## 5. Reproduction and loss ledger

Run from this directory:

    python3 -I -B verify.py > RESULT.replay.json
    cmp RESULT.replay.json RESULT.json
    python3 -I -B -O verify.py > RESULT.optimized.json
    cmp RESULT.optimized.json RESULT.json

The verifier uses only the standard library and explicit require checks,
so optimization does not remove checks. Deterministic output reports:
99,999 canonical response/state tests; 255 independent literal graph
profile tests; 29 collision identities; exact inverse comparisons through
2001; exhaustive seed graphs/roots of order <=5; and rational bounds for c.
The proofs, rather than a finite replay, justify the global claims.

Connection contract: physical gap pairs -> absolute endpoint image, map
`pair -> {H(left),H(right)}`. Preserved: actual graph realizability, the
prescribed even gap, and strongness of the compact pair outputs. Lost:
independent placement at a target. Restored sidecar: full state plus an
exact target-bounded inverse. Cheapest decisive tests: 11 for the compact
family, 13 for all one-apex graphs, 105 for all <=5 seed retained-root tails.
No subtraction, density claim about the full tournament spectrum, or
all-odd conclusion is licensed by these subfamily results.
