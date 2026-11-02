# Compressed exact descent and the minimum legal construction clock

**PROVED for the declared singleton/D grammar; FINITE-EXACT implementation
checks. The scalar-target entry theorem remains OPEN.**

## 1. Inherited language and source guards

Use the actual (F,q,g,h) state defined in [AMPLIFIERS.md](AMPLIFIERS.md), starting
at (1,0,0,1): one marked vertex with empty parent. Append either a singleton
or D, the strong four-vertex tournament with rows (8,1,3,6), after the current
parent. Each cut says precisely which new vertices are beaten by the mark.

The singleton actions are

`S0(F,q,g,h)=(F+q+g,2q,0,h),`

`S1(F,q,g,h)=(F,2q+h,h,h).`                                    (D1)

For a cut m of D, let (5,32,p,Q,s,t,F_B) be its inherited exact profile. Then

`D_m(F,q,g,h)=(5F+pq+(5-s)g+(F_B-5)h,32q+Qh,(5-t)h,5h).`     (D2)

Here s counts full D paths starting in the cut, t those ending outside it,
p counts ordered covers with first slot nonempty and ending outside the cut,
and Q counts covers with second slot nonempty and beginning inside the cut.
The profiles are computed literally by `diamond_profiles.py` in packet31b.
For masks m=0,...,15, the Q digits are, in order,

`0,5,6,11,9,14,15,20,7,12,13,18,16,21,22,27.`                 (D3)

They are distinct modulo32. This decoder and the original terminating
one-step inverse are inherited, not new results of this packet.

Every reachable state obeys h=5^r, F odd, F>=h,

`q>=F-h+g,    q-F+h-g even,`

with g in {0,h} if h=1, or in {0,h/5,...,h} otherwise. These necessary
conditions are **not** sufficient. The positive cut/interior proof of these
guards is inherited from31b. At h=1 the exact state is

`(F,q,g,h)=(f(q),q,q mod2,1),`

where f(0)=f(1)=1, f(2a+1)=f(a), and
f(2a)=f(a)+a+(a mod2) for a>=1. Binary digits of q, read most significant
first as singleton cuts, construct it. Leading S0 operations are identities.

## 2. Jumping over a whole singleton suffix

**Theorem D-jump.** Full-state membership in this language, together with
the smallest number of added vertices, is determined by the following finite
recurrence. At fixed r=log_5 h its running time is polynomial in binary input
length. No claim of polynomial time uniform in r, or for an unspecified
scalar target F, is made.

Suppose r>0. In any nontrivial word isolate its last D and the following
singleton suffix. Let its length be l and its chronological binary value y.
Writing q_D for the slope just after the last D gives

`q=2^l q_D+h y,  0<=y<2^l.`                                   (D4)

Because h is odd, the target and l determine exactly one candidate:

`y=(q h^(-1)) mod 2^l,   q_D=(q-hy)/2^l.`                     (D5)

For l=0 take y=0. Reject negative q_D. Set h0=h/5. The last D digit is

`Q=(q_D h0^(-1)) mod32,  q0=(q_D-Qh0)/32.`                    (D6)

By (D3), there is at most one matching D cut; reject a missing digit or
negative q0. Its output g_D=(5-t)h0 is now known. Replay the l singleton
bits on the formal state (0,q_D,g_D,h). This determines an F increment I,
the final q, and the final g. Require the last two to match the target.
For each allowed previous g0 recover

`F0=[F-I-pq0-(5-s)g0-(F_B-5)h0]/5.`                           (D7)

Require integrality and all source guards, then recurse at (F0,q0,g0,h0).
If h0=1, the only candidate g0 is q0 mod2. Add cost4+l to the best successful
predecessor. Reconstruct the actual chronological word, not merely its cost.

**Why the search is finite and complete.** If q=0, the guards force F=h,g=0;
D_0 repeated r times is a minimum word of cost4r. Otherwise, delete every
S0 while q=0; it is a state identity. After q first becomes positive every
singleton at least doubles it. Thus the total number of remaining singleton
letters is at most bit_length(q), so it suffices to try

`0<=l<=bit_length(q).`

If q_D=0, a suffix with leading zero digits is redundant and is discarded.
Every normalized word has a last D and one of these suffixes. Equations
(D5)-(D7) recover its genuine predecessor, while every successful recursive
branch replays to the target by the literal construction rules. The terminal
h=1 case is exact. This proves both directions by induction on r.

All words realizing a state have exactly r D atoms. The recurrence compares
every possible normalized last block, using its actual vertex cost, and the
base word has minimum singleton cost. Induction therefore proves the minimum
cost assertion. Neutral padding cannot improve a minimum. It can be retained
as a separate positional sidecar if reconstructing a larger supplied object.

At each level there are at most6(L+1) predecessor branches, L=bit_length(q).
Predecessor slopes never increase and the depth is r. A safe recursive-node
bound is (6(L+1)+1)^r. Suffix replay, exact arithmetic and modular inversion
take polynomial bit cost per node. Recursive state coordinates have
input-bounded bit length: F,q,h do not increase backward and g<=h. Temporary
suffix arithmetic also has polynomially bounded length. This proves the
fixed-r complexity assertion. It is an improvement over enumerating all
2^l suffixes, not a new complexity result for arbitrary tournaments.

## 3. An admissible clock bound and least-resistance recursion

Set R=q+h. A singleton gives R'<=2R. A D cut gives

`R'=32R-(27-Q)h<=32R.`                                       (D8)

Starting from R=1, any word with r D atoms and s singletons satisfies
R<=2^(5r+s). Its added-vertex cost therefore has the lower bound

`4r+max(0,ceil(log2(q+h))-5r).`                               (D9)

Every reachable state also has a normalized word of cost at most
4r+bit_length(q), by the proof above. The lower bound is admissible for
branch-and-bound pruning; it need not be attained.

The repaired least-resistance rule is the Bellman equation

`cost(z)=min_{legal predecessor y -> z}(cost(y)+cost(y -> z)),`

with unreachable states assigned infinity. "Legal" here is recursive
membership, not positivity, congruence satisfaction, or low current rank.
This construction is acyclic after removing identities because F+q+h
strictly decreases under every retained inverse step.

The old square25 observer has33 states under its **fixed named square-path
complements and observation**. It shows why a cost should be attached to a
well-defined operation graph. It neither supplies our state coordinates nor
proves closure under arbitrary metagraph operations. Likewise the older
Paley "least resistance" work warned that series/parallel reduction stops at
irreducible character data. Here the missing data are actual predecessor
membership and the overwritten endpoint amplitude g.

## 4. Exact greedy failures

Define the hostile algorithm precisely: at each inverse step choose the
admissible predecessor minimizing F+q+h, breaking ties lexicographically;
do not backtrack. Necessary guards from section1 are still enforced.

| Target state | Least-rank predecessor | Why it fails | Minimum word | Actual vertices |
|---|---|---|---|---:|
| (17,12,0,5) | (7,6,4,5) | forbidden H=7 | D_2 S0 | 6 |
| (49,48,0,5) | (21,24,4,5) | forbidden H=21 | D_9 S0 S0 | 7 |
| (19,14,0,5) | (9,7,3,5) | impossible joint state, attainable scalar9 | D_8 S0 | 6 |
| (15,10,0,5) | (5,5,5,5) | succeeds but spends an extra vertex | D_1 S0 | 6 |

For the last row greedy uses D_0 S1 S0, on7 vertices. Adjacency-DP verifies
all four displayed actual outputs; minimality here is within this grammar.
These are small explicit witnesses, not claimed first failures in an
exhaustive ordering of all graphs.

The third row cannot be repaired merely by blacklisting7 and21. Indeed any
parent with h=5 has one strong component D and otherwise singleton SCCs:
nontrivial SCC path counts multiply, the order-five strong floor is9, and
the only strong tournament of count5 is D. Hence the grammar covers **all**
marked parents with h=5. For q=7, a positive slope before D would give
q>=32. Thus q=0 before D. A positive Q other than7 cannot reach7 by later
doubling and additions of5; Q=0 permits singleton slopes0,5,10,15,... but
not7. Digit7 uniquely selects D_8, which gives (11,7,1,5). Every subsequent
singleton changes q to at least14. Therefore (9,7,3,5) is globally unrealizable
as a marked-parent state, even though tournaments with H=9 exist.

## 5. What the inverse does not prove

The inherited one-step inverse already terminates; this packet compresses
its singleton stretches and minimizes construction cost. It does not prove
that every odd F outside{7,21} occurs in this grammar. The inherited scalar
entry bound and the known grammar misses remain in force. Nor do these four
scalars preserve the full automorphism group, arbitrary selected ports or
literal complemented relations: the81/27 hostile still rules out such a
promotion. The actual word and graph are returned precisely so those future
questions can use the missing incidence and transport information.
