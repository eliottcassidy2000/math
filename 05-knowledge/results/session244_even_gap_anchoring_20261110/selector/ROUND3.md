# Efficient endpoint universality with its full state retained

**PROVED construction and identities; finite exact replay.** Construction
proposed by the parent Session244 agent and independently certified here.
Endpoint existence and trivial-root-stabilizer existence are inherited from
THM-5001; the quantitative cost and complete state below sharpen that route.
No claim of literature priority or global optimality is made. Chan–Pak remains
OPEN.

## 1. Inheritance and the new quantitative bridge

Sources:

- `01-canon/theorems/THM-5001-every-positive-endpoint-quotient-has-a-small-strong-realizer.md`;
- `05-knowledge/results/session233_composition_20261109/endpoint/PROOF.md`;
- `05-knowledge/results/session243_prime_gap_bridges_20261110/constructions/GAP_UNIVERSALITY.md`;
- the first-round `proof.md` in this directory for the canonical response word.

The inherited canonical source w(N), for N>=2, has state
`(g(N), epsilon, N,1)`, where epsilon=N mod2. Its graph has
`ceil(log2 N)+1` vertices. Write g=g(N). Append O and then Z, retaining the
old root a during construction. Denote the O-newborn by v and the final
Z-newborn by z. Now designate v as the root.

**Theorem (PROVED).** This strong rooted tournament has full state

    (H,E,U,d) = (g+2N, N, 2g+2N-epsilon, g+N-epsilon).       (1)

Its root stabilizer is trivial, and

    2N+1 <= H <= 3N,
    U <= 4N,
    H <= U <= (4/3)H,
    |V| = ceil(log2 N)+3.                                  (2)

For N=1 use a directed triangle with any root, whose state `(3,1,3,1)`
also agrees with (1) under g(1)=1. It has3 vertices.

### Proof of the entire profile

The total count g+2N is the inherited compact gap family's lower endpoint.
Deleting v leaves w(N).Z, so

    d = g+N-epsilon.

Inside this deletion graph, the endpoint count at z is g-epsilon by the
J transfer formula. Every vertex except z beats v, while v beats z.
Deleting a terminal v therefore gives every Hamilton path of the deletion
graph except those ending at z. Thus

    E = (g+N-epsilon)-(g-epsilon) = N.

For U, the distinguished path from v is either the singleton v, contributing
d, or starts v,z because z is v's unique outneighbor. Removing v in the
second case leaves a rooted cover of the deletion graph at z. Its response
is g+N, again by J. Therefore

    U = d+(g+N) = 2g+2N-epsilon.

The graph is strong by the same reachability argument as the compact pair:
w(N)'s old root reaches all old vertices; z reaches a, every old nonroot
reaches z, and the new sink v reaches z. The old root reaches v directly.

Any automorphism fixing v fixes its unique outneighbor z, and then fixes
z's unique outneighbor a. The remaining original vertices form a transitive
tournament, which has no nonidentity automorphism. Thus the stabilizer of v
is trivial. For N=1 this is immediate in a directed triangle. Consequently
the endpoint-orbit quotient also equals N, as in the inherited interpretation.

Finally 1<=g<=N gives the H and U bounds. The response ratio follows from
`3U<=4H`, equivalent to `2g<=2N+3epsilon`. The word-length count proves
(2). This proves every claimed coordinate, not only the endpoint sum.

### Equality and order comparisons

Within this construction, H=2N+1 iff N is a power of2, including1.
Indeed g(N)=1 iff all letters in its canonical word are O, by (1) of the
first-round proof and the strict positive increment at every Z.
For N>=2, H=3N iff N=2^k+1 for k>=1, using the inherited all-Z equality
classification. N=1 also has H=3N.

THM-5001's prior order bound was `max(5,ceil(log2(N+1))+3)`. The present
bound improves it by one vertex at powers of2 with N>=2, by two at N=1,
and agrees elsewhere. No minimal-order theorem for arbitrary endpoint
realizers is asserted.

Targeted inspection of THM-5001, its proof and response sidecar, and
Session243's compact gap proof found no preceding stated endpoint bound
H<=3N. The old proof's named endpoint examples had H=23,187,461,1573,9899
for N=7,21,37,101,257. Formula (1) now gives respectively

    (N,H,U,d) = (7,19,23,11),
                 (21,61,79,39),
                 (37,107,139,69),
                 (101,289,375,187),
                 (257,771,1027,513).

The verifier directly recounts these graphs, including the257 case.
This is a quantitative improvement in this explicit route, not a claim that
no cheaper endpoint witnesses were known by any other construction.

## 2. Controlled response is useful, but fixed-root continuation still has a trap

Keep the selected root v and append t>=1 copies of Z. The general transfer
identity from ROUND2 gives

    H_t = H+(2^t-1)U-E-(2^t-2)d.

Substituting all four coordinates of (1) simplifies exactly to

    H_t = (2^t+1)(g(N)+N)-epsilon.                       (3)

This identity also holds at N=1 with the triangle convention. It gives a
complete arithmetic inverse for these particular tails: for epsilon=0 or1,
choose t>=1 with `2^t+1` dividing `h+epsilon`, then solve

    g(N)+N = (h+epsilon)/(2^t+1),
    N mod2 = epsilon.                                 (4)

If the quotient is A, testing `ceil(A/2)<=N<=A-1` suffices and is necessary,
using1<=g(N)<=N. There are finitely many t to test since g(N)+N>=2.
The inverse is complete for this family, not universally successful.

There is a sharp arithmetic distinction:

- For even N, (3) is a nontrivial product, so every t>=1 output is composite.
- For odd N, g(N)+N is even, so `H_t = -1 mod 2(2^t+1)`.

Therefore no Mersenne prime target greater than3 can be supplied by these
nonempty fixed-root Z tails. In the odd-N branch, h+1 is a power of2 and
has no odd factor2^t+1; in the even-N branch the output is composite.
For example31 and127 are excluded from this whole tail family at all N,t,
although they are legitimate full tournament counts elsewhere. Exact inverse
testing also excludes13 and173. The explicit rerooted constructions in
ROUND2 overcome both of these latter targets by changing the root profile.

The hostiles do not refute other suffixes, free rerooting, or Chan–Pak. They
show concretely that uniformly low response-to-count ratio plus endpoint
universality is still not an arbitrary anchored translation theorem.

## 3. A linear-cost input for the synchronized endpoint-gap compiler

The sibling `../composition/PROOF.md` gives the exact synchronized pair
rows `9D-3E+11U` and `9D-E+11U`, with common endpoint and deletion coordinate
`3D-E+3U`. Feeding it (1) yields

    lower = 31g+37N-11epsilon,
    upper = 31g+39N-11epsilon,
    common endpoint = common deletion = 9g+11N-3epsilon.

The difference remains2N, both final responses equal their full counts,
and the upper count is at most70N. This is a linear quantitative cost for
that synchronized construction, not a competing claim against the compact
pair's smaller5N cost: the latter lacks the same output synchronization.
The formulas follow by substitution and the sibling compiler's proved
matrix identities. The sibling verifier owns physical checks of that
compiler; this packet independently certifies the input profiles.

## 4. Reproduction and scope ledger

From this directory:

    python3 -I -B verify_round3.py > ROUND3.replay.json
    cmp ROUND3.replay.json ROUND3.json
    python3 -I -B -O verify_round3.py > ROUND3.optimized.json
    cmp ROUND3.optimized.json ROUND3.json

The verifier directly checks129 full graph profiles (N=1,...,128 and257),
strongness, size/equality bounds,16 actual root-stabilizer enumerations,
and96 literal fixed-root Z continuations. It also computes the complete
family image (3) for target values at most2001 using its proved finite
inverse bounds; there are512 such values. No finite prefix is promoted to
an unbounded coverage conclusion. The full theorem follows from the
combinatorial profile proof, not numerical extrapolation.

Connection contract: compact gap source -> efficiently rooted endpoint
witness, map = append OZ then reroot at the O-newborn. Preserved: the physical
full count and strongness. Restored: arbitrary endpoint value with linear
bounds on H,U,d and trivial stabilizer. Lost: independent target placement
of H. Cheapest hostile to a fixed-root repair: prime31, excluded by the
factorization (3). A future construction must change the root menu or use
another operation rather than suppress this retained arithmetic sidecar.
