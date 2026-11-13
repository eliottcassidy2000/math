# Growing the determinant modules reaches four exact prime exits

**FINITE-EXACT:** the complete specified one- and two-vertex extension
universes, six response pairs and four prime exits. **PROVED:** the
induced commuting construction and the scoped automorphism statements.
**OPEN:** unbounded determinant renewal and a universal prime/7p selector.
These are new witnesses for an inherited constructor, not new global
spectrum records or an assertion of minimum tournament order.

## 1. The response equation and the changed modules

The [polarized repair](../target_collatz_transitions_20261113/POLARIZED.md)
uses arbitrary nonempty modules M,N with responses

    h=H(M), p=P_M(2), j=H(N), q=P_N(2).

Its original count, prepared count, and balanced exit are

    A=12pq+6hq+2pj−hj,
    Z=12pq+2hq+6pj−hj,
    Y=3Z+2hj=3A+2e,
    e=j(6p+h)−6hq.

Three old apex arcs change, then two root-child births give the balanced
state(Y,Z,Z,Y). The prepared root keeps its unique predecessor and full
old star. A prescribed e=±1 gives an exact3A±2 transition.

The old nontrivial pair was

    M0=(4,9,82,21,35,143,187,31), (h,p)=(181,3424),
    N0=(6,12,48,21,131,219,159,15), (j,q)=(179,3416).

Its exit was435779329=3·145259777−2. The fixed-partner and padding
obstructions did not exclude adding vertices that create new strong
modules. That is the operation tested here; no old module arc changes.

Every single-vertex extension of M0 and of N0 was tested:256 orientations
per side,242 and239 distinct response pairs. None of their ordered
pairings has e=±1. This is a statement about adding **one to each**, not
a prohibition on unequal extensions or arbitrary repairs.

With **two vertices added to each**, all131072 labelled extensions per
side yield56614 and58244 distinct profiles. Exactly six ordered response
pairs satisfy e=±1. They are:

| (h,p) | (j,q) | e | original A | exact exit Y |
|---|---|---:|---:|---:|
| (4553,101144) | (4231,94696) |−1|118358708801|355076126401|
| (4577,103880) | (3665,83792) |1|107497517119|322492551359|
| (4811,104600) | (4525,99136) |−1|128222027201|384666081601 **prime**|
| (5305,112888) | (4583,98288) |−1|137285354561|411856063681 **prime**|
| (5815,117688) | (7687,156856) |1|228758168383|686274505151 **prime**|
| (6491,129368) | (6637,133384) |−1|213936379073|641809137217 **prime**|

All twelve displayed modules are strong. All six exits have gcd(Y,Z)=1.
The first two targets have factors19 and11 respectively. Primality of the
other four is checked by trial division through their integer square
roots; no probable-prime test is being called a proof.

## 2. Complete finite search, and independently reconstructed graphs

For base rows T on0,...,7, two added vertices8,9 are specified by masks
a,b in0,...,255 and e0 in{0,1}. Bit i of a or b means i points to the
corresponding new vertex; e0 means8→9. Thus

    row_i = old_i + bit_i(a)·2^8 + bit_i(b)·2^9,
    row_8 = (255 xor a)+e0·2^9,
    row_9 = (255 xor b)+(1−e0)·2^8.

These2·256² choices enumerate every induced two-vertex extension of the
specified labelled base, without an isomorphism or strongness filter.
`pair_growth2.json` supplies one actual attachment witness for every
successful response pair. Representatives are enough because the
determinant depends only on the responses; the result counts profile
pairs, not all labelled graph pairs.

The C++ search computes H and P2 by exact subset endpoint DP. It solves
the pairing equation by the necessary guards gcd(h,p)=1 and3∤h, then
the modular inverse

    j(6p+h)≡±1 (mod6h).

Every positive j up to the maximum observed right response is tested,
and q is recovered exactly before dictionary lookup. This retains every
possible e=±1 pair; it is not a nearest-ratio or sampled search.

The Python checker independently reconstructs all twelve graphs, recounts
their H and P2, and checks the tournament and strongness conditions. It
also enumerates Hamilton permutations for the third row's two modules.
Their26-vertex sources and28-vertex outputs are certified by the inherited
universal module formula, not by a claimed full28-vertex path census.

## 3. A physical commuting square, retaining the old balanced port

Use the inherited named ordinal spine

    r → s → M → t → N → u → v,

and apex a, which beats r,t,N before preparation and r,s,N,u after it.
Extend M0,N0 by the two specified vertices each, leaving all old names
and the other incidences unchanged. This embeds the old22-vertex source
in the new26-vertex source. Applying the same three arc reversals on
both levels preserves the induced embedding.

In the prepared source add x with r→x→all other old vertices; mark x.
Then add y with x→y→all other vertices, retaining x. Identifying the old
and new x,y embeds the original24-vertex balanced output as an induced
subtournament of each28-vertex output. Every old arc and the entire old
star at x survive. Its sole predecessor remains the physical vertex r.
The newly enlarged marked state is the specified(Y,Z,Z,Y), rather than
an unspecified witness with the same scalar count.

This gives four actual induced extensions from the old balanced output
to prime counts, while retaining its controller. It does not show that
one of these new pairs has another compensated extension, or that a new
source count equals the preceding exit. Those are separate return tests.

## 4. The visible balanced interface also controls global automorphisms

The strong-parent stabilizer law is inherited from
[ANCESTRY A1–A2](../adaptive_metagraph_reformulation_20261103/ANCESTRY.md).
Here is a direct argument sufficient for these outputs, with the boundary
made explicit. It applies to any rooted T with indegree(r)≥1, even when
T is not strong.

In NR(T,r), exactly x and y have indegree one. Every old vertex other
than r receives from both new vertices; r receives from y and at least
one old predecessor. Their arc x→y distinguishes the two members of this
canonical pair. Every automorphism therefore fixes x,y, then fixes r as
the unique predecessor of x. Restriction and extension give

    Aut(NR(T,r)) ≅ Aut(T,r).

Using the universal rooted maps, the output state is

    (H+2F,H,H,H+2F).

A root-fixing automorphism acts freely on Hamilton words and on the
ordered marked two-path covers counted by F. Thus gcd(H,F)=1 makes its
stabilizer trivial and makes the balanced output primitive. In the six
prepared sources H=Z,F=Z+hj, and gcd(Y,Z)=1 proves this guard. Hence all
six outputs are globally asymmetric, including the two composite ones.

The indegree boundary is real: if r has indegree zero, the same operation
gives C3→(T−r), with group C3×Aut(T−r), and count3H(T). For example
NR on the singleton creates C3. Conversely, starting at a rooted C3
has trivial point stabilizer and gives the rigid five-vertex H9 graph,
but its balanced marked content is3. Global rigidity alone is therefore
not the primitive-interface criterion.

Nor is NR itself a repeatable primitive return: applying it immediately
to a balanced state(Y,Z,Z,Y) gives(3Y,Y,Y,3Y). A genuine repair or changed
port is needed before a subsequent useful balanced exit. The current
ports and clock lanes supply other return mechanisms; their entry must
still be established rather than inferred from asymmetry.

## 5. The current handoff test

The clock lane gives an actual indefinite changing-module rule, with an
exact finite selector for its A-seed word language. Its complete test of
these six targets, allowing every initial left transitive size, finds no
entry (380 bounded word states). See[clock PROOF §5](clock/PROOF.md) and
the supplied entry receipt. This does not exclude the arbitrary-seed rule
or another clock. It prevents a false composition of two positive results.

Reproduce the extension search with the C++17 source `pair_growth2.cpp`;
the command and exact controls are recorded in the packet audit. Use
`pair_growth.py` for the one-birth universe and `pair_check.py` for the
independent witness, embedding, prime and automorphism-boundary controls.
