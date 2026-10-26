# Independent audit of the cyclic bridge criterion and order-ten diagnostic

2026-10-26. **ACCEPTED / PROVED:** the sufficient all-size bridge criterion
below. **FINITE-EXACT:** the explicitly bounded controls and two complete
endpoint-fixed word censuses. This is not a necessary criterion, a theorem
that a qualifying word always exists, or a black-core connectivity theorem.

## 1. Exact bridge theorem and proof

Let n>=3. Let P=(0,1,...,z), z=n-1, be a directed closed blue Hamiltonian
word in T. Put k=d_T^+(z), I=N_T^-(0) minus {z}, and B=K_P(T). Choose a
directed Hamiltonian word W in B minus vertex0, with first vertex1 and last
vertexz. Define

    C=(1,2,...,z,0),   S=K_C(T),
    Q2=(W,0),         Q3=(0,W),
    g=r_Q3 r_P,       J=I union {1}.

If d_T^+(1) != k and either g has even order or g(J) != J, then C,Q2,Q3
are three successive directed black words and their K moves reach the exact
labelled target B.

**Legality and target.** C is directed because P is directed and z->0. Its
edge set differs from P exactly on the two pairs01 and0z. Consequently S
and B agree on every pair not incident to0, so W is directed in S as well.
S retains z->0, making Q2 directed. In U=K_Q2(S), every W edge is retained.
The pair01, omitted by both C and Q2, has been reversed twice and again has
direction0->1. Thus Q3 is directed. Writing E(W) for its interior path pairs,

    E(C)=E(P) XOR {01,0z},
    E(Q2)=E(W) union {0z},
    E(Q3)=E(W) union {01}

gives E(C) XOR E(Q2) XOR E(Q3)=E(P). Three complements therefore produce
exactly K_P(T), with every vertex and arc retained literally.

**First and second blackness.** Anti-reversal pairs a word's endpoint
outdegrees to a sum n-1. Original P blue gives d_T(0)=n-1-k. Therefore blue C
would imply d_T(1)+n-1-k=n-1, or d_T(1)=k. At S, the first vertex1 gains the
usual path-endpoint correction and the last vertex0 loses it:

    d_S(1)=n-d_T(1),       d_S(0)=k-1.

Blue Q2 would similarly require n-d_T(1)+k-1=n-1, giving the same excluded
equality. Hence the one degree guard certifies **both** initial words.

**Third blackness.** A positional reversal preserves its own unordered path
edge set, so blue/black status agrees across that word's K edge. Q3 may
therefore be tested at B. P remains blue in B. It follows that Q3 is blue
in B if and only if g=r_Q3 r_P is an automorphism of B: composing two anti
maps gives an automorphism, and composing that automorphism with r_P gives
the converse. Because Q3 and P have the same endpoints, g fixes0 andz.
The exact outgoing neighborhood of0 in B is I union {1,z}, so g must
preserve J=I union {1}. If it does not, Q3 is black.

An even-order tournament automorphism is impossible. An appropriate power
would be a nonidentity involution, which swaps some pair and would preserve
and reverse its directed arc simultaneously. Therefore even order of g
also certifies Q3 black. This proves the theorem.

This proof works for arbitrary endpoint degree k and imposes no clone,
module, strongness, mixedness, or disjoint-port assumption. It retains the
actual J, rather than just its cardinality. Neither existence of a
1-to-z Hamiltonian word nor satisfaction of the blackness guards is proved.

## 2. Independent finite controls and two equality/menu hostiles

`bridge_check.py` imports only the local independent `core.py`, copied from
the sealed first audit. It reconstructs all closed blue marked sources at
n=3,...,7 without an endpoint-degree restriction, then tests every
endpoint-fixed W permutation. This gives source counts 1,2,8,32,256. It checks
all legal bridge targets, actual degree identities, the exact blue/automorphism
equivalence, endpoint-port preservation, and each sufficient implication.
There are 151 sufficient-gate hits, all giving actual three-black ledgers.
These finite checks corroborate, rather than replace, the proof above.

The copied 27 original repair certificates contain26 cyclic bridges. All 26
have g(J)!=J;25 satisfy the degree guard, and ten also have even-order g.
All 26 are literally three-black, including the one at degree equality.

For the equality source

    T=(126,156,440,1008,1248,1474,1926,1793,1555,1079,15),

d_T(1)=k=4. Exhausting all 8!=40320 endpoint-fixed W permutations yields132
directed W. C and Q2 are black for every one. Exactly 127 choices also have
Q3 black; the remaining five have Q3 blue. This source is the **cyclic**
repaired source with the transient ordered clone block (5,4), not the source
whose supplied certificate was noncyclic. Thus degree equality is not a
connectivity obstruction or even an obstruction to this bridge construction.
The degree condition is sufficient, not necessary.

As a separate control, the order-eleven source whose discovered certificate
was noncyclic is

    (126,220,472,976,1952,1102,1936,1825,1571,1063,15).

The same complete 40320-permutation census finds122 directed W:119 give
three black moves and three have Q3 blue. A noncyclic *discovery certificate*
was therefore not evidence that cyclic bridges were absent. The two sources
and their different word counts are kept separate in the receipt.

## 3. Bounded independent order-ten dimensions

`n10_check.py` is an independent signed-component implementation, not the
Gaussian elimination used by the producer's diagnostic. It reuses only
the explicit word formulas and literal replay helpers. It constructs the
raw source and blue equations itself, with all actual preceding flips.

All 140 marked endpoint-four source cells at n=10 have dimension10, giving
143,360 labelled marked sources. The130 base cells and ten endpair cells
refine to151 nonempty, disjoint and exhaustive affine leaves. Their partition
cardinalities and pairwise disjointness are checked. Every old word remains
uniformly directed and every ledger retains the exact labelled target.

Exactly 142 leaves are uniformly three-black. The nine remaining leaves each
have one dimension-zero blue locus, and the nine decoded literal sources are
distinct. Their source arcs, ports, old blue flags, and literal old targets
are replayed. This diagnostic does no replacement search and asserts no
connectivity obstruction. A separate later repair certificate may settle
those nine sources; their failure under this old menu alone does not do so.

## 4. Reproduction and limits

    python3 bridge_check.py > bridge.normal.json
    python3 -O bridge_check.py > bridge.optimized.json
    cmp bridge.normal.json bridge.optimized.json
    python3 n10_check.py > n10.normal.json
    python3 -O n10_check.py > n10.optimized.json
    cmp n10.normal.json n10.optimized.json

Both normal/optimized pairs pass and are byte-identical. Every gate is an
explicit exception. The manifest pins the local independent core, both
checkers, copied literal certificates, all four receipts, and this report.
No whole black-core conclusion, optimal threshold, or all-order W-existence
claim is authorized by this packet.
