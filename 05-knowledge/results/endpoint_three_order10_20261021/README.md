# The remaining order-ten endpoint-three case

October 20, 2026. **FINITE-EXACT, complete affine certificate.** No tournament
census is used. The endpoint-three theorem's order-ten boundary is closed by
168 marked affine cells and three explicit exceptional ledgers. Together with
the inherited structural theorem for every n>=11, this gives the same
unfiltered three-black detour conclusion for every n>=10. The proof here is
finite and exact at n=10; it is not a replacement all-order structural proof.

## 1. Statement, inheritance, and limits

Let T be a tournament on ten labelled vertices and let P=(a,...,z) be a
specified directed Hamiltonian word. For a directed word Q, let K_Q(T) retain
the arcs of Q and reverse every other arc. A word is **blue** when its actual
positional reversal is an anti-automorphism of its source tournament, and
**black** otherwise. Assume P is blue, z->a, and

    d_T(a)=6,    d_T(z)=3.

Then three successive directed black K moves end at the **literal labelled**
tournament K_P(T). There are no source-mixedness, target-mixedness,
strong-partner, or isomorphism filters. In particular the two marked endpoints
belong to the same black component. This does not prove connectivity of the
whole black core, a global bridge classification, or any unrelated endpoint
layer. Reading this marked line from its open partner gives endpoint scores
(4,5), but it is still the same line.

The closest mechanism is the complete n>=11 partition in
`endpoint_three_repairs_20261020h/README.md`, with the eleven ledgers implemented
in its `root_check.py`. The hostile is the sharp order-ten cyclic-cut example,
not a hypothetical bridge. The corrected near misses are failure to preserve
physical port labels, mistaking a failed sufficient score test for actual
blue colour, and confusing K-complement with full tournament converse. The
least-used sidecar is the affine locus on which a particular intermediate
word is actually blue.

The inspected inherited commit is
`c4a7875652a5ea516e93a31c3c6410d2204773b2`. Exact source SHA-256 values are in
`source_pins.json`, including the endpoint-three notes and checker, the marked
trijection audit, and the correction ledger. This packet imports none of the
inherited Python implementations. Its small ledger selector is a frozen copy
of the predecessor's formulas, checked independently rather than granted its
n>=11 blackness conclusions at n=10.

Concept board: mirrored exceptional ports; ordered-graph chart coordinates;
path-edge XOR; actual anti-reversal equations; and singleton template failures.
The mechanism is a parameter-space reduction, not a class enumeration.

## 2. Exact source coordinates and the trijection

Relabel by P, so P=(0,1,...,9), a=0, b=1, x=2, y=7, c=8, z=9.
There is a unique two-element set I contained in {2,3,...,8} with

    N^-_T(0)={9} union I,
    N^+_T(9)={0} union J,    J={9-i : i in I}.

Port overlaps are permitted and retained. The middle on {1,...,8} has the
standard directed path and reflection anti-automorphism i->9-i. Every
nonadjacent middle pair is free once per reflection orbit; adjacent pairs
are fixed forward. There are twelve free orbits, in the following order:

    (1,3), (1,4), (1,5), (1,6), (1,7), (1,8),
    (2,4), (2,5), (2,6), (2,7), (3,5), (3,6).

Each listed pair denotes its orbit under (i,j)->(9-j,9-i). Variable x_k=1
means the listed ordered arc points forward. These are distinct representatives
of all twelve orbits. Every assignment gives a valid closed blue source;
conversely every source gives exactly one assignment and one I. Thus this is
an exact parametrization of 21*2^12=86,016 **marked labelled presentations**,
not unlabelled tournament classes. The closing arc itself makes T strong;
no additional connectivity assumption is imposed on its partner.

In the ordered-graph trijection, a backwards nonadjacent source arc j->i is
the graph edge {i,j-1} on ordered vertices {0,...,8}. Therefore the twelve
variables above are complements of twelve graph-edge orbit bits, while the
port incidences fix the remaining boundary bits. K_P is graph complement
in this chart, and blue P is invariance under the specific order reflection
s->8-s. Changing the Hamiltonian word changes the chart. The three words
below keep this chart-change sidecar; no arbitrary graph relabelling or cut
aggregate is substituted for a legal tournament move.

## 3. The complete affine reduction

For each I, split by the three actual arc values

    T(1,8), T(2,8), T(2,7).

They are three distinct free bits, so the eight cells are disjoint, nonempty,
and each has 2^9 sources. There are exactly 168 cells. Within each cell the
inherited eleven-branch selector chooses one fixed ordered triple (Q1,Q2,Q3).
One representative per cell is constructed only to select this fixed branch.
The certificate stores all three words for every cell, including every port
label. No interior orientation or port overlap is discarded.

Write any directed source-arc indicator as

    T(u,v)=linear_form(x_0,...,x_11)+constant over GF(2).

After a K word Q its constant toggles exactly when {u,v} is not a path edge.
Consequently every intermediate arc remains affine. For each of the 27 joins
in the three words, reducing its required value by the three cell equations
proves it directed throughout that cell. This gives 4,536 symbolic join
checks, not 86,016 source replays.

For each word Q in its actual intermediate source S, let r_Q be the literal
positional reversal. Q would be blue precisely when all 45 equations

    S(u,v)=S(r_Q(v),r_Q(u)),    0<=u<v<=9,

hold. Adjoin the three cell equations. There are 504 such affine systems.
Exactly 501 are inconsistent. Each inconsistent system is accompanied by a
list of equation indices whose XOR is the impossible equation 0=1. The
other three systems have rank twelve, each with exactly one solution. For
each, the certificate supplies twelve XOR proof vectors deriving
x_k=the specified bit, and verifies that assignment satisfies every equation.

All three words also satisfy

    E(Q1) XOR E(Q2) XOR E(Q3)=E(P),

where these E are **undirected labelled** path-edge sets. Three K complements
therefore give the exact target K_P(T), not merely an isomorphic object. Both
checkers independently verify the equivalent affine arc identity for every
cell. The third word is tested at its actual third source, not merely at a
score proxy or an unrelated representative.

The only feasible blue-step systems are:

| I | cell values (T(1,8),T(2,8),T(2,7)) | inherited branch | blue step | unique mask |
|---|---|---|---:|---:|
| {3,6} | (0,1,0) | closed_epsilon_one | 3 | 2202 |
| {5,6} | (0,0,0) | closed_epsilon_zero | 3 | 2446 |
| {7,8} | (0,0,0) | closed_endpair_yx | 2 | 3525 |

A mask encodes sum x_k*2^k in the explicitly listed orbit order. These are
three failed template presentations, not three failed detours or bridges.
The first two failures share the common final word (0,9,1,2,3,4,5,6,7,8).

## 4. Three literal replacements

Adjacency rows below mean that row u has bit v set exactly when u->v.
`exceptions.json` retains these rows, their literal targets, both port sets,
and all directed words. They are sufficient fixed certificates: their
discovery searches are not dependencies of the result.

### I={3,6}, J={3,6}, mask 2202

    T=(438,724,808,339,676,842,149,812,594,73).
    Q1=(0,1,2,3,4,5,8,6,7,9)
    Q2=(0,1,2,3,5,8,6,7,9,4)
    Q3=(0,1,2,3,5,6,7,8,9,4)

### I={5,6}, J={3,4}, mask 2446

    T=(414,628,616,466,420,841,657,806,582,25).
    Q1=(0,1,2,3,4,5,8,6,7,9)
    Q2=(0,1,2,3,9,8,7,4,6,5)
    Q3=(0,1,2,3,9,7,4,6,8,5)

### I={7,8}, J={1,2}, mask 3525

    T=(126,44,120,752,994,704,898,775,559,7).
    Q1=(0,1,2,3,4,6,7,8,5,9)
    Q2=(0,1,2,3,8,5,9,4,6,7)
    Q3=(0,1,2,3,8,9,4,5,6,7)

At each intermediate source, the direct replays exhibit a literal ordered
pair (u,v) whose arc value differs from its required anti-image
(r_Q(v),r_Q(u)). All nine words are therefore actually black. The full
witness tuples are in the receipts. The first source's three endpoint score
sums all equal nine, illustrating why endpoint scores alone cannot certify
blackness. Every replacement is directed and ends at the exact labelled
partner.

**Conclusion.** All but these three source points have the inherited three
black words; these three have the displayed replacement words. The partition
is exhaustive and no filter is applied, proving the stated order-ten theorem.

## 5. Controls, audit, and reproduction

`check.py` uses a small GF(2) elimination routine with equation provenance to
produce `certificates.json` (about 155 KB). It also verifies every provenance
vector, all symbolic joins, all symbolic targets, and the three fixed repairs
with directed arc sets. It imports only the Python standard library.

`verify_certificates.py` is independent of that implementation: it imports
neither `check.py` nor any inherited kernel, reconstructs the source forms by
canonical middle-pair orbits, and verifies the supplied XOR vectors directly.
It performs **no Gaussian elimination and no search**. It proves that every
one of the 168 cells occurs exactly once, checks all 4,536 joins, verifies the
501 contradictions and all 36 unit equations for the three singleton systems,
and replays the literal repairs using adjacency bit rows. A changed or omitted
case, inconsistent proof vector, wrong unit assignment, invalid directed join,
blue replacement, or wrong labelled target raises an explicit exception.

The sharp cyclic-cut hostile is also retained:

    (190,44,760,176,994,704,907,770,559,11).

Its cyclic cut (1,2,...,9,0) is genuinely directed and blue. The third
replacement ledger above repairs it as well. Thus the method's sharp n=10
blackness boundary is not mistaken for a metagraph obstruction.

Run from this directory, with Python 3.10 or later:

    python3 check.py --output receipt.normal.json
    python3 -O check.py --output receipt.optimized.json
    cmp receipt.normal.json receipt.optimized.json
    python3 verify_certificates.py --output independent.normal.json
    python3 -O verify_certificates.py --output independent.optimized.json
    cmp independent.normal.json independent.optimized.json

Both receipt comparisons were byte-identical. Optional inherited-source
verification is available with `python3 check.py --repo /path/to/math`.
That flag reads but does not execute or modify the pinned source files.
Every mathematical check uses an explicit exception and remains active
under optimization. Receipt hashes pin the checker, exception data, inherited
source manifest, and complete affine certificate.

The scope is deliberately finite: this closes n=10, extends the inherited
unfiltered endpoint-three range to n>=10, and supplies a reusable exact
residual-isolation technique. No enumeration of all order-ten tournaments,
all marked Hamiltonian words, or the general black metagraph is claimed.
