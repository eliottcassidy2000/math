# Four black moves dissolve an interior R4 block

**PROVED, all orders:** let A,W be arbitrary nonempty tournaments. If
`A -> C3 -> 1 -> W` is black-active, it is connected to `A -> R4 -> W`
by four black moves after choosing the context paths appropriately.
Every arc internal to A,W and every arc between the final components is
restored literally. Converse gives the alternative replacement
`R4 -> 1 -> C3` whenever its target is active.

In particular both replacements are available in the singleton/C3/R4
SCC language whenever their endpoints are active and both contexts are
nonempty. This is a boundary-path theorem, not a count-based inference.
The proof has been independently audited by the metagraph lane.

## 1. A literal four-K identity

Use `K_P` to complement every nonpath arc, retaining the consecutive
edges of a directed Hamiltonian path P. This is the same line as path
reversal after taking the converse-merged quotient. Label the R4 block
by outneighbor bit masks `B=(6,4,8,3)` on0,1,2,3. Keep chosen Hamiltonian
paths of A and W forward and contiguous in each of the following words:

```text
P1 = (A,2,3,0,1,W)
P2 = (W,3,0,A,2,1)
P3 = (W,3,1,0,A,2)
P4 = (A,2,0,3,1,W).
```

All four words are directed paths at their respective steps. Four K
operations give exactly `A -> D -> W`, where `D=(10,0,3,6)` consists
of the cycle `2 -> 0 -> 3 -> 2` followed by sink1. The formal two-port
check retains distinct first and last vertices of each context. It
checks path directions and the edge-XOR identity; internal context
path arcs are retained four times, and every internal nonpath arc is
complemented four times. Thus the identity holds for arbitrary internal
orientations of A and W. Singleton contexts are checked separately.

We prove P3 is always black and choose context paths to make the other
three moves black simultaneously. Throughout,
blue means that positional reversal is an anti-automorphism, so every
paired pair of vertices has total outdegree `n-1`.

## 2. The third move is always black, without a restriction on W

Put `m=|A|>=1`, `k=|W|>=1`, and let `w` be the outdegree in W of the
first vertex L of its chosen path. When k>=2, the endpoint degree
excesses at P2 and P3 are respectively

```text
d(first)+d(last)-(n-1) = m+2-w,       w-m-1.             (1)
```

When k=1 they are instead `m+1` and `-m`, both nonzero. The difference
comes from the first and last W vertices coinciding in the first move.

Assume k>=2 and that P3 is blue. Equation(1) forces `w=m+1`, hence
`k>=m+2`. Immediately before P3, vertex2 beats every W vertex and
vertex1; L also beats vertex1. Positional reversal sends vertex2 to L.
If k>=m+3, the reversed positions of the W block include vertex1:
W occupies positions0,...,k-1 and its reversed interval begins at m+4,
while vertex1 is at position k+1. Anti-symmetry would therefore make
L lose to vertex1, a contradiction.

The remaining case is `k=m+2`. Now positional reversal sends vertex3
to vertex1 and vertex0 to the last vertex R of W. Immediately before
P3 we have `3 -> 0` and `1 -> R`. The first arc would require `R -> 1`
under reversal, contradicting the second. Thus P3 is never blue.

## 3. A blue second move requires a large first component of W

**Stronger necessary condition:** if P2 is blue, the first strong
component of W has order at least `m+4`.

The k=1 case is already excluded. For k>=2 equation(1) forces `w=m+2`.
Immediately after K1 the exact outneighborhood of vertex1 is

```text
N_out(1)=V(A) union {3,L}.
```

The only inneighbors of L outside W consist of vertex1. Under the
positional reversal of P2, vertex1 maps to L; the A vertices map to
positions2,...,m+1 of W, while vertex3 maps to position m+3. All those
positions must lie in W. In particular `k>=m+4`.

The A images must be inneighbors of L after K1. These pairs are not
consecutive path edges, so originally in W, L beats positions2,...,m+1.
It also beats position1 because the chosen W word is a directed path.
Let q be the size of W's first strong component. If q=1, then L dominates
the remainder of W, so `w=k-1=m+2`, forcing `k=m+3`, a contradiction.
There is no strong tournament of order two. If `3<=q<=m+2`, L dominates
every other vertex of its first component, contradicting strongness.

It remains to exclude q=m+3. Immediately after K1, vertex3 beats all
of A. Reversal forces every W vertex at positions2,...,m+1 to beat
the W vertex at position m+3 after K1. Those pairs are nonconsecutive,
so in the original W the vertex at m+3 beats all vertices at positions
2,...,m+1. At q=m+3 it lies in a later strong component yet beats a
vertex of the first component, contradicting the ordered SCC join.
Consequently `q>=m+4>=5`, as claimed.

This proves P2 black whenever W begins in a singleton, C3, or R4. The
size-orientation criterion `m>=k` is an earlier sufficient special case;
it is no longer needed for the three-atom language.

## 4. Choosing the context paths without an SCC-size restriction

First, P2 can always be made black. If some vertex of the first SCC of W
has W-score different from m+2, choose that vertex to start the W path.
Every vertex of a strong tournament starts a Hamiltonian path; append
paths through the later SCCs. Equation(1) then makes P2 black for k>=2;
k=1 is already settled.

Suppose instead every vertex of the first SCC has W-score m+2. A blue
P2 would send each vertex a of A into that first SCC, at W positions
2,...,m+1, by section3. The score of that image in the original W is
the score of a after K1, namely

```text
m-1-d_A(a)+[a is the first vertex of the chosen A path] <= m.
```

Indeed its image is not a global endpoint of P1, so K1 complements its
score, and the anti-reversal complements it back. This contradicts the
assumed common score m+2. Hence P2 is automatically black in this case.

It remains to combine this choice with endpoint blackness.

The source and target cannot both be self-converse. In any self-converse
SCC word, at most one self-converse strong type occurs an odd number of
times: only the possible central component can be unpaired. Replacing
R4 by C3 and1 flips the parity of three distinct self-converse types.
Two vectors each having at most one odd coordinate cannot differ in
three coordinates. Other nonself-converse component types do not affect
this argument.

If neither endpoint is self-converse, P1 and P4 are automatically black.
Use the preceding choice for P2. Suppose the source is self-converse.
When m>=k, (1) makes P2 black for every path choice. If the displayed R4 is the central
component, its path `(2,3,0,1)` is black, so the full P1 cannot be blue.
If it is paired with another R4 in a context, choose the path in that
other copy not to match the displayed path under positional reversal.
R4 has five Hamiltonian paths and trivial automorphism group, so this
choice is available. A blue full word would induce the forbidden match
of the two component paths. When m<k, SCC symmetry has the precise form
`W=Z -> R4 -> A^op`, with Z self-converse, possibly empty. If Z is empty,
the first SCC of W has order four, so P2 is black for every choice by
section3. If Z is nonempty, the mirrored R4 is not the first SCC of W:
choose the starting vertex of W for P2 as above, then change the path
in the mirrored R4 independently. This does not change the W starting
vertex. The target is nonself-converse in all these source-symmetric cases.

Finally suppose the target is self-converse and active. Choose a black
Hamiltonian path of the target. A C3 automorphism aligns its triangle
segment with `(2,0,3)` while fixing A,W and the sink1. Its two context
segments then supply the desired paths for P4. Blue status is invariant
under K because positional reversal preserves the consecutive-edge set;
therefore P4 is also black just before its last K move. The source is
nonself-converse. For m>=k, P2 is black by (1). For m<k, peeling the
mirrored A and the next SCCs from the target gives exactly two forms:

```text
W=C3 -> A^op,            or W=Z -> 1 -> C3 -> A^op,
```

where Z is self-converse and may be empty. In the first case every
possible W starting vertex has score m+1. In the second case its score
is m+3 if Z is empty, and at least m+4 if Z is nonempty. None equals
m+2. Thus P2 stays black for every target-compatible path choice.
The third move is always black by section2.

This proves the four-black theorem. Converse gives its dual statement.

## 5. Adjacent swaps in the rigid language

Any adjacent singleton/C3 interchange between active rigid SCC words,
with both outside contexts nonempty, can pass through an R4 block by
two such four-black routes. At a word end the
[unrestricted-prefix interchange](prefix.md) gives four black moves;
the left end follows by converse. Empty outside words give the same
converse-merged type. Thus every adjacent swap whose rigid endpoints
are active has an authenticated black implementation.

These routes may leave the SCC language at intermediate steps. They
are paths in the ambient black metagraph, not edges of its induced
subgraph on terminal words.

## 6. Fixed-inventory normalization of the entire three-atom language

**Theorem (PROVED).** Fix n and f. All black-active SCC words on
`{1,C3,R4}` with n vertices and f feedback arcs lie in one connected
component of the ambient converse-merged black metagraph. They normalize
to the converse-merged type of

```text
1^(n-3f) -> C3^f.
```

The nonempty active strata have `1<=f<=floor((n-1)/3)`. This does not
assert that different f strata are different components, that the
intermediate tournaments retain f, or that the induced graph on these
terminal SCC words is connected. In fact distinct nonstrong terminal
words have no direct decreasing triangle move between them.

If a word has a singleton atoms, b triangles and c copies of R4, then
`n=a+3b+4c` and `f=b+c`. The latter is exact: feedback arcs minimize
independently in the ordered SCCs; C3 and R4 each require one. A
dissolution `R4 -> C3,1` or `R4 -> 1,C3` preserves n and f at its two
endpoints and decreases c by one.

First remove every R4 while keeping the word black-active. For an
interior R4, sections1–4 give both replacements whenever their targets
are active. If another R4 remains, both targets are active by the path
mass formula `h=5^c` and the blue-count bound in
[the arithmetic packet](arithmetic.md), section3. If it is the last R4,
the targets are rigid words and can fail activity only when palindromic.
The two replacement words cannot both be palindromes: their only
different positions are an adjacent pair with distinct letters1,C3.
Reflection would have to preserve that two-position set, making it the
middle pair, where a palindrome requires equal letters. Thus at least
one target is active.

For an R4 at the right end, use the
[two-K end dissolution](arithmetic.md) if `A -> C3 -> 1` is active.
If that target is pure blue, its rigid SCC palindrome forces
`A=1`, or `A=1 -> C3 -> V` with V a pure-blue rigid palindrome. In
either case A begins with a singleton. The independent
[empty-context boundary repair](moves.md), section10, gives four black
moves to `A -> 1 -> C3`, which is nonself-converse because its first
and last SCCs have different orders. The left-end case follows by
converse. For the word consisting solely of R4, its existing order-four
black line reaches `1 -> C3`. This covers every boundary, so induction
on c removes all R4 atoms.

Both active end replacements are available when needed: if the first
end-dissolution target and the desired alternative are both active,
the unrestricted-prefix interchange joins them. If the first target
is pure, the boundary repair above reaches the active alternative
directly. This justifies both orientations of every active replacement
edge in the abstract macrograph checked below.

The resulting rigid word has exactly `r=n-3f` singleton letters and
`t=f` triangle letters, with r,t positive. Its active words are precisely
the nonpalindromes. By section5, any adjacent swap with active endpoints
has a black path in the ambient graph. Let I count pairs in which a
triangle occurs before a singleton; `0<=I<=rt`. Word reversal sends I
to `rt-I`. Therefore every palindrome has `I=rt/2`.

If I<rt/2, successively swap an adjacent triangle/singleton pair to
decrease I until it reaches zero. All those words remain strictly below
the midpoint, hence are nonpalindromic. If I>rt/2, increase I to rt
instead; those words stay above the midpoint. The two sorted endpoints
are identified by converse. Finally, a nonpalindromic word at the
midpoint has `0<I<rt`, hence an adjacent inversion and a decreasing
move to I-1. That neighbor is active, and the preceding argument applies.
Thus all active rigid words with this inventory lie in one ambient
black component, proving the theorem.

The scalar f is not the directed-triangle potential: the original
word has `c3=b+2c`, whereas the normalized endpoint has `c3=f`.
The endpoint transformation decreases triangle count by the number
of dissolved R4 atoms, but its individual black paths may first
increase it. This is why distinct local minima can share a basin.

## 7. Exact controls and reproducibility

[arithmetic_interior.py](arithmetic_interior.py) checks the formal
two-port ledger and both singleton-port identifications. It then uses
every standard-Hamiltonian-path presentation of A of orders1,...,3
and W of orders1,...,6, for4,400 pairs. It independently enumerates
actual target Hamiltonian paths to determine activity, searches
context-path choices, and checks every resulting path, color and
literal final adjacency matrix. Blue second-move examples are retained
as hostile controls; the checker verifies their first-SCC size bound.
Separate R4/R4 and target-symmetric contexts test endpoint coupling.
An explicit order-five end certificate is replayed independently.

Run both normal and optimized modes:

```text
python3 05-knowledge/results/scc_basins_20261020d/arithmetic_interior.py
python3 -O 05-knowledge/results/scc_basins_20261020d/arithmetic_interior.py
```

The exact universe, source hashes, branch witnesses and labelled paths
are retained in [arithmetic_interior.json](arithmetic_interior.json).
Normal and optimized runs give identical stdout and byte-identical
receipts. The4,400 presentation pairs include4,391 active targets, six
cases needing a context-path change from the standard choices, and
three standard-path blue-second-move hostile examples. The JSON SHA256 is
`58cd0e2eccdc23d75db49ae6926218750d89b0654988e38d441a214391ed4675`.

The independent [abstract-word checker](arithmetic_wordgraph.py)
enumerates all7,920 three-atom words through vertex order18. It retains
3,957 active converse-merged nodes in45 nonempty `(n,f)` strata, builds
the graph of active adjacent swaps and active R4 replacements, and
finds exactly one component in every stratum. It checks interior,
left-end, right-end and whole-word operations separately; orders4 and5
each give the expected two-node macrograph. These are edges representing
the proved black paths, not direct edges in the induced metagraph.
Run `python3 05-knowledge/results/scc_basins_20261020d/arithmetic_wordgraph.py`
and its `-O` variant. Both produce identical stdout and JSON, whose SHA256 is
`c823eccd08828857b442a326329b388670ec4668e0e8b3222b41527d06088752`.

The unbounded conclusions follow from the path and SCC arguments above;
the controls are not an extrapolation from their finite order range.
