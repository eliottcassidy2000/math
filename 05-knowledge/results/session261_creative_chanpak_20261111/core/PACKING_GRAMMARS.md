# Packing budgets are finite resource grammars

**PROVED.** This is a structural reduction for bounded cycle packing, not a
proof of Chan–Pak or a claim that any fixed packing budget suffices.

The inheritance is Session124's triangle-core forbidden-word mechanism,
`05-knowledge/results/session124_cycle_packing_obstructions_20261031/RESULTS.md`.
The new statements below extend that mechanism to every packing budget and
replace a brute-force list of repeatable alphabets by a literal hypergraph
matching calculation. The related intercyclic classification and OCF inputs
are inherited, not re-audited. No external novelty claim is made.

## 1. Exact free-alphabet theorem

Fix a tournament C on k labelled core vertices. For a word w over
{0,1}^k, let T_C(w) have the core C and transitive spine v1→...→vm;
bit c of letter wi is 1 precisely when c→vi. Let A be an alphabet.
Build the following hypergraph R_C(A) on the **core vertices only**:

- {c} is an edge when both values of coordinate c occur in A.
- {c,d} is an edge when some letter in A forms a directed triangle
  with c,d.
- {c,d,e} is an edge when those three core vertices form a triangle.

Then

    max_{w in A*} ν(T_C(w)) = matching_number(R_C(A)),

where ν counts vertex-disjoint directed triangles. The same ν is the maximum
number of vertex-disjoint directed cycles, or directed odd cycles: every
cyclic subtournament contains a triangle.

**Proof, upper bound.** Every triangle meets C because the spine is transitive.
Its intersection with C is one of the listed hyperedges. Vertex-disjoint
triangles give disjoint core supports.

**Proof, attainment.** Take a maximum hypergraph matching. For a singleton
{c}, choose two letters with bits 1,0 at c and place the corresponding spine
vertices in that order. They form c→vi→vj→c. For a pair, choose its witness
letter and one fresh spine vertex. A triple already is a core triangle.
Concatenate the separate spine witnesses in any order. They use disjoint core
supports and distinct spine vertices, so all matching edges are realized
simultaneously. Cross edges cannot destroy those triangles. This constructs
an attaining word with at most twice the matching number letters. QED.

**Consequences.** A* respects a packing budget r if and only if its finite
resource hypergraph has matching number at most r. In particular at most r
coordinates vary, and |A|≤2^r. If A violates the budget, a word of length at
most 2r+2 already witnesses the failure. This tests an arbitrary freely
repeatable alphabet, not just individual letters or adjacent pairs.

### The dual budget is exact as well

Let τ(T) be the minimum number of vertices whose deletion makes T transitive,
and let τ(R) be the resource hypergraph's minimum vertex-cover size. Then

    max_{w in A*} τ(T_C(w)) = τ(R_C(A)).

A resource vertex cover X deletes every possible triangle from every word,
so τ(T_C(w))≤|X|. Conversely put t=τ(R). For every resource edge supply t
spine-disjoint copies of its witness from the packing proof, and concatenate
all of them into one word. A deletion of fewer than t total vertices has a
core part X of size<t. Some resource edge avoids X, and its t witnesses
cannot all be hit by fewer than t spine deletions. (A core triple witness
survives automatically.) Thus this word needs at least t deletions. The
statement is trivial for t=0. This proves exactness in both directions.

The matching and covering numbers need not agree, even in these tournament
resource hypergraphs. The alphabet {8,40} on the ordinal sum of two triangles
has matching number2 and vertex-cover number3. Indeed its seven-vertex word
[40] already needs three deletions, while every word over the alphabet has
packing at most2. Hence a packing-two grammar cannot in general be replaced
by a two-apex grammar without changing its objects. This does not say that
its Hamilton counts lack other two-apex realizers.

## 2. Every bounded-packing tournament has finitely many free phases

Fix r≥1. A tournament with ν(T)=r has a maximum set of r disjoint triangles.
Use their 3r vertices as core C. The complement is transitive, for another
triangle there would enlarge the packing. Thus it is T_C(w).

For each of the finitely many C, let L_C,r consist of words whose tournament
has ν≤r. This language is closed under deleting letters. Moreover it has a
finite forbidden-subsequence basis with every word of length≤2r+2.

Indeed, a forbidden tournament contains r+1 disjoint triangles, each meeting
C and each using at most two spine vertices. Keeping only their spine
vertices gives a forbidden subword of length≤2r+2. Minimal forbidden words
therefore have that bound. The core C itself has exactly packing r, so the
empty word is not forbidden.

A deterministic recognizer stores which proper prefixes of the finite
forbidden list have occurred as subsequences. Its bits only change 0→1.
If M is the number of distinct nonempty proper prefixes, there are at most M
strict state changes in an accepted word. Between them the self-loop letters
form a freely repeatable alphabet A. By section1 each such alphabet has at
most 2^r letters. Consequently every legal word has at most M+1 free phases
and M exceptional letters, with

    M ≤ sum_{j=1}^{2r+1} 2^(3rj).

This bound is deliberately crude. It is effective and independent of word
length; the r=1 inherited implementation has a much smaller M=83.

Writing L_m for the number of legal words of length m gives

    L_m ≤ 2^(rm) sum_{j=0}^{min(M,m)} binom(m,j) 2^(2rj)
        = O_r(m^M 2^(rm)).

To count, specify the at most M state-change positions. There are at most
2^(3r) choices at each of those, and at most 2^r at every other position.
The changing states need not be supplied separately: the recognizer is
deterministic. This is a presentation count, not an isomorphism count.

The exponential rate is sharp. Let C be the ordinal sum of r triangles.
In each triangle choose one variable core vertex and give both other
vertices bit0. The resulting alphabet has 2^r letters. Deleting the r
variable vertices makes every word transitive, so every word respects the
budget. Thus a universal smaller exponential base is impossible.

## 3. A global cut replaces the tempting tensor law

Let C=C3→...→C3 (r copies). A tempting tensor rule would combine any one
of the nine legal binary alphabets from each triangle independently. It is
false: cross-block mixed triangles consume different core resources.

The correct statement for **maximum-cardinality alphabets** is exact:

> There are (2r+1)3^r legal alphabets of size 2^r. Each varies in one
> coordinate of each triangle block; after deleting those variable vertices,
> all fixed bits, in the remaining core's transitive order, are a string of
> 1s followed by 0s.

**Proof.** If two coordinates vary in one block, their singleton resource
edges and the r−1 other core triangles form a matching of size r+1. Thus at
most one varies in each block. Size2^r forces exactly one per block, and
forces the full Cartesian product over those coordinates. Its r singleton
edges already saturate the budget. Any triangle on two remaining core
vertices and one spine vertex would be disjoint from them and violate the
budget. The remaining core is transitive. If c precedes d there, a mixed
triangle is precisely the forbidden fixed-bit pattern 0,1. Hence the fixed
bits must be a prefix of1s then0s. Conversely that pattern makes the whole
core-plus-spine complement of the r variable vertices transitive, so every
triangle hits one of them and the packing is≤r. There are3^r choices of the
variable vertices and2r+1 cut positions among the remaining2r core vertices.
They give distinct alphabets. QED.

For r=1 this recovers the nine binary cube edges

    {0,1}, {0,2}, {0,4}, {1,5}, {2,3},
    {3,7}, {4,6}, {5,7}, {6,7}.

For r=2, only45 of their81 Cartesian products are legal. For example
{0,1}×{1,5}, encoded as {8,9,40,41}, has resource matching number3.
The theorem explains why: the first block has a fixed0 before a fixed1
in the second block. This is a genuine inter-block coupling, not a rounding
or multiplicity error.

There is a further boundary obstruction. A legal smaller alphabet need not
extend to a legal alphabet of maximum size. The complete r=2 resource test
finds51 maximal alphabets:45 squares and the six edges

    {8,40}, {16,24}, {32,48}, {57,61}, {58,59}, {60,62}.

This last classification is **FINITE-EXACT**, for the stated six-vertex
ordinal core only. Completeness is justified: any legal alphabet has at
most one varying coordinate per triangle block, and its projection to each
block lies in a legal singleton or one of the nine edges. It therefore is
contained in one of81 products, all of whose subsets are tested. Thus a
smaller free alphabet can be trapped even though larger legal alphabets
exist elsewhere. Its missing coordinates cannot be added independently.

## 4. Clone revisits cost a packing budget

In a rooted tournament substitute a transitive m-vertex module for the
root. Write its inherited positive expansion as

    H_m = sum_s a_s s! S(m,s),  a_s ≥ 0.

Here a_s counts collapsed Hamilton skeletons with s nonempty module blocks.
If a_s>0, every m≥s realization has at least s−1 disjoint directed cycles.
Take the s−1 external path segments between successive module blocks. Each
starts at an out-neighbour of the root and ends at an in-neighbour, so a
fresh clone closes it into a directed cycle. These external segments are
disjoint, and there are enough distinct clones. Each cycle contains a
triangle. Thus bounded packing≤r for the whole clone ray forces a_s=0 for
s≥r+2. In particular its exponential bases are at most r+1.

This is a one-way implication. A coefficient vanishing statement alone does
not bound all cycle packings. It connects the literal packing budget to the
spectral/clone interface without identifying their full states.

## 5. Transfer contract and consequence for Chan–Pak

- **Source:** Session124's packing-one forbidden-subsequence language,
  and the random common-block/successor-loss source
  `session141_metric_aliquot_square_paths_20261102/independent/refinements/sources/CONNECTIONS.md`.
- **Target:** all-r cycle-packing tournament words.
- **Map:** each triangle becomes its occupied core support, a hyperedge.
- **Preserved:** maximum packing over a freely repeatable alphabet, exactly.
- **Lost:** the order and multiplicity of cycles in an individual word, and
  therefore its Hamilton count. Resource matching is not the OCF polynomial.
- **Sidecar:** the literal word and the positive Hamilton-count automaton
  from `../wildcard/` retain those multiplicities.
- **Decisive test:** build disjoint witnesses for a maximum resource matching;
  compare direct triangle packing for arbitrary words with the upper bound.

Intersecting the finite forbidden-word recognizer with the fixed-core
Hamilton-count automaton yields an exact bounded-packing synthesis grammar.
Every finite tournament occurs for some r, and OCF gives H(T)≥3^r because
its r disjoint triangles can be independently omitted or given one of two
colours. Thus a target N needs only r≤floor(log_3 N). This is an effective
structural reduction, not uniform arithmetic coverage. The remaining goal
is a target-independent productive return or a proof that a suitable bounded
packing grammar covers every allowed odd N.
