# Reflection-paired black edges: cuts, bridges, and an exact census

Status: general elementary deductions and constructive proofs below; finite-exact
census through n=8; explicitly labelled open extrapolations. No completeness claim
for the metagraph from counts, parity, or these new constraints.

## 1. Conventions and inheritance

The vertices are tournament isomorphism classes merged with converse. A tiling
retains the fixed Hamiltonian path. Complement C flips every free tile; reflection
R permutes tile coordinates without flipping bits. Literal lines are {t,Ct}.
They may be loops and remain distinct when parallel. Blue means Rt=t; black means
Rt!=t. R preserves merged vertices and pairs black lines freely.

The load-bearing inherited results are:

- THM-643-blue-black-line-parity-structure: odd unmerged fibres H/Aut; SC fibres
  have an odd blue count; NS fibres have no blue tilings.
- THM-796-three-sorted-recursive-tiling-line-node-incidence, section 4:
  every black multiplicity, loops included, is even, including within an exact
  reflection-defect fibre.
- THM-649-flip-parity-trilogy-affine-D-binomial-law and THM-790-blue-parity-law-proved: affine cyclic
  triangle current carried by legs and apex.
- THM-2453-palindromic-narayana-law-for-the-rigid-stratum: ordered singleton/C3
  stacks and their Hamiltonian-path/automorphism products.

THM-2477-rigid-strong-classification-through-sixty-two contains a substantially
stronger inherited classification; **none of the new all-order lower bounds here
depends on that classification**.

We use the explorer coordinate order (x,y), y increasing and x decreasing,
with bit 0 following the displayed n -> ... -> 1 path. Class IDs here are ranks
by minimum tiling mask, not the canonical labels in earlier atlas files.

## 2. The paired-black quotient and its exact charge

Write K for the literal black multigraph. Define Q by identifying the two black
lines exchanged by R. Thus each literal multiplicity, including each loop
multiplicity, is divided by two.

If A^K is instead the directed half-edge matrix (loops doubled on its diagonal),
then the number of Q loops at u is A^K_uu/4, not A^K_uu/2.

Let M_u be merged fibre mass and B_u its blue tiling count. Then

    deg_Q(u) = (M_u-B_u)/2.

In particular every NS vertex has **odd** degree in Q, because its mass is twice
an odd unmerged fibre. The black graph K is even-degree, but Q is not generally
even-degree. Define the charge at every vertex by

    chi_u = (M_u-B_u)/2 mod 2,

which is one at NS vertices and can be either zero or one at SC vertices.
For every vertex set S,

    |delta_Q(S)| = sum_(u in S) deg_Q(u) mod 2,
    |delta_K(S)| = 2 sum_(u in S) chi_u mod 4.

Consequences: each Q component contains an even number of charged vertices.
Each shore of a Q bridge contains an odd number. In particular, a component
containing only NS vertices has even order. These follow by the handshake lemma,
not by an unproved tournament classification.

### What is genuine cycle structure and what is forced doubling?

Use the same vertex set, including isolates if any, for K and Q. They have the
same c connected components. Loops count as one edge. Then

    beta(K) = |E(Q)| + beta(Q).

The first term is the canonical family of digon cycles from reflected line
pairs. Over F2, summing the two coordinates in each reflected edge pair maps
the cycle space of K onto that of Q, with kernel spanned by these |E(Q)| digons.
For paired loops, the sum of the two loop vectors is the corresponding kernel
generator. This is a generic graph-doubling identity, useful for interpreting
the census but **not** a new tournament-specific rigidity theorem.

## 3. Constructively forced bridges at every order

Let a_n count words in symbols 1 and 3 with total weight n. Set a_0=1 and
a_j=0 for j<0; then

    a_n = a_(n-1)+a_(n-3).

The number p_n of palindromic such words is

    p_(2r) = a_r,
    p_(2r+1) = a_r+a_(r-1).

### Family I: nonpalindromic rigid stacks

For a word w, replace 1 by a singleton and 3 by C3, with every arc directed
forward between blocks. Strong-component decomposition is unique, so different
words give different tournament classes. If w contains k copies of C3, then

    H(T_w) = |Aut(T_w)| = 3^k.

Every Hamiltonian path traverses the blocks consecutively, so the tiling fibre
has one element. Converse reverses w. Each nonpalindromic pair {w,reverse(w)}
therefore gives an NS merged vertex with M=2, B=0 and deg_Q=1.
There are (a_n-p_n)/2 such Q leaves.

### Family II: mixed leaves with a four-vertex core

Let S4 be the unique strong tournament on four vertices. In the fixed-path
coordinate order its tiling fibre is exactly {1,3,5,6,7}; its blue subset is
{1,6,7}. Thus F(S4)=5, B(S4)=3.

For any singleton/C3 word A of total weight r, form

    T_A = A -> S4 -> A^op.

The outer blocks have singleton tiling fibres. Reflection swaps their two
copies and acts only on the five core tilings. Hence T_A is self-converse with
F=5, B=3, and deg_Q=1. Unique ordered SCC decomposition distinguishes the a_r
choices. This produces a_((n-4)/2) mixed leaves at every even n>=4.

### Why these leaves give different bridges when n>=5

It is not enough to count leaves: two leaves could be opposite ends of the
same edge, as happens at n=4. Cyclic triangle current excludes that collision
for the constructed families when n>=5.

For Family I, let alpha,beta indicate whether the first and last block are C3.
The unique tiling has one reversed distance-two arc in each C3 block. Its
apex bit is zero for n>=4. The affine leg law gives

    c3(T_w)=k,
    c3(Ct)=k+n-2-alpha-beta >= n-2.

Here k>=alpha+beta because the first and last C3 blocks, when both present,
are distinct (the single C3 case is n=3).

For Family II with n>=6, let k be the number of C3 blocks in A, and alpha
indicate whether A's first block is C3. The core lies strictly inside the two
path endpoints. Every black core tiling therefore satisfies

    c3(T_A)=2+2k,
    c3(Ct)=n+2k-2alpha >= n.

At a fixed order n, every constructed leaf has c3 <= (n+2)/3. Since
(n+2)/3 < n-2 for n>=5, its partner cannot be another constructed leaf.
Every such degree-one edge is a bridge, and these bridges are distinct.

**All-order constructive bound, n>=5:**

    #bridges(Q_n) >= L_n
      = (a_n-p_n)/2 + [n even] a_((n-4)/2).

The first values L_5,...,L_12 are 1,3,3,6,7,14,17,30. At n=4 the two
constructed leaves are joined to each other, giving one bridge, not two.

This bound is attained in the independent full census through n=8.
Equality at larger orders is **OPEN**, as is the claim that all Q bridges are
leaf edges. The constructive proof requires neither equality nor connectivity.

## 4. A correction: pure-blue does not mean pendant

At n=5, the regular five-vertex tournament has H=15, Aut=5 and fibre

    {11,46,53}.

All three tilings are reflection-fixed. Their complements {52,17,10} land in
three distinct merged nodes. It is a pure-blue vertex of degree three, not a
pendant vertex, even in the underlying support graph.

All three neighbors have score sequence (1,2,2,2,3) and c3=4, but their
(H,Aut,F) triples are (11,1,11), (15,3,5), (13,1,13), respectively.
Thus even identical score/triangle shadows do not merge the three branches.

The claim “pure-blue is always pendant” in the historical HYP-3808 and
HYP-3809 material is false. The prose “PB ... appear only as blue leaves” in
THM-643's census discussion is also false; this does not affect that theorem's
proved parity statements. The strongest automatic survivor is: **a pure-blue
vertex of mass one is a leaf**. At n=7 a corresponding pure-blue degree-three
example has fibre {1216,9600,10560}.

## 5. Independent exact census

The enumerator uses no inherited canonical-isomorphism classifier. For each
unseen fixed-path tiling, it generates every directed Hamiltonian path of its
tournament and renormalizes it to the fixed path, obtaining the entire
isomorphism fibre. Different discovered fibres must be disjoint. Reflection
then pairs converse fibres. A brute-force permutation implementation verifies
all fibres independently through n=5.

| n | tilings | classes | merged | SC | blue lines | black lines | Q bridges |
|---:|---:|---:|---:|---:|---:|---:|---:|
|3|2|2|2|2|1|0|0|
|4|8|4|3|2|2|2|1|
|5|64|12|10|8|8|24|1|
|6|1,024|56|34|12|32|480|3|
|7|32,768|456|272|88|256|16,128|3|
|8|2,097,152|6,880|3,528|176|2,048|1,046,528|6|

| n | blue loops | black loops | Q loops | beta(Q) | beta(K) |
|---:|---:|---:|---:|---:|---:|
|3|0|0|0|0|0|
|4|1|0|0|0|1|
|5|0|4|2|6|18|
|6|2|24|12|209|449|
|7|0|114|57|7,797|15,861|
|8|4|618|309|519,740|1,043,004|

Every active blue subgraph is connected through n=8. Every nonempty active
black/Q subgraph is connected through n=8. These are **FINITE-EXACT**, not
all-order claims. Counts in the JSON include isolated vertices: for example
n=8 blue has 3,353 components (one active plus 3,352 NS isolates), while
black/Q has four (one active plus three pure-blue isolates).

Every Q bridge through n=8 is a leaf edge. The NS/mixed leaf counts at
n=4,...,8 are (1,1),(1,0),(2,1),(3,0),(5,1), exactly the two constructed
families. `structural_checks.json` contains every bridge endpoint with its
representative mask, type, mass, blue count, c3 and quotient degree.

At odd n=3,5,7, the SC blue-count spectra on even-c3 and odd-c3 vertices are
identical, a stronger finite observation than even multiplicities in the
combined spectrum. The respective common spectra are

    n3: {1:1}
    n5: {1:2,3:2}
    n7: {1:5,3:6,5:8,7:16,9:9}.

No class-level involution proving this odd-order equality is supplied. The inherited
blue triangle-current law proves blue bipartiteness at odd n, but bipartiteness
alone does not imply equality of the degree spectra on its two sides.

### Two further hostile controls at n=8

The full SC blue-count spectrum is

    {1:3,3:2,5:4,7:6,9:6,11:4,13:6,15:8,17:7,19:15,
     21:7,23:15,25:15,27:11,29:17,31:17,33:17,35:8,37:6,39:2}.

Thus the all-orders all-even-multiplicity extrapolation in
HYP-3809-merged-metagraph-bucket-conjecture-atlas is false at n=8: already
blue count 1 has odd multiplicity 3. This does not refute the separately
stated odd-order observation.

The maximum blue count 39 also exceeds the old conjectural cap
3^floor((n-2)/2)=27 in THM-643 C1. Two witnesses have explorer masks 40077
and 40109, with (H,Aut,F,B)=(637,1,637,39) and (621,3,207,39),
respectively. Both have score sequence (3,3,3,3,4,4,4,4) and c3=20.
These are refutations of conjectural extrapolations, not of THM-643's
proved parity theorems. No novelty claim is made for finding these failures.

## 6. Reproduction and scope

    python3 census.py --max-n 7 > census.out
    python3 census.py --min-n 8 --max-n 8 --output n8 > n8/census.out
    python3 family_checks.py > family_checks.out

`census_n*.json` includes every tiling-to-node/class map, full literal
blue/black/Q multiplicities, class H/Aut/triangle/score data, per-node
reflection orbit counts, and graph components/bridges/cycle ranks.
The compressed n=8 JSON contains the same data as the uncompressed file.
`family_checks.py` independently enumerates the constructed fibres and audits
the separation inequalities through n=12; the proof above is all-order.

The exact universe is all 2^binom(n-1,2) fixed-path tilings, with no pruning
or score-based merging. No n=9 full census was attempted. Integer signatures,
cycle ranks and connectedness are not claimed to characterize realizability;
count-preserving rewires remain genuine hostile controls.

The compact public packet omits the 33 MB full n=8 atlas and its compressed
copy. They were generated and retained locally; `FULL_ATLAS_SHA256SUMS` records
their reproducible uncompressed hashes. The generator, compact certificates,
complete small atlases through n=7, figures and audit report are included.
