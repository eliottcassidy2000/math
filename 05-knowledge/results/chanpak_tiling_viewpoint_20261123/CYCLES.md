# Closing the path and pairing the chords of its cycle

**PROVED.** At even order there is a family of marked Hamilton cycles for which every cut gives a literal metagraph line realized by a single arc reversal after relabeling. For every even order \(n\ge6\), all those lines are black. The same family, with its rim edges allowed to vary, generates every even self-converse tournament and every pair of disjoint arc reversals exchanged by a chosen involutory reversing symmetry. In one eight-vertex example the ordinary path counts around the square are 541, 539, 537, 539. The two middle tournaments form one converse pair with exactly 1078 tilings in the merged node.

## From a cyclic viewpoint to a physical flip

A closed path tiling has a directed Hamilton rim \(C=(0,1,\ldots,n-1,0)\). The rooted path obtained by cutting a rim edge is its viewpoint. The path fixes its consecutive edges and leaves one binary choice for every nonpath pair. The literal partner of a tiling complements *all* these positions simultaneously. A line is blue if the path-reversal reflection fixes its tiling, black otherwise. A cyclic carrier here is a tournament with distinguished directed Hamilton cycle up to relabeling, which on a fixed directed rim allows rotations. Identifying converse additionally allows reflections accompanied by reversing all arcs. The underlying tournament can have several distinct cyclic carriers. [THM-4729, open and closed tilings](../../../01-canon/theorems/THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md) · [Earlier cycle stabilizers, reflective carriers and seam counts](../session314_metagraph_constraint_algebra_20261113/CYCLE_CARRIERS.md)

Define \(K(T,C)\) by retaining \(C\) and reversing every chord off \(C\). This involution commutes with rotations and with reflection-converse. At a cut \(e\), the nonpath positions are exactly the off-rim chords and \(e\). Consequently the literal partner of the closed tiling is obtained by reversing *one* arc, \(e\), in \(K(T,C)\). A fixed cyclic carrier under \(K\) thus turns every rim cut into a physical one-arc line of its original merged tournament. This is a sufficient marked construction; other lines can have one-edit witnesses with different rims.

**Fixed-carrier classification.** At odd \(n\ge5\) there are no fixed merged cyclic carriers. At even \(n=2m\ge6\), the number of them equals
\[
 \frac{1}{2n}\sum_{\substack{d\mid n\\d\equiv2\bmod4}}
       \varphi(d)\,2^{(n/d)(m-1)} .
\]
These values at orders 6, 8 and 10 are 6, 256 and 52,432. Every fixed carrier at these even orders has the following description on the labelled rim. The half-turn \(A(i)=i+m\) preserves the rim arrows and reverses the orientation of every chord. No such carrier admits a reflection-converse symmetry of its rim, so none of its cuts is blue.

To prove the classification, encode the chords by \(x_{d,i}=1\) if \(i\to i+d\), with \(2\le d\le\lfloor n/2\rfloor\). At a diameter use the consistency condition \(x_{m,i+m}=1-x_{m,i}\). A rotation giving chord complement has even order when there is at least one chord; this excludes odd \(n\). Reflection-converse at axis \(i\mapsto k-i\) sends \(x_{d,i}\) to \(x_{d,k-i-d}\). It cannot give chord complement if \(2i=k-d\pmod n\) has a solution for at least one admissible \(d\). At odd order use distance two. At even \(n\ge6\) take distance two for even \(k\) and distance three for odd \(k\), including the diameter at six. Thus only rotations can witness a fixed carrier.

Let \(g\) be the least chord period. Its stabilizer has odd order \(n/g\) because the automorphism group of a tournament has odd order. A complementary rotation shifts by \(g/2\) modulo \(g\), so \(g\) is even and the actual half-turn \(n/2\) also complements. Conversely the half-turn condition suffices. Any ordinary reflection-converse symmetry would compose with that half-turn to give a forbidden complementing reflection-converse. This proves the claim about blue lines.

For the enumeration, sum over rotations of order \(d\) that complement the rows and put \(s=n/d\). Each ordinary row contributes \(s\) independent bits when \(d\) is even. Half the length of a rotation orbit reaches the antipode, so the diameter row is consistent precisely when \(d/2\) is odd. It then contributes another \(s\) bits. There are \(m-1\) rows altogether. Double-count pairs consisting of a word and a complementing rotation. Each fixed ordinary carrier contributes \(n\) pairs by orbit-stabilizer. Every such ordinary carrier lies in a distinct converse pair by the absence of reflection-converse symmetry. These two divisions give the formula. At order four there is a single ordinary carrier, which also admits reflection-converse. At order three there are no off-rim chords.

## Every even self-converse class from the same fixed chord data

For a labelled half-turn chord assignment, orient the rim forward initially and reverse the rim edges in a set \(S\). Write \(T(S)\) for the resulting tournament. The map
\[
                    S\longmapsto C\setminus A(S)
\]
preserves the merged class and the Hamilton count: it is converse followed by the half-turn relabeling. Exactly \(2^m\) rim orientations are fixed, one reversed edge chosen from each antipodal edge pair. Each has \(A\) as an anti-automorphism. Every assignment with \(A\) as anti-automorphism occurs exactly once in these fixed sets. There are \(2^{m(m-1)}\) off-rim assignments (choose the \(m\) half-row entries at each distance), so altogether this gives \(2^{m^2}\) labelled tournaments reversed by the specified \(A\).

Every self-converse tournament of even order occurs up to isomorphism. Its group of automorphisms and anti-automorphisms has order twice an odd number. An involution exists, necessarily reversing tournament edges. It can fix at most one vertex since the arc between two fixed vertices could not reverse. At even order it therefore pairs all vertices. Relabel it as \(A\). All reversing involutions of the same tournament are conjugate by an automorphism: their order-two subgroups are Sylow subgroups, and if a conjugating map reverses arcs, composing it with the source involution makes it an automorphism. More precisely, if \(\beta\) is one reversing involution in a given ordinary class, exactly
\[
                    2^m m! / |C_{\mathrm{Aut}(T)}(\beta)|
\]
fixed labelled assignments give this class. To count them choose one of the \(|\mathrm{Aut}(T)|/|C_{\mathrm{Aut}(T)}(\beta)|\) allowable involutions, choose one of \(2^mm!\) relabelings carrying it to \(A\), then divide by the relabelings identifying a copy of \(T\).

For any fixed rim orientation and rim edge \(e\), the two physical flips at \(e\) and \(A(e)\) commute. Starting at a fixed point of the set involution, flipping either gives converse tournaments and equal path counts. Flipping both returns to the fixed set. The construction covers *every* nonfixed arc under an involutory reversing symmetry at even order. Its endpoints belong to different vertex pairs, so arrange them consecutively among the representatives of the pairs and put their mates opposite. The given arc and its mate are then rim edges.

A local test decides when the scores alone force two distinct middle classes. Before flipping \(u\to v\), put \(a=d^+(u), b=d^+(v)\). The score multiplicity polynomial changes by \((z-1)(z^b-z^{a-1})\). Its original coefficients were symmetric around \((n-1)/2\). Comparing the increment with its reversal shows that score symmetry survives exactly if \(a=b+1\) or \(a+b=n-1\). Outside those two cases the middle tournaments necessarily give distinct ordinary classes. Their common merged node has twice the odd ordinary tiling mass. If either score condition holds, a full isomorphism test can still distinguish the two classes.

## Exact integer counts

**FINITE-EXACT.** The column of lines lists smaller endpoint heights for literal unit lines coming from the half-turn construction on a *forward directed rim*. The column of squares lists all \(h\) for which a fixed central orientation and opposite physical flips yield \(h,h+2,h+4\) in either direction. The two middle counts coincide. The square's edits are always physical arc reversals; literal-line membership is a separate condition in this second experiment.

| Order | Lower heights of literal unit lines | Starts of equally spaced physical squares |
|---:|---|---|
|4|3|1|
|6|23, 29, 41, 43|1, 29, 41|
|8|117, 197, 245, 257, 369, 409, 415, 417, 497, 507, 519, 521, 531, 545, 641, 649, 651, 657, 659|1, 29, 41, 173, 217, 537|

Here are three eight-vertex squares. Row \(i\) is an outgoing-neighbor bitmask in zero-based labels. The opposite edge to \(\{i,i+1\}\) is \(\{i+4,i+5\}\), taken modulo eight.

| Initial rows | \(i\) | Before, either single flip, both flips |
|---|---:|---|
|248, 125, 241, 212, 96, 8, 160, 50|3|173, 175, 177|
|154, 224, 227, 246, 70, 209, 129, 16|0|217, 219, 221|
|154, 100, 233, 210, 70, 153, 33, 82|3|541, 539, 537|

In all three rows the endpoint tournaments are rigid self-converse, and the two middle tournaments are rigid and mutually converse but nonisomorphic. The merged tiling masses are therefore \((173,350,177)\), \((217,438,221)\), and \((541,1078,537)\). The last two squares require at least eight vertices by the order-seven maximum of 189. [THM-329, maximum counts](../../../01-canon/theorems/THM-329-a038375-new-terms.md)

Run `python3 -B rim_cube.py` from this directory. The [exhaustive search](rim_cube.cpp) assigns every independent chord bit and every fixed central rim orientation. It tests every rim cut and every antipodal pair by subset path counting. The [Python replay](rim_cube.py) checks the three examples independently by directed vertex permutations and full score-compatible relabelings. The fixed carriers and squares supply exact target pairs and triples in a proper family. Coverage of all admissible odd targets remains open.
