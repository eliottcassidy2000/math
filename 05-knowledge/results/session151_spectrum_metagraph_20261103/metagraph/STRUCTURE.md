# Original and converse-merged line multigraphs

## Definitions and an exact algorithm
For n≥3 use vertices 0,...,n−1 and the path 0→1→...→n−1. Let P={(i,j):j≥i+2}, m=|P|. A bit vector t reverses the indicated non-path edges from the transitive tournament, while all path edges remain forward. Let F(t)=t+1 be the antipodal tiling, and let S(t) permute coordinates (i,j)→(n−1−j,n−1−i). S represents the converse tournament after reversing the marked path; F generally does NOT represent converse and does not descend to a function on tournament classes.

The **original** vertices are tournament isomorphism classes c(t). Its individual colored edges are the 2^(m−1) distinct unordered pairs {t,F(t)}, joined at c(t),c(F(t)). Color blue iff S(t)=t, black otherwise. The **merged** graph retains precisely those same distinct edges, but applies q(C)={C,C^op} to their endpoints. Thus a cross-edge between converse classes can become a loop; parallel edges and loops are not discarded. Degree counts a loop twice. An SC node is a singleton orbit of converse, an NS node a two-class orbit. One may implement this at any n by canonical labeling (or any exact isomorphism oracle), and enumerate one representative of each F-pair. This is an exact procedural classification, not a closed-form classification of the integer adjacency matrices and not a claim to scale efficiently to arbitrarily large n.

Define b(C,D),k(C,D) as the number of individual blue, black lines with endpoint classes C,D, including diagonal loops as edges counted once. S commutes with F, sends C to C^op and therefore yields b(C,D)=b(C^op,D^op), k(C,D)=k(C^op,D^op). If E_U is the set of original classes represented by merged U, then for distinct U,V the merged multiplicity is Σ_{C∈E_U,D∈E_V} b(C,D), and similarly for black; for U=V it is Σ_{unordered C,D∈E_U} b(C,D). This retains formerly cross-class lines as loops without doubling.

## All-order consequences
Let f=floor((n−1)/2); S has f fixed coordinates. Its fixed subspace has dimension d=(m+f)/2=floor((n−1)^2/4). Therefore there are 2^(m−1) lines, of which 2^(d−1) are blue. Both ends of a blue line are SC; all other lines are black. Original C has F_C=H(C)/|Aut(C)| tilings: tournament automorphisms act freely on Hamiltonian paths, with *one tiling per automorphism orbit*, not one orbit for an entire fiber. Rédei's theorem and odd order of tournament automorphism groups imply F_C odd. S acts on an SC fiber as an involution, so its number of fixed tilings is odd. Thus an original SC node has odd blue degree and even black degree, whereas an original NS node has zero blue degree and odd black degree. In the merged graph the NS node has even black degree; every SC node still has odd blue and even black degree. Each color's cross-degree has the same parity as its degree, because loops contribute two. Number of SC nodes is even.

**Stronger black bundle theorem.** For every n≥3 the number of black edges between any given merged pair U,V, including U=V, is even. Proof: S acts on F-orbits and fixes both merged endpoints. If an F-orbit is fixed by S then S(t)=t or S(t)=F(t). The latter is impossible because S has a fixed coordinate, and the former says blue. Thus black edges form free S-pairs, including loops. In particular the black subgraph is a union of parallel doubled edges and doubled loops, not just an even-degree multigraph. The graph obtained by forgetting multiplicity cannot support this claim.

A line {t,F(t)} is an original loop iff c(t)=c(F(t)); it is a merged loop iff c(t)=c(F(t)) or c(t)=c(F(t))^op. This tests loop status exactly without calling the witness an automorphism of T: an isomorphism T(t)→T(F(t)) need not be an automorphism of either fixed tournament. For an original class C, |Aut(C)| counts symmetries of a representative tournament, not automorphisms of the class graph. Graph vertex automorphisms are precisely permutations preserving the two colored integer multiplicities for every unordered pair (including diagonal). If graph automorphisms also permute actual line identities, the kernel over vertex permutations consists of independent permutations of same-color parallel lines, with an additional flip of the two incidence flags of each loop if the multigraph formalism treats flags as distinct. No equality with tournament automorphism groups is asserted.

## Reproducible n=3,4,5 test
Run `python check.py`; `results.txt` records every class and merged node, original-class automorphism orders, color degrees, cross-degree and loops. Canonical labels are minimum adjacency bitmasks under all permutations, with bit 1 indicating reverse direction on lexicographic vertex pairs. They are local labels, not external names. Summary:

| n | original classes | SC classes | merged nodes | all lines | blue | black | original blue loops | original black loops | merged blue loops | merged black loops |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | 2 | 2 | 2 | 1 | 1 | 0 | 0 | 0 | 0 | 0 |
| 4 | 4 | 2 | 3 | 4 | 2 | 2 | 1 | 0 | 1 | 0 |
| 5 | 12 | 8 | 10 | 32 | 8 | 24 | 0 | 4 | 0 | 4 |

At n=5 the regular cyclic tournament is code 76, has automorphism order 5, fiber size 3, and all three incident blue lines are cross-lines. The repository's independently checked source establishes these go to three distinct merged nodes, so it is not a leaf. The checker explicitly verifies its three distinct blue neighbors.

## Source and cautions
- `math/01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md`: §3 already proves the black sigma-twin bundle theorem; it is not novel here. See also THM-796.
- `math/01-canon/theorems/THM-643-blue-black-line-parity-structure*.md`: exact S,F frame and parity, with corrections; n=5 pure-blue non-leaf and n=8,9 counterexamples to a blue-fiber cap.
- THM-643 cites an independently audited n=5 witness; its linked results path is absent from the local snapshot, so this report relies on THM-643 and the independent checker rather than inspecting that underlying audit.
- `math/05-knowledge/hypotheses/HYP-3808*`: historic pure-blue pendant claim is false already n=5; its claimed pure-black-only loop onset at n=6 is an enumerated onset, not proved universal pattern. Original and merged allocation parities differ.
- `math/05-knowledge/hypotheses/HYP-3809*`: sentence claiming all tilings in a fiber are *one* Aut-orbit of Hamiltonian paths is wrong. They correspond to distinct Aut-orbits, their number H/|Aut|.
- n=2 is excluded: S has no fixed coordinate, so the black bundle proof needs separate treatment, and the usual blue line formula exponent is negative.
