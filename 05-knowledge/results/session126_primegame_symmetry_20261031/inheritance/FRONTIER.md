# Session126: exact inherited carriers and the n5 label boundary

2026-10-31. **Inherited definitions and proved mechanisms; independently verified finite n5 census; one elementary gauge obstruction.** Source snapshot `8f76b84d2960f9e14ad62a154b8533c77593702f`. No priority claim is made for inherited results. [SOURCE_PINS.json](SOURCE_PINS.json) records original paths, full-file hashes and pinned links; old source files are not duplicated in this portable packet.

## 1. What F=U+S means

For N>=2, let D(N)={d:d|N,1<d<N}. Then

`F=tau(N)-2`, `U=omega(N)-[N prime]`, `S=2^omega(N)-1-[N squarefree]`.

These count all members, prime members, and squarefree members of D(N). The prime subset is contained in the squarefree subset: **F=U+S is a count equality, not a disjoint partition, and S is not the divisor sum sigma.** The exact balanced shapes are p, p³ and p²qr with distinct p,q,r in the last shape. Their (F,U,S) triples are (0,0,0), (2,1,1), (10,3,7). Thus 7 and17 are balanced primes, whereas55=5*11 has (2,2,2) and is not balanced. Sources: [arithmetic definitions and classification](https://github.com/eliottcassidy2000/math/blob/8f76b84d2960f9e14ad62a154b8533c77593702f/05-knowledge/results/trijection_reset_20261020/arithmetic.md); [Session84 divisor frontier](https://github.com/eliottcassidy2000/math/blob/8f76b84d2960f9e14ad62a154b8533c77593702f/05-knowledge/results/session84_guarded_unification_20261028/FUS_DIVISOR_FRONTIER.md).

Closest arithmetic mechanism: the labelled prime-exponent box. It supports exact count updates, but balance is not conserved:4->12->24 has defects F-U-S equal to-1,-1,+1. Prime names, exponent values, and update guards are the useful sidecar for the FRACTRAN comparison—not a numeral matching alone.

## 2. The actual ten-node metagraph

Fix the directed Hamiltonian path5->4->3->2->1. The six nonpath arcs form a64-element bit cube. Ordinary tournament isomorphism classes are merged with their **all-arcs converse**. A literal line pairs tiling t with C(t), obtained by complementing only its six free bits. Grid reflection R permutes coordinates, without bit negation; it realizes full converse followed by vertex reversal. A line is blue iff R(t)=t, and a node is pure-blue iff every tiling in its fibre is reflection-fixed. Equivalently every Hamiltonian path has positional reversal as an anti-automorphism. This is not the one-arc-flip graph. Sources: [corrected metagraph explainer](https://github.com/eliottcassidy2000/math/blob/8f76b84d2960f9e14ad62a154b8533c77593702f/03-artifacts/explainers/09-metagraph.md); [exact n5 atlas](https://github.com/eliottcassidy2000/math/blob/8f76b84d2960f9e14ad62a154b8533c77593702f/05-knowledge/results/merged_metagraph_black_quotient_bridges_20261020/census_n5.json).

At n5 there are12 ordinary classes,8 self-converse classes and10 merged nodes, with32 literal lines:8 blue and24 black, including4 black loops. Retain loops and parallel edges. Complement C does not descend to a single-valued function on nodes.

| Pure-blue tournament | H | Aut | Entire tiling fibre |
|---|---:|---:|---|
| transitive TT5 |1|1|{0}|
| singleton -> C3 -> singleton |3|3|{16}|
| regular cyclic T5 |15|5|{11,46,53}|

Here H is the Hamiltonian-path count, **not F**. The mask gauge is tile order(5,1),(4,1),(3,1),(5,2),(4,2),(5,3), with bit0 denoting the descending free arc. Node IDs are ranks by minimum mask; the pure-blue node IDs are0,8,7 in table order. The regular T5 node has three distinct blue neighbors, so pure-blue does not mean pendant or rigid. T5 is not a Paley tournament. See [pure-blue pendant correction](https://github.com/eliottcassidy2000/math/blob/8f76b84d2960f9e14ad62a154b8533c77593702f/05-knowledge/results/merged_metagraph_black_quotient_bridges_20261020/README.md).

Current inherited strong pure-blue atoms are exactly K1,C3,regular cyclic T5. General pure-blue tournaments are the specified palindromic strong-component stacks, with rigid K1/C3 halves and an optional pure-blue center. The all-order count generating function is (1+z+z³+z⁵)/(1-z²-z⁶). Respect the independent odd-proof repair, not the old unequal-partner-sums test. Sources: [THM-4730 classification and proof repair](https://github.com/eliottcassidy2000/math/blob/8f76b84d2960f9e14ad62a154b8533c77593702f/01-canon/theorems/THM-4730-strong-pure-blue-alphabet-complete-points-of-symmetry-are-k1-c3-t5.md); [THM-2454 current all-order count law](https://github.com/eliottcassidy2000/math/blob/8f76b84d2960f9e14ad62a154b8533c77593702f/01-canon/theorems/THM-2454-pure-blue-alphabet-and-center-law.md).

## 3. What the numbers7,17,55 can and cannot mean

No bounded source search recovered a definition identifying this triple with the three pure-blue nodes. That is a recovery limit, not proof that no such interpretation exists. In the documented gauge, masks7,17,55 belong to mixed nodes1,9,6, respectively.

**Proved coordinate-gauge obstruction.** Any selection of one tiling from each pure-blue fibre contains0 and16 at Hamming distance1. Bit-coordinate permutations, even combined with arbitrary coordinate-wise bit complementation, preserve Hamming distances. The three pairwise distances among7,17,55 are3,2,3. Hence this triple cannot select the three pure-blue representatives in any such cube-coordinate gauge. Arbitrary numerical node labels, nonlinear encodings, and other interpretations are not excluded. The argument was independently reviewed in the five-family lane.

## 4. Ten divisors and ten nodes are different carriers

For p²qr the full exponent box{0,1,2}x{0,1}x{0,1} has12 divisors; removing1 andp²qr leaves10. Divisor complement C(a,b,c)=(2-a,1-b,1-c) and prime exchange Q(a,b,c)=(a,c,b) commute, with fixed counts0,4,2 for C,Q,CQ. The CQ-fixed divisors arepq,pr.

[THM-4733 divisor/tiling comparison](https://github.com/eliottcassidy2000/math/blob/8f76b84d2960f9e14ad62a154b8533c77593702f/01-canon/theorems/THM-4733-divisor-plane-of-p2qr-versus-the-marked-path-tiling-cube-3-4-3-and-4-2-4-reflections.md) gives the explicit fold of the **four mirror-moved bits** of the n5 cube:
`(p1,p2,q,r)->(p1+p2,q,r)`.
It intertwines complement/reflection with C/Q and maps onto the full12-divisor box; restricting away from the two extreme bit patterns gives the ten proper divisors. This is a quotient of a four-bit shadow. It is **not** a bijection with the ten merged tournament nodes and preserves neither H nor tournament isomorphism classes.

The full genuine tiling cube at every n>=3 has a mirror-fixed apex. Therefore reflection*complement has no fixed tiling, whereas CQ fixespq,pr. No equivariant injection of the ten-divisor carrier into that full cube can preserve these two involutions. The reverse assignment also fails because Q has fixed points and bit complement has none. This sharp fixed-point obstruction remains intact after the comparison is made precise.

## 5. Independent verification and portability

Run from any working directory, Python3 standard library only:

```sh
python3 /path/to/check_n5_independent.py > normal.json
python3 -O /path/to/check_n5_independent.py > optimized.json
cmp normal.json optimized.json
```

The checker imports no producer code and reads no data files. It exhaustively canonicalizes all1024 labelled tournaments over120 permutations; reconstructs64 fixed-path tilings; checks all64 inherited node assignments, every blue/black edge multiplicity and each class's H/Aut; and verifies the Hamming-distance obstruction. Compact expected atlas data are embedded and pinned to the source above. Failure gates remain active under-O. Relocation and checksum results are recorded in [VALIDATION.md](VALIDATION.md).

## 6. Concurrent credit boundary

At incoming commit `adfb2645207756de07a639e8632ac8264241026e`, [opus S70 HYP-9448](https://github.com/eliottcassidy2000/math/blob/adfb2645207756de07a639e8632ac8264241026e/05-knowledge/hypotheses/HYP-9448-primegame-phase-law-2-to-the-n-times-7-to-the-largest-proper-divisor-minus-one.md) is explicitly **OPEN/PRE-REGISTERED**. It overlaps comparison of the original tail15/2,1/7,55 with the15/14 variant, ten exponent registers, and proposed phase-entry states2^n7^(D(n)-1), D(n)=n/spf(n). It reports exploratory evidence through1000 and registers predictions through2500. Credit that concurrent investigation; this packet does not import its all-n phase law, assert that the confirmatory run finished, or merge speculative claims. Its supplied C simulator was inspected, not run here.

**Primary-source naming correction.** Conway's [1987 chapter scan](https://gwern.net/doc/cs/computable/1987-conway.pdf), printedp4 Theorem1 and printedp9's A–N table, contains the user's exact tail15/2,1/7,55. The15/14,15/2,55 tail is a different presentation; S70's standard/variant naming is not imported as historical attribution. The chapter was independently re-read here and in the PRIMEGAME lane. Also, S70's checker tests entry into a maximal state-prime-free block (`nostate && !inphase`), whereas its unqualified prose says every state-prime-free configuration. The latter wording fails already at15 after seed2's first step. Only its proposed phase-entry investigation is credited.
