# A Golay dodecad on the two growing chains

**PROVED.** The shifted height \(H+2=(2^{k+1}+1)^2\) comes from the zero word, at every \(k\ge1\). For each \(k\ge36\) the same six-vertex quotient supports 24 switchable physical arcs on one common Hamilton board. Restricting their states to the extended binary Golay code gives **4096 distinct ordinary classes**, all chiral, supported on **13 exact path counts**. The count-preserving part of the affine Mathieu contact symmetry is precisely the dodecad stabilizer \(M_{12}\), while the automorphism group of each individual tournament is \(C_3^w\), with \(w\) its contact weight.

See the [unequal-chain inverse and two-slot law](UNEQUAL-CHAINS.md), the [two-apex variant](APEX-CONTACTS.md) where all 4096 heights eventually separate, and the earlier [Golay cycle register](../../../session212/cycle-spectrum/PORTS-AND-GOLAY.md). These are three different count readings of actual marked arcs.

## Physical placement and path count

Use quotient rows \((26,52,57,50,32,1)\), with sites \(2,4\) replaced by modules \(A,B\). Keep the same external arcs as at their sites. Write \(D_A=H(A)\), and \(Q_A\) for the number of ordered two-path covers of \(A\), permitting empty slots. By the two-slot law,
\[
 H=4Q_AQ_B+3Q_AD_B+D_AQ_B-D_AD_B.
\]
Ordinal composition multiplies each of \(D,Q\): a cover restricts independently to each component, in each slot.

In each length-\(k\) ordered chain reserve 12 disjoint, consecutive triples, then \(k-36\) singletons. On each triple mark its outer arc. The two states are a transitive triple, with \((D,Q)=(1,8)\), and a directed triangle, with \((D,Q)=(3,12)\). The forward two-edge path through the triple exists in both states, so expansion of any fixed Hamilton word of the quotient along the two module orders gives a **common** Hamilton word. All 24 marks are actual off-path tiles in that board. The labels and marked positions do not move when a bit toggles.

Assign one Golay dodecad to the first 12 marks and its complement to the second 12. For split weights \(a,b\), put \(w=a+b\), \(x=k-a\), \(y=k-b\). There are \(a\) triangular SCCs in \(A\) and \(b\) in \(B\); the module reads are \(D_A=3^a\), \(Q_A=3^a2^x\) and likewise for \(B\). The exact formula is
\[
 \boxed{H=3^w(4\cdot2^{x+y}+3\cdot2^x+2^y-1).}
\]
At \(w=0\) this is the shifted-square rail. The Golay-forced split weights are even, so the parenthesis is \(2^y\bmod3\) and never zero. Thus \(v_3(H)=w\). In particular **no nonzero Golay word lands two below any square**: then \(H+2\equiv2\bmod3\). Also every height is \(-1\bmod16\), since \(x,y\ge24\) and \(w\) is divisible by four. So there is no literal height gap of two between any of these coded boards.

The formula is an exact inverse over all \(k\ge36\): read \(w=v_3(H)\); invert \(H/3^w\) by the unequal-chain inverse to get the unique ordered \((x,y)\), if it exists. Then \(k=(x+y+w)/2\), \(a=k-x\), \(b=k-y\). Accept exactly if \(k\ge36\) is integral and the split profile below occurs. Having a residue is not a membership test for this family.

## The 13 profiles

| \(w\) | ordered \((a,b)\) and multiplicities |
|---:|---|
|0| \((0,0):1\) |
|8| \((2,6):132,\ (4,4):495,\ (6,2):132\) |
|12| \((0,12):1,\ (4,8):495,\ (6,6):1584,\ (8,4):495,\ (12,0):1\) |
|16| \((6,10):132,\ (8,8):495,\ (10,6):132\) |
|24| \((12,12):1\) |

To get the split table without enumerating 4096 words, note that octads are a Steiner \(S(5,8,24)\): two octads cannot share five coordinates by minimum distance eight, and their five-subsets exhaust all five-subsets by counting. Self-duality forces even intersections with the fixed dodecad; sizes 0 and 8 would give a weight-4 word by adding the dodecad or its complement. The intersection sizes are therefore 2,4,6. The counts of incidences with one and two dodecad coordinates are \(12\cdot253\) and \(\binom{12}{2}\cdot77\). Solving gives \(132,495,132\). Complementation gives the weight-16 row; adding the dodecad gives the weight-12 entries with first coordinate 4 and 8, and the two extreme entries are the chosen dodecad and its complement. Subtraction from 2576 gives the central entry. The replay also checks the table directly using dodecad mask 0x6fdd. The unequal-chain injection and the valuation \(v_3(H)=w\) give exactly one height per displayed profile at fixed \(k\).

## Three inequivalent symmetries

1. **One tournament.** The six-vertex quotient is prime under modular decomposition (its 56 subsets of sizes 2 through 5 give no proper nontrivial module), rigid (an automorphism acts freely on its 23 paths and cannot have order 23 on six vertices), and not self-converse (its score multiset differs from its converse). In an expansion any module meeting several quotient parts projects to a module in the quotient, so meets them all. A vertex omitted from such a module would see some included quotient parts in each direction, since the quotient is strong; this contradicts the definition of a module. Thus the maximal proper modules are exactly the six parts, including singleton parts. An isomorphism between expansions preserves them individually; the ordered SCC list in each expanded module recovers the marked triangle positions. Inside each triangle the only freedom is cyclic rotation. This proves the class count and \(\operatorname{Aut}(T)\cong C_3^w\). Every codeword gives a separate all-black doubled node after converse merging, with mass \(2H/3^w\).
2. **The labelled edit space.** The affine group \(2^{12}:M_{24}\) acts faithfully on the Golay words with the marked-arc Hamming metric (one bit is one physical arc), preserving in particular the graph of octad moves. It acts *between* tournaments, not by vertex automorphisms of one tournament. The physical family has minimum labelled distance eight; it is not claimed that eight is the minimum after arbitrary relabelling.
3. **The count-decorated edit space.** A height-preserving affine action must fix the zero word, the only word of 3-valuation zero, hence has zero translation. The height determines the ordered split weight. In particular the specific dodecad word must map to itself, the unique word of profile \((12,0)\). Conversely a coordinate permutation in \(M_{24}\) preserving this dodecad preserves both split weights. The resulting subgroup is the [dodecad stabilizer \(M_{12}\)](https://web.ma.utexas.edu/users/allcock/research/algor.pdf), not \(M_{24}\). The merged mass also distinguishes the 13 profiles: the octave of \(H/3^w+1\) recovers \(x+y\), hence \(w\) at fixed \(k\).

## Why the Leech comparison is literal, and where it stops

Use the usual even-coordinate layer of the Leech lattice in units of \(1/\sqrt8\): the vector \(2c/\sqrt8\) lies there for each Golay word \(c\), since the residue pattern is \(c\) and the coordinate sum of \(2c\) is divisible by eight. Squared Euclidean distance in this copy is **half** the number of edited arcs. For an octad \(O\), restriction of the code to \(O\) is the even-weight subspace: its annihilator consists of the dual words supported on \(O\), exactly \(0\) and \(O\). Thus differences across an octad realize all \(2^7\) sign patterns of even sign parity on that support. This realizes precisely the \(759\cdot128=97{,}152\) octad-shaped norm-four Leech vectors as signed physical eight-arc differences; each comes from 32 ordered pairs of words. It does not realize all Leech vectors or all lattice isometries. [Leech–Sloane construction](https://neilsloane.com/doc/Me14.pdf) · [Leech coordinate model](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v17i1r18/pdf/).

The sporadic analogy is at the level of **a finite exceptional seed and an operation**, not a generation theorem for groups. The [classification](https://www.ams.org/journals/notices/201806/rnoti-p646.pdf) lists general families alongside the sporadic exceptions; it does not obtain the former from the latter. Here the six-site quotient really does generate an infinite family by module substitution. A sporadic group acts on the 24-coordinate layer, its stabilizer acts on a count-preserving quotient, and the group on an individual tournament is different again. This separation is why a rigid finite seed can organize infinite growth without its physical automorphism group explaining the growth.

Replay: run dodecad-chain-check.py with Python 3.
