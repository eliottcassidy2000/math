# What paths, cuts and line colors count

The three fibers here consist of full Hamilton paths, different cuts of directed Hamilton cycles, and fixed-path staircase boards. A path whose endpoints close is one of the cuts; the other paths contribute only to the first fiber. An ordinary board is an automorphism orbit of paths. Merging identifies a class with its converse. A blue board is fixed by staircase reflection; the remaining boards label black lines. Counts of boards below refer to line endpoints in the specified class.

## An exact odd and doubled distribution

**PROVED.** Among *all* tournaments whose Hamilton words have one permutation sign relative to a reference, the possible path counts are precisely
\[
 h=3^a15^b,\qquad a,b\ge0.
\]
This follows from the [viewpoint group theorem](VIEWPOINTS.md). Their words of SCCs have \(a\) directed triangles, \(b\) copies of \(Q_5\), and any number of singletons. Every such word occurs. Its automorphism order is \(3^{a+b}\), so it has \(5^b\) boards. The converse operation reverses the word. The [ordinary and palindromic word series](SECOND-PROOF.md#the-ordinary-isomorphism-classes-at-every-size) give the exact number of odd and doubled merged nodes by \((n,a,b)\).

The colors inside each such node are determined as follows. A palindrome has \(5^{\lfloor b/2\rfloor}\) blue boards among its \(5^b\) boards. Indeed the path order fixes the component intervals. Choose any presentation on each left-side block, use its forced reflection on the paired right block, and choose a fixed presentation on a central block, if present. Singletons and triangles have just one board, and \(Q_5\) has five boards with exactly one fixed by reflection. A nonpalindromic word and its reverse give \(2\cdot5^b\) boards, all black. Consequently for prescribed admissible \(h\) in this sector the two multiplicities are explicit, since \(b=v_5(h)\), \(a=v_3(h)-b\). The empty tournament can be removed at size zero.

For example at size 36 and \(h=2025=3^2 15^2\) this sector has words with 20 singletons, two triangles and two five-blocks. The total number of words is \(\binom{24}2\binom{22}2=63{,}756\); exactly \(12\cdot11=132\) are palindromes. Thus 132 merged nodes have 25 boards each, five of them blue, while 31,812 merged nodes have 50 boards each, all black. Other tournaments of height 2025 may have both path signs.

At order \(n\ge6\) every tournament in this sector is reducible. Its fixed-path board begins in the first SCC and ends in the last, with the first vertex beating the last. Complementing the nonpath tiles reverses that outermost arc. The opposite board therefore closes to a directed Hamilton cycle and has a strong tournament with viewpoint group \(S_n\). This holds for every board in the sector, of either color.

## A rigid example separating the three counts at all sizes

**PROVED.** Take the ascending transitive tournament on \(1,\ldots,n\), reverse only \(1\to n\), and use the ascending path as reference, for \(n\ge4\). The cyclic rim has a unique Hamilton cycle, up to rotations: a cycle must use \(n\to1\), and the rest must be increasing. Its \(n\) cuts generate \(C_n\); all its Hamilton viewpoints generate \(S_n\). The tournament is rigid: each outdegree occurs at most twice, so its automorphism group embeds in a product of order-two groups and also has odd order. Reversing the integer order is its unique anti-automorphism.

Write \(E-O\) for the number of even words minus the number of odd words in the chosen fiber.

| Fiber | Cardinality | \(E-O\) |
|---|---:|---:|
| All full paths (also all boards, by rigidity) | \(1+2^{n-2}\) | \(1+(-1)^{n-1}2^{\lceil(n-2)/2\rceil}\) |
| Cuts of directed Hamilton cycles | \(n\) | \(n\) for odd \(n\), zero for even \(n\) |
| Blue boards | One for odd \(n\), \(1+2^{(n-2)/2}\) for even \(n\) | One |

Besides the increasing word, every path has the form \(S^\uparrow,n,1,\overline S^\uparrow\), for \(S\subseteq\{2,\ldots,n-1\}\). Its sign is \((-1)^{n-1+I(S)}\), where \(I(S)\) counts pairs \(b<a\) with \(a\in S,b\notin S\). Pair consecutive interior labels. Switching membership in the first mixed pair cancels opposite signs, leaving \(2^{\lceil(n-2)/2\rceil}\) choices, each with even \(I\). The sign of a one-place rim rotation is \((-1)^{n-1}\), giving the second row. For the last row reflection acts by \(S\mapsto (\overline S)^*\), where the star reverses the interior labels. Its fixed subsets choose one from each mirrored pair; an interior midpoint rules them out. When there are mirrored pairs, switch the two central labels, which are consecutive and have opposite memberships. This pairs their signs. The increasing word remains. [Reproduction](viewpoints.py)

Adjoining one new source leaves the height unchanged, and now the class and its converse are distinct because their strong component words have opposite order. The ordinary board count remains \(1+2^{n-2}\); the merged mass doubles and every board is black. Thus the same exact height and the same nontrivial strong component can give either type of merged node.

## A consequence for the remaining Chan–Pak atoms

**PROVED reduction.** The multiplicative reduction leaves ordinary odd primes \(p\ne7\) and twisted primes \(7p\) with \(p\ge5, p\ne7\), the three other composite atoms already having witnesses. Each realization of one of these atoms has exactly one nonsingleton SCC. For primes \(p\ge5\) and all these twisted primes, its viewpoint group must be the full symmetric group on that SCC, since none of the exceptional strong classes has that height. In particular each such core includes paths of both permutation signs. [Prime and twisted-prime reduction](../chanpak_two_holes_atoms_20261123/README.md)

For \(p\ge5\), any strong core of height \(p\) is rigid. The [Busch floor, THM-4744](../../../01-canon/theorems/THM-4744-interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve.md) bounds its order by \(n\le1+3\log_5 p<p\), and the free path action makes its automorphism order a divisor of \(p\). A permutation on fewer than \(p\) letters cannot have order \(p\). Thus for any such prime, attaining height \(p\) is equivalent to exhibiting a doubled merged node of mass \(2p\) with rigid underlying tournament: append one singleton before a strong core to force converse asymmetry. Conversely rigidity reads its ordinary height directly from half its merged mass. For \(7p\) with \(p\ge11\), the same bound \(1+3\log_5(7p)<p\) makes the strong-core automorphism order either 1 or 7. The two elementary logarithmic inequalities hold at their indicated lower endpoints and their right sides then increase faster. The corresponding ordinary board count must be \(7p\) or \(p\), respectively. A method seeking the second channel must retain sevenfold symmetry alongside the class key.

This gives a necessary test for building the outstanding strong cores. The all-order restrictions on \(H\) and on the ordinary board count \(H/|\operatorname{Aut}|\) both exclude 7 and 21; their general equality as spectra is open. [THM-1370, height exclusions](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md) · [THM-4812, presentation exclusions](../../../01-canon/theorems/THM-4812-presentation-fibres-7-and-21-are-absent-at-every-order-orbit-product-lemma-vertex-transitive-lemma-closed-window.md)
