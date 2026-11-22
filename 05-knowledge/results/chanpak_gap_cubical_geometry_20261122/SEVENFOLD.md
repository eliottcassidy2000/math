# Why the cyclic sandwich cannot be a sevenfold amplifier

**Proved.** For every nonempty tournament \(T\), its ordered two-slot path-cover count \(Q_2(T)\) is not \(6H(T)\). Hence the cyclic sandwich \(C_3[1,1,T]\) never maps a count \(p\) to \(7p\), even though actual strong \(7p\) realizations exist. This rules out one physical mechanism, not arbitrary sevenfold constructions.

The sandwich count is \(H(T)+Q_2(T)\): there are three quotient words using the \(T\)-port once and a unique quotient word using it twice, with the number of ordered nonempty two-run fillings equal to \(Q_2(T)-2H(T)\).

Cutting each Hamilton path at its \(n+1\) positions yields distinct ordered two-slot covers on an \(n\)-vertex tournament. The extra covers are precisely those whose two nonempty occupied slots cannot be concatenated due to a backward joining edge. If \(T\) is strong and \(n\ge4\), choose a vertex \(v\) such that \(T-v\) is strong, by [Moon's vertex-pancyclicity theorem](https://doi.org/10.4153/CMB-1966-038-7) applied at length \(n-1\). For each in-neighbor \(u\) of \(v\), a Hamilton path of \(T-v\) starting at \(u\) gives a failed concatenation with the singleton \(v\) placed first. For each out-neighbor, use a path of \(T-v\) ending there and put \(v\) last. (Cut the Hamilton cycle in \(T-v\) at the appropriate position.) These give at least \(n-1\) distinct failed covers per eligible vertex, and the two-slot sizes distinguish different singleton choices since \(n\ge4\). There are at least two eligible vertices: if just \(v\) were omittable, no cycle of length \(n-1\) would contain \(v\), again contrary to Moon. Thus
\[
 Q_2(T)\ge(n+1)H(T)+2(n-1)\qquad(T\text{ strong},\ n\ge4).
\]
If \(T\) has two or more strong components, both \(Q_2\) and \(H\) multiply under ordinal sum (split each occupied slot at the component cuts). Hence
\[
 Q_2(T)/H(T)\ge \prod_i(n_i+1)>1+\sum_i n_i=n+1.
\]
In particular \(Q_2/H>6\) at every order at least 5.

For orders at most 4, transitive tournaments have \((H,Q_2)=(1,2^n)\), the cyclic triangle has \((3,12)\), and its ordinal join with a singleton has \((3,24)\). A strong four-vertex tournament necessarily has outdegrees \(1,1,2,2\), hence two cyclic triangles. Its count is \(H=5\) by odd-cycle packing. Counting two-slot covers by the size of the first slot gives \(Q_2=2H+2(4+2\cdot2)+6=32\): every 3-subset has count 1 plus twice its triangle indicator. None of the resulting ratios is six. Thus the simplest universal ratio lift is obstructed globally, not only by testing a finite bank of prime inputs.

# The two small holes are different cut geometries

Here is a short exact structural explanation alongside the wrapper no-go. A least missing admissible count being prime or seven times a prime is an inherited arithmetic reduction; it does not make the two initial holes mere copies of a mod-7 obstruction.

In a strong \(n\)-vertex tournament there are at least \(n-2\) cyclic triangles. Induct via a strong deletion, using Moon; the deleted vertex itself lies on a triangle since the full tournament is strong. So a strong tournament with \(H=7\) would have order at most 5 by positivity of the odd-cycle collection. At order 5 it also has a directed five-cycle, raising the count to at least \(1+2(3+1)=9\). Strong orders three and four give counts 3 and 5. Ordinal multiplication cannot make a prime 7 from smaller nonunits. This proof needs no small enumeration.

For \(H=21\), ordinal multiplication could make a composite 21 only as \(3\cdot7\), already impossible. By [Busch's minimum theorem](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3) a strong tournament of order at least seven has \(H\ge5^{(n-1)/3}\ge25\). Order at most five is handled below, so only strong order six remains.

Delete a vertex leaving a strong five-vertex core. Fix a directed Hamilton cycle on the core and mark each of its five cyclic diagonals by a bit: bit \(i\) is one when \(i\to i+2\bmod5\). Rotation, and reflection followed by taking the converse, identify the binary bracelets of length five. There are just eight. This uses a Hamilton cycle only as a coordinate frame, not as a uniqueness assumption.

For each edge pair \(i,j\) define its symmetric normalized cut weight
\[
 w_{ij}=H(T-\{i,j\})+
 \#\{i\to u\to v\to j,\ j\to u\to v\to i: u,v\text{ distinct outside }\{i,j\}\}.
\]
Indeed all nonempty even vertex-paths between the pair have size two or four; the two-vertex path leaves a three-vertex residual tournament and each four-vertex path leaves a singleton. By the [Grinberg–Stanley odd-cycle formula](https://arxiv.org/abs/2307.05569), an apex beating \(S\) adds twice the weighted count of even-order directed paths from \(S\) to its complement, with each path weighted by the complementary Hamilton count. To rewrite the increment as the directed cut \(\sum_{i\to j,\ i\in S,j\notin S}w_{ij}\) in half-height units, move each reversed-endpoint path contribution to the existing arc direction. Such a path closes along that arc to an even cycle; the correction is a sum of directed cycle circulations, whose net flow across any cut is zero. For a singleton this is its weighted outdegree \(a_i\); for two vertices \(i,j\) it is \(a_i+a_j-w_{ij}\). For sets of sizes three and four use their complementary two- and one-sets in the converse, i.e. use weighted indegrees. Thus only vertex and pair data of the five-core are needed.

The complete bracelet test is below. The last two columns list the union of the small weighted cuts in the core and converse (retain values at most the target plus one); the needed value is never present. Empty/full contacts would make the six-vertex tournament nonstrong, or keep its count at \(H(T)\), so are irrelevant.

| Diagonal word | Triangles | 5-cycles | \(H(T)\) | Needed half-jump | Singleton cuts | Pair cuts |
|---|---:|---:|---:|---:|---|---|
|00000|5|2|15|3|none|none|
|00001|4|2|13|4|5|none|
|00011|3|1|9|6|3, 4, 5|3, 5, 7|
|00101|4|3|15|3|none|none|
|00111|3|1|9|6|3, 4, 7|4, 5, 7|
|01011|4|2|13|4|5|none|
|01111|4|1|11|5|4, 6|6|
|11111|5|2|15|3|none|none|

The supplied bracelet replay computes just these small weights and counts. This is an exact compressed proof using the five-vertex cut geometry, rather than a census of all six-vertex tournaments. It explains the contrast: the hole 7 already follows from the number of cycle lengths and triangles in a single core, whereas 21 falls between allowed **weighted attachment jumps** to existing cores.

## Three explicit rigid sevenfold witnesses

Conversely, sevenfold arithmetic is not forced to have a sevenfold geometric symmetry. With bit \(j\) in row \(i\) encoding \(i\to j\), the following are strong:

| Count | Outgoing rows | Sorted outdegrees |
|---:|---|---|
| \(7\cdot5=35\) | \(24,117,121,114,96,65,1\) | \(1,2,2,2,4,5,5\) |
| \(7\cdot11=77\) | \(52,93,120,113,96,66,1\) | \(1,2,2,3,4,4,5\) |
| \(7\cdot13=91\) | \(116,237,248,241,226,192,128,1\) | \(1,1,2,4,4,5,5,6\) |

They have no nontrivial automorphisms. The automorphism group acts freely on Hamilton paths, so its prime divisors divide the count. A permutation of prime order \(p\) moves a \(p\)-cycle within one outdegree class, but the listed degree multiplicities are smaller than every prime dividing the respective count. This disproves an explanation of the sevenfold lane based solely on forced orbit symmetry, while the preceding sandwich no-go disproves a simple product-like source for the number 7. Neither falsifier rules out more flexible marked-port couplings.

Run `python3 holes-five-bracelet.py` for the bracelet weight tables and rigid witnesses.
