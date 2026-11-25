# How many branches force the Moore matching rows to mix?

**PROVED conditional on existence.** In a degree-57, diameter-two Moore graph choose any root. Any nine of its 57 branches give a connected matching cover on 504 sites. Every fixed seven branches have at least 43 connected extensions among the 50 ways of adjoining an eighth branch. Every nine-branch cover has at least 50 monochromatic edges under every binary coloring of its sites. The third conclusion forces odd cycles wholly inside these nine branches.

These tests concern the candidate graph before a global completion. The model is the [reciprocal Latin-row chart](../moore57_thirteen_20261125/LATIN.md). A head is a neighbor of the root. Its 56 other neighbors are its branch, also called a fiber. Any two branches are connected by a perfect matching. The composition around a closed branch walk is a permutation of the starting fiber, called its holonomy. Its orbits under all closed walks on a branch set are precisely the sheets of its connected components. A single simple base cycle of length \(k\) lifts a cycle of permutation length \(d\) to a graph cycle of length \(kd\). Thus even a triangle may lift entirely to even cycles. In the Hoffman–Singleton graph every triple at the tested root lifts to three hexagons. [Rooted transport identities](../moore57_switching_flags_20261124/TRANSPORT.md) · [Three-branch cycle decomposition](../moore_mobius_prime_packets_20261124/README.md)

## Decoding an oriented cycle or a packing

Fix the matching permutations \(\pi_{ij}\) and an odd length \(k\ge3\). Choose a closed word of branch indices \(i_0,i_1,\ldots,i_k=i_0\) with successive indices distinct, and a starting slot \(x_0\) in branch \(i_0\). Propagate \(x_{r+1}=\pi_{i_ri_{r+1}}(x_r)\). Closure means precisely \(x_k=x_0\), a fixed point of the word's holonomy. This gives a **simple** cycle exactly when the physical vertices \((i_r,x_r)\), \(0\le r<k\), are pairwise distinct. Repeated branch labels can therefore be essential. For chordlessness, test every cyclically nonconsecutive pair \(r,s\) with different branches and require \(\pi_{i_ri_s}(x_r)\ne x_s\). Vertices in the same branch have no edge. Conversely every simple odd cycle wholly on branch sites gives precisely this record after choosing its starting vertex and traversal direction.

For an undirected cycle identify cyclic rotations and reversal of the **lifted** record. For an oriented traversal identify rotations only. If the host edges carry prescribed directions, retain a traversal only when each matched edge points forward, including the closing edge. Latin rows specify adjacency; these arrows must be supplied separately. A collection of decoded cycles is vertex-disjoint exactly when their slot sets are disjoint separately in every branch. Edges between the cycles are allowed. Equivalently, in the commutative ring with one variable \(u_{i,x}\) per site and \(u_{i,x}^2=0\), form the product of \(1+z_{|C|}\prod_{(i,x)\in C}u_{i,x}\) over the desired canonical cycles. The sum of the coefficients of the surviving squarefree site monomials gives the length-enumerated packing polynomial. This applies to the [13-cycle packing counts](PACKINGS.md) once the permitted vertices and chord convention are chosen.

To include cycles using a branch head or the root, add these physical vertices to the same records and use the full adjacency rule: root to any head, each head to its own slots, and cross-branch slots by the matching test above. Canonicalization, simplicity, orientation, chords and disjointness then use all recorded vertices. A tournament on the same sites supplies additional arcs and needs its own adjacency and direction tests; its Hamilton path count does not follow from the Moore matching rows alone.

## Nine branches, the eight-branch boundary and the missing edges

Write \(B\) for adjacency on the 3,192 sites of all 57 branches. On the subspace of vectors with zero sum separately in each branch,
\[
 B^2+B=56I.
\]
This is the unique direct-or-two-step route identity from the transport model. Given \(m\) branches and a component with \(s\) sites per fiber, let \(w\) there be its indicator minus \(s/56\), extended by zero on all omitted fibers. The equal size follows from the bijective matchings. Each selected matching respects component membership, so \(\langle w,Bw\rangle=(m-1)\|w\|^2\). The largest available eigenvalue is 7. When \(m\ge9\), a proper component would give a nonzero vector with quotient at least 8, which proves connectivity.

At \(m=8\), a proper component attains equality and \(Bw=7w\). Evaluation at any omitted site shows it has precisely \(s/7\) neighbors in the component. Consequently \(s=7a\). Also \(s\ge m-1\) for any \(m\ge2\), since from one site to another selected branch the direct matching and the \(m-2\) routes through distinct selected branches have distinct endpoints in the component. If \(m=8,s=7\), adjoining its eight heads and the root gives an 8-regular graph of girth at least five on \(1+8+56=65\) vertices. It attains the Moore bound. Its two nonconstant eigenvalues would be \((-1\pm\sqrt{29})/2\) with equal multiplicities 32, violating trace zero. This excludes \(a=1\) in every component. For a disconnected eight-branch chart the sheet profile is therefore seven times a partition of eight into at least two parts, each at least two. The same neighbor-count argument applies to any union of components.

There is a design on each component of size \(7a\) per fiber. On its \(56a\) sites join two in different fibers when their unique direct-or-two-step route is missing from the selected eight fibers. Each site lacks \(7a-7\) contacts in each of seven destination fibers. This deficiency graph is \(49(a-1)\)-regular. The \(a\) neighbors in the component of any omitted site form a clique in it. Sites in a fixed omitted branch partition the component into 56 such cliques, since each matching is bijective. The 49 omitted branches partition the deficiency edges by the uniqueness of the midpoint. At \(a=2\) the deficiency graph on 112 sites must split into 49 perfect matchings.

## The integer energy of misaligned rows

Choose arbitrary sets \(A_i\) of \(s\) sites on each of \(m\ge2\) branches. Let \(b=\sum_{i<j}|\pi_{ij}(A_i)\setminus A_j|\), with \(\pi_{ij}\) the matching from branch \(i\) to \(j\). Equal sizes make the cost independent of the order. Exactly \(2b\) matching edges cross between selected and unselected sites inside the chart. At an inside site \(x\), let \(d_x\) count its incident crossing edges. At an outside site \(y\) count its selected neighbors by \(n_y\). The mean in each outside branch is \(\mu=ms/56\).

Use the same centered indicator \(w\), with \(\|w\|^2=ms(56-s)/56\). Across the matchings \(\langle w,Bw\rangle=(m-1)\|w\|^2-2b\), and \((B-7I)^2=15(7I-B)\). On the inside, \((B-7I)w=(m-8)w+e\), where \(e_x=-d_x\) on selected sites and \(e_x=d_x\) on the other sites. Each crossing edge contributes \(-1\) to \(\sum_x w_xe_x\). On the outside, \(((B-7I)w)_y=n_y-\mu\). Expanding the squared norm gives
\[
 \sum_x d_x^2+\sum_y(n_y-\mu)^2
 =(4m-2)b-(m-8)(m+7)\frac{ms(56-s)}{56}.
\]
The first sum is at least \(4b\). In each omitted fiber the second sum is at least \(56\{\mu\}(1-\{\mu\})\), by integrality at fixed mean. Hence
\[
 b\ \ge\ \left\lceil\frac{(m-8)(m+7)ms(56-s)/56
 +(57-m)56\{\mu\}(1-\{\mu\})}{4m-6}\right\rceil .
\]
At eight branches and \(r=s\bmod7\) this is \(\lceil28r(7-r)/13\rceil\). In particular 13 candidate sheets per fiber must lose at least 13 one-direction contacts on eight branches and at least 56 on nine. If the eight-branch bound is attained, the identity reduces to
\[
 \sum_x d_x(d_x-1)+\sum_y(n_y-1)(n_y-2)=2.
\]
Thus there is exactly one exceptional site. It is either internal with two mismatches or external with zero or three selected neighbors. All other internal mismatch counts are 0 or 1 and all other external counts are 1 or 2. This gives a test on a proposed 13-sheet arrangement without assuming a root-fixing order-13 action. Attainment is open.

## Two exceptional extensions of the same seven branches

Fix seven branches and add either branch \(j\) or branch \(k\). Suppose each extension is disconnected. Take a nonempty proper union \(C\) of components of the first, of size \(7a\) per branch, and a union \(D\) of the second, of size \(7b\) per branch. On each common fiber let \(u\) count the sites in the intersection. The matchings on the seven branches carry both sets and preserve \(u\). Evaluate \(Bw_C=7w_C\) at branch \(k\), where its right side vanishes. Sum over the \(7b\) sites in \(D_k\). Each of the seven shared fibers contributes \(u-7ab/8\). Branch \(j\) contributes \(h-7ab/8\), where \(h=|\pi_{jk}(C_j)\cap D_k|\). This gives the exact integer constraint
\[
 h=7(ab-u),\qquad ab-\min(a,b)\le u\le ab.
\]
It holds for all four choices of a side or its complement. In particular two extensions cannot give the same bipartition on the core. Taking the same lighter side of size \(7a\), with \(2\le a\le4\), would demand \(7a\le a^2\).

The orbits on one fiber of the common seven-branch cover have weights at least six by the distinct route bound above. Each extension cut groups these indivisible weights. This gives a small finite proof of the 43-extension bound.

1. Enumerate every integer partition of 56 into parts of size at least six. On its indexed parts enumerate the nonempty cuts whose lighter side has total weight 14, 21 or 28. Take one orientation in a balanced cut.
2. Join two cuts if the displayed intersection bounds hold for all four oriented pairs. Pick any binary component cut for each exceptional extension. They must be distinct and pairwise adjacent, giving a clique in this graph.
3. Exhaust all cliques by choosing their next vertex in increasing order and intersecting the remaining candidate mask with its neighbors. Prune a mask only when adding every remaining vertex cannot beat the incumbent. The counts of weight partitions by clique number are as follows.

| Maximum number allowed by this test | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Weight partitions | 407 | 675 | 218 | 41 | 8 | 2 | 1 | 1 |

**FINITE-EXACT** on all 1,353 partitions. It follows that at most seven of the 50 extensions can be disconnected. Double counting incidences between seven-sets and eight-sets shows at least \(43/50\) of all eight-sets at each root are connected. If seven exceptions occur at a common core, the unique admissible profile on it is eight orbits of size seven each. All maximum cliques are then of one of the following two types.

- In 240 cases, the seven chosen lighter triples form the seven lines of the Fano plane on seven orbits, leaving one orbit unused. Each pair of triples meets once, and their union has size seven. Equivalently the 21 pairs of points are covered exactly once.
- In 30 cases, the seven balanced cuts are the parallel classes of planes in affine 3-space over the two-element field. Orient the cuts away from one fixed orbit. Their seven indicator words of weight four close under pairwise distinct sums, so together with zero they form a 3-dimensional binary code. In a generator matrix every nontrivial linear functional is one on four columns and zero on four. Binary Fourier inversion gives each of the eight possible columns exactly once, which is the affine-plane description.

The Fano case fails when all seven extension constraints are applied to one site. Fix an exceptional branch corresponding to line \(L\). Its 35 sites outside its chosen component cut have every neighbor in the seven common branches on one of the four off-line Fano points or on the spare orbit. For such a site let \(n_1,\dots,n_4\) be its neighbor counts on the four off-line orbits. Each of the other six exceptional lines contains one of the six pairs of off-line points. Applying its equality eigenvector at the site counts exactly three neighbors in its cut, of which either zero or one comes from that exceptional branch itself. Each pair \(n_i+n_j\) therefore lies in \(\{2,3\}\). At most one count can be at least two, and a count of three would force three zeros, violating the pair lower bound. In particular \(\sum n_i\le5\). The total capacity of the 35 sites is at most 175. The four off-line orbits require \(4\cdot7\cdot7=196\) contacts by bijectivity on each common fiber. The other 21 sites have all seven contacts on \(L\), giving the contradiction.

**Only the affine case survives seven exceptions.** Label the eight atoms by \(v\in\mathbb F_2^3\) and label the seven exceptional branches by their distinct nontrivial characters \(\chi_j\). At a site in any ordinary outside branch count its seven core neighbors by atom, giving a vector \((n_v)\) of nonnegative integers with sum seven. For each of the seven balanced cuts the core contribution is three or four, since adding the single contact to the corresponding exceptional branch gives four. Thus every nonconstant unnormalized Fourier coefficient \(\sum_v n_v\chi_j(v)\) is \(\pm1\). Parseval gives \(\sum_v n_v^2=(49+7)/8=7\). Exactly one atom is missing and each other atom occurs once.

At a site of exceptional branch \(j\), all seven core neighbors are in its own four-atom side for character \(\chi_j\). The exceptional Fourier coefficient has absolute value seven and each of the other six has absolute value one. Thus \(\sum_v n_v^2=(49+49+6)/8=13\), the minimum possible for seven observations on four atoms. The multiplicities are \(2,2,2,1\). Label this site by the unique atom of multiplicity one. Label each ordinary site by its missing atom and each core site by its own atom. Every branch has seven sites of each label, by summing its core contacts on each atom. For an exceptional branch \(j\), **every matching from it preserves the bit \(\chi_j\) of the labels**. On edges to the seven core branches this is component membership. At any other destination site the core sum of \(\chi_j\) is minus the character of its label, by the two multiplicity profiles. The matching from exceptional branch \(j\) must cancel this sum by its eigenvector equation. In particular for any dependent triple of exceptional characters \(\chi_a\chi_b=\chi_c\), all three mutual matchings preserve both independent bits. Each seven-point orbit of the core is thus subject to simultaneous constraints across every one of the 50 remaining branches. Existence at this numerical bound is open.

## Nine branches force odd cycles at a quantitative rate

On nine branches color the 504 sites by signs \(w_x=\pm1\) and put zero on every other vertex of the full graph. Write \(A\) for its adjacency, \(e\) for the number of monochromatic matching edges in these nine branches, \(t_i\) for the sum of signs in branch \(i\), and \(T=\sum t_i\). The selected graph is 8-regular. From \(A^2+A=56I+J\) obtain
\[
 \|(A+8I)w\|^2=15\langle w,(A+8I)w\rangle+T^2=60e+T^2.
\]
At an inside site the squared residual is \(4d_x^2\) where \(d_x\) now counts same-sign neighbors. These contributions total at least \(4\sum d_x=8e\). Each of the 2,688 outside sites has nine signed neighbors, giving an odd residual. The nine heads contribute \(\sum_i t_i^2\). Therefore
\[
   52e+T^2-\sum_i t_i^2\ \ge\ 2688 .
\]
Summing the signs at both ends of the selected graph’s edges gives \(8T=2(e_+-e_-)\), with \(e_+,e_-\) its monochromatic counts by sign. Thus \(|T|\le e/4\), and \(\sum_i t_i^2\ge T^2/9\) gives \(52e+e^2/18\ge2688\). This excludes \(e\le49\). Also \(e\) is even since the graph has 2016 edges and every cut in an 8-regular graph has even size. If \(e=50\), \(T\) is even, \(|T|\le12\), and \(52\cdot50+8\cdot10^2/9\) alone would not exclude \(|T|=10\). The branch sums are even integers, so \(\sum_i t_i^2\ge 2|T|\). For \(|T|\le10\) the left side is then at most 2680. Consequently equality at 50 requires \(|T|=12\) and \(\{e_+,e_-\}=\{49,1\}\). In particular every coloring with \(|T|\le10\) has at least 52 monochromatic edges.

This is a lower bound on edge deletions to bipartiteness and forces an odd cycle in every nine-branch subsystem. Connectivity and odd cycles are distinct thresholds in the known degree-7 model. In the explicit Hoffman–Singleton chart at one root all 35 four-branch covers are connected and bipartite, while all 21 five-branch covers are connected and contain an odd cycle. The nine-branch estimate gives neither a prescribed cycle length nor a Hamilton-path count for a tournament.

Replay the cut enumeration, the 13-sheet bounds and the small-graph controls with `python3 -B 05-knowledge/results/chanpak_midpoint61_latin_20261125/latin_cycles.py` from the repository root. [Exact enumeration and controls](latin_cycles.py)
