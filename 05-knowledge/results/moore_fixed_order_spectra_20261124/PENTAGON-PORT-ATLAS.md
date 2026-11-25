# The pentagon port atlas: 249 sites, and a four-site Hoffman–Singleton control

This sharpens the [290-site colour-aware interval code](PORT-COLOUR-INTERVAL-CODE.md). The degree-57 result is conditional on an [order-11 action](https://www.sciencedirect.com/science/article/pii/S0024379509003735); we use its fixed pentagon and the three [remaining cyclic-orbit strata](ORDER11-ORBIT-PACKETS.md). The code is a marked tournament, not an induced subgraph. It may relabel moving orbits and normalize phase origins.

## 1. Two port colours name every remote vertex

Each root on the fixed pentagon has a colour class of 55 ports (5 orbits of 11). For root colours at pentagon distance two, their ports are joined by an equivariant perfect matching (the unique two-step route to the other fixed root). Ports over adjacent root colours have no direct links: a direct link would close a square through the roots. The port edges therefore live over the distance-two pentagon of colours.

Every remote site has one port neighbour in each colour, giving a five-coordinate fingerprint. For two distinct colours there are exactly 55 missing ordered pairs of symbols: the matching given by the port edge if the colours have distance two, or the matching obtained by walking two steps on the port pentagon if the colours have distance one. Indeed these pairs already have either an edge or a common port neighbour, and no other cross-colour pair has either. Every remaining pair has exactly one common remote neighbour by the Moore condition. Thus **each of the ten two-colour projections is a bijection onto the 55×55 array minus one perfect matching**. This is a near orthogonal array of size 55·54.

Normalize four of the five port matchings to identity by port relabelling along a path; the final matching uses at most \(5!\,11^5\) states. Use colours 0 and 1 as the two coordinates, phase-normalize the first contact, and sort remote orbits by their first port orbit and their second port symbol. These give a complete, canonical list of 5·54 remote orbit names in the chosen port coordinates. For each of the five first-port representatives and each of the three remaining colours, the 54 symbols of that colour occur exactly once on its allowed list: store just a permutation of 54. All port matchings and remote–port contacts now cost at most
$$
 \left\lceil\log_2(5!\,11^5(54!)^{15})\right\rceil=3581\ \text{bits}.
$$

The port winding is a separate observable. Following the five distance-two matchings around the colour pentagon gives a permutation of the 55 phase-labelled ports, defined up to conjugacy under gauge changes. A row can be added using only missing pairs exactly when its first symbol is a fixed point of this winding: these rows are forced by transport along the pentagon and the distance-one holes agree automatically. In particular the missing rows fill to a complete orthogonal array iff the winding is the identity.

## 2. Each remote neighbour row is five simultaneous permutations

For a remote vertex x, fix a port colour. Of the 55 symbols of that colour, exclude x's own contact and the two reached along port edges from its other contacts. The three are distinct: equality would make a triangle or square. Each remaining port has exactly one two-step route from x, through a remote vertex. Consequently **the 52 remote neighbours of x project bijectively onto those 52 allowed symbols in each of the five colours**.

For the code use just the colour-0 projection to index these neighbours. In a permitted symbol slot y there are initially 54 candidate remote orbits/phases (one phase per orbit in that port fibre). Adjacency to x excludes sharing any of its other four port contacts. In each other colour the pair consisting of y and x's contact identifies a unique candidate to exclude unless it is the missing pair. These exclusions are distinct: a remote other than x cannot share two port contacts with x, and x itself has the wrong colour-0 symbol. For the two colours at distance two from colour 0, the potential missing y's are already among the three banned symbols. The other two colours create at most two exceptional slots. Thus **50 options normally, 54 in at most two slots**, before further exclusions. Do not charge 49 uniformly: x's own orbit need not be in the target fibre.

Process selected remote rows in label order. Store one status bit per row and a step from 1 through 5 if cyclic. Its two internal shifts, if any, are then known; do not write them again. In a row remove those slots and any filled by known reciprocal entries from earlier selected rows. At each other slot use an index in the ordered candidate list, excluding self and earlier rows, and fill the reciprocal entry if that row will be read later. This recovers exactly every selected remote row. Encode the varying radices by successively splitting a half-open arithmetic interval, as in the earlier interval code. The product of radices, not the number of slots or the particular list sizes, is the worst-case bound; the decoder computes when to stop from the slot obligations. A dyadic interval index of that bit width suffices.

## 3. Why roughly 114 remote rows suffice, with a uniform bit budget

Use the nontrivial Fourier-mode adjacency B and projection \(P=(7I-B)/15\) from the [Gram decoder](FOURIER-OVERLAPS-AND-THREE-DELETIONS.md). All 25 port rows are already known. The port block of B has row absolute sum at most two, so \(P_{UU}\), where U is the port set, has spectrum in [1/3,3/5]. The Schur residual on remotes,
$$
 P'=P_{RR}-K,\qquad K=P_{RU}P_{UU}^{-1}P_{UR},
$$
is an orthogonal projection of rank \(u=r-25=113,114,115\). This follows from taking the orthogonal complement to the span of the port Gram vectors. A principal basis of the residual plus all ports is a principal basis of P, so the full matrix completes by the earlier exact formula.

Let \(c=5t\) be the cyclic-orbit count. Write \(p'_i=P'_{ii}\) and \(\delta_i=K_{ii}\). There are precisely five unit Fourier port contacts per remote; thus \(1/27\le\delta_i\le1/15\). The nonzero eigenvalues of K are those of \(I-P_{UU}\), from the block identity \(P^2=P\). In particular \(\operatorname{Tr}K=40/3\) and \(\operatorname{Tr}K^2\le(2/3)(40/3)\). Since \(P_{ii}\ge1/3\),
$$
 Q'=\sum_i(p'_i)^2\le Q_0=((270-c)49+268t)/225-8.
$$
Indeed subtract at least \((2/3-1/15)\operatorname{Tr}K=8\) from the original remote diagonal-square sum.

Let S be drawn from the determinantal principal basis of \(P'\), and average also over the ten primitive characters. Put \(I\) equal to the number of off-diagonal shifts with both endpoints in S, \(J\) the number of cyclic rows in S, and \(q=52u-2J-I\) the number of new row digits. The original lifted spectral inequality for the vector \(p'\) on remotes, zero elsewhere, gives
$$
 \sum_{i<j}m_{ij}p'_i p'_j\ge 11u^2/100-4Q_0-268t/225 .
$$
The last subtraction bounds the diagonal cyclic correction by its original value. As in the port-colour note, the mean weighted error \(\sum m_{ij}|P_{ij}|^2\) on remote pairs is at most \(a=(7020-c+2(675-c)-2(135-t))/225\). For K the analogous error is at most \(b=(3/2)\operatorname{Tr}K^2\le40/3\) since multiplicities are at most three. Cauchy–Schwarz in the weighted pair/character mean gives the residual error at most \((\sqrt a+\sqrt{40/3})^2\). The certified integer upper bounds are 93, 92, 90.

The mean of J is between \((36t-c)/15\) and \(36t/15-c/27\). Therefore the displayed spectral estimate minus the error plus twice the lower bound on the mean of J bounds \(\mathbb E(I+2J)\) below by M, say. Choose a basis with
$$
 q+5J/12\ \le\ 52u-M+(5/12)(36t/15-c/27).
$$
This choice exists by averaging and is a basis at every primitive mode by Galois conjugacy. Let D be the ceiling of the right side. Since \(5^{12}\le50^5\), the J cyclic-step indices and q row digits have combined radix product at most \(50^{D-2u}54^{2u}\): there are at most 2u exceptional slot digits, and a step index costs at most 5/12 of an ordinary digit in this bound. Encode the basis subset by its combinatorial rank. Exact integers give:

| c | u | D | Total bits, including 2-bit stratum |
|---:|---:|---:|---:|
| 25 | 113 | 4,760 | 30,847 |
| 100 | 114 | 4,762 | 30,860 |
| 175 | 115 | 4,763 | 30,868 |

The 249-site tournament supplies 30,876 contacts. On 250 sites one can reserve a Hamilton spine and still have all 30,876 nonspine contacts. This is a conditional improvement **290 → 249**, or **291 → 250** with a spine. The decoder does not depend on a table of possible Moore graphs or use the existence conjecture as an oracle.

The eight spare bits suffice for two disjoint directed triangles and a separate decoder-neutral marker: at order 249, both marker readings have at least nine Hamilton paths by the odd-cycle collection formula. The same reservation is possible off the spine at 250 by spacing the three triangle vertices so none of the reserved bits is a consecutive pair. The previous finite Brégman comparison replays with ceilings 2^1397 and 2^1404 respectively. These are uniform windows, not exact counts for a hypothetical graph.

## 4. Small controls and what to try next

**Petersen.** Use an order-3 action fixing a single vertex (on the standard disjoint-pairs representation, cycle three of the five ground points). Its three ports form one phase orbit. Each of the two remote orbits projects to the three port symbols, and its two remote neighbours must use both other symbols. No internal triangle is allowed. There is one candidate in the opposite remote orbit for each symbol, so these two crossed phases force the entire remote block. **Pentagon.** Rotation on the pentagon has one orbit and no fixed root: it is the internal-step case, not a port-atlas instance.

**Hoffman–Singleton.** For the [Robertson–Hafner affine model](https://www.math.auckland.ac.nz/~hafner/his.pdf), take the order-5 action with a fixed pentagon and use the same colour convention as above. There is one port orbit per colour, with zero winding, and four remote orbits. Normalize phases so the five missing pair matchings are all equality. After sorting by the second colour when the first is zero, the three remaining columns must be permutations fixing zero on the completed array. The only such permutations whose difference with the identity is also a permutation are the three slopes 2,3,4; they are mutually compatible. Thus there are **six** candidate ordered assignments of slopes to the three remaining colours. In the physical boundary colour order, the exact row-slot test leaves two possible remote neighbours at every vertex for slopes (0,1,2,3,4), and **zero** for every other slope order. Taking the forced edges in the first case gives a 50-site graph passing every Moore equation entry. Equivalently, the port-only Fourier Gram block has rank four here: its five-cycle has one eigenvalue 2 and four smaller eigenvalues in the nontrivial block. No remote-row packet is needed for this control.

Encode the six assignments by four-site marked-spine masks (2,5,3,4,0,7) in lexicographic assignment order, using tile positions (0,2),(0,3),(1,3). Masks **2 and 5** are a blue complement pair in the same tournament isomorphism class and have **five Hamilton paths each**. Only the first decodes to an extendible port assignment. The rank label is a choice of encoder, not a canonical graph invariant. The structural point is the explicit information loss on forgetting physical placement: neither count nor the unmarked metagraph node decides extension for this atlas.

**Degree-57 search target.** Track the port winding and the simultaneous permutation constraints, not just a single coordinate projection or an unmarked tournament node. Try the identity-winding subcase first as a deliberately stronger test: fill the missing rows and compare with the affine near array of symbols in Z/55, using coordinates \(a+ib\) for i=0,…,4 and b≠0. The differences of these slopes are units mod 55. The replay certifies that every required slot has a candidate passing all five port tests, so mere nonemptiness cannot rule out this synthetic table. But copying the Hoffman–Singleton rule “remote edges preserve slope” is impossible: each slope class has only 55 vertices and a triangle- and square-free induced graph on 55 vertices has average degree at most \(\sqrt{54}\), by counting wedges and using that a neighbour pair is nonadjacent. The required remote degree is 52, so any completion of this affine atlas would need, on average, **more than 44 cross-slope neighbours per remote**. The next meaningful filter is simultaneous five-colour matching with reciprocal rows and Moore equations, not a same-slope lift. The identity winding and this affine atlas are test hypotheses only.

Run the exact companion for the three bit inequalities, the Petersen and Hoffman–Singleton checks, the four-site metagraph calculation, and all ten projections plus fifteen full-size permutation-rank round trips for the synthetic 55-symbol atlas.
