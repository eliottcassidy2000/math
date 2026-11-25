# Gewirtz transport exclusions and two short-guard seam charts

## 1. A sharply transitive set cannot preserve projective points

Any two branches of the [edge-normalized degree-57 deck](TWO-COORDINATE-RIGIDITY.md) need a sharply transitive set of 56 label transports: the direct matching and the 55 two-hop routes (including the normalized hub route when appropriate) take each input label to each output exactly once. This tests **transport** on labels, not automorphisms of the proposed Moore graph.

For the Gewirtz candidate, the 56 labels are one \(\mathrm{PSL}(3,4)\) orbit of hyperovals in \(\mathrm{PG}(2,4)\), with adjacency if disjoint. The [Gewirtz graph description of Brouwer and Van Maldeghem](https://homepages.cwi.nl/~aeb/math/srg/rk3/srgw.pdf), §10.20, records the two families of maximum independent sets: hyperovals through a point and hyperovals avoiding a line, each of size 16. Their disjointness graph is the connected point–line incidence graph, so every automorphism either preserves both types or swaps them.

Let \(M\) be the 56-by-21 point-membership matrix and \(N\) the 21-by-21 point–line incidence matrix. Double counting, or the exact replay, gives \(M^\mathsf TM=12I+4J,\ M\mathbf1=6\mathbf1,\ \mathbf1^\mathsf TM=16\mathbf1^\mathsf T\). Thus \(M\) has full column rank. If every transport \(s\) in a sharp set \(S\) preserves type, its hyperoval and point permutation matrices obey \(P_sM=MQ_s\). Summing \(P_s\) gives \(J_{56}\), hence \(\sum_s Q_s=(16/6)J_{21}\), impossible for an integer matrix. In particular the natural \(\mathrm{PSL}(3,4)\) label model is ruled out; the earlier multiple-of-five character test did not rule it out.

## 2. Restrictions that still apply when dualities are allowed

The line-exterior matrix is \(J-MN/2\), since a hyperoval meets a line zero or twice. Let \(d\) be the number of type swaps in a sharper candidate set \(S\). For a fixed source point coclique, let \(a,b\) be the vectors of image point and line coclique multiplicities. Summing its incidence vector over a sharp set gives \(16\mathbf1\). Since \(M\) has full rank,
\[
 6a-3Nb=(16-d)\mathbf1.
\]
Thus \(d\equiv1\pmod3\). If \(d=1\), then every entry of \(Nb\) would be odd, but one line does not meet every point. The complementary case \(d=55\) is also impossible: multiplying the whole set by a type swap exchanges \(d\) and \(56-d\).

The explicit finite replay generates the full 40320-element point-preserving group and its 80640-element extension by duality. It uses 21 homogeneous points over \(\mathbf F_4\), six elementary transvections, the field automorphism, and the duality taking a hyperoval to its six exterior lines (correct the orbit by a diagonal projectivity if needed). The point-preserving group is complete because the faithful action on points is a collineation stabilizing one of three projective hyperoval orbits; \(\mathrm{P\Gamma L}(3,4)\) has order 120960, and the orbit of the chosen hyperoval family under it has at least three members. The generated 40320 reaches that upper bound. The type-swap coset doubles it. Alternatively these orders are in the cited graph description.

Among the point-preserving permutations that fix no hyperoval, the counts by number of fixed projective points are exactly 11520 with zero and 7560 with one. When \(d\) is even, \(Nb\equiv0\pmod2\); the parity of each individual line multiplicity need not vanish. For \(d=4\) it does, as follows.

### The d=4 boundary is impossible

If \(d=4\), the parity support of \(b\) has size at most four and \(Nb\equiv0\pmod2\). Any nonempty set of at most four distinct lines would have a point unique to one line (it has five points and at most three other intersections). Hence all entries of \(b\) are even **in this case**. There are then at least two equal-image pairs of the four correlations per source point, giving at least 42 point coincidences. But each pair gives, by quotient, a point-preserving hyperoval derangement that fixes at most one point, giving at most six total. This excludes \(d=4\), and by symmetry \(d=52\). No stronger even-\(d\) interval is inferred from even point incidence alone.

For each branch pair let \(D_{ij}=1\) for type preservation and \(-1\) for type swap, and \(D_{ii}=0\). This is symmetric and the normalized hub row is positive off the diagonal. The off-diagonal entries of \(D^2+D\) equal \(56-2d\equiv0\pmod3\); diagonal entries equal 56. So over \(\mathbf F_3\), \((D-I)^2=0\) and \(\operatorname{rank}(D-I)\le28\). The integer \(d\) at every sharp set lies in \(\{7,10,\ldots,49\}\); the nonhub row swap counts obey the same filter by testing the hub-to-branch pair. These are necessary conditions, not a general exclusion of the full Gewirtz transport group.

### A finite obstruction even allowing dualities: no orbit-coset batch

The replay rules out a sharply transitive set in the full group **if it is a union of complete left cosets of a subgroup of order seven**. All such subgroups are Sylow conjugate. Take \(C=\langle c\rangle\) where \(c\) has eight 7-cycles on the hyperovals. For a left coset \(gC\), record the ordered tuple of eight sets \(g(O_i)\), where \(O_i\) are the orbits of \(C\). Two cosets are mutually sharp exactly when their eight corresponding images are disjoint, because for any input in \(O_i\) a coset runs through exactly \(g(O_i)\). Identical tuples cannot both belong to a candidate.

Left-translate to include \(C\) itself. Exhaustive replay over the 80640-element group gives **201** distinct compatible tuples other than the identity coset's, with **252** compatible unordered pairs and **zero** compatible triples. Thus a mutually sharp union can use at most three complete cosets (and three is attained); the 56-target set would require eight. In particular the full Gewirtz group has no regular subgroup on the 56 hyperovals: a subgroup of order 56 acting regularly would be a union of eight such cosets. This still allows sharp sets that are not coset unions. The independent linear-span test confirms its limits: the all-ones matrix lies in the span of the full group's permutation matrices modulo 2 and modulo 3, so those linear relaxations alone do not rule out arbitrary sharp sets.

## 3. Two actual short-guard charts, with small twin-prime replays

The [six-site guard census](GUARDS-AND-DECK-DISTANCE.md) gave two guard signatures \((h,Q,R,J)\) and scalar moves \(111\to113\). Use lexicographic backward-arc masks on all 15 unordered vertex pairs and let the listed outcontacts be local to the six-site block:

| Chart | Masks | Apex outcontacts | \(R,J\) |
|---|---|---|---|
| A | 177 / 181 | \(0,1,4\) | 158,162 |
| B | 11524 / 11588 | \(0,1,2\) | 191,129 |

In either chart use binary payload \(B_0=2v\), then the six-site block, one guard, final zero, and common apex. All external arcs follow the established ordinal construction. The lower heights are exactly
\[
 H_A(v)=1662v-58q(v)+427,\qquad
 H_B(v)=1596v-58q(v)+493,
\]
and each mate is two higher; the move is one actual arc. Full endpoint-subset DP and trial division verify:

| Chart | Payload \(B_0\) | Order | Two prime heights |
|---|---:|---:|---:|
| A | 18 | 14 | 15269, 15271 |
| A, prefix 1010 | 82 | 16 | 67931, 67933 |
| B | 74 | 16 | 58907, 58909 |
| B, prefix 1010 | 664 | 19 | 528509, 528511 |

Each is minimum host order (and lower height) for the corresponding chart restriction, verified using the slope lower bound to cap the finite search. No minimality across all tournaments is asserted.

The one guard is needed for this transfer form. With no guard, the case where the final zero is last requires subtracting the paths of the preceding pointed block ending at the common apex (they cannot then step to the zero). This is \(e\), whose change is \(+2\), exactly canceling the change in \(C\) from the six-site move. The case where the final zero immediately precedes the apex is unchanged by the move. The replay checks the collapse in each small fixture.

### Combine the two readings, without mistaking scalar gluing for an arc walk

For the **same** payload their lower-height difference is \(66(v-1)\), so a prescribed absolute pair of chart readings has a unique candidate \(v=1+(H_A-H_B)/66\), with the full digit formula giving the exact equality test. For modular pairs, take the normalized integer readings \((H_A-427)/2,(H_B-493)/2\). Their matrix on \((v,q)\) is
\[
 \begin{pmatrix}831&-29\\798&-29\end{pmatrix},
\]
with determinant \(-957=-3\cdot11\cdot29\) and entry gcd 1. The [uniform digit-state law](TWO-COORDINATE-RIGIDITY.md) therefore gives complete independent modular control away from \(3,11,29\). For any modulus, the image index is \(\gcd(957,M)\) by the two-by-two integer normal form (all residues in the image are reachable). This is a joint local result; the absolute pair still has to pass its unique-candidate test.

There are also scalar three-level overlaps with **different** payloads. For A with \(v=1274\) and B with \(v=1339\) the two physical seams share a height and their union gives \(2104127,2104129,2104131\). With prefix 1010, \(v=5180,5463\) gives \(8574323,8574325,8574327\). These are shared **count values**, not a claim that the two middle tournaments are isomorphic or that there is a length-two physical walk in the metagraph. This distinguishes arithmetic overlap from the unresolved absolute selector.

Run `python3 gewirtz-transport.py` for the projective and C7 replay, optionally with `--spans` for the linear span checks. Run `node dual-guard-seams.js` for the explicit orientations and prime-height replays.
