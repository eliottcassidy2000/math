# From repeated counts to shared structures: fifteen, its neighbours, and the missing thirty

**Status: CITED classifications; PROVED elementary carrier lemmas below;
FINITE-EXACT complete checks on 4, 6 and 8 labels; proposed transfers OPEN.**
2026-10-18, `frontier_petersen_oct18`, cardinality continuation. The classical
Richelot, Fano, and Hamming correspondences are recovered, not claimed as new.
This is a research atlas, not a claim that all occurrences of an integer have
one cause. No Collatz, plane-colouring, or planar-Jacobian conjecture is resolved.

The useful extension of the owner's technique is to follow a numerical match
until it gives a **map, a preserved relation, or a missing kind of object**.
Here the fifteen matchings of six labels give every Lagrangian in a binary
symplectic four-space. At eight labels, 105 matchings miss thirty of 135
Lagrangians. Those thirty are precisely labelled affine three-spaces, or
extended Hamming codes, or Fano configurations after marking one point.
The failure of the first representation identifies the additional construction.

Artifacts: [exact checker](../../04-computation/experiments/cardinality_carriers_20261018m.py),
[receipt](cardinality_carriers_20261018m.json). This supplements the
[recurrent-number census](constants_atlas_20260926_recurrent_numbers.md) and
[small-integer role index](natural_number_roles_20261004.md).

## 1. Inheritance, scope, and working board

**Anchor:** classify the fifteens and build actual bridges. **Niche:** recover
the repo's method across tournament, coding, modular, and crystallographic
threads. **Wildcard:** use the first count mismatch to find the omitted code.

- Closest proved mechanism: the equivariant-map criterion, Lemma N1 in
  [natural bijections](procgen_tbij_20261001_natural_bijections.md), and the
  explicit tetrahedral representations in the
  [K4 atlas](../../07-reflections/tetrahedral-k4-representation-atlas-tournaments-xor-berggren-h1-codex-20260816.md).
- Canonical hostile: fifteen duads and fifteen matchings have different
  fixed-point counts for the same permutation; see §3.
- Corrected near miss: MISTAKE-059, the false inference that a 60-to-20 map
  has constant fibre size three. The historical reflection's header and final
  status still repeated it despite the correction in its body; synchronized here.
- Least-used relevant sidecar: a marked sixth/eighth label. It gives a lawful
  passage to A5/Fano geometry while making the lost full symmetric-group action
  explicit.

The live board has six concepts: **counted object and equivalence; group
action; incidence/pairing; marked source; extension in size; indecomposable
pieces**. Each proposed connection below is checked against this board.

## 2. The owner's five fifteens, precisely typed

| Occurrence | What has cardinality fifteen? | Scope and evidence |
|---|---|---|
| Sporadic sequences | The standard catalogue of 6 Zagier, 6 Almkvist–Zudilin, and 3 Cooper sporadic Apéry-like sequences | **CITED catalogue.** Straub says fifteen *known*, in two specified normalized recurrence families. This is not an exhaustiveness theorem for every integer recurrence. [Straub, introduction](https://arxiv.org/html/2301.12248v1). |
| Convex pentagonal tiles | The fifteen known families of convex pentagons that admit monohedral plane tilings | **CITED classification**, Rao's exhaustive search. These are families of shapes, not fifteen individual tilings or fifteen congruence classes of pentagons. [Rao](https://arxiv.org/abs/1708.00274). |
| Partitions of size four | Partitions of a **labelled four-element set**, Bell number B4 | **PROVED / FINITE-EXACT:** block shapes 4; 3+1; 2+2; 2+1+1; 1+1+1+1 have 1,4,3,6,1 members. Integer partitions of 4 number five. |
| Perfect matchings of K6 | Partitions of six labelled vertices into three unordered pairs | **PROVED:** 6!/(2^3 3!) = 5·3·1 = 15. Modulo arbitrary vertex relabelling there is one matching type. |
| Elliptic torsion | Isomorphism types of the full group E(Q)_tors for elliptic curves over Q | **CITED**, Mazur: C_n for 1≤n≤10 or n=12, and C2×C_(2m) for 1≤m≤4. Eleven plus four. It is not the number of subgroups of one curve, and changing the ground field changes the problem. [Mazur](https://pmihes.centre-mersenne.org/articles/10.1007/BF02684339/); [Barrios's constructive realization](https://arxiv.org/abs/2012.12475). |

Two further fifteens materially improve this list:

1. **Fifteen genus-two Richelot kernels:** maximal isotropic subgroups of
   J[2] over an algebraic closure of characteristic not two. For a squarefree
   sextic they are its three-quadratic splittings. This is the **same matching
   carrier**, not just the same count. Kernel rationality and distinct quotient
   isomorphism classes are additional questions. Florit–Smith work in
   characteristic not 2,3,5 and explicitly distinguish weighted edges from
   kernels: [§§1–2](https://arxiv.org/html/2101.00917v2).
2. **Fifteen Ogg primes:**
   {2,3,5,7,11,13,17,19,23,29,31,41,47,59,71}. The primes for which all
   supersingular j-invariants are defined over F_p coincide with the prime
   divisors of the Monster. The agreement is of the **actual prime set**, a
   much stronger starting datum than its size. Moonshine supplies modular
   structure; see [Duncan–Ono](https://arxiv.org/abs/1411.5354), and the repo's
   [Ogg note](ogg_triangular_chowla_20260927.md) with its
   [audit](ogg_triangular_chowla_20260927_audit.md).

No map has been supplied between Mazur's fifteen types, the fifteen tile
families, and the fifteen sporadic recurrences. Their equality remains a
candidate generator. Useful refinements are available: moduli and degeneration
graphs for tiles; recurrence operators, transformations and monodromy for
sequences; modular curves and forgetful maps for rational torsion. A numbering
of all three lists from 1 to 15 would not preserve any of these predicates.

## 3. An explicit fifteen bridge, including its action boundary

Let Ω have 2g+2 labelled elements. Put

    V_g = {even subsets of Ω}/(S ~ Ω\S),
    [S]+[T] = [S symmetric-difference T],
    <[S],[T]> = |S intersection T| mod 2.

The even subsets form a (2g+1)-dimensional binary vector space. The all-ones
vector is its radical for the dot product; quotienting it gives a
nondegenerate alternating space of dimension 2g. These statements require
the even-subset condition: it makes both the complement quotient and the
pairing compatible.

A matching M of Ω gives

    L_M = span{[P] : P is a pair in M}.

The pair vectors are mutually orthogonal and have exactly one linear
relation, the sum of all pairs. Thus L_M is g-dimensional and Lagrangian.
The number of Lagrangians is

    L_g = product_(i=1..g) (2^i+1).

**Proof of the count.** Double-count (nonzero v, Lagrangian L containing v).
The quotient v-perp/<v> is symplectic of dimension 2g−2. Hence
(2^g−1)L_g = (2^(2g)−1)L_(g−1), with L_0=1.

For g=2 every nonzero class has a unique two-element representative. Two
different pair classes are orthogonal exactly when their pairs are disjoint.
The three nonzero vectors of a Lagrangian are therefore precisely the three
pairs of a matching. This proves an S6-equivariant bijection

    K6 matchings <-> Lagrangians of V2 <-> genus-two Richelot kernels.

The last arrow is the classical algebro-geometric interpretation; the first
is proved here without using a curve. The geometry is the generalized
quadrangle W(3,2), often called the doily: 15 points (duads), 15 lines
(matchings), three points per line, three lines per point.

**A further bridge to operators.** Modulo scalar phases, the two-qubit Pauli
operators are `(X^a Z^b) tensor (X^c Z^d)`, labelled by F2^4. There are
fifteen nonidentity classes, and the commutation sign is the symplectic
pairing. Thus the fifteen doily lines are also the fifteen maximal commuting
triples of nonidentity classes. This is classical
([Planat–Saniga, Pauli graphs](https://arxiv.org/abs/quant-ph/0703154)).
The checker verifies the commutator using actual signed permutation actions
on four basis states, all 1,024 cases. The map preserves commutation, but
quotienting scalar phases loses product signs: retain those signs before
importing a contextuality or signed-current claim. A symplectic basis is
the extra choice, not a choice-free labelling of operators by branch pairs.

**Equal point and line counts do not erase their roles.** The permutation
(12) fixes seven duads but three matchings. The permutation (12)(34)(56)
fixes three duads but seven matchings. Thus no ordinary S6-equivariant
bijection identifies these two fifteen-element sets. The outer automorphism
twists the action; it does not preserve the original point meaning. The
checker verifies all 720 relabellings on all fifteen matching-to-Lagrangian
images. This strengthens the warning already in
[the six-point observer study](information_observers_collatz_20261014a.md#3-the-outer-automorphism-preserves-the-group-not-its-point-action).

**The marked sixth point reconnects the older A5 work.** Mark ∞ in Ω. A
matching has one pair {∞,v} and two pairs {a,b},{c,d}. Send it to the
involution (ab)(cd) of A5 fixing v. Equivalently send it to the unordered
pair of triples {v,a,b}, {v,c,d}. This is an invertible S5-equivariant map.
The marked label is necessary for this particular construction. It joins
the old [icosahedral fifteen](../../07-reflections/the-icosahedral-fifteen-s6.md)
to the Richelot carrier: fifteen involutions correspond to the fifteen
half-turn axes of the icosahedron, after identifying A5 with its rotation group.

**Another genuine split: 15 = 9+6.** Choose a symplectic decomposition
V2=U⊕W into two planes. Nine Lagrangians are products of one line in each
plane (3×3). Six are graphs of symplectic isomorphisms U→W (|Sp(2,2)|=6).
To prove exhaustiveness, a Lagrangian meeting U in a line also meets W in
a line; otherwise both projections are isomorphisms and it is such a graph.
For a product of two elliptic curves this is exactly product versus gluing
kernels, Florit–Smith §2.2. The sequence catalogue also splits 6+9, but the
matching arithmetic alone supplies no identification of those two splits.

## 4. The first failed extension finds the missing object

| g | labels | nonzero vectors | perfect matchings | Lagrangians | missing from matchings |
|---|---:|---:|---:|---:|---:|
| 1 | 4 | 3 | 3 | 3 | 0 |
| 2 | 6 | 15 | 15 | 15 | 0 |
| 3 | 8 | 63 | 105 | 135 | 30 |

The genus-three split 105+30 is classical: Boyer explicitly calls the last
thirty the Fano disposition ([slides, pp.9–13](https://www.normalesup.org/~iboyer/files/agct14Boyer)).
The following elementary binary-code proof both explains the split and
connects it to other repo lanes.

**Lemma (all g).** Taking the preimage of L in the even-subset space gives
a binary self-dual code C of length 2g+2, dimension g+1, containing the
all-ones vector. Conversely, every binary self-dual code of that length
gives such an L.

**Proof.** The preimage has dimension g+1 and all pairwise dot products zero,
so it equals its orthogonal complement. Conversely self-orthogonality makes
every codeword even; therefore the all-ones vector belongs to C-perp=C.
Quotient it out. Matching codes are direct sums of two-coordinate repetition
codes. This is an identification of labelled linear-algebra objects, not a
claim that every such code defines a rational isogeny.

**Lemma (eight labels).** Either C is a matching code, or it is a labelled
extended binary Hamming code, with weight enumerator

    1 + 14 z^4 + z^8.

**Proof.** If C contains a weight-two word supported on {i,j}, orthogonality
makes all codewords equal on i,j. Subtract that word when needed: C splits
as a repetition code on {i,j} and a self-dual code on the other six labels.
Every self-dual six-coordinate code has a weight-two word: a nonzero word
other than all-ones has weight two or four, and complement turns four into
two. Repeating the split gives four pairs.

Otherwise there are no weight-two words, and complementation also forbids
weight six. The fourteen words besides zero and all-ones all have weight
four. Choose a basis 1,f1,f2,f3. Give each label i the address
(f1(i),f2(i),f3(i)) in F2^3. Two equal addresses would make e_i+e_j orthogonal
to every codeword, hence a forbidden weight-two word in C. The eight addresses
are therefore distinct. C is precisely the evaluations of all affine binary
functions on these eight points. Its fourteen weight-four supports are the
affine planes. A different basis changes the addresses by an affine map.

Thus the thirty missing Lagrangians correspond to affine structures on an
eight-element labelled set:

    8!/|AGL(3,2)| = 8!/(8·168) = 30,
    |GL(3,2)| = (8−1)(8−2)(8−4) = 168.

An automorphism of the affine-plane system is affine: after fixing the
origin, its parallelograms recover addition. Each triple of distinct points
lies in exactly one affine plane, so these fourteen blocks form S(3,4,8).
Mark an origin and take the seven planes through it; removing the origin
gives the seven lines of the Fano plane on seven points. The **same thirty**
also count labelled Fano planes on a fixed seven-set, 7!/168=30.

This explains several previously separated numbers in one construction:

    105 pair-generated subgroups + 30 affine/Fano subgroups = 135;
    63 nonzero ambient vectors = 28 pair classes + 35 four-subset classes;
    16 codewords = 1 + 14 plane words + 1;
    puncture one coordinate: 16 disjoint radius-one balls of size 8
                            cover all 128 seven-bit words.

The punctured code has minimum distance three, since the unpunctured
nonzero words have weight at least four. Its balls are disjoint; their
cardinalities sum to 128, proving the perfect packing. The pairing is the
algebraic predicate, the code is a representation, and the packing is a
different invariant of that representation. None of these steps follows
merely by noticing that both sides have sixteen elements.

**Independent exact verification.** The checker builds all isotropic spans
in the quotient, obtaining dimensions counts (1,63,315,135). It separately
generates all matchings and all 40,320 permutations of one affine-plane
design. The thirty missing spans and thirty resulting designs agree exactly.
Every triple-incidence condition, code closure, self-duality, and punctured
packing is checked. Normal and `python -O` runs give byte-identical receipts.

## 5. Where the repo was already doing this

This was a targeted corpus search followed by reading mathematical claims and
their repairs, not an exhaustive re-audit of every linked paper or every line
of the repository. The following twenty routes show how the technique evolved.
Historical reflections are provenance; their sweeping claims do not outrank
later proofs and counterexamples.

| Recovered route | Number cluster | What to carry forward |
|---|---|---|
| [Icosahedral fifteen, June](../../07-reflections/the-icosahedral-fifteen-s6.md) | 5,6,10,15,20,24,60 | Shared-vertex/pair map to A5 involutions is exact; 60/20 does not give uniform fibres. |
| [Angle-defect spine, June](../../07-reflections/the-angle-defect-spine-five-solids-five-lattices-and-the-seventh-prime.md) | 5,6,7,168 | Curvature sign is useful; equal counts of solids and lattice types do not derive LRC or identify their classifications. Scope warning added. |
| [K4 representation atlas, August](../../07-reflections/tetrahedral-k4-representation-atlas-tournaments-xor-berggren-h1-codex-20260816.md) | 3,4,6,16,24,64 | One S4 acts on vertices, edges, matchings and orientations. Each representation has a different kernel or extra coordinate. |
| [Paley–Fano–octonion design](collatz_mod6_20260922_paley_fano_octonion_design.md) | 7,14,21 | The fourteen cyclic triangles of P7 split into two Fano systems. Preserve signs/orientation in the octonion comparison. |
| [Moonshine/source interplay](procgen_sources_20260923_interplay.md) | 6,11,12,23,24,196884 | Actual stabilizer embeddings, outer twists and eta products connect layers. Not every six-point observable survives the twist. |
| [Recurrent-number census](constants_atlas_20260926_recurrent_numbers.md) | 9,437-file census; entropy constants | Successful unification is an explicit entropy function, [THM-4487](../../01-canon/theorems/THM-4487-dip-spectrum-entropy-curve-and-sharpness-of-thin-divergence.md); numerical proximity alone failed elsewhere. |
| [Coincidence atlas](collatz_coincidence_atlas_20260927.md) | 10,12,24,139 | Dodecahedral/icosahedral quotients and Euler pentagonal identities give maps; the separate 139 appearances need a new one. |
| [Ogg/triangular/Chowla note](ogg_triangular_chowla_20260927.md), [audit](ogg_triangular_chowla_20260927_audit.md) | 15,30,42,168 | Same prime set; class numbers become Paley half-row imbalances. Node/conic/elliptic distinctions control triangular-number questions. |
| [Natural-bijection study](procgen_tbij_20261001_natural_bijections.md) | 3,4,12,56,… | Character/Burnside identities can explain equal orbit counts without natural maps. Stabilizers decide lawful maps. |
| [Natural-number role index](natural_number_roles_20261004.md) | 1 through 87211 | Group by what a number does: norm, clock, dimension, orbit size, or obstruction. |
| [Trees versus converse tournaments](trees5_tournaments4_converse_20261005.md) | 3 at the root; 6 versus 10 next | The small correspondence fails at the next size. A neighbouring parameter is a cheap decisive test. |
| [Seven/twenty-one/Mersenne study](seven_twentyone_mersenne_openai_math_20261006.md) | 7,18,19,21,63,128 | The mod-19 and mod-27 clocks share exponent indices because both have order 18; their Fourier kernels still differ. |
| [Ramanujan/Heegner study](ramanujan_heegner7_trivial_cycle_20261006.md) | 7,9,13,163 | Nine class-number-one maximal orders versus thirteen orders: changing the counted category changes the count. |
| [Six-point information observer](information_observers_collatz_20261014a.md) | 6,15 | Six pentads construct Out(S6); source action must be transported separately. |
| [128 atlas](one_twenty_eight_atlas_balance_gold_20261014.md) | 7,8,16,128,240 | Hamming sphere packing and code/lattice constructions explain some counts; poset balance and 2-adic claims need their own maps. |
| [18/19/11/22 merge](merge_18_19_11_22_collatz_hn_20261018.md) | 11,18,19,22,29 | A denominator norm and a multiplicative order are distinct invariants even when the integer agrees. |
| [Pythagorean/conic graph dictionary](pythagorean_form_conic_graphs_rp_margin_20261018.md) | 5,7,10,14,28 | Classical Petersen/Coxeter/Heawood identifications are real graph isomorphisms; the generated group can be smaller than the full automorphism group. |
| [Pinned shade support](shade_difference_support_20261018b.md) | 6,15,30 | All fifteen palette matchings are real charts. Their fibres overlap, so counts cannot be added as disjoint masses. |
| [Carry descent and transport](carry_descent_transport_branches_20261018l.md) | 15,375,4125 | Here 15 counts pairs of six factors, not palette matchings: these are the doily's point and line roles respectively. The branch conditions retain directed factor pairs; see the incoming integration in §9. |
| [Seventeen Vitali atoms](../../07-reflections/seventeen-vitali-atoms.md) | 17 | Seventeen is the author's list length, not a classified universe. Its existing warning restricts the old girth claim. Reuse the partition/choice question, not an asserted universal seventeen. |

The recovered trajectory is useful: early reflections notice numbers and
shape, later atlases ask for a common object, and the natural-bijection and
source-preserving work specify exactly what survives. The next extension is
to use a failed continuation to classify the missing pieces, as in §4.

## 6. A working atlas for numbers other than fifteen

These entries distinguish exact linked counts from candidates and different
roles. All imported claims are scoped to the references in §5; the elementary
formulas and the new finite carrier are checked here or proved above.

| N | Two or more appearances | Structural relation or next decisive test |
|---:|---|---|
| 3 | K4 perfect matchings; nonzero E[2] over an algebraic closure of characteristic ≠2; lines of F2² | The g=1 even-subset quotient gives the matching/line map; an elliptic 2-torsion marking supplies the geometric interpretation. |
| 4 | K4 vertices/faces; four T4 isomorphism classes; four elements of E[2] | Vertices and opposite faces have a natural incidence map. Tournament classes are a different quotient; vertex labels do not act on them as on the four vertices. |
| 5 | Integer partitions of 4; Platonic solids; planar Bravais lattice types | Counts agree, but no common carrier established. Platonic pairs satisfy 1/p+1/q>1/2; lattice types retain metric/centring data. A classification-by-constraints analogy, not an equivalence. |
| 6 | K4 edges; S3 elements; K6 pentads; nonproduct genus-two kernels after a symplectic split | Matching quotients and graph-of-isomorphism kernels explain specific sixes. Pentad action requires the outer twist. |
| 7 | Fano points and lines; nonzero F2³ vectors; imaginary octonion basis directions; P7 vertices | Fano incidence organizes the code and octonion conventions. A basis and multiplication signs are extra data. |
| 8 | Affine F2³ points; one radius-one ball in seven-bit Hamming space; order of phi in the repo's F9 model | First two share the code construction. Field order requires its own characteristic polynomial; equal eight is not a field isomorphism. |
| 9 | Rational class-number-one imaginary quadratic fields; product kernels in a split V2 | Both exact, unrelated classifications so far. Retain maximal versus nonmaximal orders in the CM count. |
| 10 | K5 edges; Petersen vertices; three-fold icosahedral axes | Petersen = Kneser graph on two-subsets of five labels is exact. The A5 axes instead correspond to cyclic subgroups of order three. |
| 11 | Nonzero residues modulo 23 up to sign; coordinate positions in the ternary Golay precursor; allowed cyclic types in Mazur's list | These are differently typed elevens. The actual repo bridge is the M24 frame shape 1²11² to eta(t)²eta(11t)²; no torsion-list bijection supplied. |
| 12 | Icosahedral vertices; T5 isomorphism classes; extended ternary Golay coordinates | Preserve actual action before transferring. The old T5 metagraph is not the icosahedral graph despite matching some totals. |
| 13 | Class-number-one imaginary quadratic orders; residue-clock/geometry appearances | Nine maximal plus four nonmaximal orders (−12,−16,−27,−28). This is an explicit forgetful-fibre example, not thirteen unrelated prime phenomena. |
| 14 | Fano incidence graph vertices; fourteen affine planes of F2³; fourteen weight-four H8 words | First count is 7+7; latter two are the same supports. Marking an affine origin produces the Fano plane. |
| 16 | H8 codewords; Hamming(7,4) codewords; K4 integral cut/cycle index in the repo | Puncturing gives an exact code bijection. The lattice index is a separate determinant invariant and needs an additional map. |
| 18 | Units mod27; nonzero elements mod19; order of 2 in both groups | The generator map 2^a mod27 ↦ 2^a mod19 is a group isomorphism. It does not preserve addition or Fourier kernels. |
| 21 | Fano point-line flags; Paley P7 arcs; Paley P7 automorphisms | Both flag and arc actions of F21 are regular. Choose the Fano lines as translates of {0,1,3} in F7 and base flag (0,{0,1,3}), base arc (0,1): transport by x↦ax+b, a∈{1,2,4}. This is an exact equivariant bijection after those choices; all 441 cases checked. It does not extend to the full Fano group of order168. |
| 24 | Tetrahedral flags; S4; labelled transitive T4; extended binary Golay coordinates | A flag or total order identifies the first three after a base choice. The Golay length arises from a different extension mechanism. |
| 28 | Pairs among eight labels; weight-two classes in V3; Coxeter graph vertices in the conic lane | First two are exactly the same labelled classes. The graph's quadratic-form construction must be retained to import adjacency. |
| 30 | Labelled affine structures on eight points; labelled Fano planes on seven; missing genus-three matching kernels | Exact bridges in §4; 8!/1344 = 7!/168 = 30. |
| 63 | Nonzero F2⁶ vectors; F64* elements; nonzero genus-three J[2] | Additive vector-space identification exists after choices. F64* is cyclic of order63, while J[2] has exponent2; multiplicative and additive laws are not interchanged by the count. |
| 105,135 | K8 matchings; all genus-three Lagrangians | The gap is thirty Fano/Hamming objects. This row is intentionally unequal: it predicts new structure. |
| 128 | Seven-bit words; 16 Hamming balls each of size8; separate Collatz clock periods | Perfect packing proves the first equality. A source-specific Collatz payment rule is additional information, not code coverage alone. |
| 168 | GL(3,2); Fano automorphisms; Klein-quartic automorphisms | A classical exceptional-group bridge. The Paley tournament automorphism subgroup has order21, so passing to the full Fano symmetry loses orientation. |

For the lattice row, [Bright–Cooper–Kurlin](https://journals.iucr.org/a/issues/2023/01/00/uv5012/index.html)
provides a useful refinement of the older five-count slogan: continuous
complete invariants distinguish lattices within symmetry classes. The finite
type is an informative quotient, but not a full description of the object.

## 7. How to extend the technique, and what it can do for active problems

Use the following order when a new N repeats.

1. **Type the count.** Labelled objects, isomorphism types, connected components,
   parameter families, or a finite catalogue? Specify field, characteristic,
   orientation, and allowed degeneracies. Record a failed example too.
2. **Compare a small fingerprint.** Dimension/rank, incidence degrees, orbit
   lengths, stabilizers, fixed-point counts, and a canonical operation. This
   often determines the appropriate comparison before any expensive search.
3. **Attempt an explicit map with a decoder.** Name the preserved predicate,
   information lost, and marked coordinate needed. An equivariant bijection
   requires matching stabilizers up to conjugacy; Lemma N1 allows containment
   only because its weaker target is a bijection of orbit sets.
4. **Extend the parameter one step.** Compare sequences, not isolated entries:
   matching counts 3,15,105 versus symplectic counts 3,15,135; Bell counts
   B4=15, B5=52 versus 2^4−1=15, 2^5−1=31; converse tournaments versus trees.
5. **Explain the discrepancy.** If a map remains injective but stops covering,
   enumerate the complement by an invariant and search for its primitive
   object. Here pair-generated codes omit the H8 block. This is a productive
   obstruction, not a discarded coincidence.
6. **Transport the actual obligation.** A source boundary relation or weighted
   current is richer than the carrier. Verify the transported operation and
   acceptance rule, including tails, before using the analogy as a theorem.

Fixed points of individual group elements are a useful hostile test but are
not a complete test for arbitrary finite G-set isomorphism. One may instead
compare fixed points of **all subgroups** (the table of marks), or explicitly
match transitive orbits by conjugate stabilizers. Thus even agreement of a
whole permutation character must not be silently upgraded to a natural map.

### Concrete follow-through and helper questions

**Six-colour charts — same carrier, additional predicate (OPEN transfer).**
The palette matching in [shade support](shade_difference_support_20261018b.md)
maps exactly to a Lagrangian of V2. A palette pair is a doily point; a chart
is a doily line. This provides a canonical incidence index for the existing
fifteen-chart search. It preserves membership of a colour pair in a chart
and S6 relabelling. It does not preserve actual unit-distance edges unless
they are separately stored, nor does one chart's failure imply all colourings
fail. The source colours, source-to-source path parity, and permitted chart
transitions remain the sidecar. The cheapest useful experiment is to record
obstruction supports on this fifteen-line geometry and check every transition
against the original pinned extension relation. At eight colours, should a
pair-only chart grammar be enlarged by H8-type charts? The complement of the
linear-algebra grammar suggests this question; it does not prove colourability
or forcing. A decoder to actual colour names is the first required test.

**Sporadic sequence transformations — compare operator geometry (OPEN).**
Build a graph whose vertices are the fifteen named recurrences and whose
edges are proved pullbacks, twists, binomial transforms, or symmetric-square
relations. Compare its invariant partitions and automorphism group with a
candidate fifteen-object geometry. A degree mismatch or stabilizer mismatch
would cheaply reject that specific model. The recent derivative/Clausen
work remains at its proved family scope; the catalogue size does not extend
its congruences to untouched families.
The incoming [reanchoring packet](reanchoring_kernel_20261018c.md) already
supplies the A/C/F component with explicit carry transport. Use those edges
as inherited data rather than rediscovering them; §9 records the stronger
marked comparison they require.

**Tilings — retain overlaps and moduli (OPEN).**
For each of the fifteen pentagon families, encode angle and edge constraints,
dimension, and common specializations. Their intersection/degeneration graph
is a more meaningful comparison object than fifteen arbitrary names. Its
incidence may suggest a familiar finite geometry, but that must be calculated.
Mark local boundary ports before importing any tile construction into
Hadwiger–Nelson: tile adjacency alone does not encode unit-distance edges.

**Collatz — universal code coverage is a model of what is missing (OPEN).**
The Hamming balls prove coverage because they are disjoint and their exact
masses sum to the whole universe. A Collatz certificate family needs an
analogous decoder plus a terminating or well-founded acceptance rule for
each integer source. Binary residues can encode finite paths exactly while
their inverse limit still admits nonterminating addresses. Therefore the
new helper question is: can a source-preserving refinement be decomposed
into local pieces whose uncovered complement has an independently decreasing
rank, as weight-two peeling did for the code? Keep the original source,
signed weight and integer clock. The present code construction supplies a
worked model of completion by missing pieces, not a Collatz rank function.

**Planar Jacobian — preserve local carrier data (OPEN).**
Moving between curve branch points, 2-torsion and code coordinates retains
a finite symplectic pairing. It does not retain a rational map's poles,
critical-point behaviour, or recovery of original polynomial observables.
A useful helper question is which of these source predicates can be encoded
as marked data on a degeneration graph before composing transfers. Equal
kernel counts alone cannot reopen an excluded polynomial carrier.

## 8. Reproduction, repairs, and boundaries

From the repository root:

```bash
python3 04-computation/experiments/cardinality_carriers_20261018m.py --output /tmp/cardinality.json
python3 -O 04-computation/experiments/cardinality_carriers_20261018m.py --output /tmp/cardinality-O.json
cmp /tmp/cardinality.json /tmp/cardinality-O.json
```

The computation has no random seed and no unreported filter. The second
construction is the affine-plane permutation orbit, independently of the
isotropic-subspace enumeration. A development check caught one representative
error: a raw duad containing the last label is not the chosen representative
of its complement class. Incidence now canonicalizes it first; this is also
recorded in the mistakes ledger. The receipt does not certify arbitrary
genera, geometric realizations over Q, or any active conjecture.

The historical fifteen note is repaired for MISTAKE-059's nonuniform fibres,
the split of the 24 five-cycles into two A5 conjugacy classes, and the sign of
the spherical defect. The five-solids/five-lattices reflection receives a
scope warning: its count analogy is not a derivation of LRC or an equality
of classification problems. The existing seventeen-Vitali warning remains;
seventeen is treated here as an editorial list size.

The genuinely reusable research move is now explicit: **count → retain an
action and relation → test the next size → identify the missing primitive →
verify the original predicate**. The archive already contains both successes
and sharp failures of each step. The genus-three Hamming complement supplies
a compact, completely checked example of the whole move.

## 9. Incoming work: two fifteens, two sixes, and a retained first response

Commit `63fab03129` arrived during close-out. Its
[Zagier parameter atlas](zagier_six_triples_20261018c.md),
[reanchoring kernel](reanchoring_kernel_20261018c.md),
[coupled branch compiler](coupled_support_branches_20261018c.md), and
[fractional geometry](coupled_fractional_geometry_20261018c.md) sharpen the
present comparisons. The first two proofs were read and their optimized
checkers independently invoked here: **9,938 and 10,704 checks**, with
byte-identical saved receipts. The large geometric census below is inherited
from its authenticated packet, not recomputed in this session.

**The two fifteens in the active colouring work are distinct roles.** The
375 local branches choose two of six direction factors and one of five cosets
for each: binom(6,2)·5². Its fifteen objects are **duads**. The palette charts
pair all six colour names: 5·3·1. Its fifteen objects are **synthemes**.
After identifying their abstract six-label sets, the doily relates them by
incidence, not by a relabelling-equivariant point/line identification. The
fixed-point witness 7 versus 3 in §3 applies immediately. A line contains
three duads, every duad belongs to three lines, and the resulting 45 flags
are a common incidence refinement. A branch's cosets and directed factor
constraints still have to be retained; these flags do not create a branch
equivalence. This is a concrete reason to preserve the two meanings of 15
in the research data schema.

**The two sixes require a common frame.** The incoming packet compares the
six K4 edges with the six generic cross-ratio expressions of four ordered
points. The former S4 action is faithful; the latter factors through S3 with
kernel V4. A double transposition fixes every cross-ratio expression but
moves four edges. Both are quotients of 24 ordered frames, by different
four-element stabilizers. Retaining a frame gives a correspondence; silently
identifying the quotients discards different information. This gives a second
repair move beside the missing-H8 completion: **lift to a shared marked
space before comparing two incompatible quotients**.

**The sequence transformation must carry its first response.** For the
Zagier convention of §2, put P(t)=t(t²−at+c). At a root d of P and nonzero
scale epsilon, reanchoring t=d+epsilon·z sends

    (a,b,c) -> ((a−3d)/epsilon, (b−d)/epsilon, P'(d)/epsilon²).

The corresponding sequence transform is

    u'_n = epsilon^(−n) sum_(k=0..n) binom(n,k)(−d)^(n−k) u_k.

For A=(7,2,−8), roots d=−1 with epsilon=1 and d=8 with epsilon=−1 give
C=(10,3,9) and F=(17,6,72). These are actual edges in the proposed recurrence
graph. The first-block divided carry transforms **affinely**, with a correction
from binomial and power carries. At p=7, r=1, A→C has actual residue5, naive
residue6 and correction6; 6+6=5 modulo7. A→F has actual2, naive6 and oriented
correction3. The correction is the retained first response, exactly the
kind of coordinate a count or zeroth-order polynomial forgets.

The golden triple (11,5,125) in that packet means
(L5,F5,5F5²), with signed norm −1. It is not the recurrence parameter triple
(a,b,c); that literal interpretation gives u2=5/2. Likewise F's accessory
b=6 cannot be changed to5 while retaining integrality, although its singular
set and auxiliary elliptic j-invariant stay unchanged. This connects the
number atlas to golden arithmetic through an explicit slot decoder.

**The geometric gate remains complete boundary support.** The new 714-point
configuration retains 94 mixed edges inside its relative six-core, but all
570 boundary patterns still extend. Its 8,550 fixed-pair chart tests all fail,
already forced by a same-anchor unit triangle. Thus the next meaningful
colour representation permits a variable ternary pair label and retains
the binary shade and both coordinates of every source pin. Its lift condition
is bipartiteness of each same-label edge graph plus pinned path parity.
An arbitrary cycle basis is insufficient: in K4, giving one vertex label1
and the other three label0 makes the three triangles through the first
vertex mixed while the remaining triangle is monochromatic. The predicate
is nonlinear even though the cycle space is linear.

The resulting research sequence has three distinct repairs: complete a
carrier by its missing primitive objects; lift incompatible quotients to
their common frame; and transport a first-response correction with the
compressed object. Each is useful because it identifies the data a proposed
analogy would otherwise erase.
