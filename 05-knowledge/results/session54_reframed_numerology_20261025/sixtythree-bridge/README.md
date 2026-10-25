# Paley-7 blockers as the punctured binary augmentation module

Status: PROVED explicit equivariant bridge conditional only on the independently reproduced FINITE-EXACT blocker classification. The census is checked exhaustively by two separate Hamiltonian-path algorithms. No Collatz or LRC conclusion follows.

## Recovery and scope

The closest prior mechanism is THM-4553, `01-canon/theorems/THM-4553-the-point-stripping-ladder-psl27-borel-torus-and-the-collatz-frobenius.md`, part (v): beta(P7)=6 and 63 minima, grouped as 21+21+14+7. MISTAKE-567 (in the mistakes ledger around the original claim) corrects the false statement that all 63 are Hall obstructions: the seven isolated-vertex stars retain a one-path-cycle factor. MISTAKE-607 warns against conflating Paley support with Mersenne denominators. The old `mod18_mod19_seven_sixtythree_fractal_20261006.md` explicitly labels its comparison with the Sierpinski tower numerology.

Targeted searches for blocker/augmentation, blocker/F8, even-subset/blocker and 63/equivariant did not locate this bridge in the supplied checkout. This is a recovery-plus-new-map claim, not a claim that binary augmentation modules or their Fourier decomposition are new. The historical referenced computation files are absent from this checkout; this package does not silently depend on them.

Concept board: (1) actual HP blockers; (2) Hall versus exotic distinction; (3) five G-orbits; (4) nonempty even subsets; (5) two F8 Fourier channels; (6) the different Sierpinski G-set. The decisive test is the orbit profile, followed by a type-aware map and hostile incidence/XOR tests.

## Objects and exact map

Let V=F7, Q={1,2,4}, and u->v when v-u is in Q. Let G consist of x->ax+b with a in Q and b in F7. Let B be the minimum six-arc HP-deletion sets. Let W be the even subsets of V, with symmetric difference addition, a six-dimensional F2-space. All neighborhood statements below refer to the original tournament unless explicitly residual.

For B in B define Phi(B) as follows:

1. If the residual digraph has two sources u,v, put Phi(B)={u,v}.
2. If it has two sinks, order them u->v in the original tournament, and put Phi(B)=V\{u,v,2v-u}.
3. If three residual vertices all have the same sole out-neighbor c, put Phi(B)={c} union N+(c).
4. If three residual vertices all have the same sole in-neighbor c, put Phi(B)={c} union N-(c).
5. If B is the six-arc star at c, put Phi(B)=V\{c}.

These cases have sizes 21,21,7,7,7. The common-neighbor center is unique. The checker certifies the classification and uniqueness is also visible from the residual neighborhoods.

### Why this is a bijection

The first family gives all 21 pairs; the last gives all seven six-subsets. The two closed-neighborhood families give fourteen four-subsets. Their complements are the seven out-neighborhood cyclic triangles and the seven in-neighborhood cyclic triangles, which are distinct: each neighborhood triple has sum 3c modulo seven, recovering its center; opposite neighborhoods at the same center are disjoint. The other 21 four-subsets complement transitive triangles.

The triple {0,1,2} is transitive. Its G-stabilizer is trivial, since an automorphism of a transitive triple fixes each vertex and a nonidentity affine map cannot fix three vertices. Its G-orbit therefore consists of all 21 transitive triangles. The map from an arc (u,v) to {u,v,2v-u} is precisely this orbit parametrization. Thus case 2 fills exactly the missing 21 subsets. These five image classes are disjoint and exhaustive.

Every formula commutes with x->ax+b; orientation and neighborhoods are preserved and 2(av+b)-(au+b)=a(2v-u)+b. Hence Phi is a G-equivariant bijection B -> W\{0}, not merely a match of cardinalities. All three seven-orbits have order-three stabilizers; the two 21-orbits have trivial stabilizer.

This choice is explicit and geometric, but not uniquely forced by equivariance. There are 3! times 2! times 21^2 = 5292 G-set bijections in total, since each seven-orbit has trivial equivariant automorphism group and each regular orbit has 21.

## Fourier channel proof

Use F8=F2[alpha]/(alpha^3+alpha+1), with alpha primitive of order seven. For an even subset S set

    F(S) = (sum_{x in S} alpha^x, sum_{x in S} alpha^(-x)).

This is an F2-linear isomorphism W -> F8 direct-sum F8. Indeed the subset polynomial f has degree at most six; if both coordinates vanish, its roots include the two disjoint cubic conjugacy classes of alpha and alpha^-1, and evenness supplies root 1. Thus x^7-1 divides f, forcing f=0.

For g(x)=2^k x+b, the action on coordinates is

    (u,v) -> (alpha^b u^(2^k), alpha^(-b) v^(2^k)).

The two nonzero axes have seven elements each. On mixed pairs, z=uv is fixed by translations and squared by the order-three multiplier. Nonzero F8 elements split under squaring into {1}, {alpha,alpha^2,alpha^4}, and {alpha^3,alpha^6,alpha^5}, giving mixed orbit sizes 7,21,21. Under the displayed Phi, common-out is the first axis, common-in the second, stars have uv=1, two-sinks have uv in the first three-element class, and two-sources have uv in the other.

This proves a structural 63 = 7+7+(7+21+21). Transporting addition through Phi after adjoining a formal zero gives an abstract binary vector-space structure. It is NOT native arc symmetric difference.

## Loss ledger and hostile probes

Preserved: affine G-action, exact object identity via the displayed inverse type table, orbit/stabilizer data, and a transported module structure. Lost without sidecar: native arc incidence, arc XOR, arc Hamming distance, and HP counts of general deletion sets. Sidecar: the five-type decoding rule. There is no claim about deletion sets outside the 63 minima.

1. The natural incidence boundary, XOR of the endpoint pairs of deleted arcs, is NOT Phi. It has only 21 images, each with exactly three preimages. Each 21-element blocker orbit collapses to a seven-element boundary orbit; the three seven-element blocker orbits have a common seven-element boundary image. This is an actual lossy map, not an obstruction to the type-aware bijection.
2. Native linearity is impossible for ANY identification of these 63 arc masks as the nonzero words of a six-dimensional binary code: total nonzero weight would be 63*6=378, whereas every active coordinate contributes 32 and hence total weight must be divisible by 32. It is not.
3. There is no isometry from native arc Hamming distance to even-subset Hamming distance. Arc masks 63 and 2883 have distance eight; even subsets of seven vertices have distance at most six. Under Phi these two masks have distance two. Bit positions are the lexicographically ordered 21 directed edges in the JSON/code.
4. Native XOR closure fails; result.json gives an explicit pair and its non-blocker XOR.
5. The original Sierpinski T6 carrier has eight seven-orbits and seven G-fixed vertices, whereas B has no G-fixed vertex and two 21-orbits. Therefore the earlier Sierpinski/blocker comparison remains non-equivariant numerology for those specified actions. The present target is the augmentation module, not that tower action.
6. No HP-count, Collatz-cycle, or individual-trajectory statement transfers through the F8 dictionary. In particular, no ensemble-to-orbit inference is made.

## Reproduction and evidence

Run from any directory:

    python /path/to/session54-sixtythree-bridge/verify.py
    python -O /path/to/session54-sixtythree-bridge/verify.py

The checker uses explicit exceptions, not assertions. It enumerates all 5040 vertex permutations to obtain 189 Hamiltonian paths, then checks each of 20349 five-arc sets and 54264 six-arc sets both by path-mask hitting and an independent endpoint-subset DP. Every candidate agrees; no five-set blocks; exactly 63 six-sets block. It verifies all 63*21 equivariance equations, all 64*21 Fourier action equations, the full coordinate bijection, orbit stabilizers, and hostile boundaries. Normal and optimized output are byte-identical. Full 63-entry map and representative/stabilizer data are in result.json. source-pins.json records input and artifact hashes.

This package makes no repository edits. Parent integration and independent review remain separate responsibilities.
