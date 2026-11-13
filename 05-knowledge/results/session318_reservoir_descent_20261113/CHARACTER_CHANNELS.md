# Reflection characters, branch-sensitive lifts, and an exact obstruction to descent

Provenance boundary: THM-4731 already has twin doubling, affine self-line Burnside, and outer-ring symmetric depth. THM-477's blue code and THM-452's doubling coset law supply the linear/affine vocabulary. Session22's converse channels already give the unrestricted self-converse enumerators via odd partition Burnside. The additional layer here is the joint reflection/translation decomposition of the *actual weighted class adjacency*, and its exact compression obstruction.

Let X be the n-vertex path tiling cube, n>=3, F=Fix(sigma), D=im(1+sigma), and let c(x) be the ordinary tournament class. Write rho for converse on classes. For v in F define a dart-count matrix A_v by counting x with c(x)=a and c(x+v)=b. A_0 is the positive fibre-size diagonal M. For v=kappa (the all-ones vector), this is the ordinary metagraph adjacency (loop has two darts); v can be any reflection-invariant direction. All A_v commute with rho and are symmetric. Matrices for v in F commute on the cube, not necessarily on classes.

## Fourier/positive decomposition retaining every class edge

For each d in D choose one x_d with (1+sigma)x_d=d; its coset x_d+F is the entire syndrome fibre. On this coset sigma translates by d, not by a universal deck element. For a real character chi:F -> {+1,-1} define the column vector indexed by ordinary classes
```
b_(chi,d)[a] = sum_(f in F, c(x_d+f)=a) chi(f).
```
Changing the representative multiplies the whole column by a sign and changes no outer product. Directly substituting v=f+g in the two sums proves
```
G_chi = sum_(v in F) chi(v) A_v = sum_(d in D) b_(chi,d) b_(chi,d)^t,
A_v = |F|^(-1) sum_chi chi(v) G_chi.
```
Every G_chi is positive semidefinite. Moreover b_(chi,d)[rho(a)]=chi(d)b_(chi,d)[a], by reflection on its own coset. So the sign block of G_chi is the Gram sum over exactly d with chi(d) equal to that sign. Its rank is at most the number of those d. A character trivial on D has no antisymmetric block at all; otherwise both block ranks are at most |D|/2. This is actual adjacency recovery, not just blue/black or degree statistics.

In particular A_kappa is a signed sum of positive matrices, not a positive operator. A useful sanity check: transitive node has one non-loop dart, so a 2-by-2 principal minor of its row and neighbour has negative determinant. The symmetric Markov version M^(-1/2) A_v M^(-1/2) is the compression of the translation unitary by the isometry lifting a class basis vector to its normalized indicator. Hence its spectrum lies in [-1,1]. The original cube is a permutation representation; compression is not a representation of F. The replay gives a minimal failure at n=4: translating the sole altitude coordinate {0,3} and translating the mirrored pair {{0,2},{1,3}} gives noncommuting normalized class operators. This is why a common upstairs Fourier basis does not diagonalize the metagraph downstairs.

## Boundary-sensitive voltage factor

For a non-self-converse merged vertex choose one lift a and its mate rho(a). Between two such vertices a,b the even block of the normalized A_kappa adds A[a,b]+A[a,rho(b)] and the odd block subtracts, then both divide by sqrt(M[a,a]M[b,b]). This applies on the diagonal too, with the dart convention. Between a two-sheet vertex and a self-converse vertex c the even entry is sqrt(2) A[a,c]/sqrt(M[a,a]M[c,c]); the odd block has no such vertex or entry. Between two fixed vertices just retain their ordinary normalized entry. Use normalized vectors (e_a±e_rho(a))/sqrt(2) and e_c to prove the formulas. The characteristic polynomial factors exactly into these two blocks. The odd block is a signed or C2-voltage graph on the two-sheet locus with the boundary deleted, not the naive double cover of the full merged graph. Its trace moments count closed walks avoiding the fixed locus with a factor +1 for returning to their chosen sheet, -1 for ending on the opposite sheet, using the indicated normalization. Gauge changes are row-and-column sign conjugations.

More precisely, if P is the nonnegative even block restricted to the two-sheet locus, the entrywise absolute value of the odd block is bounded by P. Its spectral radius is at most rho(P). In a connected component of P with some edge to the fixed locus, rho(P)<1: P is similar to the corresponding restriction of a stochastic walk on merged vertices (use positive eigenvector sqrt of aggregate fibre mass); a finite connected substochastic chain with a leaky row loses mass within a uniform finite time. This proves strict loss in the odd sector without estimating voltage cancellation. Even if there is no boundary edge, strict loss can still come from competing signed parallel edges. No assertion is made that all orders have either property.

## Outer-ring representation recurrence, with honest Bott boundary

In n->n+2 central inclusion, new coordinates consist of one fixed coordinate and n-1 reflected pairs. Therefore F gains n independent binary generators, D gains n-1, and syndrome and character pairs update by direct products. On a decomposed pair the eigenvalue multiplies:
```
chi_new(d_new)=chi_old(d_old) chi_ring(d_ring).
```
Also chi_new(kappa_new)=chi_old(kappa_old)chi_ring(kappa_ring). Since kappa has a 1 on a fixed coordinate, it never lies in D. At syndrome zero there are |F|/2 characters for each kappa sign, all reflection-even. At each nonzero syndrome, d and kappa are independent in F, so there are |F|/4 for each of the four joint signs. Thus joint multiplicities of the uncompressed presentation representation are |F|(|D|+1)/4 for either kappa sign in the reflection-even sector and |F|(|D|-1)/4 for either kappa sign in the reflection-odd sector (n=3 odd sector zero). This gives corresponding upper bounds on positive and negative inertia of the compressed normalized operator in each sector, by restriction of quadratic forms. This exact tensor recurrence does not imply periodicity of class adjacency, or that the class indicators factor along an outer ring. It is a finite C2 character decomposition and tensor stabilization, not Bott periodicity, K-theory, or a topological dimension-mod-eight theorem. A genuine category analogy would retain the class incidence as a correspondence, not pretend it is an equivariant map respecting all translations.

Replay: `python replay_characters.py` checks n=3..6 against literal path-pinned class identification, independent adjacency counts, reconstruction for all translation characters, parity columns, and noncommutativity of normalized kernels at orders 4–6. This is not a theorem of Chan–Pak coverage.

## Where a wreath product really occurs: cyclic modules

Do not call the reflected ring stabilization itself a wreath product. There is a genuine one on a distinct layer: canonical cyclic substitution modules. For nonempty tournaments A,B,C in the directed order of C3, let T=C3[A,B,C]. The three blocks are precisely the maximal proper modules. Indeed, a module meeting two blocks must meet the third (any vertex in the third distinguishes those two), then must contain every vertex in every block (a missing vertex sees the other two blocks oppositely). Therefore automorphisms rotate the quotient blocks, and anti-automorphisms reflect them. Write a_X for the automorphism count, j_X for the number of anti-involutions (permitting the identity on a singleton only, because there are no arcs to reverse).

For each of the three quotient reflections, call its fixed block F and swapped blocks B,C. It is allowed exactly if F is self-converse and C is isomorphic to B^op. Its contribution to the anti-isomorphism count is a_F a_B^2 and its contribution to the anti-involution count is j_F a_B: the two anti-maps on the exchanged pair are independent in the first count, and must be inverses in the second. A valid reflection is unique unless all three block classes coincide and are self-converse; two distinct valid reflections compose to a rotation and hence force equal classes. In the fully repeated case the automorphism group is the actual permutation wreath product Aut(A)^3 ⋊ C3, and the three counts (automorphisms, anti-maps, anti-involutions) are (3a^3, 3a^3, 3aj). In a homogeneous depth-d tower over a self-converse base with a automorphisms and j anti-involutions, the recurrence gives
```
a_d = 3^((3^d-1)/2) a^(3^d),
j_d = 3^((3^d-1+2d)/4) a^((3^d-1)/2) j.
```
This counts vertex reversal symmetries and not Hamiltonian path involution-fixed marks. A Hamiltonian path may cross the three blocks repeatedly; wreath symmetry does not justify treating its blue path count as a product of blue counts for blocks.
