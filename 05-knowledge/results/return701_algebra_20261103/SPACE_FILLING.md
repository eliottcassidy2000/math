# Certified Collatz addresses, Hilbert charts and binary algebra

**PROVED constructions and identities; FINITE-EXACT replay through depth 7.**
Local labels H1--H4. This develops the user's space-filling idea using
ancestry as cell inclusion. It does not require native Collatz edges to
be nearest-neighbour planar edges.

## H1 A source-preserving arithmetic address

Recover the proper inverse binary grammar SF5 from
[COLLATZ_LATTICE](../joint_incidence_clocked_frames_20261102/COLLATZ_LATTICE.md).
Start at 5. For an odd positive n not divisible by 3, take candidates
`k=2,4,6` if n=1 mod 3, and `k=3,5,7` if n=2 mod 3. Exactly one candidate
has `2^k n=1 mod 9`: successive candidates traverse the three lifts of
the required class modulo 3. Discard it. The remaining two candidates,
in increasing order, are the children labelled 0 and 1:

`m=(2^k n-1)/3`.

They are odd, not divisible by 3, and larger than n. The ordinary Collatz
map takes exactly k+1 steps from m to n, with k+1<=8. The valuation
`v_2(3m+1)=k` recovers both parent and label. Induction gives distinct
integer sources for all binary words at any fixed depth. Every such
source has a certified route to 5 and hence to 1.

Keep the full clock, not just the word length. If L inverse edges have
exponents k_1,...,k_L, let K=sum k_i. Then

`n_L=(2^K*5-A_L)/3^L`,
`A_0=0`, `A_(L+1)=2^(k_(L+1))*A_L+3^L`,
`ordinary clock to 5 = K+L`.

No assertion that every integer belongs to this binary subtree is used.

## H2 A recursive square chart preserves ancestry and scan order

Group two binary letters into one base-four letter. For every d, the
4^d words of length 2d index all cells of a 2^d-by-2^d grid. Use the
classical Hilbert recursion. At a step with side-half s, the four maps
of the preceding grid point (x,y) are

`(y,x), (x,y+s), (x+s,y+s), (2s-1-y,s-1-x)`.

The first base-four digit selects the map and the remaining digits
select its argument. Induction proves bijectivity and unit Manhattan
distance between consecutive scan indices: each subcurve is an isometry
of the smaller one, and the three joins connect adjacent endpoint cells.
It also gives

`floor(point_d(i)/2)=point_(d-1)(floor(i/4))`.

Thus removing the last two inverse Collatz instructions removes precisely
the final spatial refinement. Each cell retains its actual integer,
binary address, orientation and native step clock. This is a genuine
bijection of finite source sets with grids at every scale.

The nested closed squares have diameter sqrt(2)*2^(-d), so every infinite
address determines a point. Consecutive parameter cylinders map into
adjacent squares. If `4^(-d-1)<|t-u|<=4^(-d)`, the points lie in at most
two adjacent depth-d squares, giving distance at most
`2*sqrt(5)*sqrt(|t-u|)`. Consistent endpoint limits identify the two
base-four expansions at a boundary. Every spatial point has a nested
cell chain, proving surjectivity. This recovers a continuous Hilbert
space-filling map, now carrying the arithmetic labels on its finite cells.

The underlying curve is classical: [Hilbert, 1891, original paper](https://www.digizeitschriften.de/download/pdf/235181684_0038/log31.pdf).
The construction here is the specified source/clock transport, not a
claim to have invented that curve.

Scan adjacency is different from a Collatz edge. A Collatz parent is
represented by a containing cell at a different resolution. This
explicitly respects the earlier obstruction: 704 sources within 21
ordinary steps of 5 cannot inject into a honeycomb ball of 694 vertices
with each native edge a lattice edge. See
[GEOMETRY](../collatz_hexagonal_growth_20261102/GEOMETRY.md). The new chart
changes the retained geometric relation; it does not invalidate the bound.

## H3 Cayley--Dickson reduction becomes an address automaton

The prior exact result [C7, PROOFS](../cayley_dickson_integral_20261030/PROOFS.md) identifies the integral
Cayley--Dickson algebra at stage m after reduction modulo 2 with

`F2[(C2)^m] = F2[t_1,...,t_m]/(t_1^2,...,t_m^2)`.

The signed basis product becomes `e_i e_j=e_(i XOR j)`. Take m=2d and
index that basis by the binary inverse words in H1. Multiplication by a
basis element is now an exact permutation of certified arithmetic sources
and their Hilbert cells. This is address transport: XOR is not asserted
to be ordinary addition of the source integers, and the characteristic-2
algebra is not identified with the real nonassociative tower.

For the final two address bits, `i XOR 1` and `i XOR 3` are the two
edge-neighbours of i in its final 2-by-2 block; `i XOR 2` is the diagonal.
For a binary cell field v, define

`D_a(v)_i=v_i+v_(i XOR a)` over F2.

Then `D_1^2=D_3^2=0`, while
`D_1 D_3(v)_i=v_i+v_(i XOR 1)+v_(i XOR 2)+v_(i XOR 3)`.
The two individual nilpotent responses therefore leave a nonzero joint
response. Each update is local if the cells carry their two marked
compass ports. There are only eight possibilities for those oriented
perpendicular ports. Higher-bit XOR operations act on growing blocks;
there is no claim that the entire tower has a fixed geometric radius.

This supplies an explicit cellular automaton on a marked grid. The full
augmentation ideal has nilpotence index 2d+1, since the product of all
2d different t_i is nonzero and every longer monomial repeats a variable.
Moreover the scalar coefficient of `x e_i` reads the coefficient x_i.
For arbitrary linear combinations x, future basis-multiplication queries
consequently distinguish all 4^d linear coordinates, exactly as in the
inherited C7 minimality argument. A single basis address still needs only
2d bits; the lower bound for arbitrary vectors concerns a larger source set.
Square-zero individual responses do not justify throwing away joint data.

## H4 A finite orientation controller still needs arithmetic memory

The Hilbert recursion's linear orientation part is generated by swapping
coordinates and its negative; there are four orientation matrices.
Given a single length-2d address, its integer source and clock can be
computed by an explicit tape algorithm together with that finite controller.
At depth d all stored
integers have O(d) bits because each inverse step multiplies by at most
128 before division by 3. A direct implementation uses O(d) space and
O(d^2) bit operations, including computing the Hilbert coordinates.

This separates the finite controller from its unbounded source register.
Knowing a residue modulo 9 selects the next children, but computing their
new residues modulo 9 requires the old source modulo 27, and so on.
There is no finite-state claim for the complete arithmetic label process.
This is compatible with the repository's guarded FRACTRAN/Turing-machine
compilers: a visual scan does not remove the tape, priority or source guard.

The independent finite replay checks every source, its ordinary clock,
the grid bijections, all scan joins, ancestry and the XOR identities
through 16,384 sources. Full source/clock/coordinate hashes and clock
histograms are in [hilbert_collatz.json](hilbert_collatz.json).
