# The exact odd index of tournament migration

The Berggren tree for primitive Pythagorean triples admits unimodular matrices on its primitive coordinates and a ratio-selected unique inverse. The repository's earlier Stern–Brocot/Berggren ordinal work (THM-4057) and labelled local T4 work (THM-3364) provide the relevant distinctions between exact parent maps, modular shadows, and analogy. Here the rooted tournament migration has an equally explicit inverse on each actual branch, but a genuine extra odd index. This is a structural no-go for a unimodular identification of this specific migration alphabet, not a no-go for a different tournament pair construction or for Chan–Pak.

## Sharp arithmetic index

Use the rooted transfer matrix for a migration word of width k and integer u:

    M(w) = [[F(u),u+1],[2^k−rev_k(u),2^k]],

where F(0)=1, F(2x+1)=F(x), F(2x)=F(x)+x+[x odd] for positive x. On rooted states with one full terminating path and S=0, this is the physical (H,R) matrix from the branching-reader theorem. The norm statement below is purely about its matrix, so applies also to the usual leading-one terminal reads when used on compatible states.

Remove all leading zeroes and trailing ones. If nothing remains, the determinant is zero. Otherwise the remaining **active** word starts in 1 and ends in 0. If it has a ones and b zeros (a,b≥1), then

    d = det M(w_active) is odd, with d ≥ 2^a+2^b−1 ≥ 3.

Equality in the first inequality holds exactly for the two-run word 1^a0^b. Thus at fixed active length k≥2 the sharp minimum is 2^floor(k/2)+2^ceil(k/2)−1; all minimizing words are the balanced two-run words. In particular an upper bound on odd index gives a logarithmic upper bound on active width, regardless of padding.

Here is a direct proof. For any nonempty prefix let D=2^width, v its reversal, x=v/D, and ε its last bit. Under appending a 1 the normalized determinant δ=det(M)/D stays fixed; under appending a 0 it changes by ε−x/2. In both cases x becomes (x+b)/2 for appended bit b. This is obtained by substituting u'=2u+b, v'=v+bD and the recursion for F. It also holds for the empty-prefix initialization at δ=0, x=0, ε=0 (both unpadded extensions have determinant zero). A nonempty zero run of length l after ones therefore contributes exactly 1−x(1−2^(-l))>0 to δ, where x is the normalized reversal just before the run. For the first runs 1^a1 0^b1, this contribution is 2^(-a1)+(1−2^(-a1))2^(-b1). Multiply by the final D. With A=2^a, B=2^b, U=2^(a−a1), V=2^(b−b1), one obtains BU+AV−UV≥A+B−1 since it is increasing in U,V on 1≤U≤A/2, 1≤V≤B/2. If there is any further run the inequality is strict. Parity follows from even u and odd reversal in the active word.

Leading zeroes and trailing ones have exactly separate dyadic roles. For w=0^g·w_active·1^t,

    M(w) = diag(1,2^g) · M(w_active) · diag(1,2^t).

Each nondegenerate step therefore multiplies the odd part of the lattice index by at least 3; arbitrary sequences or dyadic renormalization cannot make it unimodular. Since F(u) is odd and the bottom-right is a power of two, the gcd of all matrix entries is 1. The Smith normal form is (1,det), so the full finite cokernel is cyclic. For the active word, D is a unit modulo its odd d: an output (h',R') is in the image lattice iff

    D h' − (u+1)R' ≡ 0 (mod d).

Indeed the other inverse numerator multiplied by D is d R' minus the bottom-left entry times the displayed numerator, hence automatically divisible by d. General padded words require first undoing the physical output-R scaling, then this odd test and finally the input-R scaling. These are precise branch-image conditions; physical parent reachability is a separate test.

## Dyadic redistribution is a literal graph gauge

For two consecutive migrations, the chronological words (0a,b) and (a,b1) produce isomorphic rooted tournaments with the same final order from any common initial rooted tournament. The first zero-labelled spine vertex of the older stage dominates everything else in that older stage, including its new root. Move that vertex to the end of the next spine. It already lost to every vertex on that next spine, and the next root already beat it as an old nonroot vertex; its new terminal bit 1 preserves the latter arc. It still dominates every vertex in the remaining old graph. All other arcs agree under this relocation. This proves the gauge on literal labels without any assumption on intermediate numeric root state.

For example the previous slope-three block (width,address)=(4,2),(8,84) agrees not just numerically but as a rooted tournament with (3,2),(9,169). This is precisely the graph mechanism for transferring a dyadic factor between adjacent matrices, (M(b)M(0a)=M(b1)M(a)). A unique parent language needs a normalization or a fixed alphabet; letting arbitrary subdivisions count as distinct histories can already double-count a single physical host.

## A literal failure of primitive preservation

Start with the transitive pair rooted at its last vertex, state (H,R,S)=(1,2,0). A weak migration (4,6) (word 0110) gives (19,42,0), still primitive. The weak migration (3,2) (word 010) then gives (183,450,0), whose gcd is 3. The replay counts both literal rooted graphs and independently checks that the unique full ending at the root survives both moves. Every later integer linear migration on these compatible coordinates preserves the common divisor; any subsequent positive read of these coordinates is a multiple of 3. This kills hard prime/7p values with p≠3 on that continuation and contrasts with Berggren's universal primitive preservation. It does not say all branches lose primitivity.

Run `python replay.py` for exhaustive finite controls of the determinant formula against 16,278 active padded words, literal isomorphism under four redistributed padding pairs, and the direct rooted gcd-loss witness. Normal and optimized transcripts agree.

## An exact physical sum of two squares

Let w* reverse and complement the bits of w. The same binary word defines a one-apex tournament over a transitive spine; reversing the spine and all arcs shows that its full count F is invariant under w→w*. Equivalently its inversion formula is invariant because reversing and complementing an earlier-1, later-0 pair gives another such pair at the same distance. Further,

    u*+1 = 2^k−rev_k(u),    2^k−rev_k(u*) = u+1.

Hence M(w*) is the transpose of M(w). Pad with both a leading 0 and a trailing 1, so the transpose identity holds for a=0w1 and a*=0w*1, and both are physical weak migrations preserving S=0. From the rooted transitive-pair state v=(1,2), apply a then a*, obtaining (h₂,R₂)=M(a)^T M(a)v. The first migrated state is (h₁,R₁)=M(a)v.

Use the converse of the final rooted graph as the incoming block of the single-return reader, with a singleton module and singleton outgoing vertex. The latter beats only the distinguished reader vertex among the incoming block. Its exact read is h₂+2R₂ by the two-run formula. This final tournament is strong even when the rooted stages are weak: in the converse, the incoming marked root starts a full Hamilton path, so it reaches every incoming vertex, and every vertex reaches the outgoing vertex and then the marked root.

Therefore the literal Hamilton count is

    h₂+2R₂ = v^T M(a)^T M(a)v = h₁²+R₁².

For w=10,100,1010 the actual resulting counts 1669,6709,25309 are prime by trial division. The smallest graph was independently counted by raw subset DP as well. An infinite primitive subfamily is w=1^j0, j>=1. Put A=2^j. Directly F(w)=A+1 and the symmetric-padding rule gives

    h₁=9A−3,      R₁=18A+2,       R₁−2h₁=8.

Since h₁ is odd and R₁ even, these coordinates are coprime and opposite parity; they yield the actual primitive Pythagorean triple

    (R₁²−h₁², 2h₁R₁, R₁²+h₁²),

whose hypotenuse is the literal tournament count above. This is a concrete link to a subfamily of Berggren's arithmetic object, not an identification of all ternary moves.

There is also a sharp arithmetic limit on this diagonal pairing: if 7 divides a sum of two squares, then both entries must be divisible by 7 because −1 is not a quadratic residue modulo 7. Repeating the argument gives an even valuation at 7. Thus this entire sum-of-squares chart cannot furnish 7p for an odd prime p distinct from 7, even though it does furnish actual primes. Breaking the transpose pair or changing the final subscriber leaves the scoped theorem; the obstruction is not to the full tournament spectrum.
