# Gaussian disks: eightfold symmetry, smaller encodings, and the missing quartic

## Scope and inheritance

PROVED below: exact orbit counts, full lattice-isometry obstruction, primitive-ray decoding, rational shell transport, and the quartic orbit invariant. These are elementary statements about specified carriers, not Bott periodicity or an equivalence of Gaussian integers and tournaments.

Inherited source snapshot: `5258d52ac95924d524ee816a7dffa7329ddcf095`. Relevant prior work:
- `05-knowledge/results/collatz_bott_marked_periodicity_20261007d.md`: marked stabilization retains data, while homotopy/stable classes erase arithmetic guards; it already separates actual Bott periodicity from arithmetic return times.
- `01-canon/theorems/THM-3339-fibonacci-three-ray-berggren-transplant-and-moving-owner-obstruction.md`: the Fibonacci Pell locus, three Berggren rays, and moving-owner obstruction.
- `01-canon/theorems/THM-3536-berggren-angle-languages-signed-c4-and-harmonic-density.md`: the mandatory odd/odd normalization repair, giving Fibonacci chamber word D,D,U rather than all D.
- `05-knowledge/results/session41_owned_attachments_fibonacci_20261023/additive/RESULTS.md`: golden/Gaussian modulo-two obstruction and modulo-three norm isomorphism.

These inherited theorems are not presented as new discoveries or exhaustively re-audited here.

## 1. Exact meaning of one eighth

Use integer squared-radius M>=0 and B_M={(x,y) in Z²:x²+y²<=M}; an arbitrary real radius R has M=floor(R²). Let N(M)=|B_M|, A=floor(sqrt(M)), D=floor(sqrt(M/2)). D4 consists of the signed coordinate permutations. Choose representatives 0<=y<=x.

There are I representatives with 0<y<x, A nonzero axis representatives, D nonzero diagonal representatives, and the origin. Their orbit sizes are respectively 8,4,4,1. Hence

N(M)=8I+4A+4D+1,
K(M)=I+A+D+1=(N(M)+4A+4D+7)/8.

Thus one octant reconstructs the disk, but its lattice-point count is strictly larger than N/8. K/N tends to 1/8 as M grows; the boundary correction is O(sqrt(M)). The origin and boundary stabilizers cannot be discarded from exact counts.

For a positive shell m, let r2(m) count signed ordered representations. Its orbit count is

k(m)=(r2(m)+4[ m is a square ]+4[ m is twice a square ])/8.

The two exceptional conditions are disjoint for positive m. The classical Gaussian factorization formula gives r2(m)=0 if any prime 3 mod4 has odd valuation, and otherwise r2(m)=4 product_(p=1 mod4)(v_p(m)+1). The factorization argument is ordinary Gaussian unique factorization: distribute each split prime's exponent between a Gaussian prime and its conjugate, with four unit choices; inert exponents must be even.

## 2. No extra global Euclidean lattice symmetry

Every origin-fixing Euclidean isometry preserving Z² sends each standard unit vector to an integer unit vector. Its columns are therefore a signed permutation matrix. Conversely all such matrices work. Consequently O(2,Z)=D4, and no larger global Euclidean lattice-isometry group reduces the generic eight-point orbit further.

This is a conditional lower bound on symmetry-only reconstruction, not an information-theoretic lower bound on the entire disk. If the object is known to be the complete disk, M alone plus the lattice rule already reconstructs it. If instead it is an arbitrary D4-invariant subset of B_M, each of the K(M) orbits can independently be included or excluded. There are 2^K(M) such subsets; any lossless fixed-length binary encoding needs at least K(M) bits. A subset's unknown orbit memberships cannot be replaced by the radius.

## 3. Three valid ways to reduce further, with different decoders

### Primitive rays

Every nonzero (x,y) has the unique form d(a,b), where d=gcd(|x|,|y|)>0 and gcd(|a|,|b|)=1. In the representative octant store only primitive directions (a,b), then regenerate d(a,b) for 1<=d<=floor(sqrt(M/(a²+b²))). This reduces stored point representatives, but uses integer dilation, not a norm-preserving symmetry. Origin is separately retained. For arbitrary subsets, every ray needs its actual selected set of d values; full radial closure cannot be assumed.

### Norm shells and rational rotations

For nonzero z,w in Z[i] of the same norm, u=w/z=w*conj(z)/N(z) lies in Q(i), with N(u)=1, and uz=w. Thus rational unit-norm rotations act transitively on each rational circle, and are transitive between its Gaussian points when used as partial maps on the lattice. They do not generally preserve Z[i] globally. Example: (3+4i)/5 sends 5 to 3+4i, but sends 1 outside Z[i].

One shell seed plus the allowed arithmetic rotation instructions can therefore generate the whole shell. Alternatively its integer m and prime factorization generate all representations. This is a stronger arithmetic decoder, not an enlarged O(2,Z) orbit. It is not a free compression of an arbitrary marked shell, since the membership/mark information still needs preserving.

### Algebraic generators

Z[i] is generated as a Z-algebra by i. Together with M this gives another extremely short specification of B_M. The description of the generator and closure operations replaces point data; it is not meaningful to compare its single symbol directly with a fraction of a point list without charging the decoder/model.

## 4. A precise fourth-order angular coordinate

The polynomial invariant ring is

Q[x,y]^D4 = Q[u,v],  u=x²+y², v=x²y².

Proof: sign changes force every monomial to have even x and y exponents, so an invariant is a polynomial in X=x²,Y=y². Coordinate swap then forces a symmetric polynomial in X,Y, which is a polynomial in X+Y and XY. Conversely u,v are invariant.

The pair (u,v) separates real D4 orbits: X,Y are the unordered roots of t²-u t+v. For integer points a candidate pair (u,v) is realizable exactly when these roots are nonnegative integer squares. Equivalently, require a nonnegative square discriminant u²-4v, integral nonnegative half-sums (u±sqrt(u²-4v))/2, and that both half-sums are squares. This is a complete arithmetic realizability test for the orbit coordinate.

Within a fixed norm shell, v alone separates the octant representatives. Norm alone first fails at m=25: (5,0) and (4,3) have equal u=25 but v=0 and144. Every positive shell m<25 has at most one D4 orbit (the checker exhausts these small shells; alternatively inspect 0<=y<=x<=4). The fourth-order angular harmonic is

cos(4 theta)=1-8v/u²,  u>0.

Hence the missing information is an exact quartic angular coordinate. The identity is a genuine relation of norm, angular phase and symmetry quotient. It does not identify this coordinate with the repo's tournament electrical invariants or an arbitrary Gaussian Hessian; such an identification requires a defined map.

## 4b. Separate orbit compression does not preserve addition

There is no addition on D4 point-orbits induced by adding arbitrary representatives. The orbit of 1 also contains -1, yet 1+1=2 and 1+(-1)=0 lie in different orbits. Thus individually quotienting two Gaussian inputs erases their relative alignment. For a pair (z,w), the correct symmetry quotient uses one simultaneous D4 action, not independently chosen actions on z and w. Its finite reconstruction can retain a relative frame, then add in that common frame. This is the exact analogue of retaining one shared chart in an attachment network, at the level of the group-action construction only.

A norm-based additive statistic needs the cross term N(z+w)=N(z)+N(w)+2 Re(z conjugate(w)). Neither pair of separate norms nor pair of separate quartic orbit coordinates supplies that term. The example (1,1) versus (1,-1), viewed as ordered pairs of Gaussian inputs, has identical individual orbit data and different cross terms. Consequently arithmetic uses of compressed color classes must retain joint phases or explicitly sum over them. The quartic separator solves a single-point orbit problem; it does not solve arbitrary additive coupling.

## 5. Why the Berggren and Fibonacci bridges are typed, not universal

The Gaussian squaring map

(n+i m)²=(n²-m²)+2mn i

has norm (n²+m²)² and yields Euclid's Pythagorean triple. Primitive positive triples require gcd(m,n)=1 and opposite parity, followed by conventional leg ordering. If coprime m,n are both odd, the displayed triple has content two and must first be divided by two.

Berggren's tree therefore parameterizes nondegenerate primitive first-quadrant integer-radius Gaussian points, modulo D4, through Euclid parameters; axis points are handled separately. It does not enumerate every Gaussian disk point or every shell. Among primitive rays with both coordinates nonzero, the smallest missing ray is (1,1), whose norm2 is not a square. Scaling a nonsquare norm by d² does not make it a square.

The repo already has a stronger positive Fibonacci bridge than a loose three-color resemblance: the positive solutions |n²-mn-m²|=1 are consecutive Fibonacci parameters. Descent (m,n)->(n-m,m) changes the sign of this form and decreases the maximum until (1,1) or (1,2). For 0<m<n, after primitive normalization their triples lie on the three explicitly described Berggren rays. The descent terminal (1,1) gives the degenerate triple (0,1,1) after normalization and is not a Berggren node. This retains source parameters, parity normalization and ancestry address. It does not supply a map from arbitrary square-sum graph island colors to those branches.

The Session41 ring comparison is also exact and local: modulo2 the golden algebra F2[t]/(t²-t-1) is a field while the Gaussian algebra F2[i]/(i²+1) has a nonzero nilpotent, so they are not isomorphic. Modulo3, t maps to 2+i and gives an isomorphism preserving conjugation and norm. Thus checking both ends yields a prime-dependent bridge, not a universal ring identification.

## 6. Bott periodicity is a different eight

Bott's theorem concerns stable homotopy of classical groups (complex period2, real period8), not the order of the planar lattice isometry group. The Clifford-module version has an eight-degree real algebra stabilization with a matrix factor; it is not an eight-point orbit operation. See the primary papers [Bott (1959)](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/bott4.pdf) and [Atiyah–Bott–Shapiro, Clifford Modules](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/abs.pdf).

The legitimate common design principle is to state what a quotient or stabilization preserves and retain the missing coordinate when reconstruction is desired. Here the omitted coordinate can be v, a radial multiplier, or a shell-membership flag. In native wiring it can be a shared port/chart. This is a methodological analogy, not a transfer theorem between the subjects.
