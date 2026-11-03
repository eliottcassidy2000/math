# Fano sign gauges form a Hamming coset, with an explicit E8 lift

**Status:** PROVED incidence-code identification, one-line correction,
parity extension, lattice properties and explicit isometry; FINITE-EXACT
Hamiltonian counts for all 128 Fano line orientations. The E8 identification
is classical; this note supplies the exact map from the preceding arithmetic
and octonion sidecar. No Collatz, Bott-periodicity, or integer-ring implication
is asserted.

## Inheritance and concept board

The preceding [squarefree-symmetry note](arithmetic_braids2_20260917_squarefree_symmetry.md)
constructs the seven nonzero squareclasses supported on three primes, the
labelled Fano plane, and an octonion multiplication-sign tournament. Its
128 imaginary-basis sign choices produce 16 distinct labelled tournaments.
This note identifies the linear map behind that quotient.

Closest inherited mechanism: the character kernel of the sign-gauge map.
Canonical hostile: 128 independently oriented Fano lines produce 112
tournaments outside that gauge family. Corrected near miss: seven vertices,
21 pairs, or a repeated count of 240 does not identify mathematical objects.
Least-used sidecar: the normal vector labelling each Fano line; it becomes
an error syndrome.

| Live concept | Exact operation | Preserved information / loss |
|---|---|---|
| Seven prime squareclasses | label `p,q,r` by a binary basis | XOR multiplication; exponent heights lost |
| Fano line signs | vertex gauges act by incidence | labelled multiplication signs; eight gauges have identical action |
| Hamming code | compute a three-bit syndrome | unique correction of one line; two errors may miscorrect |
| Extended affine code | append the parity coordinate | all affine evaluations on `F2^3`; chosen zero coordinate added |
| Construction A | lift parity classes to integer vectors | addition and metric; octonion multiplication is not transported |
| Existing E8 routes | compare actual maps | common target need not give a map between source objects |

The repo already contains the same extended code through a different route:
[THM-480, skew-tower row code](../../01-canon/theorems/THM-480-skew-tower-row-code-is-dplus.md)
identifies its order-eight code as `RM(1,3)`. Its higher-dimensional claims
are not dependencies here. [THM-868, E8 score lattice](../../01-canon/theorems/THM-868-e8-bridge-score-lattice.md)
uses tournament scores and a different map. We use neither its interpretive
conclusions about repeated eights nor a hypothetical identification of the
two source constructions. The general warning about lost coordinates is
essential in both comparisons.

## 1. The incidence map is exactly a Hamming code

Let `V=F2^3`, label the seven points by `x in V\{0}`, and label the seven
lines by their nonzero normals:

```text
L_a = {x != 0 : a dot x = 0},  a != 0.
```

Each line contains three points. A vertex sign choice is
`sigma_x=(-1)^{g(x)}`, with `g in F2^7`. On the line `{x,y,x+y}` the
multiplication sign changes by `sigma_x sigma_y sigma_(x+y)`. Its line-flip
vector is therefore the incidence image

```text
(Bg)(a) = sum_(x in L_a) g(x).
```

All sums in this section are in `F2`. Since the line indicator is
`1+a dot x`, set `S=sum_x g(x)` and `v=sum_x g(x)x` to obtain

```text
(Bg)(a) = S + a dot v.                                      (FC1)
```

The vectors `(1,x)`, `x!=0`, span `F2^4`: differences yield every `(0,v)`,
and one column supplies the first coordinate. Thus `B` has rank four and
its image is precisely the evaluation code

```text
C = {(S+a dot v)_(a!=0) : S in F2, v in F2^3}.
```

This is punctured `RM(1,3)`. Its weights are

```text
W_C(z) = 1 + 7 z^3 + 7 z^4 + z^7.                          (FC2)
```

Indeed a nonzero linear functional has four ones on the seven nonzero
points; its complement has three. The remaining words are zero and all
ones. Equivalently,

```text
C = {c : sum_(a!=0) c(a)a = 0}.                             (FC3)
```

To check this equality, every affine evaluation satisfies the displayed
three parity checks: the coordinate sums and pairwise coordinate products
over `F2^3` are even. The parity-check matrix has all seven nonzero columns,
so its kernel has dimension four. Hence `C` is exactly the binary
`[7,4,3]` Hamming code, with this explicitly specified coordinate order.

The kernel of `B` has dimension three. Every character exponent
`g(x)=u dot x` lies in the kernel, since its restriction to a line has zero
or two ones. These eight words exhaust the kernel, giving the
`[7,3,4]` simplex code with enumerator `1+7z^4`. This recovers the preceding
`128/8=16` gauge count as rank-nullity, not numerical coincidence.

The point and line coordinate sets are different. Identifying them in the
formulas uses the chosen dot-product polarity; this sidecar must be retained.

## 2. The octonion locus is an affine coset, with a unique one-line repair

Use the precise Cayley-Dickson convention from the companion:

```text
(a,b)(c,d) = (ac-db*, a*d+cb),   (a,b)*=(a*,-b).
```

For basis indices `x,y=1,...,7`, distinct, write
`e_x e_y=alpha(x,y)e_(x XOR y)` and orient `x -> y` exactly when
`alpha(x,y)=+1`. This convention is given by
[Baez, The Octonions, Cayley-Dickson construction](https://math.ucr.edu/home/baez/octonions/node5.html).
It preserves the XOR output as a sidecar and orients each Fano line cyclically.

Index line coordinates by normals `a=1,...,7`, and write its vertices in
increasing order `(x,y,z)`. Orientation bit one means the cycle
`x -> y -> z -> x`; zero means its converse. For this convention the base
orientation mask is `m0=96`, bits corresponding to normals six and seven.
The full gauge locus is

```text
G = m0 + C.                                                 (FC4)
```

Here addition is XOR. In particular, this is an affine Hamming coset,
not the zero-containing code in the stated orientation convention:
`syndrome(m0)=6 XOR 7=1`.

For any orientation `m`, define its relative syndrome

```text
s(m) = sum_(a!=0) (m(a)+m0(a)) a.
```

If the syndrome is zero, `m` is a gauge orientation. If it is nonzero,
reversing **the entire Fano line with normal `s(m)`** changes that syndrome
to zero. This reverses three tournament arcs, not one arc. The repair is
unique because the seven line columns are distinct and nonzero. Also,
the code's minimum distance three prevents two gauge orientations from
both being within distance one of the same word. Thus the 16 radius-one
balls, each containing eight orientations, partition all 128 possibilities.

The independent exhaustive census agrees with the earlier companion:

| Relative syndrome | Number of orientations | Directed Hamiltonian paths H |
|---|---:|---:|
| zero | 16 | 189 for each |
| each of the seven nonzero values | 16 | 171 for each |

These Hamiltonian values are **FINITE-EXACT**, checked both by subset dynamic
programming and by all `7!` vertex permutations for every orientation.
The algebraic code proof does not itself calculate `H`. Together, the two
results show that the unique one-line correction changes `H` from 171 to
189 in this explicitly bounded family.

**Failure boundary.** Starting from a gauge word, flipping lines with normals
one and two has syndrome three. The decoder then flips line three and lands
at a different gauge word, since the three-error mask is a weight-three
codeword. This is one-error correction, not recovery from arbitrary changes.
Independent line orientation is also more general than changing the signs
of seven basis units: a single basis sign changes a pencil of three lines.

## 3. Adding the missing point gives the extended Hamming code

Append coordinate `a=0` with the parity of the seven-bit word. For
`c(a)=S+a dot v`, this parity is `S`, since a nonzero linear functional
has weight four. Consequently the extension is exactly

```text
C_hat = {(S+a dot v)_(a in F2^3)} = RM(1,3),
W_(C_hat)(z) = 1 + 14 z^4 + z^8.                            (FC5)
```

The fourteen weight-four supports are the affine hyperplanes of `F2^3`.
Each word has weight divisible by four. For any two words `c,d`,
`wt(c+d)=wt(c)+wt(d)-2|supp(c) intersect supp(d)|` shows that their support
intersection is even. Thus the code is self-orthogonal, and its dimension
four in length eight makes it self-dual. It is the extended `[8,4,4]`
Hamming code. The extra coordinate is the zero vector in the **dual normal
coordinate space**; regarding it as the scalar octonion unit requires the
chosen identification of primal and dual binary spaces.

This step extends the linear difference code `C`. Extending the translated
orientation set `G` without first choosing a base produces an affine set,
which must not silently be treated as a linear code.

## 4. Construction A yields E8 with an explicit orthogonal map

Define

```text
M = {z in Z^8 : z mod2 in C_hat},
Lambda = (1/sqrt(2)) M.                                    (FC6)
```

The preimage `M` has index `2^(8-4)=16` in `Z^8`. Scaling all eight
coordinates by `1/sqrt(2)` multiplies covolume by `2^(-4)`, so `Lambda`
has covolume one. For `z in M`,

```text
sum z_i^2 = wt(z mod2) mod4 = 0 mod4.
```

Thus squared norms in `Lambda` are even integers. Polarization gives
integer inner products, and integrality together with covolume one makes
the lattice self-dual. This proves that `Lambda` is even and unimodular,
including the scaling factor that is easy to lose.

We can identify it with standard E8 directly, without appealing only to
the uniqueness theorem for even unimodular lattices in dimension eight.
Use coordinate order `a=0,...,7`. Let `Q` be the orthogonal matrix applying
`(u,v) -> ((u+v)/sqrt(2),(u-v)/sqrt(2))` to the pairs `(0,1)`, `(2,3)`,
`(4,5)`, `(6,7)`. Then

```text
Q(z/sqrt(2)) = ((z0+z1)/2,(z0-z1)/2,...,(z6+z7)/2,(z6-z7)/2).
```

The parity difference within each pair is the same binary coefficient
of the affine word. Therefore all eight output coordinates are integers,
or all eight are half-integers. Their sum is `z0+z2+z4+z6`, which is even:
the sum of an affine function over that four-point plane is even. Hence

```text
Q(Lambda) subset D8 union (D8+(1/2,...,1/2)),
D8 = {w in Z^8 : sum w_i is even}.
```

The right side has covolume one: `D8` has covolume two and its union with
the displayed coset is an index-two extension. Since `Q(Lambda)` also has
covolume one, the inclusion is equality. This is the standard E8 model,
as specified in [Kurkoski, The E8 Lattice and Error Correction in Multi-Level
Flash Memory, section II.C, equation (10)](https://arxiv.org/pdf/1009.5764).
The extended-Hamming construction itself is also recorded in
[Wilson, The Leech lattice, section 2](https://webspace.maths.qmul.ac.uk/r.a.wilson/talks_files/Leech.pdf).

**PROVED complete minimal-vector count.** A nonzero squared norm is at
least two by evenness. Norm two means `sum z_i^2=4`, and there are exactly
two possible integer shapes:

* One coordinate is `+2` or `-2`: `8*2=16` vectors.
* Four coordinates are `+1` or `-1`, on one of the fourteen weight-four
  code supports: `14*2^4=224` vectors.

Both types occur, giving `16+224=240` roots. They generate the lattice:
the coordinate roots generate `2Z^8`, and the weight-four words generate
the binary code, so their lifts together generate all of `M`.
The explicit map sends the root set onto the usual 112 integer roots
and 128 half-integer roots of E8. The script checks this set equality and
closure under all `240^2` root reflections with integer arithmetic.
This conclusion concerns the additive lattice and its Euclidean metric.
Indeed, if coordinate zero is identified with scalar octonion `1`, this
normalized embedding is not closed under the inherited Cayley-Dickson
product: it contains `sqrt(2)*1`, whose square `2*1` is absent because
`2*sqrt(2)` is not an integer coordinate of `M`. An octonion-order claim
would need a different embedding or multiplication normalization.

The number 240 here counts lattice roots. The earlier 240 counts labelled
members of an `S7` tournament orbit. No natural equivariant bijection
between these sets has been constructed; equal cardinalities alone supply
none. The code and the displayed isometry, not that equality, are the bridge.

## 5. Why this construction singles out three bits and eight coordinates

The incidence calculation extends to `F2^r`, `r>=2`: take all nonzero
points and all nonzero hyperplane normals. Equation (FC1) still holds,
so the image is punctured `RM(1,r)`, with length `2^r-1`, dimension `r+1`,
and weights zero, `2^(r-1)-1`, `2^(r-1)`, `2^r-1`. Parity extension gives
`RM(1,r)` in length `2^r`.

For this extended code to be self-dual its dimension must satisfy

```text
2(r+1)=2^r.
```

Among integers `r>=2`, this happens exactly at `r=3`: it fails at two,
holds at three, and the exponential exceeds the linear expression from
four onward. At three self-orthogonality was proved above. At four the
extended code is `[16,5,8]`, not a self-dual `[16,8]` code, and its
Construction-A covolume is eight. Thus the appearance of eight has an
exact rank-and-duality mechanism in this construction. It does not prove
Bott periodicity or make all appearances of eight equivalent.

The arithmetic entry point is equally specific. For `N=p^2qr`, the seven
proper nontrivial squarefree divisors are precisely the nonzero classes
of the three-prime squareclass space. Choosing `p,q,r` as a basis provides
the labels used above. The passage to squareclasses destroys exponent
heights and integer addition; the subsequent octonion signs add data.
This assignment of prime generators to distinct imaginary basis units
cannot preserve integer multiplication: `pq=qp`, whereas
`e_p e_q=-e_q e_p`. This route supplies no Collatz descent theorem.

## Reproduction and hostile controls

Run from the repository root:

```text
python 04-computation/experiments/arithmetic_braids2_20260917_fano_code.py
python -O 04-computation/experiments/arithmetic_braids2_20260917_fano_code.py --output C:/tmp/arithmetic_braids2_fano_code_optimized.json
```

The [script](../../04-computation/experiments/arithmetic_braids2_20260917_fano_code.py)
uses only the Python standard library, explicit `require` checks that remain
active under `-O`, and deterministic LF-terminated JSON. The
[output](../../04-computation/experiments/arithmetic_braids2_20260917_fano_code.json)
stores its source SHA-256, all codewords, the eight syndrome cosets, an
integer lattice basis, and its exact Gram matrix with determinant one.

The finite universes are all 128 vertex gauges, all 128 Fano orientations,
all 256 length-eight parity words, all integer vectors in `[-2,2]^8` of
squared length four, and all pairs of the 240 resulting roots. Independent
paths are the incidence formula versus explicit matrix action; affine
evaluation versus parity-check kernel; Hamiltonian DP versus all vertex
permutations; and root shape generation versus the full bounded integer box.
Hostile controls record the nonzero affine base syndrome, the two-line
miscorrection, and the weight-three words that invalidate a zero-bit
extension. Omitting the `1/sqrt(2)` scale instead leaves covolume sixteen
and minimum squared norm four; it cannot be called this normalized E8 lattice.
