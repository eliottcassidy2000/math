# Session130: periodic phase carriers, arithmetic norms, and the n=5 lines

**Status: PROVED scoped carrier/codec theorems, with independently recomputed FINITE-EXACT n=5 quotient witnesses.** No claim that PRIMEGAME generates the merged metagraph, and no literature-priority claim. The complete period classification is inherited, not reproved here.

## Results first

1. Raw PRIMEGAME periodic cycles cannot model the n=5 tiling cube equivariantly, even by taking any finite disjoint union of them. With complement modeled by a half-period and reflection by time reversal, there are too few reflection-fixed line orbits: n=5 has 8 blue lines and 12 pairs of black lines; each six-cycle supplies one of each, while each two-cycle supplies only a blue line. This is the first obstruction, absent at n=3 and n=4.
2. The correct repair is **binary decoration of a rooted six-cycle**. Its 64 decorations carry exactly the tiling complement/reflection action. There are precisely 16 reflection-compatible coordinate gauges. An explicit one is given below.
3. The six actual periodic states, not six arbitrary external primes, give a primitive, lossless arithmetic codec for all 64 decorations. The codec works uniformly for every admissible rescaling of the cycle. A line is blue exactly when its arithmetic reflection norm is a rational square.
4. This lift preserves literal complement lines and their colors, but not autonomous merged-node evolution. Temporal phase shift changes 6 of the 8 blue lines to black, and no coordinate gauge makes phase shift descend to the ten merged nodes.

## 1. Inherited objects and exact conventions

The inherited Session129 classification says every periodic orbit has least period 2 or 6. The six-cycle, rooted at its unique state immediately before instruction G, is

```
b = (4301,17765,16302,9867,8671,23023)
word = GBDEFA.
```

Every six-cycle is uniquely `b_j X`, with `gcd(X,51765)=1` and `38` not dividing X. Equivalently

`X=S*2^a*11^r*13^v*19^w*23^s`,

where the displayed exponents are nonnegative, `a*w=0`, and S has no prime factor in `{2,3,5,7,11,13,17,19,23,29}`. The gcd of its six states is X. The two-cycles are `11M <-> 13M` with the inherited admissibility condition on M.

For the n=5 fixed path `5 -> 4 -> 3 -> 2 -> 1`, free tiles in the inherited gauge are

`T=((5,1),(4,1),(3,1),(5,2),(4,2),(5,3))`.

A bit 0 means the larger label points to the smaller; bit 1 reverses that free arc. The four path arcs remain fixed. Complement C toggles all six free bits; R permutes the tile coordinates as

`r=(0)(1 3)(2 5)(4)`.

Literal lines are C-orbits, not merged nodes. A blue line has R-fixed endpoints. A black line belongs to a distinct pair exchanged by R. No tiling satisfies `CR(t)=t`, since R has fixed coordinates. There are 16 R-fixed tilings, hence 8 blue lines, and the remaining 48 tilings form 12 free Klein-four orbits, hence 24 black lines.

## 2. Raw-periodic-state obstruction and its sharp onset

**Theorem 1 (raw cycle carrier obstruction).** On a cycle of length m in `{2,6}`, give phase j the involutions `C(j)=j+m/2` and `R(j)=c-j` modulo m. No finite disjoint union of such cycles is equivariantly bijective to the n=5 tiling cube under C and R.

**Proof.** A valid source cannot contain a point fixed by CR, because the target has none. If m=6 and c is odd, the congruence `2j=c+3 (mod 6)` has two solutions, so that reflection choice is excluded. If m=6 and c is even, R has two fixed points forming one C-line, while the remaining four points form a free Klein-four orbit: one blue line and one black-line pair. If m=2, the only allowable R is the identity, contributing one blue line and no black pair. Consequently a union of a two-cycles and b six-cycles has `(blue lines,black pairs)=(a+b,b)`. Equating to `(8,12)` would require `a=-4`. This is impossible. QED.

At n=3, `(blue lines,black pairs)=(1,0)`, and one two-cycle realizes the action. At n=4 the pair is `(2,1)`, realized abstractly by one two-cycle and one six-cycle. These are statements about C/R-sets only, not their merged-node incidence.

More generally, if a coordinate reflection has d fixed positions and s transposed pairs, the binary cube has

`B=2^(d+s-1)` blue lines and `F=(2^(d+2s)-2^(d+s))/4` free Klein-four orbits,

provided d>=1. The raw 2/6-cycle model requires `B>=F`, equivalently `2^s<=3`. The n=5 reflection has s=2, the first failing order. This identifies the mechanism rather than a mere mismatch of total cardinalities.

## 3. Binary decoration repairs the carrier

Let phases be `Z/6Z`, with rooted reflection `iota(j)=-j`. A decoration is a function `x:Z/6Z->{0,1}`. Define

`(C x)_j=1-x_j`, `(R x)_j=x_{-j}`.

The phase root is intrinsic to the six-cycle: the unique G-state is phase 0, and increasing phase follows the program. Thus this phase reflection is specified without choosing among six possible starting points. The coordinate identification with tiles remains a declared gauge.

**Theorem 2 (explicit equivariant lift).** Set

`phi=(0,1,2,4,5,3)`,

meaning phase j is assigned tile T[phi[j]]. The map `E_tile(x)[phi[j]]=x_j` is a C/R-equivariant bijection from all binary decorations to all n=5 tilings. It induces a bijection of 32 literal lines preserving blue and black.

**Proof.** The assignment is bijective and satisfies `phi(-j)=r(phi(j))`; bitwise complementation commutes with every coordinate permutation. QED.

Exactly 16 coordinate bijections intertwine these reflections: choose the matching of their two fixed coordinates in `2!` ways, the matching of their two transposition orbits in `2!` ways, and the two orientations in each paired orbit in `2^2` ways. No particular one of these is forced by the FRACTRAN fractions.

This is not the unsupported statement “six steps equal six free bits.” The source has been enlarged from six states to the Boolean function space on six states, and the preserved operation is explicitly proved.

## 4. A primitive arithmetic codec from the periodic states

Factor the canonical states:

|j|b_j|factorization|
|---:|---:|---|
|0|4301|11*17*23|
|1|17765|5*11*17*19|
|2|16302|2*3*11*13*19|
|3|9867|3*11*13*23|
|4|8671|13*23*29|
|5|23023|7*11*13*23|

For a positive rational z, define the integer valuation decoder

`D(z)=(nu17(z)-nu5(z), nu5(z), nu2(z), nu3(z)-nu2(z), nu29(z), nu7(z))`.

**Theorem 3 (canonical primitive codec).** For every integer vector e in Z^6,

`D(product_j b_j^e_j)=e`.

Therefore the b_j are multiplicatively independent and their valuation vectors span a primitive direct summand of the valuation lattice. In particular

`Z_1(x)=product_j b_j^x_j`

is a lossless integer encoding of all 64 binary decorations.

**Proof.** Direct factorization gives `D(b_j)=unit_j`. Thus D is an integral left inverse of the six-column valuation matrix B. The endomorphism BD is an integral projection onto its image; consequently the image is a direct summand. QED.

**Cycle-relation interpretation.** Let `d_j=nu(b_(j+1)/b_j)` be the six valuation increments in GBDEFA order. The same decoder gives `D(d_j)=e_(j+1)-e_j`, with cyclic indices. Their only integer relation is therefore `sum_j d_j=0`; their span is the rank-five primitive sublattice corresponding to the coordinate-sum-zero lattice in Z^6. Indeed the latter is primitive, and the state lattice is an integral direct summand by the left inverse. This is the Session129 formal relation ABDEFG written as the boundary kernel of the actual cycle. Adding any degree-one phase state integrates the rank-five increment lattice to the rank-six state codec. It explains why the return relation and the six-coordinate codec fit together.

The useful degree character is

`kappa(z)=sum_j D_j(z)=nu3(z)+nu7(z)+nu17(z)+nu29(z)`.

Every b_j has degree 1. Moreover every one of the six executed instructions A,B,D,E,F,G has degree increment 0. These registers provide a circulating marker: `17 -> 3 -> 29 -> 7 -> 17`, with G and D retaining the marker. This statement is about the six-cycle instruction subsystem; kappa is not conserved by every instruction in the full program.

**Theorem 4 (all allowed scales; integral feedback absorption).** Fix any admissible X and let `a=nu2(X)`. Write `c_j=X*b_j`. The codec `Z_X(e)=product_j c_j^e_j` has the integral decoder

```
x0 = nu17(z)-nu5(z)
x1 = nu5(z)
x2 = nu2(z)-a*kappa(z)
x3 = nu3(z)-nu2(z)+a*kappa(z)
x4 = nu29(z)
x5 = nu7(z).
```

It is primitive and injective on all Z^6, not only on the 64 binary vectors. The scale is recoverable from the complete physical cycle as the gcd of the six c_j; it is an input to this scaled decoder, not something inferred from an arbitrary single codeword.

**Proof.** Allowed X has no factors 3,5,7,17,29, so `D(X)=w=(0,0,a,-a,0,0)` and `1^T w=0`. Applying D to `Z_X(e)` gives `(I+w 1^T)e`. Since `(w 1^T)^2=0`, the integral inverse is exactly `I-w 1^T`, yielding the formulas. This again gives an integral left inverse of the scaled valuation matrix. Finally the base states have gcd 1, so the scaled cycle has gcd X. QED.

This rank-one nilpotent absorption is the precise reusable algebra; no smallness or real convergence condition is involved.

## 5. Blue lines are exactly square reflection norms

Fix an allowed X. Let H_X be the free multiplicative subgroup generated by c_0,...,c_5 and let

`rho_X(c_j)=c_{-j}`

extend multiplicatively. This is well-defined by Theorem 4. Let `K_X=product_j c_j`; for X=1,

`K_1=2453532251388619737666330`.

On the binary code set `A_X={Z_X(x):x in {0,1}^6}`, define

`C_X(z)=K_X/z`, `N_X(z)=z*rho_X(z)`.

**Theorem 5 (arithmetic color criterion).** Under the explicit phase/tile gauge, arithmetic C_X models tiling complement, rho_X models R, and a literal line `{z,K_X/z}` is blue if and only if any of the following equivalent conditions holds:

1. `rho_X(z)=z`;
2. `x1=x5` and `x2=x4`;
3. `nu5(z)=nu7(z)` and `nu2(z)-a*kappa(z)=nu29(z)`;
4. `N_X(z)` is a square in the positive rationals (equivalently an integer square).

**Proof.** The first three conditions follow from the primitive decoder and the reflection permutation. The norm has basis exponents

`(2x0, x1+x5, x2+x4, 2x3, x2+x4, x1+x5)`.

If the bits match, it is z^2. Conversely if the rational norm is a square, every prime valuation is even. Applying the integral decoder implies every displayed basis exponent is even, hence each binary pair has equal bits. Since the norm is a positive integer, a rational square root is integral. QED.

The square test is independent of the selected one of 16 reflection-compatible tile gauges, although the projected endpoint node labels generally depend on that gauge. It is also uniform across allowed scales. The norm is not a norm of a quadratic number field: it is the multiplicative norm of an explicitly defined involution on H_X.

**Sharp binary boundary.** On the full integer-exponent lattice, a square norm says only `x1=x5 (mod 2)` and `x2=x4 (mod 2)`. It need not imply reflection-fixedness: `z=c1^2` has norm `(c1*c5)^2` but `rho_X(z)=c5^2 != z`. The square iff blue theorem must not be extended beyond the binary codewords without restoring the actual defect vector.

The defect vector `(x1-x5,x2-x4)` supplies two exact reflection-odd coordinates. Blue is its zero locus; R and C each negate it. The norm alone is not a complete orbit invariant: decorations with supports `{1,2}` and `{1,4}` have the same norm, but their C/R-orbits are distinct. The missing information is the relative alignment of the two unequal pairs. In the selected tile gauge these are masks 6 and 34, belonging to distinct merged nodes 2 and 5. Thus square/non-square detects color perfectly without determining the entire line or node.

## 6. The temporal operation does not survive the node quotient

Let tau shift phase indices by one: `(tau x)_{j+1}=x_j`. It commutes with C and satisfies `R tau R=tau^{-1}`.

**Theorem 6 (blue under a phase shift).** Exactly 4 of the 16 R-fixed decorations remain R-fixed after tau. They are the decorations constant on the even phases and constant on the odd phases. Consequently only 2 of the 8 blue lines remain blue; the other 6 become black.

**Proof.** If Rx=x, then `R tau x=tau^{-1}x`. Therefore tau x is also R-fixed precisely when `tau^2 x=x`. The latter says even phases have one common bit and odd phases have another. There are four choices. Complement pairs them into two lines. QED.

In particular, a different temporal root changes which mirror is used; it does not preserve the fixed tiling mirror. The canonical G-root avoids ambiguity, but it does not promote tau to a color-preserving operation.

**Theorem 7 (no merged-node dynamics from phase shift).** For every coordinate assignment from the six phases to the six free tiles, including all 16 reflection-compatible gauges, the induced coordinate six-cycle does not descend to a map on the ten merged n=5 nodes.

**Proof.** If a permutation of a finite set descends through a surjective quotient q to a function f, then f is a permutation. The source permutation maps each entire q-fibre bijectively to the corresponding f-fibre; hence fibre cardinality is preserved. In the fixed-path cube, the mask 16 is the sole tiling of the merged node `1 -> C3 -> 1`: all Hamiltonian paths begin at its source, end at its sink, and the middle C3 has three rotationally equivalent paths. A coordinate six-cycle moves mask 16 to another one-bit mask. Every other one-bit mask has a merged fibre of size at least two: masks 1 and 7 share a node, 2 and 8 share a node, and 4 and 32 share a node. The exhaustive 64-mask checker independently confirms these witnesses and the singleton fibre `{16}`. Thus the singleton fibre would have to map to a larger fibre, contradiction. QED.

In the selected gauge, the smallest equal-node pair with unequal shifted nodes is `(1,28)`: both lie in node 4, but tau sends them to masks `(2,49)`, in distinct nodes `(2,8)`. Node labels are exactly the inherited sorted canonical-adjacency labels. The checker also records a witness for each of the 16 compatible gauges.

A phase decoration can be fixed by zero, one, or all three vertex-axis temporal mirrors, but never exactly two. The 64 decorations split into counts 24,36,4 respectively (hence the 32 complement lines split into 12,18,2). To see the impossibility of exactly two, the product of two distinct mirrors is a shift by two or four, so invariance under two forces parity-class constancy and invariance under the third. The four common fixed decorations are exactly those from Theorem 6; subtracting their 12 incidences from the three mirrors' 48 total fixed incidences leaves 36 uniquely fixed decorations.

The rooted state reflection itself is not a global prime-renaming symmetry: it exchanges b_2, with five distinct prime factors, and b_4, with three. A permutation of primes preserves that count. Nor is it a physical forward FRACTRAN step. It is a derived operation on a complete rooted periodic orbit; rho_X is its multiplicative extension to the codec.

## 7. Connection contract and scope

|Field|Precise content|
|---|---|
|Source|A complete rooted admissible six-cycle, plus one binary decoration on its phases|
|Target|A fixed-path n=5 tiling, then its literal C-line|
|Map|Product of selected phase states; integral decoder; explicit phase-to-tile permutation phi|
|Preserved|All six decoration bits, C and R actions, 32 literal lines, exact blue/black color; all admitted scales|
|Lost on passing to a merged node|Marked-path presentation, phase-coordinate information and autonomous phase-shift evolution|
|Necessary sidecar|Fixed scale X and phase/tile gauge; the full marked tiling if operations are to be continued|
|Cheapest hostile|Raw-cycle count 8<12; shifted blue mask; singleton-fibre obstruction at mask16|
|What is not claimed|PRIMEGAME generates n=5 nodes, the six physical states are 64 tilings, an unchosen canonical tile gauge, or a prime-renaming/global program symmetry|

The earlier Session129 three-coordinate node codec and this six-coordinate literal-tiling codec have different carriers and different preservation guarantees. They are compatible as dictionaries after projection, but no inverse from the node codec to the literal tiling exists without selecting among its fibre.

## 8. Evidence and reproduction

`python3 check_primegame_lines.py > normal.json`

`python3 -O check_primegame_lines.py > optimized.json`

The checker imports no inherited producers or external packages. It verifies every one of 64 decorations and every one of 32 lines, reconstructs all 10 merged nodes from labelled adjacency matrices, checks all 16 compatible gauges, and tests 270 admissible scales with exponents 0..2 at 2,11,13,19,23 and S in {1,31}. It also exercises the integral decoder on signed exponents and distinguishes the norm-collision hostile. Failure gates remain active under `-O`. All-universe claims are justified by the proofs above; finite tests are not substituted for those proofs.

Source pins are retained separately. Inheritance: Session129 `partition/README.md` §4 (complete cycle classification), `periodic_density/README.md` §1 (coprimality form), `tournaments/REPORT.md` §§1–3 (tiling conventions and finite node census), and THM-796, `01-canon/theorems/THM-796-three-sorted-recursive-tiling-line-node-incidence.md` (three-sorted carrier discipline). No external literature was consulted or claimed.
