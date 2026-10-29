# Four ports recover an affine space; the same defect survives equal theta series

**Session97, 2026-10-29. Status:** FINITE-EXACT reconstruction of six marked
codes; PROVED Construction-A transport and equal-theta/nonisometry criterion.
This is a concrete synthesis using classical codes and lattices, not a claim
that isospectral nonisometric lattices are newly discovered.

This note is a new extension of the frozen [four-port checkpoint](REVERSAL_AND_FOUR_PORT_HOSTILE.md).
It imports its eight explicit length32 codes, their two coordinate-isomorphism
classes, and the complete four-port profiles.

## 1. Positive reconstruction, not only an obstruction

For each of the six zero-free codes, take its1240 four-subsets having octad
response7. Mark the existing extension point31 as zero. Define

* a+a=0;
* a+0=0+a=a;
* for distinct nonzero a,b, a+b is the unique c such that `{0,a,b,c}` has response7.

The checker verifies uniqueness for every pair and associativity for **every
ordered triple**. The result is an elementary abelian group of order32, hence
an F2-vector space of dimension5. It constructs a five-bit chart and verifies
the chart intertwines this addition with XOR. Every original octad becomes
a coset of a three-dimensional linear subspace. Since there are

`2^(5-3) * GaussianBinomial(5,3;2) = 4*155 = 620`

such flats, the octads are **exactly all affine three-flats of AG(5,2)**.
The response7 four-subsets are exactly the affine two-flats: there are
`2^(5-2)*GaussianBinomial(5,2;2)=8*155=1240`.

This explains the response table intrinsically. Three distinct points span
an affine plane, lying in seven affine three-flats. Four points either span
a plane (seven containing three-flats) or span a three-flat (one). Thus there
are no zero responses, and the exceptional response7 identifies the affine
addition law. No arbitrary field chart is inserted to force this outcome;
the group law is recovered from the named four-port data first.

The cyclic permutation on the31 nonzero coordinates preserves this addition
and is a Singer cycle. Its minimal polynomial is separately reconstructed
from the group operation. The six codes give exactly the six primitive
irreducible polynomials of degree5 over F2:

| Code generator bitmask | Recovered rotation polynomial bitmask |
|---:|---:|
| 73235 | 41 |
| 83903 | 55 |
| 92135 | 47 |
| 102641 | 37 |
| 118733 | 61 |
| 129989 | 59 |

Bit i again means the coefficient of x^i. The same code class therefore
supports six different marked cyclic field charts. The source's cycle
labelling is precisely the information that ordinary code isomorphism forgets.
The two-code sector has zero-response four-subsets, so this reconstruction
rule cannot define an affine F2^5 geometry there.

## 2. Construction A with the normalization explicit

For a doubly-even self-dual binary code C of length N, define

`L_C={x/sqrt(2) : x in Z^N, x mod2 in C}`.

This standard Construction A is discussed in
[Elkies's error-correcting-code lectures](https://people.math.harvard.edu/~elkies/M256.13/index.html), October17, and
[Rains–Sloane, Self-Dual Codes](https://www.mathe2.uni-bayreuth.de/axel/papers/sloane%3Aself_dual_codes.pdf).
The properties needed here also have short direct proofs:

* Self-orthogonality makes `x dot y /2` integral.
* Since integer squares are congruent to parity bits modulo4, doubly-even
  code weights make every squared lattice norm `sum x_i²/2` even.
* The unscaled integer preimage has index `2^(N-dim C)=2^(N/2)` in Z^N.
  Scaling by1/sqrt2 gives covolume1. Thus the integral lattice is unimodular.

No phase, scaling or determinant convention is implicit in this passage.

### Root-axis criterion

Assume further that C has minimum nonzero weight at least8. A norm-two
lattice vector corresponds to an integer vector with sum of squares4.
The only possibilities are ±2 in one coordinate, or ±1 in four coordinates.
The latter would require a weight-four codeword and is absent. Hence the
roots are exactly

`{±sqrt(2)e_i : 1≤i≤N}`,

a root system A1^N. These roots recover the coordinate axes intrinsically.
Any lattice isometry must permute them, so is a signed coordinate permutation.
Signs disappear modulo2. Consequently

**L_C and L_D are isometric iff C and D are coordinate-equivalent**

for this minimum-weight-at-least-eight class. This proof does not apply to
the Hamming/E8 case, where additional roots are present.

## 3. Equal weight enumerators give equal full theta series

Use a formal variable z that records the squared norm of the unscaled
integer vector, so a coefficient of z^m records lattice norm m/2. Set

`E(z)=sum_(a even) z^(a²)`, `O(z)=sum_(a odd) z^(a²)`.

If `W_C(X,Y)=sum_w A_w X^(N-w)Y^w`, then exactly

`sum_(v in L_C) z^(2||v||²) = W_C(E(z),O(z))`.

This identity follows coordinate by coordinate and is valid as a full formal
series. Thus equal weight enumerators imply equal **entire** theta series,
not just agreement to a tested cutoff.

Apply this to the n31 code generators75035 and92135. Their full weight
enumerators agree, while their four-port profiles prove that they are not
coordinate-equivalent. Both have minimum code weight8, so the root-axis
criterion proves their Construction-A lattices are **nonisometric even
unimodular lattices of rank32 with the same theta series and root system
A1^32**. The eight marked kernels yield exactly two isometry classes here.

The first coefficients, indexed by lattice norm0,2,4,6,8, are

`1, 64, 160704, 64543488, 4845725632`.

The checker independently expands the formal identity through lattice norm16.
Full equality is proved by the identity, not inferred from those coefficients.

Equivalently, the flat tori `R^32/L_C` and `R^32/L_D` are isospectral but
nonisometric. Fourier eigenfunctions are indexed by dual-lattice vectors
with eigenvalues `4*pi²*||v||²`; unimodularity identifies each lattice with its
dual. This is a concrete instance of the classical inverse-spectral
obstruction, not a resolution of an open spectral problem.

## 4. The lost joint ports are intrinsic lattice data

The same obstruction can be stated without choosing a code presentation.
Choose four of the canonically recovered root axes and a root on each.
Count norm-four lattice vectors whose inner product with each chosen root
is odd. Root signs do not affect the predicate. A code octad contributes
256 sign choices; the other norm-four vectors have all-even coordinate
parities and contribute zero. Thus this joint root-port response is

`256 * R4(S)`.

Responses involving up to three root axes are universally
`256*155`, `256*35`, `256*7`, respectively. But one rank32 lattice has11160
four-axis subsets with zero response and the other has none.

So even **the entire theta series plus every one-, two-, and three-root-axis
response** misses a distinction seen at four axes. This is exactly the
previous joint-data loss transported to an isometry-invariant lattice
construction. One class bit decides which lattice isometry class occurs;
to answer arbitrary marked future queries, the affine/code chart still
matters. The result is a carrier theorem with a precise observer family,
not a general assertion that spectral data never suffices.

## 5. Hamming, Golay and the Leech boundary

The same normalization gives these calibrations:

| Binary extended code | Rank | Number of norm-two roots |
|---|---:|---:|
| Hamming `[8,4,4]` | 8 | 240 |
| Golay `[24,12,8]` | 24 | 48 |
| either rank32 code class | 32 | 64 |

The first is the standard E8 construction. The second has root system A1^24:
it is the corresponding Niemeier lattice, **not the Leech lattice**, which
has no roots. Reaching Leech requires an additional neighbour/glue operation;
no such operation has been constructed in this packet.

Prior repo work
`05-knowledge/results/leech_deep_hole_swap_theorem_20261022.md` and
`THM-4758-leech-deep-hole-swap-theorem-and-the-deep-hole-theta-identity.md`
uses the A1^24 Niemeier lattice, marked cusp/glue data and shifted theta
series. Those older finite certificates have **not** been replayed here and
are not dependencies of our rank32 argument. The exact common object is
A1^24 in the Golay Construction-A calibration, not an identification of its
48 roots with Leech vectors of norm two.

A second recovered niche source,
`gold_budget_golden_golay_20261019.md`, §3, already describes reciprocal
**ternary** Golay syndrome modules over F243 and an order11 coordinate action.
Its full budget multiplier swaps the modules; its square gives the finite
Collatz-labelled action. Today's binary reciprocal-choice/reversal mechanism
is the same representation-theoretic pattern with a different alphabet and
carrier. The old conjugacy is not a global Collatz proof, and branch-point
permutations determine Jacobian2-torsion, not automatically3-torsion.

## Reproduction and scope

`check_affine_and_lattice.py` makes459670 explicit checks beyond its imported
cyclic certificate. It reconstructs all six affine groups, checks every
associativity triple and octad closure, recovers all six rotation polynomials,
and computes four theta calibrations. Normal and optimized output agree.
The root/isometry and full-theta arguments are proofs above; low-shell
computations are checks, not substitutes for those arguments.
