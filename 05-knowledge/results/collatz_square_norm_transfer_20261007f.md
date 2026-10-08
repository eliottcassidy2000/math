# Square-plus-one norms as marked source-guard coordinates

**PROVED:** the lattice quotient, marked local norm bijections, exact native
guard transport and original-source boundary below. **FINITE-EXACT:** the
declared controls and 756 unchanged receipts from the preceding bounded
inverse bank. **OPEN:** grounding every original source. This package adds
no guard coverage and discovers no ROOT certificate.

Artifacts: [program](../../04-computation/experiments/collatz_square_norm_transfer_20261007f.py)
and [output](collatz_square_norm_transfer_20261007f.out).

## 1. Inheritance and three distinct meanings of the sequence

The sequence `1,2,5,10,17,26,37,50,65,82,101,122` is `a^2+1` for a=0,...,11.
Its successive differences are `2a+1`. It has a precise square-lattice
interpretation, but its terms are not asserted to be optimum packing counts.

The closest existing norm mechanism is
[creation/Fano section5](creation_fano_20260925.md): `n^2+1` decreases under
a positive affine Collatz word exactly when the original source decreases.
The lattice mechanism is
[THM-3336, primitive Gaussian multiplication and content](../../01-canon/theorems/THM-3336-primitive-gaussian-multiplication-content-curved-farey-triangulation.md),
which retains the integer content removed by primitive normalization.
The packing interface in
[Eleven Squares routing, section2](eleven_squares_routing_bridges_20261007.md)
retains the frame, canonical mark, admissible family and clearance predicates;
an unmarked torus does not retain a packing certificate.

The positive new interface here is a **marked finite-local norm coordinate**
that transports an original source's exact binary and ternary guards without
changing their precision. The canonical hostile is positive sources3 and5:
their normalized norms coincide modulo8, although their next valuations are
1 and4. The missing local marker is a square-root branch. This is not global
norm noninjectivity on positive sources: the exact positive norm determines
its positive source by integer square root.

Our board is original source / lattice index / norm value / local branch /
native word / ROOT obligation. These types remain separate:

| Source | Target/map | Preserved predicate | Lost data or required marker |
|---|---|---|---|
| Integer lattice vector | multiplication by a+i | Norm similarity and exact sublattice membership | Primitive reduction requires its content |
| Positive odd integer n | q(n)=(n²+1)/2 | Exact positive source, globally recoverable | No compression claim; q has about twice as many bits |
| Odd local residue | q modulo powers of2 and3 | Native guard on each declared branch | Binary and ternary square-root branches |
| Authenticated inverse-word receipt | same receipt through a norm guard | Same source, child, word and payment | Child ROOT obligation still external |

## 2. The lattice index is exactly a²+1

Multiplication by a+i on Z[i] is the integer matrix

\[
 M_a=\begin{pmatrix}a&-1\\1&a\end{pmatrix},\qquad
 M_a^TM_a=(a^2+1)I,\qquad \det M_a=a^2+1.             \tag{1}
\]

The quotient is cyclic. Explicitly,

\[
 (x,y)\longmapsto x-ay\pmod {a^2+1}                   \tag{2}
\]

is surjective and has kernel exactly `M_a Z²`. Both matrix columns map to
zero. Conversely, if x-ay is a multiple of a²+1, both coordinates of
`M_a^{-1}(x,y)=((ax+y)/N,(-x+ay)/N)` are integers. Hence the square cell of
this sublattice has index and area a²+1. This is a square tiling/lattice
statement, not an optimization theorem about packing unit squares inside
an axis-aligned enclosing square.

For every integer a, `3` does not divide a²+1, so this full two-coordinate
map is invertible modulo every power of3. For even a it is also invertible
modulo powers of2. For odd a, `v₂(a²+1)=1`; its Smith factors are `1,a²+1`,
so the kernel modulo2^k has exactly two elements for every k>=1. This
full-lattice kernel and the scalar square-root branch below are distinct
objects, even though each has a twofold dyadic phenomenon.

Primitive normalization must retain content: the inherited example
`(2+i)(4+3i)=5(1+2i)` removes a factor5. An index or normalized direction
alone does not record this multiplication history.

## 3. Same-precision marked local norm coordinates

For odd integers n set

\[
 q(n)=\frac{n^2+1}{2},\qquad
 q(n)-q(m)=(n-m)\frac{n+m}{2}.                         \tag{3}
\]

**Dyadic theorem.** On each branch `n=b mod4`, b=1 or3, the map q is an
isometry and a bijection onto `q=1 mod4` in Z₂. At every finite k>=2,
it gives a bijection from the branch modulo2^k onto residues1mod4 modulo2^k.
Indeed `(n+m)/2` is odd on the same branch, so
`v₂(q(n)-q(m))=v₂(n-m)`. The domain and target have equal finite size.
The two branch decoders lift one bit at a time; the two preimages of a
given residue are negatives of each other modulo2^k.

**Ternary unit theorem.** On each branch `n=c mod3`, c=1 or2, q is an
isometry and a bijection onto `q=1 mod3` in Z₃. In this setting division
by2 is multiplication by the local unit2⁻¹. The factor `(n+m)/2` is a
3-unit on one branch. Again finite counts prove bijectivity, and the next
ternary digit is uniquely determined. The restriction to units matters:
3 and6 are different modulo9 but their q values agree modulo9. No unit
decoder or unchanged-precision assertion is applied to that ramified stratum.

Combining the two finite bijections by CRT proves the exact guard rule.
Let k>=2, d>=1 and let r be an odd 3-unit residue modulo `M=2^k3^d`.
For a supplied odd integer n,

\[
 n\equiv r\pmod M
 \quad\Longleftrightarrow\quad
 \begin{cases}
 n\equiv r\pmod4,\\
 n\equiv r\pmod3,\\
 q(n)\equiv q(r)\pmod M.
 \end{cases}                                          \tag{4}
\]

Thus a marked norm guard uses the **same prime-power precision** as the
original scalar source guard. The implementation independently lifts each
root and rejoins the residues by CRT. If the two branch labels are forgotten,
each admissible mixed norm residue has four source residues modulo M.
Positive representatives of these distinct local residues need not have the
same exact full norm. In contrast, the full positive integer q uniquely
recovers n through `n=isqrt(2q-1)` with an exact square check.

This precision claim uses the normalized coordinate q, not a raw norm residue
with an unlicensed division by2. To recover q modulo2^k from `N=n²+1`, one
needs N modulo2^(k+1). For example N(1) and N(3) both equal2 modulo8, while
q(1) and q(3) differ modulo8. At the prime3, division by2 has no precision cost.

The smallest displayed dyadic hostile is

\[
 q(3)=5\equiv13=q(5)\pmod8,\qquad
 v_2(3\cdot3+1)=1,\quad v_2(3\cdot5+1)=4.             \tag{5}
\]

For the current Mersenne source `n=2^E-1`, with odd E>=3, the markers are
fixed: `n=3 mod4` and `n=1 mod3`. Its norm coordinate also satisfies

\[
 q(n)=1+2^E(2^{E-1}-1),\qquad v_2(q(n)-1)=E.           \tag{6}
\]

This retains the source exponent as a valuation; it does not advance the
Collatz orbit, lower the exponent or ground the source.

## 4. An actual receipt interface, and its limit

The [bounded inverse cover](collatz_bounded_inverse_cover_20261007e.md)
retains252 labelled words with carrier `(P,Q,B)`, native ternary guard
`n=BQ⁻¹ modP`, and smaller child `h=(Qn-B)/P`. Their even terminal valuation
forces the endpoint n to be1mod3. Split the odd sources into their two mod4
branches, combine each branch with the original ternary residue, and apply
(4). `native_norm_guard` constructs these exact guards. `norm_receipt`
accepts a supplied ordinary source, tests its marked norm guard, then invokes
the existing ordinary receipt verifier. Its result is **the identical receipt**:
same n, same h, same actual word and common endpoint. All252 labels are
retained, including overlapping counting cells.

The test universe uses three ordinary native sources per label, giving756
identical receipts. It neither materializes an astronomical Mersenne nor
changes its residue family. In particular no mass can be added to the
preceding guard cover simply by changing to the norm coordinate.

The companion [CM bridge transfer](collatz_cm_bridge_transfer_20261007f.md)
uses simultaneous local constraints on a marked, source-owned bridge. The
comparison here is structural: both interfaces must preserve all local
markers and the original payment predicate. No CM theorem, Frobenius class,
packing theorem or unmarked norm is used to deduce a Collatz ROOT statement.

**The source-replacement obstruction is elementary.** For every nonzero
Gaussian integer z=a+bi, the odd part of `a²+b²` is1mod4. Divide the two
coordinates by their common power of2. If precisely one remaining coordinate
is odd, the sum is1mod4. If both are odd, the sum is2mod8 and one further
halving again gives1mod4. Therefore every such normalized norm source x>1
already has the known one-step descent

\[
 U(x)\le\frac{3x+1}{4}<x.                             \tag{7}
\]

But the unresolved Mersenne originals are3mod4 and cannot themselves be
odd Gaussian norms. Replacing n by q(n), observing (7), and declaring n
grounded would change the source and supply no common-future certificate.
The inherited exact norm identity reinforces this boundary: for a legal
positive affine word `m=(Pn+B)/Q`,

\[
 (n^2+1)-(m^2+1)
 =\frac{((Q-P)n-B)((Q+P)n+B)}{Q^2}.                   \tag{8}
\]

The second factor is positive. Norm descent is exactly original-source
descent, not a stronger criterion. A full marked coordinate is useful for
guard verification; it is not an independent payment or grounding rule.

## 5. Reproduction and finite scope

Controls include the displayed sequence; lattice parameters0..19 and
coordinates-5..5; dyadic kernel moduli through16; complete dyadic scalar
decoding for k2..9 and ternary-unit decoding for d1..6; complete mixed
residue decoding for k2..6,d1..3;500 exact positive-source norm inverses;
Gaussian coordinates0..30;756 unchanged actual receipts; Mersenne markers
for odd exponents3..29; and malformed-packet/branch/source hostiles.

```text
python -B -X utf8 04-computation/experiments/collatz_square_norm_transfer_20261007f.py
python -B -O -X utf8 04-computation/experiments/collatz_square_norm_transfer_20261007f.py
```

Normal and optimized execution must match the companion output. The local
decoders use exact integer arithmetic and bounded digit lifting. ROOT
certificates remain separately supplied proof objects throughout.
