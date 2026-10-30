# Independent operational and palette audit

2026-10-30. **PASS within the declared scope.** No mathematical correction or blocker found in `tournaments/PALETTE_SYMMETRY.md` or `operations/OPERATIONS_AND_OBSERVERS.md`. This audit does not assess literature novelty, distribution of primes, or convergence of aliquot sequences. The cited global classification of square-sigma prime powers is outside this lane; the particular values sigma(3)=4, sigma(81)=121, sigma(343)=400 suffice for its six-state example.

## Independent proof review

* **Joint Burnside.** An even vertex cycle has an antipodal edge reversed by a power, forbidding a fixed tournament regardless of colors. Odd vertex cycles of lengths a,b contribute gcd(a,b) edge orbits between them and (a−1)/2 within one cycle. A coloring on an l-cycle is determined by its first value, with consistency g^l(c)=c. Thus the product of Fix(g^l), not a power of Fix(g) or a simple group-order division, is necessary and sufficient.
* **Weyl quotient and gaps.** Coinvariants in characteristic zero preserve the equivariant attachment/deletion identity. Division by q2^n gives DE−ED=I. Its graded harmonic decomposition yields eigenvalue (n−j)/n on the jth harmonic, of dimension a_j−a_(j−1). GL has two color orbits and therefore a surviving degree-one mode; AGL has one color orbit but two ordered-pair orbits and therefore first surviving degree two. The gaps are 1/n for GL (n≥1), 2/n for AGL (n≥2). The AGL singleton at n=1 has no nonconstant mode. These are gaps of the quotient refresh chain, not arithmetic speedups.
* **Source-bearing formula.** Uniqueness of the universal source makes its color orbit an invariant. Moving its color to c leaves exactly G_c available on the rest, proving the sum of a_(n−1)^(G_c). The source color and a marked port are different stabilizer problems. In particular source-free dimension is generally not harmonic dimension.
* **Complete relation objects.** Equal rowspaces give invertibly row-equivalent matrices even with dependent/zero rows, proving GL completeness. For affine tuples, subtract column one: the normalized rowspace is the augmented rowspace intersected with the coordinate-one-zero hyperplane. Equality of augmented spaces is therefore sufficient as well as necessary. For a fixed functional, choose it as first row: the residual transformations are [1 0; b B]. If its word w is zero, only r−1 rows can carry rank; if nonzero, classify the remaining rows modulo w and use b to remove residual multiples of w. The zero-word rank guard is essential. Empty affine tuples need the stated separate convention.
* **Typed small counts.** The transitive triple is rigid; the directed triangle has cyclic order-three automorphism group. Their decoration counts are respectively 16 and 8 (GL), 5 and 3 (AGL), and 50 and 20 (marked). Hence h_2^GL=5−2=3 and h_3^GL=24−5=19. These harmonic dimensions do not identify arithmetic mechanisms associated with prime labels 3 or 19.
* **Operation signature.** Exact observation fibers permit arbitrary independent permutations: (S_(2^k))^(2^r). An affine observer-preserver has kernel component Ak+Bs+c and unchanged syndrome, giving |GL(k,2)|2^(k(r+1)). Commuting with all named basis toggles forces F(x)=x+F(0); exact observation forces F(0) in ker M. Thus the hierarchy 16/8/2 for the three-prime cube is correct. Permuting update names only normalizes the update family and is a weaker requirement than pointwise commutation.
* **Neutral grid and character guard.** Preserving both named coordinate equivalence relations forces a separate permutation of each coordinate's classes, hence S3×C2, not S6. The squareclass equality has no implication for an unguarded Legendre-symbol test: sigma(81) vanishes modulo 11, sigma(3) does not.
* **Aliquot reconstruction.** For gcd(a,n)=1, expanding sigma(a)sigma(n)−an gives s(an)=s(a)n+sigma(a)s(n). Since s(a)>0 for a≥2, the exact pair reconstructs n. At a=3 the formula is n=s(3n)−4s(n), on 3∤n. This does force injectivity of any representation answering both exact probes. The sigma-squareclass observer instead sees multiplication by the square 4 and loses this information. The coprimality guard cannot be dropped.

## Alternative finite checks

`check_independent.py` is standard-library-only and does not import or execute either producer checker. It generates binary subspaces by span extension, enumerates tournament types and their automorphism groups, and counts relation-object orbits through degree four. It compares these totals with a separately implemented palette Burnside sum, and compares directly counted source-bearing decorations with the color-stabilizer formula. This extends the relation-object finite check one degree beyond the producer relation checker.

| Group | a_0,…,a_4 | source-bearing counts, n=1,…,4 |
|---|---|---|
| GL | 1,2,5,24,188 | 2,5,16,94 |
| AGL | 1,1,2,8,48 | 1,2,5,24 |
| marked port | 1,3,11,70,756 | 3,11,50,378 |

Independent degree-two refresh matrices are built by refreshing a winner/loser color tuple, including the new edge, without U/D or predicted eigenvalues in the construction. Exact rational nullities establish complete eigenspace multiplicities:

| Group | eigenvalue 0 | eigenvalue 1/2 | eigenvalue 1 |
|---|---:|---:|---:|
| GL | 3 | 1 | 1 |
| AGL | 1 | 0 | 1 |
| marked port | 8 | 2 | 1 |

Operational enumeration uses the 16 independently selectable fiber flips rather than scanning 8! permutations. The affine test is polarization, and the update test is direct commutation. It obtains 16/8/2. A separate six-object permutation enumeration yields 12 named-grid automorphisms. A divisor-incidence sieve, rather than multiplicative factorization, checks 3,485 coprime probe pairs, the shared-prime counterexample, the exact 5/19 probes, and the sigma collision 14/23.

```
python3 check_independent.py > normal.json
python3 -O check_independent.py > optimized.json
cmp normal.json optimized.json
```

Both executions passed and their receipts are byte-identical. Guards use explicit exceptions. The all-degree statements rely on the proofs, not extrapolation from these finite tests. No publishing was performed by this audit lane.

## Focused final integration review: rank does not preserve source deletion

The subsequently added operations §5 is correct. On a transitive tournament,
vertex order from its source is intrinsic. The rank-one tuples (v,0) and (v,v)
have different source-deletion ranks 0 and 1. The affine tuples (0,v,v) and
(0,0,v) have augmented rank 2, but source deletion gives augmented ranks 1
and 2. In least-significant-bit-first coordinate notation, their spaces are
{0,1}/{0,3} and {0,6,7,1}/{0,4,7,3}; deleting coordinate zero is x >> 1.
Thus equal current ranks do not imply equal ranks after the same intrinsic
operation. The claim is correctly limited to failure of the rank-only
operation-preserving quotient, not arbitrary algorithmic lower bounds.
The independent checker now constructs these spaces from the color lists
rather than taking the supplied masks as inputs. Normal and optimized
receipts were regenerated and remain byte-identical. No further mathematical
extension is made in this integration review.
