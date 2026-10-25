# Independent modular-magic count audit and all-modulus theorem

Status: PROVED (elementary parametrization and finite Smith-form inclusion-exclusion); independently checked exact controls. No literature priority claim.

## Scope and inheritance

Count ordered, labelled 3-by-3 arrays over Z/nZ whose eight row, column, and diagonal sums agree and whose nine entries are pairwise distinct. In the all-modulus theorem entries are unrestricted: they need not be square residues. The audited prime-power packet instead imposes square-residue entries, zero allowed. Neither counts choices of roots, quotients by rotations/reflections, or proves an integer-square magic array.

Anchor: independently audit the prime-power exact formulas. Niche: recover all-modulus ordinary counts from collision kernels. Wildcard: determine the exact quasipolynomial period. Concept board: ring parametrization, collision arrangement, Smith torsion, affine unit-square shifts, valuation descent. Closest mechanism is the producer's division-free parametrization; hostile examples are n=1,2 and characteristic 3, where polynomial continuations and division by small integers fail. The corrected near miss is applying the eventual 2-/3-power polynomial at every depth. The least-used sidecar is the full integer Smith invariant of a collision intersection rather than only its rank.

## Theorem: every modulus

Let g_d=gcd(n,d), for n>=1. The labelled ordinary count is

A(n)=n^3-(4g_2+4)n^2+C(n)n,

where

C(n)=4g_4+4g_6+2g_2^2+2g_5+g_2g_4-2-4g_2.

This is a cubic quasipolynomial of minimal period exactly 60, including small moduli without exceptions.

### Parametrization without division

An array with top corners a,b and center c necessarily is

(a, 3c-a-b, b;
 c-a+b, c, c+a-b;
 2c-b, a+b-c, 2c-a).

For completeness, if its line sum is s, the diagonals first determine the bottom corners as s-c-b and s-c-a. The two outer columns then determine the middle outer entries as c+b-a and c+a-b. Their row forces s=3c. The remaining entries follow. Conversely substitution checks the eight lines. Thus these arrays are bijective with (c,a,b) over every commutative ring.

Writing a=c+u and b=c+v, the nine offsets are

0, u, -u, v, -v, u+v, -u-v, u-v, -u+v.

Translation by c does not affect distinctness and contributes the factor n.

### Collision kernels and exact intersections

Every difference of two offsets, up to sign, is one of

u, v, u+v, u-v, 2u, 2v, 2u+2v, 2u-2v,
u+2v, u-2v, 2u+v, 2u-v.

The kernels of the first four are contained in those of their doubles; each of the remaining eight forms really is a difference of offsets. Hence the forbidden set is exactly the union of the kernels of these eight integer row vectors:

(2,0), (0,2), (2,2), (2,-2), (1,2), (1,-2), (2,1), (2,-1).

The first four singleton kernels have size n*g_2 and the last four have size n. Every subset of at least two rows has rational rank two, since the eight slopes are different over Q. If its Smith invariants are d1,d2, its kernel over Z/nZ has size gcd(n,d1)*gcd(n,d2). Indeed unimodular integer row and column operations remain invertible modulo n, reducing the equations to d1*x=d2*y=0. The invariants can be computed without choosing reductions: d1 is the gcd of all entries, and d1*d2 is the gcd of all 2-by-2 minors.

Summing the inclusion-exclusion signs (-1)^|S| across subsets of cardinality at least two, the nonzero aggregate coefficients are:

- (d1,d2)=(1,1): -2
- (1,2): -4
- (1,4): 4
- (1,5): 2
- (1,6): 4
- (2,2): 2
- (2,4): 1

All 247 subsets are explicitly recorded in certificate.json and regenerated using only integer gcd and determinant arithmetic in verify.py. This finite identity, combined with the kernel formula and the singleton sum, gives the stated theorem for all n, not merely for the enumerated moduli.

### Exact period

The coefficients depend only on n modulo lcm(2,4,6,5)=60. More explicitly, for each residue r in {0,...,59}, use the polynomial

X^3-(4*gcd(r,2)+4)X^2+C(r)X

on n congruent to r modulo 60. The full 60-entry C(r) vector is in certificate.json, making this an explicit quasipolynomial description.

There is no period collapse. Polynomial coefficients on arithmetic progressions are unique: compare two proposed quasipolynomial representations on their common residue classes modulo the least common multiple; equality at infinitely many n implies equality of polynomials. Therefore any quasipolynomial period must be a period of the coefficient functions, in particular C. The least positive period of a periodic sequence divides every period, hence divides 60. Each proper divisor of 60 divides one of 12,20,30. But

C(60)=56, C(72)=48, C(80)=40, C(90)=44,

so none of 12,20,30 is a period. The minimal period is 60.

### Boundaries and connections

At n=1,2,3,4,5,6,7,8 the theorem gives zero, as required by nine-entry distinctness. At n=9 it gives 216. Specializing to powers of 2 gives n(n-4)(n-8) for n>=4, with A(1)=A(2)=0; powers of 3 give n(n-3)(n-5) for n>=3, with A(1)=0. At primes p>=7 it gives p(p-1)(p-7); p=5 is exceptional because the determinant-five intersections retain extra torsion.

The map to the collision arrangement preserves distinctness and forgets only the free center c, restored by a factor n. Smith rank alone loses torsion information; d1,d2 are the required sidecar. This explains why the coefficient changes at factors 2,3,5. Counts are not multiplicative: A(5)=A(7)=0 but A(35)=33600. Distinctness modulo a product does not require distinctness in either factor separately. The theorem imposes no square-residue constraint and supplies no global integer solution.

## Independent audit of prime-power square counts

Audited frozen source: source-snapshots/PROOF.md, verify.py, certificate.json. SHA-256 hashes are in this audit's certificate. Verdict: PASS, without required mathematical corrections.

1. Division-free parametrization and exact labelled scope are correct, including characteristics 2 and 3.
2. The ordinary polynomials and their exceptional low levels agree with the all-modulus theorem above. The producer's three-line collision clusters at p=3 and primitive strata at p=2 independently justify the same counts.
3. Opposite square entries sum to twice the center. Modulo 3 this forces all entries into one of the classes 0 or 1; nonunit squares are divisible by 9. Modulo 4 it likewise forces all entries to be even squares or odd squares; the former are divisible by 4 and the latter are 1 modulo 8 when k>=3.
4. Square-unit classes are exactly 1 modulo 3 for p=3 and 1 modulo 8 for p=2. The affine maps x->1+3x and x->1+8x preserve distinctness and equivalently preserve line-sum equalities at the reduced modulus. Thus unit counts are A(3^(k-1)) and A(2^(k-3)), with unavailable small levels zero.
5. Division of all entries by p^2 gives a bijection for the nonunit part to the square problem two levels lower. Equality modulo p^k becomes equality modulo p^(k-2), so no lifting factor is missing.
6. Therefore N_p(k)=U_p(k)+N_p(k-2), with N_p(0)=0 and low levels zero. Its exact generating functions are

G_2(z)=1536*z^7/[(1-z^2)(1-2z)(1-4z)(1-8z)],

G_3(z)=216*z^3*(1+27z)/[(1-z^2)(1-3z)(1-9z)(1-27z)].

7. The leading ordinary term is p^(3h). After shifts s_2=3 and s_3=1, the normalized parity-compatible descent sum tends to p^(-3s_p)/(1-p^(-6)). The lower-degree terms vanish after normalization. Hence the densities are exactly 1/504 and 27/728. Existence thresholds are k=7 and k=3, not claims about characteristic-p extension fields or integer squares.

## Reproduction and independence

Run `python3 -B verify.py` and `python3 -O -B verify.py` in this directory. All checks use explicit RuntimeError guards, not removable assert statements. The verifier recomputes all 247 Smith subsets, checks collision equivalence and exact ordinary counts for every modulus 1 through 120, and counts square arrays independently by bitset intersections of nine translated square-residue sets (powers 2 through 2^9 and powers 3 through 3^5). This census loops over all offset pairs and counts permitted centers, rather than the producer's square-pivot triple enumeration. Finally it multiplies the exact sequences by the GF denominators through degree 50.

The finite controls support the proof and expose scope errors; they do not substitute for all-level proof. No remote operations or shared-repository edits are performed by this audit. Parent publication is responsible for integration. The startup coordination command reported missing .machine-id; no identity was invented.
