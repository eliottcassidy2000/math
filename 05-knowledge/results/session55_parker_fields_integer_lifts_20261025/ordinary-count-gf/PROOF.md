# Ordinary modular magic counts: exact rational generating function

Status: PROVED corollary of the audited all-modulus count, with exact polynomial verification. No external novelty claim. This is a separate continuation packet; the preceding audit files are unchanged and pinned by SHA-256 in source-pins.json.

## Scope and dependency

A(n) counts ordered, labelled 3-by-3 magic arrays over Z/nZ with nine pairwise distinct entries, with all eight row/column/diagonal sums equal. Entries are unrestricted residues, not required to be squares. There is no quotient by dihedral symmetry, count of roots, identification with tournament counts, or integer-square existence claim.

The source theorem is session55-modular-count-audit/PROOF.md, “Theorem: every modulus,” independently verified by its associated script and certificate. It gives

A(n)=n^3-(4gcd(n,2)+4)n^2+C(n)n,

C(n)=4gcd(n,4)+4gcd(n,6)+2gcd(n,2)^2+2gcd(n,5)
      +gcd(n,2)gcd(n,4)-2-4gcd(n,2).

This packet takes that proved all-n count as input; its new work is the exact indicator expansion, rational GF, reduced poles, and period interpretation.

## 1. Divisibility-indicator expansion

Write I_d(n)=1 when d divides n and zero otherwise. The gcd formulas give

gcd(n,2)=1+I_2,
gcd(n,4)=1+I_2+2I_4,
gcd(n,5)=1+4I_5,
gcd(n,6)=1+I_2+2I_3+2I_6.

Using I_2^2=I_2 and I_2 I_4=I_4,

gcd(n,2)^2=1+3I_2,
gcd(n,2)gcd(n,4)=1+3I_2+4I_4.

Substitution yields exactly

C(n)=7+13I_2+8I_3+12I_4+8I_5+8I_6,

and therefore

A(n)=n^3-8n^2+7n+I_2(-4n^2+13n)
      +8n I_3+12n I_4+8n I_5+8n I_6.

No asymptotic approximation or small-n correction is involved.

## 2. Rational generating function

Let F(z)=sum_(n>=1) A(n)z^n, a formal power series over Q. Define

L_1(t)=t/(1-t)^2,
L_2(t)=t(1+t)/(1-t)^3,
L_3(t)=t(1+4t+t^2)/(1-t)^4.

These satisfy L_j(t)=sum_(k>=1) k^j t^k. Consequently sum_(n>=1) n^j I_d(n)z^n=d^j L_j(z^d), and

F(z)=L_3(z)-8L_2(z)+7L_1(z)
     -16L_2(z^2)+26L_1(z^2)
     +24L_1(z^3)+48L_1(z^4)+40L_1(z^5)+48L_1(z^6).

The reduced rational form is

F(z)=8z^9 Q(z)/D(z),

D(z)=(1-z)^4(1+z)^3(1+z+z^2)^2(1+z^2)^2
     *(1+z+z^2+z^3+z^4)^2(1-z+z^2)^2,

Q(z)=27+37z+119z^2+201z^3+346z^4+438z^5
     +639z^6+649z^7+729z^8+635z^9+576z^10
     +400z^11+300z^12+160z^13+94z^14+34z^15
     +14z^16+2z^17.

Thus the numerator has degree 26 and denominator degree 27. In cyclotomic notation D=Phi_1^4 Phi_2^3 Phi_3^2 Phi_4^2 Phi_5^2 Phi_6^2, since Phi_1=z-1 occurs to an even power. All factors are essential. In particular the initial coefficient is A(9)=216 and A(1)=...=A(8)=0.

## 3. Exact poles, including no cancellation at third roots

At z=1, L_3(z) has leading singular term 6/(1-z)^4. All other terms have order at most three. The pole order is exactly four.

At z=-1, only -16L_2(z^2) can have order three. With w=1+z, its leading singular term is -4/w^3, nonzero. The pole order is exactly three.

For a primitive d-th root zeta with d=3,4,5,6, write w=1-z/zeta. Whenever zeta^e=1, L_1(z^e) has leading term 1/(e^2 w^2). The order-two leading coefficients are consequently:

- d=3: 24/3^2 + 48/6^2 = 4.
- d=4: 48/4^2 = 3.
- d=5: 40/5^2 = 8/5.
- d=6: 48/6^2 = 4/3.

Every coefficient is nonzero. All other summands are analytic at the respective roots, except those already included in the displayed sum. This proves exact order two at each such root. There are no poles elsewhere: each summand has denominator a power of 1-z^e with e in {1,2,3,4,5,6}.

In particular there is NO pole at a primitive 60th root. A common quasipolynomial period is an lcm of the root-of-unity orders of the pole locations, not necessarily the order of any individual pole location.

## 4. Minimal quasipolynomial period is nevertheless 60

The divisibility formula immediately gives a quasipolynomial of period lcm(2,3,4,5,6)=60. Conversely, a quasipolynomial of degree three and period m has a GF whose poles are all m-th roots of unity: sum each constituent polynomial over its residue class, obtaining a denominator dividing (1-z^m)^4. Since the reduced F has primitive fourth-, fifth-, and sixth-root poles, every such m must be divisible by 4,5,6, hence by 60. Therefore the minimal quasipolynomial period is exactly 60.

This also follows from the source audit's coefficient-uniqueness argument: C(60)=56 differs from C(72)=48, C(80)=40, and C(90)=44, excluding all maximal proper divisors of 60.

## 5. Why 2, 3, and 5 appear, and what does not transfer

The source's collision arrangement has row vectors

(2,0), (0,2), (2,2), (2,-2), (1,2), (1,-2), (2,1), (2,-1).

Its inclusion-exclusion intersections retain Smith invariants among (1,1), (1,2), (1,4), (1,5), (1,6), (2,2), and (2,4). Their gcd kernels explain all divisibility coefficients and the lcm 60. More concretely, determinants 4,5,6,8 arise from pairs of collision rows; determinant 8 in the doubled-diagonal pair has Smith type (2,4), so the largest relevant 2-primary invariant is 4 rather than 8. Examples are determinants 5 from (1,2),(2,-1), 6 from (2,2),(1,-2), and 8 from (2,2),(2,-2). These are rank-two linear collision torsion data.

This mechanism is not a Parker quadratic-discriminant or square-residue prime exception. Square membership is not imposed anywhere in A(n); its divisibility dependence already occurs in unrestricted ordinary arrays. Identifying a shared numerical prime set would not supply a predicate-preserving map. No transfer to a Parker square construction, to a tournament enumeration, or to an integer-square solution is claimed.

The exact connection is: ordinary arrays -> offset-pair complement of a linear collision arrangement -> Smith-kernel inclusion-exclusion -> divisibility-indicator counts -> cyclotomic rational GF. The center is forgotten by the offset representation and restored by a factor n; square membership is not preserved because it is not a predicate in this carrier.

## Verification

Run `python3 -B verify.py` and `python3 -O -B verify.py` in this directory. Only the Python standard library is used. The checker constructs the denominator, combines the nine rational summands using exact polynomial division over Q, checks the complete numerator polynomial, and computes nonzero remainders modulo every cyclotomic factor to certify that no factor cancels. It independently verifies the C-indicator identity on all 60 residue classes and convolves the exact all-modulus count with D through degree 360. The persisted certificate contains all exact polynomial coefficients and cyclotomic remainders. Load-bearing tests raise explicit exceptions, so optimized replay retains them.

The original audit packet is immutable; source-pins.json records its exact dependency hashes. No external publication-priority claim is made.
