# Addendum: why the ternary projective carrier is exceptional

Status: PROVED elementary consequence of the audited Gaussian allocation bijection. This is separate from the frozen squarefree support note. No new F1 descent or hourglass relation is claimed.

Fix an odd prime ell, put a=(ell-1)/2, and choose distinct integer split primes p_1,...,p_k with p_i=1 mod 4, k>=1. The auxiliary prime ell need not equal any p_i. Set e=product p_i^a and fix an order and Gaussian conjugate orientation for those factors.

The inherited Gaussian theorem identifies D_e with the nonzero integer box [-a,a]^k modulo simultaneous negation. Reduction modulo ell is a bijection from that box to F_ell^k, and it commutes with negation. Consequently the exact indexing object is

(F_ell^k minus {0})/{+1,-1},

of cardinality (ell^k-1)/2. It is NOT in general P^{k-1}(F_ell).

There is a natural surjection from this signed quotient onto projective space, obtained by quotienting by all nonzero field scalars rather than just the two signs. Every projective line has ell-1 nonzero vectors, and signs act freely, so every fiber has exactly (ell-1)/2=a points. Thus this surjection is bijective exactly when ell=3. For k=0 both sets are empty; that vacuous case is excluded from the 'exactly when' statement.

The distinction is structural, not a coordinate-permutation convention. On a fixed nonempty support S of size r, the signed quotient has (ell-1)^r/2 elements, while the projective stratum has (ell-1)^(r-1). Their ratio is again a. Supports survive the quotient, but absolute integer exponent magnitudes generally do not: a scalar in F_ell acts by multiplication followed by taking balanced representatives in [-a,a]. This can change |j_i|.

The exact valuation statement makes the arithmetic loss visible. For every index j,

gcd(e²,d(j))=product_i p_i^(2(a-|j_i|)).

For j_i!=0 this follows from the inherited exact local valuation; for j_i=0 the inherited lower bound is at least 2a and the gcd truncates it to 2a. Thus unlike the squarefree ternary case, the support alone need not determine this gcd, and projectivization need not preserve it either. A coordinate-zero deletion still permits exact division by the corresponding p_i^(2a), with the remaining Gaussian indices unchanged, but arbitrary projective scalar identification does not.

## Minimal hostile

Take ell=5, a=2, k=1, actual split prime p_1=5, and e=25. Using pi_1=2+i gives:

- j=1: z=5(2+i)²=15+20i, gap d=600, v_5(d)=2, gcd(625,d)=25.
- j=2: z=(2+i)^4=-7+24i, gap d=336, v_5(d)=0, gcd(625,d)=1.

These are the two distinct signed classes in F_5^*/{±1}, but both map to the unique point of P^0(F_5). Endpoint squares are 625±600=5²,35² and 625±336=17²,31². So the projective quotient already destroys exact divisibility for one coordinate.

For k=2 and ell=5, using center e=5²·13²=4225 gives 12 gaps, whereas P^1(F_5) has six points; every point represents two gaps. It would therefore be false to extend the squarefree formula by replacing 3 with ell while still claiming a gap/projective-point bijection.

## Map and scope ledger

Source: bounded integer Gaussian indices modulo sign, with fixed conjugate choices. Target 1: nonzero finite-field vectors modulo sign, via reduction; this is a genuine set bijection retaining coordinates through balanced representatives. Target 2: projective points, via full scalar quotient; this identifies a distinct signed classes per point, and can lose |j_i|, gcd data and exact gap values. Necessary sidecar: the scalar class in F_ell^*/{±1}, or equivalently the original signed index class. No map identifies integer exponent addition, Gaussian phase multiplication, or numerical gap addition with field addition. The special ell=3 case is special because the full scalar group already consists of signs, independently of any S_3 symmetry of projective-line coordinates.

Reproduce: `python3 verify_uniform_exponents.py` and `python3 -O verify_uniform_exponents.py`. All gates use explicit exceptions; both runs agree exactly. Tests cover ell=3,5,7, k=1,2, the full gcd formula, every projective fiber, explicit endpoints, and independent endpoint enumeration for tested centers <=5000. These finite tests support, rather than replace, the proof above. The only inherited input is the Gaussian index/valuation theorem in `session58-hourglass-prime-support/RESULT.md`, §§3–4; no new external theorem is required.
