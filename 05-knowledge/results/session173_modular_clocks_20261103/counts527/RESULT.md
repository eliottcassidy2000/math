# The 527 composition polynomial and the normalized-44 contrast

Status: PROVED algebraic corollaries conditional on the explicitly inherited finite inventories below. Those inventories are FINITE-EXACT in the repository and have prior independent audits; this packet does not freshly reproduce their exhaustive enumeration. No external priority claim.

## Sources and conventions

All paths refer to eliottcassidy2000/math, inspected at commit 1b91742ce55a202b73c80d56561ba56d00b503a6.

- `05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md`, Sections1–2: SCC generating series, Busch bound and eventual onset.
- Same directory, `low_H_formulas.json`: exact strong17 and31 inventories.
- `04-computation/experiments/h_spectrum_multiplicities_20261021/audit_S50/inventory_375.json`: independently produced strong17,31,89 inventories, including self-converse flags.
- `04-computation/experiments/h_spectrum_multiplicities_20261021/strong_ears_625.out`: 125 total strong counts for odd375..623. Index76 is527 and equals184.
- Same experiment, `audit_S50/b3_ear_census_625.py` and `.out`: a second, independently implemented ear census confirms all125 totals entry by entry. It does not print the order histogram for527.
- `05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities/n8_classes.csv`:22 strong order8 classes with527 paths.
- `01-canon/theorems/THM-4744-interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve.md`: inventory completeness and strong minimum f(13)=625.
- `01-canon/theorems/THM-4749-the-repunit-degree-law-phi-2-of-2-is-the-triangle-phi-3-of-2-is-the-gap-and-fermat-products-reverse-the-order.md`: inherited factorization-degree law.

N(n,h) counts ordinary isomorphism classes with actual directed Hamiltonian-path count h. It is neither labelled count nor H/|Aut|. Let A_h(z) count nontrivial strong ordinary classes by order.

## 1. Exact symbolic answer for527

Inherited inventories give

    A_17(z)=2z^6,
    A_31(z)=2z^6+2z^7.

Since527=17*31 with distinct primes, the only nontrivial strong core words are one527 atom or the two orders of a17 atom and a31 atom. Ordered SCC decomposition is unique. Singleton SCCs fill the gaps freely. Thus the all-order series is exactly

    F_527(z)=A_527(z)/(1-z)^2 + 8(z^12+z^13)/(1-z)^3.       (1)

Write a_m=[z^m]A_527(z), S=sum a_m, W=sum m*a_m. The inherited625 census supplies S=184. The exact maximum order is at most12 because f(13)=625>527. The n<=7 census has maximum189; a_8=22 from the order8 census. Therefore

    A_527(z)=22z^8+a_9z^9+a_10z^10+a_11z^11+a_12z^12,
    a_9+a_10+a_11+a_12=162.                              (2)

The individual four coefficients, and hence W, are NOT recovered here. The existing printed source gives only their sum. They must not be guessed from total184.

For every n>=11, coefficient extraction from(1) gives

    N(n,527)=8(n-11)^2+184(n+1)-W
            =8n^2+8n+1152-W.                            (3)

Indeed the two-component contribution is
8*C(n-10,2)+8*C(n-11,2)=8(n-11)^2. At n=11 both binomials vanish as required, although the second formal upper argument is zero; for n>=12 this is ordinary stars-and-bars. The single-component contribution is sum a_m(n-m+1); at n=11 the only absent order12 term is exactly zero. This proves the onset11 directly. It is a guaranteed onset, NOT asserted minimal: cancellations could allow earlier polynomial agreement.

Useful exact consequences independent of W:

- The quadratic coefficient is8 and the linear coefficient is8.
- N(n+1,527)-N(n,527)=16(n+1), n>=11.
- The second difference is16.
- 1634<=W<=2120, hence -968<=1152-W<=-482, using only(2). These are bounds on admissible coefficient data, not claims that extremal distributions occur.
- The irreducible information still missing is one first moment W, not the entire census.

If one does not assume the inherited S=184, the formula remains valid with
8(n-11)^2+S(n+1)-W. Thus the census dependency is isolated cleanly.

## 2. Normalized44 means actual89, with a complete answer

Even44 cannot be H(T), by the oddness theorem. If the intended observable is (H-1)/2, then44 corresponds to89. The inherited audited inventory gives

    A_89(z)=14z^7+7z^8+2z^9,
    C_89(z)=z^8,

where C counts strong self-converse classes. Since89 is prime, only one nontrivial SCC is possible. Therefore

    F_89=A_89/(1-z)^2,
    N(n,89)=23n-149 for every n>=8.                      (4)

This onset is minimal: at n=7 the actual count is14, while the polynomial gives12. For n>=8, self-converse counts are1 for even n and0 for odd n, because the squarefree self-converse core formula gives

    J_89(z)=z^8/(1-z^2).

Consequently converse-merged counts are (23n-149+1_even(n))/2 for n>=8. These formulas distinguish three natural meanings of a tournament count rather than identifying them.

## 3. A genuine mod8/mod9/mod72 connection

Put C=1152-W. Formula(3) gives, for n>=11,

    N(n,527)=C mod8,
    N(n,527)=C-n(n+1) mod9.

Its eventual residue sequence modulo72 has EXACT least period9, regardless of the unknown moment. A period d would require the increment16dn+8d(d+1) to be divisible by72 for every sufficiently large n. Consecutive n force72|16d, hence9|d. Conversely d=9 works because the increment is144n+720.

In contrast, N(n,89)=23n-149 has exact least period72 modulo72, since gcd(23,72)=1. This is a genuine compositional contrast: semiprime527 yields a quadratic polynomial with a rigid9-step residue clock; prime89 yields a linear permutation of all72 residues. No identification with the proposed alphabet {-1,0,1,2} is asserted.

## Reproduction and limitation

`python check.py` and `python -O check.py` use standard-library arithmetic. They verify source extraction, all-order series versus polynomial identities for each allowed source order, the four extremal placements of the162 remaining strong atoms, period claims, the89 boundary and converse parity. They do not independently enumerate tournament isomorphism classes. `recover_inventory.py` is an unfinished recovery instrument adapted from the prior auditor. Recovery was attempted but not completed: its launch stopped at missing numba before computing any inventory. Nothing from the subsequently installed isolated dependency directory was run.
