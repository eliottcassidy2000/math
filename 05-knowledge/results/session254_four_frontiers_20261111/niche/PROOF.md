# Four further Laurent collision strips, and a universal forced factor

**PROVED IN THIS NOTE**, relative to the inherited exact carried-row formulas and small-block Smith theorem. The new sign certificates are exact symbolic identities in an unbounded parameter, not sampled extrapolations. General complementary width and positive physical-parameter noncancellation remain **OPEN**.

## 1. Inheritance and native target

At the pinned Session253 commit `44b3a8d8e17e3eae3ee4665ac3cc10734d341131`, `05-knowledge/results/session253_cross_problem_frontiers_20261111/moments/PROOF.md` proves the local small-block Smith exponents `(0,1,...,1)` at the carried Laurent collision x=−r and closes complementary widths ell=0,1,2,3. Its literal source is `05-knowledge/results/third_20260906_laurent.md`; no historical row proof is re-audited here.

We advance the specified next theorem through ell=7 and extract the all-width forced factor. This is a native boundary-module result, not a planar-JC proof or an integer-spectrum conclusion.

Write h=r+ell, z=2r+1, and use the inherited monic complementary polynomials

    P(t)=sum_(i=0)^ell ell!/[i!(3ell−3i)!]
                   product_(k=2i+1)^(2ell)(z+k) t^i,
    Q(t)=sum_(i=0)^(2ell) (2ell)!/[i!(6ell−3i)!]
                   product_(k=2i+1)^(4ell)(2z+k) t^i.

The complementary multiplication block is invertible iff R_ell(r)=Res_t(P,Q) is nonzero. The physical r are integers r≥2, but the certificates below hold for every real r≥0.

## 2. Universal forced-factor theorem

For every integer ell≥1, define

    F_ell(r)=product_(j=1)^ell (r+j)^j (2r+2j+1)^j.

Then F_ell divides R_ell in Q[r], and degree R_ell≤4ell². Consequently

    B_ell(r)=(-1)^ell R_ell(r)/F_ell(r)

is a polynomial of degree at most 3ell²−ell.

**Proof of divisibility.** At r=−j, z=1−2j. Every coefficient of P below t^j contains the zero factor z+(2j−1), whereas P is monic. Thus t^j divides P. Every coefficient of Q below t^(2j−1) contains 2z+(4j−2), hence Q is divisible by t^(2j−1). In particular P and Q have a common factor of degree at least j.

At r=−j−1/2, z=−2j. The corresponding zero factors z+2j and 2z+4j imply t^j|P and t^(2j)|Q. Again their gcd has degree at least j.

To upgrade common-root multiplicity to parameter divisibility, represent R_ell as the determinant of multiplication by Q on the free Q[r]-module Q[r,t]/(P), of rank ell. Over the residue field at either parameter value, its nullity equals the degree of gcd(P,Q), hence is at least j. Elementary row/column operations with unit determinant isolate a matrix of corank≥j over the parameter DVR; its determinant is divisible by the local parameter to power j. The two families of parameter values are pairwise distinct. Their product therefore divides R_ell. This argument proves a LOWER bound on those factor multiplicities, not equality in arbitrary width.

**Proof of degree bound.** The coefficient of t^i in P has degree≤2(ell−i) in r, and the corresponding coefficient of Q has degree≤2(2ell−i). Put r=1/s and t=u/s². The scaled monic polynomials s^(2ell)P(u/s²) and s^(4ell)Q(u/s²) have polynomial coefficients in s. Their u-resultant is s^(4ell²)R_ell(1/s), also polynomial in s. Thus degree R_ell≤4ell². Subtracting degree F_ell=ell(ell+1) gives the residual bound. QED.

This sharp weighted degree bound replaces the earlier convenient but looser bound8ell². It does not require a root-separation theorem.

## 3. Exact positive certificates at ell=4,5,6,7

The new symbolic computations give respectively

    ell                 4    5    6    7
    degree R_ell       64  100  144  196
    degree B_ell       44   70  102  140
    positive coefficients45   71  103  141.

Every coefficient of each B_ell is strictly positive. Complete rational coefficients in descending order are stored in `certificate.json`. The script constructs P,Q from the displayed formula, computes the exact resultant, divides by F_ell exactly, verifies that the quotient is a polynomial, and tests every coefficient sign. Thus the assertion is a finite algebraic identity whose parameter r remains unbounded. Six independent literal Sylvester-determinant evaluations per width guard the implementation; they are controls, not the proof of the polynomial identity or positivity.

Therefore for every r≥0 and ell∈{4,5,6,7}, R_ell(r) has sign (−1)^ell and is nonzero. Combining with the inherited small-block theorem, for every integer r≥2 and h=r+ell with ell∈{0,...,7}, the COMPLETE h-dimensional carried-response Smith exponents over Q[[x+r]] are

    (0 repeated ell+1 times, 1 repeated r−1 times).

Its determinant has exact collision order r−1, inverse precision loss exactly1, and cokernel (Q[[x+r]]/(x+r))^(r−1). The zero-width through three-width strips are inherited; widths4 through7 are new here. No claim about arbitrary ell follows from this finite list of unbounded-r strips.

## 4. Exchange ledger and next experiment

The polynomial residual B_ell is the exact missing sidecar between a controlled small collision block and the full matrix. This complements the planar-JC row16 ledger: in both cases an explicitly computed response controls some directions and leaves a precise residual condition, but there is no asserted map between their coefficients. The literal centralizer/Chan–Pak bridge remains the one established in Session253; extending these strips does not extend its physical positivity or tournament-realization scope.

The natural all-width conjecture is coefficientwise positivity of B_ell for every ell≥1. The stronger useful theorem would be structural complementary coprimality for every r≥2, perhaps through a separation or positive determinant identity. Coefficientwise positivity is sufficient, not necessary. The cheap next probes are ell8 and extracting a recurrence in ell; endlessly increasing width without a mechanism has diminishing value.

Run `python3 probe.py` from this directory. It requires SymPy, uses exact rational arithmetic and writes the complete certificate and concise log. This task does not replay any inherited census or outside claim.
