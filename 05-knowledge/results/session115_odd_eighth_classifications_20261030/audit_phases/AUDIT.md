# Independent audit: Boolean phases and word-address transport

**Session 115. PASS for the scoped written proofs below.** These are independent algebraic reviews and separately implemented exact checks, not formal proof-assistant certification or literature-priority assessments. No major open conjecture is resolved.

## 1. Coherent norm-table phases

Reviewed `quadratic/ODD_EIGHTH_PHASES.md` against the inherited Session112 norm-sign conventions. The gauge values s(3)=A, s(5)=B, s(6)=D, s(7)=C+D give precisely the polynomial with quadratic coefficients A,B,D and cubic coefficient A+B+C. The finite arithmetic selection map and its squareclass30 promise use the earlier exact bank; they do not extend to arbitrary integers or repeated prime multiplication.

The complete spectrum proof is correct. A nonzero normalized quadratic has alternating polar rank2, giving four Walsh coefficients of magnitude4 and four zeros. Fourier inversion at zero fixes three positive and one negative nonzero coefficient. Each cubic is an affine phase with one point toggled. Its exceptional Walsh coefficient is +6 except for the zero-point defect class, where it is −6. This yields the two distinct odd-eighth profiles, not all four probabilities in every table.

The relative GL(3,2) quotient action has exactly four orbits of sizes1,7,7,1. The point defect transforms with the coordinate change, and the base multiplication table must be transported too. The result is not a bare fixed-base table-relabeling statement and is not an abstract octonion-isomorphism invariant. These qualifications are load-bearing and are correctly retained.

For arbitrary d>=3, the extremal test |W_f(r)|=2^d−2 is equivalent to exactly one disagreement with one affine phase (possibly after its constant flip). The exceptional frequency and its sign identify the affine part; the remaining signs identify the defect. Uniqueness also follows because distinct affine words differ at least2^(d−1)>=4 places. Thus there are exactly2·(2^d)^2 such functions. Exhaustion of all odd-weight functions is special to d=3; the stated d4 weight-three hostile is valid.

`check_boolean_independent.py` starts from integer truth-table masks and uses an independent fast Walsh transform. It exhausts every Boolean function at dimensions3 and4, compares the extremal test against independently enumerated affine-radius-one words, and checks the16 normalized classes, literal GL actions and arithmetic selection values. **66,103 gates pass normally and under optimization with identical output.** There are128 defects among256 functions at d3, and512 among65,536 at d4. The d4 odd-weight hostile has maximum Walsh magnitude10, not14.

## 2. Word-address algebra classification

Reviewed `unitary/WORD_UNITARY_CLASSIFICATION.md`. The domain is the entire finite character algebra through depthA, not just primitive frequencies. The defining map reindexes a Fourier basis by valuation-word labels; it is not automatically composition on the source variable.

The necessity proof is valid: the one-letter word(A) is a primitive character and forces any multiplicative map to scale all frequencies by q/p. The word(1,A−1) then has defect2(p−q)/(p²q), exactly requiring p=q modulo2^(A−1). For sufficiency the first word-sum term cancels, while every later term contains both p−q and an additional factor2. All depths at mostA therefore obey the same scalar law.

The odd-residue algebra has the same necessity: the primitive character generates it with its half-period power equal to−1. Distinct characters are still distinguished at source point1. The inherited transport respects the signed pairing of frequencies differing by half a period, so it descends to this algebra rather than merely acting on redundant symbols.

At first failure A=v2(p−q)+2, only the second word term can survive. It contributes the half-modulus precisely when the first valuation is1, equivalently the source numerator k equals q modulo4. This proves the stated defect formula. Opposite odd numerators receive complementary defects. A real cosine consequently becomes nonzero purely imaginary on the odd coset, proving failure of real preservation. Adding a constant makes the input nonnegative without removing that nonreal image. Thus positivity and source-permutation pullback fail at the same first depth. Earlier-depth algebras embed into later ones, so these properties do not return. The sufficiency direction is literal odd scalar pullback.

The global statement is properly scoped to the union of finite character algebras and extensions where the named operations are defined. It neither rejects ordinary parity-coding conjugacies nor proves integer-orbit convergence.

`check_unitary_independent.py` uses literal shortcut dynamics to decode every source residue to its stopped valuation word, with no producer imports. It verifies permutation, signed odd-coset descent, multiplicativity, real preservation, explicit scalar maps and first-defect formulas for every odd multiplier residue pair through depth7; selected full product tables are also checked. The receipt records **44,313 gates**; normal and optimized outputs agree byte for byte.

## Reproduction

Run both audit scripts in this directory with ordinary Python and with `python -O`. They use the standard library and explicit exception gates. Their receipts separate finite universes from the all-parameter proof statements.
