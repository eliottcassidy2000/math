# Independent audit of the odd-eighth synthesis

**PASS after the character-label correction, verified in the final draft. Scoped proof audit and separate exact checks.** The two actual cross-field maps, Lucas equality representation, and typed parameter-Haar interpretation are valid with the character-label qualification below. This is not formal proof-assistant certification or a literature-priority assessment.

## 1. Affine-order events and named Boolean readings

For h of odd order l acting on a binary vector space, N_h=1+h+...+h^(l−1) is the projection onto Fix(h). The affine map(v,h) has even order exactly when N_hv is nonzero: its l-th power is the corresponding translation, whose order is either1 or2. Therefore the Boolean event f_h(v)=1[N_hv≠0] preserves the actual predicate before averaging, and the even-order mass is the average of its named zero-character readings.

In dimension3, odd matrix orders are1,3,7. Their event functions are respectively raw OR=1−delta_0, a nonzero linear form, and zero; their readings are7/8,1/2,0. The complement groups of orders1,3,7,21 give the displayed four mixtures. This is a genuine event-function map, not equality of unrelated measures.

**Character-label correction before publication:** the raw identity-event function OR is not the uniquely normalized Fano gauge representative. That representative is s=OR+x1+x2+x3, corresponding to A=B=C=D=1. Thus p_OR(0)=p_s(111)=7/8, whereas p_s(0)=3/8. More generally, subtracting the linear normalization r·v transports the observer from0 to r. For order-three matrices, normalization produces zero but the retained nonzero character has reading1/2. Discarding this named character is exactly the demonstrated loss. The relative gauge class or unordered profile alone cannot recover group mass.

## 2. The physical LRC weight through word-address transport

On odd j modulo8, chi_2(j)=sqrt(2)cos(2pi j/8), with values(+1,−1,−1,+1). The word-address permutation3→5 sends primitive numerator frequencies1,3,5,7 to7,1,3,5. Applying it to the two-character cosine expression gives

    U chi_2 = −i chi_(−2),

where chi_(−2)=(+1,+1,−1,−1) on the four odd inputs. Consequently the positive reciprocal-owner weight88/105+(32/105)chi_2 becomes nonreal. The obstruction is exact on the actual phase labels. It is not merely a failed floating-point comparison or a claim about full-space Dirichlet-character extensions.

The finite-algebra congruence gate p=q mod2^(A−1) correctly supplies positivity-preserving pullback structure for this particular operator. Even there, no map between the actual speed sets or physical LRC dynamics has been constructed. The source packet states this boundary.

## 3. Lucas phases and parameter-Haar means

For one-reset controllers, the inequalities j mod(2*5^k)<=L mod(2*5^k) separate into a parity inequality and all lower-digit comparisons. A first offending base5 digit cannot be cancelled by lower digits, giving the digitwise criterion. Frobenius factorization then identifies exactly the nonzero binomial coefficients modulo5. The word completion and phase remain needed to reconstruct the entry; counting alone loses them.

The constant-controller3/8 is an integral of the stated maximal block-mean expression over a profinite exponent parameter space, conditional on an even exponent. At finite precision E the truncated extended mean is3/8−1/(16*5^(E−1)). Passing to the Haar integral is justified by the geometric tail. This is not an ordinary Cesaro average over initial positive integer exponents, not an actual controller with noninteger valuation, and not the clock-normalized rate of one source.

## 4. Independent verification

`check_cross_maps.py` imports no producer implementation. Its **5,923 explicit gates** check:

- Exact symbolic eighth-root transport and the nonreal image of each of the four positive weights.
- All105 odd-order matrices of GL3(2), literal affine orders, norm projections, truth-table weights, and normalization/character transport.
- All four complement mixtures by actual element orders.
- Every phase j for0<=L<=100, against direct binomial residues and the digit-count product.
- Exact finite parameter-Haar layer proportions through precision5 and their truncated mean formula.

Normal and `python -O` output match byte for byte. Reproduce with `python check_cross_maps.py` and `python -O check_cross_maps.py`; `receipt.json` records the finite universes. Python3 and SymPy are required. The written proofs, not these finite tests alone, establish the all-parameter claims.

## 5. Final scope refinement (wording audit)

The final synthesis correctly separates two gates on the odd-eighth algebra: positivity-preserving pullback requires p=q modulo4, whereas invariance of the specific physical weight w8 requires p=q modulo8. Inside the first gate, c=q/p is1 or5 modulo8; c=5 reverses chi2 and exchanges the two positive weight levels. Outside it the weight image is nonreal. The producer's expanded exact checker covers all16 odd residue pairs and was rerun successfully (1,136 gates); these are additional producer checks, not part of this audit's independent5,923-gate receipt.

The atlas now states the all-a>1, L>=0 one-reset classification, with the full intersection modulus g_e and its mixed-radix equality rule. The formulas transcribed in SYNTHESIS.md agree with GENERAL_RESET.md. Ordinary Lucas/binomial support is explicitly restricted to the2,5,5,... radix case; the other cases are not claimed to be single-prime binomial support. No universal general-reset clock bound is asserted. This final addition is a transcription and scope audit; the separate general-reset auditor supplies its independent proof and exact-check review.
