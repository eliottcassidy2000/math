# Independent audit: all one-reset letters

**PASS, scoped.** Reviewed the finalized `../collatz/GENERAL_RESET.md` independently of the producer's finite-check implementation. This supplements, and does not replace or weaken, the earlier `AUDIT.md` pass. No arbitrary-controller bound or Collatz convergence statement is certified.

## Proof review

Let A=F_1, M=F_(1^L a), a>1. Their affine commutator against A^j has valuation h+z_j, where h=v5(2^a−2), z_j=0 for odd j and z_j=1+v5(j) for even positive j. The stated coefficients P,Q,B and affine order T_e are correct.

For t>1 the multiplier of M minus1 is a unit. Conjugating M to its multiplier shows that membership of A^j in its cyclic subgroup is equivalent to the simultaneous slope-subgroup and fixed-point conditions. The slope's residue is always admissible because t=2 or4; its principal-unit condition is 1+v5(j)>=min(r,e). The fixed-point condition is exactly h+z_j>=e. When e<=h the latter disappears; otherwise it forces even j and the additional exponent threshold. These give precisely both stated t>1 formulas for g_e, including the nontrivial odd-j intersection at shallow depths.

For t=1 write s=v5(P−Q), b=v5(B), r=min(s,b). A power M^n has slope and translation depths s+v5(n) and b+v5(n) (until trivial modulo5^e). If s differs from b then h=r and the intersection is trivial. If s=b=r, matching the two coordinates of A^j gives the additional threshold z_j+h−r>=e along with z_j>=min(r,e). This yields J=max(min(r,e),e+r−h). In actual reset words J>=1: h=0 forces r=0, and h>0,t=1 forces L odd and r>=1. Therefore the stated even-j formula is correct. Also r<=h follows directly from the valuation of (P−Q)−B.

The resulting g_e divide successive g_(e+1), tend to infinity, and give exactly the three stated distinct-value radix chains. In particular the jump is10, not2 or5, when t>1 and r>h>0. Repeated g_e do not introduce extra digits.

A fixed entry selects at most one coset of the cyclic intersection at each depth. Empty entry selects residue0 and simultaneously attains every layer cap1+floor(L/g_e). Summation using the inherited ordinary all-layer theorem gives the exact maximal mean; eventual growth makes the correction sum finite. This is not an interchange of unproved truncation limits.

An equality entry must have an infinite-depth phase: otherwise all sufficiently high layers vanish while the maximum has positive mass. That phase commutes rationally with M; faithfulness of positive affine words and the primitive one-reset word give exactly1^j u=w^b. Such an entry selects residuej at every depth, whose count is1+floor(j/g_e)+floor((L−j)/g_e). Equality is therefore exactly the stated residue inequalities. These are genuinely equivalent to digitwise inequalities in the *same mixed radix*, because any excessive digit contributes more than all possible lower-digit deficits. Their product count is correct. No prime-characteristic Lucas identity is claimed for the mixed radix10 branch.

## Independent finite controls

`check_general_independent.py` uses only standard-library code and literal affine cyclic-subgroup enumeration, with no producer or predecessor imports. It checks orders, the least positive intersection exponent, all small-entry layer bounds, simultaneous empty-entry maxima, every completion-phase count, divisibility chains, digit predicates and phase counts. The receipt records the finite universes and total. Normal and optimized Python outputs are byte-identical; every gate raises an explicit exception on failure.

These finite controls corroborate the proof but are not substitutes for its universal arguments. They preserve the distinct scopes of original-clock rates, per-controller-block means, and the separate parameter-Haar statistics.
