# Independent audit: second-band owner atlases

Session 116, 2026-10-30. Audited `lrc/SECOND_BAND.md` and inspected the separate producer checker. **PASS within the explicitly defined family and open upper chamber.** No global LRC or lower-chamber classification is established.

## Proof review

The family is exactly `C(q,A)=[1,q-1] union (q+A)`, with `1 in A subset [1,q-1]`. Its name is a definition, not an assertion that this is every historical meaning of “second-band.”

The inherited Session114 atlas for A={1} is the necessary global confinement step. Inside it, a new speed q+b has phase c/q at a primitive a/q. A phase c=1 occurs exactly for b=a^-1 modulo q, and supplies exactly the promoted lower wall. For c>=2 and old owner p>=2, the inequality

`2(v-p) <= p(c-1)(q+2)`

follows from `v<=q+cp` and `p(c-1)>=2`. With epsilon<2/[q(q+2)], it keeps the entire raw phase above r. Its upper bound follows because the phase decreases from c/q<=1-r. For p=1 the already-present owner q+1 is used; the displayed special estimate is valid for all q>=3. Reflection supplies the other side. Hence there are no hidden wrap components or untested additional walls. The proof genuinely covers every q>=3, not only the finite tests.

Only unit b can have phase +/-1 at a primitive center, so nonunit attachments are invisible in this chamber. Promoted owners retain their residues; all atlas bits are distinguishable. The forced promotion of 1 leaves phi(q)-1 independent bits, and the nonunit fiber has size 2^(q-1-phi(q)). A minimum u>1 forces precisely all lower units promoted and u unpromoted. If every unit is promoted the smallest owner is q+1. Thus the stated minimum multiplicities are complete. In particular composite minimum 4 for q7,A={1,2,3} is not an error.

For harmonic transfer, `2m>s` is exactly what makes m/(qm+s)>1/(q+2); also 1/(k+1)>=1/(q+1)>1/(q+2). Thus the claimed optimum is in the proven chamber. The least owner is unique, giving exactly two endpoints in the strict first-tail regime. The fully promoted case s=q+1 satisfies precisely the same argument.

For q tails, equality at 1/(q+1) requires a primitive tail phase. The owner-window bound is exactly ej<=m. Each directional owner occurs once, and disjoint windows and their two nonzero directions ensure the listed points cannot duplicate. When a core endpoint ties, ej=m implies e<=m<qm: the faster tail, not the core, determines terminal relaxation. This validates the weighted formula even at threshold m=s. Complete-period counting gives both stated affine quasiperiod increments; no asymptotic equidistribution assumption is needed.

Promotion surgery preserves owner residue/center and retains exactly the shorter j range. Therefore lost phase pairs and weighted loss are precisely the stated interval of j. The limiting losses are the difference `1/u-1/(q+u)=q/[u(q+u)]` times the universal density factors. The joint-bank criterion also follows: among all windows, the minimum owner supplies the widest available range for the common phase x=q|t-a/q|; the nonempty closed first-arrival set attains its infimum when finite.

The q7 character identity follows by substituting `w8(j)=88/105+(32/105)chi2(j)` into B and using N=2 times the number of admitted j. All four odd-eighth phases first occur at j=3 on a minimum-owner side, hence m>=3s. No connection to unrelated eighths is inferred merely from equal denominators.

## Independent exact tests

`check_independent.py` was written separately, before inspecting the producer checker. It imports no producer code. It intersects literal safe intervals for every individual integer speed using a two-pointer interval algorithm, and independently enumerates all opposite-slope affine intersection candidates to certify selected global maxima. All comparisons use exact rational arithmetic and explicit exceptions, not Python assertions.

- Every one of the 254 actual cores for q=3 through 9: full survivor intervals at three radii, including the top endpoint and a radius close to the excluded floor.
- Exact owner-atlas count, fiber sizes, and minimum-owner multiplicities for those moduli.
- 65 harmonic preserving/nonpreserving cases for every core at q=3,4,5, including fully promoted owners.
- Independent global maxima, complete strict-branch maximizers, and fastest-active-owner terminal slopes.
- Surgery survivor inclusion, count loss, and weighted loss.
- Composite-minimum and fully-promoted q7 eighth-phase thresholds and character formula.
- Exact count and mass quasiperiod increments for q=3,4,5, minimal and fully promoted cores.
- Literal hostile source {1,2,3,4,5,6,9} has global maximum 1/7, strictly below the old floor 1/5.

Normal and optimized runs must produce byte-identical receipts. The independent receipt records 1,152 passing gates. Finite exact tests supplement the proof; they do not themselves certify all q.

## Scope guardrails

No inherited positive floor is available when 2m<=s. Nonunit invisibility is not asserted below the open chamber, and arbitrary higher-band speeds are not authorized by this proof. The theorem includes the top endpoint but excludes r=1/(q+2). No defect requiring a producer change was found.
