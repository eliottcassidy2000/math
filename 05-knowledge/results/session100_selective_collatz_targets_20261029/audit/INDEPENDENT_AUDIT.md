# Independent audit: selective source lifts, clocks, and controllers

Source repository resolved HEAD at audit start: `92514b16e2f6a13b7a1f194a6e17d5e5c5fc9ee7`.
Targeted reads: the three Session100 lane notes and their checkers; no broad historical crawl.

**Verdict:** the central three-lane theorems withstand independent proof review and exact finite controls. They are reframes/obstructions, not a Collatz proof. Two minor wording qualifications below should be incorporated before publication.

## Proof review

1. **Bijection domain.** Excluding the backward orbit of `-1/b` is necessary; arbitrary positive-integer valuation words correspond precisely to the remaining odd 2-adic sources. Oddness of each tail proves exact valuations and excludes finite termination. This is a 2-adic conjugacy for both parameters, not a 5-adic conjugacy.
2. **Exact cylinders.** Backwards substitution from an odd terminal gives modulus `2^(D+1)` and residue `(2^D-C)b^(-J)`. Modulus `2^D` alone does not specify an odd endpoint and is insufficient. Compatible least representatives are nondecreasing; boundedness therefore means eventual constancy, and only a positive odd ordinary integer can have a bounded sequence. This criterion does not presuppose reaching 1.
3. **Selective target.** `Y_(3->5)=Phi_5(W_3^+)` retains every positive 3-source, while its target attractor is `-1`, the image of 1. This is an exact restatement of Collatz, not a simpler theorem already supplied by symbolic conjugacy.
4. **Denominators.** For reduced `p/q`, the numerator `2^a p-q` is coprime to odd `q`; this is true even when `gcd(b,q)>1`. Hence the reduced denominator after inversion remains a multiple of `q`. The hostile denominator 101 persists under every finite prehistory and all rotations. The proof excludes the entire specified eventual-periodic basin, not arbitrary nonconvergent itineraries.
5. **No finite safety filter.** Each exact cylinder contains infinitely many positive sources. Every alleged rejection triggered by a finite prefix also rejects such sources. This excludes safety-style forbidden-prefix rules, not all infinite-word automata acceptance conditions. The distinction in the note is correct.
6. **Clock barrier.** Multiplication of stepwise ratios gives a conditional endpoint bound if all intermediate states exceed the threshold. Without that condition the strict source barrier gives *some* earlier descent, not necessarily descent at the final endpoint. Independent witness: `b=3,n0=165,k=17,S=27`, endpoint 167; the strict barrier holds but only prior descent is guaranteed.
7. **Limsup and Haar guard.** If 1 is avoided, odd positive states are at least 3, proving the `log_2(10/3)` upper bound; reaching 1 makes the mean clock exactly 2 asymptotically. Thus the strict limsup equivalence is sound. Chernoff bounds are ambient-density controls; no pointwise integer conclusion follows from them. The finite-certified-basin variant is conditional on an actual complete certificate, as stated.
8. **Controller trap.** The shifted potential has the stated ratios. For `(5,5)`, 19 terminates from any positive queue in exactly `Q+3` steps; 83 preserves every positive queue over its 3-cycle. Avoiding both yields uniform branch contraction, infinitely many branches are necessary to remain active, and the positive lower potential bound contradicts that contraction. Thus the classification is genuinely arbitrary-queue, not limited to queue 1.
9. **Finite guard.** When `Delta>0`, all branch inputs outside the finite set E admit a uniform contraction by monotonicity in input and valuation. Any infinite run must revisit E infinitely often. Universal positive-queue certification on E is essential; queue-1 certification alone would not justify the conclusion. Empty E is sufficient. Controller return termination still does not settle source-marked descent after resets.

## Wording corrections

- Bijections §2: replace “negative 5-adic rational sources” by “negative rational sources in the 2-adic 5x+1 system.”
- Controllers §1: qualify uniqueness of the shift *from valuation-1 leaves* by `r>=2`. At `r=1`, valuation 1 is a branch, so the proposed characterization is vacuous. The formula `h=1/(b-2)` remains well defined and later cases/Delta-positive finite-guard theorem are unaffected.

## Reproduction and scope

`python3 -B audit/check_independent.py` and `python3 -B -O audit/check_independent.py` both PASS, with identical output. No assertions are optimized away.

- 117 exhaustive exact cylinders for `b=3,5,9`, word lengths 1..3, digits 1..3, enumerating every odd residue; each includes a half-modulus countercheck.
- 1,452 prehistories: `b=3,5,9`; terminals `19/101,-1/3,5/9,7/15`; lengths 0..4, digits 1,2,7.
- 5,000 native controller runs: both `(3,3)` and `(5,5)`; odd sources 1..999; starting queues 1,2,3,7,19. Repeated native `(n,Q)` certifies finite experiment traps independently of the proposed special-input diagnostic.
- 15,716 triggered source-clock barriers for `b=3,5,7`, odd sources 1..501, first 60 steps; all have prior descent.

These are finite implementation checks. Universal statements rest on the reviewed arguments above. Audit owns only `audit/`; lane authors/parent incorporate wording corrections and handle publication.

## Atlas checker pre-publication review

The atlas checker also passes (`atlas_replay.json`). Three orders `(1,1,2)`, `(1,2,1)`, `(2,1,1)` have carries `b²+2b+4`, `b²+2b+8`, `b²+4b+8`, giving successive-difference ratio `b/2`. Two-letter permutation differences are `2^a-2^c`, independent of multiplier. Exact guard cylinders remain distinct, so differences compare *affine charts*, not alternative trajectories of one source.

For a branch-preserving affine candidate `h(x)=Ax+B`, `A!=0`, the identity `h((bx+1)/2^a)=(c h(x)+1)/2^a` first implies `b=c` from x-coefficients; the a=1 and a=2 constant identities then give `B=0,A=1`. On actual branches, the infinite exact-cylinder point sets justify extending each equality to an affine identity. The checker samples identities off-domain at x=1,3, which is legitimate only as polynomial-identity controls, not as native-orbit examples. This rigidity is limited to branch-preserving affine maps; the nonaffine 2-adic conjugacies in the bijection lane are not contradicted. The final atlas prose was subsequently reviewed as recorded below.


## Final atlas proof review

`atlas/AFFINE_ATLAS_FINGERPRINT.md` has now been reviewed in full. The stronger **rational-function** rigidity theorem is valid: characteristic zero makes the eigenvalues `2^(-j)` distinct for different integer Laurent degrees; `R(y/2)=R(y)/2` therefore leaves only degree 1 in the Laurent expansion. The formal Laurent embedding of rational functions is injective, so this proves `R(y)=Ay` as a rational identity, not merely locally. Nonconstancy gives `A!=0`; the original constant and linear coefficients force `A=1` and `b=c`.

The passage from native-cylinder equality to rational identity is justified: clearing denominators gives a polynomial vanishing at infinitely many inputs, and a nonzero rational function has only finitely many poles. The two branch domains need not overlap. Finite affine samples do not independently test the arbitrary-rational theorem; the theorem rests on this valid formal proof.

The three-row observer is accurately scoped to permutation-chart differences at common input, with branch guards retained. Ratio `b/2` is invariant under a common invertible affine output-coordinate change. The minimal-length claim concerns precisely this observer: it neither claims all two-slot information is blind nor claims arbitrary observers require three slots. Final verdict: PASS, with earlier wording corrections incorporated by parent. No further mathematical correction requested.
