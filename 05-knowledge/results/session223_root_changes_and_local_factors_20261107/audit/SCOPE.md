# Session223 convention and implication audit

Status: scope audit; inherited facts CITED or FINITE-EXACT. No all-target coverage proved.

## Source map
- Canon THM-4904 (`01-canon/theorems/THM-4904-root-change-word-fibers-carmichael-chamber-parity-law-sqrt-density.md`) and `05-knowledge/results/root_change_carmichael_chamber_20261107/{PROOF.md,audit/AUDIT.md}`.
- Original physical source and C=1 readout: `session199/physical-adjoint/PROOF.md`, sections 1–3; height: `session200/root-weighted-height/PROOF.md`, section 1.
- Complete C=1 grammar: `session197/lowcycle/C1_COMPLETE.md`, section 3.
- Chan–Pak precise formulation is CITED from `05-knowledge/results/chanpak_increments_thm013_thm017_merge_20261104.md`, section 0: Conjecture 6.2, Computational complexity of counting coincidences, TCS 2024, arXiv:2308.10214: every positive odd number other than 7 and 21 is H(T) for a tournament. This audit has not independently fetched the paper.

## Conventions
H counts directed Hamiltonian paths (vertex orders), not modulo reversal. C counts directed Hamiltonian cycles modulo cyclic rotation. The root t=(D,E,U,d)=(1,0,2,1) is a two-vertex rooted edge. Each chronological O/Z/J is one new vertex. D=H of current source, F=D-E+U is the H after appending Z or J. Only J changes the marked root.

M_K(N) counts terminal-Z words with <=K Js and terminal-J words with <=K Js, plus the empty word at N=1; trailing O words are omitted as they preserve D. For N>1, its support equals {F(t.w): #J(w)<=K}. This support equivalence is not a word-by-word bijection: a terminal-J counted word may use only K-1 Js in its prefix. To obtain strong C=1, append J to a prefix witnessing F=N. This needs at most K+1 physical root changes; the extra closing J is excluded from the SOURCE budget K. N=1 is exceptional and is not a C=1 readout.

The dagger reverses the source prefix and swaps O/Z. Its invariant is F, not the raw D of that prefix. Word multiplicities are not distinct numerical-profile, labelled-graph, or isomorphism-class multiplicities.

## Implications and missing steps
1. An effective theorem giving every odd N>B in F_{<=2} supplies a strong C=1 realization for that tail. It does not supply the small-target repairs.
2. If B<=10^9, the inherited census covers the intervening tail down to2647. Then only the 29 C=1 gaps need unrestricted/strong witnesses, apart from genuine exclusions7,21 and N=1. The list is 7,13,21,41,43,55,67,77,83,163,173,231,233,247,251,277,299,421,449,503,533,569,583,637,821,1067,1339,1861,2647. Existing strong finite spectra cover all other27, but cite or retain their witness certificates when closing the implication.
3. If B>10^9, additional verified coverage up to B is needed. A merely existential finite B yields only a cofinite statement, not exact {7,21} coverage.
4. Strong realization of every odd N>=3 except7,21 is stronger than the published unrestricted formulation. C=1 realization for every such target is FALSE (13 already obstructs it). Strong C=1 eventual coverage is compatible with the observed29 permanent C=1 holes.
5. Omission of7,21 for all tournaments is inherited from THM-1370 / THM-079 / THM-115, not newly proved by a restricted grammar missing them. General parity and H=1 transitive examples are separate boundary facts.
6. Singleton density zero, growing means, local residue laws and finite no-gap windows do not imply cofinite support. A modular witness is not an exact integer witness unless size bounds prevent wraparound and a shared physical word supplies all coordinates.

## Size efficiency
For a source prefix with k O/Z letters and r<=K Js, W=U-d>=2^k and F>=W+2d>=2^k+2. Thus F=N>=3 gives k<=floor(log2(N-2)). Its closing-J C=1 graph has exactly k+r+3 vertices, at most floor(log2(N-2))+K+3. This is a bound on every such witness, inherited from the W mechanism, not a search-time bound. For K=2 the bound is floor(log2(N-2))+5. The strongest global O(log N) statement requires finite repairs; any finite collection adds only a fixed bound. Enumeration of all words of logarithmic length is not automatically polynomial in bit length.

## Existing cutoff evidence
FINITE-EXACT, independently reproduced in the inherited audit: F_{<=2} has140 odd gaps up to10^9, largest3,780,967. F_{<=3} has44 gaps through4,000,000, largest27,781. F_{<=7} has29 gaps through10^6, largest2647. Overlap plus monotonicity in K gives every odd2647<N<=10^9 in F_{<=7}. These are target bounds across arbitrary word length (via growth pruning), not merely bounded graph orders. Depth7 examples1877,1897 have14-vertex J-closed witnesses and8 physical Js. No all-order conductor follows.

## Hostile tests for proposed all-order proofs
- Empty source; N=1; first targets3,5,7,13,21; all29 C=1 holes; largest F2 hole3,780,967; near-minimum multiplicities277,438,727 and301,189,927; exact-depth7 value1897.
- Every identity on arbitrary physical states versus only the root, and all zero-length block boundaries. In particular F(Z^a O^b)=(2^a+1)(2^b+1) requires a,b>=1. Zero blocks instead give F(O^b)=1+2^(b+1) and F(Z^a)=1+2^(a+1).
- Inspect whether inverse branches preserve a full physical predecessor recipe, root budget, parity and exact value, and strictly descend in a well-founded measure. Modular reachability alone fails this gate.
- A residue Markov law for fixed-length independently mixed O/Z blocks with inserted Js is a different probability space from F<=X or X<=F<2X weighted enumeration. Prove the change of measure, renewal/counting transfer and error bound before claiming the observed height-weighted 3-adic local factor or density. Positive limiting mass in every residue does not prove every integer is represented.
- Do not quotient words before counting multiplicities; JZ=OJ merges distinct words. Verify child D, parent F and closure cycle count with a literal builder, independently of transfer formulas.

## Follow-up proof audits

### Three-adic fixed-length law: PASS after periodicity repair
Reviewed `../threeadic/PROOF.md`, sections1–8. The d=0,U!=0 chain is periodic with period2, not mixing. The revised resolvent argument correctly handles this: the first row starts in the d=1 aperiodic class and has no pole at -1. Consequently the three-factor resolvent has order3 only at1, with coefficient pi_1 R Pi R Pi, and at most order2 at-1. Coefficient normalization gives O(1/L), including zero-length block boundaries. This is a rigorous finite-dimensional rational generating-function proof, not a fit.

Independent `audit_mod3.py` checks all81 states, constructs the stationary projector Pi and the signed -1 projector Q, verifies Pi^2=Pi, Q^2=Q, Pi Q=Q Pi=0, P Pi=Pi P=Pi, and P Q=Q P=-Q exactly. Independent recurrent-component discovery finds sizes1,1,1,6,18,18. Exact projection arithmetic reproduces one-J law(41/108,8/27,35/108) and two-J law(34/81,47/162,47/162). Output in MOD3_VERIFY.txt. The proof of aperiodicity in unit-d classes and period2 in the six-state class is valid. At this stage of review the higher-adic statement was only an exact transfer identity; the subsequent all-k proof is independently audited in HIGHER_AUDIT.md. Height-truncated HYP-9494(E) remains open.

### Transformation bridges: PASS, with a wording clarification requested
The full-state JZ=OJ identity and its iterated block identity are valid and credited as inherited. The two-run product rectangle, next-O defect2^b-2^a and positive-exponent domain are valid. The infinite graph on distinct integers2^a+1, a>=1, has exactly edge3–33: mod4 excludes exponents both>=2; factoring (u-1)(u+1)=2^(b-2) leaves only2,4. Prime-family distinctions and the composite Fermat counterexample are correct. Reflection composition and fixed-point statements are correct; the explicit closed nonsimple walk is an effective hostile control.

Ran the bridge worker's normal and optimized verification paths, identical outputs (3,280 source words,1,600 pairs,19,900 exponent pairs,16,380 reflection checks). The all-order arguments above were checked independently by symbolic substitution and elementary factorization, not inferred from those runs.

Applied clarification: the sentence saying a two-edge switch preserves degrees 'precisely when all four indicated edges are legal' mixes degree balance with graph legality. Degree balance is automatic for four distinct endpoints, regardless of square-sum legality. Write instead: 'For four distinct vertices, replacing edges(a,b),(c,d) by(a,c),(b,d) preserves every degree. To be a legal simple-graph square-sum switch, the removed edges must exist and the inserted edges must be permitted and absent from the remaining graph.' Connectivity still needs a separate check. This does not affect the product or one-edge theorem.

### Digital cutoff reformulation and selector obstruction: PASS
Reviewed `../cutoff/REPORT.md`. The block identity uses chronological bits with Z=1; the reversed digit weight in a is correct, the adjacent-pair correction in c is necessary, and the leading-bit initial E term is correct. Substituting J then C gives precisely the displayed alpha,beta,gamma, including C empty. Substituting middle B gives equation(3) with the correct plus beta*e*d contribution. B empty is handled separately and must not be fed to the nonempty bit definitions. A,C remain arbitrary words; thus this is an exact membership condition, not coverage.

The mixed-period drift proof is valid. The virtual orbit has all U/d strictly between0 and1 because a mixed binary cyclic string is neither all0 nor all1. The Z-run contribution is >2^(-l)>0, so delta_P>0. The first period sets E to its periodic value; the displayed formula is accordingly asserted only for m>=1. After any fixed suffix, the D-direction coefficient stays positive: the displayed(x,E,W,d) updates preserve the nonnegative E<=d cone, no nonzero vector can map to zero, and F=x+W+d>0 on its nonzero elements. The perturbation is correctly not called a physical source.

The transducer argument remains valid when input/output lengths differ or epsilon labels occur. Any loop on an accepting path is pumpable and therefore contains no J under a fixed budget. A mixed O/Z loop would give an output exponential-plus-constant and an input readout exponential-plus-nonzero-linear, impossible even when their exponential bases differ or an exponential coefficient vanishes. Trimming and label expansion reduce to finitely many SCCs, each monochromatic or source-empty, yielding bounded source runs. The inherited height bound then gives only polylogarithmically many source words at height X. The result excludes this literal-word encoding class, not output-only automata or general arithmetic algorithms, nor cofinite grammar support.

Normal and optimized verifier outputs agree:494 mixed periods,19,760 exact second-difference checks and28,830 digital two-J equations. The proof was checked independently of these regression runs. Prior-art wording must credit Session220's already asynchronous selector obstruction; only extension to literal bounded-J witnesses is new here. This attribution correction is applied in the published report.

### Higher-adic completion supersedes the earlier limited scope

The subsequent `../threeadic/ALL_K.md` has now been independently audited: PASS for the all-k projection and explicit valuation factors in the FIXED-LENGTH two-J ensemble. See HIGHER_AUDIT.md for independent group connectivity, branch derivation, exact factors, normalization and zero-residue tails. The earlier statement above that only a higher-adic transfer identity was proved described the first draft; it is superseded by this scoped analytic extension. The height-window bridge remains unproved, and no cutoff/coverage conclusion changes.
