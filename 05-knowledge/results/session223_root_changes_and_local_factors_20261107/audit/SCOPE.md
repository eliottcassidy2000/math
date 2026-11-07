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
