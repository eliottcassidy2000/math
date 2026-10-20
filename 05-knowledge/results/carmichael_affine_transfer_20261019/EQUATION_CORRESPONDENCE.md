# Equation-by-equation transfer ledger

Source read: entire supplied 9-page PDF and its full extracted text. Equations on printed pages 3–6 were additionally rendered and visually checked. Source: Howe–Huffman–Kohlmeyer–Kruer, *The exponent of growth of the Carmichael counting function*, October 8, 2026. The paper credits Larsen's pruning/Fourier strategy; the present notes make no external novelty claim for classical Fourier, coding, Krawtchouk, or coset-counting ingredients.

Scope: a mathematical reuse ledger, not a portable Lean audit. Source claims in (S),(A), nonvanishing, and the final formal verification remain source claims here. Elementary transfer results below have standalone proofs in the companion files.

## Introduction and input interface

- (1.1), exponent-one limit: preserved as the *conditional* asymptotic consumer of the finite compiler. No independent analytic proof is supplied.
- (1.2), finer Pomerance asymptotic: no transfer. Even the exact-cardinality entropy improvement supplies neither its secondary term nor effective y-thresholds.
- (2.1), smooth prime reservoir: the actual role is to provide N large distinct labelled multiplicatively independent integer atoms with locally compatible divisibility. Distinct labels are indispensable for injectivity; arbitrary repeated group elements are allowed only because the labels remain distinct primes. The group lemmas require no primality.
- (2.2), modulus size: keep the logarithm visible, log L<=y^delta log(5y), rather than replace it by y^(delta+sigma). This sharpens the sufficient Bernoulli and exact-k mixing inequality from eta>2delta+4sigma to eta>2delta+2sigma. The coordinating lane independently identified this sharpening.
- (2.3), ambient character estimate: replaced abstractly by (A*) in FIXED_CARDINALITY_TRANSFER.md. It has two distinct consumers: translated short-arc bounds and a quotient-index bound. Neither consumer requires equidistribution of the selected P or T.
- (2.4), global zero-free half-plane: source claim and source-dependency boundary only. It is not used directly by any finite-group theorem here.
- (2.5), weighted prime-character sum: source claim and analytic bridge only. A proposed replacement estimate must retain modulus uniformity before it can feed the finite transfer.

## Modulus and pruning

- (3.1), L,G and labelled P,A: abstract source is an ambient finite positive counting measure on a finite abelian group. Concrete targets are the binary strip control columns, tournament incidence columns, and arithmetic correction columns. These maps are displayed, not inferred from shared vocabulary.
- (3.2), log-lcm budget: the useful entropy debt is log|G|, not the huge value |G|. Keeping this logarithmic debt exact is what makes the sharper parameter condition work. A generic finite integer matrix uses log of its finite image/cokernel size instead.
- (3.3), D,U,R,A: retain D=ceil(log_2|G|); use U0=y^delta log(5y); keep R proportional to ambient-to-selected density Q/N. The R^2 phase penalty is a quantitative density cost, not an unexplained universal spectral gap.
- (3.4), size bounds: the coarse U=y^(delta+sigma) is sufficient but needlessly loses sigma twice in the mixing-scale inequality. It is not used in the sharpened compiler.
- (3.5), beta+2delta+2sigma<1: retained because it makes A0=o(N), separately from the mixing exponent condition.
- (3.6), (R-1)A0<=Q: retained as one convenient sufficient bound. The weighted-power version uses Q+2 sum_(j<R)(1-j/R)A_j directly.
- (3.7), nontrivial support after pruning: strengthened to non-concentration on *every character fiber* by affine pruning. This is necessary for exact cardinalities; pruning only the identity phase leaves all-equal-label hostiles.

## Character phases

- (4.1), character energy: independent selection uses E_chi=n-Re sum chi. Exact-size selection requires Delta_chi=n-|sum chi|. Their difference is mathematically consequential: all equal nonidentity labels have E_chi>0 and Delta_chi=0.
- (4.2), Fejer polynomial: use its translate F_R(exp(-i phi)xi(g)) for every phase phi. Absolute ambient Fourier bounds survive the phase factors. This transfers phase anti-concentration to the direction of the selected sum.
- (4.3), ambient Fejer mass: the same positivity argument works for arbitrary positive ambient weights. Its second independent use, ordinary quotient orthogonality, yields [G:H]<=(Q-A0)/(n-A0). This bounds the otherwise uncontrolled cardinality period.

The improved finite dispersion constant is

    Delta_chi >= (N/R^2) min(2/D,3/(4pi^2)),

rather than coarsening both branches to N/(100DR^2). The low-order branch uses unequal-pair counts and the exact squared-distance identity, while the high-order branch uses translated Fejer arcs. The asymptotic scale improves because of logarithm retention, not because of these constants.

## Random products and counting

- (5.1), 2delta+4sigma<eta: sharpened to 2delta+2sigma<eta after retaining logarithms. The upper bound eta<1-sigma ensures sampling remains sparse.
- (5.2), size-window count: replaced by a binomial lower bound binom(n,k)/(2|H|) for each sufficiently mixed exact k compatible with the affine coset. Its main missing datum is h=ord(aH), and the index bound proves h=O(Q/N).
- Unnumbered independence identity before (5.3): literal transforms (A1),(T3),(L4) in APPLICATIONS_AND_EXACT_CODES.md. The invariant is a dual-code weight in exponent-two targets, a cut size for tournament parity, and a character evaluation on correction columns for lifts.
- (5.3), Bernoulli modulus identity: retain a common phase theta before extracting the kth coefficient. The resulting exact-size character bound is b_(n,k)^(-1)exp[-rho(1-rho)Delta_chi]. Claiming the old independent product after conditioning would be false.
- (5.4), identity probability: ordinary character inversion holds for every target, not just 1. The affine compiler restricts possible targets to a^kH. In exponent-two groups, Fourier positivity additionally proves identity mass at least uniform without any energy assumption, although nonidentity targets and exact sizes need not have that property.
- (5.5), Chernoff tail: unnecessary for the exact-k compiler. It remains useful for identity-positive exponent-two or inversion-paired reservoirs, where it alone can produce a large size-window family.
- (5.6), upper integer size: with k<=y^eta, use X(y)=(5y)^(y^eta), eliminating the (1+nu) collar.
- (5.7), entropy: use log binom(n,k)>=k log(n/k), retain log n>=(1-sigma)log y-log2, and subtract log|H|. The resulting limiting exponent is 1-sigma-eta.

## Parameters, limit, and formal scope

- (6.1), parameter choice: any fixed choices satisfying beta+2delta+2sigma<1 and 2delta+2sigma<eta<1-sigma suffice. Parameters can shrink for the conditional exponent-one result. There is no claimed effective dependence of the threshold on them.
- (6.2), collar loss: the fixed-size route removes (1-nu)/(1+nu); its exact entropy loss is sigma+eta at the exponent level.
- (6.3), liminf-to-all-X step: retain the source's continuous strictly increasing X(y), but with no cardinality collar. The integer k(y) can jump because the upper size X(y) is defined independently from it. Thus the argument still reaches all sufficiently large X, not only a sequence.
- Section 7 and Appendix A: provenance and verification scope only. The selected 320-module revision and author-side axiom report are not substituted for a clean portable build. Nothing here promotes the analytic source assertions to independently established theorems.

## New hypotheses with explicit completion gates

1. **Fixed-weight Pascal mixing:** quantify the least t,k giving near-uniformity on the parity-compatible affine hyperplane. The exact Krawtchouk transform is proved; a uniform sharp bound is open here.
2. **Cut-spectrum fixed-k mixing:** for connected edge-reversal graphs, determine a structural criterion beyond non-bipartiteness that bounds every relevant centered cut energy. Bipartiteness is exactly the additional affine obstruction; it does not measure bottlenecks.
3. **Sparse first-lift multiplicity:** find natural integral matrix families whose correction-column reservoirs satisfy strong phase dispersion. First-lift solvability and binary correction supply are different obligations.
4. **LRC ambient phase reservoir:** an actual positive measure and target-preserving map would be needed to use translated Fejer anti-concentration. No such map is established; no theorem about LRC follows from the shared kernel.

## Reproduction and scope

Run from the containing workspace:

    python3 session8-character-transfer/verify_transfer.py > session8-character-transfer/verification.json
    python3 -O session8-character-transfer/verify_transfer.py > session8-character-transfer/verification_optimized.json
    cmp session8-character-transfer/verification.json session8-character-transfer/verification_optimized.json

The exact universe is frozen in the script: cyclic sampled label lists; exact affine pruning of labelled subfamilies of uniform cyclic reservoirs; all binary characters of twelve Pascal codes; all input words in the specified finite strips; all labelled tournaments through six vertices; all connected simple labelled graphs through five vertices; and 300 admissible integer first-lift instances. It includes independent direct enumeration and dynamic-programming paths. The Fourier inequalities and translated Fejer evaluations are numerical checks with explicit tolerance, separately typed from the integer certificates.
