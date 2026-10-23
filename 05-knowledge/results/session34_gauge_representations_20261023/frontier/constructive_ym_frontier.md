# Session 34: constructive Yang–Mills frontier audit

Checked 2026-10-23 UTC. Read-only public-source inheritance audit, not proof certification. No user-supplied PDF was opened. All working downloads are public repository files; no Git mutations occurred.

## Sources and status

- Official current repository: https://github.com/ethanwsussman/Yang--Mills-LLM-experiment . README describes a continuing LLM-assisted construction project, with QED3/YM3 expositions and an apparent YM4 argument. This is author-reported progress, not a validated solution.
- Original proposal: https://github.com/ethanwsussman/Yang--Mills-LLM-experiment/blob/main/README.pdf . Downloaded and read. It separates existence, axioms, and mass gap, and explicitly warns that raw output contains gaps/mistakes. Its stated snapshot date is October 6, 2026.
- Theorem 0: https://github.com/ethanwsussman/Yang--Mills-LLM-experiment/blob/main/theorem_0_complete_account_raw_I.tex . Downloaded public 10,554,720-byte source; inspected hypotheses, theorem statement, and explicit status discussion, not its full proof or Parts II–III. Master file confirms the three-part title. Chapter 0 explains status vocabulary.
- Current rotation manuscript: https://github.com/ethanwsussman/Yang--Mills-LLM-experiment/blob/main/RP%20O(d)%20trick%20(RP1)/RPM1.pdf . Downloaded and inspected Theorem A and scope qualifications, pp. 1–3. Source: https://github.com/ethanwsussman/Yang--Mills-LLM-experiment/blob/main/RP%20O(d)%20trick%20(RP1)/RPM1.tex . RP1 README expressly says correctness has not been satisfactorily verified: https://github.com/ethanwsussman/Yang--Mills-LLM-experiment/blob/main/RP%20O(d)%20trick%20(RP1)/README.md . This is the current linked target, not verified identical to the older filename YM_Rotation_Invariance_simplified.pdf.

## Theorem 0: exact constructive scope

Part I, H1–H7 and theorem label 3.1:thm-main: four-dimensional Wilson SU(N) lattice gauge theory on periodic tori; odd sufficiently large blocking factor; ordered auxiliary choices; sufficiently small reference and bare couplings; sufficiently large torus; first-passage RG cutoff; additional source-radius/coupling thresholds. Clauses (0)–(i) address coupling flow/UV stability for N≥2; substantive curvature-limit clauses are stated for SU(2).

The announced targets include cutoff-uniform separated action-density correlations, subsequential continuum/infinite-volume limits, coordinate reflection positivity, a restricted reconstruction, a two-point nontriviality witness, and tuned finite-torus uniqueness up to a coupling parameter. These are NOT all established unconditionally: the explicit OS-status discussion makes correlation/limit claims implications from named incomplete or unproved inputs. Full rotation covariance in clause (vi) additionally assumes an unproved vanishing far-sphere rotation-current flux. The source does not assert a mass gap, clustering, general field operators/Wightman functions, or coincident-point/tempered-distribution control. Its two-point witness is expressly not a proof excluding generalized-free behavior.

## RPM1: different conditional bridge

Theorem A assumes an existing Borel probability law on real Schwartz distributions with coefficients in a finite-dimensional orthogonal O(d) representation, d≥2, planar weights in −2,…,2. It requires translation invariance; E exp(F(χ)) ≤ exp(C(||χ||₁+||χ||₂²)) for every real test; and reflection positivity for all coordinate and diagonal mirrors x_j=0 and x_j±x_k=0. It concludes O(d) covariance only of moments with pairwise disjoint supports. Noncovariant contact terms are explicitly allowed. Lattice symmetry alone is explicitly not RP.

This source presents a proof, but the repository's warning precludes describing it as independently verified here. Importantly, Theorem 0's separated-support, locally finite-order moment bounds do not supply RPM1's full Borel law/exponential bound, nor does its coordinate RP alone supply diagonal RP. One cannot simply substitute this manuscript to remove Theorem 0's infrared hypothesis.

## One next obligation informed by finite positivity

Prove an actual reflection-compatible characteristic-functional bridge for a specified Yang–Mills lattice observable family: for every coordinate AND diagonal mirror H and every finite collection of supported tests, identify

K^(a,V)_{ij} = E[exp(i F_(a,V)(χ_i − Θ_H*χ_j))]

as PSD, with the correct reflected tensor/gauge action, and show entrywise convergence to a continuum characteristic functional as a→0 and V→∞. PSD then passes to the limit by the elementary closedness of finite PSD cones. This is one precise bridge obligation, not a claim that convergence exists.

Our prior graph-extension Gram positivity and finite gauge/representation calculations inform only the algebraic factorization/checking step. They do not identify these matrices with the Yang–Mills characteristic-functional matrices. Gauge-group representation symmetry is not spacetime O(4) symmetry. No finite rank, eigenvalue, or graph-gluing calculation supplies tightness, uniform exponential bounds, continuum convergence, or a physical mass gap. Those analytic inputs remain separate; even the bridge above would not complete the construction.

Stopping condition reached: exact source scopes and a single non-automatic transfer obligation identified. No new proof hunt undertaken.
