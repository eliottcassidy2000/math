# Three external preprints merged: MUB(6) = 3, local smoothing with Kakeya in every dimension, and "ETH is false" — typed hooks to the Eisenstein/Gauss-sum lattice, the 4D Kakeya program, the classification atlas, and the dipole law

**Session:** klein-2026-10-27 (worktree `codex/session-overnight-20261027`), 2026-10-27.
**Status of everything below:** the three papers are **EXTERNAL PREPRINT CLAIMS** (author line "OpenAI"; the openai/math
collection that [seven_twentyone_mersenne_openai_math_20261006.md](seven_twentyone_mersenne_openai_math_20261006.md) read
on 2026-10-06). Nothing in them was re-derived here. Three parallel readers produced page-cited digests; a web search on
2026-10-27 found **no independent public confirmation or refutation** of any of the three. Every connection to the
repository is typed (IDENTITY / TOOL / ANALOGY / NO TRANSFER / COINCIDENCE), as the precedence rule requires.

**Owner's prompt (verbatim, with three uploaded PDFs):** "merge in recent results on The maximum number of mutually
unbiased bases in dimension six. merge in prior results on the kakeya conjecture and quantum chromodynamics. pursue the
new reflection argues that every distance lower bound in the repo comes from one identity: a move negates the
tournament's matrix and adds the path. Where the bound fails to be exact, a forced symmetry explains why. also the
Order-11 prime cycle-index types and the order-12 ledger".

**A mismatch, stated plainly.** The uploads are `paper-6.pdf` = *Local smoothing for the Euclidean wave equation*
(2026-10-25, 59 pp.), `paper-7.pdf` = *The Exponential Time Hypothesis is false* (2026-10-17, rev. 10-25, 67 pp.) and
*The maximum number of mutually unbiased bases in dimension six* (2026-09-24, 76 pp.). No quantum-chromodynamics paper
was attached. "Prior results on QCD" is read as the repository's own QCD/root-lattice analogy
([metagraph_hall_bott_20261023 §6](metagraph_hall_bott_20261023/README.md),
[metagraph_five_port §5](metagraph_five_port_phase_frontier_20261023.md)); section 4 merges it with the dipole law.
The ETH paper is merged as what it is (section 3). The dipole census and the order-11/12 ledger are in the companion
packet [order11_prime_types_n12_ledger_dipole_census_20261027/](order11_prime_types_n12_ledger_dipole_census_20261027/README.md).

## 1. Mutually unbiased bases in dimension six (N(6) = 3)

**Claim (Theorem 1.1, computer-assisted).** Under the arithmetic and compiler contract of its Appendix A (IEEE 754
binary64, round-to-nearest, no contraction, a disassembly check that two `pow(.,2)` calls were lowered to multiplication
in the one recorded binary), a documented complete execution (4 October 2026 UTC, 100 native invocations, 74 minutes
wall) certifies that no four pairwise mutually unbiased orthonormal bases exist in `C^6`; with the classical lower bound
of three (Lemma 1.2) this gives `N(6) = 3`, Zauner's 1991 conjecture. Corollary 1.3: the maximum number of pairwise
Killing-orthogonal, `*`-closed Cartan subalgebras of `sl_6(C)` is three.

**What is proof and what is certificate.** All lemmas of its sections 2–6 (equivalences and phase gauges, the nine-phase
arrangement with seed grid `N = 112`, the vanishing-sum filter, seed balls of radius `0.032`, completion over a seed ball
through `U(2)` cells with Taylor bounds and separating inequalities, Root certificates and phase balls, the atlas
refinement, exclusion of two further bases through the common Hadamard displacement, closed vector balls and an exact
clique search) are mathematical proofs; the certificate is only that the finite traversal terminated with zero unresolved
records (seed cover 151,744,086 bin pairs → 249,878 centers; 3,083 → 600 → 22,379 records; 163 initially unresolved,
all closed by a fallback; final C shards 282/282/282/281 with 0 unresolved). Its own caveats: the hardware/compiler
semantics "remain assumptions"; an earlier 7 September execution "has not been authenticated"; every rejection uses a
necessary condition with explicit margin, so caps never exclude.

**Hooks, typed.**

- **IDENTITY (number ring).** Its Lemma 1.2 builds the three bases as `A_r ⊗ B_r` on `C^2 ⊗ C^3` with the qutrit bases
  `I, F, DF` (`F` the `3×3` Fourier matrix, `D = diag(1, ω, ω)`); the qutrit overlaps are `(F*DF)_{jk} = (1+2ω)/3` or
  `(2+ω²)/3`, squared modulus `1/3`. The element `1+2ω` with `(1+2ω)² = −3` is exactly the skew operator of the directed
  3-cycle on the centred triangular lattice ([THM-4813 (B)](../../01-canon/theorems/THM-4813-collatz-numerator-coordinate-forcing-cocycle-gauss-sum-intertwiner-and-sharp-spectral-exception-budgets.md):
  Paley skew = Gauss-sum multiplication; [THM-4814 (E)](../../01-canon/theorems/THM-4814-order-ten-prime-cycle-index-types-order-eleven-ledger-skew-hadamard-nuclear-norm-extremality-and-black-core-distance-of-regular-tournaments.md):
  the root lattice `A_{q−1}` of a doubly regular tournament is a `Z[(1+√−q)/2]`-module through `(I+S)/2`). The qutrit
  unbiasedness `|1+2ω|² = 3` and the Paley identity `S² = J − 3I` on `Z/3` are the same equation `g² = −p` for the
  quadratic Gauss sum, `p = 3`. In every prime dimension `p ≡ 3 (mod 4)` the standard (Ivanović / Wootters–Fields)
  overlaps are quadratic Gauss sums `g/p` with `g² = −p`, the equation satisfied by the Paley skew matrix `S_p`.
  Typed IDENTITY for the algebra; ANALOGY for the objects (bases are not tournaments).
- **COINCIDENCE (no tournament shadow of the dimension-6 obstruction).** Six is the first non-prime-power dimension.
  On the tournament side there is no obstruction at six: the nuclear-norm maximum `6√5` at order 6 (THM-4814 (C)) is
  attained, by the skew-conference tournament `S Sᵀ = 5I` (Paley conference matrix of order 6), and skew-Hadamard
  tournaments exist at order 7. Do not promote "6 is special" across the analogy.
- **TOOL.** Vanishing sums of six unit exponentials (its §2.3, Lemma 2.7 with the bound `‖c + JΔ‖ ≤ 5πw²/4`, Lemma 6.5
  `|Σ e(z_j)| ≥ 6 − 2π√6‖z‖`) are the same object as the roots-of-unity coverage frame of
  [LEM-020](../../01-canon/theorems/LEM-020-coverage-spectrum-frame.md); only the tool is shared.
- **TANGENT.** Corollary 1.3 reads "unbiased" as Killing-orthogonality (`12 tr(XY)` on `sl_6`) of Cartan subalgebras
  spanned by the projectors `bb* − I/6`; the repository's remark that every strong tournament generates `sl_n`
  (metagraph_hall_bott §6) concerns matrix units graded by arcs, a different grading. Nothing transfers.
- **Not in the paper:** real or skew-Hadamard matrices, Paley constructions, tournaments, circulants, Gröbner bases, SAT,
  interval-arithmetic libraries (it uses closed bins with explicit one-sided binary64 margins).
- **Prior read in this repository:** [paper_reset2_transfers_20261007.md](paper_reset2_transfers_20261007.md) used the
  paper's finite-recursion interface ("a processed node must retain all its unproved children") as a template for the
  Mersenne reset-two ledger; it did not read the mathematics. Section 1 above is the first mathematical merge.

## 2. Local smoothing for the Euclidean wave equation (Kakeya and Nikodym in every dimension)

**Claims.** Theorem 1.1: for `n ≥ 2`, `2 ≤ p < ∞`, `ε > 0`, `‖U(T)f‖_{L^p(R^n×[1,2])} ≤ C ‖J^{s_+(p)+ε} f‖_{L^p}` with
`s_+(p) = max{0, (n−1)/2 − n/p}` (Sogge's conjecture; `ε = 0` at `p = 2`). Engine: the critical annular estimate at
`p_c = 2n/(n−1)` with loss `[log(2+N)]^{2/p_c}` (Theorem 1.2). Corollary 9.1: the spherical extension estimate on
`B(z,R)` with `[log(2+R)]^{2/p_c}`, the **Kakeya maximal estimate** `‖K_δ f‖_{L^n(S^{n−1})} ≤ C_n [log(2/δ)]^{2−2/n} ‖f‖_{L^n(R^n)}`,
and **every set containing a unit segment in each direction has Hausdorff dimension `n`**; (9.7) the Nikodym maximal
estimate and dimension statement; Theorems 10.1/11.1 Bochner–Riesz and Stein's square function; Corollary 12.1
paraboloid restriction and Schrödinger local smoothing. The proof is analytic throughout: a conserved current for a
superposition of waves, quaternionic coordinates adapted to the wave, a Moore-cofactor determinant identity with
divergence-free rows, a **determinant energy** built from a quaternionic Monge–Ampère Dirichlet potential, and a
thin-cylinder supersolution argument costing one logarithm per passage. Its own caveats: `ε`-loss for `p > 2`; the
optimal log power is open between `1/p_c` and `2/p_c`; the Kakeya log power `2 − 2/n` versus the necessary `1 − 1/n`;
restriction only in the interior of the range; Nikodym on manifolds only locally with constant curvature. Its references
[62]–[64] are three further OpenAI preprints (4D Kakeya dimension, 3D Kakeya maximal, elliptic-capacity sphere
restriction), not read here.

**What it changes here, typed.**

- The reference file [CORE-PAPERS-KAKEYA-4D-2026-08-24.md](../reference/CORE-PAPERS-KAKEYA-4D-2026-08-24.md) kept the
  Euclidean linear Kakeya conjecture in `R^4` **OPEN** (Katz–Zahl `3.0598…`, Borges et al. preprint). A dated block now
  records the external claim: every Euclidean Kakeya set in every dimension has full Hausdorff dimension, with the
  maximal estimate above; **UNVERIFIED here**; the file's warning that multilinear, sticky, plany, finite-field and
  Arithmetic-Kakeya statements must not be merged stands.
- **NO TRANSFER to the repository's finite objects.** The 4D program's theorems are finite-field and combinatorial
  ([THM-4235](../../01-canon/theorems/THM-4235-finite-four-dimensional-torus-kakeya-quadratic-carrier.md) quadratic
  carrier over `F_61`, [THM-4236](../../01-canon/theorems/THM-4236-four-dimensional-kakeya-monic-cubic-focal-selector.md)
  focal selector), the Kakeya numbers `K(A_5) = 15`, `K(2I) = 30`
  ([THM-870](../../01-canon/theorems/THM-870-kakeya-numbers-icosahedral.md)), the `F_p^2` one-line-per-direction minima
  `7, 17, 31` ([HYP-2235](../hypotheses/HYP-2235-finite-field-kakeya-falconer-carrier.md)) and the LRC(14) "needle"
  and "comb" theorems named Kakeya by analogy (THM-1128, THM-1154, THM-1236, …). The paper contains no finite-field,
  polynomial-method, decoupling or incidence content (its digest checked Dvir, Bourgain–Guth, polynomial partitioning:
  absent). None of those objects is touched, in either direction.
- **TOOL RESEMBLANCE only.** The determinant energy's comparison `(n−1)/(n+1) X_0 ≤ F ≤ n/(n+1) X_0` and the Moore
  cofactor identity `det(A + t vv*) = det A + t v*C(A)v` are the same elementary facts the repository uses for
  Hadamard-type and nuclear-norm bounds (THM-4814 (C)); nothing beyond the identities is shared.

## 3. "The Exponential Time Hypothesis is false"

**Claims.** Theorem 1.1: for each fixed `k ≥ 2` a deterministic algorithm counts the satisfying assignments of a
`k`-CNF on `n ≥ 20` variables and `m` clauses in `(n+m)^{O_k(1)} 2^{O_k(n log log n / √log n)}` bit operations and
polynomial space (for `m ≤ Dn`: `2^{O_{k,D}(n/√log n)}`); a rational-interpolation route gives exponent
`ρ(n) = n (log log n)² / √log n`. Consequences stated: ETH, SETH, Gap-ETH, NSETH and their randomized forms refuted;
`#ETH` refuted; subexponential counting of independent sets, proper `q`-colourings (`2^{o(v)}`) and perfect matchings
(`2^{o(v+e)}`); Corollary 1.2 (Max-`k`-SAT within `εm` in randomized subexponential time); Corollary 1.4 / C.1 circuit
lower bounds `C_L(n) > c n √log n` for some `L ∈ E^NP`; Corollary C.2 inversion of constant-locality maps. Technique:
clause-respecting partitions of the assignment cube, **shared linear summaries** (an `(r, τ, D)`-design with
`r = Θ(√log n /(log log n)²)` groups and `≤ 2τn` integer measurement rows), **low-degree replacement forms** built from
swap matrices and a rigidity-style `(−1)^j` interpolant, Hermite/Taylor recovery of a univariate `H(T)` of degree
`≤ 3nR`, counting common zeros of sparse `F_2` systems by Beigel–Tarui modulus amplification mod `2^ℓ`, a counting-
preserving sparsification, exact fraction-free arithmetic, and (Appendix D) a tensor-contraction route with 2-adically
tapered Hermite recovery. Its own caveats: `2^{o(n)}` bounds do not transfer under general polynomial-time reductions,
so no general NP problem is shown subexponential; OVC-based and unrestricted-width CNF-SAT lower bounds stand; the
circuit bound is not a superpolynomial bound for NP.

**What it changes here, typed.**

- **NO DEPENDENCE.** No canon or result file uses ETH or SETH as a hypothesis (grep over the knowledge, canon and
  reflection trees: zero occurrences of either). The problem classification by arithmetical level and busy-beaver size
  ([THM-4622](../../01-canon/theorems/THM-4622-logical-status-dictionary-absoluteness-choice-and-automorphisms-for-jc2-collatz-hadwiger-nelson.md),
  [procgen_atlas_20260924](procgen_atlas_20260924_collatz_implication_atlas.md),
  [analogy_reduction_atlas_20261025](analogy_reduction_atlas_20261025/README.md)) is logical, not fine-grained; it is
  unaffected either way.
- **NO TRANSFER to `H(T)`.** The repository's Hamiltonian-path count `H(T)` and its symmetry quotient `tc = H/|Aut|`
  (THM-1370, THM-4812) are computed by subset dynamic programming in `O(2^n n^2)`. A `k`-CNF encoding of Hamiltonian
  paths needs `Θ(n^2)` variables (position × vertex), so Theorem 1.1 gives `2^{o(n^2)}`, weaker than the DP; Corollary
  1.3's graph counts (independent sets, colourings, matchings) have `O(v)` variables and are the genuine beneficiaries.
- **TOOL RESEMBLANCE only.** Modulus amplification mod `2^ℓ`, 2-adic tapering `q_t = h·max(0, u − v_2(t))` and Hermite
  recovery over `Z/2^ℓ` resemble the repository's 2-adic clocks (THM-4556 and the Collatz numerator coordinate THM-4813);
  the resemblance is in the ring, not in any statement.
- The machine–Collatz dictionary (a `k`-state machine is a modular loop with `O(k)` branches), the three levels of
  invariant (regular, semilinear, beyond) and the two tool identities (the replacement interpolant as a Bonferroni
  truncation; modulus amplification as the Hensel idempotent) are in the reflection
  [eth-false-busy-beaver-sizes-and-modular-loops](../../07-reflections/eth-false-busy-beaver-sizes-and-modular-loops-the-machine-collatz-dictionary-klein-20261027.md).
- The paper cites "the Unique Games Theorem [Ope26]" (another preprint of the same collection) only to contrast
  NP-hardness of approximation with its own subexponential approximation; that claim is not examined here.

## 4. The QCD / root-lattice analogy merged with the dipole law (prior results on "quantum chromodynamics")

The repository's QCD content is one typed analogy (metagraph_hall_bott §6; five_port §5): an arc `i → j` maps to the
root `e_i − e_j` of `A_{n−1}`; the sum over arcs is the score-deviation vector `2s − (n−1)`; matrix units
`[E_ij, E_jk] = E_ik` retain composable paths but every strong tournament generates all of `sl_n`, so the grading must be
kept. SU(3) enters only as the Lie algebra whose root system `A_2` is the Eisenstein lattice; quark colours and gluons
are representation dimensions, not graph colours (the note says so). What the dipole law
([THM-4811 (A)](../../01-canon/theorems/THM-4811-arc-strong-connectivity-is-the-exact-complement-distance-to-the-nonstrong-locus-black-refinement-and-the-dipole-law.md))
adds, PROVED and elementary:

- **The dipole is a root.** A native move along `Q = (a, …, b)` sends the deviation vector `d = s − (n−1)/2 ∈ (1/2)Z^n`
  (sum zero, so a point of the `A_{n−1}` weight space; integral for odd `n`) to `−d + α`, `α = e_a − e_b` a root. A word
  of length `ℓ` gives `d_ℓ = (−1)^ℓ d_0 + Σ_i (−1)^{ℓ−i} α_i`, a signed sum of `ℓ` roots.
- **`L` is the root-word metric.** The minimal number of roots whose signed sum is a sum-zero integer vector `v` is
  `‖v‖_1/2`; minimising over relabellings is the sorted matching. So for any target `T'`,
  `d(T, T') ≥ min(ceil_even (1/2)Σ|s_(i)(T) − s_(i)(T')|, ceil_odd (1/2)Σ|s_(i)(T^op) − s_(i)(T')|)` (the parity of `ℓ`
  decides which of `±d_0` is compared); for `T' = TT_n` both terms are THM-4810's `L = (1/2)Σ|s_(i) − i|`. The score
  bound is thus the `A_{n−1}` root-word distance between sorted score vectors, which is the exact content of
  "arcs are roots".
- **The `A_2` case.** For `n = 3` the roots are the six Eisenstein units; the directed 3-cycle's skew operator is
  `1+2ω` (THM-4813 (B)), the qutrit MUB overlap of section 1; the three ROOT rows of the owner's directive embed in
  the same lattice (THM-4813 (B): `−RK`, the ideal `(√−3)` times a unit).
- **Typing.** IDENTITY for the lattice algebra (arcs ↦ roots, dipole ↦ root, `L` ↦ root-word metric); ANALOGY for the
  physics, exactly as before. No confinement, Wilson-loop or gauge-theoretic statement of the repository is touched
  (the only other QCD mentions are a six-gluon-amplitude Zenodo deposit read for its symbols, and a critical review of a
  "transformers are lattice gauge theories" manuscript).

## 5. What was and was not done

Done: three complete page-cited digests (readers in parallel, no code from the papers executed); the hooks above;
the dated external-claim block in the 4D Kakeya reference file; the companion packet's census and ledger. Not done:
any re-derivation or audit of the three papers; reading their cited sibling preprints (Unique Games; 4D/3D Kakeya;
elliptic-capacity restriction); any change to the status of a repository theorem on their account — external
preprints do not close repository OPEN markers, they annotate them.
