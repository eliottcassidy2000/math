# Session29 synthesis: interfaces, spectra, and what symmetry actually preserves

This session began with `3^k+1`, Collatz resets, and a proposed Dehn-invariant analogy for rank88 Gaussian Hessians. It developed several exact bridges, and rejected several tempting ones. The results below distinguish proved all-order statements, finite exact censuses, and still-open extrapolations. All linked packets include proofs, executable controls, and separate audit records.

## 1. Main outcomes

1. **Every positive odd integer is a blue-fibre count.** An explicit binary construction realizes every odd h as the number of mirror-fixed spine tiles of a strong self-converse tournament, at odd order at most `2 ceil(log2(h+1))+1`. For h>=3 it also gives a strong even-order realization at `2 ceil(log2(h+1))`. Blue7 has an order6 witness with ordinary Hamiltonian count37; blue21 has an order10 witness with ordinary count913. This does not fill the forbidden ordinary Hamiltonian counts7 or21.
2. **An infinite metagraph boundary has exact wiring addresses.** Intrinsic one-outgoing-arc cuts reconstruct each seamed target from rooted blocks. The alphabet has13 outer rooted classes,17 inner doubly rooted classes, and5 mirror-fixed inner classes. This yields exact target-class counts, edge multiplicities, colours, automorphism products, and rational Hamiltonian-mass functions for every block length.
3. **All324 selected order10 lines are global nonbridges.** The initial boundary graph had eight leaf bridges. Six disappear in a larger64-letter boundary; the remaining two have explicit4+6 detours. The distinction between a bridge in a restricted picture and a bridge in the full metagraph is substantive, not cosmetic.
4. **The rank88 extremal-theta scalar temperature profile is globally certified.** It has exactly two reciprocal global minima, and the self-dual point is a strict local maximum. Arithmetic transcendence also excludes a self-dual Hessian zero for the specified extremal-lattice family. The required angular eigenvalues remain undetermined; no new rank88 lattice is constructed and no lattice optimality is proved.
5. **The odd-exponent Collatz family is a full finite-pattern encoder, but a poor source budget.** It realizes every finite accelerated valuation word and has density-one finite descent within every retained modular guard. The family remains too sparse in original-source intervals to pay a153N budget. Neither finite descent nor exponent density proves convergence.

## 2. A vocabulary that retains the needed information

- A **spine tile** is a tournament with a fixed labelled Hamiltonian path, modulo the repository's framing convention. Its class fibre has size H(T)/|Aut(T)|.
- A **blue fibre** counts mirror-fixed tiles. It is not H(T), not a fixed-path count before symmetry division, and not necessarily the number of distinct incident lines when loops are present.
- A **seam** is the unique forward arc across a specified target cut. Other cross-arcs point backward.
- A **rooted address** records an exit-rooted first block, doubly rooted internal blocks, and an entry-rooted last block. It retains isomorphism data that a scalar path count forgets.
- A **path-port state** `(H,S,R)` retains total paths, paths beginning at the active portal, and a sum of endpoint and two-path-cover counts.
- A **signed controller** is the tournament on anti-involution pairs that specifies which pair controls each allowed sign. Reversing traversal changes the sign requirement; forgetting that change produces false bounds.
- A **boundary-incidence deficit** is full merged-node tile degree minus the incidence already represented by a chosen boundary family. It counts missing incidences, not missing edges, since an internal unrepresented edge contributes twice.

There are two different bridge notions throughout: a directed strong bridge inside a target tournament, and a bridge in the undirected metagraph. The former exposes intrinsic seams; the latter is an adjacency property of isomorphism-class nodes.

## 3. Exact rooted reconstruction and its counts

For a strong tournament let delta(X) count arcs leaving X. The identity

`delta(A)+delta(B)=delta(A∩B)+delta(A∪B)+|A\B||B\A|`

forces all one-outgoing-arc cuts with both shores of size at least2 to be nested. Crossing cuts contradict strong connectivity; disjoint or covering incomparable cuts have too many cross-arcs. Consequently each prescribed cut size has at most one such cut.

In an r-block target, the cuts at5,10,...,5(r−1) are intrinsic. Every isomorphism preserves their rooted blocks, and converse reverses the address. For r>=2 the number of ordinary target classes is

`169 * 17^(r−2)`.

The generating function for converse-merged target nodes is

`1/2 * [169 z²/(1−17z) + 13 z²(1+5z)/(1−17z²)]`.

This is an exact all-order classification of this boundary's targets, not of every tournament of order5r. It generalizes to any fixed block size at least2 and any independent mirror-closed tile alphabet. Intrinsic cuts can be recognized in O(n³) time using the union of an in- and out-arborescence, followed by SCC tests on candidate strong bridges. No n! target labelling search is needed.

The total outgoing line count is18^r. The blue line count is18^(r/2) for evenr and4*18^((r−1)/2) for oddr. These count genuine lines because the source nodes are nonstrong and every target is strong.

Three distinct rational functions must not be conflated:

- The sum of target H over all source lines, with parallel multiplicity, is generated by
  `(172z+45144z²+50400z³)/(1−874z−43344z²−50400z³)`, minus172z to retain r>=2.
- The sum of H over ordinary target classes, each once, is
  `(93313z²+3359832z³+4057812z⁴)/(1−823z−37140z²−45180z³)`.
- The unweighted target-class count is the address formula above.

At r2 there are169 ordinary target classes,91 merged target nodes, and324 lines. The sum of full merged target degrees is84097; only324 incidences are captured by this family. The remaining83773 incidences are a concrete budget for the missing surroundings, not evidence that the small drawing determines the entire metagraph.

A subtle self-converse family survives without any incoming blue line from this boundary: the fifth fixed inner address consists of the two asymmetric letters755/949. At odd r>=3 it yields `13*17^((r−3)/2)` such target nodes. They still have blue tiles globally. Thus “no blue line visible in a selected family” does not mean “not self-converse.”

See [rooted reconstruction](session29-target-wiring/REPORT.md), [detour completion](session29-bridge-extension/REPORT.md), and the offline [wiring atlas](session29-wiring-atlas/atlas.html). The atlas's data and JavaScript interactions were checked; actual browser rendering was unavailable and is explicitly unverified.

## 4. Reflection gives a real Lorentzian observer

The18 letter transfers satisfy

`K_mirror(t)=J K_t^T J`, where `J=diag(1,−1,1)`.

This is a structure-preserving map: it respects word composition, reflection, and the target Hamiltonian observable. An even mirrored word has count `H²−S²+R²`. These physical port states satisfy H>=S, so the expression is nonnegative despite the indefinite form.

The aggregate transfer has three states. Every shifted3-by-3 Hankel matrix of its count sequence has determinant

`−19416360960000 * 50400^s`

and signature(2,1). The polynomial `q(x)=x²−874x−45144` has square moment `L(q²)=−116640000`. Therefore even the tail of this scalar sequence cannot be a positive Hilbert-space moment sequence. The operator itself has real eigenvalues and can be made self-adjoint in another positive inner product; the obstruction is the specified scalar observable, not the bare pole locations.

The average target H per outgoing line grows as `c*(rho/18)^r`, where rho is the largest root of `x³−874x²−43344x−50400`, approximately921.115395615, and c is approximately0.22830918837. This is not the growth of A000568 or a global Hamiltonian maximum.

See [exact moment corollaries](LORENTZ_MOMENT_COROLLARIES.md).

## 5. Signed pairs, the full blue spectrum, and a faster algorithm

For a fixed anti-involution, exact path counting needs only a used-pair mask and the current signed endpoint. Previous signs can be summed out. This gives O(m²2^m) arithmetic operations and O(m2^m) storage on2m or2m+1 vertices, improving the initial ternary-state implementation. Fixed odd cycles use the same state space with different initial and terminal conditions.

If f_sigma counts fixed Hamiltonian paths and c_sigma is the centralizer order inside Aut(T), then

`blue(T)=f_sigma/c_sigma`.

The quotient is essential. P7 has H189, automorphism order21, centralizer3 and fixed-path count9, hence blue3 among9 total tiles. It is mixed, despite appearing in the separate converse-Gluck exception story.

The signed controller formula consists of internal factors `1−b_vu*b_vw` and an optional fold factor `1+b_vu`. It is an exact finite antiferromagnetic interaction model, not a metaphor imported from physics. Its mean over independent labels is m! for every fixed controller. After choosing local outgoing bipartitions, all remaining star phases optimize independently as `A+sum|B_v|`.

For the binary construction, add a universal controller whose outgoing labels are all epsilon=±1. Internal placements vanish; endpoint placements give

`f_new=f+(1+epsilon)A`, `A_new=2A`.

Starting at f=A=1 gives every odd value from1 to2^m−1, uniquely encoded by the row-sign bits. The controller is transitive and the anti-involution centralizer is trivial. A final positive row makes the construction strong; h1 uses C3 on the odd side. Thus31,63 and127 are exact all-positive-row blue counts, without any assertion that these values are rare or unrestricted in the ordinary H spectrum.

Unused odd-order mate arcs can break the centralizer without changing f. Gluck's odd-order regular-subset theorem therefore shows that the maximum blue count in order2m+1 equals the maximum fixed-path count F_m. This is a genuine use of the earlier group-theoretic thread. [Gluck, Corollary1, p66](https://doi.org/10.4153/CJM-1983-005-2).

See [signed spectrum](session29-coherent-halfpaths/REPORT.md) and [independent audit and faster algorithm](session29-signed-audit/AUDIT.md).

## 6. Average versus maximum: two different asymptotic questions

The classical odd-partition formula counts self-converse classes without enumerating tournaments. If b_n is that count and

`a_m=(m)_3*2^(8−4m)/3`,
`c_m=[(m)_5*2^24/5+(m)_6*2^28/18]*2^(−8m)`,

then the proved expansion is

`b_(2m)=2^(m(m−1))/m! * [1+a_m+c_m+O(m^9*2^(−12m))]`,
`b_(2m+1)=2^(m²)/m! * [1+a_m/4+c_m/16+O(m^9*2^(−12m))]`.

Consequently the mean blue fibre, divided by m!, tends to1 from below. The remainder proof controls all odd partitions, not merely finitely many defect types. The exact census formulas are classical and already present in the repository; their provenance is retained.

The maximum is a separate question. Exact finite fixed-path maxima through seven pairs are

`1,3,9,39,195,1185,8295`.

They match the established sequence

`C_m=m! * sum_(j<=floor(m/2)) 1/(2^j j!)`.

For C_m itself, the EGF is `exp(z²/2)/(1−z)` and C_m/m! tends to sqrt(e). It counts partial matchings with ordered unmatched vertices; equivalently it is a moment of an independent exponential-plus-Gaussian variable. NONE of this proves F_m=C_m for arbitrarym. A valid matching injection or another upper-bound proof, together with an attaining family, is missing. The first unsigned coherence guess was false because reversing an edge reverses its sign requirement. Even transitive controllers cease to maximize the genuine signed problem at five pairs.

These statistics are explicitly defined here; the older unlocated symbol M_k is not silently redefined.

## 7. Collatz: full local language is not global payment

For odd k, `u_k=(3^k+1)/4` is a bijective residue encoder on odd2-adic classes. In fact `v2(u_k−u_l)=v2(k−l)` for distinct odd exponents. Every finite accelerated valuation word occurs infinitely often, with exact geometric cylinder density. The ordinary Collatz clock is the sum of valuations plus the number of accelerated steps.

A uniform carry estimate yields an exact finite-depth endpoint-descent density and an exponentially decaying exceptional bound. The eventual-first-descent set has density1 among odd exponents, separately in every retained mod15/mod35 register. This is a positive theorem, but not convergence: descent may leave the family, and an endpoint-descent probability is not a first-descent CDF.

The contrasting source bound is decisive for the debt problem: an interval `[4N,256N)` contains at most two family sources. Infinite exponent density cannot be relabelled as153N distinct original-source credit. No residual R' bound, uniform coverage, or cross-term debt closure is claimed.

See [finite language](session29-collatz/REPORT.md) and [density-one descent](session29-collatz/phase2/DENSITY_ONE_DESCENT.md).

## 8. Dehn, angular information, and the completed scalar certificate

A harmonic angular tensor is rotation-covariant. Ordinary scissors congruence forgets precisely that frame. An additive equivariant map into the harmonic-quartic space cannot become a nontrivial rotation-invariant Dehn valuation by simply taking a quotient. Spectra also lose relative alignment: two angular operators can have the same spectrum but different spectra after combination.

The valid replacement is a direction-enriched Dehn object on actual polyhedral edges, retaining edge directions and angle classes. A regular tetrahedron supplies a nonzero example. No map from a rank88 shell to such a three-dimensional polyhedron has been supplied.

A concrete lattice test separates E8⊕E8 and D16+ despite identical scalar theta series: angular trace powers agree through degree3 and first separate at degree4. E8⊕E8 has an all-temperature saddle with exact inertia(71+,64−). This is an actual information hierarchy, not an equality inferred from the number24 or a shared modular weight.

For rank88, the certified larger minimizing temperature lies between1.098604471215707 and1.098604471215708; the minimum ratio lies between11135.8763212 and11135.8763213. The reciprocal is the only other global minimum, and1 is the only other stationary temperature. Nesterenko's theorem makes the self-dual scalar ratio transcendental, while the actual angular operator has algebraic eigenvalues, excluding a zero there. That arithmetic nondegeneracy extends to the stated extremal even-unimodular dimensions24m+16. It does not determine the sign of every eigenvalue or prove global energy optimality.

See [Dehn report](session29-dehn/REPORT.md), [rank16 test](session29-rank16/REPORT.md), [rank88 certificate](session29-gaussian/REPORT.md), and [arithmetic extension](session29-gaussian/ARITHMETIC_ADDENDUM.md).

## 9. Highest-value remaining obligations

- Prove or refute the signed maximum/matching envelope beyond seven pairs; do not promote the sqrt(e) extrapolation.
- Classify the13 remaining full64 target-boundary leaves globally, and continue the black-component quotient analysis. Raw multigraph bridges differ from bridges after suppressing parallel lines.
- Extend rooted-address reconstruction past the one-forward-seam setting, especially strongly2-edge-connected targets.
- Pay Collatz residual debt with distinct original sources and actual clocks; finite symbolic universality alone cannot do it.
- Determine the actual rank88 angular minimum eigenvalue and shell data needed for stability.
- Recover the intended definitions of the historical D(z) and M_k before claiming they have been solved.

No proof of Collatz convergence, Hadwiger–Nelson's exact chromatic number, the lonely runner conjecture, RH, the global ordinary Hamiltonian spectrum, or complete tournament metagraph classification is claimed. The progress is in exact constructions, obstructions, refined counts, and algorithms that make those larger questions better posed.
