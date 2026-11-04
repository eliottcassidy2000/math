# Session 179 synthesis: exact arithmetic state and the limits of finite shadows

**Status.** Five independently audited research packets, two elementary interface corollaries independently reviewed in [the synthesis audit](AUDIT_CROSSFIELD.md), and explicit next proof obligations. No major open problem is claimed solved, and historical priority is not established. The shared advance is to retain the arithmetic information needed to verify a claim, rather than infer an exact theorem from a finite quotient, a scalar resemblance, or a local rank.

## 1. The five deliverables and their separate scopes

1. **[Periodic Collatz template atlas](../collatz/RESULT.md).** An arbitrary positive valuation word has a rational cycle center and a phase-corrected proper integer energy. For a noncenter source, one valuation of its displacement from the center gives the exact first mismatch against the periodic word, without orbit search. Matching steps preserve energy while strictly decreasing precision. This extends the older valuation-two rank; changing templates is not proved to decrease any common rank. [Independent audit](../collatz/AUDIT_REPORT.md).

2. **[Collatz grammar boundary and arithmetic certificates](../algorithms/RESEARCH.md).** Every context-free sublanguage of the sharp slope-contracting language has a uniform exponential contraction margin. Consequently the complete sharp language is not context-free and its count set is not Presburger-definable. Exact affine triples, parity legality, and a strict height threshold nevertheless give compositional certificates for individual descents. This is not a lower bound on proving Collatz: a conservative context-free family could still cover every integer. [Independent audit](../algorithms/AUDIT.md).

3. **[Crossfield transfer and valuation duals](../crossfield/REPORT.md).** Prime valuations turn divisor-sum product identities into exact integer-column equations. A replayable separator excludes a finite catalogue of 30 blocks even under arbitrary nonnegative multiplicities. Separate inherited max-product and transfer-operator controls distinguish pathwise descent from averaged spectral contraction. Neither a finite catalogue nor a selected path language gives universal integer coverage. [Independent audit covering this packet](../algorithms_audit/AUDIT.md).

4. **[Primitive-order domains for odd perfect supports](../oddperfect/PROOF.md).** Classical Bang–Zsigmondy, with every exception checked against Euler parity, forces each exponent length to be one of at most s−1 in-support multiplicative orders and demands all its nontrivial divisor levels. Exact closed rows and column balance give a complete search of at most (s−1)^s row combinations for a specified support and Euler prime. There is no global bound on the support. [Independent audit](../oddperfect/AUDIT_REPORT.md).

5. **[One-apex tournament modular synthesis](../spectrum/REPORT.md).** An exact three-register append interface yields, for every odd q and every residue r, a strong one-apex tournament with H(T)≡r mod q and at most 3q² vertices. This is an all-order modular realization theorem, not exact realization of arbitrary odd counts. Permanent exact gaps such as 7 and 21 are compatible with it. The packet's direct subset-DP controls independently check its interface. [Separate independent audit](../spectrum/AUDIT_REPORT.md).

The five interfaces retain different information. None substitutes for another simply because each can be described using words, flows, ranks, or prime factors.

## 2. What is inherited, and what this session actually adds

### Inherited facts and mechanisms

- Euler parity, multiplicativity of σ, multiplicative orders, LTE, and Bang–Zsigmondy are classical inputs.
- Semilinearity of context-free Parikh images is classical. The algorithms packet explicitly cites that input; the strict uniform-gap deduction is proved there.
- Exact affine branch composition, finite max-cycle/potential duality, and linear separation are established mechanisms. The crossfield packet identifies where the repository already used them and does not rename them as discoveries.
- The old Collatz rank `3^K t²` at n−1=2^K t is explicitly recovered, not rediscovered.
- Tournament SCC factorization, endpoint interfaces, one-apex pair counts, and the exact global gaps 7 and 21 have prior repository authorities.

### New deductions or newly delivered integrations

- The Collatz phase correction works for every positive valuation word, including centers with negative D, and gives an exact mismatch formula and proper local rank. Repeated words are proved redundant within a starting phase.
- The algorithms packet proves the uniform-rate obstruction for every context-free subset of the sharp slope language and provides an authenticated strict-descent checker, with an explicit first positive certified source.
- The crossfield packet supplies a specific exact 23-coordinate separator for a 30-block multiplicative catalogue, retaining outside primes and multiplicities. The method is inherited; the displayed certificate and its hostile controls are the delivered work.
- The odd-perfect packet imports primitive-divisor existence into the S177/S178 order/valuation architecture, replacing a dense exponent interval by sparse divisor-complete candidate domains. This is a classical-theorem corollary and a new integration relative to those packets, not a historical novelty claim.
- The tournament packet proves transitivity of the append-gate action modulo every odd q, including multiples of 3, and turns it into strong tournaments with the explicit quadratic order bound.

## 3. Two rigorous interfaces between packets

These corollaries explain how the results can be used together without silently strengthening their conclusions.

### 3.1 Compile a template match directly into a descent certificate

Fix a Collatz valuation template w=(a_0,...,a_(L−1)) and a positive odd source n which is not its exact center. The atlas gives the exact number J of matching odd steps. A valuation letter a is the shortcut parity block `1` followed by a−1 zeros. Therefore the J-letter matched prefix corresponds to a legal shortcut word of length S_J, and its affine triple is

    (A,B,C)=(3^J, 2^(S_J), C).

**Corollary.** This matched prefix gives a strict descent if and only if

    (2^(S_J)−3^J)n > C.

Its triple can be computed by powering the full-template triple and appending the remaining prefix, without evaluating all intermediate trajectory states.

**Proof.** The mismatch theorem authenticates every valuation letter before J. Translating each accelerated odd step into its shortcut block authenticates the shortcut word. Exact affine composition gives its endpoint (3^J n+C)/2^(S_J). Comparing this positive endpoint to n proves the criterion. Binary powering of the composition operation gives the stated compiler. If J=0, the empty triple gives equality, so no descent is falsely certified. QED.

This is a concrete local certificate pipeline. A template may mismatch before providing descent. Choosing another template still requires new work: no common decreasing rank across template changes follows from this corollary. Arithmetic values can have large bit length despite a short expression.

### 3.2 A one-row-per-prime separator after primitive-domain reduction

Fix support S and Euler prime q. Compute the primitive-complete domains and retain only rows whose σ value has no outside factor and the correct 2-adic valuation. If any source has no retained row, reject immediately. Otherwise each source p has a finite nonempty set R_p. A row corresponding to length n has odd-prime balance vector

    d_(p,n)[r] = v_r(σ(p^(n−1))) − (n−1)[r=p], r∈S.

An odd perfect number in this branch requires one row per source with sum of these vectors equal to zero.

**Corollary (sound normalized separation).** For any rational vector y indexed by S, if

    sum_(p∈S) min_(n∈R_p) y·d_(p,n) > 0,

the branch contains no odd perfect number.

**Proof.** Every selection of one row per source has y-pairing at least the displayed positive sum. Its vector sum therefore cannot be zero. QED.

This combines the sparse-domain theorem with the crossfield dual mechanism while retaining the one-row-per-prime constraint. It is a different relaxation from the crossfield homogeneous cone, which allows arbitrary multiplicities of catalogue blocks. Neither relaxation is asserted always to separate an integer-infeasible branch. This corollary is not accompanied by a new numerical separator or a claimed benchmark improvement.

It suggests a rigorous next experiment: build the sparse closed-row catalogue, propagate S178 capacity/cut bounds, and seek exact rational normalized separators on the remaining branches. Every proposed separator must be replayed against **every** remaining row with exact arithmetic.

## 4. Twelve precise older source leads

Paths are relative to the repository root `math/`. The packet reports provide their fuller surveys. The entries below name concrete operations recovered, not endorsements of all prose in the sources.

1. `05-knowledge/results/collatz_branch_toll_rank_20261004.md`, §4 and §6: root-centered `3^K t²` rank and the larger-child family 223+62208t → dependency 233+65536t. Imported into the arbitrary-template rank; neutral inverse towers add no coverage.
2. `05-knowledge/results/collatz_run_block_family_20261005.md`, §2: exact inverse phases for maximal valuation-one runs. Generalized to periodic-template parsing, without importing its strong-reset coverage beyond its guard.
3. `05-knowledge/results/collatz_weak_reset_family_20261005.md`, Proposition R: positive rational cycle centers behind weak first descents. Supplies the correct center object, not a proof that all trajectories eventually leave every obstruction.
4. `01-canon/theorems/THM-2211-carry-regime-root-transducer-and-infinite-autonomous-index.md`: an external unbounded carry stream is a resource, not a finite autonomous state. Explains why three scalar registers or an affine triple need not have bounded information content.
5. `05-knowledge/results/crossroads223_20260926_automata.md`, §§3–5: complete finite residue lifts preserve cycle-mean and weighted-growth data. More residue fibers alone do not remove the expanding language obstruction.
6. `01-canon/theorems/THM-4594-the-maximal-class-decided-collatz-sieve-and-the-sign-barrier.md`, §4: positive cycles can occupy slope-contraction thresholds. The new checker retains strict height inequality instead of checking the slope alone.
7. `01-canon/theorems/THM-3078-z234-direct-farkas-normalization-and-four-two-high-boundary.md`: explicit rational duals and independent certificate replay. Transferred to exact divisor-sum valuation columns.
8. `05-knowledge/results/procgen_tension_20260926_ranks_christoffel_duality.md`, Theorems R and D: cycle/potential alternatives and optimal rank defects. This lane is inherited; the new work does not present max-plus vocabulary as a solution to Collatz.
9. `01-canon/theorems/THM-4640-the-repair-algebra-at-minus-one-collisions-are-increasing-two-three-representations-ladders-and-gersonides-carries.md`, §5(a): a Bang–Zsigmondy primitive-divisor argument in the October 14 Wilson/Collatz repair algebra. The recovered primitive-factor tool yields sparse odd-perfect exponent domains.
10. `01-canon/theorems/THM-4146-rational-three-cycle-order-six-lift-horizontal-divisor-fibre-firewall.md`, order-six finite-ring discussion: exact composite-modulus period can occur without a new primitive prime. Retained as a hostile control, not confused with the primitive-prime theorem.
11. `07-reflections/the-apex-recursion-modular-tree-recursive-H-and-recursive-carving-s531.md`: nested reverse-arc blocks couple instead of simply multiplying. The append gates retain an explicit state for one tractable nested family.
12. `01-canon/theorems/THM-3134-tournament-endpoint-jet-and-c3-newton-profile-transform.md`: general substitution needs a full endpoint interface. The small one-apex interface is operation-specific and cannot be generalized to arbitrary ears or substitutions without extra data.

Further verified boundaries worth keeping beside these leads are THM-4507 (finite valuation corrections do not provide a global Collatz logarithmic potential), THM-4094 (exact tournament prime lanes), THM-1370 (global gaps 7 and 21), and the S177/S178 support and valuation packets.

## 5. Actual next proof priorities

### Priority A: arithmetic catalogue integration

The most bounded and immediately testable continuation is the odd-perfect pipeline in §3.2. First exhibit a support/branch where sparse primitive domains and normalized separation save work beyond either original filter. Then ask whether an infinite support family admits a symbolic domain or dual certificate. Stop short of claiming global progress until there is an argument controlling support, not just exponents on a supplied support. The concrete five-prime support already rejected in S178 is a mechanism control, not a new number-theoretic exclusion record.

### Priority B: cross-template Collatz descent, not more local ranks

Use the parser and authenticated affine checker to search for a precisely specified mismatch transition between two primitive templates. The target theorem must provide a common proper rank which strictly decreases at that transition and prove every intended source satisfies its arithmetic guards. A finite bank of observed transitions is not a coverage proof. Keep cycle-center cases explicit; a center is not automatically the terminal 1.

A distinct viable formulation is a conservative context-free certificate family with a provable source-covering argument. The grammar theorem does not rule that out. It rules out representing the entire sharp slope language by such a grammar. Establishing conservative coverage remains equivalent in strength to the relevant Collatz descent assertion, so renaming the family is not progress by itself.

### Priority C: bounded exact lifts for tournament counts

The modular state graph is already connected for every odd modulus. The unresolved task is to control the exact count, or its excess above a target, while lifting a residue witness. A useful theorem must bound this excess tightly enough to address an exact ordinary-prime or exceptional-7p lane. More modulus censuses cannot exclude an exact odd value, because this construction realizes every odd-modulus residue. No bounded-lift claim is made here.

### Avoid spending another session on these false upgrades

- Adding powers of one template does not create a new cycle center or phase-zero ordering.
- Neutral inverse towers yield new child heights but not new source grounding.
- Primitive-level Hall matching is redundant once all corresponding edges are forced; missing level existence is the useful condition.
- Averaged operator contraction does not imply pathwise descent.
- Finite modular reachability does not imply exact integer realization.
- A general undecidable programming class does not establish undecidability of the fixed Collatz map.

## 6. Verification map and scope ledger

- `../collatz/periodic_template_rank.py`: 1,364 words, 51,702 nonfixed controls and five fixed controls. Independent exact-precision audit: 3,060 cylinders and long boundaries through 820 matched letters. Normal and optimized outputs match.
- `../algorithms/check_certificates.py`: 16,382 binary words, 49,146 source lifts, all 13,522 contracting words in the declared universe, and first-descent/conservative controls through 10,000. The proof, not this finite range, establishes the grammar result.
- `../crossfield/valuation_dual.py`: exact replay of all 30 catalogue columns, hostile `7^3` control, and positive controls from 6 and 28. `../crossfield/verify.py` checks its affine-state and averaged-spectrum controls.
- `../oddperfect/check_primitive_domains.py`: 189 support/Euler cases, 124 directly closed rows retained, and 129 primitive levels. Its independent audit confirms the arithmetic and the external theorem's exceptions.
- `../spectrum/modular_apex.py`: all 1,023 binary words through length nine against independent Hamiltonian-path DP, every odd modulus through 101, and the affine group identities in the declared smaller range. The all-modulus theorem is algebraic, not extrapolated from those checks.

The two corollaries in §3 are elementary deductions with proofs, newly assembled here; their cross-packet implementation and any efficiency claim remain future work. The independent synthesis audit approved both interface corollaries and their scope boundaries.
