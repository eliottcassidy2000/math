# Session138 verification and checkpoint history

## Results and provenance

**All current headline claims are independently accepted in their stated scopes.** The all-order asymptotic coefficient is proved; exact finite-order coefficients are proved through eight. Arbitrary strong-quotient scalar zero-rigidity and unrestricted spectrum completeness remain open.

| Checkpoint | Mathematical scope | Evidence retained |
|---|---|---|
| First, published at `6fff918603cec8a2b9e233bbb627d703214d6953` | Multistar coefficient asymptotic to 4/n, matching-support sharpness, exact unrestricted coefficients through seven, reflector/resource theorems, restricted profile rigidity | [Frozen scope](CHECKPOINT1.md), original lane selectors and audit pins, [first portable replay](REPLAY_CHECKPOINT1.json); 118 published blobs verified |
| Second continuation | Sharp order eight; `K_n=8/n+O(n^(-2))` and finite-order/near-candidate bounds; exact all-parameter sign chambers for three actual modules | Separately curated source files, independent proofs, exact receipts, and continuation replay |

The first-checkpoint lane bytes are not rewritten to claim later results. The first source-integrity and replay files retain their original counts. Updated root summaries describe the combined release; historical through-seven statements remain explicitly scoped.

## 1. All-order asymptotic theorem: proof and exact controls

[Main theorem](asymptotic_stability/PROOF.md) · [Independent audit](independent/asymptotic_audit/AUDIT.md) · [Reviewed source pins](independent/asymptotic_audit/REVIEWED_SOURCES.json)

For `n>=4`, the accepted bound is

```text
kappa_n/[1+n*kappa_n^2/(2*c0)] <= K_n <= kappa_n,
c0=(1-exp(-4))/4,
K_n=8/n+O(n^(-2)).
```

The audit checks Pfaffian normalization, one-edge completed-square bounds, independent endpoint interpolation, finite-interval integration, weighted greedy matching, and the retained binary-baseline square completion. It verifies the explicit small-variance and large-variance branches uniformly, without assuming an extremizer exists. The empty determinant remainder and prohibited positive order-two gap extension are retained as boundaries.

The self-contained independent `check.py` imports no producer kernel. Its exact Fraction/Hamilton-DP controls comprise:

- 796 matrices: all ternary-edge matrices through order four, plus fixed-seed rational cases and a sharp fair-matching control at each order five through eight;
- 4,566 one-edge quadratic/interpolation checks;
- 796 full-margin product bounds;
- 1,584 rational telescoping gap checks in two edge orders;
- 765 baseline/triangle/wedge/matching/square-completion checks, including ten near-candidate controls;
- 384 arbitrary nonnegative rational weighted graphs and 2,467 greedy groups.

Normal and optimized receipts match. Exponential integration and asymptotics are mathematical arguments, not floating-point tests. The near-candidate conclusions concern expected binary gap, probability, and variance mass; no finite-n all-supported-transitive or all-pairs-fair classification is claimed.

The [resolvent cross-check](stability/RESOLVENT_ROBUSTNESS_CROSSCHECK.md) supplies a second analytic robustness mechanism with weaker constants: `u-|v|>=exp(-12V)` and an integrated nonlinear gap. The same baseline-preserving matching argument yields the same asymptotic rate, while the sharper displayed finite-order bound uses the Pfaffian constant four. This is a second mechanism, not merely another execution of the first program.

## 2. Sharp order eight: the final zero-wedge certificate

[Theorem](stability/SHARP_EIGHT_GRAPH_CERTIFICATE.md) · [Cycle-reduction audit](independent/ORDER_EIGHT_CONTINUATION_AUDIT.md) · [Independent graph audit](independent/n8_audit/FINAL_GRAPH_AUDIT.md)

The accepted graph inequality is

```text
41m-32m2+8m3-4m4+8p21 >= 0.
```

It has **no extra wedge term**. Earlier `+2w` variants are not substituted for this stronger load-bearing certificate.

The independent graph implementation constructs seven-vertex isomorphism orbits by adjacent-transposition BFS, rather than importing producer representatives or a graph atlas. It covers all 2,097,152 labeled seven-vertex graphs in 1,044 orbits. Adding all 128 neighborhoods for a new vertex gives 133,632 augmentations and covers every eight-vertex graph up to isomorphism. This is not a literal scan of all `2^28` eight-vertex masks.

Matching counts use vertex-subset dynamic programming; wedge-plus-edge patterns are counted wedge-first. These differ from the producer's permutation-orbit construction and pattern-mask counting. The exact minimum is zero, maximum 10,808, with two zero augmentations: empty and perfect matching. Normal and UBSan outputs match. All arithmetic is integral; separate-affinity interpolation and the independently reviewed tournament-cycle reduction establish the theorem over the whole real variance cube.

Positive-variance equality in `Delta>=87V/128` follows from the retained positive triangle reserve and H=1 classification, not merely from interpreting a graph zero count. Exact equality beyond order eight remains open.

## 3. All-parameter mixed-profile sign chambers

[Producer theorem](zero_balance_later/NO_ZERO_THEOREM.md) · [Independent proof audit](independent_zero_balance_later/AUDIT.md) · [Source pins](independent_zero_balance_later/SOURCE_PINS.json)

The accepted theorem covers the three specified actual Q4 modules for every positive integer side multiplicity. It also includes arbitrary-side positivity for seed 996 with sharp constant 135618, arbitrary other-side thin positivity for seeds 867/995, and the sufficient `p^2 R` negative chamber. It does not classify arbitrary side tournaments for the latter two seeds.

The proof uses an exact tridiagonal run-word response, a log-concave binomial-profile likelihood-ratio lemma, and a one-crossing comparison. The all-parameter quantifiers do not come from finite scans.

The independent checker imports only the previously frozen independent profile kernel, not producer algorithms. It reconstructs:

- all three actual twelve-vertex profile pairs;
- all 507 entries of the three response kernels by forced run words;
- 1,331 literal three-module contexts;
- 64 positive Stirling side profiles;
- 852 positive log-concave sequences, 5,112 lift checks, and 2,556 nonnegative multiplier checks;
- the non-log-concave hostile `(2,1,2)` and the sharp singleton-side control;
- 12,288 bounded parameter sign controls, explicitly diagnostic rather than all-order premises.

Normal and `-O` receipts are byte-identical, SHA256
`2510edaf1ee4c59cf1d96619ab71c5bea57376a539e464423c0f2e2f8fa5d770`.
The frozen checkpoint-one independent kernel remains unchanged. None of these three seeds gives an actual scalar-zero counterexample in the proved family; arbitrary strong-quotient zero-rigidity stays open.

## 4. First-checkpoint verification retained

The [stability audit](independent/STABILITY_AUDIT.md), [extended through-seven audit](independent/EXTENDED_STABILITY_AUDIT.md), [structural/resource audit](independent/OTHER_LANES_AUDIT.md), and [profile audit](independent/PROFILE_DICHOTOMY_AUDIT.md) keep their original scopes.

- Independent order-five and order-six Bernstein routes cover respectively 24 signed tensors across twelve orientation orbits and 112 signed tensors across 56 orbits. These are exact tensor certificates over the real cube, not weighted-matrix samples. The full order-six optimized audit is not claimed as an additional replay.
- Independent six/seven graph censuses cover all 32,768 and 2,097,152 labeled graphs, with 76 and 421 zero residuals. Different producer and Gray-code kernels support the finite-to-continuous reduction.
- The reflector verifier checks all 1,099 contexts through order five, all 48 minimal one-edge isolators, 162 actual full-profile substitutions, and sixteen all-order constructor-image values through 200. The resource checker tests 10,000 exact core budgets and 1,331 affine-ray periods; universal scope comes from the source-normal-form proof.
- Profile reconstruction checks fifteen cores, seventeen full profiles, 228 block-profile cases, 836 literal cuts, 108 literal cyclic contexts, 669,924 bounded-core case-reduction states, and three actual mixed-profile hostiles. The producer's separate 3,633,840-permutation leaf check is corroborated by a different independent partition kernel. These are not an expanded strong-quotient census.

The first assembled-packet replay passed twelve scientific executions and two builds, with every output matching and all 108 curated lane inputs unchanged. [Original receipt](REPLAY_CHECKPOINT1.json) · [Original source-integrity record](SOURCE_INTEGRITY.json)

## 5. Fresh continuation replay: terminal PASS

The relocated continuation replay passed **ten scientific executions and four builds**, with all ten outputs byte-matching their frozen receipts and no stderr or UBSan diagnostics. It covers the independent asymptotic checker, producer and independent zero-balance checkers, and both final order-eight graph implementations; each Python checker ran normally and under `-O`, and each C++ implementation ran optimized and under UBSan. Total elapsed time was 19.644 seconds.

All **143 curated lane files** in the relocated tree remained unchanged. The 108 first-checkpoint files remain byte-identical, all 115 prior packet paths are preserved, and the continuation adds 35 explicitly selected lane files. The new source-integrity pass checks 27 manifest entries and seventeen reviewed-source pin entries. These are checks of the stated retained scopes, not a new independent proof or an expansion of any census.

[Continuation replay receipt](REPLAY_CHECKPOINT2.json) · [Continuation source integrity](SOURCE_INTEGRITY_CHECKPOINT2.json)

## 6. Curated inventory and replay discipline

First-checkpoint selectors remain unchanged: `stability/PUBLIC_FILES.txt`, `exact_resources/PUBLIC_FILES.txt`, `zero_rigidity/public_manifest.json`, and `independent/FIRST_CHECKPOINT_FILES.json`, plus the dedicated profile audit packet.

Continuation selectors are:

- [N8_PUBLIC_FILES.txt](stability/N8_PUBLIC_FILES.txt): four producer files, plus this selector and `N8_SHA256SUMS`;
- [SECOND_CHECKPOINT_FILES.json](independent/SECOND_CHECKPOINT_FILES.json): eleven independent audit/source/receipt files, plus the selector and `SECOND_CHECKPOINT_SHA256SUMS`;
- the final asymptotic proof and separately pinned resolvent cross-check, added explicitly;
- [zero-balance manifest](zero_balance_later/public_manifest.json): six payloads plus both manifests;
- the six-file `independent_zero_balance_later` audit packet, including its own `SHA256SUMS`.

No executable, cache, superseded graph variant, exploratory parameter scan, or additional exact-equality candidate is included. Retained finite diagnostic controls have their scopes stated above. A repeated program is a same-source replay, not a newly independent proof.

Copy the assembled packet before replaying programs that may create sidecars; retain sibling directories. Python tensor programs require NumPy; the new continuation Python checks use the standard library. C++ graph certificates require a C++17 compiler. All gates remain active under optimization.

```sh
sha256sum -c SHA256SUMS
(cd stability && sha256sum -c N8_SHA256SUMS)
(cd zero_balance_later && sha256sum -c SHA256SUMS)
(cd independent_zero_balance_later && sha256sum -c SHA256SUMS)

W=$(mktemp -d /tmp/s138-replay.XXXXXX)
cp -R . "$W/packet"
cd "$W/packet"
```

Representative continuation routes:

```sh
(cd independent/asymptotic_audit &&
 python3 -B check.py > replay.normal.json &&
 python3 -B -O check.py > replay.optimized.json &&
 cmp replay.normal.json normal.json &&
 cmp replay.optimized.json optimized.json)

(cd independent_zero_balance_later &&
 python3 -B audit_zero_balance.py > replay.normal.json &&
 python3 -B -O audit_zero_balance.py > replay.optimized.json &&
 cmp replay.normal.json normal.json &&
 cmp replay.optimized.json optimized.json)

(cd independent/n8_audit &&
 g++ -std=c++17 -O2 -Wall -Wextra -Werror graph8_no_wedge.cpp -o graph8_replay &&
 ./graph8_replay > replay.normal.json &&
 cmp replay.normal.json no_wedge.normal.json)
```

No formal proof-assistant certification, external-priority claim, efficient exact Hamilton counter, unrestricted spectrum theorem, or general zero-rigidity theorem follows from these receipts.
