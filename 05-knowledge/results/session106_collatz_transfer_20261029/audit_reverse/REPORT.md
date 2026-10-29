# Independent reverse/descent audit

**Verdict: PASS within the explicitly stated scopes.** No producer code is imported. Proofs and exact native enumeration were inspected independently; no historical74-return search was rerun. The all-order claims below are checked mathematically, not inferred from the finite receipts.

## Reverse forest

Inspected `reverse/REVERSE_FOREST.md` and `check_reverse.py`.

1. The inverse formula follows by solving each affine branch; source guard and positivity are retained. On81Z, nonzero branches force exact3-adic valuations4,6,5,4 respectively, because their affine parenthesis is a3-adic unit. Their adjusted unit residues61,13,70,34 are distinct. Thus there is at most one nonzero inverse.
2. A zero inverse exists iff v3(y)>=5, preserves adjusted unit chi and lowers valuation by one. At the two possible splits, the zero alternative lasts respectively one or two edges before leaving the allowed valuation/chi combinations. It cannot branch again. Choosing the nonzero edge yields a spine with terminal zero side paths of length<=2. This proves width<=3 at **every depth and endpoint**, without a valuation cutoff.
3. Each later split on the spine occurs after a positive defect increment. Its zero side path retains the defect before that increment. Therefore the spine and distinct contemporaneous side paths have strictly different total D-J. The identity5J-3D=2k is correct for all five branch clocks. Together these prove the stated endpoint/native-clock uniqueness.
4. The argument is scoped to the five-return language on positive81-multiples. It supplies neither a forward escape bound nor payment below an original source. The27-layer collision does not itself disprove all possible width theorems outside81Z; the document correctly limits its conclusion to failure of this proof/interface.

The first inspected producer checker used removable assertions. This was reported and the producer replaced them with explicit RuntimeError checks before finalization. The independent checker always uses explicit exceptions.

**Independent native census:** enumerate each positive81-multiple source through300000; run the literal ternary queue until its first completed return; identify growing-return words; store incoming native edges. For every81-multiple endpoint through300000, compare that dictionary with inverses independently obtained from the full odd-native-word affine formula and literal valuation replay. Because every accepted edge grows, no source above the endpoint cap can be a missing inverse. Then explore every such endpoint to reverse depth15, checking width, distinct defects, and clock recovery. This is a complete bounded endpoint census, not sampling reverse edges.

## Source-paid atlas and rank barriers

Inspected `descent/SOURCE_PAID_ATLAS_AND_STATE_BARRIER.md`, `check_descent_transfer.py`, and `atlas16.json`.

* The final-oddness source guard does imply every prefix's exact valuation: the suffix carry is odd, so reducing the concatenation identity modulo2^(D_j+1) forces the prefix numerator to have valuation exactly D_j. Prefix quotients are therefore authenticated native endpoints.
* Payment at the smallest nonterminal cylinder member plus B>A proves payment for the entire positive ray. Slope alone is not used. The checker propagates both constant and coefficient to verify exact divisibility and oddness for every parameter value.
* The whole-cylinder first-payment decision is prefix-dependent but independent of a sampled representative; accepted cylinders therefore form a prefix-free family. The finite covered/unresolved partition is not an all-integer convergence statement.
* The finite-sidecar lower bound is correct: k unpaid all-one returns require k+1 different sidecar states at a fixed original source. Choosing k=floor(log2(N+1))-1 gives S>=floor(log2(N+1)). The document correctly excludes unbounded current integer/carry/clock sidecars from its no-go claim.
* The residue-polynomial obstruction is correct. On arbitrarily large n=-1 modulo2^(M+1), one native a=1 step expands n while preserving its residue modulo2^M. Nonnegative polynomial restrictions on an unbounded ray have nonnegative leading behavior; any nonconstant such polynomial is eventually increasing. Constants cannot strictly decrease. This does not rule out general oscillatory or source-dependent ranks.

The independent checker recomputes every delivered depth16 certificate from its native word, checks infinite-ray payment by exact slope and least-source endpoint, verifies disjoint coverage, and separately follows each of32768 odd residues until first individual descent or precision exhaustion. Its words and unresolved list agree exactly with the delivered atlas. The latter equality is checked only at the declared precision, not asserted for every possible atlas design.

## Positive generalization of the interval carrier

Added `recurrence/GUARDED_MACRO_INTERVALS.md`: the inherited source/endpoint lift works for any authenticated guarded odd-affine macro. Every conjunctive affine observation of the common original source and endpoint reduces to a scalar half-line; historical conjunctions therefore collapse to a largest lower and smallest upper bound. This is a substantive positive sufficient-state result, not a claim of global minimality.

A disjunction q<=1 OR q>=3 is a decisive nonconvex hostile; congruence conditions require explicit cylinder splitting. Disjoint macro guards give additive source weights, and exhaustiveness is separately needed for conservation. A growing-only finite menu cannot silently replace the complete native transition domain.

The checker verifies the macro lift and mixed source/endpoint inequalities for all five growing words after a nonempty prefix, at12 parameter values each. Algebraic substitution proves its unbounded scope.

## Reproduction

    python3 check_independent.py > normal.json
    python3 -O check_independent.py > optimized.json
    cmp normal.json optimized.json

Receipts agree byte-for-byte. Only Python's standard library is required; the checker reads the relative `../descent/atlas16.json` artifact but imports no producer implementation. No global termination, least-survivor growth, or tail-contraction claim is certified.
