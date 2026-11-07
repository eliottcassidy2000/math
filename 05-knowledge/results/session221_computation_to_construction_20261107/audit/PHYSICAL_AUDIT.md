# Session221 single-arc construction audit

**Verdict: PASS.** Independently audited `session221-cellular/PROOF.md`, its stated physical model, all quantifiers and formulas. The main theorem is a physical difference constructor: for every integer `N>=1`, two strongly connected tournaments on `5+ceil(log_2 N)` vertices differ in exactly one arc and their Hamiltonian-path counts differ by `2N`. This is not absolute spectrum coverage. No publication was performed by this auditor.

## Correctness and scope

The inherited one-apex weighted inversion formula gives the adjacent `10 -> 01` response `2(1+L+R)`. Its immediate exterior neighbors are response-null. Choosing the complemented binary digits of `N-1`, then the padding bit 1, realizes the prescribed half-response `N`. Direct expansion gives the ordinary absolute counts `4N-3-f(N-1)` and `6N-3-f(N-1)`.

Reversing the chain arc between those adjacent vertices, and then swapping their names, is exactly the tape exchange. All other chain incidences agree because the positions are adjacent. This independently proves the single-arc, rather than merely two-apex-arc, claim.

Appending the common final bit 1 preserves the response. For `N>=2`, both words begin 0 and end 1, so both are strong by the displayed spanning cycle. Their absolute counts are `12N-9-f(N-1)` and `14N-9-f(N-1)`. The `N=1` strong exception is `0101 -> 0011`, with counts 9 and 11 on five vertices. The ordinary `N=1` exception uses three vertices and counts 1 and 3.

At a fixed site the left and right binary registers independently range over finite integer intervals, proving the exact unrestricted response range and its convolution multiplicity. The maximum half-response at chain length `n` is `2^(n-3)` ordinarily and `2^(n-4)` for a move preserving strongness. Thus the displayed order bounds are optimal **within adjacent-chain one-apex flips**, not among arbitrary single-arc tournament moves.

The strong fixed-site intervals have nonzero lower bounds. In particular chain length 5 gives only half-gap 2, and length 6 only 3 and 4. The final proof correctly avoids claiming a full strong fixed-order initial interval.

The quadratic mixed term for simultaneous swaps at sites distance `d>=2` is `w(d-1)`, including the exceptional small distances. Hence the parallel interference formula is correct; higher interactions vanish. Rule184 conservation and termination follow from ordinary inversion decrease. The displayed suffix-observer identity also proves the scoped obstruction to a nontrivial, equal-length, context-independent additive rewrite.

## Independent physical computation

The audit imports no producer code for these checks. Its graph labels put the apex last, unlike the producer's apex-first convention. It uses a separately written backwards endpoint recursion on vertex subsets, rather than the producer's forward DP or weighted formula.

`audit_physical.py` completed:

- 1,023 actual tape tournaments, every tape of length 0 through 9;
- 1,793 local physical response comparisons;
- 968 parallel physical comparisons;
- 36 complete unrestricted fixed-site ranges;
- the ordinary selector for every `N=1,...,128`, each pair counted on literal graphs and checked for exactly two changed apex arcs.

`audit_strong.py` completed:

- 321 actual strong local moves through chain length 9;
- 321 exact one-chain-edge-flip / adjacent-relabel adjacency identities;
- 21 complete strong fixed-site ranges;
- every strong selector for `N=2,...,64`, with independent exact absolute counts.

The two scripts are independent reconstructions from the mathematical statement. Finite checks supplement the all-N proof rather than supply the infinite quantifier.

## Inheritance and originality boundary

Read from commit `38fdf9f234159aaedb2d0ba1e0f7bcdefb49e817`, the October20 `feedback_word_hamilton_spectrum_20261020/proofs.md`, sections 2–4, already proves the weighted inversion formula, exact adjacent response, strict rearrangement, strong interior increment-two exclusion, and fixed-content extrema. These are explicitly marked inherited in the candidate proof.

The assessed contribution is the all-even physical single-arc inverse, its exact base counts and encoding-optimal size, the strong lift, the complete marked response intervals, and the parallel interference identity. No claim of external literature priority was checked. General arc-flip formulas and evaluation automata are inherited. Neither signed responses nor even differences themselves are Hamiltonian-path counts.

## Compact release and reproduction

Release these audit files together:

- `PHYSICAL_AUDIT.md`
- `audit_physical.py`
- `audit_strong.py`
- `physical_receipt.json`
- `strong_receipt.json`
- `PHYSICAL_SHA256SUMS`

From their directory run:

    python3 -B audit_physical.py
    python3 -B audit_strong.py

There are no external dependencies and no dependency on the separate fixed-apex compiler packet. No unresolved correctness blocker remains in the audited physical construction.
