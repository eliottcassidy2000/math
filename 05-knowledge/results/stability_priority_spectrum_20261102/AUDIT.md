# Verification and scope audit

Status: same-agent proof review with separate native computational paths.
No new independent-agent or formal verification is claimed. Session139's
526 proof and its audited local inequality are credited dependencies; the
new320 combination has the proof and checks supplied here.

## Reproduce

Run in this directory, with Python's standard library only:

    python -B verify_stability.py
    python -B verify_stability_320.py
    python -B verify_clocks.py
    python -B verify_carriers.py
    python -B seal.py --check

All gates use explicit exceptions and remain active under `python -O`.
Normal and optimized runs have byte-identical JSON receipts. The total is
45,029 local checks plus seven inherited81/27 checks. Quantifiers over all
orders, all moduli, and all horizons follow from the proofs; finite replay
does not establish them by itself.

| Checker | Universe and independent path | Receipt |
|---|---|---|
| Stability466 | Rational endpoints; symbolic deletion identity;24 native weighted matrices at orders4..6 with full independent rounding; fair pairs4..12 |3,475 gates |
| Stability320 | Rational dual and tail endpoints; polynomial identity; orders320..1499 finite replay; even/odd recurrences |3,683 gates |
| Arithmetic clocks |132 ordered programs of length1 or2 over reduced a/b with1<=a,b<=4, each modulo1..12; first9 trace horizons; native inputs1..10000; compatible CRT charts |37,806 gates |
| Carriers | Actual six-vertex spectral witnesses and source extensions; exhaustive square-sum path existence1..25; all169 quadratic-ring states; literal five-vertex cycles; weighted order-five zeta witness |65 local plus seven inherited gates |

## Hostile controls and failed implications

- Ordinal fair pairs attain the proposed coefficient; they defeat a
  dimension-free lower bound and the unproved normalized five-vertex
  local-to-global transfer. At order6 the latter fails by3/32.
- The old cycle-support quotient changes15 to11 at five vertices. The
  literal cycle list and native Hamiltonian DP agree on15.
- Binary cospectral tournaments have H=13,17; weighted cospectral
  tournaments already differ by1/3456 at order5. These are exact matrices,
  not approximate eigenvalue computations.
- The twelve-vertex pair has67 identical labelled cyclic triples and equal
  H=462915, but automorphisms81,27 and ordered-path presentation fibres
  5715,17145. The reused automorphism routine's source is pinned.
- [1/2,2/3] at horizon2 has minimal residue period12 versus the naive36.
  At horizon9 its124416 residue classes compress to35 trace responses.
  This does not give a finite carrier for all complete traces.
- The repeated1,8 square-sum walk passes each edge test but fails vertex
  uniqueness. A transitive tournament order exists even for an empty graph;
  adjacency markings decide whether it represents a Hamiltonian witness.

## All-size proof audit

For320, the retained matching has at least eight edges because
g(160)=98867/12500>7 and g increases. The geometric tail is bounded per
deleted edge by a discrete endpoint argument, and both endpoints decrease
with R. Every inequality keeps B until the final contradiction. The local
dual uses W_i from disjoint groups; W is spent once. If1-lambda R is
nonpositive, the positive D coefficient already contradicts the budget;
otherwise substituting B<=kappa D preserves the inequality direction.
The L=0 case uses t>kappa, forcing D=B=0 and hence the inherited matching
equality locus. The V=0 boundary is handled by the binary equality theorem.

For FRACTRAN, fractions are reduced, inactive rules are identified by
earlier denominator divisibility, and original indices survive deletion.
HALT erases the input for the specified observation. Least period is
proved by testing N=b_i, and the all-modulus simplification is primewise.
The trace-period necessity proof separates the2-adic and3-adic coordinates;
its response-class count has an independently enumerated valuation model.

The general spectrum, realization of all prime-power context states,
weighted orders10..319, LRC(14), RH and Collatz remain open. The session
does not claim a new complexity theorem from an invariant collision.

## Repository checks

Owned-path diffs pass `git diff --check`. The documentation checker is run
through a Windows extended-length path because an existing canon filename
exceeds the ordinary path limit. The initial long frontier line was
shortened to respect the bounded startup-surface budget. Historical proof
packets remain frozen; current navigation points to the320 refinement.
Source pins and artifact hashes provide the exact reproduction boundary.
