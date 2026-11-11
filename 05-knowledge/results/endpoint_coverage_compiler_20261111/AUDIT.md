# Exact replay and mathematical scope

**PASS for all packet computations in normal and optimized Python.**
This is same-agent review using independent counting implementations, not
an independent-person proof audit. `audit.json` records each command and
byte-for-byte agreement of its normal and `-O` result. `manifest.json`
pins packet files; the two imported prior modules are separately pinned
inside `audit.json`. Run `python3 audit.py` to reproduce the entire replay.

## Universe and positive controls

- Ordinary tournament counts 1,1,2,4,12,56,456 for orders 1..7, generated
  by all one-vertex extensions and exact degree-cell canonicalization.
- Every nonempty selected set and unselected port of every unit-offset
  base in that universe: 1,407 presentations, of which 634 are stable.
- 1,276 actual inflated-graph endpoint controls; coefficient inversion
  through the full required size bound for each of ten named families.
- All 763 separately stable port pairs: literal sizes (2,2), complete
  two-variable coefficient certificates and 127 mixed-size controls.
- All 4,054 one-vertex children of ordinary bases through six vertices
  independently check the response polynomial. The content-pair census
  through seven uses exact endpoint coefficients and all Boolean cuts.
- 60 actual constructed tails, 30 direct path-backtracking controls,
  ten nontransitive C_3 substitutions, and six requested 99-ray points.
- 5,000 direct prime-wheel exponent controls, 96 binary-language controls,
  exact CRT inclusion-exclusion and a prime certificate for 419 by trial
  division. No probable-prime assumption is used.

## Hostile controls and retained failures

The full d vector is necessary: testing only block size two misses the
(1,0,-4) defect; the (1,2,-2) example has nonmonotone successful sizes.
Their success sets {1,2} and {1,3} are proved exactly, not inferred from
the six checked values. An arbitrary module at a stable port has offset
H(M), verified with C_3 and proved by contraction.

The converse pair with H=37 and W=14 has U=99 versus 117. Equal-H and
equal-defect projections therefore fail as tail states. The ternary30/90
overlap, by contrast, is a verified selected-set isomorphism and is counted
as one coefficient orbit in the summary.

The wheel accepts infinitely many composite values: exponent
2+8140892760*t, t>=1, always gives a proper divisor 419. Tournaments
realizing those counts are still constructed by the same tail rule.

## Scope reviewed

- E1's D_(n-1)=2 has n>=2; the n=1 boundary is separately stated.
- E2's finite basis bound depends on an unselected port and the specified
  endpoint conditions. The arbitrary-module formula is not a formal
  replacement of P_j(M) by a single scalar.
- E3's s>=1 cases are complete lines. Higher starts are reported as tails.
- E4 counts marked presentations over ordinary bases, not automorphism
  orbits of marks. The pair-stability conclusion is limited to bases of
  order at most seven, though its block sizes are unbounded.
- E5 preserves original arcs and both endpoint roles. Its rational-series
  formulation does not prove coverage of all integers.
- E8 excludes only the stated inflation-only additive-counter compiler.
  It does not exclude all finite-state or Zeckendorf-based constructions.
- The inherited all-c<=4097 bank includes different ray mechanisms;
  this packet constructs selected-set tail seeds for its named families.
- No new universal seed theorem, Chan–Pak proof, primality criterion,
  improved A000568 formula, or global optimal graph-order bound is claimed.

The bounded navigation check is rerun after integrating incoming commits;
the initial run reported an existing tight startup-byte budget plus the
proposed frontier insertion. The insertion was removed pending integration.
