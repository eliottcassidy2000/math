# Proof and computation audit

**Status: PASS, same-agent audit with alternate computational paths.
Not independently refereed.** Proof labels are local; no canon theorem IDs were
reserved or promoted. All listed finite universes are complete unless a bound
or chosen family is explicitly stated.

## Universes and positive / hostile controls

| Checker | Universe and inherited filters | Alternate path / controls |
|---|---|---|
| sign_spectrum.py | Full S_n cycle-sign sums n=2..7; literal support and completion cubes n<=5; bigons k=2..8 with literal cubes k<=3; theta p=3,5,7 and k=2,3,4 | Actual weighted flip-matrix moments when d<=6; full-group two-arm tree and cycles-plus-isolate as positive controls; every C2 target type for each bigon |
| paired_oddness.py | All 2^(n^2) labelled aligned inputs and all S_n actions for n=2,3,4 | Literal affine class orbits versus signed support orbits; all eight flips; cancellation witnesses; classical counts 8,88,2752 |
| lrc_reflection.py | Primitive subsets of {1,...,10} of sizes1..4, thresholds the distinct members of {1/(k+1),1/4,2/5,1/2}; 1274 interval cases and50960 grid cases q=1..40 | Exact Fraction interval intersections versus direct integer residue tests;840 rational phase checks; separate AP13 boundary example |
| unit_distance.py | All45 pairs of the explicit ten-point algebraic configuration; seven rational circle neighbors; eight-vertex rigid orientation | Symbolic zero/nonzero checks, not tolerance comparisons; K_(2,3) two-circle obstruction; eight distinct directed degree pairs |
| inherited hostile | The supplied 12-vertex AAA/BAA substitution pair | Exact isomorphism search and path DP: groups81/27,67 identical labelled cyclic triples,462915 paths each, named guard true/false |

Normal and optimized checker counts are respectively the same:
470,25,109177,136. These counts measure active checks, not independent theorems
or statistical evidence. No random sampling is used in the new finite receipts.

The sign and paired checkers share the primitive coordinate-action helper.
Their literal orbit enumerations are alternate paths, not independent
implementations of every convention. The weighted matrix uses literal
coordinate flips and independently verifies spectral moments through degree
d+1. The imported 81/27 control is inherited code, not a new implementation.

## Proof review

- **Types:** source A-orbits are marked completion classes, not unmarked S_N
  tournament classes. Target X carries the specified same A action.
  The sign-even test is relative to the joint stabilizer.
- **Complements:** S1-S2 use free-pair complement. P2-P3 specify exactly which
  of T,W,mask is flipped. No equality with native Hamilton-path complement
  is inferred without a codec.
- **Necessity:** the spectrum obstruction is necessary for an operator-preserving
  isomorphism, not a sufficient classification invariant. The full-A Euler
  equivalence is inherited with that hypothesis intact.
- **Weak counts:** every nontrivial C2 action on the bigon target is free on
  tournaments. The r=0 case is handled separately. Only the natural edge-set
  target has the stated exponential complement-fixed count.
- **Oddness:** tournament automorphism groups are odd by Cauchy's theorem and
  the arc on an involution's swapped pair. This is used only after establishing
  that one full tournament is fixed. Reversing both cores needs the separate
  anti-automorphism argument.
- **LRC endpoints:** closed intervals and isolated points count as components;
  the origin is excluded because delta>0 and V is nonempty. Dilation and
  primitive normalization are distinguished.
- **Geometry:** an embedding means distinct Euclidean points, with every
  nonedge checked for the faithful example. It need not be a crossing-free
  drawing. The all-size length-three construction deletes finitely many
  forbidden parameters at each step, leaving rational choices.
- **Capacity:** the decoder sees the abstract faithful graph, with no real
  coordinate or uncharged metadata channel. The lower bound uses a labelled
  output upper bound, hence is valid also for unlabelled output.
- **Novelty:** original equinumerosity, ordinary Walsh duality, tree maps,
  bigon/theta tier separation, and Eplett counts are credited. The new packet
  adds scoped operation spectra, explicit flip actions/censuses, LRC trace
  limits and geometric comparisons. No external priority is established.

## Repairs and rejected inferences

Two historical reflections now carry maintained qualifications. The degree-six
claim is refuted by an exact degree-seven witness. The LRC/self-converse
identification is marked unsupported, rather than falsely declaring that
an undefined correspondence has been disproved. The old lattice construction
and tournament parity facts survive within their actual scopes.

During proof writing, the length-three geometric argument was stated using
finite forbidden parameter sets: b=(1,0) would require a_x=1/2, which the chosen
interval excludes. A blanket assertion that all unwanted terminal coincidences
occur only at parameter zero would be false and is not used.

The paired census's separate-sign filter differs from Session90's global
even-graph filter; their counts must not be compared as if they were the same
predicate. The explicit change of variables in P1 resolves the apparent mismatch.

## Reproduction and source pins

Run reproduce.py from any directory. It copies the four checkers to a temporary
directory, runs each with and without -O, compares all mathematical receipt
fields, and checks the inherited hostile. It does not rewrite repository
receipts. Only SymPy version and a floating nonedge-gap diagnostic are excluded
from comparison; all coordinates, edges, exact check outcomes and counts agree.

manifest.json pins this packet and the actual inherited computational / theorem
dependencies using SHA256 after CRLF-to-LF normalization. make_manifest.py is an
explicit maintenance operation, not called by reproduction. A changed dependency
must be reviewed before refreshing its pin.

The repository documentation check passed (bounded startup surface,
frontmatter-derived status, correction sentinels and router smoke). All local
links in the five packet notes resolve after path normalization; the git diff
whitespace check passed. Normal and optimized reproduction passed, including
the inherited 81/27 hostile. No major-open-problem closure,
unrestricted carrier minimum, or all-order algorithmic speedup is claimed.
