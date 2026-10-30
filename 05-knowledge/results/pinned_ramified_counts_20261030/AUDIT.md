# Audit and reproduction

**Status:** displayed all-size proofs reviewed by the same agent. Exact
computations compare literal permutations with symbolic cycle inventories.
No independent human/subagent audit, Lean build, or literature-priority
claim. The Hexagon theorem is a cited preprint claim; its application in
A3 remains conditional. Arithmetic identities and fibre classifications
from Sessions78/81/82/108/111 are credited as inherited.

## Validity review

- Domain: p>=5, k>=1, v_p(D)=1; positive integers in arithmetic. S is
  nonempty and pointwise named. Empty S, setwise pins, unmarked roots and
  arbitrary complements are outside the new cycle-index theorem.
- Directed-color relation and actual ring vertices are retained. The
  clock is the specified p-power unit, not field Frobenius. No eigenvalue,
  cyclic-triple statistic or abstract group order substitutes for it.
- In P1 the p-fold duplicate parametrization is divided out exactly once.
  At k=1 the whole root group consists of fibre shifts. The earlier false
  unrooted cycle extrapolation is not used.
- Fixed points of all powers determine cycle partitions, not just ranks.
  The r>=2 argument uses h(0)=0 and the strict socle valuation gap.
- P3 first takes an order-p power, then uses disjoint supports. Without
  the root this reduction was not proved. The probability statement needs
  IID across vertices; deterministic testing needs no randomness.
- LCA and j classify the displayed cycle index, not full pin orbits.
  Two same-size pin configurations with tau1/tau3 have different orbit
  counts despite equal rigidity probability. A frame orbit is not an
  ordinary tournament class.
- A1 strips only primes at which the scale is a unit. A height bound on
  the actual target count remains compulsory. Small moduli do not remove
  physical integral realization guards or prove a smaller carrier.
- A3 retains fixed K, all exponent-profile candidates, integer and divisor
  guards, and the preprint dependency. No claim of practical speedup or
  of a generic counting engine for tournaments is made.

## Exact finite universes

1. Complete literal root groups at (p,k,D)=(5,1,5),(5,2,5),(5,2,10),
   (7,1,7):625,3125,3125,117649 permutations. For each, every listed pin
   stabilizer is filtered from those permutations and its ENTIRE cycle
   distribution compared to P1. Check every element's order-p descendant;
   compare the order-dividing-p set to all fibre shifts. Root, every
   residue occupancy2..p, and available tau1/2/3 cases are included.
2. Additional selected maps at (5,3,10) and (7,2,21) check actual higher
   clocks and socle perturbations. These are selected controls, not a
   full group census at those larger orders.
3. Exhaust all1024 binary words on two free five-point fibres, all25
   translations:64 orbits,900 rigid words,36 rigid orbits. This independently
   tests the consequence of the probability theorem. Sparse constructions
   are checked against every element of each full root group above.
4. Symbolic inventories and recurrence evaluations cover p=5,7,11,13,
   k=1..4, all ell and all j. Large evaluations use exact residues modulo
   1000000007 WITH the group-order factor retained before division; no
   floating point or giant claimed explicit enumeration is involved.
5. Regenerate the complete19-source343 fibre by proper-divisor sieve
   through29241, and independently evaluate all prime2/7 queries by
   factorized sigma. The sharper finite bound and the inverse formulas
   are inherited. The finite minimal modular separator search uses the
   pigeonhole lower bound19, then exhausts moduli19..71.
6. Scaled-zero tests exhaust C=0..2^b for b=0..10, scale bases2 and10,
   E=0,1,17,1000000. Check unit inverses, exact lcm bounds, and both
   missing-height/nonunit hostiles. No huge scale integer is constructed.
7. Replay the native81/27 transporter control:67 identical labelled
   cyclic triples,1620/1728 regular binary masks. No new three-run or
   complement-recovery assertion for that pair is made.
8. Incoming Session112's reflection-axis theorem supplies P5's source
   hypothesis. Independently filter250 actual unsigned two-pin maps at625;
   verify all125 reflections' cycles and1000 subgroup orbit counts using
   union-find, then compare directed/unsigned palette counts. The integer
   and tilted pin sets control phase presence and absence. All1024 binary
   words on a smaller diagonal-reflection action test the difference
   between900 directed-rigid and zero unsigned-rigid words.

## Incoming work and attribution correction before publication

Fetched `fea3b07542`/`ecd0185d34` before publishing. Their modular observer
already has exactly A1's cutoff, so A1 is marked as a credited replay.
No mathematical statement was false and no prior packet publication needed
retraction. P5 instead transfers the incoming reflection theorem to the
new joint pin/color observer, with the source theorem explicitly pinned.

## Reproduction

From the repository root (Python standard library):

    python -X utf8 -B 05-knowledge/results/pinned_ramified_counts_20261030/check.py > receipt.json
    python -X utf8 -O -B 05-knowledge/results/pinned_ramified_counts_20261030/check.py > optimized.json

Compare these with the committed packet receipt. All gates use explicit
exceptions; `-O` does not disable them. Native group/cycle checking and
symbolic enumeration are distinct paths; both reuse the pinned Ring
arithmetic implementation. Thus this is not a fully independent source
implementation or independent proof review.

    python -X utf8 -B 05-knowledge/results/pinned_ramified_counts_20261030/seal.py

The manifest checks LF-normalized artifact bytes, immutable source pins
and the two live runtime dependencies. It is an integrity check, not a
proof checker. Exact final gate counts and replay hash are recorded below
after both runs finish.

Final replay: **3,423,143 local gates plus7 inherited hostile gates**.
Normal and optimized receipts agree byte-for-byte; receipt SHA256:

    95cf132c6629576eb4917bd8e350f2e5097b02b2b745980361ed04e2ea3d56f8

The complete groups comprise124524 literal root permutations, with45
additional selected higher-depth maps. The phase check has250 literal
maps and1000 generated-subgroup checks. The zero-test replay has16464
scaled-count cases. These are finite validation universes, not the proof
of the all-size statements. After integrating the credited incoming
checkpoints, `agents/check_docs.py` passed, `git diff --check` passed,
and the manifest verified7 artifacts,12 immutable source pins and2 live
runtime inputs. No full experimental replay was repeated after this
documentation-only integration; the pinned runtime bytes are unchanged.
