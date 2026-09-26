# Audit record: recursive families, square sums, and pairing attainment

**Status: independently audited proofs plus FINITE-EXACT computation.**
Ordinary Collatz remains OPEN. This is a current audit record; the proof
notes and canonical statements, not this narrative, define the results.

## Proof reviews

- **Flow / THM-4500:** root independently derived the local primal/dual
  telescope, L1 contraction, finite-policy vertex characterization, and
  the slow-phase attainment proof. Geometry re-derived them separately.
  Its independent controls covered 16 modulus-four policies, 64 state
  bounds, 2592 child minimizations from 81 potentials, 6561 L1 potential
  pairs, and 240 exact telescopes. For attainment it checked 3660 safe
  rational selector cases at three scales and 100000 actual ancestor
  iterates. Verdict SOUND. A finer-dyadic-conditioning justification was
  made explicit. No supremum-norm Bellman contraction is asserted.
- **Orbits / THM-4501:** root and geometry independently checked exact
  affine lifts, triangle scope, valuation disjointness, infinite-union
  density with its tail error, phase order, and dyadic closure. Geometry
  also checked 72 direct-membership/floor-count cases on small bases.
  Verdict SOUND. Fixed-depth growth wording was tightened. The 27 motif
  is exclusive within its selected triple, not the full prefix.
- **Fractal:** automata independently verified the block-IFS lower
  dimension, critical Hausdorff-measure upper bound, fixed-horizon
  natural density, and the two-completion hostile. It checked 2046
  residue words at lengths 1..10 above exact carry thresholds and 100
  direct rational compositions. Verdict SOUND. Dimension h is inherited;
  the newer sharp cover asymptotic supplies the critical-measure corollary.
- **Squares:** root independently checked the unicyclic deletion proof,
  additive C4 identity, degree-preserving switch boundary, explicit cycle
  witness, and signless kernel. Automata accepted the lift/operator map.
  An omitted degree cap in the first wording of the unicyclic criterion
  was repaired: the two deleted-edge endpoints must each have degree at
  most three. The structural equivalent, proof, and enumeration already
  imposed this. The degree-four two-tail witness is retained.

The independent checks described above were separate audit probes; the
following retained programs are the persistent reproduction paths.

## Retained finite universes and hostile controls

| Lane | Universe / retained exact computation | Controls and boundary |
|---|---|---|
| Flow | All 2^24 final residual states; rational modulus-32 primal/dual; all 16 modulus-four policies; 255 independent fixed-root subtree cases; exact selector through120000; actual modified sources through160000;479996 ancestor checks | Overflow bounds; mean1; L1 envelope; normal/-O; v3 policy fails6->9->14; L-infinity factor4/3; direct patch fails26/39 |
| Orbits | All2047 parity cylinders through depth10; odd sources<=10000 at inverse depths1..4;248 actual motif lifts; exact phase orders; floor counts at2^40,2^80,2^120 for four K | Literal27 inverse obstruction; common-tail growth distinction; selected27 versus full4347 support; explicit59-step descent; odd CRT loss |
| Squares | All paths up to reversal in Q_n,1<=n<=25;3898 labelled connected unicyclic graphs on3..6 vertices; triangles/C4 through46; two full46-vertex cycles | Self-loops excluded; degree-two endpoint22; disconnected2-factor after an otherwise valid switch; unused edge blocks lifted path |
| Fractal | Exact DP through512; independent brute words through16; Spitzer recurrence through128; all allowed blocks through8; dyadic quotients through2^12;80 exact rational prefixes | Integer barrier comparisons; every occupied prefix extendible; negative periodic points; one rational sequence has different real and dyadic limits |

All floating values are display approximations. The deep flow computation
uses NumPy integer arrays with explicit range guards; all certificate
fractions and sums are exact. The fractal DP does not use floating values
to decide positivity. Every assertion needed for the checks survives
optimized Python. Computations validate the stated finite universes;
the infinite conclusions rest on the accompanying proofs.

## Reproduction and hashes

    python -B 04-computation/experiments/crossroads_family_20260926_audit.py
    python agents/check_docs.py

The runner executes normal and optimized versions, then compares stdout
to the retained file (allowing only Git's CRLF/LF normalization). Its
`--record` mode is reserved for intentional output refreshes.
The manifest is generated from exact staged Git blobs and excludes itself
and shared navigation. Verify against committed blobs, not platform-specific
working-tree line endings.

Namespace reservation and first checkpoint: `6d5865d4e`. Scope promotions
followed independent acceptance. No result depends on an unproved reserved
file, on a numerical extrapolation, or on external Collatz proof claims.
