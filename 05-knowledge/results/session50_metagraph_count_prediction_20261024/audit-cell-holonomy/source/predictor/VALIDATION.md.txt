# Validation and prediction boundary

## Exact results

Normal and optimized Python runs produce byte-identical stdout and receipts.
The blind result is 54,256 pointed states, refinement upper trace 55,068,
strong order-seven error 14, total order-eight correction 28, corrected trace
55,040, and A000568(8)=6,880.

The recorded operation counts are:

- 54,256 child extensions and 54,256 child refinement calls.
- 532 nonempty input-parent refinement calls.
- 61,944 parent-subset action evaluations.
- 34,801 parent-action arc comparisons.
- 132,603 refinement rounds across parent diagnostics and child calculations.
- Zero child canonicalization calls and zero deletion/re-rooting isomorphism
  calls, by the inspected source boundary and absence of such dependencies.

These are operation counts, not measured speed or complexity comparisons.

## Independent input validation, outside blind prediction

`freeze_inputs.py` executes the frozen growth generator only through order
seven, then serializes representatives and their full actual actions.
`validate_inputs.py` independently rechecks each supplied group against that
source's exhaustive degree-cell permutation search, checks pairwise
nonisomorphism of all representatives, and sums n!/|Aut(T)|. For orders zero
through seven the class counts are 1,1,1,2,4,12,56,456 and the labelled masses
are 1,1,2,8,64,1024,32768,2097152, respectively. These are exactly
2^(n(n-1)/2), so distinctness and complete actions prove complete atlases.
No validation function is imported by the blind entrypoint; its largest
canonicalized graph order is seven.

This validation uses lower-order isomorphism computations openly. It is not
part of a claim that no isomorphism information is needed anywhere.

## Hard blind boundary

`predictor.py` imports only seven named standard-library modules. Its source
AST guard rejects unapproved imports, indirect function calls, named
canonicalization/deletion/re-rooting operations, dynamic execution, and target
answer or trace constants outside the guard's own blacklist. A runtime audit
hook permits reads only of the executable source and the SHA-256-pinned input
atlas, permits writes only to the requested receipt, and forbids subprocess,
network and native dynamic-library activity. Every mathematical check uses an
explicit exception and remains active with Python -O.

The supplied input groups are validated as permutation groups preserving
parent arcs. Full group/atlas completeness is separately certified, not
claimed from closure or from the hash alone.

The output's zero-call fields document this inspected implementation; their
basis is not merely an initialized counter. The audit receipt lists its full
call/import surface. An independent audit is additionally required before
publication and is owned by the surrounding session's audit packet.

## Hostile controls

`hostile_controls.py` proves rejection of:

- an import of the generation module;
- a child canonicalizer call;
- a deletion lookup call;
- dynamic execution;
- injected target and corrected-trace literals;
- a changed/incomplete parent atlas;
- unsupported target order nine;
- an unexpected runtime file read.

A minimal directory containing only the predictor, atlas, poisoned answer
receipt and poison generator still reproduces the correct baseline without
using either poison file. This establishes that the historical source and
answer receipts are not required blind inputs.

## Replay

`replay.py --manifest-only` verifies literal hashes and archived source pins,
then copies only manifest-listed files into a fresh temporary directory. It
runs the predictor, the separate input validator, the separate input generator,
and hostile controls under normal and optimized Python. Each deterministic
output is compared byte for byte; generated input bytes are compared with the
frozen atlas. It then runs the blind predictor from a second minimal directory
containing only the two required files. `replay.out` records the completed
replay; the manifest excludes its own digest and this replay output to avoid
self-reference.
