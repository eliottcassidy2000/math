# Validation and provenance

Producer replays: `normal.json` and `optimized.json` are identical JSON, generated
by separate ordinary and `python3 -O` processes. They contain literal witness
rows, orbit/error histograms, universe hashes, and named check counts.
All runtime tests use explicit exceptions. No external packages are needed.

Independent checks within the program:

- Regular-seven classification: residual-degree generation, followed by all7!
  vertex relabelings and separately backtracked full automorphism actions.
- Five-cell crossings: complete binary regular-matrix generation, literal
  C5×C5 orbit partition, separate fixed-matrix counts, analytic bipartite-cycle
  count2040, and full vertex automorphism backtracking of each class.
- Seven-cell substitutions: all64 quotient orientations, all root-preserving
  S3 quotient relabelings, all three internal seven-types, actual refinement,
  strong connectivity, and independent full vertex automorphisms.
- Triangle/singleton crossings: all13824 presentations; literal rotation/block
  permutation orbits on the10368 strong presentations; direct automorphisms;
  all3456 constant-attachment exclusions. Cyclic constants must collapse9+1;
  transitive constants retain3+3+3+1 but are not strong.

The resulting strong-class counts are96 in the5+5 sector,12 in7+1+1+1,
and128 in3+3+3+1 with all crossings nonuniform. Their positive-error class
counts are84,4,32, and error contributions3360,56,576 respectively.

This is not unrestricted target enumeration. The exhaustiveness bridge is the
proved stable-cell partition classification in THEOREM.md. Other order-ten
classes exist, but the proof shows their error is zero or they are not strong.

The total-order-ten coefficient uses inherited E_s,9=716. The Session50 source
statement and result receipt are byte-pinned under `sources/`; runtime checks
of3992 do not depend on them. No target a_10 value, order-ten atlas, full
order-nine atlas, or bulk U_10 trace is accessed.

`verify_manifest.py` checks exact file membership, byte hashes, archived source
pins, and normal/-O equality. It works after a plain directory copy to an
unrelated location. The manifest excludes itself and optional verification
receipts; it does not exclude any source file or mathematical deliverable.

Distribution extension: `distribution-normal.json` and
`distribution-optimized.json` agree byte-for-byte. The checker independently
reconstructs all765 order-eight reroot destinations from85 inherited parents,
with15 regular-nine components. Including canonicalization of the85 supplied
parents,850 eight-vertex canonicalizations are performed. None is performed on
nine or ten vertices. The resulting error histogram is0:3,18:3,36:2,72:7;
combining it with the family certificates and SCC padding yields the160-class
order-ten positive-error distribution and raw second moment.
