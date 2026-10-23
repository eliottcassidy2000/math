# Session40 first independent audit

**Status: PASS.** The source packet is `hopf/MANIFEST.json` (or its original working folder `session40-hopf-incidence/MANIFEST.json`) with SHA256 775ebd76ac73d78082793af063a5905b838bcb55fe1d56098ed2d539242207aa. The note SHA256 is503684db4fb36ef00571ab022bae1b05beee7233fa4789ab2a69edd393a682c1. No claim of external novelty is made.

## Proof review

The chordal gluing proof is valid for actual bijective edge transports, without requiring commutativity. Triangle identity reduces every cycle by chords. In the F3 translation restriction, the hidden phase dimension is cycle-space dimension minus actual triangle-boundary rank. One marked owner does not force nonflat cycles to become flat.

The fixed-root imaginary-Bargmann tournament encoder has a genuine local-transitivity restriction: projected outneighbors lie in one open semicircle and their angle order is transitive. The source-over-C3 and P7 counterexamples therefore exclude these tournaments from this specific encoder. Simultaneous complex conjugation of all spinors, including the root, flips imaginary products and hence all arrows. This does not imply that the native electrical linking form admits an anti-isometry.

A graph base has H²=0, so connection holonomy must not be identified with a nonzero Hopf first Chern class. The note states this boundary correctly. The Bargmann identity and handedness hostile retain actual shared triple data rather than independent overlap magnitudes.

The five-tap modulo8 binomial identity is correct for all dyadic depths2^r with r≥2. Its quarter taps are parity data, middle tap modulo4, endpoints modulo8. The two quarter parities may be aggregated, so coordinate-sampling minimality must not be called unrestricted encoding minimality. The exact binomial valuation proves the surviving support. A homomorphism C3→a finite2-group is zero, but nonlinear encodings are not excluded.

Square-shell reflections and their iterated affine endpoint formula retain partial-domain restrictions. They do not characterize Hamiltonian feasibility by the endpoint formula alone. The inherited strict-summand corrections are stated correctly.

## Independent controls and additional exact carrier

`check_hopf_guards.py` independently performs8,833 checks: every binomial coefficient through dyadic depth4096, valuation controls, all64 four-vertex graph supports and4,096 F3 edge-phase assignments, and the Bargmann handedness hostile. The producer's larger S3 and shell-word checks are separate.

`MATRIX_TREE_PATH_PROJECTION.md` proves the exact Hamiltonian degree-sector projection of the classical weighted matrix-tree polynomial. It includes endpoint pair coefficients and mixed edge derivatives. The checker independently verifies all74 simple graphs of orders2,3,4, every endpoint pair, and every first/mixed-second edge derivative:3,057 controls. This is an exact all-size expression, not an efficient Hamiltonicity algorithm or a tournament wiring classification.

Run both Python scripts normally and with `python3 -O`. The matrix-tree checker requires SymPy; the Hopf guard checker uses only the Python standard library. Both scripts use explicit failure conditions, so optimized Python retains the checks.
