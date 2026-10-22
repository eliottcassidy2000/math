# Independent audit: Session29 seamed metagraph transfer

Verdict: PASS, within the stated edge-weighted family scope. No counterexample found. The all-word result rests on a constructive path-cover transfer, not extrapolation from sampled word lengths. This audit independently reconstructs the alphabet and matrices without importing any producer code.

## Independent checks

`independent.py` uses endpoint subset dynamic programming, independently regenerates the 5+13 seed spine tilings, constructs all 90 state/letter atom graphs, and derives their five-state path-cover signatures. It verifies:

- All 90 rows agree with the producer's permutation/cut enumerator.
- The DE/DB collapse and E/Q compression hold exactly for each of the 18 letters.
- Each compressed first column equals the actual one-block port signature.
- All 18 per-letter mirror identities K_tau(t)=J K_t^T J hold, J=diag(1,-1,1).
- All 324 full ten-vertex target Hamiltonian counts agree with matrix products.
- The six specified general length-three full 15-vertex targets agree.
- All 72 mirror-fixed length-three targets agree individually; their sum is 2290528.
- Aggregate K and the degree-three rational GF agree exactly.
- An independently formed six-dimensional quadratic-moment operator gives the stated degree-six even and odd blue rational functions directly by symbolic inversion, without recurrence fitting.

Normal run is recorded in `run.out`; `evidence.json` contains all independently computed raw and compressed matrices. Tests use explicit exceptions, not removable assertions.

## Constructive proof audit

A prefix of a final Hamiltonian path has at most two path components: at most one can leave the prefix through its unique outgoing seam; every other component must terminate at the final Hamiltonian endpoint. The same restriction applies to an intermediate two-component cover whose active component ends at the new portal: at most one nonportal endpoint is permitted.

The old components are maximal prefix runs of the new cover. Consequently old-atom to old-atom arcs must be forbidden in the contracted transition graph. Allowing them would count nonmaximal cuts and destroy uniqueness. The producer correctly forbids them. New-to-old arcs depend only on whether the old portal is a component start; old-to-new arcs depend only on whether it is a component end. These are precisely the five recorded states. Thus contraction and expansion are inverse maps, with multiplicities equal to old path-cover counts, for every word length.

The DE/DB equality also has a direct explanation. With two old components, atom1 has no outgoing seam and must furnish the sole nonportal endpoint. Atom0 must therefore leave through the seam to the new first vertex. An attempted edge from that new first vertex into atom0 would close a directed cycle, irrespective of whether atom0 represents a singleton. So the one apparent DE/DB adjacency difference never participates in a valid cover.

The further three-coordinate reduction is a finite exact identity for this specified alphabet, not claimed for arbitrary block families. This is sufficient: the exhaustive 18-letter certificates establish a common quotient of the already proved all-word transfer.

## Graph and colour conventions

Prior THM-4729 L3 proves the nonstrong-world independence and is unaffected by its separate L4 pure-blue proof gap. For every source word with r>=2, first-to-last is a nonspine chord, and complementing it closes the preserved spine to a Hamiltonian cycle. Every target is therefore strong and outside this nonstrong source family.

Each source tiling contributes one distinct outward line. There is no factor of two: no line has both endpoints in the source family. Converse merging changes node assignments but does not identify distinct tilings or parallel lines. Blue means a mirror-fixed tiling, not an edge merely fixed setwise by a reflection. Complement and mirror commute, so both endpoints have the same colour.

Hence the line totals 18^r, blue totals 18^(r/2) for even r and 4*18^((r-1)/2) for odd r, and their difference for black are correct. These are multiplicity counts, not distinct strong-target node counts.

THM-2183 concerns unrestricted arc-reversal distance. Its Hamming-product consequence is a separate lens and supplies no fixed-spine-complement adjacency assertion. The addendum preserves this distinction.

## Certified formulas and scope

K=[[172,-24,408],[24,-4,68],[408,-68,706]], h_r=e0^T K^r e0.

For r>=1, the GF is (172z+45144z^2+50400z^3)/(1-874z-43344z^2-50400z^3). Actual outward-edge H mass begins at r=2, so subtract 172z. The r=2 total is195472 and its blue contribution11522; black is183950.

The even blue resolvent includes the empty-word term1; the odd blue resolvent includes the one-letter term42z. Remove both when interpreting actual outward edges. The addendum's subtraction of130z from total-minus-blue for actual black outward edges is correct.

For mirror-fixed even words the exact identity H_target=H^2-S^2+R^2 is a real structural connection between mirror and an indefinite adjoint. It does not identify an integral Lorentzian lattice, a Leech cusp, or a Dehn invariant, nor does indefinite self-adjointness alone establish a general real spectrum. The displayed aggregate's three distinct real eigenvalues can separately be verified from its cubic.

Remaining open: target canonical classes, target automorphism orders, collisions of different lines at the same target, and complete metagraph wiring. This audit does not promote those stronger claims.

## Producer sources reviewed

`/workspace/shared/session29-profiles/ADDENDUM.md`, `seamed.py`, `transfer.py`, `check_compression.py`, `blue_transfer.py`, `lorentzian_transfer.py`, and their JSON outputs. Hashes of the reviewed source snapshot are in `SOURCE_SHA256SUMS`.

## Portability amendment, 2026-10-22 07:07 UTC

PASS. `check_compression.py` changed exactly two source lines: the JSON input and output absolute paths became paths relative to the script via `Path(__file__).parent` (with the Path import on the input line). Reversing these two substitutions reconstructs the previously reviewed SHA256 exactly, 4c18cba563e34cd268c45939e46c1bbfc6bd65d2b297e5ab7f187f60372c5760. The amended hash is792d8b451e98b62b1b520d029f3dc06324ad9fc12f1016e0e34b2bcef46b51f1. Every other reviewed source hash remains unchanged.

An independent isolated scratch-copy replay from working directory `/tmp`, in both normal and optimized Python, reproduces identical stdout and compressed-matrix bytes. Output SHA256 remains14279c5ef98948cd67fbf2aacc86a377fe6cb3c53d3eba8cf4c025b096955d68. No mathematical change. `portability_audit.py` and `portability.out` preserve this verification, and the reviewed-source and deliverable manifests are refreshed.
