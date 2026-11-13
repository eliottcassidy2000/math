# Final audit, scope, and reproduction

**PASS for the scoped statements in this packet.** All-order claims are proved
in the notes; exhaustive computations retain their finite universe. Primary
inputs remain cited, and prior constructions remain inherited. The unrestricted
Hamilton spectrum, compact class-wiring recurrence and black-core connectivity
are not promoted.

The clean research base is `f90cd879d156a37f0027d1719a25b8731918db6a`.
Session315 is separately pinned at
`ae6c9620d8718e3ce5dd710ef6ecb1e00b93436a`. Source documents and byte lengths
are bound in the packet's source records; incoming proofs were read but their
larger paired-site census was not rerun. No chronological log was used as a
proved dependency.

## Proof review

**Joint fibres, P1–P4.** Root checked the finite-core formulas, the uniform
Vandermonde valuation bound, the half-core parity boundary, cancellation at
x=1 over F2, and the least-period argument including all possible integer
periods. A separate agent checked P1/P2 against synthetic binomial and rational
series families. The low-H completeness upgrade uses Busch's exact strong
minimum: from order seven it is at least25. Root reopened Theorem1 and
Corollary1 of the primary paper. Thus the all-order H<=21 claim has a stated
cutoff proof, not an extrapolation from six vertices. The fixed-count formula
and H19 rank use the inherited optimal pairing lemma; no new general theorem
on polynomial sequences is claimed.

**Memory and score closure.** Root checked the decomposition I=LE+R,
elimination of the hidden state, D-self-adjointness and the Gram equality
`M(F_2-Q^2)=B^T D B`. The zero/nonzero implication is therefore iff, with a
positive diagonal witness whenever closure fails. The staircase/triangle-flip
argument was checked for all n, including the distinct n5 incidence pattern
and n>=6 maximum. The finite characteristic polynomial, scalar Hankel
determinant `1/2059408867230`, and full Krylov ranks were recomputed in an
independent implementation. The minimal dimension statement is explicitly a
moment representation including F_0, not an unqualified transfer-function
degree statement.

**Pairing and wiring.** Root's labelled-orbit enumeration and the producer's
fixed-spine enumeration independently give the same minima and histograms.
[The independent review](metagraph/PAIRING_REVIEW.md) checks every finite
group implication and the physical converse/reflection lift. It also records
the final repaired code transport: all ten metadata rows, twelve ordinary
classes, and all100 entries of the order-five adjacency matrix agree after
the two canonical gauges are reconciled. Vertex rigidity leaves dart/parallel
edge permutations fixing vertices outside the conclusion.

The later general lift lemma is separately reviewed in
[PAIRING_LIFT_REVIEW.md](metagraph/PAIRING_LIFT_REVIEW.md). Its semidirect
product needs a choice of names in each double fibre; its polynomial counts
identity as an involution and restricts the base group to preserve fibre
sizes and any requested H/colors. Independent brute-force controls include
exchanged double fibres, where the two coherent lift choices are essential.

**Signed actions and eta products.** [Independent proof review](audit_sporadic/REPORT.md)
checks the signed pair-stabilizer iff, the composite odd-cycle controls,
integral-root necessity/sufficiency, support depth and valuation equality
boundary, and the modular-cusp obstruction. Root checked the same arguments,
including the difference between a divisibility filter and a coprimality
filter. The title was narrowed to fixed finite-matrix recurrences. Primary
eta inputs were reopened by the independent reviewer and are linked in that
report. No moonshine theorem is claimed as new.

The later recovered lattice section is inherited, and was separately checked
by root against its referenced construction: the Golay sum guard is retained;
the null vectors have Gram matrix `[[0,-1],[-1,0]]`; their unit pairing gives
primitive integral splitting; the printed projection kills both inner
products. The Leech identification is cited. The 25-vertex score lift is null
iff transitive; cyclic score lifts do not define the proposed null quotient.
The review script also independently checks the displayed Gram arithmetic.

## Exact executable universes

| Program | Complete declared universe and independent controls |
|---|---|
| `pairings/check.py` | All labelled orientations through n6, partitioned into ordinary classes; DP versus literal Hamilton permutations; finite-core versus first-SCC convolution through n160 for H<=21; modular periods and explicit optimal involutions. |
| `metagraph/replay.py` | All fixed-spine tilings n3..6; physical path/orbit reconstruction of the same matrix; exact weighted walk moments through length8; triangle-twist controls through n32; 868 always-active checks. |
| `audit/check_small.py` | Independent adjacent-transposition orbit partition of all labelled orientations n3..6; all literal Hamilton paths; native matrix, colored incidence, Gram defect, minimal scalar cubic and Krylov ranks; every H-preserving converse-equivariant involution and every degree-compatible graph vertex permutation n3..5; 137384 always-active checks. |
| `sporadic/check.py` | Signed orbital and independent GF2 constraints for odd m3..13; all512 cyclic pair completions at m5; eta products versus derivative recurrence through degree128, including roots; complete Paley7/23 row-code word sets. |
| `audit_sporadic/check.py` | Independent eta Euler-division versus recurrence through degree80; permitted roots and first forbidden coefficient; literal signed action checks at m3,5 and stabilizers through m9; all512 cyclic pair completions; cusp controls. |
| `sporadic/check_other_lanes.py` | P1 binomial controls R1..48,k1..5; 2040 synthetic reduced-pole cases for P2; Lorentzian Gram arithmetic; 382887 always-active checks. Its prose root-code review predates the final pairing block, which is separately reviewed above. |
| `metagraph/pairing_lift_check.py` | General lift formula versus all ordinary permutations in eleven synthetic cases with at most six ordinary objects; subgroup and double-orbit controls; n3..5 graph/abstract pairing polynomials; 22 always-active checks. |

Each program has matching normal and optimized-Python receipts. Explicit
exceptions retain the gates under `-O`; producer code is not imported by the
independent checker. Matching receipts are finite implementation controls,
not proofs of unbounded quantifiers. The figure has a separate standard-library
generator which derives and validates the ordinary converse pairs; the final
render was visually inspected for legible labels and unclipped content.

Run from the repository root:

```sh
python3 -I -B 05-knowledge/results/joint_metagraph_dyadic_20261113/verify_packet.py
python3 -I -B 05-knowledge/results/joint_metagraph_dyadic_20261113/verify_packet.py --replay
```

The first command verifies the closed packet manifest, promotion inverses,
review/source bindings and stored normal/optimized equality. The second also
reruns all seven programs in both modes and compares their exact output bytes.
It does not recreate the historical external source documents.

## Repairs and status lineage

Two prepublication mistakes were repaired and logged in `01-canon/MISTAKES.md`.
Raw dart cardinalities need automorphism/converse factors, but those common
factors cancel in the normalized merged-node law from uniform physical path
sampling. Separately, a copied score-sorted canonical code was captioned as an
unrestricted minimum. The first pairing review missed that gauge mismatch;
independent figure construction exposed it and every affected label was
transported. Neither repair changes the final incidence counts or theorem.

[PROMOTIONS.json](PROMOTIONS.json) records exact before/after hashes and
reversible replacements for all three proof candidates. Only scoped status
lines and the overbroad eta title were changed at promotion; their substantive
proofs remain the reviewed texts. Older review hashes and pre-audit manifests
are preserved as historical receipts, not asserted to be final file hashes.
The closed root manifest is the final byte authority. Any later substantive
edit requires a new audit and manifest.
