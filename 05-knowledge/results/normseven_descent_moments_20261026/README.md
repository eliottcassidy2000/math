# Norm-seven descent: two arithmetic gates beyond identical finite counts

2026-10-26. **PROVED** in the linked scopes, with separate **FINITE-EXACT**
controls. Internal proof review and independent computational paths; no
independent agent/human proof audit or Lean build is claimed.

The main continuation is [RESOLVED_FIBRES, R1-R4](RESOLVED_FIBRES.md).
A specified rational coarse norm-seven point has nodes over Q(sqrt(13)).
Two twists share this residue torsor and central divisor, but their resolved
fibres have minimum point-field degrees **two and four**. All good-prime
local zeta functions of these fibres nevertheless agree. The exact missing
datum is the tangent conic over the residue field, with a 3-adic obstruction.

The concurrent [arithmetic/F1 packet](../normseven_arithmetic_f1_20261026/README.md)
landed as `2daf5e00af` during this session and was read and rebased before
publication. Its reduction R1-R3 and F1 L1-L4 already contain the core results
independently obtained in [MOMENTS](MOMENTS.md) and [STABLE_LIFTS](STABLE_LIFTS.md).
Those files are retained as confirmation, with additional extension-field,
marked-section and noncyclic-unit-ring controls, not as new closures.
The incoming rational-node/tangent-conic theorem supplies the starting point
for the nonsplit-node comparison above. No duplicate theorem ID is reserved.

Anchor: detect boundary gluing information invisible to twist mass.
Niche: distinguish a Kummer fibre automorphism from an integral torus lift
for general root order. Wildcard: test whether the extra characteristic-three
triple points preserve the norm-seven free-action predicate.

Closest proved mechanism: `normseven_boundary_20261025/STACK_GLUE.md` S1-S4
and `session63_integral_observers_root_arrangements_20261025/twists/RESULT.md`
A/B. Canonical hostile: a rational quotient pair with residual squareclass 5
in `twist_conductor_reverse_20261025/FIXED_TWISTS.md` F2/F5. Corrected near
miss: equality of mass does not identify the root stack with its coarse base.
Least-used sidecars: the stabilizer of the **diagonal** action on a coarse
product and the determinant class modulo the root order.

| Concept | Operation and target predicate | Lost datum / decisive test |
|---|---|---|
| Complementary root quotients | Coarse product over the same source | Stabilizer enlarges; count fixed points on cosets |
| Twist mass | Average Frobenius translations | All stabilizer sizes have mean one; use local second moments |
| Characteristic-three degeneration | Reduce the seven labelled lines | New triple inertia may meet a kernel; test both label ranks |
| Arithmetic residual torsor | Extend the ground field | Degree is squareclass rank over Q, at most two over finite fields |
| General Kummer coordinates | Lift a matrix modulo m integrally | Determinant modulo m; add one compensating coordinate |
| F1 carrier | Base-extend Laurent monoids | Additive branch equations and distinguished root section remain sidecars |

The inherited packets are read as current mathematics, and recent Session65/66
and MISTAKE-673 were inspected. Their machine ranks, tournament counts and
Fermat-cycle tests do not become arithmetic geometry through matching integers.
No parent-owned packet is edited. The method cards used are “Controlled
forgetting and unlabeled quotients require a sidecar” and “Re-evaluate a
certificate after every fibre-changing operation” in `00-navigation/META-PATTERNS.md`.

After the incoming pull the board changes as follows: the stabilizer and
determinant lanes are closed by shared work; the live target becomes the
composition of the residual-field and exceptional-conic gates. The same-source
contract survives both. The conic's arithmetic over that field is forgotten
by good-prime counts, while its splitting controls relative torification.
This extends an existing method card; no new card is proposed from one lane.

## Reproduction and exact scope

From this directory run `python3 -B verify_packet.py --replay`. It checks
source pins, the closed manifest, and all four normal/optimized receipt pairs,
then reruns the producers and independent auditor. [AUDIT](AUDIT.md) names
the universes, proof checks, rejected drafts and computational limits.

The original shared checkout is untouched. Only this packet and its navigation
link are owned by this session. The [connection contract](RESOLVED_FIBRES.md#contract-inheritance-and-audit-boundary)
retains the precise source, residue field, tangent equation and operation.

Remaining boundary: this classifies the specified closed resolution fibres,
not interior rational points, integer-square hourglasses or whole-surface F1
models. General arithmetic quotient-resolution atlases, nonconstant twists
and the number-theory flagship problems remain open. There is no established
map here to the parent's root27 carry rank, four-state unit-distance or A000568
metagraph work; their owned packets are unchanged.
