# Norm-seven descent: hidden sheets, arithmetic moments, and one-coordinate stabilization

2026-10-26. Research continuation from `983768efde84`. This packet is under
construction; results are not promoted until the included audit and replays pass.

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
