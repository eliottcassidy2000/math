# Dyadic conductor tower and the marked blue-fibre bridge

2026-10-21. **PROVED within the stated dyadic and marked local scopes;
independent review FINAL PASS.** Read the [proof](CONDUCTOR_TOWER.md),
[resolved corrections](CORRECTIONS.md), and [independent audit](review/evidence/INDEPENDENT_AUDIT.md).

## Results

The actual Eisenstein-order class groups at conductors 2,4,8 have orders
1,2,4, corresponding to idoneal integers 3,12,48. At f=2^a, a>=2, the group
is C_(2^(a-2)) x C2 with explicit generators and ideal-extension maps.
At conductor sixteen the sequence to conductor eight is nonsplit: every
lift of the specified class A has order four. The marked shore exchange at
conductor eight does not lift through this next extension.

For a finite blue parent path, a chosen zero-bit origin and shore marking
give the explicit map (e,h) -> B^e C^h. Ideal extension detects e XOR h;
the two blue lifts form its kernel. Genus characters recover the individual
bits. Two exceptional size-two strata retain distinct inactive-shore
stabilizers; neither is the diagonal conductor kernel.

This is not an identification of the full tournament metagraph or its
path-flip incidence. It supplies no electrical conductances or response map.
Conjugation is inversion, not the shore exchange. The complete global list
of 65 idoneal numbers is not assumed.

## Inheritance and correction history

The incoming [marked-fibre network precursor](../metagraph_fibre_network_20261020k/README.md)
already proves the restricted idoneal classification and marked regular
4/2/1 torsor identifications. The present additions are the literal all-level
dyadic tower, genus/blue-predicate interface, nonsplitting and symmetry-lift
obstructions. The [balanced two-port companion](../balanced_deletion_groups_ports_20261021/README.md)
retains the incidence data not supplied by the class-group quotient.

An f=40 hostile refutes extending the dyadic square-kernel statement to all
multiples of eight. The final proof is narrowed accordingly. Only the final
corrected proof, scripts, receipts and current_inputs/ snapshot are packaged
here. Superseded original inputs and their old execution receipts are not
published. The independently reviewed final-only audit records the concise correction
history and depends only on its current_inputs/ snapshot.
The final verdict and FINAL_INPUT_PINS.json control the portable replay.

## Reproduction

From this directory, with standard-library Python only:

    python3 -B verify_conductor_tower.py
    python3 -B -O verify_conductor_tower.py
    cd review/evidence
    python3 -B independent_checker.py
    python3 -B -O independent_checker.py

Producer normal/optimized output is identical across 28,754 exact checks;
independent output is identical across 46,719 checks. The finite conductor
range is 2,4,8,16,32,64,128; the all-level conclusion rests on the written
induction, not extrapolation. Independent inputs and verified primary-source
locations are retained in FINAL_INPUT_PINS.json and SOURCE_CHECKS.json.
