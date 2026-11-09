# Session231 dense rooted fibers: deliverable map

## Main result

Start with **SHORT_INTERVALS.md**. Independently audited: for every k≥0 a single physical full interface supports k+1 consecutive odd Hamilton counts, with strong witnesses on12k+13 vertices. The explicit baseline is

    D_k=(245051264·4096^k+1602510k+1118023)/24843,
    counts D_k−2j, 0≤j≤k.

This closes the inherited missing *unbounded consecutive repair width* mechanism. It does NOT prove baseline occupancy or Chan–Pak completeness.

- `SHORT12_UNIT_LOOP_WITNESS.json`: eight-vertex core, cut124, twelve-bit X/Y words.
- `verify_short12_unit.py`, `SHORT12_CHECK*.json`: primary exact control.
- `verify_short12_independent.py`, `SHORT12_INDEPENDENT*.json`: independent six-state scalar control throughk15.
- `PLACEMENT.md`, `verify_contexts.py`, `CONTEXT_CHECK*.json`: later context recurrence, residue, lower bound and content filter. These later refinements are algebraic proof candidates pending independent audit unless separately promoted by the parent.

## Earlier positive constructions and lineage

- `UNIT_INTERVALS.md`, `UNIT_LOOP_WITNESS.json`, `UNIT_LOOP_SEEDS.json`, `verify_unit_family.py`: independently audited1188-bit unit loop, now superseded only in efficiency by SHORT12.
- `POSITIVE_FAMILY.md`, `LOOP_WITNESS.json`, `verify_family.py`: independently audited unsaturated h95 fixed-step190 unbounded fiber, plus general cycle-gain quantization.
- `baseline_and_gcd.py`, `BASELINE_GCD*.json`: all-k H95 baseline residues149,169,119 mod190 and gcd-one control; not a prime claim.
- `PROOF.md`: earlier movable-core automaton, H29 endpoint-erasure pair, H45 five-value hole, bounded extension controls, complete census and residue hostile. See status labels in that file rather than treating all discoveries as independent audit promotions.

## Complete finite censuses

`ERASE9.stats` freezes the final order9 total:191536 classes,193541720 erased interfaces, maximal run3 overall and2 unsaturated. Corrected order7/8 totals are104432/3385144. `ERASE_CLASSES7/8.*`, `ERASE7.*` and `RESIDUE6/7/8.*` retain the exact receipts. The earlier phantom-empty-map diagnostic affected only reported group totals and has been fixed in source.

The simple mod8 missing-residue conjecture is false at order8; the explicit five-value hostile is in PROOF.md. Complete small-core maximum3 does not contradict the all-order moving construction, whose deletion cores become larger than9.

## Reproduce compact positive certificates

Run inside this directory (standard-library Python only):

    python3 -B verify_short12_unit.py > SHORT12_CHECK.json
    python3 -B -O verify_short12_unit.py > SHORT12_CHECK.optimized.json
    cmp SHORT12_CHECK.json SHORT12_CHECK.optimized.json
    python3 -B verify_short12_independent.py > SHORT12_INDEPENDENT.json
    python3 -B -O verify_short12_independent.py > SHORT12_INDEPENDENT.optimized.json
    cmp SHORT12_INDEPENDENT.json SHORT12_INDEPENDENT.optimized.json
    python3 -B verify_contexts.py > CONTEXT_CHECK.json
    python3 -B baseline_and_gcd.py > BASELINE_GCD.json

`core.py` and `verify.py` provide shared exact routines. `verify.py` has a main guard, so importing its literal enumerator does not run or print its separate census.

## Regenerate large discovery inputs (do not need to publish generated files)

The inherited `classgen.c` is copied unchanged. From this directory:

    gcc -O3 classgen.c -o /tmp/dense-classgen
    /tmp/dense-classgen 9 classes9
    g++ -O3 -std=c++17 erase9.cpp -o /tmp/dense-erase9
    /tmp/dense-erase9 classes9_9.txt 64 > ERASE9.jsonl 2> ERASE9.stats
    g++ -O3 -std=c++17 saturation_bank.cpp -o /tmp/dense-saturation
    for n in 3 4 5 6 7 8 9; do /tmp/dense-saturation $n classes9_$n.txt > SATURATED$n.jsonl 2> SATURATED$n.stats; done

Generated `classes9_9.txt` (~30MB), `SATURATED9.jsonl` (~14MB), duplicated class files, binaries, temporary search traces, and Python bytecode need not be published. Preserve the generator, census/search sources, small physical witnesses, exact final receipts, and this reproduction route. The all-order short theorem depends on neither the large census files nor discovery search completeness.

## Precise next question

The intervals have real unbounded width but exponentially separated unpadded baselines. Arbitrary common singleton contexts retain the unit spacing. Prove, or find an obstruction to, placing one of these real windows over each sufficiently large target odd integer. Gcd-one, residue coverage, entropy, and the mere existence of long intervals do not discharge that selector.

## Compact bounded-minimum and diameter refinements

- `MINIMALITY_AND_DIAMETER.md`: precise core-order≤9 shortest-loop scope and an all-order O(L) one-core diameter bound; these refinements await their separate audit.
- `MINIMALITY_SEEDS.json` (~compact filtered seed list), `check_shorter_loops.py`, `SHORTER_LOOP_CHECK*.json`: default check needs no large generated bank. Regenerate the bank with the commands above, then run `python3 -B check_shorter_loops.py --regenerate` to independently recreate the seed list and repeat the check.
- `CENSUS_INPUT_HASHES.json`, `ERASE9.stats`, `ERASE9.jsonl`, `SATURATED3..9.stats`, and `ORBIT_MASS.json`: retain small hashes/count receipts while omitting generated large data.
- `short12_all_positions.py`, `SHORT12_ALL_POSITIONS.json`: bank-relative intermediate-position control. This discovery control needs regenerated SATURATED files, unlike the main standalone positive verifiers.

`FINAL_SHA256SUMS` pins the final compact proof-and-certificate packet. The older checkpoint manifests document intermediate milestones and are not substitutes for the final manifest.
