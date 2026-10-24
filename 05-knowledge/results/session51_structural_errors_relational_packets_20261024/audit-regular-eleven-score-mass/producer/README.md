# Atlas-free regular-eleven score-mass correction

Two exact fixed-labelled-score recurrences and the independently proved18
nonrigid regular-eleven classes give:

- 1205 rigid and1223 total regular-eleven classes.
- 13333 almost-regular-ten ordinary classes.
- Regular correction133210; total order-eleven correction152690.
- a11=(U11−152690)/11, with U11 not computed here.

No almost-regular10 or regular11 atlas or target class count is an input.
This is a score-count/orbit-mass sidecar, not bare metagraph sufficiency.

Reproduce:

    python3 check.py > normal.json
    python3 -O check.py > optimized.json
    cmp normal.json optimized.json
    python3 verify_manifest.py
    python3 replay.py

The standalone proof is THEOREM.md; inherited inputs are source-pinned under
sources/. The physical-copy replay repeats normal/-O and12 semantic hostile
rejections, each also rejected by the manifest.
