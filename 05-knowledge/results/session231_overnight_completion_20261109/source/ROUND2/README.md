# Round2: two-clock bank rigidity

`PROOF.md` gives the all-order candidate theorem for arbitrary fixed source context U beta^s W gamma^s S and every physical odd Hamilton coefficient a: a genuinely separated two-clock pure-dyadic ray forces a=1. This is scoped beyond Session191's one-clock theorem, not a completeness result. Independent audit pending.

The positive mirrored-word ray family is inherited from `square_sum_spectrum_overnight_20261031e/CENTERED-FANS.md`. The closed half-word formulas and `decode.py` give an exact scalar-only square-root candidate test and literal strong graph constructor for that existing family. Their evidence is separated from global F2 occupancy.

Reproduce:

    python3 -I -B verify.py > replay.json
    cmp replay.json CHECK.json
    python3 -I -B -O verify.py > replay.optimized.json
    cmp replay.optimized.json CHECK.optimized.json

Discovery only: `rays.py`, `RAYS.json`, `DISCOVERY.txt` describe the 61,504-template exploratory search; its positive rows mostly reproduce elementary one-block rays. `three_clock_probe.py` and `THREE_CLOCK_PROBE.json` test 10,976 bounded three-clock templates against actual a>1 bank rows and find no genuinely separated pure-dyadic ray. This is a finite negative probe, NOT an all-order three-clock obstruction. The main theorem and decoder do not depend on either search.
