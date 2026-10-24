# Structural regular-eleven nonrigid correction

The complete nonrigid regular-eleven spectrum is4 zero-error,7 error22,5 error66,
and2 error88 classes. Error44 is impossible. If P10 counts almost-regular
ten-vertex ordinary classes, then

    a_regular11=(P10+120)/11
    E_regular11=10P10-120
    a11=(U11-10P10-19360)/11

P10 and U11 are not computed here. The proof exhausts full automorphism orbit
partitions, not arbitrary equitable presentations; the full ordinary groups
and pairwise isomorphisms are checked on bounded candidate families.

Run python3 check.py and python3 -O check.py; their JSON outputs agree exactly.
THEOREM.md supplies the exhaustive coverage proof. SOURCE_PINS.json preserves
the independent tournament primitives and inherited nonregular/SCC constants.

Reproducibility:

    python3 check.py > normal.json
    python3 -O check.py > optimized.json
    cmp normal.json optimized.json
    python3 verify_manifest.py
    python3 replay.py

The manifest-only replay runs from an unrelated directory and repeats12
semantic mutation rejections under normal/-O, each also rejected by the byte
manifest. No assertion-dependent validation is used.
