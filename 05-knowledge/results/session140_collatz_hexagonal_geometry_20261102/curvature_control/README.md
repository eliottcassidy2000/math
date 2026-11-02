# Curvature control for binary branching

Read `PROOF.md` for the all-level argument and scope restrictions. This
package is a bounded abstract-tree control; it does not construct a Collatz
or regular hyperbolic-hexagonal model.

The exact rational checker passes in normal and optimized Python modes:

```sh
python verify.py
python -O verify.py
```

Both saved outputs are identical. The checked universe comprises 255
vertices, 254 edges, 32,385 vertex pairs, and 32,131 edge pairs. The infinite
conclusion comes from the proof, not this finite enumeration.

`draw.py` uses installed Matplotlib only to produce the scientific SVG and
PNG. The final PNG was visually inspected on 2026-11-02: no crossing, clipped
label, or text overlap was found. Exact vertex coordinates and circle centers
are in the source; display samples are not mathematical certificates.

`PUBLIC_FILES.json` is the intended portable file list; `SHA256SUMS` covers
those files. Generated `.mplconfig` font caches are not part of the deliverable.
No Git or remote changes were performed by the construction agent.
