# Exact validation record

Both frozen-receipt checkers pass in ordinary and optimized Python.

- native_check.py:703,231 explicit conditions per replay, including receipt equality. Covers all1,140 cheap triples against8 real repair sources;703,230 conditions precede receipt comparison.
- Native relational route controls:672 routes,6,188 events and6,860 arrivals. The largest observed interior event is14 slots/38 clocks; the theorem's20/52 bound is proved, not inferred from that sample.
- Native boundary controls:160 direct accepted events,136 nonidentity, observed maxima20 slots/66 clocks. The universal22/76 boundary bound follows from exact source costs and the fixed-grade proof.
- theory_check.py:2,584,306 explicit conditions per replay, including receipt equality. Covers every255 nonempty ternary allowed table for the generic grade/cost tests, exact small clipped interval gaps and grade distances, full-boundary descent, actual native equal-target interpretation, and hostile relational controls.
- Every replay input is pinned to commit c09b642a629b7d12b529b88c2697abaee78e73e7. Inputs are portable copies; historical Markdown is deliberately archived as text.
- This is not an independent external audit. Independent reviews can use reproduce.py and the included source closure.

The first native run wrote its deterministic receipt; subsequent normal and -O runs compared it exactly. Neither checker uses Python assertions. No Git write, commit or push was performed by this worker.
