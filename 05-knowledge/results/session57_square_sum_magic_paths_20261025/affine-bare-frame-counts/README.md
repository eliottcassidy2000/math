# All unique-sum frames over one affine square plane

The exact classification is 30(q−9) disjoint unordered frames: 6(q−9) AP/AP and 24(q−9) AP/non-AP. Every frame has at least one AP part. The 30 gauge-orbit types split six plus 24. Each orbit has q members and exactly nine forbidden gauges.

THEOREM.md gives the proof, orientation convention, square-convention limits, and corrected path interpretation: 1080(q−9) distinct unoriented six-vertex paths, and the same number of unordered two-path covers, within this prescribed-cross-sum class. This does not count every path or biclique of the ambient square-sum graph.

Run with the standard library:

    python verify.py
    python -O verify.py

Both runs verify the q81, q729 and q2187 pinned witness planes. The q81 run exhausts every ambient first triple and every viable second triple, explicitly checks all 30 gauge orbits, and deduplicates all 77,760 six-vertex paths. The two larger fields have 21,600 and 65,340 disjoint frames. All guards survive optimized Python.

The original AP-only packet remains frozen. This is a separate supplement; section 5 corrects its overly cautious path-count wording and clarifies that both q81 and q729 witnesses contain zero. SOURCE_PINS.json and SHA256SUMS identify the exact local inputs and outputs. No repository or remote was modified.
