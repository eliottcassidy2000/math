# The H533 fixed-port obstruction is escaped at exactly three cycles

Session197. **Explicit positive construction and a sharp cycle threshold within a stated class.** This packet is separate from the reviewed all-order C<=2 classification in `FEEDBACK533.md`.

**Theorem.** Among strong tournaments with H=533 and feedback vertex number at most two, the minimum Hamilton-cycle count is exactly three. A ten-vertex witness has feedback number exactly two and is built while retaining the same two initial ports throughout.

Use outgoing row masks

    [2,268,905,881,999,967,903,779,513,3].

Its three Hamilton cycles, written from vertex0, are exactly

    (0,1,2,3,4,5,6,7,8,9),
    (0,1,2,7,3,4,5,6,8,9),
    (0,1,3,4,5,6,2,7,8,9).

In particular the graph is strong. The adjacent verifier obtains H=533 independently by subset endpoint DP and by literal directed-word recursion, and obtains the three displayed cycles both from anchored cycle DP and from literal recursion. These are exact counts of the displayed finite graph, not extrapolations.

Delete vertices1,3. The remaining vertices are transitive in the order

    4,5,6,2,7,8,9,0.

Every single-vertex deletion leaves a directed triangle; the receipt gives ten explicit triangle certificates. Thus the feedback number is exactly two.

For a physical fixed-port construction, begin with the edge1->3, keep the ordered ports[1,3], and append vertices in the displayed transitive order. At each birth, a mask indicates which current ports the newborn beats; it loses to every older nonport. The sequence of masks is

    1,1,1,2,3,0,1,1.

All newborns are forgotten as ports. The verifier checks all45 arcs of the resulting graph against the displayed rows. No root replacement, equal-count rebuilding, or change of the two exposed vertices is used.

Finally, the separately proved all-order theorem `FEEDBACK533.md` states that every strong H533 tournament with C<=2 has feedback number exactly three. It excludes every smaller cycle count in the class considered here, at every order. Combined with this explicit C=3 witness, it proves the sharp minimum. This argument does not classify all H533 tournaments with C=3, does not minimize their order, and makes no claim about a universal fixed-port scalar compiler.

## Reproduction

Run `python3 verify_feedback533_c3.py` and `python3 -O verify_feedback533_c3.py`. Both write the identical `feedback533_c3_certificate.json`. Only the Python standard library is used. The exact hash of the separate negative theorem is bound in the receipt; its published proof and independent review are at commit231bda1d4b2f57506bfcd6a2b2b3decb5f5f089d, under `session197/lowcycle/feedback533/`.
