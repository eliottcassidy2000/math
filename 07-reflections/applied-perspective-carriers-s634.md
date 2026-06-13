# Applied Perspective Carriers: S634

The S629/S630 lesson becomes more useful once it is phrased as a compression
rule.  Do not collapse a cross-problem echo to one scalar.  Keep the perspective
or carrier channel that explains why the scalar is legal in one world and
forbidden in another.

For LRC, THM-381 says the observer is lonely exactly when the marked observer
is a source in the threshold tournament.  The converse tournament changes this
to the marked observer being a sink.  Those are not two unrelated target lists:
they are one source/sink pair under reversal.  Through `n=6`, S634 confirms the
source and sink target counts are both `U(n-1)`, and the merged source/sink
target count is again `U(n-1)`.  This is small, but it is the right small:
the quotient preserves the marked source predicate and forgets only a duplicate
obstruction orientation.

For tournament class counting, the same move gives a clean Burnside invariant.
Let `P_n` be the number of rooted perspectives, i.e. unlabelled tournaments
with a vertex root modulo automorphism.  Let `F_n` be the number of rooted
perspectives fixed by the self-converse perspective flip.  Then rooted
perspectives modulo all-edge reversal are

```text
Q_n = (P_n + F_n)/2.
```

The computed values through `n=6` are

```text
n:   1  2  3   4   5    6
P:   1  2  4  12  48  296
F:   1  0  2   0   8    0
Q:   1  1  3   6  28  148
```

This is a practical counting target.  We should compute `Q_7` next by combining
a rooted-perspective count with the S630 self-converse fixed-term audit.  It is
also a warning: counting rooted classes without quotienting by converse may be
paying twice for the same obstruction in any source/sink-symmetric problem.

For unit distance, S634 finished the exact-core audit requested by HYP-2206.
All five stored exact `n=21`, `57`-edge graph6 cores are traceable.  Therefore
the split

```text
57 = 20 spine edges + 37 bulk edges = 20 + C_hex(3)
```

is graph-real for these cores, not just inherited from a displayed Moser slab.
The first core has a very different degree profile (`18` vertices of degree
`5`, three of degree `8`), while the other four have a degree-`3` vertex and
the found Hamiltonian spines often use it as an endpoint.  That points toward
the next real unit-distance computation: endpoint-compatible ears and faithful
one-vertex extensions, not another raw edge-count beam.

Tournament Analysis in this session used proof routes as vertices:

- LRC source/sink pairing;
- rooted-converse quotient counts;
- `n=21` unit-spine bulk audit;
- `n=22` frontier extension constraints;
- raw `Phi_3` scalar matching.

The majority route tournament was transitive and ranked LRC source/sink pairing
first.  That ranking feels right: it is the least speculative, preserves the
exact LRC predicate, and immediately halves a duplicated obstruction list.  The
raw scalar route came last; S627-S630 have now made that warning stable.

The creative next step is to treat "perspective fixed under reversal" as a
general proof-obligation datatype.  In LRC it is observer-source plus threshold
labels.  In unit distance it is a Hamiltonian spine plus bulk shell.  In
tournament enumeration it is a rooted perspective fixed by a converse action.
Those are not the same object, but they rhyme in the precise way algorithms can
use: quotient only after naming the side channel that survives the quotient.
