# LRC14 n=4 Einheit/Erdos-870 tournament models

Representation abundance is not minimality. That is the transferable lesson
from the Erdos-870 negative answer, and the n=4 tournament toy model makes the
same distinction visible without asymptotics.

After origin claimed HYP-3145 for the broader Erdos-870 filler-core
interface, HYP-3146 for the shift-package/canary policy, HYP-3147 for the
n=3 edge-flip kernel, and HYP-3148 for the live-core deletability audit, this
packet became the prompt-exact HYP-3199 refinement that keeps the named
`a,b,c,x,y` section and deletion audit visible.

The fixed-Hamiltonian-path chart is a cover. It is useful because it exposes
the tiling coordinates `a,b,c`, but the score-class quotient is not a group
quotient: the `S` class has five representatives, and products such as `S*S`
depend on the hidden representative. The table looks algebraic only because
we are reading a slice of the cube.

The partial-score `(0,1,1,2)` chart is different. Freeing `x=a` and `y=b`
after fixing `c` gives an exact Einheit/Klein-four section:

```text
T = E
+ = x
- = y
S = x+y
```

This is the chart that deserves to feed proof packets. It preserves the two
essential order-two coordinates, while the fixed-path `c` coordinate is
deletable for class reachability. The right proof-sidecar is therefore not
"how many ways can I represent this class?", but "which coordinates are
essential after the quotient has been chosen?".

The analogy to the Erdos-870 proof architecture is useful as a discipline, not
as a theorem transfer. In the formalized additive-basis result, a large
representation count is compatible with no minimal subbasis; the construction
uses an order-two source and filler/clustered structure depending on `k`. In
this n=4 chart, `x,y` play the order-two core, `c` plays filler, and the large
`S` fiber is the warning sign. If LRC14 packets only record the abundant fiber,
they can promote the wrong coordinate.

Operational rule for the LRC proof attempt: never promote score class,
fixed-path tiling count, or A000568-like class abundance without a
minimality/deletion sidecar and a quotient-congruence audit. Add these fields
to the edge-witness packet before using the signal:

```text
n4_fiber_multiplicity_by_class
n4_quotient_congruence_defect
n4_einheit_torsor_status
n4_deletable_arc_coordinate
n4_minimal_chart_status
erdos870_minimality_sidecar_status
```

Tournament Analysis should use model/proof carriers as vertices, not runners,
raw arcs, or score labels. The n=4 scout's selected carrier path starts with
`minimality_deletion_sidecar`, then `einheit_xy_exact_chart`, then the
`erdos870_abundance_nonminimality_gate`, and only later the raw tiling cover
or score-sequence shadow. That is the ordering the larger LRC packet grammar
should copy.

Continuation note: the next useful compression is explicit:

```text
x = a OR c
y = b OR c
S = c OR (a and b)
```

This tiny monotone circuit maps the three-flip tiling cover onto the two-bit
Einheit chart.  It is class-preserving and lossy: it forgets whether `S` came
from the apex/canary edge `c` or from the live pair `a,b`.  That makes the
projected `c` action a transformation monoid rather than a `V4` action on
score classes; on the canonical table it sends `T,+,-` to `S` and sends the
canonical `S` representative back to `T`.

The n=3 edge kernel is the smaller model to keep beside it.  With `C3` cyclic
and `T3` transitive, the one-edge flip kernel is
`C3 -> 3*T3`, `T3 -> 1*C3 + 2*T3`, so the random transition matrix is
`[[0,1],[1/3,2/3]]`, with stationary distribution `(1/4,3/4)` and
eigenvalue `-1/3`.  Keep the Worpitzky `1,4,1` degree-three descent word as a
separate sidecar, not as the kernel itself.

For the broader LRC route, the new signal ledger is: Lee-Yang/Pascal cap with
de Moivre-Laplace bulk, `phi4` off-circle dip, reciprocal radius balance
`q0=q6*R^6`, the k=8 bimodality functional
`L_y=p0+p6+(1/10)*p3`, bounded-core degree ceiling `4`, the Abel-Ruffini
quintic-wall warning, ear/Omega recursive witnesses, and a
Newton/Maclaurin quartic moment test at the arithmetic progression section.
Each item is only useful if it preserves a witness sidecar rather than
becoming another scalar abundance score.
