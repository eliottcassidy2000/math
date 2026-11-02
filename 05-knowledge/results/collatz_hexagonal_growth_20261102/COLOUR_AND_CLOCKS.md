# Functional graphs, palette transport, and arithmetic clocks

**PROVED elementary deductions and recovered results, as labelled.**
No Hadwiger--Nelson bound or universal Collatz convergence result follows.

## 1. The ordinary Collatz arrow graph is exactly three-colourable

For any function f on a countable set, forget orientations, delete loops,
and merge opposite arrows into a single undirected edge. Every component
of this simple graph has at most one cycle. Indeed, a finite connected
subgraph with two independent cycles would have more edges than vertices,
but its edge orientations would require more than one outgoing edge at
some vertex. Any two cycles in one infinite component would already lie
in such a finite connected subgraph.

Trees are two-colourable; a unicyclic component is three-colourable by
colouring its cycle and extending alternately along its attached trees.
Thus every such functional simple graph is three-colourable. For ordinary
positive Collatz the triangle1,4,2 is present, so its chromatic number is
exactly3, unconditionally and independently of convergence.

For a finite tree component of v vertices its chromatic polynomial is
q(q-1)^(v-1). For a component with a unique simple cycle of length ell>=3,
it is

    (q-1)^(v-ell)[(q-1)^ell+(-1)^ell(q-1)].           (7)

Multiply across components. This gives a complete enumeration formula
for this particular colouring problem. verify.py independently checks
all total maps on1..5 labelled vertices, evaluating(7) against literal
three-colourings. Loops are deleted and opposite arrows merged here;
retaining loops as colouring edges would change the question.

**Hadwiger--Nelson boundary.** Its graph includes EVERY pair of embedded
points at Euclidean distance one. Such edges may cross. A three-colourable
selected arrow graph therefore cannot by itself force a larger plane
chromatic lower bound. Adding actual distance-one edges may change the
graph, but each added edge needs a geometric certificate. The old
[cube/flower/spindle comparison](../hadwiger_six_seven_eight_20261025/metric_spindle/RESULT.md)
already illustrates this distinction: the seven-image cube fold has
chromatic number3, while the Moser spindle on seven points has4.

## 2. Complete palette gluing on a functional component

This is an application of the recovered
[Session62 local-coefficient theorem](../session62_agent_theme_bridges_20261025/holonomy/TWISTED_PAIR_LIFTS.md)
and [the C2 times S4 palette reduction](../hadwiger_six_seven_eight_20261025/fano_eight/RESULT.md).
It is not a new general holonomy theorem.

Now keep one edge-cell for each native arrow, including loops and opposite
arrows as distinct cells. Equip each edge with an invertible transport on
a finite frame set X; traversal backwards uses its inverse. This is a
different carrier from the simple colouring graph in section1. Each weak
component still has at most one independent cycle, including a loop or
two-edge cycle when present.

On a component with no cycle, choosing a frame at any reference vertex
extends uniquely along every finite path: there are exactly |X| sections,
even for an infinite component. On a component with one cycle, let h be
the product of transports around it, carried to the reference vertex.
There are exactly |Fix(h)| sections. Trees attached to the cycle introduce
no additional obstruction. These statements follow from uniqueness of
paths after removing one cycle edge; no averaging is involved.

For three pairs of six palette labels, X=F2^3 and h(x)=Px+v, P in S3.
The full classification is

|Pair permutation P|Required gain v|Number of sections|
|---|---|---:|
|identity|v=0|8|
|transposition(i j)|v in{0,e_i+e_j}|4|
|three-cycle|v_1+v_2+v_3=0|2|
|any other gain in the same row|not satisfied|0|

This is the equation(I-P)x=v: its kernel has dimension equal to the
number of cycles of P, and its image gives the displayed conditions.
Over all48 affine transports, the fixed-frame histogram is

    zero:33, two:8, four:6, eight:1.                   (8)

Completion parity accepts24 transports, nine of them falsely. In the
recovered decomposition W=C2 times S4, parity discards the C2 factor;
the nine false positives are exactly the nine derangements of the four
even frames. The other15 parity-preserving permutations fix a frame.
This explains that count without identifying it with an unrelated family
of fifteen objects. Only in the three-cycle row does parity alone suffice.

This completely solves a fixed transport problem on a functional graph.
It says nothing about whether each arithmetic component reaches the
designated cycle: an infinite directed ray with identity transports has
all frames and never returns. The5n+1 cycle in GEOMETRY.md supplies a
second concrete hostile to interpreting frame compatibility as convergence.

## 3. Finite residue pictures lose a digit per inverse step

For the recovered SF5 inverse maps j_a(n)=(2^a n-1)/3, current legal
branch choices depend on n mod9. Nevertheless n mod9 alone does not
determine the next mod9 state. At the same legal exponent a=2,

    19=37 mod9,       j_2(19)=25, j_2(37)=49,
    25=7 mod9,        49=4 mod9.

Both parents and children are valid positive odd units, and both children
are larger than their parents. More generally for every k>=2,

    n=19, n'=19+2*3^k
    j_2(n')-j_2(n)=8*3^(k-1),

which is nonzero modulo3^k, although the parents agree modulo3^k.
The same legality guard is preserved because they agree modulo9.
Thus no fixed quotient n mod3^k, k>=2, is closed under even this one
legal inverse operation. This is an exact observer limitation, not a
complexity lower bound or an obstruction to all finite controllers.

For a fixed legal exponent word of length h, composition has the form

    j_word(n)=(2^A n-b)/3^h,                          (9)

where A is the sum of exponents and b an explicit integer. A source
congruence modulo3^(h+k) therefore determines its terminal state modulo
3^k. Among two legal starts, terminal congruence modulo3^k is equivalent
to source congruence modulo3^(h+k), since2^A is a3-adic unit.
This is an exact precision statement on each legal word domain, not a
claim that every word is legal or every source belongs to SF5.

The SF5 choice of the word itself is determined for h steps by the root
modulo3^(h+1), together with branch selections: the j-th choice tests a
mod9 state after j divisions by3. The labels and growing precision are
the arithmetic data a six-state residue hexagon omits.

This is the same denominator mechanism as the prior
[guarded FRACTRAN precision theorem](../stability_priority_moduli_20261102/FRACTRAN_PRECISION.md),
applied to this inverse family. The native time is separately preserved:
an exponent a costs a+1 ordinary steps and a shortcut steps. For example
5 has SF5 children13(a3) and53(a5); their ordinary returns cost4 and6.
Macro-level parity alone therefore cannot recover native-clock parity.

## 4. What a decreasing rank does prove

For a map f and a forward-invariant terminal set A, every orbit reaches
A iff there exists rho:X->N with rho=0 on A and rho(f(x))<rho(x) off A.
Necessity uses the first-hit time; sufficiency uses well-foundedness of N.
This is an exact equivalence, not a construction of a global Collatz rank.

SF5 supplies such a rank on a decidable proper subtree: its address depth
decreases by one under U until5, and its accumulated native clock decreases
by one along each expanded C edge. The inherited recognizer terminates
because accepted macro-parents are smaller positive integers. The known
convergent source7 is rejected by that recognizer. Thus a sound rank plus
unboundedly many certified sources still leaves a genuine coverage problem.

The distinction also appears in energy. On a binary tree, split unit flow
equally: the total unit-resistance energy is sum_(h>=1)2^h(2^-h)^2=1.
But under its fair boundary measure, the evaluation norm for a chosen
length-h cylinder is2^h and diverges along a single infinite address.
Finite transport energy is not finite source-evaluation energy. This is
the precise boundary highlighted by SF4 and the deleted-atom hostile in
[the source-floor packet](../collatz_source_refinement_floor_20261005.md).
