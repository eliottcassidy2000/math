# Uniform guarded sibling transport

**PROVED:** the all-clock templates below, with six native complements,
peak order two above the source, and common path frames at every untouched
old vertex. **FINITE-EXACT:** literal replay at clocks 1 through 32, and fresh
endpoint/phase decoders at clocks 1 through 6. **OPEN:** next-letter 24/44
realization, output-language closure, and the full Hamiltonian spectrum.

## Inheritance and the precise improvement

The closest mechanism is the one-helper arbitrary-arc compiler in
[`adaptive_metagraph_reformulation_20261103/METAGRAPH.md`, §5](../../adaptive_metagraph_reformulation_20261103/METAGRAPH.md).
It searches either orientation's tournament for an adjacent guarded arc and
uses an inverse gadget when needed. Its working application is
[`recursive_seed_switching_20261103/transport/NOTE.md`](../../recursive_seed_switching_20261103/transport/NOTE.md).
We do not claim a new arbitrary-arc reachability theorem.

The earlier entry family and its three literal sibling flips are in
[`TWO_DIGIT_ENTRY.md`](../../multiport_defect_grammar_20261103/transport/TWO_DIGIT_ENTRY.md).
The new result supplies a bounded native word valid uniformly in the
unbounded prefix length, with a stronger common-frame statement. The hostile
example is the long arc 0→2 of the transitive triangle: no Hamilton word in
that source contains it adjacently. Thus searching only the original graph
does not prove arbitrary-arc guarding. The corrected near miss is that the
ordinary append-helper gadget requires the arc's second endpoint not to be
the last path vertex, and may change the old last vertex's frame. A terminal
pad resolves both issues. The least-used sidecar is the actual named pair
of guarded words, before passing to any class or count quotient.

The live concepts are marked component signatures / homogeneous prefix
inflation / native guarded words / joint endpoint response / phase decoding.
The map is an actual named graph rewrite; its invariant is every old arc
except the specified pair at each completed gadget. Scalar H, strongness,
and phase are not required to survive its native intermediate states.

## A two-helper common-frame refinement

Let P be an actual directed Hamilton word of T containing a→b consecutively.
Add a fresh terminal pad z, with every old vertex beating z. Add a fresh
helper h, which beats only a among all currently present vertices. Write
P=L,a,b,R. Use the two native guards

    P1 = (L,a,b,R,z,h),
    P2 = (L,a,h,b,R,z).

A native move retains the consecutive path arcs and reverses every other
arc. P1 is directed in the birth graph. After its move, a→h and h→b are
directed because both corresponding arcs were initially reversed and were
not consecutive in P1; every other consecutive arc of P2 was retained by
P1. Hence P2 is a valid guard in the actual intermediate graph.

On pairs of old vertices, the path edge sets differ only at {a,b}. After
both complements and deletion of h,z, precisely that old pair reverses.
For every old vertex other than a,b, its predecessor and successor agree
between P1 and P2. This includes both old path endpoints and works also
at order two. No old name is deleted. Peak order is |T|+2.

The pad itself is not part of a strong intermediate graph; no strongness
claim is made for helper births or native intermediate states. This is
stabilized native reachability, not a fixed-order black-only detour.

## Why a finite quotient guard works for every clock

Suppose a source has a nonempty transitive block A=(a1→…→at) whose exterior
incidences are homogeneous. Contract A to one marker α. Any directed
Hamilton word of the quotient lifts to one in the source by replacing α
with (a1,…,at). Every old exterior adjacency of the word remains valid,
and the inserted internal adjacencies follow the transitive order.

If the requested arc is exterior to A, the two-helper construction above
therefore lifts too. This argument needs homogeneous incidence only in
the source of the gadget. The first native move can alter nonconsecutive
internal prefix arcs; the explicit second guard remains legal, and the
second move restores every prefix arc. Every prefix frame is common.

Use distinct formal names a1,a2,… for the all-order family, disjoint from
0,…,4,r and the helpers. Finite receipts use names 10,…,10+t−1 for t≤32
and r=99; this finite numeric convention is not extended across a collision
with r or a helper.

## The three explicit all-clock repairs

Let r→A→B, where A=TT_t and t≥1. A vertex of B beats r exactly when it lies
in W. The initial suffix B is TT2[0,1]→C3[2,3,4], with 2→3→4→2 and
W={0,1,2}. The three repairs and their quotient guards are

| Old arc to reverse | Guard in the current contracted source |
|---|---|
| r→4 | (α,0,1,r,4,2,3) |
| 1→r | (α,0,1,r,3,4,2) |
| 4→2 | (α,0,1,4,2,r,3) |

Inspection of the declared arcs proves each guard directed. Inflate α and
apply the padded gadget. After all three repairs, B is the transitive word
0→1→2→3→4 and W={0,2,4}, exactly the previously constructed digit-4 target.
There are six native complements and six temporary births/deletions,
independent of t, with peak order t+8. Every original name survives.

For the seven-coordinate component signature, write h=H(B), f=F_W and
j=J_W. The inherited exact attachment formula is

    H(M(TT_t;B,W)) = h−j+2^t f,
    E_r(M,2)−H(M) = j.

Literal five-vertex signatures, independently recomputed by the checker,
give every completed stage:

| Stage | h | f | j | Full H |
|---|---:|---:|---:|---|
| Initial digit-2 suffix | 3 | 21 | 16 | −13+21·2^t |
| Add 4 to W | 3 | 33 | 28 | −25+33·2^t |
| Remove 1 from W | 3 | 27 | 14 | −11+27·2^t |
| Break the terminal triangle | 1 | 21 | 12 | −11+21·2^t |

All completed graphs are strong: the terminal strong component of B meets
W at every stage. Every B vertex can reach a port, then r, then A; A reaches
all B. Thus the sibling map increases H by 2 and decreases E_r(2) by 2,
while retaining f=21 at its endpoints. It does not multiply H by three.

## Named sidecars and finite validation

The receipts save complete physical-name lists and matrices at all completed
states, and both actual native guards. For t=1,…,6 they freshly compute all
30 ordered two-component endpoint joins J(a,b) on ports {r,0,1,2,3,4} and
the golden/F49/F343 charts at each port. This pair matrix is a useful
readout; it is not claimed to replace all higher joint-port coordinates.
The full named graph remains available for those queries.

The exact chart definitions and inherited literal 81/27 control are in
[`PHASE_CONTROLS.md`](../../multiport_defect_grammar_20261103/transport/PHASE_CONTROLS.md).
Here X=(S_r(2)−e_r)/2, Y=(S_r(2)+e_r)/2. The decoders are

    F49 plus charge: (a,b,c)=(X+4(H−Y), Y, H−Y) modulo 7,
    Frobenius:       (a,b,c) ↦ (a+4b,−b,c);
    F343 vector:    (H,X,Y) modulo 7,
    Frobenius:       (H,X,Y) ↦ (H+4X+6Y,2X,3X+4Y) modulo 7.

Changing the actual graph requires fresh decoding; neither native operation
is asserted to implement field Frobenius. These new sources are not the old
81/27 pair, and we do not transfer its symmetry labels or numerical responses
to them. That inherited control was unchanged and not rerun here.

Reproduction from repository root:

    python 05-knowledge/results/c4_port_returns_20261103/native/check_native.py
    python -O 05-knowledge/results/c4_port_returns_20261103/native/check_native.py --output 05-knowledge/results/c4_port_returns_20261103/native/receipt_optimized.json

The normal and optimized receipts must be byte-identical. The finite universe
contains 96 arc gadgets / 192 native moves, plus the transitive-triangle
non-guard and the order-two/last-arc positive controls. The literal native
interpreter acts on unordered pairs; it does not infer legality from the
count identities. The all-clock result follows from the templates and the
inflation proof, not extrapolation from these 32 clocks.

## Exact remaining boundary

This gives an all-clock native bridge between the two known sibling outputs.
It does not produce the missing two-letter targets with (c,f)=(-35,63) and
(-29,63), certify a larger component-bank hit, or show either output returns
to a language carrying both ternary digits. The component-bank inverse and
the native guard compiler discharge different obligations: an accepted
actual target supplies something to compile; the compiler cannot create a
target with a requested count by itself.
