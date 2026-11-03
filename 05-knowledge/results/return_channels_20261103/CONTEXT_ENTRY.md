# A context can enter the observer image without realizing the original source

Status: **REFUTED** right-cancellation of observer-state membership;
**FINITE-EXACT** complete tests in the declared box. No accepted prior theorem
is retracted; the failed implication was a candidate next-step simplification.

Let S be the states realized by65-bank words. Forward invariance says
`v in S => vP_l in S`. The converse is false even if v satisfies every
necessary physical cone condition. An exact witness is

    v=(3,0,3) not in S,
    vP_a=(3,3,3)=v(q1) in S.                           (E1)

The guarded inverse proves the first exclusion at every word length and
constructs q1 for the second state. Thus the inverse search must consider all
legal predecessors; it cannot cancel a named final matrix solely because the
output state happens to have some bank presentation.

There is an actual graph behind this example. The rows `(10,12,9,0)` describe
`C3 -> singleton`, with port3. Its physical state is(3,0,3). A singleton
one-backward-bridge join gives rows `(26,28,25,0,8)`, port4, state(3,3,3).
The bank representative q1 has that same state but **four** vertices; the
supplied extension has **five**. Literal permutation counts recover all three
states, including the marked one-backward-step defect used to obtain C.

The complete bounded-H3 word census has12 words, and q1 is the ONLY word
with state(3,3,3). Thus no five-vertex bank presentation has that state. More
intrinsically the five-vertex source has no Hamilton path starting at its
port, whereas q1 has one. [The port-guard theorem](PORT_GUARD.md) proves this
zero cannot be repaired by ANY further bank append. Equal observers must not
substitute the four-vertex object for the actual five-vertex extension.

The complete new test universe is every absent state satisfying the physical
cone with odd h<=15 and0<=x,y<=2h, followed by each of the65 letters. There
are942 such states and61230 pairs. Exactly49 outputs enter S; every output
has an explicit returned source word in `context_entry.json`. These are
observer-image entries. Most formal inputs have not been supplied with a
literal graph, and none is declared realized merely by the output's success.

The strongest survivor is forward image closure and the full guarded inverse.
A useful next question is which contexts enlarge the observer image in a
controlled way while preserving a separately supplied graph and clock. The
four ternary closure failures remain absent after our chosen bridge, so(E1)
does not supply the missing global return operation.
