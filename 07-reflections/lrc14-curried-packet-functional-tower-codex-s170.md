# LRC14 Curried Packet Functional Tower

codex-2026-06-24-S170

The phrase "think of everything as functions; curry them" gives a surprisingly
strict version of the current LRC14 proof discipline.

The old failure mode was:

```text
row -> number
```

and then arguing as if the number remembered the row.  HYP-2990 through
HYP-2997 already pushed against this by asking what a quotient may forget.
HYP-3002 makes that operational: every time a proof route fixes one argument,
sums one coordinate, or evaluates one lane, the remaining object should still
be a function on the unforgotten fiber.

The master type is:

```text
E : S -> P(S) -> root -> lane -> fiber -> certificate -> verdict.
```

This makes several older themes click into one grammar.

The Fibonacci row is not `F_n`; it is `A(n)(k)=C(n-k-1,k)` until we decide to
sum over `k`.  Zeckendorf is not a count; it is a section choosing one legal
support after the path-conflict function is known.  Goldbach is not a count; it
is a high-entropy prime-pair fiber with local residue functions.  Fermat
polygonal coverage is not a yes/no theorem shadow; it is a bounded-budget
function whose local residues are absorbed by a fixed arity.  Farey `p+q`,
`p*q`, and powers are not equivalent mutations; they are different partial
evaluations of `F(p)(q)(e)(lane)`.

For LRC the same structure was already present:

```text
danger(S)(t)(v)
C(S)(t)
Phi(S)(k)(v|k)(atom_bank)
R(packet)(section)(grid_exit)
lost_Q(c)(y)(x,x')
```

The proof idea is to never let these collapse silently.  If a route uses a
tournament shadow, the shadow is a function of a gauge and a time/packet
choice.  If it uses a Fejer certificate, the certificate is a function of a
packet family, harmonic degree, and atom bank.  If it uses a Farey product, the
product is a function inside the exact root fiber, not a replacement for the
root.

This is a useful way to read AP/GW.  They are not merely rows with zero open
safe mass; they are rows where all currently known curried functions close
early at a boundary equality return.  K33 behaves differently: the product
partial evaluation exposes an incidence/state-lift handoff before the row can
be scalarized.  Covering rows behave differently again: the all-covered chart
is only one partial evaluation of a multi-chart boundary-moment function.

The practical next step is boring and important: packet records should carry a
`curried_call_signature`.  A future classifier should not only say
`route=K33-STATE-LIFT`; it should say something like:

```text
E(row)(K33_packet)(p=3,e=1)(product)(cross_handoff)(Fejer d159)
```

That string is not decoration.  It is the proof audit trail: every parenthesis
marks a coordinate that has been fixed, and every missing coordinate becomes a
lost-coordinate function that must be zero, reconstructed, exact, annihilated,
descended, boundary-equal, or named as residual debt.

The main creative upshot is that "functions" and "cocycles" are almost the same
warning in this repo.  A cocycle is what appears when a curried function refuses
to be independent of the coordinate we tried to forget.  The LRC14 proof should
try to show that all such refusals are already named before a strict
counterexample can be formed.
