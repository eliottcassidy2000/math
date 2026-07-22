# Coherent subset lift versus orbit incidence across LRC, DvdK, and planar JC

The most useful connection from THM-2101--2103 is not a shared formula. It is
a fork in what one may lawfully do after replacing a structured family by one
scalar relation.

## The fork

Let a fine object contain a distinguished subset or fiber `S`, and let a
quotient remember only a sum, measure, bracket layer, or pair weight. There are
two legitimate continuations.

1. **Coherent-lift route:** retain the coordinate that decides whether all
   local equalities lift simultaneously.
2. **Orbit-incidence route:** if a transitive group moves `S`, every translate
   has the same scalar value, and a global sum is independently known, sum over
   the group. Uniform incidence may turn the subset identity into a
   contradiction.

Confusing the two routes is dangerous. A scalar equality rarely reconstructs
its subset. Conversely, when equivariance and a full-object invariant are
available, refusing to average can hide a short proof.

## LRC: pair averaging needs a coherent-fiber sidecar

THM-2126 computes every transverse pair weight from the primitive relation

```text
A g+B f+C f'=0.
```

The Fourier coefficient vanishes whenever one relation coefficient is
divisible by seven. Eight terminal characters can therefore give all 28 pair
weights `5/343` and tree total `5/49`. The exact graph is not wrong; it is
incomplete. At the guard boundary the eight terminal phases occupy one common
mod-seven fiber, and the graph has integrated that alignment away.

The hostile CRT family shows this is not repaired by scalar divisor
completeness or one-deletion primitivity: those filters coexist with a tree
margin tending to zero. Yet THM-2104 immediately escapes the same family
because all quotient coefficients lie in one `2`-, `3`-, and `5`-adic
valuation layer, while THM-2114's content blocker excludes its primitive
coefficient plane. The next carrier must remember finite-ring needles, full
guard-fiber centers, clock residues, all maximum bases, or higher
intersection multiplicity. A
tournament tie-break on equal edges destroys the missing coordinate rather
than restoring it.

## Planar JC: first-layer synchronization needs divisibility

For a proper-power top pair

```text
F=A H^m,       G=B H^n,
```

THM-2102 makes the first lower layer satisfy

```text
L=A m H^(m-1)Q-B n H^(n-1)P.
```

In a nonresonant early defect, `J(H,L)=0` forces `L=0`. That scalar
synchronization still does not lift to a common approximate root. The lift
exists exactly when

```text
H^(m-1) divides P,       H^(n-1) divides Q.
```

Thus the missing sidecar is the quotient class `[P] mod H^(m-1)` (equivalently
the matching class for `Q`). The explicit `(m,n)=(2,3)` control shows `L=0`
with nonzero quotient. Here averaging would only erase the obstruction; later
Keller bracket equations must kill or classify it.

## DvdK: transitivity makes averaging decisive

THM-2101 has the opposite structure. Vanishing positive constant terms makes
the distinguished subset `S` of small roots satisfy

```text
sum_(alpha in S) alpha^(M-1)/Phi'(alpha)=1.
```

Irreducibility makes the Galois action on all roots transitive, and the subset
identity is equivariant with invariant right side. Summing every translate of
`S` counts each root uniformly. Lagrange interpolation independently gives

```text
sum_(all roots alpha) alpha^(M-1)/Phi'(alpha)=0.
```

Uniform incidence therefore turns the information loss into the contradiction
`|Gal|=0` in characteristic zero. This works because the three necessary
ingredients are present: equivariant subset weights, a transitive finite
action, and a separately controlled full-root sum.

## Reusable decision test

Whenever a new quotient looks unexpectedly sharp, ask:

```text
What is the distinguished subset/fiber?
What simultaneous compatibility did the scalar forget?
Can that compatibility be retained as a divisibility/phase/multiplicity sidecar?
Or does a transitive action make every translate identical and the full sum impossible?
```

THM-2126 answers “retain the fiber,” THM-2102 answers “retain the divisibility
class,” and THM-2101 answers “sum the orbit.” The distinction is structural;
no LRC-to-GMC-to-JC implication is claimed.
