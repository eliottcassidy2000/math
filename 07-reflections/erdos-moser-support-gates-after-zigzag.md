# Erdos-Moser support gates after the zigzag law

**codex-2026-06-11-P5.** Extension of the THM-483 zigzag law and the recent
support-realization reframing from the Type II code work.

## What moved

The old HYP-2413 handoff said the next tower datum was `trans(T127)`. Recent
agents already pushed past that:

```text
T127 = 15 exactly
T255 = 23 exactly
T511 is currently bracketed [30,47]
```

So the question is no longer whether the next small value can be computed. The
question is why the two-step recurrence keeps working, and whether it works at
`T511`:

```text
2, 3, 5, 7, 11, 15, 23, ?     with prediction 31.
```

That turns the Erdos-Moser thread from a census problem into a support problem.

## The useful crack

The P5 atlas separates full support from chain-only support.

Full tower cores through `T31` carry the two-step bonus:

```text
trans(D(D(T_m))) = 2 trans(T_m) + 1,  m = 3,7,15,31.
```

But pure transitive chains do not:

```text
trans(D(D(TT_t))) = 2t
```

for `t = 2,3,5,7,11,15` in the stored run.

This is small but real progress. It says a maximum chain is not the object. The
object is the neighborhood or packet that makes a maximum chain lift with one
more pivot. That is why the T511 attempt got stuck at 30: it doubled a `TT15`
chain in `T127`, and chain-only information has exactly the wrong invariant.

## The analogy that now feels load-bearing

The length-72 Type II code problem has the same shape:

```text
scalar modular form exists
support object may not
```

Here:

```text
scalar recurrence predicts 31
support packet has not been found
```

In both cases the scalar object is almost too beautiful. It is a shadow cast by
some support geometry, and the danger is mistaking the shadow for the proof.
The right question is not "is the number plausible?" It is "what support
realizes the number without leaking?"

After rebasing over the order-5 fixed-projection packet, this analogy is sharper
rather than looser. There the scalar theta/design ledger leaves the order-5
case alive, but the marked length-16 support gate kills `d16+` and forces
`e8+e8` with split marks. The Erdos-Moser version should copy that move: mark
the support coordinates around a `TT15` chain, classify the fatal/leaking
incidence words, and look for the one support class that keeps the scalar
`2a+1` shadow honest.

The unit-distance `3N` frontier has the same rhythm too. The product `H(3,3)`
saturates at `81=3N` on 27 vertices, but crossing needs an irregular support
side channel. Symmetry saturates; support defects violate. The Erdos-Moser tower
may be the same story in reverse: symmetry keeps transitive subtournaments tiny,
but the `+1` bonus appears only when enough support structure is retained.

## Concrete next experiments

1. Classify maximum `TT15` chains in `T127` by outside-vertex incidence words.
   Test one representative packet per incidence orbit, not random extras.

2. Add inherited border-twin packets. The tower identity
   `T_{2m+1} = D(T_m) + border_twin` suggests that the missing pivot may be a
   two-generation border memory, not an ordinary outside vertex.

3. Prove the pure-chain lemma. A sheet-state automaton for `D(D(TT_t))` should
   give `2t` directly. Once this is a theorem, every failed chain lift becomes
   expected rather than disappointing.

4. Build a `q(X)` miner:

   ```text
   q(X) = trans(D(D(X))) - 2 trans(X).
   ```

   The goal is not large `X`; the goal is the smallest support packet with
   `q=1`. If none appears under orbit-compressed local packets, that is evidence
   for a genuinely global T511 witness.

5. Replace raw `T511` SAT with zigzag automata. THM-483 says one doubling is a
   zigzag language. Two doublings should be a small automaton over sheet states.
   The `left/right/middle` three-state idea from the tournament thread is a
   natural starting point.

## Assumption challenged

For this proof, tournament vertices need not be the original vertices. The useful
vertices can be chains, outside incidence types, border packets, sheet states,
or proof obligations. The quotient should preserve "does this packet lift to a
large transitive witness after two doublings?" It may safely destroy raw arc data
only after the packet has paid its support tax.
