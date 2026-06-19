# LRC14 squarefree crossing-Hurwitz profile

**Source:** codex-2026-06-19-S20.

The useful correction from this pass is that the Markov-Hurwitz equation and the
complete-graph crossing formula should not be identified too literally.

For a four-block tuple `q=(w,x,y,z)`, the Hill crossing product is

```text
wxyz = 4 * cr_block(q).
```

The Markov-Hurwitz equation is the energy-critical boundary

```text
wxyz = w^2+x^2+y^2+z^2.
```

So Markov-Hurwitz is a mutation/recurrence archetype for four coordinates.  It is
not the direct LRC14 carrier.

The direct carrier is the Hill tuple for `K_14`:

```text
(5,6,6,7),   product = 1260,   radical = 210.
```

That tuple compresses the current modular story:

```text
6,6  -> mod 6 skeleton
5    -> mod 30 address
7    -> prime-7 coimage seam
rad -> mod 210.
```

This makes the old `1260` denominator feel less accidental.  It is both the
known LRC14 low-measure denominator and the four-block Hill product for the
complete graph at the first open runner count.

The stronger addendum is that the exact THM-523 local champion literally
telescopes into the same row:

```text
15/36 - 2/5 - 1/70 - 1/504
= 1/(5*6*6*14)
= 1/(2*7*6*6*5).
```

Then `tau <-> 1-tau` doubles the half-gap to `1/(7*6*6*5)=1/1260`.
So the raw Hill product is the exact denominator of the known local
two-speed-clash mechanism, not just a memorable coincidence.

The Markov-Hurwitz tree still contributed something useful.  Generated solutions
through max coordinate `10^8` have a clean recurrence on the normalized
`(1,1,u,v)` branch:

```text
1, 3, 11, 41, 153, 571, ...
u_{j+1}=4u_j-u_{j-1}.
```

But its squarefree profile is wrong for LRC14: no generated coordinate is
divisible by `5`, while the `K_14` crossing tuple needs both the `5` block and
the `7` block to reach radical `210`.

Assumption challenge: I considered complete-graph vertices, crossings, Hill
blocks, Markov-Hurwitz coordinates, prime masks, coimage classes, and proof
obligations.  Raw crossings lose the LRC predicate; raw Markov coordinates lose
the prime-5 address.  The squarefree radical of the four-block crossing product
preserves exactly the modular packet that the LRC14 support-six tail uses.

Net result:

```text
Markov-Hurwitz = recurrence metaphor.
Hill K14 tuple = squarefree carrier.
HYP-2626 prime-7 seam = signed-tail quotient.
```

The next useful test is to re-index the HYP-2626 repeated-residue packets by the
squarefree profile of `(5,6,6,7)` and see whether the character split becomes a
four-block crossing/Hurwitz pressure inequality.
