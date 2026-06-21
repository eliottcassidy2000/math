# Half-arc-transitivity is the converse Z_2 made rigid

*kind-pasteur-2026-06-21, Thread C (the user's Doyle-Holt lead)*

The user pointed at the Doyle-Holt graph — "any edge maps to any other but only in one of two
ways" — and asked if that Z_2 orientation ambiguity is the tournament converse Z_2 (reverse all
arcs; half-tiling fundamental domain, THM-549). The answer is yes, and tighter than expected: the
two notions are the **same Z_2** seen one categorical level apart.

## The single sentence

> A graph is half-arc-transitive **iff** it carries an automorphism-invariant orientation D whose
> converse D^op is **not** realized by any automorphism.

This is Tutte's invariant-orientation characterization. Now read it as tournament theory:
- D is a tournament-like object (an oriented graph, here even a regular 2-in/2-out one).
- D^op is its converse.
- "Aut cannot realize D^op" is **exactly** "D is non-self-converse (NS)."

So **half-arc-transitive = vertex-transitive + NS**. The Doyle-Holt graph is the smallest carrier
of a vertex-transitive NS orientation; on the pure-tournament side, the F_21 Cayley tournament at
n=21 is the smallest vertex-transitive NS *complete* orientation. Two "smallest" objects (27 and
21), same phenomenon, both forced off the abelian/circulant world.

## Why it must leave the abelian world

The deepest part is a one-line obstruction that the project has been circling for months from the
other side. On any abelian group, **inversion x -> -x is an automorphism**, and it reverses every
orientation through the identity. So:
- every edge-transitive Cayley graph of an abelian group is arc-transitive (Chen-Quimpo) — no
  half-arc graphs there;
- every **circulant tournament is self-converse** — the converse Z_2 is realized by i -> -i.

These are the *same fact*. It is why the project's SC SPINE (circulants, Paley) is "too symmetric"
and why the genuine Z_2-rigidity — the NS SEA of the merged metagraph G_n/Z_2 — only appears once
you leave the cyclic world. The Doyle-Holt graph is the canonical external witness that this
boundary is real and that 27 (resp. 21) is exactly where it first bites.

## The resonance

Three times now the same Z_2 has surfaced: (a) the half-tiling as fundamental domain of T<->T^op
(THM-549); (b) the SPINE/SEA split of G_n/Z_2 (SC = converse realized, NS = converse rigid); (c)
the half-arc-transitive graph (arcs split into the two converse-orbits). In every incarnation the
*half* is the fundamental domain of an unrealizable involution, and the obstruction to halving is
abelian-ness. The mathematics keeps drawing the same line — "is the converse a relabelling?" — and
the answer partitions every world we look at into a symmetric spine and a rigid sea.

The open thread: the Holt graph's invariant orientation D is a **regular partial tournament** with
its own OCF / Hamiltonian-path count. If D's H is computable and its half (under the unrealized
converse) tiles the way THM-549 predicts, the analogy would stop being a dictionary and start being
a theorem.
