# LRC14 Seven-Sector Net Cap

The useful move this session was to stop staring directly at `mu_{1/7}`.
The endpoint asks for a cap on the set where the cluster phases form a
`1/7`-net.  But any `1/7`-net must hit all seven fixed sectors.  That gives a
larger, cruder object, `S7(E)`, and a sufficient inequality:

`meas(S7(E)) <= cap_k`.

This is less elegant than "consecutive minimizes `mu`", but it may be easier
to prove.  The exact bounded bank had no violations, and the sector margins
were large.  The first four maxima were consecutive; `k=12,13` preferred small
perforations, but still stayed well below the caps.

The negative result matters just as much.  Local compression is false in the
simple form one wants.  It can increase `mu_{1/7}`, and it can decrease sector
cover.  So the proof cannot simply push every shape toward consecutive by
one-step gap smoothing.  If compression survives, it has to be a global
potential, a multi-step sorted-bank argument, or a relation-height descent,
not a local monotone move.

This connects cleanly to the subtorus relation-lattice work.  Sector cover is
a piecewise-constant occupancy observable on the orbit `x -> E*x`.  High
relation height should force something close to random sector words; low
height should expose short affine relations and become finite.  That is the
same shape as HYP-2599/HYP-2601, just with a harder discontinuous indicator and
a more direct sufficient cap.

Assumption challenged: I did not use runners as tournament vertices.  I tried
fixed sectors, proof routes, relation height, gap compression moves, and net
sets.  The seven-sector quotient preserves the LRC-free predicate needed here:
if the cluster leaves no empty arc longer than `1/7`, it must meet every fixed
sector.  It loses cyclic adjacency and boundary ownership, so it is probably a
coarse upper-bound machine rather than the final explanatory geometry.

The new proof shape is:

1. Replace exact net measure by fixed-sector cover measure.
2. Prove the sector cap globally, perhaps by relation-height split.
3. Feed that back into HYP-2602's union closure.
4. Return to the colored CRT layer only after the continuous cap is secure.

LRC(14) is not proved, but the target is now one step less brittle.
