# HYP-2748 — The Doyle-Holt half-arc-transitive graph IS the converse Z_2 made rigid; its tournament-side incarnation is the vertex-transitive NON-self-converse tournament (smallest at n=21, F_21)

(Renumbered from HYP-2747 after a concurrent same-prompt collision: mac-mini-S14 pushed
HYP-2747 = "LRC(14) vs CJJ LP hierarchy / Lovász-theta" first. Different topic; deferred to
first-pusher per the collision protocol.)

**Status:** CONFIRMED (computationally verified; the "no half-arc Cayley graph on an abelian
group" half is a classical theorem, Chen–Quimpo)
**Source:** kind-pasteur-2026-06-21-kpswf5 (Thread C — the user's Doyle-Holt arc-flip lead)
**Thread:** C (Doyle-Holt half-arc-transitive graph <-> tournament arc-flip / converse Z_2)

---

## The question (user's lead)

The Doyle-Holt graph (27 vertices, 4-regular) is the smallest **half-arc-transitive** graph:
vertex- and edge-transitive but NOT arc-transitive — "any edge maps to any other but only in
one of two ways" = a Z_2 orientation ambiguity. The user asked whether this is the tournament
**converse Z_2** (reverse all arcs; fundamental domain = the half-tiling, THM-549) and whether
any tournament-derived graph is half-arc-transitive.

## The answer (precise)

**Half-arc-transitivity IS exactly "the converse Z_2 is not realized by the automorphism group."**
By Tutte's theorem, a half-arc-transitive graph X carries an Aut(X)-**invariant orientation** D
(an oriented graph / "partial tournament"): Aut(X) preserves D, the two arc-orbits are
{arcs of D} and {arcs of D^op}, and **no automorphism sends D -> D^op**. The map D <-> D^op is
the converse involution; "Aut cannot realize it" is the literal content of half-arc-transitivity.

This is the SAME Z_2 as the tournament converse T <-> T^op (THM-549/550), one categorical level
across (graph vs digraph).

### Three verified facts

**(I) The converse Z_2 IS realized on every circulant tournament => circulants are never
half-arc carriers.** On Z_n (abelian), inversion mu: i -> -i is always a relabelling, and it
sends a circulant tournament's arc difference (j-i) to -(j-i), i.e. T -> T^op. So **every
circulant tournament (Paley included) is SELF-CONVERSE**. Verified p=7,11,19,23: inversion maps
the QR connection set S to its complement Z\{0}\S exactly. This is the digraph shadow of the
classical theorem (**no half-arc-transitive Cayley graph on an abelian group**; equivalently
every edge-transitive Cayley graph of an abelian group is arc-transitive — Chen–Quimpo). Confirmed
computationally: Paley graphs P_5, P_13 and the wiggly metagraph = hypercube Q_m are all
**arc-transitive** (arc-orbits = 1; too symmetric). The hypercube Q_6 (n=5 wiggly) has
|Aut|=46080, arc-transitive. So the wiggly metagraph is NOT the half-arc object.

**(II) The Holt graph's invariant orientation is a rigid regular partial tournament.** Built the
Holt graph exactly as the metacirculant M(3,9) (vertices Z_3 x Z_9, layer-i intra diffs scaled by
4^i, inter-layer diff): **|V|=27, |E|=54, 4-regular, girth 5, |Aut|=54, vertex-transitive,
edge-transitive, arc-orbits=2** (half-arc-transitive). Its invariant orientation D (one arc-orbit
as "forward") is a **2-in/2-out Eulerian oriented graph** (a 4-regular tournament-like object).
Exact orbit computation: **54 automorphisms PRESERVE D, 0 REVERSE it (D -> D^op)** — the converse
Z_2 is provably unrealizable. This is the exact dictionary:
  - D  <-> a tournament T
  - D^op <-> the converse T^op
  - "Aut(H) cannot reverse D"  <->  "T is NON-self-converse (NS)"
  - the two arc-orbits  <->  the converse Z_2 fundamental-domain split = the HALF (THM-549).

**(III) The genuine tournament-side Doyle-Holt = the vertex-transitive NON-self-converse
tournament; smallest at n=21 on F_21.** A VT tournament is the digraph cousin of a half-arc
carrier when it is NON-self-converse (no automorphism realizes T^op). On Z_n (abelian) inversion
always realizes the converse, so circulants are out; the carrier must be NON-ABELIAN
(exactly as Holt needed the metacyclic group of order 27, not a cyclic group). Built an explicit
**F_21 = Z_7 rtimes Z_3 Cayley tournament** (multiplier 2, order 3 mod 7): it is a genuine
tournament, vertex-transitive (21 left translations act regularly), and
**`nx.is_isomorphic(T, T^op)` = False — NON-self-converse**, verified directly. This matches canon
(THM-052, MISTAKE-013): at n=21 there are 88 circulant VT tournaments (all SC) and 22 non-circulant
F_21 VT tournaments (all NS). **n=21 is the smallest order with a vertex-transitive non-self-converse
tournament — the digraph analog of "smallest half-arc-transitive graph = Holt at 27."**

## Where the analogy is tight vs. loose

- **TIGHT:** half-arc-transitivity = "converse Z_2 unrealizable" = the NS (non-self-converse)
  property = the SEA of the merged metagraph G_n/Z_2 (CLAUDE.md). Both require a non-abelian /
  non-circulant carrier. Both split arcs/tilings into TWO classes (= the half, THM-549).
- **LOOSE:** Holt is a *uniform* object (vertex-transitive, |Aut| large), whereas a generic
  NS tournament has trivial Aut. The right tournament match is not a generic NS tournament but a
  **vertex-transitive NS** one (F_21), which is rare. The metagraphs G_n, E_n, G_n/Z_2 are NOT
  vertex-transitive (iso classes have different fibers/H), so they are NOT half-arc-transitive —
  the half-arc structure lives on the *carrier* (F_21), not on the *quotient metagraph*.

## Engineering / structural payoff

- Gives a crisp invariant-theoretic restatement of self-converse vs. NS: **T is self-converse iff
  its "arc-orientation" Z_2 is realized by a relabelling.** For Cayley tournaments this is exactly
  "the group has an anti-automorphism acting as inversion on the connection set."
- Connects the converse-quotient half (half-tiling, THM-549; merged metagraph) to a classical
  algebraic-graph-theory dichotomy (Chen–Quimpo: abelian => arc-transitive), giving an external
  citation for why the SC spine (circulants/Paley) is "too symmetric" and the NS sea is where the
  genuine Z_2-rigidity lives.

## Files
- `04-computation/holt_metacirculant_kpswf5.py` (+ `holt_builder_kpswf5.py` — abelian search, found
  nothing, the expected null result) -> builds/verifies the Holt graph; edge list
  `05-knowledge/results/holt_graph_edges_kpswf5.txt`.
- `04-computation/half_arc_transitive_threadC_kpswf5.py` -> candidate sweep
  (wiggly = hypercube Q_m arc-transitive; Paley arc-transitive).
- `04-computation/half_arc_converse_z2_kpswf5.py` -> the definitive (I)/(II)/(III) analysis.
- `04-computation/f21_ns_tournament_kpswf5.py` -> explicit F_21 VT NS tournament, verified NS.
- Outputs in `05-knowledge/results/*_kpswf5.out`.

## Cross-refs
-> THM-549 (half-tiling = converse-quotient fundamental domain), THM-550, THM-052 (scalar M for
self-converse VT tournaments; F_21 at n=21), MISTAKE-013 (VT does NOT imply SC), CLAUDE.md
(merged metagraph SPINE=SC / SEA=NS), HYP-2657 (doubling order 3 mod 7), HYP-2740 (codex's
Doyle-Holt *sign-rigidity* metaphor on the LRC Delsarte Tanner graph — a DIFFERENT, coding-side
use of the same lead; this HYP is the literal tournament-graph realization).
External: Doyle (1976), Holt (1981); Chen & Quimpo (no half-arc Cayley graph on abelian groups);
Tutte (invariant orientation of a half-transitive graph).
```
```
