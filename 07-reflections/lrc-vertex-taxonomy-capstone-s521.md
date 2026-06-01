# A taxonomy of LRC-to-tournament vertex-definitions, and the conservation of difficulty (S521 capstone)

*claudebox-2026-06-01-S521, capstone of a long multi-turn exploration of "what can
the vertices of the tournament be?" for the Lonely Runner Conjecture. Organizes
every encoding tried, by faithfulness, restriction, and the deep tool it exposes,
and states the recurring meta-lesson.*

## The question

Encode LRC as "which iso-classes / structures can a speed set exhibit, and is a
designated *lonely* structure always among them?" The art is choosing the VERTICES
so the exhibitable set is meaningfully restricted yet still FAITHFUL (loneliness is
a function of the structure). Two failure modes of restriction recur:
**observer-blindness** (forget the marked vertex) and **single-modulus** (see only
`t=k/n`). Faithful = determine every `||v_i t||` with the observer marked = the
full gap structure.

## The vertex-definitions explored (this session + S521)

| vertices | faithful? | realizable size | deep tool exposed |
|---|---|---|---|
| runners, half-turn (round tournament) | NO (obs-blind) | `A000016(m)` | round/locally-transitive tournaments |
| residues mod q, multiplicative walk | partial | orbit `<= phi(q)` | THM-369 sieve, doubling resonance |
| residue-orbit at n | NO (single-modulus) | `<= phi(n)` | — |
| **sectors (occupancy), chip-firing** | strict | `~` compositions | sandpile / abelian dynamics |
| **sector-VECTOR (cutting sequence)** | strict | poly `~ n*Sum v` in `n^m` | symbolic dynamics, discrepancy |
| **danger graph** | yes (obs-aware) | poly, `<<2^C(n,2)` | perfect unit-circular-arc graph; `chi=omega=`congestion |
| gap-type necklace | yes | `2^n-2`, `#tight<=n-1` | binary necklaces |
| two-gap `{LL,LS,SL,SS}` | yes (minimal) | 4 | THM-384/386 dynamics |
| marked tournament (observer) | yes | round + safety bits | observer = NZ-flow source |
| witness-times `k/(v_i+v_j)` / pairs | yes | `~C(m,2)*`sum | balanced apex pair (Thm A) |
| **bad-arc nerve / winding** | yes (strict) | `Sum v` sub-arcs | nerve lemma / circular covering |
| crossing word (factors) | measure | bounded (periodic) | Rauzy / factor complexity |
| relation lattice `Sum a_i v_i=0` | dual | rank `m-1` | nowhere-zero flow of the speed system |

## What each thread concluded

- **Flow / nowhere-zero flow.** Occupancy = divergence of the crossing-flow;
  **observer lonely <=> observer is a NZ-flow SOURCE** (`indeg=N(t)=0`). The runner
  tournament always carries circulations; LRC = expel the observer from all of
  them. Crisp, but **Menger adds nothing** (source = local degree) and **Tutte
  flow-coloring duality fails** (the danger graph is non-planar at congestion
  `omega>=5`).
- **Coloring.** The danger graph is a **perfect unit-circular-arc graph**,
  `chi=omega=`max congestion in a `1/n`-arc; `chi in {1,..,n}`; **LRC = drive
  `chi` so the observer's color is free**, and the EXTREMAL sets make the danger
  graph **edgeless** at `t=1/n`. But `chi<=n` is automatic, so no chromatic BOUND;
  the content is **arc-covering / congestion-minimization**.
- **Topology (nerve).** The danger sub-arcs cover `S^1` <=> strict-LRC fails;
  covering <=> a winding cycle in the overlap nerve. A clean topological criterion.
- **Symbolic dynamics.** The sector-vector is the **cutting sequence** of the
  torus line (polynomially thin in the `n^m` grid); LRC = the central box appears;
  the crossing word's factor complexity measures orbit complexity (bounded for
  rational lines). Opens discrepancy / bounded-remainder-set tools.
- **Sectors / chip-firing.** `m` tokens hop on an `n`-cycle (fixed circulant
  arena = multiplication; tokens = addition); LRC = clear the observer's two
  forbidden nodes; sandpile invariants are the candidate tool.
- **Witness/pair.** Faithful, restricted to `k/(v_i+v_j)`; the tournament on top
  is inert; the crisp content is `tight <=> ` binding pair sums to exactly `n`.

## The conservation of difficulty (the meta-lesson)

Across every vertex-definition the same pattern holds:

1. The exhibitable set is **polynomially thin** in an exponentially large ambient
   (cutting sequence: `~n*Sum v` of `n^m`; danger graph: `<<2^C(n,2)`; etc.).
2. The LRC criterion is always clean: **a single central object appears** (central
   box / observer-color-free / source / uncovered gap / `LL`).
3. The faithful encodings are all quotients of one finest invariant (the full gap
   structure) and refinements of the minimal one (`{LL,LS,SL,SS}`).
4. **The difficulty is conserved.** Restricting the static structure pushes the
   difficulty into the *dynamics* (how the structure is generated as `t` varies) or
   the *Diophantine data* (the continued-fraction / discrepancy of the line). No
   encoding makes the difficulty vanish; it only changes address.

## Where the difficulty actually lives, and the two surviving tools

Every faithful, restricted encoding reduces LRC to the same core — *the observer's
two forbidden arcs are not perpetually covered* — and the second-moment /
anti-concentration route is RULED OUT (S521). The two routes that survive, both
pointing at the same object:

- **Covering systems.** The discretized covering of `Z/N_*` by generalized APs
  (with the doubling / `2`-resonance structure); Hough / BBMS minimum-modulus
  density bounds are the tool that could forbid it and bound the speeds.
- **Discrepancy / symbolic dynamics.** The cutting sequence of the rational torus
  line and whether its central factor appears; bounded-remainder-set / discrepancy
  bounds for the finite orbit `{v_i a/q}`.

The capstone recommendation of the whole S521 arc: stop hunting for a magic
encoding (difficulty is conserved) and attack the **covering-system reformulation**
directly, using the exact resonance structure (doubling pairs, `n=2*odd` first-even
bridge, binding pair sums to `n`) as the description of the over-correlation a
counterexample must have — the one thing a density theorem can forbid.
