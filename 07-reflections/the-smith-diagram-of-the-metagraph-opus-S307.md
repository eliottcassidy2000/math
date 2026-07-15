# The Smith diagram of the metagraph

*opus-2026-07-14-S307. Owner directive: consider Smith diagrams and squaring
the square against the tournament/metagraph corpus; sweep the niche threads for
connections. Companion data:
04-computation/smith_diagram_of_the_metagraph_opus_S307.py (+ .out).*

## 1. The move: the metagraph IS the network

BSST read a squared rectangle as a resistor network (nodes = horizontal
segments, squares = unit resistors, side = current, KCL = widths fitting, KVL =
heights fitting, aspect ratio = effective resistance, tree counts = the
denominators). S300 applied that dictionary to the LRC deck. The natural home-
project move is one level up: read **G_n itself** as the network — nodes = iso
classes, edges = wiggly (d=1) tile-flips with conductance = flip multiplicity,
source = the transitive class, sink rail = the classes at minimal score-
variance x (the distributed/circulant end), unit current injected. This is the
Smith diagram of the flow of transitivity.

## 2. The exact results (Fractions through n = 6; float64 at n = 7)

| n | nodes | R_n (transitivity resistance) | κ(G_n) (weighted tree complexity) |
|---|---|---|---|
| 4 | 4 | **3/7** ≈ 0.4286 | **21** |
| 5 | 12 | 112209451/437566662 ≈ 0.2564 | 12,993,979,594,752 |
| 6 | 56 | exact (~50-digit fraction) ≈ 0.1090 | exact 111-digit integer (stored) |
| 7 | 456 | ≈ 0.0699 (float) | — |

- R_4 = 3/7 and κ(G_4) = 21 = 3·7: the first rungs are clean; the resistance
  roughly halves per step (0.43, 0.26, 0.11, 0.07) — the metagraph becomes an
  ever-better conductor of transitivity as parallel structure proliferates.
- Total edge-current grows (1.57, 3.72, 4.67, 7.92): the flow spreads over
  more parallel paths even as R_n drops.

## 3. THE CONCORDANCE LAW (the headline)

> **DECIDED NEXT SESSION (S308): the law is a small-n artifact in its strong
> forms.** The 132 n=7 pairs are REAL (exact refinement) — all at adjacent
> levels; at n=8 discordances spread to level-distance 32 (65,921 pairs,
> 0.31%). What survives: LEVEL-MEAN/MEDIAN monotonicity (all computed n) — the
> axis is a mean-field harmonic coordinate. See
> the-concordance-decision-and-the-en-dual-opus-S308.md. The text below is the
> S307 state, kept for the record.

> **Through n = 6, in exact arithmetic, the harmonic potential φ and the
> score-variance axis x induce the SAME strict order on every comparable pair
> of classes: x(A) > x(B) ⟹ φ(A) > φ(B), with zero exceptions
> (5/5, 58/58, 1345/1345 pairs).** At n = 7 (float solve): 92,502/92,634 =
> 99.86%, the 132 discordant pairs pending an exact/high-precision recheck.

The axis the owner proposed for the flow study (S304) — a purely LOCAL score
quantity — turns out to be, empirically exactly, a monotone reparametrization
of the GLOBAL electrical potential of the wiggly layer. The "flow of
transitivity" is not a metaphor: the score-variance spectrum IS the voltage
ladder of G_n. (Named check: decide the 132 n=7 pairs exactly; if real, they
are the first classes where local score data and global network position
disagree — worth naming individually.)

## 4. Where the current flows (spine/ribs/sea, and the Paley bottleneck)

Current split by the class-level edge types:

| n | spine (SC–SC) | ribs (SC–NS) | sea (NS–NS) |
|---|---|---|---|
| 4 | 0.43 | 1.14 | 0.00 |
| 5 | **2.09** | 1.53 | 0.09 |
| 6 | 0.43 | 1.84 | **2.39** |
| 7 | 1.36 | 2.86 | **3.71** |

- The bulk transport migrates from the spine (n = 5) to the SEA as n grows —
  the electrical confirmation of the geometric-alignment picture ("the sea
  dominates at large n"), now with exact currents.
- **But the maximum-current edge is a SPINE edge at every n, and it is the
  last edge into the sink**: at n = 5 the (x=8, H=13)–(x=0, H=15) edge
  (multiplicity 9) carries 0.51 — half the unit current; at n = 7 the
  (x=8, H=159)–(x=0, H=189) edge (multiplicity 135) carries 0.099. And
  H = 15, 189 are the MAXIMUM Hamiltonian-path counts at n = 5, 7 — the sink
  bottleneck is the quasi-regular/**Paley** class. **However diffusely the sea
  carries the current mid-axis, the final approach to the circulant funnels
  through the self-complementary backbone into the Paley node.** (This joins
  the flip-rank/Paley thread: the same object that is the LRC extremal and the
  flip-rank convergent is the electrical sink bottleneck of G_n.)

## 5. The BSST dictionary, entry by entry

| squared square | metagraph |
|---|---|
| rectangle dissected | G_n with source/sink rails (transitive / distributed) |
| unit-resistance square | a wiggly class-edge; conductance = flip multiplicity |
| current through a square | transitivity current on the edge (exact Fractions) |
| aspect ratio = R_eff | the transitivity resistance R_n (3/7 at n = 4) |
| Smith diagram's bus bars | the transitive node (top) and the x_min rail (bottom) |
| matrix-tree denominators | κ(G_n): 21; 12,993,979,594,752; … |
| KVL potentials | φ ≡ the score-variance axis (the concordance law) |
| "perfect" (distinct squares) | distinct currents per edge — holds for the max-current bottlenecks; full distinctness pattern unexamined (named) |
| electrical planar duality | the G_n / E_n duality candidate: E_n (the even-graph metagraph) as the DUAL network — the sharpest untested Smith connection (named for the next session; E_n data exists at n ≤ 7 in the even-graph thread) |

## 6. The archaeology: what existed, what this session fills, and the niche threads

The comprehensive sweep (full report in the S307 session log) found:

- **This computation fills a NAMED gap.** Effective resistance on tournament
  carriers was proposed in the classic-results atlas (item 108,
  tournament-matrix-classic-results-atlas-codex-s210) and executed only on LRC
  objects (the 7-sector cycle with Foster's rule; the Green-current trap
  graphs) — "no effective-resistance / harmonic-potential computation on G_n
  or E_n exists" until now. Spanning trees HAD been computed — the arc-flip
  metagraph's closed form κ(n) = (1/N)Π(2k)^mult (metagraph-spectral-toolkit,
  with THM-588's algebraic connectivity 4) and the merged metagraph's
  unweighted counts (gn_merged_tutte_s219, n ≤ 7) — but not the wiggly
  multiplicity-WEIGHTED complexity computed here, and never with potentials.
- **The repo's existing electrical duality is Ihara/Bass**: cut ⊕ cycle =
  sandpile ⊕ even-graph (three-avenues reflection + HYP-3728: the Bass
  factorization's determinant half is the Laplacian/spanning-tree/cut side —
  our G_n network — and the (1−u²)^{r−1} half is the cycle side = E_n). So the
  BSST "planar duality" row of the dictionary lands on an existing theorem-
  grade structure: **E_n is the cycle half; the Smith network built here is
  the cut half; the missing computation is the E_n-side potential/resistance**
  (E_n data exists: A002854 counts, flip-rank ρ_E = 1,2,6,9, chordality loss
  at E_7, spectrum scripts) — the sharpest next pull.
- **Guardrails inherited**: the tournament↔even-graph map is a cube-iso but
  NOT a structure-iso (HYP-3807); the cut/cycle projection is odd-n only; and
  the deck's KCL story (THM-767 → MISTAKE-146 → THM-771's F = Q + Ω − σ)
  teaches that incidence identities outlive conservation heuristics.
- **Niche threads the Smith frame re-energizes** (from the 20-thread sweep +
  the MPA-01..30 card index): the sandpile/chip-firing observer thread (the
  critical group n^{n−2} as the K_n-level tree count under our class-level
  κ's); the H-saddle/waggly duality and the waggly Schurian coherent-
  configuration framing (association-scheme spectra ⟹ network spectra); the
  line-metagraph rigidity (klein-S163) as the edge-space network; the
  iso-classes-as-free-monoid grading (does R_n factor over SC primes?); the
  flip-rank/Paley thread — now directly connected, since the electrical sink
  bottleneck IS the Paley class; and the Mode-B tower (THM-793): does the
  bundle structure give a recursion R_n = f(R_{n−2}, leg data)? Each is one
  session from a decisive computation; the resistance-through-the-tower and
  the E_n-side potentials are the two I would run first.
