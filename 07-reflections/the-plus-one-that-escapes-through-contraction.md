# The +1 that escapes through contraction

**Source:** mac-mini-2026-06-15-S6. Dispatch (human): connect the Alcuin number (= vertex
cover number, or +1) to tournaments via "conflict graph as tournament," Hamiltonicity, and
Kuratowski/Robertson–Seymour. Canon: THM-520, HYP-2553..2555, OPEN-Q-107, T829 (renumbered off a same-dispatch collision with kind-pasteur's THM-519/HYP-2550-2552, whose complementary OCF-Ω-specific Alcuin result this sits beside).

## The shape

The vertex cover number τ is one of the most well-behaved graph parameters under the minor
order: H ≼ G ⟹ τ(H) ≤ τ(G), cleanly, in all three minor cases. So {G : τ(G) ≤ k} is
minor-closed, and Robertson–Seymour hands you a *finite* forbidden-minor obstruction set for
free — the exact same theorem that gives Kuratowski/Wagner's {K₅, K₃,₃} for planarity. τ
sits comfortably inside the well-quasi-ordered universe.

The Alcuin number is τ "plus at most one." That +1 is a tiny correction — a single boat seat.
And yet it is precisely the thing that **breaks minor-closure**. Exhaustively over all 208
graphs on ≤6 vertices: Alcuin is monotone under vertex- and edge-deletion (0 failures), but
it *increases* under edge **contraction** (8 failures). The smallest witness is as clean as
it gets: K₂,₄ has Alcuin 2; contract any edge and Alcuin jumps to 3.

So {G : Alcuin(G) ≤ k} is not minor-closed. The same Robertson–Seymour machinery that
domesticates τ gives *nothing* for Alcuin — no finite minor-obstruction set. The +1 escapes
the WQO world entirely, through the one door (contraction) that deletion can't reach.

## Why contraction, specifically

The mechanism is visible in the witness. In K₂,₄ the two hub vertices live in the same part —
they are non-adjacent. A minimum vertex cover is exactly those two hubs, and because they
don't conflict with each other you can park them together and shuttle the four leaves with
two seats. Contract an edge and the two hubs *merge their adjacencies and become adjacent to
each other*: now the cover is internally conflicting AND each cover vertex fights every leaf.
That is the star's pathology — an over-committed center — and it costs the extra seat.

**Contraction is the only minor operation that can add an edge inside a minimum vertex
cover.** Deletion removes constraints (Alcuin can only ease); contraction can *densify the
cover itself*. The Alcuin +1 is a measure of cover-internal conflict, and cover-internal
conflict is exactly what contraction manufactures. The +1 isn't noise riding on τ — it is
the part of the river-crossing difficulty that is invisible to the minor order's deletion
half and lives entirely in its contraction half.

## The recurring motif

This is the third time this project has met the same skeleton in a month, and naming it is
worth more than any one instance:

> A robust, well-structured invariant, plus a stubborn ±1 correction — and the correction is
> exactly the hard, non-closed, non-generic part.

- **Rédei.** Every tournament has Hamiltonian paths (clean existence); the content is that
  their count is *odd* — a parity the existence proof doesn't see.
- **LRC(14)** (last session). The lonely measure is structurally positive and the best 2025
  methods reach 96.6% of the bound; the conjecture lives in the final ~1–4% sliver — the
  exact value the asymptotic machinery can't close.
- **Alcuin** (this session). τ is minor-closed and finitely obstructed; Alcuin = τ + (the
  bit that breaks minor-closure).

In each case the skeleton is well-quasi-ordered / spectral / asymptotic and *tractable*, and
the correction is parity-like / exact-value / contraction-born and *carries the real
difficulty*. The project keeps rediscovering that the interesting mathematics is not the
clean structure but the one-unit correction sitting on top of it — and that the correction
characteristically escapes whichever closure (minor, absolute-convergence, continuous
relaxation) tames the skeleton.

## What the tournament bridge did and didn't give

The conflict-graph→tournament map (forward arc = edge) is *faithful* on the structure that is
order-and-density combinatorics: independent sets become reverse-transitive runs, cliques
become forward runs, induced cherries become 3-cycles, and Rédei's odd parity falls on every
conflict graph for free. All exact, all verified. But the two creative bridges I hoped for —
that the Alcuin +1 or graph Hamiltonicity would show up as *strong-connectivity* of T_G —
are both false. The +1 and Hamiltonicity live on the cover/scheduling axis, not the cycle
axis; the tournament sees the clique/independent geometry perfectly and the +1 not at all.
That negative is itself the lesson: the +1 is not a tournament-isomorphism invariant. It is
genuinely a property of how a cover sits in the graph — which is why contraction, not any
spectral or cyclic feature, is what moves it.
