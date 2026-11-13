# Selective returns, Fano wiring, and exact boundary channels

2026-11-13. **Chan–Pak, Hadwiger–Nelson and BSD remain open.**
This session gets a repeatable strong exact return by changing a named
side module, not by applying another common wrap to the whole graph.

## The return

- [Four visits](FOUR_VISIT.md): both Session 320 return pairs have
  strictly ordered cover responses at every nontrivial coordinate.
  This rules out not just four visits, but **every** common strong
  homogeneous continuation of the full returned graphs. The exact
  factored differences make the obstruction transparent.
- [Selective rays](SELECTIVE_RAYS.md): open the quotient and attach
  a named selective module. Two finite bivariate run tables show
  that equality persists for any common tournament in the new
  module. Taking TT_m gives nested strong equal-count pairs
  with H=534717·2^m−1047 and H=4595601·2^m−1915, for every m>=1.
  The second starts from the earlier prime-count return; every term
  is 2 modulo 3. The pairs remain nonisomorphic up to converse
  by Q(2). One-apex exhaustive controls show why simply asking
  for the whole two-slot response to become equal was too strong.
- [Selective side return](SELECTIVE_SECOND_RETURN.md): a different
  common cut triples the prior prime count by replacing a
  transitive side pair with a directed triangle. That *specific*
  side port has precisely those two equalizing modules; the
  unbounded selective rays use different contacts.

## What the merged metagraph must obey

- [Fano geometry](FANO_METAGRAPH.md): the order-21 Cayley block
  has seven canonical directed triangles, an exact Paley–Fano
  majority quotient and a unique dissent arc that carries the
  Frobenius port. Reversing the quotient forces the wrong
  local port action to reverse the triangles. The seven
  independent twists give six converse-merged classes;
  all six are regular, form an independent set in the line
  metagraph, and have six different Hamilton counts.
- [Wedge wiring](WEDGE_WIRING.md): the difference of two aligned
  two-slot wedge codes recovers their agreement support up to
  complement. A literal line is exactly the candidate support
  that is a spanning path and is directed by the local wedges.
  The spine is unique from order five in a given alignment.
  This constrains, but does not classify, all alignment orbits.
- [BSD stencil](BSD_WIRING_STENCIL.md): every aligned line changes
  the binary descent matrix by an explicit alternating form
  built from the marked path Laplacian. Away from the ends and
  short diagonals it reduces to a four-neighbor stencil.
  [BSD parity](BSD_PARITY.md) independently eliminates a selective
  transitive port over F4: the descent nullity sees only the parity
  of the port length. The two equal-H selective rays have distinct
  descent nullities at every length. These are binary/first-descent
  results, not elliptic-curve rank claims.

## Square-sum and coloring boundary relations

[Square overlap ports](SQUARE_OVERLAP_PORTS.md) tracks the native
size-15 to size-25 patterns. The missing size-24 path asks for
one empty virtual-edge channel; size 25 changes the available
ear and splits into channels of sizes four and six. The three
size-23 endpoint pairs cost respectively two, two and three
old-edge deletions to preserve. In the last case the only
two-cut pairing closes a cycle instead of joining the pieces.
This is the same reason the older Hadwiger–Nelson work retains
joint labeled port-color relations rather than a color count.

The remaining Chan–Pak question is not whether a common exact
return can be repeated: it can, on these rays. It is whether
one can choose the side reservoir and its named cut **from the
prescribed target** without losing the exact run relation.
