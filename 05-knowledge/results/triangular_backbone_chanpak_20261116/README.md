# Exact Chan–Pak repair: triangle budgets and reflection carries

## Results

- [Triangle budget and exact exclusions](CORE.md). If an induced
  extension strictly increases the count from H, it increases by at
  least 2⌈log₃ H⌉. The proof charges old directed triangles to forced
  odd-cycle packs, then injects Hamiltonian paths into the first
  vertex of each triangle. A fixed half-gap r has a finite induced-pair
  kernel on at most 6r vertices. Exact first three positive gaps:
  +2: 1→3, 3→5; +4: 1→5, 5→9; +6: 3→9, 5→11, 9→15.
  Paired completions of a single unknown arc give all-order hereditary
  obstructions for particular templates, even at the realized count 35.
- [Physical reflection compiler](BALANCE.md). Reflection-even and
  reflection-odd choices split the head/tape count exactly, with no
  mixed response. At tape length 2R+4 there are exactly 2^(R+1)
  masks where a marked pivot increases the count by two. A selected
  subfamily gives 2^(R−1) distinct, ordered, physical pairs plus an
  exact signed-mask inverse. A rigidity lemma explains why matching the simple inner
  signature cannot produce a second constant-scale anchor at
  that same pivot. Its absolute anchors are separated by
  order 2^R and supply Θ(√X) pairs up to X, not a filled interval.
- [Arithmetic inverse and crossovers](NUMBER.md). Linear block
  balances are equivalent to odd-hypotenuse Pythagorean triples.
  The same reflection centers the linear and squared balances at
  2T_n and 4T_n. The 8→9 chain is the first case of a square-layer
  rotation with exact carry on the consecutive-leg Pell ray.
  The Schur microcosm and the 7825 proof isolate the use of a
  forced clash rather than a raw number of constraints.

## What this does and does not solve

This excludes the entire strategy of a uniformly bounded positive
change with the old graph kept induced at arbitrarily high counts.
The physical compiler sidesteps that barrier by changing one
old arc. Its exact low step is solved; its selection from an
arbitrary absolute target is not. The next local-to-global move is
to find another physical bank with fine steps over a radius
comparable to its anchor spacing, or to replace the high bank from
the target. Chan–Pak remains open.

Reproduce with Python 3: check_core.py (all labelled orders through
6, arc-completion examples, seven triples and 6494 active variables),
check_balance.py (literal Hamiltonian path DP, both scaffolds and
binary carry masks), and packing_certificate.py (weighted split).
