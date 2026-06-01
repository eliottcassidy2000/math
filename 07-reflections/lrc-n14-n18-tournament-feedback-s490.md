# LRC n=14 / n=18 Tournament Feedback S490

The session alternated between `n=14` and `n=18` exactly as requested: press
one until it stuck, inject a noisy outside idea, then carry that noise back to
the other side as a new test.

The outside noise was deliberately modest:

- circular-arc independent sets: no-lonely static masks are independent sets in
  a cycle, hence unsafe gaps are edge covers;
- strong tournaments: a real pressure SCC should carry Hamiltonian-path/cycle
  abundance, not merely one local obstruction;
- zonotope covering radius: LRC can be read as fundamental-domain covering,
  so quotient ladders should be wall-crossing objects.

References used for the noise pass:

- Arthur H. Busch, "A Note on the Number of Hamiltonian Paths in Strong
  Tournaments", Electronic Journal of Combinatorics, 2006:
  <https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3>
- Hsu and Tsai, "Independent Sets in Circular-Arc Graphs", Journal of
  Algorithms, 1995:
  <https://doi.org/10.1006/jagm.1995.1031>
- Henze and Malikiosis, "On the covering radius of lattice zonotopes and its
  relation to view-obstructions and the lonely runner conjecture":
  <https://arxiv.org/abs/1609.01939>

## n=14

The `n=14` hard rows remain the most proof-like.  The S490 scale table
reproduces the familiar debt growth:

```text
d=2:  gap/th=5/264,  unprotected=26
d=7:  gap/th=5/924,  unprotected=84
d=14: gap/th=5/1848, unprotected=168
```

This is the product-depth export story again: shrinking the visible gap does not
erase debt.  It moves debt to deeper rows.

The new Tournament Analysis observation is that the sampled mobile pressure
graphs are still leaf-like:

```text
n14 d=7,d=14: pressure triangles=0, largest pressure SCC=1
```

So the hard rows are not showing a cyclic mobile dependency core.  They are
large-debt, positive-gap, pressure-peelable rows.

That suggests a proof route:

```text
endpoint-private row + pressure leaf -> peel
no pressure leaf -> labelled pressure SCC
labelled pressure SCC -> arithmetic endpoint-cycle contradiction
```

The stuck point is all-times certification.  Sampling selected endpoint/gap
times is not a proof that every chamber peels.  But it gives a sharp target for
the next script: implement pressure-core peeling over the exact endpoint rows.

## n=18

The `n=18` side is more exciting as a discovery machine.  It is `2*3^2`, so it
has both the first-even seam and a square odd-prime payload.  S490 tested
scales `2,3,6,9,18` with the structural skip `8`:

```text
d=2:  gap/th=9/352, unprotected=42
d=3:  gap/th=3/176, unprotected=56
d=6:  gap/th=3/352, unprotected=120
d=9:  gap/th=1/176, unprotected=176
d=18: gap/th=1/352, unprotected=352
```

The `d=18` row matches the n=14 pattern in spirit: shrinking the visible gap by
a gate multiple doubles the endpoint debt.  But the debt now lives in a mixed
`2` and `3^2` product tree.

The surprise is that n=18 also looked pressure-peelable in the sampled
tournament lift:

```text
n18 d=3,d=9,d=18: pressure triangles=0, largest pressure SCC=1
```

This matters for disproof hunting.  If even the richer mixed-torsion row fails
to produce a nontrivial pressure SCC, then n=18 is probably not a near-disproof
machine in the naive sense.  It is more likely a proof laboratory for the
general mixed-prime mechanism.

## Safe-Gap Masks

The circular-arc noise clarified a useful language.  At a static time, a runner
is lonely iff the two adjacent circular gaps around it are both safe.  Thus a
no-lonely configuration has safe gaps forming an independent set in `C_n`, or
equivalently unsafe gaps forming an edge cover.

At witness times S490 saw:

```text
n14 rows: 7 safe gaps, 1-2 adjacent-safe pairs
n18 rows: 9 safe gaps, 1-2 adjacent-safe pairs
```

The numbers are exactly half the circle in both cases.  The current examples
are not close to a no-lonely mask.  They visibly contain adjacent safe pairs.

This suggests a second proof track:

```text
endpoint debt says which walls are unpaid
safe-gap mask says which cycle edges are uncovered
pressure leaves say which blocker can be removed
```

The proof should use all three, not just one scalar gap.

## Current Judgment

`n=14` is still the best proof candidate.  It has the most exact machinery, and
the Tournament Analysis lift keeps showing leafiness rather than cyclic
pressure.

`n=18` is the best forced-random discovery candidate.  It is the smallest
frontier denominator in this range with a genuinely mixed `2*3^2` product-tree
shape.  But after S490, I would no longer call it the best disproof candidate.
The next disproof-like event to search for is very concrete:

```text
pressure_largest_scc > 1
```

If bounded perturbations around the `d=9` and `d=18` rows still never produce
that, then the right n=18 project is a proof certificate, not a counterexample
construction.
