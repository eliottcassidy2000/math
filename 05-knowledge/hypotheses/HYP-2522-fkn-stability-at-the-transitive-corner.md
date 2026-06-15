# HYP-2522: FKN Stability at the Transitive Corner Is Local, Not Global

**Status:** OPEN synthesis with exact staircase-cube evidence through `n=7`
(`fkn_near_transitive_tiling_codex.py`, codex-2026-06-15).

## Claim

The right tournament analogue of Friedgut-Kalai-Naor near the staircase tiling
model is **not** that `H(T)` is globally close to a degree-1 dictator on the
full cube. The correct statement is local:

1. the transitive tiling is the free corner;
2. the unique strongest degree-1 escape direction is the extreme tile `(n,1)`,
   i.e. reversing the min-max arc;
3. that flip produces the near-transitive family `H = 2^(n-2) + 1`;
4. genuine interaction is then measured by the top Möbius coefficient on the
   vertex-deletion lattice, not by degree-1 mass alone.

So the dictator shadow exists, but only as a **corner instability direction**.

## Exact evidence

`04-computation/fkn_near_transitive_tiling_codex.py` exhaustively analyzes the
tiling cube for `n=4,5,6,7` using three metrics:

- shell distance from the transitive corner (`d = popcount(bits)`),
- the degree-1 Walsh coefficients of `H` on the staircase cube,
- the top Möbius interaction
  `mu_H(T) = Σ_{U subset V} (-1)^(n-|U|) H(T[U])`.

### 1. Global degree-1 domination fails

The degree-1 Walsh energy share of `H` drops rapidly:

```text
n=4: 0.750000
n=5: 0.490887
n=6: 0.323670
n=7: 0.214557
```

This is anti-FKN in the naive global sense: the first Fourier layer is not
capturing most of the action as `n` grows.

### 2. The extreme tile is the unique dictator-shadow coordinate

The largest degree-1 coefficient is always the extreme tile `(n,1)`:

```text
n=4: tile (4,1), coeff -1.000000
n=5: tile (5,1), coeff -2.000000
n=6: tile (6,1), coeff -4.000000
n=7: tile (7,1), coeff -9.937500
```

And the distance-1 shell is rigid:

```text
flip strip 1  -> H = 3
flip strip 2  -> H = 5
flip strip 3  -> H = 9
flip strip 4  -> H = 17
flip strip 5  -> H = 33
```

More precisely, flipping tile `(a,b)` with strip `a-b-1 = s` produces
`H = 2^s + 1` in the single-flip shell. The unique maximal shell-1 move is the
largest strip `s = n-2`, i.e. `(n,1)`, exactly the near-transitive tournament.

### 3. The interaction defect is parity-sensitive

For the best single flip `(n,1)`, the top Möbius interaction is:

```text
n=4: mu_H = 0
n=5: mu_H = 2
n=6: mu_H = 0
n=7: mu_H = 2
```

So the corner dictator-shadow does not inject the same irreducible interaction
at every `n`; the first visible defect appears to depend on parity.

## Interpretation

The staircase cube has two distinct regimes:

- **global cube:** `H` is a high-interaction observable; degree 1 is too small;
- **transitive corner:** one coordinate dominates the first escape direction,
  namely the min-max flip.

This matches the user's reframe:

```text
board      = mean-field / free transitive corner,
tournament = interacting perturbation,
mu_H       = irreducible deviation from the free assembly of subtournaments.
```

The right FKN-style proof target is therefore:

> among small perturbations of the transitive corner, the near-transitive
> direction is the unique maximizer of first-order Hamiltonian growth.

That is much sharper than saying the whole cube is dictator-like.

## Recursive / deletion-lattice reformulation

The user's `A+B+C-D-E-F+G` language is the Möbius transform on the deletion
lattice. If `f(T)` were assembled freely from lower-order subtournaments, then
its top Möbius coefficient would vanish. The quantity

```text
mu_f(T) = Σ_{U subset V} (-1)^(n-|U|) f(T[U])
```

is therefore the clean recursive interaction defect.

For `f=H`, the transitive corner has `mu_H=0`, while the near-transitive corner
direction shows the first parity-sensitive defect. This is the recursive place
to look for a proof.

## Tournament Analysis / Assumption Challenge

The natural vertices here are **not** runners or cycle lengths. They are:

- tiling coordinates `(a,b)` in the staircase gauge,
- subtournament faces `T[U]` in the deletion lattice,
- and the irreducible interaction coefficient `mu_H`.

This quotient preserves shell distance, degree-1 response, and recursive defect,
while destroying the full labeled tournament. That destruction is intentional:
it isolates the transitive-corner stability mechanism.

Challenged assumption: "FKN for tournaments" should not mean a global junta or
dictator theorem for all tilings. The evidence says the useful analogue is a
**local corner theorem** about the unique dominant single-flip direction.

## Next route

1. Prove the shell-1 law: single flip in strip `s` gives `H = 2^s + 1`.
2. Prove that `(n,1)` uniquely maximizes the absolute degree-1 Walsh coefficient
   of `H`.
3. Explain the parity pattern `mu_H(near-transitive) = 0/2/0/2` and determine
   the exact formula.
4. Extend from one-flip perturbations to low-support perturbations and test
   whether a junta statement holds inside a bounded-radius ball around the
   transitive corner.
