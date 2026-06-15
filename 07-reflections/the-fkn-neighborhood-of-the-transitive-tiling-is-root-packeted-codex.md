# The FKN neighborhood of the transitive tiling is root-packeted: the dictator atoms are interval roots, not anonymous bits

**Source:** codex-2026-06-15-S11. Extends the live FKN tiling-cube thread
(`T823`, THM-511, upstream `fkn_tiling_cube_macmini_0615s3.py`) by resolving the
first two Hamming shells around the transitive ground state.

## The sharpened local picture

The new upstream checkpoint framed the tiling cube as an FKN world:

- the all-aligned point is the **transitive ground state**;
- Hamming distance from it is the backward-arc **energy**;
- level-1 / single-flip content is the **dictator** end of the story.

That is right, but locally it is still too coarse. The new exact correction is:

> **The free dictator atoms are the non-simple interval roots of `A_{n-1}`.**
> A one-flip tiling is not just "one bit turned on"; it is a specific interval
> upset `(x,y)`, i.e. a score transfer `e_y-e_x` along the chain.

So the radius-1 shell is not a homogeneous Boolean sphere. It is a triangular
atlas of interval roots.

## Shell 1: all atoms are different, but H and c3 see only height

THM-513 makes the shell-1 geometry exact.

For a single flipped tile `(x,y)` with `g=x-y-1` intermediate vertices:

```text
score defect = e_y - e_x,
c3           = g,
H            = 2^g + 1.
```

Three things happen at once.

1. **The full iso class remembers the interval address.**
   The one-flip shell contains exactly `m=C(n-1,2)` pairwise non-isomorphic
   tournaments. At `n=6`, the energy-1 shell has `10` distinct iso classes, not
   `4` or `1`.

2. **The score remembers the root exactly.**
   A one-flip move is a unit transfer from `x` down to `y`.

3. **`H` and `c3` forget the address and remember only the gap.**
   All shell-1 atoms with the same `x-y-1` share the same `(H,c3)`.

This is a precise version of "the naive `2^m` count ignores structure." Even at
radius `1`, different invariants collapse different parts of the root packet.

## Shell 2: the first nonlinear layer is packet incidence

At radius `2`, the first cyclic interaction is already exact and clean.

Let the two active tiles have gaps `g1,g2`.

```text
same-end  -> c3 = g1 + g2 - 1
cross-end -> c3 = g1 + g2 + 1
disjoint  -> c3 = g1 + g2
```

This is the packet law:

- **same-end** packets lose one relay triangle because the second reversal
  destroys exactly one of the first packet's triangles;
- **cross-end** packets create one extra relay triangle;
- **disjoint** packets simply add.

So shell `2` is not organized by energy alone either. It is organized by the
incidence graph of the two interval roots.

The Hamiltonian-path side reads the same way at second order. Existing quadratic
theorems already say:

- same-end packets are subtractive for `H`;
- cross-end packets are constructive;
- disjoint packets are constructive in the stored data through `n=8`.

That is the first packet-interaction correction to the pure dictator picture.

## The recursive answer to the user's prompt

The user asked for an abstract recursive view, not just the raw
`2^(n choose 2)` count. The right recursion is now visible.

- Going from `n` to `n+1` adds one new top row of interval roots `(n+1,y)`.
- Old coefficient data is preserved on old tiles (THM-299).
- The local neighborhood therefore grows by **one triangular root layer at a
  time**.

So a tournament near the transitive state is best read as a packet of interval
upsets layered onto the ranking, not as an arbitrary Boolean word.

## The FKN reading after the correction

FKN says "most mass on level 1" means "close to a dictator." In this chamber the
correct refinement is:

```text
close to level 1
  means
close to a sparse packet of interval-root upsets on the transitive ranking.
```

The dictators are the interval roots. The first question is not only "how many
bits are on?" but also:

- which intervals were flipped,
- what their gap heights are,
- and how those intervals meet.

That gives the local metric package:

1. energy `k`,
2. height mass `sum g_i`,
3. same-end / cross-end / disjoint packet counts,
4. higher relay / overlap faces.

Those are the metrics that illuminate the pattern structure of tournaments near
the ranking end.

## The bridge back to Möbius / inclusion-exclusion

The user also asked for the `A+B+C-D-E-F+G` viewpoint. Upstream T823 already
identified that as the Möbius extraction of irreducible content. The local
packet picture says what that irreducible content is near the ground state:

- shell `1`: a single interval root;
- shell `2`: pair incidence;
- shell `3+`: the first genuine relay / overlap packets.

So the Möbius sieve and the Walsh/FKN grading are two views of the same local
hierarchy:

```text
root atom -> root packet -> higher packet.
```

That is the concrete extension of the user's intuition.
