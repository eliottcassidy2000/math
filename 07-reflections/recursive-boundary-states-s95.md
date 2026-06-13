# Recursive Boundary States S95

**Session:** codex-2026-05-29-S95
**Script:** `04-computation/recursive_insertion_exposure_s95.py`
**Stored output:** `05-knowledge/results/recursive_insertion_exposure_s95.out`

## The Pattern

A recurring failure mode in the repo is trying to recurse with a scalar invariant.
The recursion usually becomes simple after replacing the scalar by the right
boundary object.

```text
wrong state:  H(T), score(T), endpoint matrix, first invariant
right state:  the interface statistic seen by the next operation
```

This is the same lesson as the tiling-Hamiltonian ratio result:

```text
tiling count is not H alone;
tiling count is H divided by the automorphic boundary tax.
```

and the same lesson as the Cartan/Omega result:

```text
raw symmetric part is inert;
the useful symmetric sector appears after the right lift.
```

## Vertex Insertion Clarification

The strip-recursion tangent suggested that adding a vertex should be computable
from a smaller tournament. The endpoint matrix `M[a,b]` is not enough: it only
remembers where Hamiltonian paths start and end.

The right boundary state is:

```text
start_T[b] = # Hamiltonian paths of T starting at b
end_T[a]   = # Hamiltonian paths of T ending at a
Q_T[a,b]   = # old-vertex permutations where a is immediately followed by b
             and every other consecutive step is a valid edge of T.
```

Call `Q_T` the **bridge-exposure matrix**. It counts broken Hamiltonian spines
with exactly one slot left to repair.

For a new vertex `v`, write:

```text
sig[x] = 1 if v -> x
sig[x] = 0 if x -> v
```

Then:

```text
H(T+v)
  = Σ_{b: sig[b]=1} start_T[b]
  + Σ_{a: sig[a]=0} end_T[a]
  + Σ_{a: sig[a]=0, b: sig[b]=1} Q_T[a,b].
```

This was exhaustively verified for all base tournaments with `n<=5` and every
new-vertex signature.

## Why The Naive Formula Failed

The tempting formula was:

```text
H(T+v) = insert v into Hamiltonian paths of T.
```

That misses paths where `v` repairs an invalid old adjacency. If a path in
`T+v` contains

```text
a -> v -> b
```

then deleting `v` leaves the old consecutive pair `a _ b`. This need not be an
edge `a -> b` of `T`. The new vertex can bridge a non-edge in the old direction.

So the recursive state must count not only exposed valid arcs in Hamiltonian
paths, but also one-break spines. This is the boundary-state correction.

## Exhaustive Data

For base sizes `n=2..5`:

```text
formula failures:
  0, 0, 0, 0

distinct H:
  1, 2, 3, 7

distinct score:
  1, 2, 4, 9

distinct endpoint_M:
  2, 8, 44, 664

distinct hp_exposure_E:
  2, 8, 64, 1024

distinct boundary_Q:
  2, 8, 64, 1024

distinct insertion_response:
  2, 7, 56, 932
```

`H`, score, and endpoint matrix all become insufficient. The bridge-exposure
state gives a direct formula.

The valid-arc HP exposure matrix `E_T[a,b]` is also extremely strong at small
`n`; by `n=4` it distinguishes every labeled tournament in the test range. That
is useful as a fingerprint, but the insertion formula itself wants `Q_T`.

## Cross-Thread Translation

### Deletion-Contraction

THM-082 says:

```text
H(D) = H(D\e) + H(D/e).
```

The contraction term is exactly a boundary state: it records paths where the
chosen edge is forced as a glued interface. The theorem is not "H is recursive"
in the scalar sense; it is "H is recursive after the edge interface is exposed."

### Tiling Fibers

The explorer fixed-path fiber obeys:

```text
F(C) = H(C)/|Aut(C)|.
```

Again, `H` alone is not the recursive quotient mass. The missing boundary object
is the stabilizer of the class under relabeling. Automorphisms are the boundary
tax paid when moving from labeled Hamiltonian paths to one fixed Szele slice.

### Omega / OCF

`H = I(Omega,2)` is another boundary-state move. The directed tournament is hard
to recurse on directly; the undirected conflict graph exposes the intersections
between odd cycles. The scalar `t3` fails inside residual fibers; `Omega`
separates them.

### Path Homology

The LES/deletion proofs for Betti constraints keep saying: vertex deletions are
not enough unless the relative group is controlled. The relative group is the
boundary state. The proof succeeds when the interface `H_*(T,T\v)` is small or
structured.

### Transfer Matrices

Endpoint transfer matrices are terminal boundary states. For vertex insertion,
they are too terminal: they forget the interior slots where a new vertex can be
inserted. The next transfer object should be bridge exposure, not just endpoints.

## Meta-Hypothesis

**Boundary-State Principle.**
Whenever a tournament invariant appears not to recurse, the obstruction is
usually that the proposed state forgot the boundary through which the next
operation acts. Add the smallest interface statistic that sees that operation,
and the recursion becomes linear or a cut sum.

Examples:

```text
fixed-path tilings:     H              -> H/Aut
vertex insertion:       endpoint M     -> bridge exposure Q
OCF residual fibers:    t3             -> Omega independent sets
path homology:          deletion Betti  -> relative LES group
flip differences:       H(T),H(T')      -> contracted minors
```

## Next Tests

- Try to update `Q_T` itself recursively under vertex insertion. If possible,
  the strip recursion becomes a true dynamic program.
- Compare `Q_T` against transfer matrix Walsh layers: the one-break spines may
  be the missing degree-one boundary of `M`.
- Search whether `Q_T` descends well to isomorphism classes or whether it is a
  labeled-state object like the tiling staircase.
- Test whether `Q_T` or valid-arc exposure `E_T` predicts `H` residuals inside
  fixed score or fixed `t3` fibers.

## Working Slogan

Recursion is not scalar. Recursion lives on the boundary.
