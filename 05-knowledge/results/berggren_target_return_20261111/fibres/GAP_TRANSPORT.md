# Named gap transport: a repeatable neutral repair with two exact costs

**PROVED for all triangle-chain sources. FINITE-EXACT sharp first-source
peak under the stated floor.** These are actual old-arc changes, not claims
that the macro is a native path-complement edge. The source and final
target have the same H; intermediate values are stated explicitly.

## A local macro on actual vertices

Suppose consecutive SCC blocks are a directed triangle A and singleton v,
with A->v. Keep every name, all other arcs, and the triangle's internal
orientation. Reverse the three arcs between A and v in any order. The
four-vertex block has count3,5,5,3 through the operation. Its final form
is v->A. Thus in any larger ordinal triangle-chain source with r triangles,
the full trace is

```
3^r -> 5*3^(r-1) -> 5*3^(r-1) -> 3^r.                  (1)
```

The mixed intermediate block is S4, the strong four-vertex H5 core. No
count formula for arbitrary overlapping modules is being assumed: the
entire four-vertex block is uniformly seen by every outside SCC, so H
multiplies over this ordinal join at each step.

There is a second path to the SAME named endpoint. With A labelled
a0->a1->a2->a0, reverse a2->a0 to make TT3; move v left past a2,a1,a0
by three adjacent transitive swaps; restore a2->a0. This five-flip trace is

```
3^r -> 3^(r-1) -> 3^(r-1) -> 3^(r-1) -> 3^(r-1) -> 3^r. (2)
```

Every original vertex survives, the final triangle has its original arcs,
and both paths end at the identical adjacency matrix. The first path
preserves the lower bound H>=H(source); the second never exceeds
H(source). Confusing those two constraints would create a false barrier.

## Complete normalization of every gap vector

For T(g) of TRIANGLE_FILTRATION.md, repeatedly interchange any singleton
with the triangle immediately before it. This bubble process terminates
at the actual named graph

```
TT_(sum g_i) -> C3^r.
```

It performs exactly `sum_i i*g_i` macros, because a singleton originally
in gap i crosses precisely i triangles. Equations (1) or (2) can be
chosen independently at every crossing. This proves a uniform,
repeatable same-count repair of the whole infinite family, with complete
vertex retention and a known count envelope throughout.

The endpoint has the full known extension filtration
`C_k=3^max(r-2k,0)` and lifetime ceil(r/2). Thus the macro converts every
gap pattern to the smallest possible lifetime allowed by the untouched-
triangle lower bound. The current H and every triangle's original names
survive. Root indegree, cover response and adjacency to other components
need not survive, so this is not automatically a prepared-root return for
the separate `3H+2` / `3H+4` program.

## Minimal source and the exact floor boundary

At the seven-vertex source C3->K1->C3, H9 and C_1=3. Moving its middle
singleton left gives K1->C3->C3, also H9 but C_1=1. The three-flip route
has trace9,15,15,9.

Under the constraint that EVERY intermediate graph have H>=9, peak15 is
optimal for any nontrivial arc-flip path from this exact source. Complete
enumeration of its21 single-arc neighbours shows that every neighbour
of H>=9 has H>=15. The displayed path attains that bound.

Without the floor, the five-flip path has trace9,3,3,3,3,9 and peak9.
The theorem therefore proves no unqualified peak15 obstruction. This
explicit alternative is a hostile control for the scope of the optimum.

## Reproduction and composition contract

Run `python -B gap_transport.py`. It checks all56 gap vectors with r1..3,
each gap<=2 and total singleton count<=3, including both complete named
paths for every crossing, the resulting adjacency matrix, and the exact
number of macros. The21-neighbour minimum-source certificate is stored
in `gap_transport.json`. The universal statements follow from the local
four-vertex block and the inversion count, not from those finite tests.

This links the small H3/H5 cores to an infinite family of neutral source
repairs and to the extension-divisibility invariant. It changes a known
arithmetic obstruction while keeping the target count. It cannot on its
own leave the scalar spectrum {3^r}; global Chan–Pak coverage remains OPEN.
