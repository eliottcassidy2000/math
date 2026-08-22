# Visible zero modes as edge-parameter origin-affine circuits

**Status:** structural synthesis for
[THM-3307](../01-canon/theorems/THM-3307-aggregate-origin-affine-criterion-and-six-minimal-visible-kernel-circuits.md).
The proof source is the
[aggregate reduction](../04-computation/gmc_visible_kernel_aggregate_circuit_reduction_20260803.py),
with [frozen output](../05-knowledge/results/gmc_visible_kernel_aggregate_circuit_reduction_20260803.out).
The earlier [independent audit](multi-edge-visible-kernel-and-bordered-resolvent-independent-audit-codex-20260803.md)
remains a genuinely different `2^11` reconstruction.

Everything here concerns finite static symmetric relation-walk graphs.  It is
not a tournament, chronological transition system, FC/GMC mechanism, LRC
mechanism, or asymptotic bit-complexity theorem.

## What changed

The previous exact atlas found six minimal observer-visible zero-mode supports
by exhausting all `2048` subsets of eleven optional row edges.  That result
was reliable but did not explain the number six.

The selector bipartition supplies the missing representation.  In fixed
small/full order the adjacency is

```text
A=[0 N; N^T 0],             N is 14 x 12.                  (1)
```

The all-ones visibility question splits into `ker(N^T)` on the small side and
`ker(N)` on the full side.  An explicit row identity puts the active full-side
observer in `row(N)`, so every full-side kernel vector is invisible.  The
small-side condition is

```text
N^T xi=0,                  1^T xi != 0.                    (2)
```

Normalizing the second relation to `1^T xi=1` says that the origin belongs to
the affine hull of the active small-row neighbourhood vectors.  This is the
right geometric object.  The optional edges are parameters that change those
vectors; a minimal visible edge set is a minimal parameter support enabling
the origin-affine relation.  It is not automatically an ordinary matroid
circuit on the edge set.

## The aggregate compiler

The fourteen coefficients group by their five small selector rows.  Four row
sums and a few distinguished coordinates reduce the twelve column equations
to

```text
b=0,        j a1=0,       d4=-x,       a4=e x,
w c1+D=0,  A+(p-q)x=0,
(1-s)x+rC=0,              (1-u)x=0,
C+(1-v)x=0,               S=x+A+C+D.                       (3)
```

The letters `p,e,w,h,q,j,r,s,k,u,v` are the eleven edge indicators in core
order.  Exact lift formulas prove that `(3)` loses no solution.

When `w=1`, a free `D` aggregate supplies a visible vector for every choice of
the other indicators.  When `w=0`, visibility forces `x != 0`, and division by
`x` leaves the complete Boolean law

```text
u=1,             1-s=r(1-v),             q-p+v != 0.       (4)
```

The six circuits are now a two-line case split:

```text
w;
p,r,u;   p,s,u;   q,r,u;   q,s,u;   s,u,v.                 (5)
```

The indicators `e,h,j,k` can alter births, kernel dimension, mass, and
recurrence order, but they do not enter first visibility.  That is a proved
irrelevance only for the predicate in `(4)`, not for the full sequence.

## Persistence and fragility become algebraic

The `w` circuit has primitive direction of sum/norm `(3,15)`.  Every remaining
edge delta annihilates it, so its Rayleigh contribution gives mass at least
`3/5` on every superset.  This is not merely monotonicity of the predicate;
it is persistence of one fixed vector.

For the fragile circuit `{p,r,u}`, equation `(4)` explains three different
kill moves:

```text
add q: q-p+v becomes 0;
add s: 1-s=r(1-v) fails;
add v: 1-s=r(1-v) fails.                                  (6)
```

Other additions can preserve the primitive vector, rotate the visible
subspace, or introduce the persistent `w` direction.  Thus the state that a
sequence compiler must retain is

```text
minimal support + primitive direction + stabilizer,                       (7)
```

not just the scalar mass or recurrence denominator.

## A two-stage closed-form sequence compiler

THM-3307 and the earlier bordered-resolvent audit now separate two tasks that
were previously mixed:

1. **Head compiler.**  Use `(3)--(4)` to decide whether a zero-eigenvalue
   observer head exists and to recover its primitive direction and mass.
2. **Tail compiler.**  Embed in the fixed 26-node universe and update the base
   resolvent through the local border

   ```text
   G_new(x)-G_base(x)=b+x s(x)^T(C-xK(x))^{-1}s(x).          (8)
   ```

The first stage prevents a rational denominator from silently deleting a
finite zero-mode atom.  The second reuses base powers and inverts only the
changing local interface; every two-edge update in the audited universe has
local dimension at most ten.

This is a concrete answer to “compute the sequence efficiently via a closed
form” for this finite family: compute the exceptional head through the
aggregate incidence presentation, then compile the infinite tail through a
small bordered inverse.  No asymptotic complexity claim follows until a
uniform bound on the border is proved for a growing family.

## Cross-frontier lesson without object identification

The new LRC threshold-layer multicovers and these visible circuits share a
procedure, not a mathematical reduction:

```text
retain one common joint object;
choose the smallest feature support on which failure appears;
derive a normalized dual/affine circuit;
prove every smaller support feasible or invisible;
retain the lost physical or chronological sidecar.                         (9)
```

In LRC the support elements are nested load tails and the witness is a modular
majorant contradiction.  Here the support elements are optional symmetric
relations and the witness is an origin-affine kernel vector.  Their target
predicates and quotient losses are different.  What transfers is the search
discipline: minimize support in the representation where the common object is
still present.

## Next operations

1. Symbolically compile all visibility strata, not only the minimal supports,
   from `(3)` and attach exact projection-mass formulas on each Boolean cell.
2. Extend the bordered inverse from pairs to arbitrary subsets, quotienting a
   persistent kernel direction only after its finite observer head is booked.
3. Search other selector-paired relation systems for a row-span certificate
   like the full-side identity.  Without it, both halves of `(1)` contribute
   and the present reduction does not transfer unchanged.
4. Keep the tournament boundary strict.  An intrinsic orientation has not
   been supplied; the static relation is symmetric and must not be
   tournamentized merely to reuse tournament vocabulary.
