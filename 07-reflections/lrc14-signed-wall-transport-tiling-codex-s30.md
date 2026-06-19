# LRC14 Signed Wall Transport Tiling

The prompt's "weighted positive and weighted negative" clue feels right only
after refusing to treat those quantities as primitive.

HYP-2642 gives the exact scalar ledger for AP9 versus the endpoint defect:

```text
positive = 9749/158760
negative = 2659/39690
negative - positive = 887/158760.
```

HYP-2643 then says the visible fold count has not changed at all.  AP9 and
`(0,1,2,3,4,5,6,7,9)` both have `12` nontrivial folds; only the target address
changes, with three folds transported from `8` to `9`.

So the picture I now trust is not "there are more bad signs than good signs."
It is:

```text
a common wall tiling has off-diagonal transport;
each transported atom keeps an address;
the signed L_y value is the shadow of that addressed transport.
```

This is exactly the same warning as the tournament arc-flip theorems.  When an
arc flips, the delta is not the number of affected cycles.  The delta lives in
the lost/gained odd cycles through the arc, plus complement/contraction
correction.  When AP9 becomes nearAP9, the wall flip is not just a sign change.
It moves local mass to a neighboring address, and the contraction address is
the fold target/residue/missed-sector bucket.

The tiling model supplies the conservation law.  A common refinement of the two
wall arrangements gives a finite set of atoms.  Every atom starts in one bucket
and ends in another.  That is a transport matrix.  The positive and negative
weights are just what happens after applying the dual valuation `g(N)` and
forgetting the address.  If we forget too early, the proof becomes a fragile
comparison of two scalar sums.  If we keep the address, positive mass can be
paired against negative mass by wall, target, and sign type.

The S30 scout made this concrete.  For the moving endpoint map `8 -> 9`, the
old-sector and new-sector marginals are both exactly uniform:

```text
old sectors: 1/7, 1/7, ..., 1/7
new sectors: 1/7, 1/7, ..., 1/7.
```

The signed effect appears only after this balanced sector transport is filtered
through missed-sector state.  Before weighting, the wall atoms are
`274/2205` positive, `2269/17640` negative, and `4393/5880` neutral.  After
weighting, the scalar shadow returns the exact HYP-2642 surplus
`887/158760`.

The geometry/topology analogy I would keep is a transported hole worldline.  LRC
is not merely asking whether a hole exists; the reduced problem already has
pigeonhole holes.  It asks whether a useful hole reaches the observer after
projection.  AP is a critical confinement state.  The endpoint near-AP is a
small boundary twist that keeps the number of fold events but changes the
address of the transported hole.  The weighted surplus is the boundary integral
of that twist.

The sharp next lemma is therefore not another raw fold-count lemma:

```text
For one-gap endpoint rows F_s=(0,1,2,3,4,5,6,7,7+s),
the addressed wall-transport matrix from AP9 has maximal L_y at s=2.
```

Then the finite near-AP check becomes a structural statement about a transport
matrix.  Rows that do not live in this near-AP transport class should be handed
to the existing Freiman/signed-shell machinery.

Tournament Analysis also changes its vertex set.  I would not use runners or
arcs as the main vertices.  The tournament should rank proof quotients:

```text
signed wall transport
> fold target transport
> arc-flip contraction address
> tiling checksum
> signed pair-sum gauge
> exact-period phi packet
> raw positive/negative totals
> raw runners.
```

The pairwise observable is preservation of LRC-relevant information under
projection.  The winning quotient is the one that still has a finite transport
checksum and still projects to the exact `L_y` difference.  By that standard,
raw positive/negative totals are late-stage shadows, not the object to prove
about first.
